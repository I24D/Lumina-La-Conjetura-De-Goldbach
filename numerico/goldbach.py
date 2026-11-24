#!/usr/bin/env python3
"""
goldbach.py — herramienta numérica del proyecto Goldbach-LUMINA.

ESTADO EPISTEMOLÓGICO: todo lo que produce este script es EVIDENCIA NUMÉRICA,
NO demostración. La conjetura fuerte ya está verificada hasta 4·10^18
(Oliveira e Silva, Herzog y Pardi, 2014). Este script no compite con ese
récord: es una herramienta reproducible para generar datos (p_min(n), r(n),
comparación con la predicción de Hardy–Littlewood) y poner a prueba heurísticas.

Subcomandos:
  selftest               pruebas internas contra fuerza bruta
  verify --hasta B       comprueba todos los pares 4 <= n <= B y registra
                         los récords de p_min(n), el menor primo p con n-p primo
  reps --hasta M         r(n) exacto para n <= M (por FFT) y comparación
                         con Hardy–Littlewood

Cada ejecución guarda un JSON con su procedencia (SHA-256 del script,
versiones, parámetros, tiempos) en numerico/resultados/.
"""
from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
import math
import platform
import sys
import time
from pathlib import Path

import numpy as np

AQUI = Path(__file__).resolve().parent
RESULTADOS = AQUI / "resultados"
C2 = 0.66016181584686957392781211  # constante de los primos gemelos


# ---------------------------------------------------------------- cribas
def criba_simple(n: int) -> np.ndarray:
    """Array booleano es_primo[0..n]."""
    s = np.ones(n + 1, dtype=bool)
    s[:2] = False
    for i in range(2, math.isqrt(n) + 1):
        if s[i]:
            s[i * i :: i] = False
    return s


def primos_hasta(n: int) -> np.ndarray:
    return np.nonzero(criba_simple(n))[0].astype(np.int64)


def criba_segmento(lo: int, hi: int, base: np.ndarray) -> np.ndarray:
    """es_primo para los enteros de [lo, hi). 'base' debe contener todos
    los primos <= sqrt(hi-1)."""
    seg = np.ones(hi - lo, dtype=bool)
    for p in base:
        p = int(p)
        pp = p * p
        if pp >= hi:
            break
        inicio = max(pp, -(-lo // p) * p)
        seg[inicio - lo :: p] = False
    for k in (0, 1):
        if lo <= k < hi:
            seg[k - lo] = False
    return seg


# ------------------------------------------- Miller–Rabin determinista
_BASES_MR = (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37)  # válido para n < 3.3e24


def es_primo(n: int) -> bool:
    if n < 2:
        return False
    for p in _BASES_MR:
        if n % p == 0:
            return n == p
    d, s = n - 1, 0
    while d % 2 == 0:
        d //= 2
        s += 1
    for a in _BASES_MR:
        x = pow(a, d, n)
        if x in (1, n - 1):
            continue
        for _ in range(s - 1):
            x = x * x % n
            if x == n - 1:
                break
        else:
            return False
    return True


# ------------------------------------------------------- núcleo verify
def pmin_bloque(a0: int, a1: int, base: np.ndarray, impares: np.ndarray, pmax: int):
    """Para cada par n de [a0, a1) (a0 par y >= 6) calcula p_min(n), el menor
    primo impar p con n-p primo, buscando entre p <= pmax.
    Devuelve (n_arr, pmin, n_no_cubiertos); pmin vale 0 donde no se encontró."""
    lo = max(0, a0 - pmax)
    isp = criba_segmento(lo, a1, base)
    n_arr = np.arange(a0, a1, 2, dtype=np.int64)
    pmin = np.zeros(n_arr.size, dtype=np.int32)
    resto = np.arange(n_arr.size)
    for p in impares:
        if resto.size == 0:
            break
        v = n_arr[resto] - p
        ok = v >= 2
        acierto = np.zeros(resto.size, dtype=bool)
        acierto[ok] = isp[v[ok] - lo]
        pmin[resto[acierto]] = p
        resto = resto[~acierto]
    return n_arr, pmin, n_arr[resto]


def pmin_ingenuo(n: int) -> int:
    for p in range(2, n // 2 + 1):
        if es_primo(p) and es_primo(n - p):
            return p
    return 0


# --------------------------------------------------------- utilidades
def sha256_de(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def procedencia(params: dict) -> dict:
    return {
        "script": str(Path(__file__).name),
        "script_sha256": sha256_de(Path(__file__)),
        "python": sys.version.split()[0],
        "numpy": np.__version__,
        "plataforma": platform.platform(),
        "params": params,
        "inicio_utc": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
        "estado_epistemologico": "EVIDENCIA NUMÉRICA — no es demostración",
    }


def guardar(nombre: str, datos: dict) -> Path:
    RESULTADOS.mkdir(parents=True, exist_ok=True)
    sello = dt.datetime.now().strftime("%Y%m%d_%H%M%S")
    ruta = RESULTADOS / f"{sello}_{nombre}.json"
    ruta.write_text(json.dumps(datos, indent=2, ensure_ascii=False), encoding="utf-8")
    return ruta


# ------------------------------------------------------------ comandos
def cmd_selftest(_args) -> int:
    fallos = 0
    rng = np.random.default_rng(12345)

    # 1. Criba segmentada frente a criba simple
    ref = criba_simple(300_000)
    base = primos_hasta(600)
    for _ in range(200):
        lo = int(rng.integers(0, 299_000))
        hi = lo + int(rng.integers(1, 1000))
        if not np.array_equal(criba_segmento(lo, hi, base), ref[lo:hi]):
            fallos += 1
            print(f"FALLO criba_segmento [{lo},{hi})")

    # 2. Miller–Rabin frente a la criba
    for n in range(0, 300_001):
        if es_primo(n) != bool(ref[n]):
            fallos += 1
            print(f"FALLO es_primo({n})")
            break

    # 3. p_min por bloques frente a fuerza bruta, con bloques pequeños
    #    para ejercitar los bordes entre segmentos
    impares = primos_hasta(5000)[1:]
    base = primos_hasta(200)
    a = 6
    while a < 20_002:
        b = min(a + 998, 20_002)
        n_arr, pmin, nc = pmin_bloque(a, b, base, impares, 5000)
        if nc.size:
            fallos += 1
            print(f"FALLO: sin cubrir {nc[:5]}")
        for n, pm in zip(n_arr.tolist(), pmin.tolist()):
            if pm != pmin_ingenuo(n):
                fallos += 1
                print(f"FALLO p_min({n}) = {pm}, ingenuo = {pmin_ingenuo(n)}")
                break
        a = b

    # 4. r(100) = 6 (3+97, 11+89, 17+83, 29+71, 41+59, 47+53)
    r = representaciones(1000)
    if r[100] != 6:
        fallos += 1
        print(f"FALLO r(100) = {r[100]}")
    # 5. r(n) por FFT frente a fuerza bruta para n <= 1000
    for n in range(4, 1001, 2):
        bruto = sum(1 for p in range(2, n // 2 + 1) if es_primo(p) and es_primo(n - p))
        if r[n] != bruto:
            fallos += 1
            print(f"FALLO r({n}) = {r[n]}, bruto = {bruto}")
            break

    print("SELFTEST:", "OK" if fallos == 0 else f"{fallos} FALLOS")
    return 0 if fallos == 0 else 1


def cmd_verify(args) -> int:
    B = int(args.hasta)
    L = int(args.bloque) // 2 * 2
    pmax = int(args.pmax)
    prov = procedencia({"hasta": B, "bloque": L, "pmax": pmax})
    t0 = time.time()
    base = primos_hasta(math.isqrt(B) + 2)
    impares = primos_hasta(pmax)[1:]

    records = [{"n": 4, "p_min": 2}]  # 4 = 2 + 2 es el único caso con p = 2
    max_global = 2
    contados = 1
    suma_pmin = 2
    hist = {}
    contraejemplos = []
    a = 6
    ultimo_aviso = t0
    while a <= B:
        b = min(a + L, B + 1)
        b += b % 2  # extremo derecho par, así el bloque incluye B si es par
        b = min(b, B + 2)
        n_arr, pmin, no_cub = pmin_bloque(a, b, base, impares, pmax)
        # los pares no cubiertos con p <= pmax se revisan uno a uno sin límite
        for n in no_cub.tolist():
            pm = pmin_ingenuo(n)
            idx = (n - a) // 2
            pmin[idx] = pm
            if pm == 0:
                contraejemplos.append(n)
                print(f"!!! CONTRAEJEMPLO CANDIDATO: n = {n}")
        cm = np.maximum.accumulate(pmin)
        previo = np.empty_like(cm)
        previo[0] = max_global
        previo[1:] = np.maximum(cm[:-1], max_global)
        for i in np.nonzero(pmin > previo)[0].tolist():
            records.append({"n": int(n_arr[i]), "p_min": int(pmin[i])})
        max_global = max(max_global, int(cm[-1]))
        contados += n_arr.size
        suma_pmin += int(pmin.astype(np.int64).sum())
        vals, cnts = np.unique(pmin, return_counts=True)
        for v, c in zip(vals.tolist(), cnts.tolist()):
            hist[v] = hist.get(v, 0) + c
        a = b
        if time.time() - ultimo_aviso > 30:
            ultimo_aviso = time.time()
            print(f"  ... {a:,} ({100 * a / B:.1f} %), max p_min = {max_global}, "
                  f"{time.time() - t0:.0f} s", flush=True)

    t = time.time() - t0
    esperado = (B - 4) // 2 + 1
    datos = prov | {
        "fin_utc": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
        "segundos": round(t, 2),
        "pares_comprobados": contados,
        "pares_esperados": esperado,
        "contraejemplos": contraejemplos,
        "max_p_min": max_global,
        "huella_suma_p_min": suma_pmin,
        "records_p_min": records,
        "histograma_p_min": {str(k): v for k, v in sorted(hist.items())},
    }
    ruta = guardar(f"verify_{B}", datos)
    ok = (contados == esperado) and not contraejemplos
    print(f"verify hasta {B:,}: {contados:,} pares ({'OK' if ok else 'REVISAR'}), "
          f"max p_min = {max_global}, {t:.1f} s -> {ruta.name}")
    print("records:", ", ".join(f"{r['n']}:{r['p_min']}" for r in records))
    return 0 if ok else 2


def representaciones(M: int) -> np.ndarray:
    """r[n] = número de pares no ordenados {p, q} de primos con p + q = n,
    para 0 <= n <= M."""
    isp = criba_simple(M)
    f = isp.astype(np.float64)
    tam = 1 << (2 * M + 1).bit_length()
    F = np.fft.rfft(f, tam)
    conv = np.fft.irfft(F * F, tam)[: M + 1]
    r_ord = np.rint(conv).astype(np.int64)
    mitad = np.zeros(M + 1, dtype=np.int64)
    mitad[0::2] = isp[: M // 2 + 1].astype(np.int64)[: mitad[0::2].size]
    return (r_ord + mitad) // 2


def serie_singular(n: int) -> float:
    """S(n) = 2*C2 * prod_{p | n, p > 2} (p-1)/(p-2), para n par."""
    s = 2 * C2
    m = n
    while m % 2 == 0:
        m //= 2
    p = 3
    while p * p <= m:
        if m % p == 0:
            s *= (p - 1) / (p - 2)
            while m % p == 0:
                m //= p
        p += 2
    if m > 1:
        s *= (m - 1) / (m - 2)
    return s


def prediccion_hl(n: int) -> float:
    """Predicción de Hardy–Littlewood para los pares NO ordenados:
    (1/2) * S(n) * integral_2^{n-2} dt / (log t * log(n-t))."""
    t = np.linspace(2.0, n / 2.0, 200_001)
    integ = 2.0 * float(np.trapezoid(1.0 / (np.log(t) * np.log(n - t)), t))
    return 0.5 * serie_singular(n) * integ


def cmd_reps(args) -> int:
    M = int(args.hasta)
    prov = procedencia({"hasta": M})
    t0 = time.time()
    r = representaciones(M)
    pares = np.arange(4, M + 1, 2)
    cero = pares[r[4::2] == 0]
    # muestra logarítmica de n para comparar con Hardy–Littlewood
    muestra = sorted({int(x) // 2 * 2 for x in np.geomspace(100, M, 60)} | {M // 2 * 2})
    comp = []
    for n in muestra:
        pred = prediccion_hl(n)
        comp.append({"n": n, "r": int(r[n]), "pred_HL": round(pred, 3),
                     "cociente": round(r[n] / pred, 5)})
    # el "margen" por tramos: min r(n)/pred_HL en [2^k, 2^(k+1)), muestreado
    margenes = []
    k = 10
    while (1 << (k + 1)) <= M:
        lo, hi = 1 << k, 1 << (k + 1)
        ns = np.arange(lo, hi, 2)
        if ns.size > 400:
            ns = ns[np.linspace(0, ns.size - 1, 400).astype(int)]
        cocs = [(r[n] / prediccion_hl(int(n)), int(n)) for n in ns]
        mn = min(cocs)
        mx = max(cocs)
        margenes.append({"tramo": f"[2^{k}, 2^{k+1})", "min_cociente": round(mn[0], 4),
                         "en_n": mn[1], "max_cociente": round(mx[0], 4)})
        k += 1
    t = time.time() - t0
    datos = prov | {
        "fin_utc": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
        "segundos": round(t, 2),
        "pares_con_r_cero": cero.tolist(),
        "comparacion_HL": comp,
        "margen_por_tramos_muestreado": margenes,
    }
    ruta = guardar(f"reps_{M}", datos)
    print(f"reps hasta {M:,}: pares con r(n)=0: {cero.tolist()} ; {t:.1f} s -> {ruta.name}")
    for c in comp[-8:]:
        print(f"  n={c['n']:>12,}  r={c['r']:>9,}  HL={c['pred_HL']:>12,.1f}  r/HL={c['cociente']}")
    for m in margenes[-5:]:
        print(f"  {m['tramo']}: min r/HL = {m['min_cociente']} (n={m['en_n']}), max = {m['max_cociente']}")
    return 0 if cero.size == 0 else 2


def cmd_contrastar(args) -> int:
    """Compara los récords de p_min de un JSON de 'verify' con OEIS A025018/A025019."""
    datos = json.loads(Path(args.json).read_text(encoding="utf-8"))
    ref = json.loads((AQUI / "referencia_oeis_pmin.json").read_text(encoding="utf-8"))["records"]
    B = int(datos["params"]["hasta"])
    esperados = [tuple(r) for r in ref if r[0] <= B]
    obtenidos = [(r["n"], r["p_min"]) for r in datos["records_p_min"]]
    if max(r[0] for r in ref) < B:
        print(f"AVISO: la referencia solo llega a n = {max(r[0] for r in ref):,}; "
              f"los récords por encima no se pueden contrastar")
        obtenidos_cmp = [r for r in obtenidos if r[0] <= max(x[0] for x in ref)]
    else:
        obtenidos_cmp = obtenidos
    ok = obtenidos_cmp == esperados
    print(f"{len(obtenidos_cmp)} récords obtenidos frente a {len(esperados)} de OEIS hasta {B:,}: "
          f"{'COINCIDEN' if ok else 'DIFIEREN'}")
    if not ok:
        for a, b in zip(obtenidos_cmp, esperados):
            if a != b:
                print(f"  primera diferencia: obtenido {a}, OEIS {b}")
                break
    return 0 if ok else 3


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("selftest").set_defaults(f=cmd_selftest)
    v = sub.add_parser("verify")
    v.add_argument("--hasta", type=float, required=True)
    v.add_argument("--bloque", type=float, default=2e7)
    v.add_argument("--pmax", type=int, default=20000)
    v.set_defaults(f=cmd_verify)
    rp = sub.add_parser("reps")
    rp.add_argument("--hasta", type=float, required=True)
    rp.set_defaults(f=cmd_reps)
    c = sub.add_parser("contrastar")
    c.add_argument("json", help="archivo JSON producido por 'verify'")
    c.set_defaults(f=cmd_contrastar)
    args = ap.parse_args()
    return args.f(args)


if __name__ == "__main__":
    sys.exit(main())
