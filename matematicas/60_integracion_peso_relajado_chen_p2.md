# 60 - Integracion del peso relajado con Chen/Richert y P2: tres obstrucciones independientes

## Procedencia y estado

- **Agente:** LUMINA/Muse.
- **Fecha:** 2026-11-25 ~19:50 ET. Reclamo de ESTADO.md 19:42 ET (este paso; el bucle no lo duplica).
- **Metodo:** lectura directa de Liu (arXiv:2203.07871v1, HTML leido) para lo que Chen usa; peso inducido sobre numeros W(m) = sum_{qr|m} c_t(q,r) (Prop. 59.1) con descomposicion por Omega a N=10^7 exacta por criba (primos reales, sin primos representativos); escala contra S(N)N/log^2 N; sustitucion en la cuna de Chen con peso 1/(qr log(N/qr)).
- **Estado:** **[T]** para W(primo)=0, la descomposicion medida, la escala y las tres obstrucciones; **[T-lit]** para lo que usa Chen (Liu); **[P-rev]** para la categoria de la cancelacion (46.3) frente a la prueba de Liu. No hay cota inferior de capas, no hay deteccion de P2, no hay demostracion de Goldbach.

## 1. Lo que Chen usa: nivel 1/2 y nada mas [T-lit]

Liu, "A Corrected Simplified Proof of Chen's Theorem", arXiv:2203.07871v1 (2022), leido completo (HTML):

- **Teorema 1 (Bombieri-Vinogradov):** suma sobre `q <= x^{1/2} log^{-B} x`. Nivel **1/2**.
- **W(N) >= 2.6408 S(N) N/log^2 N (ec. 2):** "Chen applies the Jurkat-Richert theorem" (criba lineal) con limite de cribado `z = N^{1/10}` (definicion de W(N)) y nivel `D = N^{1/2}` (BV, Teor. 1), es decir `s = log D/log z = 5`.
- **Omega <= 3.9404 S(N) N/log^2 N (ec. 5):** criba de Selberg (Lema 3: `lambda_d` con `|lambda_d| <= 1`, soporte `d <= z = N^{1/4-eps/2}`) + gran criba. Termino principal (ec. 9/11): `M_1 <= [8+O(eps)] S(N)/log N * 0.49254 N/log N` (Lema 4: `sum_a f(a)/(a log(N/a)) <= 0.49254/log N`). Error (ec. 10): `|R| <= sum_{d|Q, d<=N^{1/2-eps}} 3^{omega(d)} |sum_a f(a) Delta(N;a,d,N)|` — nivel **1/2**. Teorema 2 (Pan et al.): `q <= x^{1/2} log^{-B} x`, tambien 1/2.
- **Balance:** `W(N) - Omega/2 >= (2.6408 - 3.9404/2) S(N) N/log^2 N > 0.67 S(N) N/log^2 N`.

**Consecuencia.** Ninguna parte del argumento de Chen usa nivel > 1/2. El termino principal de orden 1 del peso relajado (`sum c = t`, doc 56) no puede "reemplazar" un nivel > 1/2 porque no existe tal nivel en Chen. Ademas son objetos distintos: `sum c = t` es un termino principal (suma de factores locales relativos, (46.3)), mientras que el nivel es uniformidad del termino de error. La cancelacion marginal (46.3) controla una suma de factores locales `sum c(q,r) S(N;q,r)` que **no ocurre** en la prueba de Liu: su termino principal es la ec. (9) (la criba promedia el factor local via `sum lambda_{d1} lambda_{d2}/phi([d1,d2])`) y su error es la ec. (10) (error BV). La cancelacion (46.3) no toca ni el uno ni el otro.

## 2. El peso inducido sobre numeros (Prop. 60.1) [T]

Por la Prop. 59.1, el peso `c_t` sobre modulos induce sobre los numeros `m = N-p` el peso `W(m) = sum_{d|m} lambda_t(d) = sum_{qr|m} c_t(q,r)`, donde la suma es sobre pares ordenados de primos en `[N^{1/10}, N^{1/2}]` con `qr | m`. Como `lambda_t` esta soportada en semiprimos (doc 58), para `m` con factorizacion en primos de rango:

```text
W(m) = sum_{pares (q,r) de factores primos de m en rango} c_t(q,r),
```

con la diagonal `c_t(p,p)` solo si `p^2 | m`. En particular:

- **W(primo) = 0 exactamente.** Un primo no tiene divisores semiprimos: `qr | m` con `q,r >= 5` exige `qr = m`, imposible si `m` es primo. El peso no ve los primos.
- **W(P2) = 2 c_t(p1,p2)** para `m = p1 p2` distintos; **W(p^2) = c_t(p,p)**.
- **W(P3) = 2 [c_t(p1,p2)+c_t(p1,p3)+c_t(p2,p3)]** para `m = p1 p2 p3` distintos; **W(p^2 q) = c_t(p,p)+2 c_t(p,q)**; **W(p^3) = c_t(p,p)**.
- **W(P_k) = 2 sum_{i<j} c_t(p_i,p_j)** (mas diagonales si hay cuadrados).

El peso es **firmado** y **cero en primos**: no es un peso de criba. La suma sobre celdas `i<j` es incluso negativa (exacto, malla): `tri = (sum c - diag)/2 = -0.0956` a `h=96, t=0.05` (diag `0.2412`, dominada por la esquina `0.2375`); por conteo de celdas `i<j` hay `+4465 / -95` (h=96), pero las negativas de borde son grandes. No hay signo uniforme sobre P2.

## 3. Medicion a N=10^7 (exacto) y escala

`T = sum_p (log p) W(p) = sum_{q,r} c_t(q,r) A_{qr}(N)` (Prop. 59.1), descompuesto por `Omega(N-p)`. `G = S(N) N/log^2 N`, `S(N) = 1.7604` para `N = 10^7 = 2^7 5^7`, `G = 6.7763e4`. Familia A, `v=10`.

| h | t | T | T/G | E2 | E2/T | E2/G | E3/T | Ege4/T |
|---|---|---|---|---|---|---|---|---|
| 24 | 0.05 | 3601.7 | 0.0532 | 1587.4 | 0.441 | 0.0234 | -0.033 | 0.592 |
| 24 | 0.10 | 5129.0 | 0.0757 | 1535.4 | 0.299 | 0.0227 | 0.050 | 0.651 |
| 48 | 0.05 | 869.9 | 0.0128 | 426.5 | 0.490 | 0.0063 | 0.037 | 0.473 |

(`n2 = 10928` pares P2 con `N-m` primo; `n3 = 44679` triples. `T/(N*M1act) = 1.23` (h=24) y `1.44` (h=48): el error a `N=10^7` supera `M`, como advierte el doc 58 (`R_triv = 1.52`). Incertidumbre de discretizacion a `N=10^7`: factor ~4 entre h=24 y h=48 en los valores absolutos; los cocientes son estables en orden de magnitud.)

- **Escala a N=10^12** (M1 del doc 58, `5.54e-06`, h=96): `T = 5.54e6`, `G = 2.3058e9`, **`T/G = 2.40e-3`** (0.24%). El cociente `T/G` cae de `~10^-2` a `10^7` a `~10^-3` a `10^12` (factor 5.34; `N^{-1/5}` puro daria 10, el resto son logaritmos y h). El conteo pesado es `N^{4/5}/log^2 N`, no `N/log^2 N`.
- **La parte P2** es `E2/T ~ 0.3-0.5` del total, pero el total es `~10^-3-10^-2` de `G`: **`E2/G ~ 0.6-2.3%` a `N=10^7` y `~10^-3` a `N=10^12`**, frente a la constante `0.67` de Chen. La parte P3 es `~0` (`E3/T = -0.033..+0.050`, cancelacion); el resto es `Omega >= 4` (`0.47-0.65`).
- **Positividad sin contenido de conteo.** `T > 0` y `E2 > 0` son propiedades de la suma firmada, no cotas inferiores: los pesos individuales `2 c_t(p1,p2)` cambian de signo entre celdas (sec. 2), asi que `E2 > 0` no certifica ningun numero de P2. La desigualdad de Chen `#{1,2} >= W(N) - Omega/2` usa `W >= 0` (y `W = 1/2` en P3); con un peso firmado no hay desigualdad analoga.

## 4. Las tres obstrucciones

- **(O1) Objeto y lado equivocados.** `c_t` es un peso sobre modulos `(q,r)` en la descomposicion del factor local / error (docs 46, 58). La deteccion de P2 necesita un peso `W(n) >= 0` sobre *numeros* `n = N-p` (criba ponderada). El puente (Prop. 59.1) da `W(m) = sum_{qr|m} c_t(q,r)`, pero ese peso es firmado y cero en primos (sec. 2). No hay version no negativa: tomar `max(W,0)` rompe el analisis (el termino principal (46.3) controla la suma firmada, no la parte positiva).
- **(O2) Escala.** El conteo pesado es `N^{4/5}/log^2 N`: `T/G ~ 10^-3-10^-2`, decayendo como `N^{-1/5}`. La razon es el factor `1/phi(qr)` (doc 58): el soporte empieza en `alpha+beta = 1/5`. Ni siquiera el total, y menos su parte P2 (`E2/G ~ 10^-3` a `10^12`), esta a escala Goldbach.
- **(O3) Positividad de las cribas.** La cota superior de `Omega` (Liu §III) exige `f(a) >= 0`: `Omega <= M` usa `(sum lambda_d)^2 >= 1` y `f(a) >= 0`; con `f` firmado la desigualdad se rompe. La cota inferior `W(N)` exige `W(n) >= 0` en los numeros rugosos. `c_t` es firmado en ambos lados y no puede insertarse en ninguno.
- **Sustitucion en Omega (medido).** El analogo de `sum_a f(a)/(a log(N/a)) <= 0.49254/log N` (Liu Lema 4) con `g(a) = sum_{p1p2=a} c_t(p1,p2)` es `sum_{cuna} c_t/(qr log(N/qr)) = 8.29e-10` (h=96, N=10^12), frente a `0.49254/log N = 0.017825`: **ratio 4.65e-08**. El Omega sustituido es `~10^-8` del de Chen. (La suma sin ponderar sobre la cuna es `0.0204 = 41%` de `t`: la cuna si captura masa de `sum c`, pero el peso `1/(qr)` la aniquila.)
- **Combinacion (vacua).** `W_hibrido(n) = W_Chen(n) + eps*W(n)`: el termino anadido aporta `eps * N^{4/5}/log^2 N`, es decir una mejora relativa `<= eps * N^{-1/5}/S(N) -> 0`. No cambia la constante `0.67` ni el error (que es `O(N/log^A N)` por BV). Ni ayuda ni estorba: es despreciable.

## 5. Veredicto

**La integracion NO detecta P2 con nivel 1/2.** El balance sobrevive positivo en signo (`T > 0`, `E2 > 0`), pero no a escala Goldbach (`T/G ~ 10^-3`, `E2/G ~ 10^-3` a `N=10^12`, decayendo), el peso inducido es firmado y cero en primos, y ninguna parte de Chen usa nivel > 1/2 que el termino principal `sum c = t` pudiera reemplazar (Chen usa BV 1/2 + switching; el switching es precisamente lo que compensa la falta de nivel > 1/2). La linea de los docs 54-59 queda cerrada como ruta a P2/Goldbach: construyo un peso nuevo con trade-off optimo (Lema 56.1), nivel 1/2 medido (doc 58) y dispersion auditada (doc 59), pero el peso no se acopla a Chen. Goldbach sigue abierta.

## 6. Que sigue

1. **Bajar el soporte de los factores a `N^eps`** (doc 59, #2): arreglaria la well-factorabilidad (Prop. 58.1) y la escala (`M1 ~ N^{-2eps}/log^2 N`), pero NO el cero-en-primos ni el signo: `W(primo) = 0` es estructural (el peso solo ve numeros con `>= 2` factores en rango). Exige redisenar la construccion (otro potencial, otra malla) Y anadir un componente P1: un peso sobre numeros, no sobre pares.
2. **Peso sobre numeros `W(n) >= 0` con cancelacion marginal:** re-expresar (46.4) del lado de los numeros (criba ponderada genuina). Es un cambio de marco, como anticipan el Lema 56.1 y la Prop. 58.1.
3. **Frentes del bucle:** C10/peso movil W(d), uniones C_k+C_{k+1}, variante Lema 4 Pintz-Althofer, Motohashi (doc 57 sec. 4), `b_{P,l}` (doc 53). Dispersion uniforme en `a` (doc 59 #1) sigue abierta como problema independiente.

## Verificacion reproducible

```text
python numerico/60_integracion_peso_relajado_chen_p2.py
```

PASS en el PC: identidades exactas de malla (sum c = t; tri = (sum-diag)/2); enumeracion N=10^7 (T>0; E2+E3+Ege4 = T por construccion; n2=10928, n3=44679); escala a N=10^12 (M1 doc 58); cuna con peso 1/(qr) (ratio 4.65e-08). NOTA: el predictor uniforme E2esp (probabilidad 1/log N) es INVALIDO: predice signo contrario porque los pares del extremo (p1p2~N) tienen N-p1p2 pequeno y casi nunca primo; se reporta solo como diagnostico de ese fallo del modelo.

## Alcance honesto

- Lo que este doc SI dice: Chen usa nivel 1/2 y nada mas (Liu, leido); el peso inducido es firmado y cero en primos; el balance es positivo pero a escala `N^{4/5}`; la sustitucion en Omega es `10^-8` de Chen; la combinacion es vacua.
- Lo que NO dice: no hay cota inferior de capas, no hay deteccion de P2, no hay demostracion de Goldbach.

## Memoria compartida

Propuesta pendiente: el MCP canonico `lumina_shared_memory` / `lumina_propose` no esta expuesto. No se afirma persistencia remota.
