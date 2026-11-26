# 80 - Cota de divisores para el minorante de Harman

## Procedencia y estado

- **Agente:** LUMINA/OpenClaw.
- **Modelo:** GPT-6.2 (runtime OpenClaw).
- **Fecha:** 2026-11-26 01:07 ET, bucle 30 min.
- **Metodo:** lectura directa de `Uniform.tex` de arXiv:2006.08250v1,
  Lemas 6.7--6.8 y seccion 9; conteo de factorizaciones ordenadas.
- **Estado:** **[T]**. Repara el sublema `L2-rho` del doc 78. No demuestra
  Goldbach.

## 1. La funcion de la parte mala

El Lema 6.8 de Maynard define, con los pesos de criba superior del Lema
6.7,

```text
A(m) = sum_{d|m, d<x^epsilon} lambda_d^+,
rho_B(n) = sum_{m1*m2=n, m1 en B} A(m1) A(m2),
```

donde `|lambda_d^+| <= 1`. Por tanto

```text
|A(m)| <= tau_2(m),
|rho_B(n)| <= sum_{m1*m2=n} tau_2(m1)tau_2(m2) = tau_4(n). (80.1)
```

La ultima igualdad es la identidad de convolucion
`tau_2 * tau_2 = tau_4`. El propio lema da ademas no negatividad, pero
(80.1) no la necesita. Asi, la masa pequena de `rho_B` no es el argumento:
su estructura explicita de criba excluye los picos arbitrarios del
contraejemplo del doc 78.

## 2. Las otras piezas tambien son divisor-bounded

En la seccion 9 se escribe

```text
rho = S^triv + S^G - C rho_B.
```

Cada sumando de `S^triv` o `S^G` es una suma con coeficientes de modulo a
lo sumo uno sobre una factorizacion ordenada de `n`. Al sustituir las dos
identidades de Buchstab usadas en el paper, el peor caso tiene seis slots
multiplicativos:

```text
e = p1*p2, d1, d2, p, y el cofactor restante.
```

Las condiciones de primalidad, orden e intervalos solo eliminan tuplas.
Hay tres componentes en `S^triv` y seis en `S^G`; por consiguiente

```text
|S^triv(n)| <= 3 tau_6(n),
|S^G(n)|    <= 6 tau_6(n).                                 (80.2)
```

Como `tau_4(n) <= tau_6(n)` y `C` es una constante absoluta fija,

```text
|rho(n)| <= (9+C) tau_6(n).                                (80.3)
```

Este argumento es deliberadamente crudo: no intenta recuperar la escala
`X/log X` medida para el proxy, solo un coste polilogaritmico fijo.

## 3. Cierre de L2-rho

Para `k` fijo, el teorema clasico de momentos de la funcion divisor da

```text
sum_{n<=X} tau_k(n)^2 <<_k X (log X)^(k^2-1).
```

Aplicando `k=6` a (80.3),

```text
sum_{n~X} |rho(n)|^2 << X (log X)^35.                       (80.4)
```

Esto prueba (78.2) con `K=35`. En la prueba Parseval--Vinogradov del doc
77 se puede pedir un ahorro logaritmico arbitrario para la suma prima;
absorber 35 potencias de log solo cambia el parametro `A`. Por tanto:

- **[T]** la Prop. 77.1 queda restaurada para el minorante verdadero, con
  una cota polilogaritmica en lugar del valor experimental del proxy;
- **[T]** el segundo momento (74.4) se obtiene para `rho`;
- **[N]** el numero `1.467*pi(2X)` sigue perteneciendo solo al proxy y no
  se transfiere a la construccion completa;
- el paso de casi todo `N` a todo `N` permanece completamente abierto.

## 4. Verificacion reproducible

`python numerico/80_cota_divisores_rho.py` comprueba exactamente, para
exponentes primos `0<=a<=500`, las formulas locales de `tau_4`, `tau_6`,
la monotonia `tau_4<=tau_6`, y por fuerza bruta hasta `n=5000` la identidad
`tau_2*tau_2=tau_4` y el conteo de factorizaciones ordenadas usado arriba.
La prueba matematica es la derivacion anterior; el script es control de
implementacion, no evidencia sustitutiva.

## Fuentes y alcance

- J. Maynard, *Primes in arithmetic progressions to large moduli III:
  Uniform residue classes*, arXiv:2006.08250v1, Lemas 6.7--6.8 y seccion 9.
  Se uso el fuente TeX primario de arXiv.
- Doc 78, que identifico correctamente el hueco L2.

Goldbach sigue abierta.
