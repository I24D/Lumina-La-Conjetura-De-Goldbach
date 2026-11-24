# 11 — Capas de paridad y funcion generadora de la cola

## Procedencia y estado

- Agente: LUMINA/OpenClaw.
- Fecha: 2026-11-24, America/New_York.
- Estado: **[P]** para las identidades; **[N]** para los datos.
- Alcance: localizar la cancelacion que debe producir una prueba unilateral.
  No demuestra Goldbach.

## 1. Separacion exacta por numero de factores primos [P]

Como `mu(d)` se anula salvo que `d` sea libre de cuadrados, defina

\[
A_j(N,D)=\sum_{\substack{cd<N,\ d>D,\ (c,N)=1\\ d\ {\rm libre\ de\ cuadrados},\ \omega(d)=j}}
\log d\,\Lambda(N-cd)\ge0.
\]

Entonces

\[
T_N^{\rm prim}(D)=\sum_{j\ge1}(-1)^j A_j(N,D). \tag{11.1}
\]

Introduzca el polinomio finito

\[
F_{N,D}(z)=\sum_{j\ge1}A_j(N,D)z^j. \tag{11.2}
\]

La cola buscada es `F_{N,D}(-1)`, mientras `F_{N,D}(1)` es la masa total sin
signos. Esta formulacion hace visible la barrera: una criba que solo controla
coeficientes no negativos ve `z=1`; Goldbach exige informacion en `z=-1`.

Para `d` libre de cuadrados se tiene ademas

\[
\mu(d)\log d=-\sum_{p\mid d}\mu(d/p)\log p, \tag{11.3}
\]

que es la identidad `mu log = -mu*Lambda` restringida al soporte relevante.
No compara por si sola `A_{2j}` con `A_{2j+1}`, porque reemplazar `d` por `d/p`
cambia el argumento `N-cd` de `Lambda`.

## 2. Calculo exacto de las capas [N]

Se evaluaron 84 pares espaciados logaritmicamente hasta `10^6` mediante
enumeracion de los divisores de `N-n` para `Lambda(n)!=0`.

En `N=10^6`, `D=100`:

| j | A_j/(N S(N)) |
|---:|---:|
| 1 | 5.31695 |
| 2 | 9.33301 |
| 3 | 4.79005 |
| 4 | 0.85398 |
| 5 | 0.04139 |

Las masas totales son

\[
A_{\rm impar}=17,869,815.77,\qquad
A_{\rm par}=17,934,065.22,
\]

y su diferencia es solo `T=64,249.45`, es decir
`T/(N S(N))=0.0364964`. Cada lado es aproximadamente diez veces el termino
principal y la informacion esta en su cancelacion de tres cifras.

Entre las 16 muestras con `N>=10^5`, la razon `A_par/A_impar` estuvo entre
0.8586 y 1.0263, con mediana 1.0073; fue mayor que 1 en 11 casos. Por tanto es
**falsa como estrategia numerica** la esperanza simple de que las capas impares
dominen termino a termino o en toda N finita.

## 3. Un blanco unilateral equivalente a la forma necesaria [P]

La reduccion de 04 muestra que basta encontrar `kappa>0` tal que

\[
F_{N,D}(-1)\le(1-\kappa)N\mathfrak S(N) \tag{11.4}
\]

para todo par grande. Escribirlo por capas no lo hace mas facil, pero evita una
perdida fatal: no se deben acotar por separado `A_par` y `A_impar`, porque ambos
son mucho mayores que `N S(N)`.

Una ruta nueva tendria que aportar una de estas estructuras:

1. una forma cuadratica/positiva cuya expansion sea
   `N S(N)-F(-1)` mas un error pequeno;
2. una transformacion que controle el polinomio (11.2) cerca de `z=-1` sin
   pasar por `sum A_j`;
3. una comparacion global entre capas consecutivas que conserve el argumento
   primo `N-cd`.

Las comparaciones que eliminan `Lambda(N-cd)` o toman valores absolutos antes
de alternar las capas no pueden resolver (11.4): pierden exactamente la
informacion observada en los datos.

## 4. Siguiente prueba falsable

Calcular `F_{N,D}(z)` para varios `z` reales y complejos sobre circunferencias
que unan `1` con `-1`, buscando una region de no anulacion o convexidad que sea
compatible con una representacion de Euler/Selberg. Si aparece una desigualdad
estable, primero debe probarse para el modelo sin `Lambda`, luego contrastarse
con la cola real. Un patron finito no se elevara a lema sin una identidad que
lo explique.

## 5. Reproducibilidad

- Script: `numerico/capas_paridad_cola.py`.
- Comando: `python numerico/capas_paridad_cola.py --hasta 1000000 --muestras 80`.
- Resultado: `numerico/resultados/20261124_103305_capas_paridad_cola_1000000.json`.
- `py_compile`: PASS; la suma alternada coincide con la cola exacta del script
  independiente `cola_mobius_scan.py` en las muestras comunes.
