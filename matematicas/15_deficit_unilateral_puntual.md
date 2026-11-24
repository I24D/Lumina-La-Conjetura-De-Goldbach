# 15 — Déficit unilateral puntual: calibración de ρ(N)

## 0. Procedencia y estado

- Agente: LUMINA/OpenClaw (sesión telegram).
- Fecha: 2026-11-24, America/New_York.
- Estado: **[P]** para las definiciones e identificaciones elementales de las
  secciones 1–2 (reformulaciones de 04/11); **[N]** para la medición de la
  sección 4. No es una demostración de Goldbach.
- Umbral de pre-registro fijado **antes** de leer los resultados (sección 3).

## 1. El objeto puntual

De `04` (Prop. 4.1, verificada en Lean) y su versión primitiva (4.15):

\[
R(N)=S_N(D)-T_N(D),\qquad
R(N)=N\mathfrak S(N)-T_N^{\rm prim}(D)+O_A(N/(\log N)^A),
\]
con `D=floor(N^{1/3})`, `T^prim` la cola restringida a `(c,N)=1` y
`|T_bad|<=K_N` (Prop. 4.2, también en Lean). Definimos, para `N` par:

\[
\rho(N)=\frac{T_N^{\rm prim}(D)}{N\mathfrak S(N)},\qquad
r(N)=\frac{R(N)}{N\mathfrak S(N)},\qquad
\mathrm{gap}(N)=1-r(N)-\rho(N).
\]

Elementalmente, `gap(N)=(N𝔖(N)-S_N(D))/(N𝔖(N))+T_bad/(N𝔖(N))`: absorbe el
error del término principal y la parte mala. Ojo: `T_bad` solo está acotada
por `K_N` (Prop. 4.2), **no** es pequeña punto a punto a priori; la
sección 4 mide el `gap` real por `N`.

La lección de `02` (Li–Liu) y la hipótesis `H_ε` de `04` §5 dicen que
Goldbach para pares grandes solo necesita la cota **de un solo lado**:

\[
\rho(N)\le 1-\varepsilon \quad\text{eventualmente.}
\]

No se necesita cota inferior de la cola. El objetivo de `06` R-2
("mejorar la constante 1 a 1−ε") es exactamente esta desigualdad.

## 2. Identificación del supertrazo Z/2

`11` define las capas por paridad de `omega(d)`:

\[
A_j(N,D)=\sum_{\substack{cd<N,\ d>D,\ (c,N)=1\\ d{\rm\ libre\ de\ cuadrados}\\
\omega(d)=j}}(\ldots)\mu(d)\log d\,\Lambda(N-cd),
\]
y su identidad (11.1):

\[
T_N^{\rm prim}(D)=\sum_{j\ge1}(-1)^jA_j(N,D).
\]

Es decir: **el supertrazo de la graduación Z/2 de `13` §6 es, por
definición, la cola primitiva de `04`.** Las masas par/impar de `11`
(17.934M / 17.870M en `N=10^6`) difieren exactamente en `T^prim`
(64 249, igual al cruce de `06`). Por tanto la desigualdad de pendientes
buscada en `13` §6 y la `H_ε` de `04` son la misma condición, ahora
escrita de forma puntual y normalizada:

\[
\boxed{\ \rho(N)\le 1-\varepsilon\ }
\]

La pregunta empírica que este documento calibra: **¿qué valor tiene
`sup_{N par} rho(N)` y cómo deriva con la escala?** La medición es
puntual en `N` (no un promedio, como el Teorema 7.1 de `07`), tal como
exige `14` §5.

## 3. Diseño de la medición y umbral pre-registrado

`numerico/deficit_unilateral.py` reutiliza `exact_row` de
`cola_mobius_scan.py` (validado por el cruce independiente de `06`:
`T^prim/N = 0.0642495` en `N=10^6`):

- `D=floor(N^(1/3))` como en `04`; cola primitiva literal
  (`d>D`, `mu(d)!=0`, `gcd(m/d,N)=1` con `m=N-n`);
- `𝔖(N)` con `C2` truncada en primos <= 200 000 (suficiente para
  `N<=10^6`);
- lotes: **todos** los `N` pares en `[4,2000]`; muestreo con semilla
  20261124: 300 pares en `(2000,20000]`, 150 en `(20000,200000]`,
  40 en `(200000,10^6]`. Ejecución 2 (densificación, mismos umbrales
  pre-registrados): 600/400/300 en las mismas bandas.
- por `N`: `rho`, `r`, `gap`; por banda y global: supremo, `N` del
  supremo, top-10, mediana, mínimo (lado dos caras de `06` R-2),
  `max gap` y `min r`.

**Umbral pre-registrado** (fijado antes de ejecutar):

- si `sup rho >= 0.9` en alguna banda con deriva creciente entre
  escalas: la ruta unilateral está empíricamente moribunda y se
  registra como descarte;
- si `sup rho <= 0.5` estable (sin deriva mayor que `+0.05` por
  década): `H_ε` es plausible con margen amplio y el problema analítico
  es la *uniformidad* en `N` (los `N` raros), no la constante;
- caso intermedio: registrar y refinar.

## 4. Resultados [N]

**Ejecución 1** (lote pre-registrado; `numerico/resultados/20261124_224341_
deficit_unilateral_1000000.json`; 1489 `N` pares; 21 s; rev `9ae94d2`):

| banda | n | sup ρ | N* | mediana ρ | min ρ |
|---|---|---|---|---|---|
| 4–2000 (exhaustivo) | 999 | 0.2755 | 1556 | −0.400 | −0.946 (N=1260) |
| 2000–20000 | 300 | 0.3601 | 9224 | −0.253 | −0.792 (N=3510) |
| 20000–200000 | 150 | 0.3685 | 30658 | −0.064 | −0.712 (N=22950) |
| 200000–10^6 | 40 | 0.2464 | 599948 | +0.034 | −0.561 (N=443190) |

**Veredicto contra el umbral pre-registrado: `sup rho = 0.3685 <= 0.5`,
sin deriva ascendente entre escalas.** Para todos los `N` pares probados,
`rho(N) <= 0.37`; el margen unilateral disponible es `eps >= 0.63` en el
rango medido. La mediana migra −0.400 → −0.253 → −0.064 → +0.034: la cola
relativa se encoge al crecer `N`. El lado dos caras de `06` R-2 no se
aproxima: el mínimo observado es −0.95 (en `N=1260`, asintótica inválida
ahí) y −0.71 ya en banda adulta; la cota inferior `−(3+o(1))N𝔖` queda
lejos.

**Estructura de los peores `N`** (factorización verificada):
`30658=2·15329`, `9224=2³·1153`, `29632=2⁶·463`, `7702=2·3851`,
`119398=2·59699`, `7442=2·61²`, `25076=2²·6269`, `109958=2·54979`,
`9028=2²·37·61`, `14876=2²·3719`, `1556=2²·389`, `599948=2²·127·1181`.
Todos tienen **un solo primo impar grande** (o `61²`, o dos medianos): son
exactamente los `N` con serie singular mínima `𝔖(N)≈2C₂`. El supremo vive
donde el término principal es más pequeño — los mismos `N` "difíciles" de
Goldbach.

**Cadena (4.1)+(4.15) validada por `N`:** `gap = 1−r−ρ` es pequeño en `N`
grandes típicos (|gap| ≲ 0.1); los gaps grandes aparecen co-localizados con
`ρ` muy negativo (p. ej. `N=443190`: `ρ=−0.561`, `r≈1.000`, `gap≈0.561`):
la parte mala `T_bad` solo está acotada por `K_N` — resultado honesto
sobre dónde la cadena asintótica sigue floja punto a punto. En banda alta,
`min r = 0.976` (`N=348224`).

**Ejecución 2 (densificación, mismos umbrales pre-registrados):**
`numerico/resultados/20261124_224653_deficit_unilateral_dense_1000000.json`
(2299 `N` pares; 117 s):

| banda | n | sup ρ | N* | mediana ρ | min ρ |
|---|---|---|---|---|---|
| 4–2000 (exhaustivo) | 999 | 0.2755 | 1556 | −0.400 | −0.946 (N=1260) |
| 2000–20000 | 600 | 0.3853 | 2284 | −0.168 | −0.830 (N=11220) |
| 20000–200000 | 400 | 0.4201 | 30916 | −0.007 | −0.792 (N=39900) |
| 200000–10^6 | 300 | 0.2535 | 606722 | −0.013 | −0.649 (N=289170) |

La densificación elevó el supremo global de 0.3685 a **0.4201**
(`N=30916`) pero **sigue <= 0.5 y sin deriva ascendente**: la banda alta
(300 muestras en `(2·10^5, 10^6]`) se estabiliza en 0.25 con mediana
≈ −0.01. El peor grupo se agrupa en `N≈30–33·10^3`
(30916=2²·59·131, 31718=2·15859, 30482=2·15241, 33172=2²·8293),
todos con serie singular casi mínima. **Veredicto robusto: `H_eps`
plausible con margen amplio; el margen unilateral medido disponible es
`eps >= 0.58` en el rango probado.**

## 5. Consecuencia para el programa

- El margen unilateral es amplio y estable (`ρ ≤ 0.37` medido); lo que
  exige la prueba es *cualquier* `ε>0` **uniforme** en `N`. El problema de
  uniformidad tiene ahora forma concreta: los peores `ρ` se concentran en
  las familias de `𝔖` mínima (`N=2^a·p`, `N=2·p²`, `N=2^a`).
- **Siguiente experimento pre-registrable:** barrido exhaustivo de las
  familias de `𝔖` mínima (`N=2^a ≤ 10^6`, `N=2p`, `4p`, `2p² ≤ 10^5`) y
  contraste de `sup ρ` familiar contra el global. Si el supremo global lo
  determinan estas familias delgadas, la uniformidad pendiente se reduce
  a dos sucesiones con estructura multiplicitiva explícita — un blanco
  mucho más pequeño para un argumento maximal/paridad (familia `08` §4)
  que el "todos los N".

## 5-bis. Resultado del barrido de familias [N]

Ejecutado el mismo día (`numerico/deficit_familias.py`, 15 458 filas
exactas, 396 s; salida
`numerico/resultados/20261124_225515_deficit_familias.json`).
Pregunta pre-registrada en §5: ¿`sup ρ` familiar ≥ 0.4201 (sup global
denso)? **Veredicto: SÍ — las familias delgadas capturan y superan el
supremo global.**

| familia | n | sup ρ | @ N | mediana ρ | captura |
|---|---|---|---|---|---|
| `2^a` (a≤19) | 18 | 0.3990 | 8192 | +0.044 | no |
| `2p` (N≤10^5) | 5132 | **0.4816** | 2218 | **+0.231** | sí |
| `4p` (N≤10^5) | 2761 | **0.4798** | 34324 | **+0.230** | sí |
| `2p²` (N≤10^5) | 47 | 0.3945 | 32258 | +0.186 | no |
| ventana (25000,40000] | 7500 | **0.4841** | 32728 | ≈0.001 | sí |

La ventana (25000,40000] es exhaustiva y **post-observacional**
(añadida tras ver que el sup global denso estaba en 30916=2²·59·131,
que no pertenece a ninguna familia estricta); se reporta separada
por higiene del pre-registro.

Observaciones:

1. **Nuevo supremo medido de `ρ`: 0.4841 en `N=32728=2³·4091`** — otra
   vez `2^a·(primo grande)`. El sup familiar no crece con `N` dentro
de la familia: en `2p` ocurre en `N=2218`, en `4p` en `N=34324`;
en ningún brazo el máximo está en el extremo superior del rango.
2. **Las medianas de `2p` y `4p` son +0.23** — es decir, en las
   familias de `𝔖` casi mínima la cola primitiva es sistemáticamente
   *más grande* que la media (frente a mediana global ≈ −0.013). Esto
   es exactamente lo que predice la heurística: `𝔖(N)` mínimo →
   denominador pequeño → `ρ=T_prim/(N𝔖)` grande. La cola primitiva
   es un *déficit crónico* en estas familias, no un fenómeno puntual.
3. El mínimo absoluto del brazo ventana es `ρ=−0.878` en
   `N=30030=2·3·5·7·11·13` (primorial): los primoriales tienen cola
   primitiva fuertemente negativa. Extremos positivos y negativos de
   `ρ` viven en sucesiones de estructura multiplicitiva explícita.
4. Consecuencia para `H_ε`: el `ε` medido se ajusta a
   `ε ≥ 1−0.4841 ≈ 0.52`, aún amplio, pero el sup es sensible al
   muestreo (0.37 global → 0.48 en familias delgadas). La uniformidad
   pendiente queda acotada al problema: ¿puede `sup ρ` de `2p`/`4p`
   crecer con `N` más allá de ~0.5? Siguiente experimento
   pre-registrable: `2p` y `4p` densos hasta `10^6` (no solo `10^5`).
- `07` da control promedio (Parseval sobre `N`); el conjunto excepcional
  observado tiene justo la forma de familias delgadas de `𝔖` mínima.
- No se presupone positividad de `R(N)`: `r(N)` se registra como dato;
  la desigualdad objetivo solo acota la cola por arriba.

## 6. Verificación reproducible

- Comando: `python numerico/deficit_unilateral.py --limit 1000000`
  (ejecución 1) y `--dense` (ejecución 2); semilla 20261124.
- Barrido de familias: `python numerico/deficit_familias.py` (sin RNG;
  familias deterministas + ventana post-observacional etiquetada).
- Salida: `numerico/resultados/20261124_*_deficit_unilateral_*.json`.
- Dependencias: `cola_mobius_scan.py` (mismo directorio).
