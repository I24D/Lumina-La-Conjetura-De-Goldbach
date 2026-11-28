# Doc 127 — Errata y reparaciones del doc 125 (auditoría adversarial del doc 126)

**Fecha:** 2026-11-28
**Agente:** LUMINA Muse
**Estado:** errata verificada; no es una demostración de Goldbach.
**Base:** doc 125 (125a_prime_window.md, 125b_erem_pos.md y sus
scripts) y doc 126 (auditoría adversarial de LUMINA Chat,
126_auditoria_lumina.md). Cada corrección se verificó contra el
texto del 126 y contra los docs/scripts 125 antes de aceptarla.
**Reproducibilidad:** script `127_reparaciones.py`, salida
`salida_doc127_local.txt`: chequeos R1–R5.

Etiquetas: [T] = analítico; [M] = medición/enumeración;
[P] = pendiente.

---

## 0. Estado: qué sobrevive del doc 125 tras la auditoría

Según la propia auditora (doc 126, §0 y §5), y confirmado por la
verificación propia de este doc (§§1–5), el núcleo del doc 125
sobrevive íntegro:

- **(PRIME-WINDOW) refutada y (OUT-MASS) absoluta caída [T]**
  en el régimen u_L → 2: el término principal de la doble suma
  es F(2) = log(3/2)·log(4/3) = 0,1166450426… > 0 mientras
  Q_c(z) ~ c_Q/log z → 0; el cociente diverge como
  (F(2)/c_Q)·log z con F(2)/c_Q ≈ 0,0678.
- **Lemas de paridad de 125b y umbrales 343 y 1331 [T]:**
  E_Rem = 0 para L/r ≤ 343; E_Rem ≤ 0 para L/r ≤ 1331;
  (EREM-SIGN) abierto solo para L/r > 1331.
- **Cierre analítico del estrato 5 vía (POS) [T]:** c ∈ {1,5,7}.
- **La anulación uniforme E_Rem = 0 es falsa [T-instancia]:**
  contraejemplo exacto reforzado por la auditora y verificado
  con las funciones de la serie en la §5 de este doc.
- La ruta con signo (NET-LAYER) + (P-MASS) sigue sin refutar
  y sin demostrar: es el muro restante de la banda alta.

Las cuatro correcciones del doc 126 son reales y se aceptan
todas; ninguna cambia un veredicto del doc 125. Una precisión
de atribución: la corrección 3 (signo de ε_Π) recae sobre el
doc **125a** (§3), no sobre el 125b — la §3 lo documenta con
evidencia de grep. Detalle por corrección en §§1–4.

---

## 1. Corrección 1 — El dominio de la fórmula PW es (3/2, 9/4), no (3/2, 3)

**Señalamiento (126, §1.2).** El 125a imprime la asintótica

$$PW(L) = \log\frac{3}{u_L}\cdot
\log\frac{u_L/3}{u_L - 3/2} + O\!\left(\tfrac1{\log z}\right)$$

«para u_L ∈ (3/2, 3) fijo» (125a, líneas 127, 145 y 151). Pero
la suma exterior corre sobre los primos permitidos en
(L z^{−3/2}, L^{1/3}], y ese intervalo solo es no vacío si

$$L z^{-3/2} < L^{1/3} \iff \tfrac23 \log L < \tfrac32 \log z
\iff u_L < \tfrac94.$$

Para u_L = 9/4 los extremos coinciden (intervalo vacío); para
9/4 < u_L < 3 el cociente (u_L/3)/(u_L − 3/2) es < 1 y el
segundo logaritmo sería negativo, aunque PW(L) ≥ 0 por
definición: la fórmula impresa es imposible en (9/4, 3).

**Verificación propia [T + M, chequeo R1].** El álgebra de
extremos se rehízo: con L = z^u, L z^{−3/2} < L^{1/3} ⟺
u − 3/2 < u/3 ⟺ (2/3)u < 3/2 ⟺ u < 9/4. El script R1 barre
u ∈ {1,6; 2; 2,2; 9/4; 2,5; 2,9} y comprueba en cada caso el
signo de log((u/3)/(u − 3/2)) y la vacuidad del intervalo:
positivo y no vacío exactamente en (3/2, 9/4); cero en 9/4;
negativo (fórmula imposible) en (9/4, 3). Además
F(2) = log(3/2)·log(4/3) = 0,1166450426…, el valor que la
auditora cita.

**Corrección aplicada.** El enunciado correcto es: para
u_L → u fijo con **3/2 < u < 9/4**,

$$PW(L) = \log\frac{3}{u}\,
\log\!\left(\frac{u/3}{u-3/2}\right) + o(1) \qquad [T].$$

Quedan supersedidos los pasajes del 125a que extienden el
análisis a u_L ∈ (3/2, 3), incluido el comentario «si u_L → 3
el factor log(3/u_L) → 0 y el análisis cambia» (125a, §3):
para u_L ≥ 9/4 no hay capa₁ con las dos condiciones y no hay
fórmula que extender.

**Efecto en los veredictos: ninguno.** El régimen de la
refutación es u_L ≈ 2 (los marcos de la serie viven ahí y la
sucesión admisible de la auditora da u_L = 2 + O(1/log z)),
estrictamente dentro de (3/2, 9/4). La refutación de
(PRIME-WINDOW) y la caída de (OUT-MASS) absoluta no se tocan.

---

## 2. Corrección 2 — La constante de Mertens M′ está corrida por 2/3

**Señalamiento (126, §1.3).** El 125a define
S(x) = Σ_{5≤p<x} 1/(p−1) y escribe (línea 86)

$$M' = M - \tfrac32 + \sum_{q\ge5} \frac1{q(q-1)},$$

justificándolo con «quitar los primos 2 y 3 resta 1 + 1/2».
Pero la descomposición usada es 1/(p−1) = 1/p + 1/(p(p−1)),
y quitar p = 2, 3 de Σ 1/p resta 1/2 + 1/3 = **5/6**, no
1 + 1/2 = 3/2. (El 1 + 1/2 impreso son los valores de
1/(p−1) en p = 2, 3: se restaron las contribuciones de la
suma equivocada.) La constante correcta es

$$\boxed{\,M'_{corr} = M - \tfrac56 +
\sum_{p\ge5}\frac1{p(p-1)}\,}$$

y la impresa es menor por exactamente 3/2 − 5/6 = **2/3**.
Con la constante impresa, E(x) := S(x) − log log x − M′ no
tiende a 0: tiende a 2/3.

**Verificación propia [T + M, chequeo R2].** (i) Aritmética
racional exacta (fracciones): 1/2 + 1/3 = 5/6 y
3/2 − 5/6 = 2/3. (ii) Identidad finita exacta a precisión de
máquina para x = 10^5: Σ_{5≤p<x} 1/(p−1) =
Σ_{5≤p<x} 1/p + Σ_{5≤p<x} 1/(p(p−1)), y
Σ_{5≤p<x} 1/p = (Σ_{p<x} 1/p) − 5/6. (iii) Numérico:
M′_corr = −0,465346… (con B_1 = 0,26149721… y la serie
Σ_{p≥5} 1/(p(p−1)) = 0,10649… sumada hasta 10^7 con cola
< 10⁻⁷; la constante impresa da −1,132013…), y a x = 10^7:
E_corr(x) = S(x) − log log x − M′_corr = +0,0000 (a 4
decimales), mientras E_impresa(x) = +0,6667: con la
constante impresa E(x) no va a 0, va a 2/3, como dijo la
auditora. Coherente con la medición del propio 125a
(§4, K1): E(z) = S(z) − log log z ∈ [−0,452; −0,446] para
z ≤ 71, en camino hacia M′_corr ≈ −0,465, no hacia la
constante impresa ≈ −1,132.

**Corrección aplicada.** Sustituir M′ por M′_corr en la
definición de la §2 del 125a y en la afirmación «E(x) → 0»:
con M′_corr, E(x) = O(1/log x) y la derivación queda limpia.

**Efecto en los veredictos: ninguno.** En la identidad de
ventana M_near⁽¹⁾ = S(z) − S(L^{1/3}) la constante aditiva
entra como E(z) − E(L^{1/3}) y se cancela con cualquiera de
las dos constantes; el término principal log(3/u_L) y la
refutación no dependen de M′. El error era de normalización
de la §2, no de la §3.

---

## 3. Corrección 3 — El signo de ε_Π es −, y está en el 125a (no en el 125b)

**Señalamiento (126, §1.4).** En la suma exterior, el 125a
define ε_Π ≥ 0 como la masa de los primos del intervalo
completo excluidos por Π(p₀) ≥ 1/2 o degeneración, pero
escribe (125a, §3, líneas 114–117)

$$\sum_{p_0\in\mathrm{capa}_1}\frac1{p_0-1}
= S\!\left(L^{1/3}\right) - S\!\left(L z^{-3/2}\right)
+ \varepsilon_\Pi.$$

Si ε_Π es masa excluida no negativa, el signo correcto es
**menos**: suma filtrada = suma del intervalo − ε_Π.

**Verificación propia [T + M, chequeo R3].** (i) Localización
exacta: ε_Π aparece únicamente en el 125a (md: líneas
115–122; la definición «ε_Π ≥ 0 es la masa de los primos
excluidos» está en la línea 119). **No aparece en el 125b**
ni en su script (grep sin resultados en 125b_erem_pos.md y
125b_erem_pos.py). La atribución de esta corrección al 125b
que circuló en el encargo de este doc no se sostiene: la
corrección es del doc 126 §1.4 y recae sobre el 125a §3.
(ii) El signo es contabilidad pura: para cualquier conjunto
finito, Σ_filtrada = Σ_completa − Σ_excluida; el script R3 lo
comprueba sobre el intervalo real (L z^{−3/2}, L^{1/3}] de un
nivel de ventana con un filtro artificial de masa ε > 0: la
fórmula impresa (con +) falla y la corregida (con −) es
exacta.

**Corrección aplicada.** En el 125a §3 debe leerse

$$\sum_{p_0\in\mathrm{capa}_1}\frac1{p_0-1}
= S\!\left(L^{1/3}\right) - S\!\left(L z^{-3/2}\right)
- \varepsilon_\Pi
= \log\frac{u_L/3}{u_L-3/2} + \Delta E' - \varepsilon_\Pi.$$

**Efecto en los veredictos: ninguno, y ninguna afirmación
del 125b depende de ε_Π** (no ocurre en ese doc). En el 125a
el propio texto acota el daño: ε_Π = 0 en toda la ventana
medida y ε_Π = 0 para z grande (Π(p₀) ≍ c/log z < 1/2 y
L/p₀ > 7³ en la capa), así que el límite en u → 2 y la
refutación quedan intactos. Era un error algebraico de
escritura, no de la cadena.

---

## 4. Corrección 4 — El chequeo K3 del 125b era vacuo; reparación y re-ejecución

**Señalamiento (126, §3).** En `125b_erem_pos.py`, función
`check_K3`, subfunción `sqfree_cs`:

    def rec(start, cur):
        ...
        for p in pr:
            if p < start or cur * p >= B:
                break

La lista `pr = primes_upto(100000)` empieza por 2, y la
llamada inicial es `rec(5, 1)`: al primer paso p = 2 < 5
dispara el `break` y la recursión muere. `sqfree_cs(125)`
devuelve la lista vacía, el conjunto de candidatos que
sobreviven queda vacío **por construcción**, y esa subprueba
de K3 pasaba de forma vacua. El PASS global del 125b no
validaba ese barrido.

**Verificación propia [M, chequeo R4].** Se reprodujo la
subfunción original verbatim: `sqfree_cs(125) = []`,
confirmado. El defecto es exactamente el diagnosticado.

**Reparación aplicada.** La que sugiere la auditora:
separar las dos condiciones — `if p < start: continue`,
conservando `break` solo para `cur * p >= B`. Re-ejecución
real del barrido reparado:

- Candidatos generados: **37** squarefree c < 125 con todos
  sus factores primos ≥ 5 — exactamente la cifra del doc
  126 — y un generador por fuerza bruta independiente
  (factorización por división tentativa, sin recursión)
  produce el mismo conjunto de 37.
- Filtro de ramas (POS) del propio K3 (cota 625/49 si
  5 | c; cotas 875/121 y 1375/169 si 5 ∤ c): supervivientes
  fuera de {1, 5, 7}: **ninguno**. El barrido, ahora sí
  ejecutado, confirma el cierre.
- La aritmética restante de K3 no dependía del bug y se
  reverificó: 875/121 = 7,2314 < 7,24; 1375/169 = 8,1361
  < 8,14; sup 125p/q² sobre pares consecutivos no excede
  1375/169; la igualdad solo en p = 11; 11p² − 125p + 44 > 0
  para todo primo p ≥ 13 hasta 10^5.

**Re-ejecución de la enumeración extendida.** El chequeo K4
del 125b (evaluador del doc 123a; no usa `sqfree_cs`, no
estaba afectado por el bug) se re-ejecutó desde el script
127 y reproduce su resultado en los tres marcos declarados:

| marco | c observados | N máx | violaciones |
|---|---|---|---|
| A, z = 61 | {1, 5, 7} | 3 | ninguna |
| A, z = 67 | {1, 5, 7} | 3 | ninguna |
| B, z = 61 | {1, 5, 7} | 3 | ninguna |

**Efecto en los veredictos.** El cierre del estrato 5 ya era
[T] por la derivación analítica vía (POS) (confirmada por la
auditora y por la aritmética de K3); lo que estaba roto era
solo el respaldo [M] del barrido de candidatos. Tras la
reparación, el [M] queda **restaurado de verdad**: la
enumeración corre, genera los 37 candidatos y no deja
supervivientes. Ningún veredicto cambia; una afirmación de
reproducibilidad que era falsa («el barrido pasó») pasa a
ser verdadera.

---

## 5. El contraejemplo exacto de la auditora, verificado con las funciones de la serie

**Aporte (126, §2.2).** La auditora refuerza la falsedad de
la anulación uniforme E_Rem = 0 con una instancia exacta que
no depende de la enumeración exhaustiva del 125b: z = 61,
s = 215077 (primo admisible), L = ⌊61⁵/s⌋ = 3926, r = 11,
d′ = 35 ∈ Rem(11), y como L/r ≤ 1331 el Lema 4 fuerza el
signo: E_Rem(11) ≤ −1/φ(35) = −1/24 < 0.

**Verificación propia [M sobre funciones de la serie,
chequeo R5]** — todo con la maquinaria de los docs 120a/125b
(`Level110`, `band_instances`, `lados`, `rem_info`,
`merged_failure`), sin aceptar nada de palabra:

- 61⁵ = 844596301; 215077 es primo; ⌊61⁵/215077⌋ = 3926. ✓
- r = 11 es visible y es instancia de banda en el nivel:
  u_r = log(3926/11)/log 61 = 1,42974 ∈ [1, 3/2);
  Π(11) = 0,393814 < 1/2 (el cociente exacto del doc 126,
  151963425/385875968, coincide a < 10⁻⁹); n̂ = 4 > 1, con
  testigo no degenerado 7: 11·7³ = 3773 < 3926. ✓
- `lados(11)` clasifica d′ = 35 como **Rem**; en el modelo
  V̂⁽¹¹⁾ la única condición impar es 11·7³ = 3773 < 3926,
  así que 35 ∈ V̂⁽¹¹⁾. ✓
- El fallo en la lista fusionada de 11·35 = 385 (primos
  11, 7, 5) es exactamente el anunciado:
  `merged_failure` devuelve (m, j, p_m) = (3, 2, 5), con
  5²·385 = 9625 ≥ 3926: la condición de Nathanson falla en
  la tercera posición, con j = 2 par ≥ 2, como manda el
  Lema 2. Luego 35 ∈ Rem(11). ✓
- `rem_info(11)`: |Rem(11)| = 1, Rem(11) = {35} (la forma
  fuerte de enumeración del 125b también se reproduce en
  esta instancia), y

$$E_{Rem}(11) = -\frac1{\varphi(35)} = -\frac1{24}
= -0,0416666\ldots < 0. ✓$$

**Incorporación.** La falsedad de la anulación uniforme
queda registrada como **[T-instancia]**: una sola instancia
exacta, verificable a mano y por la maquinaria, basta para
refutar el enunciado uniforme — no hace falta confiar en la
enumeración del 125b (que de todos modos la reproduce). El
rescate de signo del 125b no se toca: esta instancia vive en
L/r = 356,9 ≤ 1331, donde el teorema dice E_Rem ≤ 0, y el
contraejemplo lo satura con igualdad en su único elemento.

---

## 6. Estado final del doc 125 tras la errata

**Veredictos del doc 125: sin cambios.** Sobreviven tal como
se publicaron: (PRIME-WINDOW) refutada [T] y (OUT-MASS)
absoluta caída [T] en u_L → 2; lemas de paridad con
E_Rem = 0 (L/r ≤ 343) y E_Rem ≤ 0 (L/r ≤ 1331) [T];
anulación uniforme falsa [T-instancia, §5]; estrato 5
cerrado [T vía (POS)] con respaldo [M] restaurado (§4);
(EREM-SIGN) abierto solo para L/r > 1331.

**Texto del doc 125 que debe leerse corregido:**

| lugar | impreso | corregido |
|---|---|---|
| 125a §2 | M′ = M − 3/2 + Σ_{q≥5} 1/(q(q−1)); «resta 1 + 1/2» | M′ = M − 5/6 + Σ_{p≥5} 1/(p(p−1)) (§2) |
| 125a §3 | PW(L) = … + O(1/log z), u_L ∈ (3/2, 3) | mismo enunciado solo para 3/2 < u_L < 9/4 (§1) |
| 125a §3 | Σ capa₁ = S(L^{1/3}) − S(L z^{−3/2}) + ε_Π | … − ε_Π (§3) |
| 125b, script K3 | `if p < start or cur*p >= B: break` | `if p < start: continue` + `break` por tamaño (§4) |

**Inventario (H2) tras la errata: idéntico al del doc 125**
— (FC-2) [P]; (NET-LAYER) + (P-MASS) [P], la ruta con signo
de la banda alta; (EREM-SIGN) [P] solo z > 121;
(BAND-LOW-NONDEG) [P]; (SLAB-V′) [P] con familia (POS) y
estrato 5 [T] — salvo dos precisiones de escritura: el
dominio de la fórmula PW es (3/2, 9/4) y las constantes de
la §2 del 125a son las corregidas en la tabla. La cadena del
doc 115, (DIST-𝒟) + (H2) ⇒ UNIF-U-COND, no se toca.

**Goldbach sigue abierto.** Este doc corrige la escritura y
la reproducibilidad del 125; no demuestra ningún residuo
nuevo ni la conjetura.
