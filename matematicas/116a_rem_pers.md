# Doc 116a — (REM-PERS-small) es vacía: cierre por paridad y por tamaño

Sección 116a del doc 116. Ataca únicamente el residuo **(REM-PERS-small)**
del doc 114 §4–§5: los x ∈ Rem⁻ cuyos fallos fusionados de r·x caen
*todos* en primos ≤ q_min(x), incluido el caso ω(x) = 1. Notación y
definiciones: las de los docs 113–114 (posiciones contadas desde el
primo mayor: el primo propio q_m ocupa la posición propia m; la
condición del modelo en la posición m es r·q_m²·P_m < L con P_m el
producto de los primos propios ≥ q_m; la condición fusionada en p es
p²·(producto de los primos fusionados ≥ p) < L).

## 0. Veredicto en breve

**(REM-PERS-small) está CERRADO [T]: la familia residual es vacía,
siempre, no solo en ventana.** Para todo x ∈ V̂⁽ʳ⁾ con ω(x) impar
(en particular todo x ∈ Rem⁻), *ningún* fallo fusionado de r·x puede
caer en un primo ≤ q_min(x):

- en q_min(x) mismo: imposible por **paridad** (Lema A, §1);
- en un primo p | r con p < q_min(x): imposible por **tamaño**
  (Lema B, §2);
- en un primo p | r con p = q_min(x): imposible por coprimalidad.

Por tanto todo fallo de todo x ∈ Rem⁻ cae en un primo > q_min(x), la
persistencia del doc 114 [T] aplica **sin excepciones**, ψ(x) ∈ Rem⁺
para todo x ∈ Rem⁻, y la contribución de la familia residual a
E_Rem es exactamente 0 — no hace falta ni cota ni apareamiento
alternativo. Corolario [T]: **Rem no contiene ningún primo**
(Lema C, §4). (DESP-SIGN) queda reducido, tras los docs 114+116a,
únicamente a (R-FLOW-71+); en todo marco con z_A ≤ 71, (DESP-SIGN)
es ahora [T] completo. (H2) no se declara probada; Goldbach tampoco.

Verificación numérica (script `116a_rem_pers.py`, PASS, J0–J5,
2,4 s): las 4298 instancias del doc 114 (4113 de marcos + 185 de
estrés z=53), 5568 + 498 elementos x ∈ V̂⁻ examinados, 0 fallos en
primos ≤ q_min, 0 elementos residuales, 0 primos en Rem, aterrizaje
de ψ 39/39 + 36/36 — cruce exacto con K4/K7 del doc 114.

## 1. Lema A [T] — q_min(x) no puede fallar (paridad)

Sea x ∈ V̂⁽ʳ⁾ con ω(x) impar. El único primo propio de x que es
≤ q_min(x) es q_min(x) mismo, y su posición propia es ω(x): la
última, **impar** por hipótesis. Por el teorema R1 del doc 114
§1.2 [T], todo primo propio que falla la condición fusionada lo
hace en posición propia **par** (si fallara en posición propia
impar m, su lado izquierdo fusionado sería ≤ r·q_m²·P_m y el
modelo fallaría en m, contra x ∈ V̂⁽ʳ⁾). Luego q_min(x) no puede
ser un fallo fusionado de r·x. ∎

El argumento vale también para ω(x) = 1: la posición propia de x
primo es 1, impar, y el mismo teorema lo excluye como fallo R1.

## 2. Lema B [T] — ningún primo de r bajo q_min(x) puede fallar (tamaño)

Sea x ∈ V̂⁽ʳ⁾ con ω(x) impar, y sea p | r con p < q_min(x).
Afirmación: p no falla la condición fusionada de r·x.

*Prueba.* Como ω(x) es impar, la posición propia ω(x) —la de
q_min(x)— lleva condición del modelo, y como P_{ω(x)} = x (el
producto de *todos* los primos propios), que x pase el modelo da

$$r\cdot x\cdot q_{\min}(x)^{2} \;<\; L. \tag{II}$$

Los primos fusionados ≥ p son: todos los primos de x (todos ellos
son ≥ q_min(x) > p) y los primos de r que son ≥ p. El producto
parcial fusionado en p es por tanto

$$x\cdot\prod_{\substack{q\mid r\\ q\ge p}} q
\;=\; \frac{x\cdot r}{\mathrm{rad}_{<p}(r)},\qquad
\mathrm{rad}_{<p}(r):=\prod_{\substack{q\mid r\\ q<p}} q \;\ge\; 1.$$

Si p fallara la condición fusionada, se tendría
L ≤ p²·x·r/rad_{<p}(r) ≤ p²·x·r. Combinando con (II):

$$r\cdot x\cdot q_{\min}(x)^{2} \;<\; L \;\le\; p^{2}\cdot r\cdot x
\quad\Longrightarrow\quad q_{\min}(x)^{2} < p^{2}
\quad\Longrightarrow\quad q_{\min}(x) < p,$$

contra la hipótesis p < q_min(x). ∎

Tres observaciones de alcance:

1. El argumento es puramente de tamaño: no usa la paridad de la
   posición fusionada de p ni que x ∈ Rem — vale para todo
   x ∈ V̂⁽ʳ⁾ con ω(x) impar, esté en Rem o en Both.
2. El factor rad_{<p}(r) solo refuerza la desigualdad (hace el
   lado fusionado más pequeño); la prueba ya cierra con rad = 1,
   es decir, aunque r no tenga ningún primo bajo p salvo p mismo…
   en cuyo caso p < q_min(x) sigue siendo contradictorio con (II).
3. Es el espejo exacto, del lado Rem, del lema de holgura del
   doc 114 §1.1: allí un fallo A2 exigía un primo de r *bajo* el
   primo que falla; aquí un fallo R2 bajo q_min queda ahogado por
   la condición del modelo en la última posición propia, que el
   modelo sí ve (posición impar) aunque no vea los primos de r.

## 3. Teorema [T] — la familia residual es vacía

Sea x ∈ Rem⁻ (luego x ∈ V̂⁽ʳ⁾ y ω(x) impar). Todo fallo fusionado
de r·x es de tipo R1 (primo propio) o R2 (primo de r), sin mixtos
(doc 114 §1.2). Si algún fallo cayera en un primo p ≤ q_min(x):

- R1: p es propio y ≤ q_min(x), luego p = q_min(x) — excluido por
  el Lema A.
- R2: p | r y p ≤ q_min(x); p ≠ q_min(x) por coprimalidad
  ((x, r) = 1), luego p < q_min(x) — excluido por el Lema B.

Por tanto **todo** fallo fusionado de r·x cae en un primo
> q_min(x) (y hay al menos uno, pues x ∈ Rem). La persistencia
del doc 114 §4 [T] —al quitar q_min(x), la posición fusionada y
el producto parcial de cualquier primo > q_min(x) no cambian—
aplica entonces a todo x ∈ Rem⁻ sin excepción:

$$\psi(x) = x/q_{\min}(x) \in \mathrm{Rem},\qquad
\mu(\psi(x)) = +1,\qquad \psi(x) \in \mathrm{Rem}^{+}.$$

La familia (REM-PERS-small) es **vacía [T]**. Su contribución a
E_Rem(r) es exactamente 0: no es una subfamilia pequeña que haya
que acotar, es que no existe. La identidad de flujo del lado Rem
del doc 114,

$$E_{Rem}(r) = -\sum_{y \in Rem^+} a(y)\,\bigl(1 - R_R(y)\bigr),$$

vale ahora [T] módulo únicamente la cota R_R(y) ≤ 1.

## 4. Lema C [T] — Rem no contiene primos (la cláusula ω = 1)

Sea x = q un primo en V̂⁽ʳ⁾, con r visible (r ∈ V⁺(L): r divide a
un elemento de V⁺(L) y V⁺(L) es cerrado por divisores, Lema 2 del
doc 114 [T]). Los fallos fusionados posibles de r·q:

- R1 en q: posición propia 1, impar — excluido por el teorema R1
  del doc 114 §1.2 [T].
- R2 en p | r, p < q: excluido por el Lema B (ω(x) = 1 es impar).
- R2 en p | r, p > q: el único primo de x es < p, así que no entra
  ni en el conteo de primos fusionados ≥ p ni en el producto
  parcial en p: la posición fusionada de p y su lado izquierdo en
  r·q coinciden *exactamente* con los de r solo. Si r·q fallara
  en p, r fallaría su propia condición, contra r ∈ V⁺(L).

Luego r·q no tiene ningún fallo fusionado: q ∈ Both. En
particular **ningún primo está en Rem** [T], y la cláusula
«ningún primo x ∈ Rem⁻» de (REM-PERS-small) queda vacía por
teorema, no por la observación constructiva del doc 114
(ψ(q) = 1 ∈ Both): el apareamiento ya no necesita ese caso
especial porque el caso no ocurre.

## 5. Verificación numérica (script, J0–J5, PASS)

Reutilizando la maquinaria del doc 114 por import (sin tocarla),
sobre las mismas instancias del doc 114:

- **J0 (cruce):** 4113 instancias de marcos + 185 de estrés z=53
  (construcción del K7); Rem⁻ = 39 (marcos) y 36 (estrés) —
  idénticos a K4/K7 del doc 114.
- **J1 (hipótesis):** r ∈ V⁺(L) en las 4298 instancias; la
  condición (II) r·q_min²·x < L se cumple para los 5568 (marcos)
  + 498 (estrés) elementos x ∈ V̂⁻ examinados, 0 violaciones.
- **J2 (núcleo):** los fallos fusionados de todos los x ∈ V̂⁻ son
  39 + 36 en total (exactamente uno por cada x ∈ Rem⁻, ninguno en
  Both⁻); **0 en primos ≤ q_min(x)**: 0 en q_min mismo (Lema A),
  0 de tipo R2 con p < q_min (Lema B).
- **J3 (Lema B directo):** para todo x ∈ V̂⁻ y todo p | r con
  p < q_min(x) —2423 pares en marcos, 221 en el estrés— el lado
  fusionado p²·x·r/rad_{<p}(r) queda < L siempre; la razón máxima
  contra L es 0,7581 (marcos) y 0,4916 (estrés): el ahogo del
  Lema B tiene margen, no es un empate.
- **J4 (familia):** elementos residuales (x ∈ Rem⁻ con todos sus
  fallos ≤ q_min): **0**. Primos en Rem: **0**. Todo x ∈ Rem⁻
  tiene ≥ 1 fallo, todos > q_min.
- **J5 (aterrizaje):** ψ(x) ∈ Rem⁺ para 39/39 (marcos) y 36/36
  (estrés), ahora explicado por el teorema de §3, no solo medido.

## 6. Estado exacto tras 116a

- **(REM-PERS-small): CERRADO [T]** — familia vacía, contribución
  idénticamente 0.
- (DESP-SIGN), tras los docs 114 + 116a, queda reducido [T] a un
  único residuo: **(R-FLOW-71+)**, las cotas R_A(y) ≤ 1 y
  R_R(y) ≤ 1 para y con min(y) ≥ 71. En todo marco con z_A ≤ 71,
  (DESP-SIGN) es [T] completo por esta vía.
- Nada de lo anterior toca (DIST-𝒟) ni cierra (H2) en el marco
  del puente (z_A → ∞): (H2) sigue abierta exactamente por
  (R-FLOW-71+). Goldbach no se declara demostrada.

## 7. Reproducibilidad

`python3 116a_rem_pers.py | tee salida_doc116a_local.txt` — solo
stdlib; importa `114_desp_sign_banda.py` (y a través de él los
docs 113/112/111/110/109) sin modificarlos. Salida verificada:
PASS, chequeos J0–J5, 2,4 s.
