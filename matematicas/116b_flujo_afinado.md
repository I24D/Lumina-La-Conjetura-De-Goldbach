# Doc 116b — (R-FLOW-71+): cotas afiladas, y el enunciado puntual es falso

Sección acotada del doc 116: solo el residuo **(R-FLOW-71+)** del
doc 114 §5 — las razones de flujo
R_A(y) = Σ_{q : yq ∈ Add⁺} 1/(q−1) y R_R(y) análoga debían ser ≤ 1
para todo y con min(y) ≥ 71, para cerrar (DESP-SIGN) por el
apareamiento φ(x) = x/q_min(x). Notación y definiciones: docs
109/113/114 (𝒫 = primos impares < z_A, ∤ A; A = 3, luego 𝒫 parte
de 5; V⁽ʳ⁾, V̂⁽ʳ⁾, Add, Rem como en el doc 114).

## 0. Veredicto en breve

- **(R-FLOW-71+) es FALSO como enunciado puntual [T].** Hay
  contraejemplos *dentro del marco* (nivel L ∈ ℒ_A, r visible,
  verificados con los predicados de la serie): R_A(y) = 1,00220
  con min(y) = 251 (J2), y R_R(y) = 1,25620 con min(y) = 4099
  (J3). Más aún, R_A(y) → ∞ como log log Y en la familia del
  contraejemplo: ninguna cota puntual de este flujo puede valer.
- **(DESP-SIGN) NO queda refutado.** La identidad de flujo da
  E_Add(r) = −Σ_y |a(y)|(1 − R_A(y)): un y con R_A(y) > 1 aporta
  un término positivo, pero el signo total puede salvarse por
  compensación con los demás y. Lo que muere es la *suficiencia
  puntual*; el residuo exacto pasa a ser una compensación global
  (FLOW-COMP), §6.
- **Las cotas afiladas pedidas sí existen y son teoremas [T]**
  (§1): la condición fusionada en la posición del propio q da la
  cota cúbica (★); los primos de r volteados por la inserción de q
  dan (★★), con su caso limpio (★★★). Verificadas sin una sola
  violación en los 1640 pares admisibles del doc 114 (J1). Ellas
  *caracterizan* la obstrucción (§2): el contraejemplo es la única
  configuración que las esquiva a todas en rango.
- **Medición [M] (J4):** en z = 101, sobre la muestra de niveles
  y r con ω(r) ≤ 2, el máximo de R_A crece con min(y):
  0,8982 (71) → 0,9125 (73) → 0,9263 (79) → 0,9392 (83) →
  0,9514 (89, 97) — siempre < 1 en ese rango, como debe ser: la
  Mertens disponible a esa escala no alcanza. La violación es un
  fenómeno de escala mayor, no de ventana.

## 1. Admisibilidad exacta de q y las tres cotas afiladas [T]

Sea y ∈ Add⁻ (μ(y) = −1, ω(y) impar), r > 1, (y, r) = 1, y sea
q ∈ 𝒫 con yq ∈ Add⁺. Como φ(yq) = y exige q = q_min(yq):

(i) q < min(y), (q, ry) = 1; (ii) ryq < L; (iii) r·y·q ∈ V⁺(L)
(condiciones fusionadas); (iv) yq ∉ V̂⁽ʳ⁾ — automático: las
condiciones del modelo de yq en las posiciones ≤ ω(y) son las de
y (q ocupa la posición ω(y)+1, par), y y ya falla el modelo.

En la lista fusionada descendente de ryq, los primos > q ocupan
las mismas posiciones y productos parciales que en ry: sus
condiciones se heredan [T]. Las únicas condiciones *nuevas* son
la de la posición de q y las de los primos de r menores que q,
cuya posición sube en 1 y cuyo producto parcial gana el factor q.
Con t(q) = #{p | r : p > q} y s_q = ∏_{p|r, p<q} p:

**(★) [T]** Si t(q) es impar, la posición fusionada de q,
1 + ω(y) + t(q), es impar (ω(y) impar), y su condición da
q² · (ryq/s_q) < L, es decir
**q³ < L·s_q/(r·y) ≤ L/y**.

**(★★) [T]** Si p | r, p < q, y #{p' | r : p' > p} es par, entonces
p ocupa posición par en ry (su posición es 1 + ω(y) + #{p' > p})
y posición impar en ryq, con producto parcial nuevo q·y·r_{≥p}
(todos los primos de y superan a p; r_{≥p} = ∏_{p'|r, p'≥p} p').
Su condición da **q < L/(p²·y·r_{≥p})**.

**(★★★) [T]** En particular, el primo mayor de r, ρ = max(r),
cumple #{p' > ρ} = 0 (par) siempre: todo q admisible con
ρ < q < min(y) verifica **q < L/(ρ³·y)**.

Verificación [M] (J1): en las 4113 instancias del doc 114, los
1640 pares (y, q) admisibles cumplen (★), (★★) y (★★★) sin una
sola violación, y el máximo de R_A en ventana se reproduce:
0,6000. Las tres cotas son condiciones necesarias reales, no
artefactos.

## 2. La configuración libre: la obstrucción exacta [T]

Un q admisible escapa de (★) si t(q) es par, y escapa de (★★)
dentro del rango q < min(y) si el primo volteado que lo ata da
una cota mayor que min(y). La configuración extrema es:

**r = r₀ primo fijo, y = Y primo.** Para q ∈ (r₀, Y): t(q) = 0
(par: sin (★)); el único primo de r bajo q es r₀, que sí ata por
(★★) — su posición en ry es 1 + ω(y) = 2, par — pero la cota es
q < L/(r₀³·Y), y con L en la ventana de abajo vale
L/(r₀³·Y) ≈ Y²/(2r₀²) > Y: el atado existe y es vacío en rango.
Para q < r₀ no hay primos en 𝒫 (𝒫 parte de 5 ≤ r₀ si r₀ = 5;
en general q = r₀ queda excluido por (q, r) = 1 y los q < r₀
suman ≤ Mertens(r₀ − 1), fijo).

La ventana de niveles que realiza la configuración, para r₀ = 5:

  max(Y³, 125·Y·q_max) < L ≤ 5Y³,   q_max = mayor primo < Y.

- L ≤ 5Y³: el modelo de y (e de yq) falla en la posición 1:
  r·Y²·Y = 5Y³ ≥ L ⟹ y, yq ∉ V̂⁽ʳ⁾.
- L > Y³: la condición fusionada de ry (e de ryq) en la posición
  1 se cumple: Y³ < L.
- L > 125·Y·q: la condición fusionada de ryq en la posición 3
  (el primo r₀ = 5, volteado) se cumple para todo q ≤ q_max:
  25·(Y·q·5) < L. La posición de q es 2 (par): sin condición.
- Tamaño: ryq = 5Yq < L, heredado de la anterior.

Luego **todo** q ∈ 𝒫 con 5 < q < Y es admisible, y

  R_A(y) = Σ_{5<q<Y} 1/(q−1) = Mertens(Y−1) − 1/4 → ∞

por el teorema de Mertens para Σ 1/(p−1). La desigualdad puntual
R_A(y) ≤ 1 es por tanto falsa, y no por poco: es falsa con
margen creciente. El fallo no está en el apareamiento (la
identidad de flujo es exacta [T], doc 114) sino en pedir la
cota término a término.

## 3. Contraejemplo R_A en marco [T] (J2)

Marco A del todo concreto: z = 331, D = z⁵, A = 3, nivel
L = ⌊D/s⌋ ∈ ℒ_A con s primo ≤ z^{10/3}; r = 5 (visible en cada
nivel: 5Y ∈ V⁺(L), verificado); y = Y primo. Membresías por los
predicados de la serie (in_Vplus doc 109, in_Vhat doc 113,
convención de lados() doc 114); R_A por suma exacta de
fracciones, y coincide con la forma cerrada Σ_{5<q<Y} 1/(q−1)
en cada fila:

| Y | s | L | R_A(y) |
|---|---|---|---|
| 73 | 2042687 | 1945083 | 0,76346 |
| 101 | 771269 | 5151504 | 0,82415 |
| 151 | 230807 | 17214364 | 0,90856 |
| 211 | 84629 | 46948396 | 0,97144 |
| **251** | **50261** | **79051268** | **1,00220** |
| 283 | 35069 | 113296524 | 1,02855 |
| 307 | 27479 | 144590262 | 1,03552 |

Primer cruce: **Y = 251**, R_A(y) =
22633668515439323928857875/22584054598420241438576064 = 1,00220,
con los 50 primos de (5, 251) todos admisibles. Como la familia
continúa (R_A ~ log log Y), (R-FLOW-71+) no se repara subiendo
el umbral 71: para todo umbral existe Y que lo viola.

## 4. Contraejemplo R_R en marco [T] (J3)

El lado Rem admite la misma enfermedad por otra puerta: el fallo
fusionado no necesita estar en la posición 1. Sea r = ρ = 4999
(primo), y = Y₁·Y₂ con Y₁ = 5009, Y₂ = 4099 (μ(y) = +1),
y = 20531891, min(y) = 4099. Con L en la ventana
(ρ·y·251³, ρ·Y₁·Y₂³]:

- El modelo de y vale: posición 1, ρ·Y₁³ ≈ 6,3·10¹⁴ < L.
- El fusionado de ry falla en la posición 3 (el primo Y₂):
  Y₂²·(Y₁·ρ·Y₂) = ρ·Y₁·Y₂³ ≥ L. Luego y ∈ Rem⁺.
- Para cada q ∈ 𝒫, q ≤ 251: yq sigue en el modelo (posición 3:
  ρ·q³·y < L por la ventana) y el fusionado de ryq sigue
  fallando en la misma posición 3 (el producto parcial hasta Y₂
  no incluye a q). Luego yq ∈ Rem⁻: **todo** q ≤ 251 es
  admisible, sin exclusiones ((q, ρy) = 1 automático en el rango).

Instancia en marco: z = 6052, D = z⁵, s = 5,
L = ⌊D/5⌋ = 1623770295721736806 ∈ ℒ_A; ρ visible
(Y₁·ρ ∈ V⁺(L), verificado). Resultado exacto:

  R_R(y) = Σ_{5≤q≤251} 1/(q−1) = Mertens(251) = 1,25620 > 1,

con 52 primos admisibles. R_R(y) ≤ 1 también es falso como
enunciado puntual en marco.

## 5. Medición extendida [M] (J4)

Con los predicados locales (sin construir niveles completos):
z = 101, 10 niveles (8 muestreados de ℒ_A + 2 dirigidos a la
familia), todos los r con ω(r) ≤ 2 y todos los y libres de
cuadrados con min(y) ≥ 71 (completo sobre esa muestra):
1542 instancias Add⁻ y 242 Rem⁺. Máximos por cubeta:

| min(y) | máx R_A | máx R_R |
|---|---|---|
| 71 | 0,8982 | 0,4167 |
| 73 | 0,9125 | 0,4167 |
| 79 | 0,9263 | 0,4167 |
| 83 | 0,9392 | 0,4167 |
| 89 | 0,9514 | 0,4167 |
| 97 | 0,9514 | 0,4167 |

(los máximos de R_A en r = 437 = 19·23, L = 104060401, y primo).
En z = 71 (6 niveles, r con ω(r) ≤ 2, enumeración de y con
ω(y) ≤ 6 por nivel y r; el presupuesto DFS no se agotó): máx R_A = 0,8493
(L = 25411681, r = 323, y = 61), máx R_R = 0,4167.

Lectura: la tendencia con min(y) es creciente y compatible con
el mecanismo de Mertens de §2; en z ≤ 101 nadie viola, porque
la suma disponible (≤ Mertens(97) en la familia pura) aún no
alcanza 1 tras las exclusiones. La medición del doc 114
(máx 0,6625 en z = 53) y esta forman una sola curva que el
contraejemplo de §3 continúa hasta cruzar en min(y) = 251.

## 6. Estado exacto y consecuencia

- **(R-FLOW-71+): cerrado — como FALSO [T].** No es un residuo
  abierto: es un enunciado refutado, con la obstrucción exacta
  identificada (la configuración libre de §2: y primo, r primo
  pequeño, q en posición fusionada par y atado de volteo vacío
  en rango; en Rem, fallo fusionado heredado en posición ≥ 3).
- **(DESP-SIGN): sigue abierto**, y su ruta por apareamiento
  necesita ahora la compensación global. El enunciado que
  hereda el residuo es, para cada r visible:
  **(FLOW-COMP)** Σ_{y ∈ Add⁻} |a(y)|·(R_A(y) − 1)⁺ ≤
  Σ_{y ∈ Add⁻} |a(y)|·(1 − R_A(y))⁺, y su análogo en Rem —
  es decir, que cada término de flujo de la descomposición del
  doc 114 sea ≤ 0 sin pasar por la cota puntual. En ventana la compensación tiene holgura
  estructural (doc 114: masa Add⁺ 18,41 contra 93,29 socia),
  pero eso es [M], no prueba. (FLOW-COMP) es estrictamente más
  débil que (R-FLOW) y no es comparable con (REM-PERS-small)
  ni con (BAND), que quedan como estaban.
- Nada de esto toca (DIST-𝒟) ni la otra mitad de (H2).
  Goldbach no se declara demostrada.

## 7. Reproducibilidad

`116b_flujo_afinado.py` (solo stdlib; importa los docs 109/114
sin tocarlos): J0 reproduce el umbral de Mertens del doc 114 con
fracciones exactas (Σ_{5≤p≤67} = 160107799/160240080 = 0,99917);
J1 verifica (★)/(★★)/(★★★) en los 1640 pares del doc 114
(0 violaciones; máx R_A = 0,6000 reproducido); J2 escanea la
familia de §2 en marco y halla el primer cruce; J3 verifica la
instancia de §4; J4 es la medición de §5. Salida:
`salida_doc116b_local.txt` — PASS, 1,1 s.
