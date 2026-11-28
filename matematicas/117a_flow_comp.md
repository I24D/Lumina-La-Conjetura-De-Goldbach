# Doc 117a — (FLOW-COMP): la masa mala contra la masa buena, por conductor

Sección acotada del doc 117: solo el residuo **(FLOW-COMP)** del
doc 116b §6 — la compensación global que hereda (DESP-SIGN) tras
la refutación puntual de (R-FLOW-71+). Notación: docs 113/114/116b
(a(y) = μ(y)/φ(y); R_A(y) = Σ_{q : yq ∈ Add⁺} 1/(q−1); la identidad
de flujo E_Add(r) = −Σ_{y∈Add⁻} |a(y)|(1 − R_A(y)) es [T], doc 114).

Descomposición por conductor r visible, y análoga en Rem:

  BAD(r)  = Σ_{y∈Add⁻} |a(y)|·(R_A(y) − 1)⁺,
  GOOD(r) = Σ_{y∈Add⁻} |a(y)|·(1 − R_A(y))⁺,
  GOOD(r) − BAD(r) = −E_Add(r)   [identidad, sin calcular R_A].

(FLOW-COMP)_Add es BAD(r) ≤ GOOD(r) para todo r; idem Rem.

## 0. Veredicto en breve

- **(FLOW-COMP) NO se cierra [T] en general, pero queda reducido
  [T] a un sub-enunciado explícito y estrictamente más débil,
  (FC-2)** (§4): una comparación cabeza/cola sobre la masa de
  Add⁻, sin flujos. (FC-2) usa solo ~15 % de la masa buena real
  en el marco decisivo y aun así sobra.
- **Prueba decisiva [M]: en el marco del contraejemplo del
  doc 116b (z = 331, r = 5, L = 79051268 ∈ ℒ_A), donde BAD > 0
  de verdad, (FLOW-COMP) vale con margen ×28**: BAD = 1,086·10⁻³,
  GOOD = 0,030399, E_Add(5) = −0,029314 < 0. Los 13 y malos son
  exactamente los primos 251 ≤ Y ≤ 317 de la configuración libre.
  El lado Rem en el mismo nivel: E_Rem(5) = −0,000556 (convención
  doc 114), BAD_R = 0.
- **La masa mala por r está cuantificada [T]** (§1): malo ⟹
  min(y) ≥ 73; en la configuración libre prima, Y vive en la
  ventana [(L/r₀)^{1/3}, L^{1/3}) y cada término decae como
  r₀^{1/3}·log log L / L^{1/3}; el total de la familia prima
  decae como (log log L)/log L (con PNT; sin PNT el conteo
  trivial no da decaimiento del total, solo del término).
- **La masa buena tiene margen uniforme [T]** (§2): todo
  y ∈ Add⁻ con min(y) ≤ 37 aporta (1 − R_A(y)) ≥ 1 − Mertens(37)
  = 2581/18480 = 0,139665. En los marcos A/B el 100 % de la masa
  de Add⁻ es cabeza ≤ 37; en z = 331 la cabeza (0,0321) supera
  a BAD incluso descontada por el margen (holgura 4,13).

## 1. Masa mala: cuantificación

**(1a) Umbral [T].** Si R_A(y) > 1 entonces min(y) ≥ 73: la cota
cruda del doc 114 da R_A(y) ≤ Σ_{5≤q<min(y)} 1/(q−1) ≤
Mertens(67) = 0,999174 < 1 si min(y) ≤ 67 (fracciones exactas,
K0). Lo mismo para R_R. Luego BAD(r) solo recibe masa de la
cola min(y) ≥ 73 de Add⁻.

**(1b) Configuración libre prima [T].** Sea r = r₀ primo e
y = Y primo, Y > r₀, en Add⁻. El modelo falla en la posición 1:
r₀·Y³ ≥ L; el fusionado cumple su posición 1: Y³ < L. Luego

  (L/r₀)^{1/3} ≤ Y < L^{1/3},

una ventana de ancho multiplicativo r₀^{1/3}. Cada término malo
vale |a(Y)|(R_A(Y) − 1) ≤ (Mertens(Y−1) − 1)/(Y−1) =
O(log log Y / Y) (teorema de Mertens), es decir
O(r₀^{1/3}·log log L / L^{1/3}): **decaimiento polinómico por
término**. El número de primos malos es ≤ π(L^{1/3}) −
π((L/r₀)^{1/3}); con PNT el total es
O((log log L)/log L) → 0. Sin PNT (π(x) ≤ x) el total no decae:
el decaimiento del total es exactamente tan fuerte como PNT en
la ventana cúbica. En el nivel decisivo (L = 79051268, r₀ = 5):
la ventana es Y ∈ [251,0…, 429,2…) ∩ {Y < z = 331}, y malo ⟺
R_A(Y) = Σ_{5<q<Y} 1/(q−1) > 1 ⟺ Y ≥ 251: los malos son
{251, 257, …, 317}, 13 primos, BAD = 1,086·10⁻³ [M, K2].

**(1c) Compuestos [T + M].** Si ω(y) ≥ 3 y min(y) ≥ 73,
φ(y) ≥ (min(y)−1)²·(max(y)−1) ≥ 72²·(max(y)−1), y
R_A(y) − 1 ≤ Mertens(min(y)−1): cada término malo compuesto es
≤ (log log z_A + O(1))/(72²·(max(y)−1)). En el marco decisivo
no hay **ningún** y compuesto malo: la identificación exacta
(K2) da el soporte de BAD = los 13 primos de (1b), sin resto.

## 2. Masa buena: margen uniforme en la cabeza

**(2a) Margen [T].** Si y ∈ Add⁻ y min(y) ≤ 37, la cota cruda del
doc 114 da R_A(y) ≤ Σ_{5≤q<min(y)} 1/(q−1) ≤ Mertens(37) =
15899/18480 = 0,860335 (fracción exacta, K0). Luego

  GOOD(r) ≥ (1 − Mertens(37)) · HEAD37(r) = 0,139665 · HEAD37(r),

con HEAD37(r) = Σ_{y∈Add⁻, min(y)≤37} |a(y)|. Es un margen
uniforme de verdad (no el 0,00083 que daría min(y) ≤ 67): la
cabeza pequeña de Add⁻ compensa con holgura constante.

**(2b) ¿Cuánta masa es cabeza? [M].** Marcos A/B (4113 instancias
del doc 114): la fracción mínima de masa de Add⁻ con min(y) ≤ 37
es 1,0000 — toda la masa es cabeza, y BAD = 0 idéntico (z ≤ 53).
En z = 331, r = 5, L = 79051268: HEAD37 = 0,0321 y la cola
min(y) ≥ 73 pesa 0,0462 — **la cola NO es pequeña**; lo que salva
la compensación no es que la cola pese poco, sino que su exceso
(R_A − 1)⁺ es cero salvo en los 13 primos libres y ≤ 0,0452 en
ellos. Aun usando solo la cabeza descontada:
0,139665 · 0,0321 = 0,00448 ≥ BAD = 0,001086 (holgura 4,13).

## 3. Medición decisiva y de estrés [M] (K1–K4, script PASS)

**K2 — el marco donde (FLOW-COMP) puede fallar de verdad**
(z = 331, r = 5; enumeración DFS exhaustiva de {y : 5y ∈ V⁺(L)}
y de V̂⁽⁵⁾, R_A exacto en fracciones sobre min(y) ≥ 73):

| L | E_Add(5) | BAD | GOOD | BAD/GOOD | holgura (FC-2) |
|---|---|---|---|---|---|
| 46948396 (control) | −0,036949 | 1,086·10⁻³ | 0,038035 | 2,85·10⁻² | 4,98 |
| 79051268 (cruce Y=251) | −0,029314 | 1,086·10⁻³ | 0,030399 | 3,57·10⁻² | 4,13 |
| 113296524 (cruce Y=283) | −0,022764 | 7,30·10⁻⁴ | 0,023494 | 3,11·10⁻² | 4,89 |

R_A(251) = 1,00220 reproducido; máx R_A = 1,04522 en y = 317.
Lado Rem en L = 79051268: E_Rem(5) = −0,000556 < 0, BAD_R = 0
(ningún y ∈ Rem⁺ con min(y) ≥ 73 en este marco), GOOD_R =
5,56·10⁻⁴. **(DESP-SIGN) sobrevive a su peor marco conocido,
por compensación, con margen ×28 en Add.**

**K1 — marcos A/B:** 4113 instancias (nivel, r): máx E_Add = 0,
máx E_Rem = −0,0 (ambos ≤ 0, [M]); BAD = BAD_R = 0 idéntico;
el mínimo de GOOD es 0 y se alcanza solo en instancias
degeneradas (Add⁻ vacío: no hay nada que compensar).

**K3 — estrés z = 101** (3 niveles de ℒ_A, 40 conductores con
ω(r) ≤ 2 por nivel): máx E_Add = 0, BAD = 0 en todo r;
máx R_A sobre min(y) ≥ 73: 0,8024 (L = 1050274) y 0,8476
(L = 10375222) — coherente con el J4 del doc 116b.

**K4 — el único término malo Rem construido (doc 116b, J3):**
|a(y)|(R_R − 1) = (Mertens(251) − 1)/φ(20531891) = 1,248·10⁻⁸
exacto. En Rem, BAD_R = 0 en todo marco medido; la enfermedad
puntual de R_R existe pero su masa es despreciable donde se ha
podido construir.

## 4. Veredicto: (FLOW-COMP) reducido a (FC-2)

**(FC-2)_A [P]:** para todo r visible,
Σ_{y∈Add⁻} |a(y)|·(R_A(y) − 1)⁺ ≤ (1 − Mertens(37)) ·
Σ_{y∈Add⁻, min(y)≤37} |a(y)|; y **(FC-2)_R**, el análogo en Rem.

- **[T] (FC-2) ⟹ (FLOW-COMP).** Por §2a, GOOD(r) ≥
  (1 − Mertens(37))·HEAD37(r) ≥ BAD(r); idem Rem con GOOD_R.
  Con las identidades de flujo del doc 114 y el cierre de
  (REM-PERS-small) del doc 116a, (FLOW-COMP) es exactamente el
  ingrediente que la ruta de apareamiento necesita para
  (DESP-SIGN): la cadena queda (FC-2) ⟹ (FLOW-COMP) ⟹
  (DESP-SIGN) por esta ruta.
- **(FC-2) es estrictamente más débil que (FLOW-COMP):** ignora
  toda la masa buena con 37 < min(y) < 73 y la de min(y) ≥ 73
  con R_A < 1 — en el nivel decisivo, el miembro derecho usa
  0,00448 de los 0,030399 de GOOD real (~15 %). Es un enunciado
  solo sobre la distribución de masa de Add⁻ (cabeza contra
  exceso de cola), sin flujos ni apareamientos.
- **[M] (FC-2) vale en todo lo medido:** trivialmente donde
  BAD = 0 (marcos A/B, z = 101) y con holgura 4,13–4,98 en los
  tres niveles de z = 331, los únicos con BAD > 0.
- **Obstrucción exacta para cerrarlo [T]:** haría falta una
  cota inferior uniforme de HEAD37(r) cuando L → ∞ en el marco
  del puente, comparable al decaimiento (log log L)/log L del
  total de la familia prima mala (§1b). Ninguna de las dos
  cantidades está probada en ese régimen: HEAD37(r) > 0 es
  plausible (la cabeza es el grueso de Add⁻ en todo marco
  medido) pero no es teorema. Ahí, y solo ahí, vive lo que
  falta de (FLOW-COMP).

Estado: (DESP-SIGN) depende por esta ruta de (FC-2) solamente;
(BAND-δ), el otro residuo promediado de (H2) según el doc 116,
no se toca aquí. Goldbach no se declara demostrada.

## 5. Reproducibilidad

`117a_flow_comp.py` (solo stdlib; importa los docs 109/114/116b
sin tocarlos): K0 Mertens exactos; K1 descomposición BAD/GOOD
sobre las 4113 instancias del doc 114 vía `instancias()`; K2
enumeración DFS exhaustiva en z = 331 con R_A/R_R exactos
(fracciones) sobre la cola min(y) ≥ 73 e identificación del
soporte de BAD; K3 estrés z = 101; K4 término exacto del J3 del
doc 116b. Salida: `salida_doc117a_local.txt` — **PASS**, 0,3 s.
