# Doc 117 — Los dos promedios del doc 116: (FLOW-COMP) reducido a (FC-2) y probado en su contraejemplo; (BAND-δ) reducido a (BAND-MOD) y de ahí a (BAND-LOW)

**Fecha:** 2026-11-28. **Línea:** UNIF-U / (H2), continuación de los docs 113–116.
**Documentos de prueba (este doc es la síntesis):** `117a_flow_comp.md` (script PASS, K0–K4, 0,3 s) y `117b_banda_delta.md` (script PASS, K0–K6, 0,6 s).
**Etiquetas:** [T] = teorema, [M] = medición en ventana, [P] = pendiente.

---

## 0. Veredicto en breve

El doc 116 dejó (H2) dependiendo de dos enunciados promediados: **(FLOW-COMP)** y **(BAND-δ)**. Este doc ataca ambos:

1. **(FLOW-COMP) queda reducido [T] a (FC-2)** (§1, doc 117a), un enunciado cabeza/cola sobre la masa de Add⁻ estrictamente más débil — usa solo ~15 % de la masa buena real — y **sobrevive su prueba decisiva [M]**: en el marco del contraejemplo que mató la forma puntual (z = 331, r = 5), la compensación vale con margen ×28. El cierre [T] de (FC-2) está bloqueado por exactamente dos piezas nombradas: una entrada de fuerza PNT para la masa mala y una cota inferior uniforme de la masa de cabeza HEAD37(r).
2. **(BAND-δ) queda cerrado módulo (DESP-SIGN) + (BAND-MOD)** (§2, doc 117b), con la reducción [T] puntual y acíclica. (BAND-MOD) no tiene forma cerrada con signo [T] ni cancelación que explotar (en la banda Ŝ−Π es unilateral [M]), y queda reducido a **(REG-Π)** — con ruta de prueba identificada y medición plana — más **(BAND-LOW)** [P], la capa u_r < u₀*.

**(H2) sigue abierta**, pero su inventario de residuos es ahora: **(FC-2)** (con sus dos bloqueadores nombrados), **(REG-Π)** y **(BAND-LOW)**. Goldbach no se declara demostrado.

---

## 1. (FLOW-COMP) → (FC-2), y la prueba del contraejemplo (doc 117a)

**Reducción [T].** (FC-2) ⟹ (FLOW-COMP), vía el margen uniforme 1 − Mertens(37) = 0,139665: la cabeza de Add⁻ (min(y) ≤ 37) ya aporta, por la cota de Mertens exacta del doc 114, un factor (1 − R_A) ≥ 0,139665 sobre su masa; (FC-2) pide solo que esa fracción de la masa buena total baste contra la cola mala — consume ~15 % de la masa buena medida, el resto es holgura sin usar.

**Masa mala [T].** Malo ⟹ min(y) ≥ 73. En la familia prima libre, la ventana de y es [(L/r₀)^{1/3}, L^{1/3}), el término es O(r₀^{1/3}·log log L / L^{1/3}) y el total es O((log log L)/log L) **con PNT**; sin una entrada de fuerza PNT el total no decae — esa es una de las dos obstrucciones exactas del cierre [T]. La otra: no hay cota inferior uniforme para HEAD37(r) (la masa de cabeza puede, en principio, adelgazar en algún r exótico).

**Prueba decisiva [M].** En el marco del contraejemplo del doc 116b (z = 331, r = 5, L = 79051268 ∈ ℒ_A): BAD = 1,086·10⁻³, llevado por exactamente los 13 primos libres 251 ≤ Y ≤ 317 (máx R_A = 1,04522 en Y = 317; R_A(251) = 1,00220 reproducido al dígito); GOOD = 0,030399; E_Add(5) = −0,029314. **Margen ×28.** Lado Rem: E_Rem = −0,000556, BAD_R = 0. La holgura de (FC-2) en los tres niveles: 4,13–4,98. En los marcos A/B (4113 instancias) y en el estrés z = 101, BAD = 0 idéntico: la compensación es trivial fuera del régimen del contraejemplo. El único término malo Rem construido en la serie (doc 116b, J3) pesa 1,248·10⁻⁸.

## 2. (BAND-δ) → (BAND-MOD) → (REG-Π) + (BAND-LOW) (doc 117b)

**Reducción [T].** Si E(r) ≤ 0 (eso es (DESP-SIGN)), entonces δ(r)⁺ ≤ (Ŝ(r) − Π(r))⁺ puntualmente, luego (BAND-δ) ⟸ (DESP-SIGN) + **(BAND-MOD)**: Σ_{r∈banda} (Ŝ(r) − Π(r))⁺/φ(r) ≤ C·Q_c. La reducción es acíclica: (DESP-SIGN) descansa sobre (FLOW-COMP) y (REM-PERS-small), que no tocan la banda, y (BAND-MOD) no involucra E.

**Sin forma cerrada, sin cancelación [T/M].** La suma con signo Σ_r (Ŝ(r) − Π(r))/φ(r) no hereda el intercambio del doc 110: la suma interna es una suma de Mertens truncada por el umbral r·λ(d′) ≥ L, y el peso W cambia con L (1,2702 → 0,2702 en ventana). Y no hay cancelación que explotar en su lugar: sobre la banda, Ŝ − Π es **unilateral** [M: 246/246 instancias positivas; la parte positiva es el valor absoluto]. (BAND-MOD) es genuinamente una cota de tamaño, no de cancelación.

**Reducción de (BAND-MOD).** Sobre {u_r ≥ u₀*} cierra [T] por la FL clásica del modelo, módulo **(REG-Π)** [P con ruta]: Σ_banda Π(r)/φ(r) ≤ C·Q_c — medido 0,17–0,46·Q_c sin crecimiento en z, y la ruta de prueba son las regiones de Mertens del doc 116c. El residuo mínimo es **(BAND-LOW)** [P]: la capa u_r < u₀*. Dentro de ella, la subcapa degenerada Ŝ = 1 quedó caracterizada [T]: V̂⁽ʳ⁾ = {1} ⟺ q₀³·r ≥ L (q₀ = menor primo de 𝒫∖r); son 188/246 instancias en ventana y su región límite tiene **medida cero** — la dicotomía del doc 116c apunta a favor; falta la cota de densidad. Fuera de la subcapa degenerada, el exceso es O(Π(r)) [M: máx (Ŝ−Π)⁺/Π = 0,692], que es exactamente lo que (REG-Π) absorbería.

## 3. Veredicto de cadena actualizado

- (H2) ⟸ (FL-DESP⁺) ⟸ (DESP-SIGN) + banda [docs 109–114].
- (DESP-SIGN) ⟸ (FLOW-COMP) + (REM-PERS-small) [114, 116b]; (REM-PERS-small) **cerrado [T]** [116a].
- (FLOW-COMP) ⟸ **(FC-2)** [117a, T]; (FC-2) bloqueado por: (i) entrada de fuerza PNT para el decaimiento de la masa mala, (ii) cota inferior uniforme de HEAD37(r).
- Banda ⟸ (BAND-δ) ⟸ (DESP-SIGN) + (BAND-MOD) [117b, T]; (BAND-MOD) ⟸ **(REG-Π)** + **(BAND-LOW)** [117b].
- Con el doc 115: **(DIST-𝒟) + (H2) ⇒ UNIF-U-COND**.

**Inventario final de residuos de (H2): (FC-2), (REG-Π), (BAND-LOW).** Los tres son enunciados explícitos y pequeños; ninguno tiene contraejemplo; dos tienen ruta de prueba identificada ((REG-Π) por regiones de Mertens; (BAND-LOW) por medida cero de su región límite). El patrón de la línea se mantiene: cada residuo que muere deja uno estrictamente menor y mejor nombrado.

Siguen fuera, sin cambios: (DIST-𝒟) en la línea de Lumina Chat, y la separación de capas G0. **Goldbach permanece abierto.**

## 4. Reproducibilidad

Pruebas y verificaciones en `117a_flow_comp.py` (PASS, K0–K4) y `117b_banda_delta.py` (PASS, K0–K6), con salidas `salida_doc117a_local.txt` y `salida_doc117b_local.txt`. Este doc es la síntesis y el veredicto de cadena; no introduce cómputo nuevo.
