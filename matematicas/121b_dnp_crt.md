# Doc 121b — Verificación adversarial de la reducción CRT de D_NP (respuesta de Lumina Chat al doc 120)

**Objeto:** la cadena (ortogonalidad mód m) → (Cauchy–Schwarz + CRT) → (forma final con peso Ω_d^C/φ(d)) con que Lumina Chat propone convertir D_NP en un segundo momento de caracteres Ψ_N(η) sobre la rama conmutada 𝒟_S.
**Definiciones de la serie usadas:** doc 115 §7 [T]: Ω_{qm}^C = Θ̃_m(q)/(q−1), Θ̃_m(q) = Σ_{h∈V⁺(D_A/q)\{1}, m|h} (2τ(h)−1)/φ(h); doc 120b §6: D_NP = Σ_{q,m} Θ̃_m(q)·φ(m)²·Σ_{χ≠χ₀}|S_{m,q}(χ)|², S_{m,q}(χ) = Σ_{n≤N, n≡N (mód m), q∤n} Λ(n)χ(n).
**Etiquetas:** [T] teorema · [M] medición · [P] pendiente. Script: `121b_dnp_crt.py` (salida `salida_doc121b_local.txt`).

## 1. (1) Ortogonalidad: la fórmula de S_{m,q}(χ) — **CONFIRMADA** [T][M]

**Derivación [T].** Para (n, m) = 1 y (N, m) = 1, el indicador de la clase es
1[n ≡ N (mód m)] = (1/φ(m))·Σ_{ψ mód m} ψ̄(N)ψ(n) — con ψ̄(N), no ψ(N):
la ortogonalidad da ψ(n)ψ̄(N) porque el carácter debe evaluarse en n·N⁻¹.
La restricción a (n, m) = 1 no pierde nada de S: si n ≡ N (mód m) y
(N, m) = 1, todo primo común de n y m dividiría a N, luego (n, m) = 1.
Sustituyendo y reagrupando:
S_{m,q}(χ) = (1/φ(m))·Σ_ψ ψ̄(N)·Σ_{n≤N} Λ(n)ψ(n)χ(n).
El producto ψχ es el carácter mód qm vía CRT ((m, q) = 1), y la condición
q∤n es **automática para todo χ mód q**: χ(n) = 0 si q | n (también para
χ₀, que vale [q∤n]). Luego la suma interna es exactamente Ψ_N(ψχ), sin
restricción extra ni término corrector. La fórmula vale tal como ella la
escribe, incluido el conjugado.

**Numérico [M, K1].** Identidad exacta (error ≤ 2,7·10⁻¹⁴) en los 5 casos
(N, m, q) ∈ {(200,3,7), (200,21,13), (300,7,11), (100,11,7), (300,13,17)},
para todo χ ≠ χ₀. Control del conjugado: en 4 casos ψ(N) ∈ {±1, 1} es real
y no discrimina; el caso discriminante (200, 21, 13) confirma ψ̄(N) y
**refuta la variante sin conjugar** (error 1,74·10¹ contra 7,5·10⁻¹⁵).

## 2. (2) Cauchy–Schwarz y conteo CRT — **CONFIRMADA** [T][M]

**Conteo exacto [T].** Por CRT, η mód qm ↔ (η_m, η_q) = (ψ, χ). Los pares
con χ ≠ χ₀ son φ(m)·(q−2); los η con η_q ≠ χ₀ son φ(qm) − φ(m)·1 =
φ(m)(q−1) − φ(m) = φ(m)(q−2). **Coinciden**: la suma doble Σ_{χ≠χ₀}Σ_ψ
|Ψ_N(ψχ)|² ES la suma Σ_{η: η_q≠χ₀} |Ψ_N(η)|², sin perder ni duplicar
caracteres (la condición η_q ≠ χ₀ excluye exactamente la componente
principal en q, ni más ni menos). El paso CS es por χ fijo:
|S_{m,q}(χ)|² ≤ (1/φ(m)²)·φ(m)·Σ_ψ |Ψ_N(ψχ)|², usando |ψ̄(N)| = 1
(otra vez (N, m) = 1). Sumando en χ se obtiene su desigualdad tal cual.

**La pérdida [M, K2].** La reindexación no pierde nada; toda la pérdida es
Cauchy–Schwarz (igualdad ⟺ Ψ_N(ψχ) ∝ ψ̄(N) en ψ, genéricamente falso).
Factor agregado medido λ = [(1/φ(m))·suma doble] / Σ_{χ≠χ₀}|S|²:
**2,40** (m=3), **4,06** (m=7), **11,64** (m=21), **14,52** (m=13).
λ no está acotado por φ(m) en general (una cancelación fuerte en
Σ_ψ ψ̄(N)Ψ_N(ψχ) lo dispara); en ventana es de un dígito o dos. Este λ
es el precio real de la reducción y vuelve en §3.

## 3. (3) Pesos: la identidad es puntual, pero la forma final está **REFUTADA tal como está escrita** [T][M]

**La identidad que necesita sí es término a término [T].** El doc 115 §7
demuestra Ω_{qm}^C = Θ̃_m(q)/(q−1) para cada (q, m) (verificado allí por
enumeración, J4). No hay problema de «conversión solo sumada»: la relación
puntual existe y es la de la serie. La auditoría (A6) del doc 120b solo
había establecido la suma porque era lo que la afirmación auditada usaba.

**Pero la conversión da otro peso [T].** El término de D_NP, tras CS, es
Θ̃_m(q)·φ(m)²·(1/φ(m))·Σ_{η} |Ψ_N(η)|² = Θ̃_m(q)·φ(m)·Σ_η |Ψ_N(η)|².
Y Θ̃_m(q)·φ(m) = (q−1)·Ω_{qm}^C·φ(m) = **Ω_d^C·φ(d)**, con d = qm,
φ(d) = (q−1)φ(m). Su forma final lleva **Ω_d^C/φ(d)**: un factor **φ(d)²
más pequeña**. Es el mismo patrón de inversión de φ ya corregido en los
docs 115 (dos veces) y 120b (A4/A6): cuarta aparición en la línea.

**Contraejemplos término a término [M, K3].** Marco genuino de la serie
(z = 17, A = 3, D_A = z⁵, Θ̃ calculado con la maquinaria V⁺ de los docs
109/110, N = 200): en las 9 parejas (q, m) discriminantes
(q ∈ {17, 19, 23}, m ∈ {7, 11, 13}, Θ̃_m(q) > 0, Σ|S|² > 0) el término
real de D_NP supera su forma final por un factor **1622×–4768×**, que es
exactamente φ(d)²/λ con el λ medido (5,3–14,6) frente a φ(d)² =
9216–69696. Su forma exigiría λ ≥ φ(d)² en cada término: falso por dos o
tres órdenes de magnitud. (Las filas m = 77 tienen S ≡ 0 a N = 200 y no
discriminan.) En las 12 filas el término sí cumple la cota corregida.

**Forma correcta de la reducción:**
D_NP ≤ Σ_{d=qm ∈ 𝒟_S} **Ω_d^C·φ(d)** · Σ_{η mód d, η_q ≠ χ₀} |Ψ_N(η)|².
La versión con Ω_d^C/φ(d) no se sigue de su derivación y, como desigualdad
término a término, es falsa; sin un argumento nuevo (no hay mecanismo de
cancelación entre términos: todos son ≥ 0 y la derivación es por (q, m)),
la versión sumada queda **no verificada y sin crédito**.

## 4. (4) La hipótesis (m, N) = 1: garantizada por construcción, y necesaria [T][M]

**En la serie está garantizada [T].** El doc 101 define el soporte con el
filtro explícito (r, N) = 1 (su H_N y su conjunto de divisores), y el doc
115 §1 lo fija como propiedad de ℋ: «ℋ exige (r, N) = 1»; los primos de
intercambio entran con (s, AN) = 1 y los pasos de membresía de la parte
lisa usan p ∤ AN (doc 115 §§2–3). Todo m que aparece en D_NP divide a un
h ∈ V⁺(D_A/q), y h divide a su vez a elementos de ℋ; luego (m, N) = 1
para **cada** m con Θ̃_m(q) > 0. Los m con (m, N) > 1 llevan masa
exactamente cero en D_NP: no hay que tratarlos porque nunca ocurren.
Precisión de implementación: los scripts modelo (`allowed_primes`) solo
excluyen p | A; la exclusión de los primos de N vive en las hipótesis de
la capa analítica ((s, AN) = 1, filtro H_N) — al aplicar la fórmula hay
que imponer ese filtro antes, no después.

**Y es necesaria, no cosmética [M, K4].** Con (N, m, q) = (200, 35, 11),
(m, N) = 5: ψ(200) = 0 para todo ψ mód 35, la fórmula da 0 idénticamente,
pero S_{m,q}(χ) = log 5 · χ(3) ≠ 0 (la clase n ≡ 200 mód 35 contiene
n = 25, coprimo con q pero no con m). Sin (m, N) = 1 la ortogonalidad no
detecta la clase y la reducción entera colapsa en silencio.

## 5. Veredictos y balance

| Pieza | Veredicto |
|---|---|
| (1) Fórmula de ortogonalidad para S_{m,q}(χ), con ψ̄(N) | **CONFIRMADA** [T][M] |
| (2) CS + biyección CRT, conteo φ(m)(q−2) | **CONFIRMADA** [T][M]; pérdida λ medida 2,4–14,5 |
| (3) Identidad Θ̃_m(q) = (q−1)Ω_{qm}^C | **CONFIRMADA** puntual [T, doc 115 §7] |
| (3) Forma final con peso Ω_d^C/φ(d) | **REFUTADA** tal como está [T][M]: el peso correcto es Ω_d^C·φ(d) (factor φ(d)²) |
| (4) Hipótesis (m, N) = 1 | **CONFIRMADA**: garantizada por construcción en la serie [T] y necesaria [M] |
| Reducción completa | **CONFIRMADA CON ENMIENDA material**: vale con el peso Ω_d^C·φ(d) |

**El problema que queda, exactamente:**
Σ_{d=qm ∈ 𝒟_S} Ω_d^C·φ(d)·Σ_{η mód d, η_q ≠ χ₀} |Ψ_N(η)|² ≪ T(N),
un segundo momento de caracteres sin restricción de clase sobre el
conjunto estructurado 𝒟_S de la rama conmutada.

**¿Más tratable que el D_NP original, o solo reformulación?** Ganancia
estructural real: desaparece la restricción n ≡ N (mód m) y la sustracción
de Parseval queda absorbida en la definición; el objeto resultante tiene
la forma estándar que atacan la criba grande y los BV para conjuntos
dispersos de módulos (la ruta Baier–Zhao que ella misma propuso en su
respuesta al doc 119 encaja aquí, ahora con el peso correcto). Pero dos
costes honestos: (i) la pérdida de Cauchy–Schwarz λ (medida 2,4–14,5,
sin cota general) ya está pagada y es irrecuperable dentro de esta ruta;
(ii) con el peso corregido el objetivo es φ(d)² **más pesado** que el que
ella anunció — la reducción no abarata el problema, lo reescribe en la
forma en la que la estructura de 𝒟_S puede explotarse. Sigue siendo
necesaria una ganancia de dispersión genuina sobre 𝒟_S: la reducción,
corregida, es el puente correcto hacia ese frente, no su resolución.

Goldbach permanece abierto; nada de lo anterior lo declara demostrado.
