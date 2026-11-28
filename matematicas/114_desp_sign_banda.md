# Doc 114 — (DESP-SIGN) por apareamiento y la banda degenerada por conteo

**Fecha:** 2026-11-28. **Línea:** UNIF-U / (H2), continuación del doc 113.
**Script:** `114_desp_sign_banda.py` (stdlib, importa las definiciones de los docs 109–113).
**Salida:** `salida_doc114_local.txt` — chequeos K0–K7, **todos PASS** (4,1 s).
**Etiquetas:** [T] = teorema (prueba completa en el texto), [M] = medición por enumeración exacta en ventana, [P] = pendiente.

---

## 0. Veredicto en breve

El doc 113 dejó (FL-DESP⁺) reducida a la FL clásica sobre el modelo más dos residuos:
**(DESP-SIGN)** — el lema de signo E(r) ≤ 0 — y la **banda** u_r < u₀. Este doc ataca ambos.

1. **(DESP-SIGN) queda reducido exactamente [T]** a dos enunciados residuales pequeños y nombrados:
   **(R-FLOW-71+)** — la cota de flujo R_A(y) ≤ 1 (y R_R(y) ≤ 1) solo para y con min(y) ≥ 71 —
   y **(REM-PERS-small)** — la persistencia del fallo fusionado cuando *todos* los fallos de r·x
   caen en primos ≤ q_min(x). Fuera de esos residuos, todo es teorema: dos lemas de herencia
   (§2, §3) y las identidades de flujo (§3, §4). Ambos residuos están **vacíos en las 4298
   instancias medidas** (4113 de los marcos del doc 113 + 185 de estrés con z = 53), y
   (R-FLOW) está **probado [T] en todo marco con z_A ≤ 71** por una suma de Mertens exacta
   (Σ_{5≤p≤67} 1/(p−1) = 0,99917 < 1). (DESP-SIGN) **no** se declara demostrado en el marco
   del puente (z_A → ∞ con N): ahí los dos residuos siguen vivos.
2. **La banda queda reducida exactamente** a **(BAND)**: la sumabilidad de la cota cruda
   (1−Π(r)²)/φ(r) sobre la sub-banda {Π(r) < 1/2}. La sub-banda {Π(r) ≥ c₀} está **cerrada [T]**
   (≤ Q_completo/c₀²). En ventana la banda entera suma ≤ 0,6991·Q_c con la cota cruda a
   u₀ = 1,5 [M] y su masa real es ≤ 0,3806·Q_c [M]; en u_r < 1 la cota cruda es *exacta*
   (δ > 0 ⟹ S(r) = 1, degenerado) [M, 349 instancias].
3. **(H2) NO queda demostrada.** El inventario exacto de lo que falta está en §8. Si
   (R-FLOW-71+), (REM-PERS-small) y (BAND) cierran, (H2) se sigue por el álgebra del doc 113
   con constante ≈ 6,4 para (H2⁗) más la constante de banda.

Goldbach no se declara demostrado en ningún punto de este documento.

---

## 1. Estructura de E(r): la diferencia simétrica, caracterizada

Con las definiciones del doc 113: V̂^{(r)} es el soporte modelo (V⁺(L/r; 𝒫∖r)) y
V^{(r)} = {d′ : (d′,r)=1, r·d′ ∈ V⁺(L)} el soporte desplazado real. Sean

- **Add(r)** = V^{(r)} ∖ V̂^{(r)} (entran al imponer la condición sobre el producto),
- **Rem(r)** = V̂^{(r)} ∖ V^{(r)} (salen).

**Identidad [T]:** E(r) = S(r) − Ŝ(r) = Σ_{d′∈Add} a(d′) − Σ_{d′∈Rem} a(d′), con
a(d′) = μ(d′)/φ(d′). Verificación por enumeración (K0/K1): error máximo 5,3·10⁻¹⁶ sobre
las 4113 instancias; |Add| total = 6509, |Rem| total = 1049 (los conteos del doc 113, K6);
E(r) ≤ 0 en las 4113 (máx E/Π = 0,0 a precisión de máquina).

### 1.1 Caracterización de Add (K1)

Todo d′ ∈ Add falla el modelo por al menos una posición impar m (condición de Nathanson
r·q_m²·P_m ≥ L), y el fallo es de exactamente dos tipos:

- **A1 — evaporación de paridad** (2527 elementos): la posición fusionada de q_m en r·d′,
  j = m + #{p | r, p > q_m}, es **par**: la condición desaparece al fusionar.
- **A2 — holgura** (3958 elementos): j sigue impar, pero la condición fusionada es más débil:
  su lado izquierdo es el del modelo dividido por ∏_{p|r, p<q_m} p.
- **Mixtos** (24 elementos): ambos tipos a la vez. Ningún elemento de Add carece de fallo del modelo.

**Lema de holgura [T]:** si q_m falla en posición fusionada impar (tipo A2), entonces r tiene
un primo < q_m. *Prueba:* si ∏_{p|r, p<q_m} p = 1, el lado izquierdo fusionado en q_m coincide
con el del modelo (los primos de r mayores que q_m entran igual en ambos productos parciales),
luego sería ≥ L, contra d′ ∈ V^{(r)}. ∎ Verificado: 0 contraejemplos en ventana (K1).

### 1.2 Caracterización de Rem (K1)

Todo d′ ∈ Rem pasa el modelo pero r·d′ falla la condición fusionada en alguna posición impar:

- **R1 — flip de paridad** (707 elementos): el primo que falla es un primo propio q de d′.
- **R2** (342 elementos): el primo que falla es un primo de r. Sin mixtos; ninguno sin clasificar.

**Teorema [T]:** todo primo R1 falla en posición propia **par** de d′. *Prueba:* si q = q_m
con m impar fallara fusionado, su lado izquierdo fusionado q²·(∏_{p|r,p>q} p)·P_m ≤ r·q²·P_m,
y el fallo fusionado da L ≤ r·q²·P_m: el modelo fallaría en m, contra d′ ∈ V̂^{(r)}. ∎
Verificado: 707/707 en posición propia par (K1). Es el espejo exacto de A1: los elementos de
Rem⁻ peligrosos (μ = −1, ω impar... véase §4) pierden su condición por paridad.

---

## 2. Lema 2 [T]: V⁺(L) es cerrado por divisores

**Enunciado:** si d ∈ V⁺(L) y p | d es primo, entonces d/p ∈ V⁺(L). Por inducción,
e | d ⟹ e ∈ V⁺(L).

**Prueba.** Sean p_1 > … > p_k los primos de d y p = p_j. Tamaño: d/p < d < L, y el primo
mayor de d/p es ≤ p_1 < L. Condiciones de Nathanson de d/p, por primo p_i:

- *i < j* (p_i mayor que p): su posición y su producto parcial en d/p son los de d
  (p no entra en ninguno de los dos). La condición es idéntica a la de d en i. ✓
- *i > j* (p_i menor que p): su posición pasa de i a i−1 y su producto parcial de P_i a
  P_i/p. Si i−1 es par no hay condición. Si i−1 es impar (i par), hay que ver
  p_i²·P_i/p < L. La condición de d en la posición impar i−1 da
  p_{i−1}²·P_{i−1} < L con P_{i−1} = P_i/p_i, es decir P_i < L·p_i/p_{i−1}². Luego
  p_i²·P_i/p < L · p_i³/(p·p_{i−1}²). Si i−1 = j (p_{i−1} = p): el cociente es
  p_i³/p³ < 1 pues p_i < p. Si i−1 > j: p_{i−1} < p y p_i < p_{i−1}, y el cociente es
  < p_{i−1}/p < 1. En ambos casos < L. ∎

**Corolarios [T]:** (a) r visible ⟹ r ∈ V⁺(L) (r divide a algún elemento de V) — licencia
retroactivamente el uso del doc 113 §4; (b) V^{(r)} es cerrado por divisores en d′
(e | d′, r·d′ ∈ V ⟹ r·e | r·d′ ⟹ r·e ∈ V); (c) 1 ∈ V^{(r)}.
Verificación [M] (K2): 0 violaciones sobre los 4353 elementos de V de todos los marcos,
probando d/p para *cada* primo p | d.

---

## 3. El apareamiento en Add: Lema 1, flujo y el umbral de Mertens

El apareamiento es **φ(x) = x / q_min(x)**: quitar el primo *menor* (no el mayor: dividir
por el primo frontera manda el socio a Both en 1626/1640 casos [M, exploración previa]).

**Lema 1 [T]:** si ω(x) es par y x ∉ V̂^{(r)}, entonces x/q_min(x) ∉ V̂^{(r)}.
*Prueba:* las condiciones del modelo de x/q_min son las de x en las posiciones
1, …, ω(x)−1, con idénticos productos parciales (el primo menor ocupa la última posición,
ω(x), par: no lleva condición). Si x falla el modelo, falla en alguna posición impar
m ≤ ω(x)−1 (la posición ω(x) es par y no puede fallar; el chequeo de tamaño se hereda
pues x/q_min < x). Esa condición es también condición de x/q_min. ∎
Además, el paso al modelo es hereditario sin condiciones: x ∈ V̂^{(r)} ⟹ x/q ∈ V̂^{(r)}
para todo primo q | x (mismo argumento de posiciones). Y por el Lema 2,
x ∈ V^{(r)} ⟹ x/q_min ∈ V^{(r)}.

**Consecuencia [T]:** φ manda Add⁺ = {x ∈ Add : μ(x) = +1, x > 1} dentro de Add⁻ =
{y ∈ Add : μ(y) = −1}: el socio conserva el fallo del modelo (Lema 1, ω(x) par porque
μ(x) = +1) y conserva la pertenencia a V^{(r)} (Lema 2). Verificado [M] (K3): 1640/1640.

**Identidad de flujo [T]:** agrupando por la imagen y ∈ Add⁻, con
R_A(y) = Σ_{q : yq ∈ Add⁺} 1/(q−1) (la razón exacta de términos apareados es
a(yq)/|a(y)| = 1/(q−1), pues φ(yq) = φ(y)(q−1)):

$$E_{Add}(r) = -\sum_{y \in Add^-} |a(y)|\,\bigl(1 - R_A(y)\bigr).$$

Verificada por instancia con error ≤ 5,6·10⁻¹⁷ (K3). Luego **E_Add(r) ≤ 0 ⟸ R_A(y) ≤ 1**
para todo y. La masa peligrosa total de Add⁺ es 18,41 contra 93,29 de masa socia en Add⁻:
el flujo tiene holgura estructural, no es un empate.

**Cota de Mertens [T]:** R_A(y) ≤ Σ_{q ∈ 𝒫, q < min(y)} 1/(q−1) ≤ Σ_{5≤p<min(y)} 1/(p−1)
(en la línea, 3 | A luego 𝒫 ⊆ primos ≥ 5). Con fracciones exactas (K3):

$$\sum_{5 \le p \le 67} \frac{1}{p-1} = 0{,}99917 < 1, \qquad
\text{añadiendo } p=71:\ 1{,}01346 > 1.$$

Por tanto **R_A(y) ≤ 1 es [T] si min(y) ≤ 67, y en todo marco con z_A ≤ 71**. Medición:
máx R_A = 0,6000 en los marcos del doc 113 (en z=29, L=2261685, r=323, y=23:
1/4+1/6+1/10+1/12); en el estrés z=53, máx R_A = 0,6625 = suma parcial de Mertens hasta 17
exactamente (K7). La cota cruda deja de bastar solo cuando min(y) ≥ 71, lo que exige
z_A ≥ 73: cualquier prueba en ese régimen debe usar la admisibilidad (r·y·q ∈ V), no
solo q < min(y).

---

## 4. El lado Rem: herencia, persistencia y su residuo

Apareamiento ψ(x) = x/q_min(x) sobre Rem⁻ = {x ∈ Rem : μ(x) = −1}.

- **Mitad modelo [T]:** ψ(x) ∈ V̂^{(r)} sin condiciones (herencia del modelo, §3).
- **Persistencia [T, parcial]:** si r·x falla la condición fusionada en un primo
  p* > q_min(x), entonces r·ψ(x) ∉ V: al quitar q < p*, la posición fusionada de p* y su
  producto parcial (primos ≥ p*) no cambian, luego la misma condición sigue fallando.
- **Residuo (REM-PERS-small):** x ∈ Rem⁻ cuyos fallos fusionados caen *todos* en primos
  ≤ q_min(x) — incluye el caso ω(x) = 1 (x primo: ψ(x) = 1 ∈ Both, descubierto por
  construcción). Para esos x no hay prueba de que ψ(x) ∈ Rem.

**Medición (K4):** 39/39 socios caen en Rem en los marcos del doc 113; **ningún** x ∈ Rem⁻
tiene ω(x) = 1; **todos** tienen un fallo fusionado en un primo > q_min(x) (el mecanismo [T]
cubre el 100% de los casos medidos); en el estrés z=53, 36/36 (K7). La familia residual
está **vacía en ventana**. La identidad de flujo análoga vale [T dado el aterrizaje]:

$$E_{Rem}(r) = -\sum_{y \in Rem^+} a(y)\,\bigl(1 - R_R(y)\bigr),$$

error ≤ 3,5·10⁻¹⁸ (K4); máx R_R = 0,2500 en ventana (0,4167 en el estrés z=53), con la
misma cota de Mertens de §3. La masa de Rem⁻ es minúscula: 0,0288 en total, contra 18,41
la de Add⁺.

Medición conjunta (K0, K7): E_Add(r) ≤ 0 y E_Rem(r) ≤ 0 **por separado** en las 4298
instancias; el peor E observado es +1,1·10⁻¹⁶ (ruido de coma flotante, z=53).

---

## 5. Estado exacto de (DESP-SIGN)

**(DESP-SIGN) queda reducido [T] a:**

- **(R-FLOW-71+):** R_A(y) ≤ 1 y R_R(y) ≤ 1 para y con min(y) ≥ 71 (solo posible si
  z_A ≥ 73). Probado [T] fuera de ese régimen (Mertens, §3). La cota cruda no basta ahí:
  Σ_{5≤p≤71} 1/(p−1) > 1; hace falta explotar la admisibilidad r·y·q ∈ V⁺(L).
- **(REM-PERS-small):** ningún x ∈ Rem⁻ tiene todos sus fallos fusionados en primos
  ≤ q_min(x), y ningún primo x ∈ Rem⁻. Vacío en las 4298 instancias medidas.

En todo marco con z_A ≤ 71, (DESP-SIGN) es [T] módulo (REM-PERS-small) únicamente.
En el marco del puente (H2) — z_A = z_A(N) → ∞ — ambos residuos son enunciados vivos,
aunque estrictamente más pequeños que el original: (R-FLOW-71+) es una desigualdad tipo
Mertens sobre conjuntos de primos admisibles, y (REM-PERS-small) es un enunciado de
posición de fallos. **(DESP-SIGN) no se declara demostrado.**

---

## 6. Corolario [T/M]: S(r) ≤ 1, luego δ(r) ≤ 1 − Π(r)

El mismo flujo cierra la cota del soporte completo: V^{(r)} es cerrado por divisores
(Lema 2), contiene a 1, y agrupando V^{(r)+}∖{1} sobre V^{(r)−} vía x ↦ x/q_min(x):

$$S(r) = 1 - \sum_{y \in V^{(r)-}} |a(y)|\,\bigl(1 - R_{V^{(r)}}(y)\bigr) \le 1
\quad\text{si } R_{V^{(r)}}(y) \le 1,$$

con la misma cota de Mertens: **S(r) ≤ 1 es [T] si z_A ≤ 71** (y en general módulo
(R-FLOW)). Medición (K5): máx S(r) = 1,000000 sobre las 4113 instancias, y
δ(r) ≤ 1 − Π(r) sin violaciones. Esto legitima la cota cruda de la banda:
t(r)⁺ ≤ (2Π(r)+δ(r))·δ(r)/φ(r) ≤ (1−Π(r)²)/φ(r).

---

## 7. La banda u_r < u₀ por conteo directo (K6)

Cota cruda sumada sobre la banda, por nivel (máximos sobre los 240 niveles):

| u₀ | máx_L Σ_banda (1−Π²)/φ / Q_c | de la cual Π(r) < 1/2 | máx_L masa real Σt⁺/Q_c | máx #r en banda/nivel |
|---|---|---|---|---|
| 1,5 | 0,6991 | 0,4408 | 0,3806 | 77 |
| 2,0 | 1,1345 | 0,5971 | 0,3806 | 136 |

- **Sub-banda {Π(r) ≥ c₀}: CERRADA [T].** Ahí (1−Π²)/φ ≤ (Π²/φ)/Π² ≤ (Π²/φ)/c₀², y
  Σ_banda Π(r)²/φ(r) ≤ Q_completo por definición del marco (doc 110). Cota: Q_c/c₀²
  (con c₀ = 1/2: 4·Q_c; la medición da mucho menos).
- **Sub-banda u_r < 1: la cota cruda es exacta [M].** En las 349 instancias con u_r < 1 y
  δ(r) > 0, S(r) = 1 sin excepción: soporte degenerado V^{(r)} = {1} en efecto,
  δ(r) = 1 − Π(r), y la "cota" es la identidad.
- **Residuo (BAND):** Σ_{r visible, u_r<u₀, Π(r)<1/2} (1−Π(r)²)/φ(r) ≤ C·Q_completo
  uniformemente en el marco. En ventana vale con C = 0,60 [M, u₀ = 2]. Las rutas de
  conteo disponibles no lo dan en el marco: acotar 1/φ vía Π_full pierde (ln z)², y
  contar por visibilidad + cofactor pierde factores de escala |V|. (BAND) es estrictamente
  más débil que (FL-DESP⁺) en la banda (solo pide sumabilidad), pero sigue **[P]**.

**La banda no se declara cerrada [T]** en el marco; está cerrada en ventana [M] y
reducida al enunciado (BAND) sobre la sub-banda Π(r) < 1/2.

---

## 8. Veredicto de cadena (tarea 4)

**(H2) NO queda demostrada por este doc.** Lo que el doc 114 cambia es el tamaño y la
forma de lo que falta. Inventario exacto, con la cadena de dependencias:

- Doc 109: (H2) ⟺ forma corregida sup Q_s ≪ Q_completo ⟸ (H2′), cota L² de defectos.
- Doc 110: (H2) ⟺ (H2′) salvo constantes ±1; queda (H2″), exceso de Gram.
- Doc 111: (H2″) ⟺ (H2‴): masa positiva P^(L) ≤ C·Q_completo.
- Doc 112: bloque r = 1 cerrado [T] por la FL clásica; queda (H2⁗): masa positiva sobre
  conductores visibles r > 1.
- Doc 113: (H2⁗) ⟸ FL clásica sobre el modelo + (DESP-SIGN) + banda, con constante
  ≈ 6,4 si (FL-DESP⁺) vale con el ρ de la FL clásica.
- **Doc 114:** (DESP-SIGN) reducido [T] a **(R-FLOW-71+)** + **(REM-PERS-small)**
  (ambos vacíos en 4298 instancias; (R-FLOW) probado si z_A ≤ 71); banda reducida a
  **(BAND)** sobre {Π(r) < 1/2} (sub-banda Π ≥ c₀ cerrada [T]; ventana ≤ 0,70·Q_c).

**Teorema condicional actualizado:** (R-FLOW-71+) + (REM-PERS-small) + (BAND)
⟹ (DESP-SIGN) + banda ⟹ (FL-DESP⁺) ⟹ (H2⁗) ⟹ (H2), con la constante del doc 113
(≈ 6,4) más la constante de (BAND). Los tres residuos son enunciados de conteo/signo
explícitos, sin BV/BDH/criba/GRH, y ninguno tiene contraejemplo en ventana.

Siguen fuera de esta cadena, sin cambios: (DIST-𝒟⁺) (docs 107–108, la pieza de
distribución de L105(b)) y la separación de capas G0. Goldbach permanece abierto.

---

## 9. Reproducibilidad

`python3 114_desp_sign_banda.py | tee salida_doc114_local.txt` — stdlib únicamente;
importa `113_fl_desplazado_positivo.py` (que encadena 112→111→110→109). Chequeos:
K0 cruce con doc 113 (6509/1049, E ≤ 0); K1 caracterización A1/A2/R1/R2 con los dos
lemas de posición [T] verificados sin excepción; K2 Lema 2 por enumeración (0/4353·ω);
K3 flujo Add exacto + Mertens con fracciones; K4 flujo Rem + persistencia medida 100%;
K5 S(r) ≤ 1; K6 banda (u₀ = 1,5 y 2, splits, degeneración u_r < 1); K7 estrés z = 53
(185 instancias fuera de los marcos del doc 113). **PASS global.**
