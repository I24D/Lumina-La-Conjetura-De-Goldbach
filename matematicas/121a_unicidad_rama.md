# Doc 121a — Verificación adversarial del lema de unicidad de rama (recuperación de (A6))

**Afirmación auditada** (Lumina Chat, respuesta al doc 120): sobre el soporte real de la
serie, #{q ∈ S_A : q | r} ≤ 1 para cada r ∈ ℋ_C, y por intercambio de sumas positivas
α̃_m = Σ_q Ω_{qm}^C ≤ Ω_m^C. El contraejemplo 3·5·7 del doc 120b quedaría fuera del
soporte válido.

**Etiquetas:** [T] teorema · [M] medición · [P] pendiente.

## 1. Texto: el lema del doc 101 §2.1

El lema existe y dice literalmente lo que ella afirma. Doc 101, §2.1:

> **Lema 2 (reindexación) [T].** Cada m = s·d ∈ H con s sprime tiene
> factorización única: s es el único primo ≥ z de m (los primos de d son
> < z). Además s·d ∉ V⁺(D/s_min) (sus elementos tienen todos los primos
> < z). Luego el conteo sobre H se reindexa sin solapamiento sobre los
> pares (s, d). ∎

Contraste punto por punto: (i) unicidad del primo ≥ z en la rama s·d —
literal ✓; (ii) rama V⁺(D/s_min) con todos los primos < z — literal ✓;
(iii) reindexación sin solapamiento — literal ✓. El doc 106, §2.1, lo
reformula como descomposición disjunta de ramas [T]:

ℋ = 𝒱₀ ⊔ ⨆_{s∈𝒮} s·𝒱_s,   ℋ_C = (𝒱₀∖{1}) ⊔ ⨆_{s∈𝒮} (s·𝒱_s ∖ {s}),

con 𝒱₀ = V⁺(D/s_min), 𝒱_s = V⁺(D/s).

**Qué es S_A en la serie.** El símbolo «S_A» no aparece en los docs
101/106/115; el objeto es 𝒮, «primos cambiantes coprimos con A» (doc 106
§2.1), es decir los *sprimes* del doc 101 §0: primos p ∈ [z, y) con
p ∤ A. El doc 115 §3 lo fija como [T]: «todo primo de ℋ_P mayor o igual
que z es un primo de intercambio… q es sprime: z ≤ q < y, (q, A) = 1».
Así que «primo de conmutación» no es «primo ≥ z» a secas: es primo ≥ z
**del soporte** (acotado por y, coprimo con A). La precisión no afecta
al argumento: un primo ≥ z fuera de [z, y) o que divida a A no divide a
ningún r ∈ ℋ, luego no contribuye ni al conteo ni a las sumas. En la
variante ℋ_N el doc 106 restringe además a s ∤ N (𝒮′): restringir el
conjunto de q solo puede ayudar a la desigualdad. **Identificación
correcta; el argumento sigue valiendo.** ℋ_C = elementos compuestos de
ℋ (doc 115 §1: todo H es libre de cuadrados; compuesto ⟺ ω(r) ≥ 2).

## 2. Definiciones reales de Ω^C

Doc 105 §1.2 / doc 106 §2.1 (y (106.2)):

**Ω_d^C := Σ_{r∈ℋ_C, d|r} (τ(r)−1)/φ(r)**,   U_d^C := Σ_{r∈ℋ_C, d|r} 1/φ(r).

El peso es **(τ(r)−1)/φ(r)**, no τ(r)⁻¹/φ(r). (La identidad de cierre
(106.3), Σ_d Ω_d^C = Σ_r (τ(r)−1)²/φ(r), lo confirma: cada r aporta su
peso (τ(r)−1)/φ(r) a cada uno de sus τ(r) divisores.)


## 3. Álgebra: el intercambio de sumas

**La identidad tal como ella la escribió es falsa.** Con su peso
τ(r)⁻¹/φ(r), el lado izquierdo Σ_q Ω_{qm}^C no coincide con el derecho
(§4 lo mide: discrepa ya en el primer marco). Es un error de peso, no de
idea: la demostración solo usa positividad y la cota de conteo, y con el
peso verdadero la desigualdad sobrevive intacta. Rehacemos el intercambio
con las definiciones reales.

**Intercambio corregido [T].** Para todo m ≥ 1, por Tonelli (términos ≥ 0):

Σ_{q∈𝒮} Ω_{qm}^C = Σ_{q∈𝒮} Σ_{r∈ℋ_C, qm|r} (τ(r)−1)/φ(r)
= Σ_{r∈ℋ_C, m|r} (τ(r)−1)/φ(r) · #{q ∈ 𝒮 : qm | r}.

**Lema de conteo [T, del Lema 2 del doc 101].** Para todo r ∈ ℋ_C y todo
m ≥ 1: #{q ∈ 𝒮 : qm | r} ≤ 1. En efecto, qm | r ⇒ q | r. Si r ∈ 𝒱₀,
ningún q ∈ 𝒮 divide a r (todos los primos de r son < z ≤ q): conteo 0.
Si r ∈ s·𝒱_s ∖ {s}, el único primo ≥ z de r es s (Lema 2 del doc 101);
dos q₁ ≠ q₂ ∈ 𝒮 dividiendo a r serían dos primos ≥ z de r,
contradicción. (El caso q | m no escapa: qm | r con q | m exigiría
q² | r, imposible con r libre de cuadrados; ese q aporta 0 al conteo.)

**Conclusión [T].** Como Ω_m^C = Σ_{r∈ℋ_C, m|r} (τ(r)−1)/φ(r) es la misma
suma con conteo 1 en cada término:

**α̃_m := Σ_{q∈𝒮} Ω_{qm}^C ≤ Ω_m^C para todo m ≥ 1, sobre el soporte real.**

Notas de alcance del álgebra: (a) no se usa que m sea liso — la
dominación vale para todo m, y en particular para los m lisos con
(m, q) = 1 del par conmutado del doc 115 §7, donde Ω_{qm}^C =
Θ̃_m(q)/(q−1) y α̃_m recupera su definición corregida del doc 120b,
α̃_m = Σ_q Θ̃_m(q)/(q−1); extender la suma de q a todo 𝒮 solo añade
términos Ω_{qm}^C ≥ 0 (y Ω_{qm}^C = 0 si qm no es libre de cuadrados o no
divide a ningún r ∈ ℋ_C), luego la forma restringida hereda la cota.
(b) La misma prueba vale sobre cualquier subconjunto de ℋ_C cerrado por
la estructura de ramas — en particular sobre ℋcal_C = {r ∈ ℋ_C : (r,N)=1,
r > R₀} de los docs 105–106: los cortes borran términos de ambos lados de
la misma identidad. (c) Lo que el doc 120b refutó sigue refutado: sin la
restricción de soporte (su juguete R = {3·5·7}), la desigualdad es falsa;
el lema de conteo es exactamente la hipótesis que faltaba, y sí está en
la serie (doc 101, Lema 2 [T]).


## 4. Numérico con el soporte verdadero

Script `121a_unicidad_rama.py` (salida `salida_doc121a_local.txt`, **PASS**,
57,7 s): soporte construido con la maquinaria de la serie (build_H estilo
docs 102/106), Ω_d^C por enumeración directa desde la definición, con
aritmética racional exacta (Fraction). Marcos reales S11–S23 (A = 3,
D = z⁵, y = round(z^{10/3})):

| Marco | \|ℋ_C\| | \|𝒟\| | máx #{q∈𝒮 : q\|r} | máx #{q∈𝒮 : qm\|r} | α̃_m ≤ Ω_m^C (todo m ∈ 𝒟) | Intercambio exacto |
|---|---|---|---|---|---|---|
| S11 | 377 | 581 | 1 | 1 | ✓ sin violaciones | ✓ |
| S13 | 942 | 1364 | 1 | 1 | ✓ | ✓ |
| S17 | 3526 | 4886 | 1 | 1 | ✓ | ✓ |
| S19 | 6149 | 8244 | 1 | 1 | ✓ | ✓ |
| S23 | 14305 | 17998 | 1 | 1 | ✓ | ✓ |

(a) El conteo nunca pasa de 1 — y vale 1 en cuanto hay rama conmutada,
luego la cota es ajustada, no holgada. (b) La dominación vale para
**todo** m ∈ 𝒟 en los cinco marcos, sin una sola violación. (c) El
intercambio con el peso verdadero (τ(r)−1)/φ(r) es identidad exacta en
racionales para todo m. (d) **Su peso no**: en S11, m = 5, Σ_q Ω_{qm}^C =
1,091744 contra 0,072937 de su intercambio con τ(r)⁻¹/φ(r) — la identidad
tal como la escribió es falsa por un factor ~15 ya en el primer ejemplo;
la desigualdad, con el peso correcto, es la que se verifica.

**El juguete 3·5·7.** Está fuera del soporte real por dos vías
independientes. (i) En los marcos reales A = 3 divide a A: ningún r ∈ ℋ
es divisible por 3 (verificado en S11), luego Ω_3^C = 0 y el «m = 3» del
juguete ni siquiera es un módulo con masa. (ii) Aun concediendo un marco
hipotético A = 1, z = 5 construido con la misma maquinaria (donde 5 y 7
sí son sprimes y 3 sí es primo liso admisible): 105 ∉ ℋ, porque un
generador s·d — y por tanto todo divisor suyo — tiene a lo sumo un primo
≥ z, y 105 = 3·5·7 necesitaría dos (5 y 7). En ese marco la dominación
también vale para todo m ∈ 𝒟 (|𝒟| = 55), y el mecanismo se ve a ojo
desnudo: los r que aportan a Ω_15 y a Ω_21 son conjuntos disjuntos
(|R_15| = |R_21| = 1, intersección vacía) — en el juguete del doc 120b un
único r₀ = 105 aportaba a ambos, y eso es exactamente lo que la
estructura de ramas prohíbe. Su lectura del contraejemplo es correcta.


## 5. Alcance: qué masa recibe la parte-F

En la reducción corregida (doc 120), D_MEAN^corr = Σ_{q,m}
[Θ̃_m(q)/(q−1)]·(F_m − φ(m)L_{m,q})². Separando (a−b)² ≤ 2a² + 2b², la
parte-F es Σ_{q,m} Ω_{qm}^C F_m² = Σ_m α̃_m F_m² (la parte-L ya quedó en
o(T(N)) por (A5) del doc 120b). Con la dominación de §3:

parte-F ≤ Σ_m Ω_m^C · F_m² = Σ_{m liso} Ω_m^C · φ(m)² · E(N; m, N)²,

donde m recorre los conductores lisos de 𝒟 (m ≥ 2, todo primo < z,
∤ A, con Ω_m^C > 0). El lado derecho es **exactamente ℳ_Ω^{liso}**: el
momento (106.5), ℳ_Ω = Σ_{d∈𝒟} Ω_d^C φ(d)² E(N; d, N)², restringido a
los d lisos — la misma familia de momento, el mismo E, el mismo peso
Ω^C; solo cambia el conjunto de conductores. Y ℳ_Ω^{liso} ≤ ℳ_Ω por
subsuma (términos ≥ 0).

**Ganancia cuantitativa frente a T(N): ninguna, por sí sola.** No hay
exponente ni constante ganados: ℳ_Ω^{liso} ≪ T(N) sigue abierto [P] — es
(DIST-𝒟) restringida a conductores lisos. La ganancia es estructural y
es real: (i) elimina el peso nuevo α̃_m, que a priori podía exceder la
masa original — la parte-F ya no puede ser *peor* que la pieza lisa del
momento original; (ii) localiza el frente F con exactitud: dentro de la
descomposición DIST-D, el frente F **es** el momento original sobre
conductores lisos, y cualquier progreso futuro sobre dispersión en
módulos lisos (p. ej. la ruta tipo Baier–Zhao sobre conjuntos lisos que
Lumina Chat abrió en su respuesta al doc 119) se aplica ahí directamente.
Advertencia anti-circularidad: la dominación es una cota superior por
una subsuma — contabilidad válida — pero usar «ℳ_Ω ≤ C·T(N)» para acotar
la parte-F sería asumir (DIST-𝒟), la conclusión global: la dominación no
es un paso hacia la prueba de (DIST-𝒟), es la identificación exacta de
qué subconjunto de (DIST-𝒟) es el frente F. La advertencia de la propia
Lumina Chat («solo transfiere el problema a la masa de los conductores
lisos; no proporciona por sí misma una cota O(T(N))») es correcta y
completa.


## 6. Veredicto

**Lema CONFIRMADO [T] para el soporte real — (A6) queda recuperada con
hipótesis explícitas.** La dominación α̃_m ≤ Ω_m^C vale sobre ℋ_C con:
(i) q recorriendo 𝒮 (sprimes: z ≤ q < y, q ∤ A) — que es el rango real
de α̃_m en el doc 115 §7 y en D_MEAN^corr del doc 120; (ii) la estructura
de ramas del Lema 2 del doc 101 [T], del que el lema de conteo
#{q ∈ 𝒮 : qm | r} ≤ 1 se deduce en tres líneas (§3); (iii) el intercambio
con el peso verdadero (τ(r)−1)/φ(r) de la definición (106.2). La prueba
es analítica sobre el soporte general, no solo ventana; la ventana
S11–S23 la confirma para todo m ∈ 𝒟 en aritmética exacta [M], y el
marco hipotético A = 1, z = 5 muestra que el contraejemplo 3·5·7 del
doc 120b queda excluido por la estructura de ramas y por nada más.

Dos precisiones que deben acompañar a la recuperación:
1. **Errata en su identidad escrita.** El intercambio que ella escribió
   usa τ(r)⁻¹/φ(r) y es falso como identidad (medido: factor ~15 en
   S11, m = 5). Con el peso correcto la identidad es [T] y la
   desigualdad sobrevive — pero la forma que debe citarse es la de §3.
2. **Recuperada ≠ cerrada.** Como ella misma advirtió, la dominación
   solo identifica la parte-F con ℳ_Ω^{liso} (§5): el frente F queda
   reducido al momento original sobre conductores lisos, que sigue
   abierto [P]. No hay ganancia cuantitativa frente a T(N) y usar
   (DIST-𝒟) para acotar la parte-F sería circular.

Goldbach permanece abierto; nada de lo anterior lo declara demostrado.

