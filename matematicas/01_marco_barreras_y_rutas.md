# 01 — Marco, resultados conocidos, dónde se rompe cada ataque y rutas posibles

**Procedencia:** Claude Code (`claude-opus-5-5`), 2026-11-23. Documento de trabajo, estado PROPUESTA.
**Fuente del contenido:** conocimiento matemático estándar del modelo. Las referencias marcadas *(verificar)* deben comprobarse en la fuente original antes de citarlas.

## 0. Etiquetas de estado (obligatorias en todo el proyecto)

| Etiqueta | Significado |
|---|---|
| **[T]** | Teorema publicado y aceptado por la comunidad |
| **[T-Lean]** | Verificado formalmente en Lean o Mathlib |
| **[C: H]** | Demostrado suponiendo la hipótesis H |
| **[N]** | Evidencia numérica. **Nunca** es demostración |
| **[H]** | Heurística o conjetura |
| **[P]** | Propuesta del proyecto, sin revisar |
| **[P-rev]** | Propuesta revisada por otro agente. Ojo: no es revisión independiente si los agentes comparten usuario de Windows |
| **[FALSO]** | Refutado; se conserva en el registro para no repetirlo |

---

## 1. Enunciado y reducciones elementales

**Conjetura de Goldbach (fuerte o binaria).** Todo entero par n ≥ 4 es suma de dos primos.

Notación: r(n) = #{(p, q) primos : p + q = n}, contando pares ordenados, y la versión con pesos

  R(n) = Σ_{a+b=n} Λ(a) Λ(b),  donde Λ es la función de von Mangoldt.

**Reducción 1.1 [T, elemental].** Sea R*(n) = Σ_{p+q=n} log p · log q, solo con primos. En R(n) − R*(n) quedan los términos en que a o b es potencia p^k con k ≥ 2. Cada uno de esos términos vale como mucho Λ(a)·log n, y Σ_{a = p^k ≤ n, k≥2} Λ(a) = ψ(n) − θ(n). Por tanto

  0 ≤ R(n) − R*(n) ≤ 2 (ψ(n) − θ(n)) log n ≤ 2,9 · √n · log n,

usando ψ(x) − θ(x) < 1,43 √x (Rosser–Schoenfeld; *verificar constante*).
**En consecuencia:** basta probar que **R(n) > 2,9 √n log n** para todo n par mayor que el límite verificado por ordenador. Por debajo de 4·10¹⁸ la conjetura ya está verificada **[N]**: Oliveira e Silva, Herzog y Pardi (2014).

**Conjetura de Hardy–Littlewood (1923) [H].**
  R(n) ∼ 𝔖(n) · n,  con 𝔖(n) = 2C₂ ∏_{p | n, p > 2} (p − 1)/(p − 2),  C₂ = ∏_{p>2} (1 − 1/(p−1)²) ≈ 0,6601618.
Para n par, 𝔖(n) ≥ 2C₂ ≈ 1,32. **Lo que hay que probar es "mucho menos" que esto:** basta con que R(n) sea mayor que una cantidad de tamaño √n log n, frente al valor esperado ≈ 1,32 n. La dificultad no está en la precisión. Está en conseguir una cota inferior **para cada n**, no en promedio.

---

## 2. El método del círculo: dónde se rompe exactamente

Sea S(α) = Σ_{m ≤ n} Λ(m) e(mα), con e(x) = e^{2πix}. Entonces

  R(n) = ∫₀¹ S(α)² e(−nα) dα = ∫_𝔐 + ∫_𝔪   (arcos mayores 𝔐 y arcos menores 𝔪).

- **Arcos mayores** (α cerca de a/q con q pequeño): dan el término principal 𝔖(n)·n·(1 + o(1)). Esto se controla con los ceros de las funciones L de Dirichlet [T], y mejor aún con GRH.
- **Arcos menores, el cuello de botella.** La mejor cota trivial es Parseval:
  |∫_𝔪 S² e(−nα)| ≤ ∫₀¹ |S|² = Σ_{m≤n} Λ(m)² = n log n − n + o(n)  [T].
  **El término de error es log n veces más grande que el principal (≈ 1,32 n).** La pérdida es solo un factor logarítmico, pero no se arregla mejorando las cotas de |S(α)|: la masa L² ya es demasiado grande. Haría falta **cancelación de fase** entre S(α)² y e(−nα) para un n concreto.
- **Por qué el caso ternario sí funciona (Vinogradov 1937, Helfgott 2013) [T].** En ese caso ∫_𝔪 |S|³ ≤ sup_𝔪 |S| · ∫|S|², y sup_𝔪 |S| ≪ n (log n)^{−A}. Un factor de |S| "sobra" para ganar. En el caso binario no sobra ninguno.
- **Circularidad.** Probar |∫_𝔪 S² e(−nα)| < (1 − δ)·𝔖(n)·n para cada n *es* esencialmente el problema binario reformulado. Además, **los arcos menores no están gobernados por los ceros de las funciones L de Dirichlet**, así que GRH no los controla punto a punto.

### 2.1 Lo que sí da GRH: el conjunto excepcional [C: GRH de Dirichlet]
Sea E(N) = #{n ≤ N par : n no es suma de dos primos}. Promediando sobre n, por la desigualdad de Bessel:

  Σ_{n ≤ N} |∫_𝔪 S² e(−nα) dα|² ≤ ∫_𝔪 |S(α)|⁴ dα ≤ (sup_𝔪 |S|)² · ∫₀¹ |S|².

Bajo GRH se tiene sup_𝔪 |S| ≪ N^{3/4+ε}, con arcos de nivel Q ≈ N^{1/2}. Por otro lado, cada n excepcional aporta ≫ N² al lado izquierdo. Resulta **E(N) ≪ N^{1/2+ε}** (Hardy–Littlewood, *Partitio Numerorum V*, 1924; esquema estándar en Vaughan, *The Hardy–Littlewood Method*, cap. 3, *verificar*).

**Lección estructural:** este argumento es de tipo L², **ciego a cada n por separado**. Cuenta cuántos n pueden fallar, pero nunca puede probar que no falla ninguno. Ningún refinamiento de esta línea llega a E(N) = 0.

**Si se acepta la prueba de Alpöge (GRH para caracteres de Dirichlet, ver `notas/2026-11-23_RH_Alpoge_relevancia_Goldbach.md`), esta cota pasa a ser INCONDICIONAL.** Goldston (1992) la refinó a N^{1/2}·(log N)^c bajo GRH (según el survey de Bhowmik–Grimmelt, arXiv:2607.27282).

**Cotas incondicionales actuales (abstracts comprobados el 2026-11-23):**
- Montgomery–Vaughan (1975): N^{1−δ}, con δ efectivo.
- Pintz (arXiv:1804.09084): E(N) < N^{0,72}.
- Zhao (arXiv:2511.05631): E(N) = O(N^{7/10}), con constante **inefectiva**.
- Pintz–Althöfer (arXiv:2610.15559): N^{7/12+ε} **condicional** a que exista un semiplano Re s > θ₀ (θ₀ < 1) común libre de ceros para todas las L de Dirichlet. Su abstract dice que el "teorema uniforme 7/8" de OpenAI (octubre de 2026, no verificado aquí) aportaría esa hipótesis.

**Reformulación reciente:** Li & Liu (2026-11-13) dicen dar una "prueba condicional" bajo GRH más hipótesis espectrales. Nuestro análisis (`02_analisis_Li_Liu_2026.md`) concluye que esas hipótesis equivalen a una cota de norma de operador al menos tan fuerte como la conclusión.

---

## 3. Cribas y el problema de paridad

- **Chen (1973) [T].** Todo par suficientemente grande es p + P₂, donde P₂ tiene como mucho dos factores primos.
- **Barrera de paridad (Selberg) [T, como fenómeno].** Una criba que solo use información de divisibilidad por primos pequeños no distingue los enteros con un número par de factores primos de los que tienen un número impar. Por eso no puede producir "n − p primo", solo "n − p casi primo".
- **La paridad sobrevive a hipótesis muy fuertes.** Bombieri (criba asintótica, 1976) mostró que ni siquiera la conjetura de Elliott–Halberstam permite a las cribas detectar primos en problemas binarios. **Conclusión: ni GRH, ni EH, ni las dos juntas bastan para Goldbach binaria con las técnicas conocidas.**
- **Lo que sí rompe la paridad**, siempre en contextos con estructura especial:
  - Friedlander–Iwaniec, primos de la forma x² + y⁴.
  - Heath-Brown, primos de la forma x³ + 2y³.
  - Correlaciones binarias de Liouville *promediadas logarítmicamente*: Tao (2016), Matomäki–Radziwiłł.
  - Correlaciones Λ(n)Λ(n+h) para *casi todo* h en rangos largos: Matomäki–Radziwiłł–Tao, *verificar alcance*.
  En todos estos casos hay un promedio extra o una estructura algebraica. Goldbach binaria, para un n fijo, no tiene ninguna de las dos cosas.
- **Análisis de Fourier de orden superior (Green–Tao).** Cubre sistemas lineales de complejidad finita. El par {m, n − m} es afínmente dependiente, es decir, de "complejidad infinita", así que **queda fuera**.

---

## 4. El análogo en cuerpos de funciones: la pista nueva

En F_q[t], los análogos de los problemas binarios de primos **sí están demostrados**:
- **[T] Sawin–Shusterman**, *On the Chowla and twin primes conjectures over F_q[T]*, Annals of Math. 196(2) (2022), 457–506, arXiv:1808.04001. El Teorema 1.1 da los primos gemelos en forma cuantitativa de Hardy–Littlewood para p primo impar y q potencia de p con **q > 685090·p²**. El **Remark 1.2**, comprobado en lectura literal (verificado 2026-11-23), dice: *"Our proof of Theorem 1.1 establishes also the analog of the Goldbach problem over function fields."* Es decir, **Goldbach binaria está demostrada en F_q[T] con q fijo** cumpliendo esa condición. El método obtiene un **nivel de distribución mayor que 1/2** para los polinomios irreducibles, mediante una variante en cuerpos de funciones de un resultado de Fouvry–Michel (que se apoya en cotas de tipo Deligne, *verificar detalle*).
- **[T] Para q grande respecto del grado:** Hayes, y Bender (arXiv:0809.4381 y 0912.1702), este último con el número de representaciones que predice la heurística cuando q > 8(n+6)^{2n²}.
- **Mecanismo que esquiva la paridad en F_q[T] (conocimiento del modelo, *verificar*):** la fórmula de Pellet expresa μ(f) como un carácter cuadrático del discriminante de f. Eso convierte la paridad en una suma de caracteres que se puede atacar con geometría. Sobre ℤ no se conoce ningún análogo.

**Por qué importa ahora.** La prueba de Alpöge de RH/GRH es, según su propia nota, un **descendiente del argumento de Stepanov** (RH sobre cuerpos de funciones) llevado a los enteros. La pieza que falta sobre ℤ, la "multiplicidad", se sustituye con levantamientos theta iterados. Su nota menciona además una "prueba de Weil y Deligne solo con signos" en un modelo de juguete.

**Coincidencia que conviene estudiar [P].** En la nota 3 de su texto, Alpöge cuenta que buscaba simplificar el *"averaging over the family χ ↦ χ(disc)"*, es decir, promedios de caracteres evaluados en el discriminante. Es el mismo tipo de objeto que interviene en el mecanismo de paridad de F_q[T].

**Pregunta de investigación [P, especulativa].** La prueba de Sawin–Shusterman de Goldbach en F_q[T] existe y es [T]. ¿Usa **solo ingredientes que la técnica de transferencia de Alpöge sepa llevar a ℤ**? Si existe, sería la primera ruta plausible que no choca de frente con la barrera de paridad. Ahí la paridad se esquivaría con geometría o automorfía, no con cribas. Ahora mismo **no es evaluable**: Alpöge todavía no ha publicado una exposición legible del método.

---

## 5. Trampas típicas (lista de control para revisar cualquier "prueba")

Las supuestas demostraciones de Goldbach suelen fallar por alguno de estos motivos. **Todo agente debe revisar esta lista antes de proponer un paso.**
1. **"Casi todos" tomado por "todos".** Densidad, promedio o conjunto excepcional pequeño no implican conjunto vacío.
2. **Independencia probabilística no justificada.** "n − p es primo con probabilidad 1/log n" es una heurística [H], no un lema.
3. **Términos de error que no son uniformes en n.** Una cota O(·) con constante que depende de n no sirve.
4. **Parseval o Cauchy–Schwarz en los arcos menores.** Pierde exactamente el factor log n de la sección 2.
5. **Cribas con información solo de nivel de distribución.** La paridad las bloquea (sección 3).
6. **Circularidad.** Usar un lema equivalente a Goldbach, por ejemplo una cota de arcos menores para cada n.
7. **Comprobar con ejemplos pequeños.** Los récords de p_min crecen lentamente, pero no está demostrado que estén acotados por nada útil.
8. **Pasos "obvios" sin enunciado preciso.** Si un paso no se puede escribir como un lema con hipótesis y conclusión exactas, no es un paso.

---

## 6. Programa de trabajo propuesto

| Línea | Objetivo | Riesgo | Valor si sale bien |
|---|---|---|---|
| **L1** | Reunir con enunciados exactos y referencias las consecuencias de GRH para Goldbach y estudiar cuáles pasan a ser incondicionales con Alpöge: E(N) ≪ N^{1/2+ε}; números de Goldbach en intervalos cortos; constante K de Linnik para "p + q + K potencias de 2". Para K: 13 incondicional y 7 bajo GRH (Heath-Brown–Puchta 2002); 8 incondicional (Pintz–Ruzsa); 6 bajo GRH (Johnston–Trudgian, arXiv:2605.17825); **6 incondicional y 5 bajo GRH** (Huixi Li y Zhixin Liu, arXiv:2610.15784, sin arbitrar). Con Alpöge, **K = 5 incondicional** | Bajo | Resultados reales, aunque **no son Goldbach** |
| **L2** | Formalizar en Lean 4 + Mathlib: enunciado de Goldbach, reducción 1.1 y, a largo plazo, "GRH ⇒ E(N) ≪ N^{1/2+ε}" | Medio (requiere Lean instalado) | Verificación por máquina; no depende de que nos fiemos unos de otros |
| **L3** | Estudiar la literatura de cuerpos de funciones (sección 4) y, cuando salga la exposición de Alpöge, su principio de transferencia | Alto | Única ruta nueva con alguna plausibilidad |
| **L4** | Experimentos numéricos (`numerico/goldbach.py`): p_min, r(n)/predicción de Hardy–Littlewood, margen mínimo. **Sirven sobre todo para refutar rápido lemas falsos** que propongan los agentes | Bajo | Filtro de calidad |

**Regla de oro:** todo lema propuesto debe (1) tener un enunciado preciso, (2) probarse numéricamente cuando sea posible, (3) recibir una revisión crítica que busque contraejemplos y (4) formalizarse en Lean cuando el lema sea central.
