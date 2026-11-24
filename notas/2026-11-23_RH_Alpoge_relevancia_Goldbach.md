# Noticia: demostración anunciada de la hipótesis de Riemann (nov. 2026) y su relevancia para Goldbach

## Procedencia
- **Autor de la nota:** Claude Code (modelo `claude-opus-5-5`), sesión local del usuario.
- **Fecha:** 2026-11-23.
- **Método:** búsqueda web y lectura directa de fuentes. Nature y Science no se pudieron leer (muro de pago). La nota de Alpöge (`alpo.ge/remark.pdf`) se leyó extrayendo el texto del PDF.
- **Estado:** PROPUESTA, no verificada de forma independiente. No es un hecho consolidado del proyecto.
- **Clasificación sugerida:** `internal` (son noticias públicas y no contienen datos personales).
- **Registro en `lumina-mcp`:** PENDIENTE. El MCP `lumina_shared_memory` dio "Connection closed" en esta sesión. Esta nota es el respaldo local que pide el punto 8 del contrato.
- **Posible conflicto de interés:** Claude es de Anthropic, la misma empresa del autor del anuncio. Por eso esta nota se apoya en fuentes externas: Lean FRO/Mathlib, madrimasd, TechCrunch, el blog de Vakil y un panel independiente en GitHub.

---

## 1. Hechos con fuente verificable

| Hecho | Fuente |
|---|---|
| El 20-11-2026, Levent Alpöge (matemático de Anthropic) publicó una nota de 4 páginas. Dice: *"I believe I have, working with Claude, discovered a proof of the Riemann Hypothesis"*. Añade que el método resuelve también la **hipótesis generalizada (GRH)** y que espera llegar a la extendida (ERH) "con más tiempo". | alpo.ge/remark.pdf (leído) |
| La prueba **solo existe en Lean**: unos 17 millones de líneas, repartidas en `alpo.ge/lean.tar.gz.part01…part14`. No hay prueba escrita para humanos. Alpöge dice que la versión completa en texto pasaría de 250 000 páginas y que piensa escribir una exposición legible. | Nota de Alpöge; madrimasd |
| El 21-11-2026, el **Lean FRO y la Mathlib Initiative** publicaron una declaración conjunta: *"We hereby certify that both the statement and the proof of the Riemann Hypothesis are correct."* El enunciado es el `RiemannHypothesis` de Mathlib. Una "inert export" de la prueba pasó por el kernel oficial de Lean, por los verificadores `con-leche` y `con-ron` y por otros verificadores independientes. | Blog de leanprover-community |
| La declaración **no** certifica que la matemática se entienda: *"We look forward to all the activity in the mathematical community that will help explain the new ideas in this proof."* | Mismo blog |
| Según el panel independiente `nasqret/riemann-rev`, la prueba solo usa los 3 axiomas estándar de Lean (`propext`, `Classical.choice`, `Quot.sound`). Dos kernels independientes (Lean y NanoDA, en Rust) aceptaron la prueba. Recompilaron con éxito Sutherland (MIT), DeDeo, Arledge y Teżewski, entre otros. | github.com/nasqret/riemann-rev |
| Según ese mismo panel, el teorema formal cubre **GRH para todos los caracteres de Dirichlet χ mod N**, además de RH para ζ. | Mismo panel. Ojo: parte del contenido se obtuvo de un resumen automático, no de lectura literal. |
| Kevin Buzzard planteó si el enunciado en Lean es realmente RH "o uno sutilmente distinto". La comparación de enunciados del Lean FRO responde a esa duda en el caso de `RiemannHypothesis`. | madrimasd 22-11 |
| El Instituto Clay exige publicación y 2 años de aceptación general, así que el premio no podría concederse antes de finales de 2028. | madrimasd 22-11 |

### Demostración rival "toda humana" (es otro caso, no mezclar)
- El 20-11-2026, Ravi Vakil (presidente de la AMS) publicó un "Hurried announcement": Alessandro Masullo (doctor por Stanford, dirigido por Brian Conrad) propone una prueba de GRH, *"all human, with no AI"*. El manuscrito se titula *The Music of the Primes* y tiene unas 960 páginas (dato de madrimasd).
- El 22-11, Vakil añadió: *"Several people have contacted me to report problems in the proposed proof, identified with the assistance of different AI models"*. También: *"I make no claim as to the correctness of the proof strategy."*
- **Estado:** abierto. No está formalizada en Lean.

### Contexto previo (según madrimasd, sin verificar aquí)
- Agosto de 2026: una versión de investigación de Claude elevó la proporción demostrada de ceros en la línea crítica (del 41,6 % al 67,2 %).
- Octubre de 2026: OpenAI anunció un resultado "cuasi-Riemann", con ceros de parte real entre 1/8 y 7/8.

## 2. Estrategia de la prueba de Alpöge (según su nota; él mismo dice que es un resumen "con pérdida")
- Es un descendiente del **método de Stepanov** (con el que se prueba RH sobre cuerpos de funciones). Sobre los enteros falta la "multiplicidad", y se sustituye con **cadenas iteradas de levantamientos theta** entre grupos unitarios totalmente definidos, U(X₀) → … → U(X_m). Cada levantamiento "imita" una multiplicidad, y se toma m ≈ √X.
- La prueba trabaja con una suma de von Mangoldt con peso en rampa, Σ_{n≤X} Λ(n)(X/n − 1)₊ − X log X + (1+γ)X, y con la positividad de Guinand–Weil.
- También usa: períodos GGP, el seesaw de Rallis, Chowla–Selberg, la desigualdad de pendientes de Bost, una torre de Iwasawa (elegida con ayuda de Behrend), una forma debilitada de Vandiver y una fórmula de traza relativa.
- La cronología del propio autor va del 06-10-2026 (el "determinante de Stepanov") al 16-11-2026 ("GRH").

## 3. Relevancia para Goldbach (lo más importante para este proyecto)

> Esta sección se basa en **conocimiento matemático estándar del modelo**, no en las fuentes de arriba. Antes de citarla en algo formal, hay que comprobar cada referencia.

1. **La prueba de RH/GRH NO demuestra la conjetura fuerte (binaria) de Goldbach.** No se conoce ninguna implicación GRH ⇒ Goldbach binaria. El método del círculo aplicado a problemas binarios no controla los "arcos menores" punto a punto, ni siquiera suponiendo GRH; solo los controla en media cuadrática. Además, las cribas siguen bloqueadas por el **problema de paridad**.
2. **Goldbach débil (ternaria):** ya estaba demostrada sin condiciones (Helfgott, 2013). Antes se había probado bajo GRH (Deshouillers–Effinger–te Riele–Zinoviev, 1997). La noticia no añade nada aquí.
3. **Lo que SÍ cambiaría si GRH de Dirichlet queda aceptada:**
   - Hardy–Littlewood (1924, *Partitio Numerorum V*) probaron, suponiendo GRH para las L de Dirichlet, que el conjunto excepcional E(N) (pares ≤ N que no son suma de dos primos) cumple **E(N) ≪ N^{1/2+ε}**. Con el teorema de Alpöge, esto pasaría a ser **incondicional**. La mejor cota incondicional conocida está en torno a N^{0,72} (Pintz; verificar referencia).
   - Varios resultados de "números de Goldbach en intervalos cortos" de longitud polilogarítmica (Linnik y sucesores; verificar) también pasarían a ser incondicionales.
   - Pero "casi todos los pares" **no es lo mismo que** "todos los pares". Aunque E(N) sea pequeño, eso no prueba que sea vacío.
4. **Verificación numérica conocida:** Goldbach se ha comprobado hasta 4·10¹⁸ (Oliveira e Silva, Herzog y Pardi, 2014). Es evidencia, no demostración.

## 4. Líneas de trabajo que se proponen al equipo (para discutir; nada decidido)
- **A. Formalizar en Lean 4 + Mathlib la implicación "GRH de Dirichlet ⇒ E(N) ≪ N^{1/2+ε}"**, usando como hipótesis el teorema de Alpöge. Es concreto y verificable por máquina, y daría un avance real y formal hacia Goldbach. Aun así, **no sería la conjetura**.
- **B. Infraestructura numérica:** un verificador de Goldbach reproducible y estadísticas de representaciones r(n) frente a la predicción de Hardy–Littlewood. Sirve para probar ideas rápido; su valor demostrativo es nulo.
- **C. Estudio del método de Alpöge:** ver si la técnica de "multiplicidad por levantamientos theta" dice algo sobre problemas aditivos binarios. Es **muy especulativo**. Hay que esperar a que exista una exposición legible.
- **Regla común:** separar siempre *hipótesis*, *evidencia numérica* y *demostración*. Ninguna revisión entre agentes que compartan usuario de Windows cuenta como revisión independiente (contrato, punto 5).

## Fuentes
- Nota de Alpöge: https://alpo.ge/remark.pdf
- Declaración Lean FRO / Mathlib: https://leanprover-community.github.io/blog/posts/riemann-hypothesis/
- Panel independiente: https://github.com/nasqret/riemann-rev
- Vakil sobre Masullo: https://math216.wordpress.com/2026/11/20/hurried-announcement-alessandro-masullo-grh/
- madrimasd (22-11): https://www.madrimasd.org/blogs/matematicas/2026/11/22/151406
- madrimasd (23-11): https://www.madrimasd.org/blogs/matematicas/2026/11/23/151425
- TechCrunch (20-11): https://techcrunch.com/2026/11/20/ai-claimed-a-solution-to-the-167-year-old-riemann-hypothesis-and-a-human-did-too/
- No leídas (muro de pago): Nature https://www.nature.com/articles/d41586-026-03812-7 · Science https://www.science.org/content/article/ai-may-have-solved-math-s-greatest-mystery-can-humans-understand-answer
- No leída: página de verificación de Sutherland, https://math.mit.edu/~drew/grh/
