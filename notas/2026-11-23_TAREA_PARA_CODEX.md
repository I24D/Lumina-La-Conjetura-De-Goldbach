# Tarea para Codex: continuar el proyecto Goldbach (relevo de Claude Code)

**Procedencia:** redactado por Claude Code (`claude-opus-5-5`) el 2026-11-23, a petición de Dal. Es un encargo, no una afirmación matemática.
**Destinatario:** Codex.

## 0. Antes de empezar (obligatorio)
1. Lee `README.md` y `ESTADO.md`, en especial las entradas fechadas de Claude Code al final, y los documentos `matematicas/06_descomposicion_bilineal.md` y `matematicas/07_sesgo_mobius_goldbach_y_casi_todo_N.md`.
2. **Reglas del proyecto, sin excepción:**
   - no sobrescribas archivos de otros agentes; añade entradas nuevas o crea archivos nuevos;
   - pon cabecera de procedencia (agente, modelo, fecha, método, estado);
   - etiqueta cada afirmación con [T], [T-Lean], [P], [P-rev], [N], [C: H], [H] o [FALSO];
   - la evidencia numérica no es demostración;
   - si encuentras un error en otro documento, **repórtalo explícitamente**: no lo corrijas en silencio.
3. **Lean 4 + Mathlib ya está instalado.**
   - Proyecto: `C:\Users\dal_n\lean\goldbach_lumina`, fuera de OneDrive a propósito.
   - elan está en `C:\Users\dal_n\.elan\bin`. Si `lake` no aparece en el PATH, abre una terminal nueva o antepón esa ruta.
   - Compilar: `cd C:\Users\dal_n\lean\goldbach_lumina` y luego `lake env lean GoldbachLumina/<Archivo>.lean`.
   - Usa `ArithmeticFunction.moebius` y `ArithmeticFunction.vonMangoldt` con nombre completo; μ/Λ dan ambigüedad.
   - Para comprobar que no hay huecos, añade `#print axioms <teorema>`: deben salir solo `propext`, `Classical.choice` y `Quot.sound`.
   - Si escribes archivos desde PowerShell 5.1, hazlo **sin BOM**; Lean falla con BOM.
4. **Memoria compartida:** si el MCP `lumina_shared_memory` funciona en tu sesión, registra como **propuesta** los aportes pendientes de Claude (docs 06 y 07; ver ESTADO). Si no funciona, déjalo anotado como PENDIENTE.

## 1. Estado en una frase
La reducción de 04 está verificada en Lean: `R_eq_S_sub_T` y 5 lemas más, en `lean/README.md`. Todo el obstáculo está en la **parte de Tipo II / arcos menores** de la cola, con d > √N (06, F1 y F2). Allí hay que distinguir los d con ω(d) par de los d con ω(d) impar: es la barrera de paridad. En la caja c = 1 la cola contiene un bloque diádico de la propia suma de Goldbach (06, Lema 6.0 y R-3).

## 2. Tus tareas, por orden de prioridad

### T1. Revisión adversarial independiente (prioridad máxima)
Revisa e intenta **refutar**:
- **06, Lema 6.1.** ¿Es correcto el uso de BV, con máximo en y y en la clase, para d ≤ √N/(log N)^B? ¿Y el tratamiento del término principal mediante μ·1_{(·,N)=1}·id/φ = μ * h_N?
- **06, Lema 6.2** (Tipo I de Vaughan). ¿Es correcto que (q, c) = 1? ¿Está bien la reducción de clases no reducidas en BV para μ?
- **06, R-2 y R-3.**
  - ¿Vale la Prop. 4.6 de **tu** documento 04 con D = √N/(log N)^B?
  - ¿Es correcto que T' (la parte con m = cd compuesto) es igual a N𝔖(N) + o(N) sin condiciones?
- **07, Teorema 7.1** ((HB) para casi todo N). Revisa:
  - Parseval sobre N;
  - la aplicación de Davenport con θ = cα;
  - que la condición (c, N) = 1 se elimine correctamente;
  - el recuento de excepciones.
- **Formato:** un documento nuevo, `matematicas/08_revision_codex_06_07.md`. Para cada punto, veredicto: correcto / error / hueco, con el detalle. Recuerda que tu revisión es [P-rev], no independiente, porque compartimos usuario de Windows.

### T2. Búsqueda bibliográfica sobre el fenómeno de 07 §§2–3
¿Está publicado el sesgo de clase de M0(N) = Σ_{n<N} μ(n)Λ(N−n) (media de Z = −0,22 si 3 | N; correlación −0,49 con 𝔖(N), N ≤ 10⁶)? ¿Y su explicación por la fórmula de arcos mayores P_Q(N) = Σ_{q≤Q} μ(q)/φ(q)·Σ_{n<N} μ(n)c_q(N−n)?
- Busca correlaciones μ×Λ, desarrollos de Ramanujan de correlaciones (Coppola, Murty, Saha y otros), sesgos de tipo Chebyshev en problemas binarios y Möbius–Goldbach.
- Informa con fuentes primarias. Si ya está publicado, márcalo como conocido y cita la referencia. **No afirmes que es nuevo sin evidencia.**

### T3. Reproducir y ampliar el experimento de 07
1. Ejecuta `python numerico/mobius_goldbach_aleatoriedad.py --hasta 1e6` y `python numerico/arcos_mayores_mobius.py --hasta 1e6 --Q 1 6 30`, y comprueba que reproduces las tablas de 07.
2. **Pregunta abierta clave:** ¿el residuo M0 − P_Q (parte de arcos menores) tiene algún **sesgo de signo** o asimetría que pueda servir para la cota unilateral (HB+) de 05, es decir, T ≤ (1−κ)N𝔖(N)? Mide la asimetría y las colas del residuo normalizado, por separado para N par, y su dependencia en N (10⁵, 10⁶, 10⁷).
3. Haz lo mismo con la **cola completa** T_N^prim (todas las c, no solo c = 1), usando `numerico/cajas_bilineales.py` o `numerico/cola_mobius_scan.py`.
4. Guarda los JSON con procedencia en `numerico/resultados/`.

### T4. Formalización en Lean (si te queda tiempo)
Siguientes piezas elementales, en un archivo nuevo `GoldbachLumina/Reduccion2.lean`. No modifiques `Reduccion.lean`.
- La reindexación (c, d) de 04: T_N(D) en forma de divisores = Σ_c Σ_{D<d≤(N−1)/c} μ(d) log d Λ(N−cd).
- La cota elemental (4.12) de 04 sobre potencias propias.
- La cota de la Prop. 4.2 de 04: |T^bad| ≤ K_N.

Copia las fuentes que compilen a `lean/` en OneDrive y añade una fila a `lean/README.md` con el resultado de `#print axioms`.

## 3. Lo que NO hay que hacer
- No persigas (HB) caja por caja en la región K > √N. En las cajas con C pequeño es Goldbach por bloques (06, R-3).
- No uses Bombieri–Vinogradov por encima de su nivel (√N/(log N)^B) ni cancelación "esperada".
- No declares Goldbach demostrada, ni total ni condicionalmente, sin una demostración completa revisada.

## 4. Al terminar
Añade una entrada fechada al final de `ESTADO.md` con lo demostrado, lo refutado, lo pendiente y los archivos, y una línea en `README.md`.
