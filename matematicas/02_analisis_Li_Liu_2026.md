# 02 — Análisis crítico de Li & Liu, *A Conditional Proof of the Goldbach Conjecture* (13-11-2026)

**Procedencia:** Claude Code (`claude-opus-5-5`), 2026-11-23. Estado: **[P]**, análisis propio sin revisar.
**Fuente:** Jiamin Li y Jianya Liu (Universidad de Shandong), PDF de 63 páginas publicado en https://mrc.sdu.edu.cn/english/info/1018/1514.htm el 2026-11-13. No está en arXiv ni ha pasado revisión por pares. Se leyeron el texto extraído de las págs. 1–15 (resumen, §1, §2 completa y el inicio de §3) y los enunciados localizados en el resto. **No se han revisado las demostraciones de §§3–7** (construcción del operador, Teorema 6.2, Lema 7.2, Proposición 7.4).

---

## 1. Qué afirman

Fijan un peso suave ϕ con soporte en [1/3, 2/3], simétrico y con 0 ≤ ϕ ≤ 1. Definen
  G(N) = (1/N) Σ_{n+m=N, (nm,N)=1} Λ(n)Λ(m) ϕ(n/N) ϕ(m/N),  c_ϕ = ∫₀¹ ϕ².

- **GRH(N):** todos los ceros no triviales de las L de Dirichlet primitivas con conductor que divide a N tienen parte real 1/2.
- **Construcción** (para N > e¹⁰, suponiendo GRH(N)). Con los caracteres impares mod N obtienen una proyección exacta de G(N). Restan la diagonal y un "modelo Beta" y queda un residuo E_N(s, z). Ese residuo define un operador **compacto, autoadjunto y de Hilbert–Schmidt** T_N sobre un espacio de Hardy.
- **Identidad exacta (2.39):** c_ϕ 𝔖(N) − G(N) = ⟨T_N v_N, v_N⟩, con ‖v_N‖² = Q_N ≪_ϕ 1.
- **Cota maestra (2.19):** |G(N) − c_ϕ 𝔖(N)| ≪_ϕ ‖T_N‖, por Cauchy–Schwarz.
- **Lo que da GRH sola (2.24, Teorema 6.2):** ‖T_N‖ ≤ ‖T_N‖_HS ≪ (log N)². No basta, porque el término principal es c_ϕ 𝔖(N) ≍ 1 (como mucho del orden de log log N).
- **"Función espectral" (1.6):** X_N(q) = ∏_j (1 − 4 μ_j² δ_N^{−2} q(1−q)), con q cuaternión, μ_j los autovalores de T_N y δ_N² = 𝔖(N)(log N)³/N.
- **Hipótesis S:** X_N no tiene ceros con Re q > 1 − η_N, donde η_N = h(N)²(log N)³/(N 𝔖(N)). Da el Teorema 1.1, G(N) ∼ c_ϕ 𝔖(N).
- **Hipótesis R:** todos los ceros de X_N tienen Re q = 1/2. Da el Teorema 1.2, con error O(δ_N).

## 2. Observación central: las hipótesis son la conclusión reescrita [P]

Los ceros de X_N se obtienen **algebraicamente** a partir de los autovalores. Cada factor 1 − 4μ²δ^{−2} q(1−q) se anula en q(1−q) = δ²/(4μ²). Esa raíz está en Re q = 1/2 si |μ| ≤ δ, y es real (fuera de la recta) si |μ| > δ. **El propio artículo lo dice** en su "diccionario exacto" (§2.4.2):

  X_N(q) ≠ 0 en Re q > 1 − η  ⟺  ‖T_N‖ ≤ δ_N / (2√(η(1−η))),

  **(2.28)** Hipótesis S(N; h) ⟺ ‖T_N‖ ≤ 𝔖(N) / (2h(N)√(1 − η_N)),
  **(2.30)** Hipótesis R(N) ⟺ ‖T_N‖ ≤ δ_N.

En consecuencia:
1. **Las "hipótesis espectrales" son cotas de la norma de operador de T_N con otro nombre.** La función cuaterniónica X_N no aporta información: es un cambio de notación. Sus "ceros" no son objetos independientes, como sí lo son los ceros de ζ.
2. **La hipótesis es al menos tan fuerte como la conclusión.** S pide ‖T_N‖ = o(𝔖(N)), y la conclusión solo necesita |⟨T_N v_N, v_N⟩| = o(𝔖(N)) para **un** vector v_N. Como la norma es un supremo sobre todos los vectores, la hipótesis es *más fuerte* que lo que se quiere demostrar.
3. **La distancia hasta lo conocido es grande.** GRH da ‖T_N‖ ≪ (log N)², y S pide o(1) por un factor h(N) → ∞. El salto es un factor h(N)(log N)² / 𝔖(N): **peor** que el salto clásico de los arcos menores del método del círculo, que es un factor log N (ver `01`, §2). El propio artículo cuantifica el salto en (2.29): η_N/η_N^GRH ≍ h²(log N)⁴/𝔖² → ∞.
4. **Ni siquiera está claro que la hipótesis sea cierta.** ‖T_N‖ pequeña exige que la forma cuadrática sea pequeña en *todos* los vectores del espacio de Hardy, no solo en el vector asociado al peso de Goldbach. Si algún vector "malo" diera un valor grande, S y R serían falsas aunque Goldbach fuera cierta. **No hemos comprobado si pasa.**

**Conclusión:** el trabajo no es una reducción de Goldbach a una hipótesis más débil o independiente. Es una **reformulación**: Goldbach (para N grande, bajo GRH) se sigue de una cota de norma de operador que es al menos igual de difícil. El título "A conditional proof" es técnicamente correcto, pero puede inducir a error, porque la condición contiene la dificultad entera.

## 3. Lo que sí puede servirnos [P]

1. **La identidad exacta (2.39)**, si es correcta, escribe el error de Goldbach como una forma cuadrática de un operador autoadjunto compacto. Con la GRH de Alpöge, la construcción sería **incondicional**. El artículo dice además que su Apéndice B la comprueba numéricamente en N = 30030.
2. **Basta una cota de un solo lado (observación nuestra).** Para que haya representaciones solo hace falta que G(N) > 0, es decir, ⟨T_N v_N, v_N⟩ < c_ϕ 𝔖(N) − (aporte de potencias de primos, ≪ N^{−1/2} log² N). Como T_N es autoadjunto, ⟨T v, v⟩ ≤ λ⁺_max(T) ‖v‖², donde λ⁺_max = max(0, mayor autovalor). Por tanto

  **Goldbach(N) ⇐ λ⁺_max(T_N) · Q_N < c_ϕ 𝔖(N) − O(N^{−1/2} log² N).**

  **No hace falta controlar los autovalores negativos.** Una cota superior para G(N), del tipo de la criba de Selberg, controla en v_N la dirección contraria y es fácil. La dificultad, es decir, el **problema de paridad, queda localizada en la parte positiva del espectro de T_N**. Esto no hace fácil el problema, pero lo plantea en forma de **positividad**, que es justo el tipo de enunciado que atacan los argumentos de signo de Weil, Guinand–Weil y Stepanov. Según su nota, son la base del método de Alpöge.
3. **Pregunta concreta derivada [P]:** ¿se puede acotar λ⁺_max(T_N) con un argumento de positividad (determinantes, núcleos definidos positivos) que no pase por ‖T_N‖? Es la versión de la línea L3 que más se puede atacar dentro de este marco.

## 4. Tareas que se proponen
- **T-LL1:** revisar las demostraciones de §§3–7, sobre todo la construcción de T_N y la identidad (2.39), y si la constante de Q_N es de verdad uniforme en N.
- **T-LL2 (numérica): HECHA en parte (Claude, 2026-11-23) [N].** Script `numerico/li_liu_check.py`; resultado en `numerico/resultados/20261123_204322_li_liu_check_10000000.json`.
  - Peso: bump C^∞ con soporte en [1/3, 2/3]; c_ϕ = 0,163897.
  - Muestra: 40 valores de N por tramo diádico, desde 2¹⁰ hasta 10⁷.
  - Resultado: **max |G(N) − c_ϕ𝔖(N)|/δ_N entre 0,14 y 0,23, sin tendencia creciente**, con media entre 0,05 y 0,07. El error relativo baja de 0,56 a 0,015.
  - Lectura: la *conclusión* del Teorema 1.2 encaja con los datos, como cabía esperar, porque es Hardy–Littlewood con error de raíz cuadrada. **Esto no contrasta la hipótesis R** (‖T_N‖ ≤ δ_N), solo su consecuencia.
  - Pendiente: reproducir el Apéndice B (N = 30030).
- **T-LL3:** buscar si un vector de prueba sencillo da ⟨T_N w, w⟩ grande. Si lo da, refutaría S y R tal como están enunciadas.
- **T-LL4:** estudiar el signo de la parte positiva del espectro de T_N (§3, punto 2).
