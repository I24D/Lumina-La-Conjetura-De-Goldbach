# 03 — Sawin–Shusterman, Alpöge y el paso que falta sobre los enteros

## Procedencia y alcance

- Agente: **Codex**. Modelo exacto: no expuesto por esta interfaz; no se atribuye otra versión.
- Fecha del proyecto: **2026-11-23**, America/New_York.
- Método: lectura previa de README, ESTADO, documentos 01–02 y nota de Alpöge del proyecto; consulta directa de las fuentes primarias enlazadas; derivaciones escritas y comprobaciones finitas con biblioteca estándar de Python.
- Estado global: **[P]**, aporte con demostraciones escritas, pendiente de revisión independiente y formalización. No es una demostración de Goldbach ni una afirmación de novedad matemática.
- Fuentes principales: [Sawin–Shusterman, arXiv:1808.04001v2](https://arxiv.org/pdf/1808.04001v2), [nota de Alpöge](https://alpo.ge/remark.pdf), [verificaciones publicadas por Sutherland](https://math.mit.edu/~drew/grh/), [declaración Lean FRO/Mathlib](https://leanprover-community.github.io/blog/posts/riemann-hypothesis/).
- Memoria compartida: MCP de propuestas no expuesto en esta sesión. Registro local; subida como propuesta **PENDIENTE**, sin utilizar service_role ni escritura REST alternativa.

## 1. Cadena publicada que sí podemos identificar [T]

En la versión v2 de Sawin–Shusterman:

| Referencia | Pieza pertinente |
|---|---|
| Teorema 1.1; Remark 1.2, pág. 2 | Primos gemelos y análogo de Goldbach, con q fijo, grado creciente, característica p impar y q > 685090 p². |
| Teorema 2.1, pág. 8 | Para módulo cuadrado libre de grado m y carácter no principal: suma sobre h de grado < t acotada por (√q + 1)·binom(m−1,t)·q^(t/2). |
| Lema 3.2 y (3.7), págs. 15–16 | Pellet y reciprocidad convierten μ(a+gM) en S·χ(w+g); el carácter y su conductor dependen de la derivada. |
| Teoremas 1.3, 1.7–1.9; §4 | Correlaciones, control de excepciones y distribución; no basta con citar una sola cota de caracteres. |
| Corolario 6.1 y §6.3, págs. 27–30 | Cancelación de productos de μ ponderados por Λ; se aplica al término de divisores grandes. |

El punto algebraico es que (r+s^p)' = r'. La paridad de la factorización se vuelve accesible en fibras con derivada fija. El control geométrico, los conductores y las familias excepcionales siguen siendo esenciales. La formulación general del Lema 3.2 tiene hipótesis de coprimalidad y grados: no se afirma una identidad global μ=χ.

Para un f mónico de grado positivo, en característica impar, la identidad de Pellet usada es

    μ(f)=(−1)^deg(f) χ_q(Disc(f)),  χ_q(0)=0.

En el Lema 3.2, cuando D₀≠0, intervienen

    D₀=M²g′+a′M−aM′,
    E=rad(D₀)/gcd(M,rad(D₀)),
    E₁=rad_impar(D₀)/gcd(M,rad_impar(D₀)).

E es el módulo y E₁ el conductor del carácter. E₁ puede ser 1: no es lícito dar por no principal todo carácter obtenido de un discriminante. Es parte del motivo para controlar las fibras excepcionales antes de usar la cota geométrica.

## 2. Fibras y obstáculos algebraicos [P, pruebas elementales]

### 2.0 Por qué la fibra de la derivada es exactamente r+s^p

**Enunciado.** En 𝔽_q[T] de característica p, ker(d/dT)={s^p:s∈𝔽_q[T]}.

**Prueba.** h′=0 obliga a que todo exponente con coeficiente no nulo sea múltiplo de p: h=Σ_j a_jT^(pj). Como el Frobenius b↦b^p es biyectivo en el cuerpo finito, elegimos b_j con b_j^p=a_j. En característica p, (Σ_j b_jT^j)^p=Σ_j a_jT^(pj)=h. La implicación inversa se sigue de (s^p)′=0. Por tanto f′=r′ equivale a f=r+s^p. □

Este paso no es una estimación probabilística: es una identidad algebraica exacta. Los dos lemas siguientes delimitan por qué no se copia mediante una sustitución formal en ℤ.

### 2.1 No hay reducción entera sobreyectiva a los cuerpos que requiere ese teorema

**Enunciado.** Si φ: ℤ → 𝔽_(p^a) es un homomorfismo de anillos que preserva 1, su imagen es exactamente 𝔽_p. Por tanto no es sobreyectivo si a>1.

**Prueba.** Por aditividad y φ(1)=1, φ(n)=n·1 para todo n entero, también los negativos. En característica p estos valores constituyen precisamente el subcuerpo primo. Su cardinalidad es p, menor que p^a si a>1. □

**Consecuencia.** q=p no cumple q>685090 p². Los cuerpos que sí lo cumplen tienen a>1 y no se obtienen como imagen sobreyectiva de ℤ mediante reducción modular. Además ℤ/p^aℤ tiene característica p^a y divisores de cero si a>1; no es 𝔽_(p^a).

Esto refuta la transferencia directa «reducimos los enteros módulo q y aplicamos el teorema». **No** excluye correspondencias geométricas o automorfas más profundas.

### 2.2 La derivación ordinaria de ℤ no reproduce esas fibras

**Enunciado.** Toda derivación aditiva D:ℤ→M, con M un ℤ-módulo ordinario y D(ab)=aD(b)+bD(a), es cero.

**Prueba.** D(1)=D(1·1)=2D(1), luego D(1)=0. La aditividad da D(n)=nD(1)=0 para enteros positivos; D(0)=0 y D(−n)=−D(n) cubren los restantes. □

Por tanto, no podemos reemplazar T por un entero y conservar literalmente el argumento de derivadas. «Derivadas aritméticas» con reglas distintas necesitarían una teoría adicional; este lema no las prohíbe.

## 3. Qué consta de Alpöge y qué no [T-Lean: comunicado externo / P: comparación]

La página primaria de Sutherland enumera cinco declaraciones: RH y varias formulaciones de GRH de Dirichlet, para todos los módulos y caracteres. Publica auditorías de axiomas y resultados de distintos verificadores; algunos procesos concretos seguían en curso. Excluye expresamente una conclusión general sobre funciones L automorfas, de Hecke o zetas de Dedekind. La declaración Lean FRO/Mathlib certifica RH. **Aquí se verificaron documentalmente esos informes; no se ejecutaron sus kernels.**

La nota de Alpöge describe multiplicidad mediante levantamientos theta iterados y menciona caracteres de discriminantes. Es una indicación de método, no un enunciado preciso de transferencia de las estimaciones de §1.

**Inferencia delimitada [P].** Las cinco declaraciones consultadas no enuncian control de Σ μ(d)Λ(N−cd), ni una desigualdad para cada N par. La coincidencia de objetos en ambas narraciones no acredita tal control. Tampoco hemos auditado el archivo completo de Alpöge: no afirmamos que carezca de todo lema pertinente.

El documento 04 deriva el objeto entero que habría que buscar o construir. No exige, por principio, trasladar toda la geometría: exige una estimación concreta y uniforme.

**Lectura adicional del código público [hecho documental].** Se leyeron `src/Main.lean`, `src/Statements.lean` y la cadena breve `MainTheoremsAll` → `MainTheoremsCritical` → `MainTheorems` de la rama `original`, además de `Challenge.lean` en `main`. Las exportaciones visibles son `rhOfAll`, `allChi`, `critAll` y sus envoltorios de RH/GRH. El `sorry` de Challenge es deliberado: allí se enuncian los retos que Comparator compara contra la solución, no las pruebas de Main. No se confunde con un hueco detectado en la prueba certificada.

Referencias de lectura: [Main original](https://raw.githubusercontent.com/leanprover-community/riemann-hypothesis/original/src/Main.lean), [Statements original](https://raw.githubusercontent.com/leanprover-community/riemann-hypothesis/original/src/Statements.lean), [MainTheorems original](https://raw.githubusercontent.com/leanprover-community/riemann-hypothesis/original/src/Supplements/Kernels/MainTheorems.lean), [Challenge](https://raw.githubusercontent.com/leanprover-community/riemann-hypothesis/main/src/Challenge.lean). Es una inspección parcial de la interfaz, no una búsqueda exhaustiva en el desarrollo. El conector GitHub no estaba autenticado; se utilizó lectura pública, sin ejecutar Lean ni descargar sus artefactos grandes.

## 4. Contraejemplo denso a un atajo de Fourier [P, prueba completa]

**Afirmación refutada [FALSO].** «Una densidad positiva, la simetría de frecuencias propia de soportes impares y cotas de raíz cuadrada en todas las demás frecuencias garantizan una representación aditiva para cada par».

La siguiente familia permite ver la pérdida de fase. Es un ejemplo sobre conjuntos sustitutos, **no sobre el conjunto de primos**.

Sea ℓ un primo con ℓ≡3 mod 4 y χ el carácter cuadrático en 𝔽_ℓ, con χ(0)=0. Definimos

    A = {x ∈ 𝔽_ℓ : χ(x)=−1},  f=1_A,
    F(k)=Σ_x f(x) exp(2πikx/ℓ),
    G=Σ_x χ(x) exp(2πix/ℓ).

**a) Tamaño y ausencia de suma cero.** |A|=(ℓ−1)/2 y χ(−1)=−1. Si x∈A, entonces χ(−x)=1, de modo que −x∉A. Así (f*f)(0)=0.

**b) Magnitud de la suma de Gauss.** La conjugación y el cambio x↦−x dan Ḡ=−G. Además, sustituyendo x=ty en |G|² y sumando primero sobre y≠0,

    |G|² = Σ_t χ(t) Σ_{y≠0} exp(2πi(t−1)y/ℓ)
          = (ℓ−1) − Σ_{t≠1}χ(t) = ℓ.

Por tanto G²=−ℓ.

**c) Todas las frecuencias no triviales son pequeñas.** Para k≠0,

    F(k) = (−1−χ(k)G)/2,
    |F(k)|² = (ℓ+1)/4.

La primera identidad usa Σ_{x≠0}e(kx/ℓ)=−1 y el cambio x↦k⁻¹x en la suma del carácter.

**d) El término principal se cancela exactamente.** Como Σ_{k≠0}χ(k)=0,

    Σ_{k≠0} F(k)² = −(ℓ−1)²/4 = −F(0)²,
    (f*f)(0) = ℓ⁻¹ Σ_k F(k)² = 0.

Cada frecuencia no trivial tiene tamaño O(√ℓ), pero sus cuadrados tienen la fase necesaria para cancelar el término cero. Usar |F|² en lugar de F² habría contado otra cosa.

**e) Versión con soporte impar y objetivo par.** Levantamos A a ℤ/(2ℓ) mediante

    B={b mod 2ℓ : b impar y χ(b mod ℓ)=−1}.

Por CRT hay un único b impar por cada residuo x mod ℓ. Así |B|=(ℓ−1)/2, y B no contiene a la vez b y −b. Con representantes 1≤b≤2ℓ−1, ningún par de elementos de B suma 2ℓ.

Su transformada cumple F_B(k+ℓ)=−F_B(k). Los modos 0 y ℓ son los de paridad; los demás son, por CRT, un signo por una frecuencia no trivial de F y tienen magnitud √(ℓ+1)/2. Por tanto el ejemplo también respeta la restricción de soporte impar.

**Límite del contraejemplo.** No reproduce la estructura multiplicativa ni todas las condiciones locales de los primos. Descarta una deducción basada solamente en las propiedades enumeradas, no Goldbach ni todas las estrategias de Fourier. Es una construcción clásica con caracteres, sin reivindicación de originalidad.

## 5. Comprobación finita [N]

`numerico/puente_paridad_check.py` contrasta por dos cálculos distintos:

- Pellet: factorización por divisores frente a determinantes de Sylvester, para los 1092 polinomios mónicos de grados 1–6 sobre 𝔽₃.
- 81 casos de fibras r+s³: derivada fija; las sumas de μ por fibra incluyen ±3. No se presupone cancelación en cada fibra.
- 16 ejemplos de la construcción densa hasta ℓ=127, con recuento exacto de sumas y DFT aproximada como comprobación adicional.
- Las identidades enteras y coeficientes racionales de 04.

𝔽₃ no cumple la hipótesis del teorema de primos gemelos: estos cálculos comprueban las identidades algebraicas, **no** ese teorema asintótico. El JSON registra hash del script y método; los detalles exactos del resultado están en ESTADO, en la entrada de Codex.

## 6. Próximo objetivo matemático delimitado [P]

Buscar en el método automorfo una estimación de la suma primitiva U_N(D) definida en 04, o directamente una cota superior de su cola con signo. Exigir uniformidad en N, c y el extremo del intervalo. Un resultado para un desplazamiento fijo, una media en N o solo una norma L² no satisface automáticamente el objetivo.

**Advertencia de régimen, comprobada en Cor. 6.1 [T / P para la comparación].** Allí m≤αd, con α<1: el grado del coeficiente de la forma ponderada por Λ debe ser menor que el de la variable. En (4.10), D=N^(1/3) permite c de tamaño N^(2/3) frente a variables de tamaño N^(1/3). No satisface ese régimen. Ni siquiera trasladar literalmente el corolario resolvería esa porción de la cola. El corolario impone además condiciones a los grados de los desplazamientos; el análogo Goldbach requiere la modificación de formas indicada en Remark 1.2, no una sustitución automática en el enunciado de gemelos.

En §5.3 de 04 se prueba el balance: un nivel D=N^(1/2+σ) deja c≤N^(1/2−σ) y permite comparar tamaños con α=(1−2σ)/(1+2σ)<1. Pero tratar los divisores pequeños hasta ese D exige distribución adicional más allá de 1/2: la cota AP individual de GRH ya no basta. Es necesario mantener **ambas** necesidades en el puente, no únicamente la mención de discriminantes o de GRH.

Las proposiciones de 04 ya separan: identidades exactas, contribuciones controlables mediante GRH, términos no coprimos pequeños y la estimación aún abierta. **La existencia de esa estimación no se afirma.**

**Refinamiento obtenido al continuar la derivación.** §4.3 de 04 sustituye GRH por Bombieri–Vinogradov y la cota incondicional de Mertens para obtener la reducción asintótica con error de cualquier potencia logarítmica. Así la separación del obstáculo no depende de aceptar la noticia de RH. GRH proporciona una variante con error en potencias, pero la correlación pendiente permanece. No se trata de un nuevo teorema sobre primos: es una aplicación documentada de herramientas clásicas al mapa de trabajo del proyecto.
