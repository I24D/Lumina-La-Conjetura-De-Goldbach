# Medición decisiva del exponente de conteo de 𝒟_S (rama conmutada de D-NP)

**Agente:** Claude Code (Anthropic, modelo `claude-opus-5-5`), brazo verificador de la línea Goldbach.
**Fecha:** 2026-11-28.
**Encargo:** de Muse, por orden de Dal (encargo 1 de 4, prioridad máxima). Medir el exponente de conteo de 𝒟_S con los soportes más grandes que aguante la PC: exponentes locales por ventanas, ajuste global con intervalo honesto y veredicto frente a 0,420232.
**Método:**
1. Definiciones leídas de los docs 86 §1, 94 §0, 101 §2.4, 106 (106.2), 115 §7, 120b §1, 121 §2 y 121b §3.
2. Reimplementación propia con un algoritmo de conteo nuevo y exacto (§1), sin importar código de nadie.
3. Script: `numerico/verificacion_conteo_ds_claude_code.py`. Salida: `numerico/resultados/20261128_verificacion_conteo_ds_claude_code.txt`. **6/6 PASS**, 280 s.

**Estado:** PROPUESTA hasta el contraste de Muse.
**Etiquetas:** [T] teorema, [M] medición, [P] pendiente. Una medición no demuestra un exponente asintótico; el §4 da además el teorema.

**Goldbach sigue abierta.**

---

## §0. Veredicto en una línea

En la ventana superior medida, z = 307..709 (X = z¹⁰ hasta 3,2·10²⁸, |𝒟_S| hasta 8,6·10¹⁰), el exponente de conteo de 𝒟_S es **α = 0,4659, con IC95% [0,4649; 0,4669] entero sobre 0,420232** [M]. Asintóticamente es **1/2** [T, §4].

La condición "exponente de conteo < 0,420232" de la ruta Baier–Zhao **no se cumple**: ni en la medición ni en el límite. El 0,4211 que se citaba es el α de |H| (no de 𝒟_S) en la ventana pequeña z = 11..37, donde todavía hay saturación preasintótica.

---

## §1. Qué se cuenta, exactamente

**Definiciones** (docs 115 §7, 120b §1, 121 §2):
- Ω_{qm}^C = Θ̃_m(q)/(q−1), con Θ̃_m(q) = Σ_{h∈V⁺(D/q)∖{1}, m|h} (2τ(h)−1)/φ(h).
- 𝒟_S = {d = qm : Ω_d^C > 0}.
- W_d = Ω_d^C·φ(d).
- Familia del doc 101: A = 3, D = z⁵, y = round(z^{10/3}), V⁺ de Nathanson con nivel ⌊D/q⌋.

**Consecuencia** [T]. Por el cierre por divisores de V⁺ (demostrado en `verificacion-101-claude-code.md` §2.4), Θ̃_m(q) > 0 ⟺ m ∈ V⁺(L_q), excepto m = 1 cuando V⁺(L_q) = {1}. Por tanto:
- |𝒟_S| = Σ_q |V⁺(L_q)| − #{q : V⁺(L_q) = {1}}.
- Σ_{d∈𝒟_S} W_d = Σ_q Σ_{h∈V⁺(L_q)∖{1}} (2τ(h)−1)·h/φ(h).

La última fórmula sale de Σ_{m|h} φ(m) = h.

**Algoritmo nuevo, exacto.** Sea T(d) = máx_{m impar} p_m²·(p₁⋯p_m). Para d > 1 se cumple d ∈ V⁺(L) ⟺ L > T(d). Luego

Σ_q |V⁺(L_q)| = Σ_{d∈V⁺(L_max)} #{q ∈ 𝒮 : q ≤ ⌊D/(T(d)+1)⌋}.

Basta enumerar V⁺ de **un solo** nivel (11 millones de elementos en z = 709) y contar primos con una criba segmentada propia, en vez de recorrer los 8,6·10¹⁰ elementos de 𝒟_S.

**Validaciones** (V1, V2):
- |H| coincide exactamente con el doc 101 (z = 11..29) y con mi enumeración del doc 101 (z = 59..151).
- |𝒟_S| y ΣW_d coinciden exactamente con la **definición literal** (106.2): enumerando ℋ_C y calculando cada Ω_d^C, en z = 11, 13, 17, 19.

---

## §2. Los datos

Tabla resumida; la salida tiene los 41 puntos z = 11..709:

| z | X = z¹⁰ | primos pequeños | #sprimes | \|𝒟_S\| | ΣW_d |
|---|---|---|---|---|---|
| 37 | 4,8·10¹⁵ | 9 | 15 386 | 132 588 | 1,08·10⁶ |
| 151 | 6,2·10²¹ | 33 | 1 171 294 | 6,80·10⁷ | 1,22·10⁹ |
| 401 | 1,1·10²⁶ | 76 | 25 130 477 | 6,02·10⁹ | 1,65·10¹¹ |
| 709 | 3,2·10²⁸ | 124 | 152 595 370 | 8,61·10¹⁰ | 2,88·10¹² |

**Exponentes locales de |𝒟_S|** (pendiente entre z consecutivos respecto de ln X) [M]:
- Para z ≤ 37 oscilan entre 0,41 y 0,51.
- Desde z = 151 están **todos** entre 0,445 y 0,472 (E1), con tendencia ascendente: 0,454 (151), 0,462 (199), 0,455 (307), 0,466 (401), 0,470 (607), 0,472 (659), 0,467 (709).
- Los de ΣW_d rondan **0,50** desde z ≈ 60.

**Ajustes por ventanas** (mínimos cuadrados sobre ln X):

| ventana | α(\|H\|) | α(\|𝒟_S\|) | α(ΣW_d) | veredicto para 𝒟_S |
|---|---|---|---|---|
| z = 11..37 (ventana del doc 102) | **0,4211** [0,4185; 0,4238] | 0,4444 [0,4268; 0,4619] | 0,5177 [0,5018; 0,5336] | entero sobre |
| z = 11..59 | 0,4225 [0,4208; 0,4241] | 0,4359 [0,4276; 0,4441] | 0,5078 | entero sobre |
| z = 59..151 | 0,4486 [0,4472; 0,4500] | 0,4487 [0,4473; 0,4501] | 0,5008 | entero sobre |
| z = 151..401 | 0,4582 [0,4566; 0,4598] | 0,4583 [0,4566; 0,4599] | 0,5017 | entero sobre |
| **z = 307..709** | 0,4659 [0,4649; 0,4668] | **0,4659 [0,4649; 0,4669]** | 0,5020 [0,5013; 0,5028] | **entero sobre** |
| z = 11..709 (global) | 0,4461 [0,4432; 0,4491] | 0,4484 [0,4461; 0,4508] | 0,5018 | entero sobre |

**Sobre la honestidad del intervalo.** Los puntos son deterministas, no una muestra. El IC95% de mínimos cuadrados mide la dispersión del ajuste en la ventana, no una incertidumbre sobre el límite. La evidencia honesta tiene tres piezas que coinciden:
1. Todos los exponentes locales desde z = 151 superan 0,4452.
2. Esos exponentes crecen con z.
3. El teorema del §4 fija el límite en 1/2.

---

## §3. De dónde salía el 0,4211

- Es el α de **|H|**, el soporte entero (parte lisa más rama conmutada), en z = 11..37 (doc 102 §3). Lo reproduzco exactamente (F1) y efectivamente queda **a caballo** de 0,420232.
- Para **𝒟_S** en esa misma ventana, el α es 0,4444 [0,4268; 0,4619], ya **entero sobre** el umbral. |𝒟_S| no es |H|: en z pequeños quitar los q con V⁺(L_q) = {1} pesa mucho (en z = 11, |𝒟_S| = 578 frente a |H| = 802).
- Aun así, la ventana pequeña es preasintótica. Con z ≤ 37 hay a lo sumo 9 primos pequeños, así que |V⁺(L)| ≤ 2⁹ en cada nivel, y el exponente medido ahí está deprimido por saturación.

---

## §4. El límite asintótico [T]

**Teorema.** En la familia del doc 94 (X = N/A, z = X^{1/10}, y = X^{1/3}, D = X^{1/2}/(log X)^{C0}, cualquier C0 ≥ 0 fijo), para todo ε > 0 y X suficientemente grande:

X^{1/2−ε} ≤ |𝒟_S| ≤ 2X^{1/2},  y  |𝒟_S| ≤ Σ_{d∈𝒟_S} W_d ≤ |𝒟_S|·X^{o(1)}.

Por tanto, el exponente de conteo de 𝒟_S, tanto sin pesos como con los pesos W_d, es **1/2**.

*Prueba.*
- **Cota inferior.** Es la del Teorema 101-CC (`verificacion-101-claude-code.md` §2.7): por el Lema A, |V⁺(L)| ≥ L^{1−4δ−o(1)} para cada nivel, y entonces Σ_q |V⁺(L_q)| ≥ X^{1/2−2δ/3−o(1)}. Para X grande ningún V⁺(L_q) es {1}, porque L_q ≥ L_min ≍ X^{1/6} > 5³.
- **Cota superior.** |𝒟_S| ≤ Σ_q L_q ≤ 2D.
- **Pesos.** Cada término de W cumple 3 ≤ (2τ(h)−1)·h/φ(h) ≤ 2τ(h)·O(log log X), y τ(h) ≤ X^{o(1)}. ∎

**Comparación con el umbral, en las dos convenciones posibles:**
- **Respecto de X** (la convención del doc 101 y del 0,4211): el límite es α_X = 1/2 > 0,420232.
- **Respecto de N con A = N^{c\*}** (la del objetivo de Chen del doc 100, X = N^{1−c\*}): α_N = (1−c\*)·α_X.
  - En la ventana superior, α_N = 0,4473 [0,4464; 0,4482], entero sobre.
  - En la ventana del doc 102, 0,4267 [0,4098; 0,4435], a caballo.
  - En el límite, α_N = 0,480058 > 0,420232.

**El veredicto asintótico no depende de la convención.** El de ventana sí, pero solo en la ventana pequeña.

---

## §5. Qué significa para la ruta

1. La ruta Baier–Zhao sobre 𝒟_S, **tal como se formuló** (exige un exponente de conteo < 0,420232), **no es viable** [M + T]. El conjunto tiene exponente 1/2 y la ventana grande ya lo muestra por encima de 0,465. No es una diferencia de milésimas: en el límite es de 0,08.
2. Si la versión de la ruta usa el conteo **pesado** ΣW_d (el peso Ω_d^C·φ(d) del doc 121), el exponente es ≈ 0,50 ya en todas las ventanas medidas.
3. **Lo que NO dice este resultado.** No dice que D-NP sea falso ni que el momento de caracteres sobre 𝒟_S sea grande. Dice que cualquier cota que lo controle **contando módulos** (tipo criba grande o Baier–Zhao con el tamaño de 𝒟_S) necesita ahorrar al menos X^{0,08} por otra vía.
4. **Precisión de alcance.** La condición "0,420232 exigido" viene del doc 121 §3 y de tu mensaje. El texto original de la ruta Baier–Zhao (respuesta de Lumina Chat al doc 119) no está en el repositorio, así que no he podido auditar de dónde sale exactamente ese umbral. Mido las dos cantidades candidatas: el conteo puro y el pesado.

---

## §6. Cómo reproducir

```
cd "C:\Users\dal_n\OneDrive\conjetura de golbach Lumina\numerico"
python -B verificacion_conteo_ds_claude_code.py
```

Solo biblioteca estándar. Tarda unos 5 minutos: el último punto, z = 709, criba los primos hasta 3,2·10⁹ en segmentos y enumera 1,1·10⁷ elementos de V⁺. Memoria de alrededor de 1 GB.

*Goldbach fuerte sigue abierta.*
