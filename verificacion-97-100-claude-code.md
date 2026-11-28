# Verificación de los puntos [P-rev] del doc 97 §8 y del doc 100 §7

**Agente:** Claude Code (Anthropic, modelo `claude-opus-5-5`), tercer brazo: verificación y cómputo en la PC.
**Fecha:** 2026-11-28.
**Encargo:** de Muse, por orden de Dal. Contrastar cada [P-rev] contra los docs 94, 96, 98 y 99 íntegros y los scripts 94d, 94e, 97 y 99.
**Método:** (1) lectura íntegra de las fuentes; (2) re-ejecución sin cambios de los cuatro scripts de Lumina Chat; (3) reimplementación propia desde las definiciones escritas, sin importar código de nadie: `numerico/verificacion_97_100_claude_code.py`, salida en `numerico/resultados/20261128_verificacion_97_100_claude_code.txt`, **26/26 PASS**, 0,9 s.
**Estado:** PROPUESTA. Muse debe contrastarla antes de aceptarla.
**Publicación:** subido a GitHub el 2026-11-28 por orden expresa de Dal, para que Muse y Lumina Chat puedan leerlo. Publicarlo no lo convierte en verificado: sigue siendo PROPUESTA hasta el contraste de Muse.
**Etiquetas:** [T] teorema (prueba escrita o identidad exacta comprobada), [M] medición (no es prueba), [P] pendiente. Veredicto por punto: **confirmado / refutado / no verificable**.

**Goldbach sigue abierta.** Nada de este informe la demuestra ni la acerca más allá de precisar qué falta.

---

## §0. Fuentes y advertencias (leer primero)

**Dónde está cada archivo.** La carpeta de OneDrive y el repositorio público han divergido, y ninguno de los dos tiene todo:

| Fuente | Dónde la encontré | SHA-256 (16) |
|---|---|---|
| doc 94, 96, 98, 99 | OneDrive `matematicas/` (no están en GitHub) | 94 `a016c73545c95c04` · 96 `625714891ea26693` · 98 `ecfc3da6b0f9413c` · 99 `18967167f270aba7` |
| scripts 94, 94d, 94e, 97, 99 | OneDrive `numerico/` (no están en GitHub) | 94 `522de0f145858048` · 94d `de31002c40be812b` · 94e `0b92899d66a23fb1` · 97 `e0dd2eb8382291e6` · 99 `d97a402837cc9db1` |
| doc 95, 97, 100 | **Solo en GitHub** (clon anónimo, `3a82328`, 2026-11-28 08:06). **No están en OneDrive.** | 95 `61a9b229970a2ee1` · 97 `66e6deefde26ec89` · 100 `28f74e2f2518d9c5` |
| clon local de Lumina Chat | `C:\Users\dal_n\dev\Lumina-La-Conjetura-De-Goldbach`, rama `main` **7 commits** por delante de su `origin/main` (incluye `e51a4c5`, `4147a6a` y `d27bdb3`, este de 2026-11-28 06:06) | sus docs 94/96/98/99 y scripts son **idénticos byte a byte** a los de OneDrive |

En el historial git del clon no hay ninguna versión del doc 99 ni del script 99 con N=999983 (búsqueda `git log -S999983`).

**Reproducibilidad.** Los scripts 94d, 94e, 97 y 99, ejecutados sin cambios (Python 3.14.10), reproducen **byte a byte** sus salidas guardadas en `numerico/resultados/20261127_*.txt` (diferencias solo de BOM y CRLF). [M]

**Independencia: dos niveles distintos.**
- *De método:* mi código no importa ni ejecuta el de Lumina ni el de Muse. Uso algoritmos distintos: enumeración recursiva de V⁺, aritmética exacta con `Fraction`, una ruta para S_r **sin caracteres** y un test de primitividad propio.
- *De contrato:* **esto no es revisión independiente** según la regla 4 del README y el punto 5 del contrato. Comparto usuario de Windows con las corridas de Lumina en esta PC. Lo que aporto es control de reproducibilidad y comprobación de identidades, no una revisión independiente.

**Memoria compartida.** El MCP `lumina_shared_memory` estaba caído en esta sesión, así que no registré la propuesta allí (regla 6). Tampoco anoté nada en `ESTADO.md`, para no adelantarme a tu contraste.

---

## §1. Resumen

| Punto | Veredicto | Evidencia clave |
|---|---|---|
| 97-§8.1 re-verificación global | **Confirmado** el núcleo; **refutadas** tres afirmaciones accesorias del doc 97 | ver 97-§8.3e, 8.3b y 8.4b |
| 97-§8.2 ¿el doc 96 ya tenía Δ_r? | **Confirmado que NO lo tenía** | doc 96 §5 l.135-142; la errata §7 l.170-187 es posterior |
| 97-§8.3 L_s y rangos de 94d | L_s **recuperado**; rango Q_s·log z **confirmado**; mínimo de Q_s/P **refutado** (es 1,00186); hipótesis del factor local **refutada** | chequeos A1, A5-A8 |
| 97-§8.4 regla de atribución 94e | Regla **confirmada**; la hipótesis "a = A" del doc 97 **refutada** | chequeos B2, B3 |
| 97-§8.5 desigualdad de Cauchy del doc 94 | **Confirmada** (forma exacta y álgebra) | doc 94 §6 l.279-296 |
| 97-§8.6 los seis [P-rev] del doc 95 | 4 **confirmados**, 1 **confirmado parcialmente**, 1 **refutado en su letra** | §2.6 |
| 97-§8.7 enunciado de sup Q_s en el doc 96 | **Confirmado**: límite fijo y ya normalizado por P_{A,N} | doc 96 l.53-64 |
| 97-§8.8 enviar el veredicto a Lumina Chat | **No verificable** (es una acción, y la hace Muse) | — |
| 100-§7.1 modelo exacto del doc 99 | **Resuelto.** Las cifras del chat salen del mismo script con **N=999983**; la "hipótesis z=3 no confirmada" queda **refutada** | chequeos D1, D2, D3 |
| 100-§7.2 los 16 conductores del script 97 | **Resuelto**; las coincidencias se deben a r=7 y r=9, comunes a ambos conjuntos | chequeo C6; script 97 l.10 |
| 100-§7.3 prueba de Δ_r del doc 98 y la condición p∥r | Prueba **confirmada** [T]; p∥r explícita; el "espurio" del doc 100 §2.2 queda **refutado tal como se plantea** | chequeos C4, C5, C9, C10 |
| 100-§7.4 ¿cambia alguna cifra de 94d/94e por Δ_r? | **Confirmado que no** (estructural) | scripts 94d/94e |
| 100-§7.5 | = 97-§8 | §2 |
| 100-§7.6 | = 97-§8.6 | §2.6 |
| 100-§7.7 experimento E100.1 | **No verificable** con las fuentes del encargo. Muse ya lo ejecutó en el doc 101; no lo he re-verificado [P] | — |

---

## §2. Doc 97 §8, punto por punto

### 97-§8.1 Re-verificar el doc 97 contra el doc 96 íntegro y los scripts 94d/94e

**Confirmado en lo sustancial.**
- La identidad diagonal, el producto completo, la cancelación del término principal y la necesidad de Δ_r se sostienen frente a las fuentes íntegras.
- Tres afirmaciones accesorias del doc 97 **no sobreviven**: la lectura "r > L_s" (8.3b), la explicación del extremo 2,829 por el factor local de N (8.3e) y la regla "a = A" con su cifra de ~16× (8.4b).

El detalle va en cada subpunto.

### 97-§8.2 ¿Contenía el doc 96 el término Δ_r?

**Confirmado que NO lo contenía.**
- El doc 96 §5 (l.135-142) enuncia en un recuadro, como "identidad exacta", S_r = Σ_{d|r} μ(r/d)φ(d)E(N;d,N), sin Δ_r.
- La errata §7 (l.170-187) se añadió después (cita los docs 97 y 98) y retira esa identidad.
- Por tanto, la corrección del doc 97 §3 era **al documento mismo, no solo a la transcripción del chat**.

Numéricamente [M], la tabla del doc 97 §3.2 se reproduce entera, incluido r=6 (chequeo C7, desviación máxima 4,9·10⁻⁵ por redondeo a 4 decimales).

### 97-§8.3 Definición de L_s, rangos de 94d y cociente Q_s/P

**(a) Definición de L_s: recuperada (confirmado).**
- Doc 96 l.17: L_s = D_A/s.
- En el proxy del script 94d (l.37-40): L = ⌊D/s⌋, con D = ⌊√X⌋, X = ⌊N/A⌋ y z = ⌊X^{1/10}⌋, **sin** la pérdida (log X)^{C_0} de D_A.
- La partición del script (l.82-85) va por el **índice divisor r** de la identidad positiva, no por el conductor ah (lo dice bien el doc 98 §5, l.155-170).

**(b) Los umbrales son r > L_s/z y r > √L_s. La lectura "r > L_s" del doc 97 queda refutada.**
- El doc 96 l.84 y el script (l.85) usan **√L_s**.
- El doc 97 §4.2 (l.316-329) comparó su masa en r > z² con una supuesta "masa en r > L_s de 18,15 %". Esa cantidad **es idénticamente 0**: todo r de la diagonal divide a algún d ∈ V⁺(L_s) con d < L_s, así que r < L_s. [T] Medido además en los 20 escenarios (chequeo A8).
- Valores correctos [M]: máximo en r > L_s/z = **9,8375 %** (e=60, u≈2); máximo en r > √L_s = **18,1459 %** (e=60, u≈3).

**(c) Rango de Q_s·log z ∈ [1,55603; 2,82873]: confirmado.** [M] Lo reproduzco escenario por escenario en los 20 casos, con Q_s **racional exacto** por la diagonal y contrastado con la doble suma de Gram (diferencia ≤ 5,6·10⁻¹⁷) (chequeos A1, A3).

**(d) Rango de Q_s/P que da el doc 96, [1,00266; 1,70054]: mínimo refutado.**
- El mínimo real es **1,00186** (e=50, u≈4); el máximo 1,70054 sí se confirma (chequeo A5).
- Es una errata menor del doc 96 l.81, que tomó el valor de e=40. El doc 94 l.391 ya daba "z=28: 1.0019" correctamente. El "[1,003; 1,701]" del doc 97 hereda esa errata.

**(e) La hipótesis del doc 97 (§4.3, l.336-342, y §5.3.2, l.416-422) de que "el extremo 2,829 lo explica el factor local de N" queda refutada.** [M]
- En los 20 escenarios de 94d, N = 2^e y A = 3: los primos excluidos son siempre {2, 3}, así que el factor local es **el mismo en todos** (chequeo A6).
- Descomposición Q_s·log z = (P·log z) × (Q_s/P):
  - P·log z ∈ [1,552; 1,688] solo varía por la convergencia de Mertens en z;
  - Q_s/P ∈ [1,0019; 1,7005] **lleva todo el rango**, ordenado por u:

| u objetivo | Q_s/P | Q_s·log z | tamaño de V⁺(L_s) |
|---|---|---|---|
| 5/3 | 1,534 – 1,701 | 2,434 – 2,829 | 1 – 4 |
| 2 | 1,334 – 1,454 | 2,096 – 2,454 | 2 – 11 |
| 3 | 1,032 – 1,101 | 1,602 – 1,823 | 11 – 216 |
| 4 | 1,002 – 1,007 | 1,556 – 1,700 | 14 – 995 |

- El extremo 2,82873 está en e=60, u≈5/3 (chequeo A7). La cercanía con el 2,874 del modelo E1 del doc 97 es **casual**.

**Observación nueva [M], para quien use este perfil.** En u≈5/3 el proxy es **degenerado**: V⁺(L_s) tiene 1, 2 o 4 elementos. En e=40, V⁺ = {1}, de modo que Q_s = 1 y Q_s/P = 1/P exactamente. El extremo superior del rango mide la truncación extrema de un proxy pequeño (D = √X, z ≤ 57). No dice nada sobre la constante asintótica.

### 97-§8.4 Regla exacta de atribución de masa en 94e ("ah ≤ R0") y tabla κ

**(a) Regla: confirmada.**
- Script 94e l.41-48 y doc 98 l.174-186: para cada a | A y cada h, la masa es m(a,h) = φ(a)φ(h)|K_A(h)|²; se cuenta como "baja" si ah ≤ R0 = (log N)^κ.
- La masa total es A·G_A. Lo verifiqué: (a,h) = 1 y Σ_{a|A} φ(a) = A. [T]

**(b) La hipótesis del doc 97 §4.1 (l.300-308), "atribuye la masa a r = A·h, con el mayor a", queda refutada.**
- La regla reparte la masa entre **todos** los a | A con peso φ(a).
- Para A = 3 es exactamente la mezcla **(1/3)·[corte h ≤ R0] + (2/3)·[corte 3h ≤ R0]**, como identidad racional exacta (chequeo B3).
- Las dos "lecturas" del doc 97 son los extremos de esa mezcla, no la regla.

**(c) Tabla κ: confirmada** con el soporte reconstruido desde cero (1444 pares (s,d), 1451 partes nuevas, G_A = 0,450566952208, chequeos B1 y B2; desviación máxima 4,9·10⁻⁷ puntos). [M]

| κ | masa baja (regla exacta) | masa alta | mejora de amplitud √(total/alta) | alta con lectura "h ≤ R0" | alta con lectura "3h ≤ R0" |
|---|---|---|---|---|---|
| 1,0 | 49,169430 % | 50,83 % | ×1,40 | 46,05 % | 53,22 % |
| 2,0 | 73,278523 % | 26,72 % | ×1,93 | 20,59 % | 29,79 % |
| 2,5 | 87,245823 % | 12,75 % | ×2,80 | 5,94 % | 16,16 % |
| 2,7 | 92,878861 % | 7,12 % | **×3,75** (energía ×14,04) | 1,06 % (×9,71) | 10,15 % (×3,14) |
| 2,8 | 95,249831 % | 4,75 % | **×4,59** (energía ×21,05) | 0,50 % (×14,10) | 6,87 % (×3,81) |
| 3,0 | 99,111943 % | 0,89 % | ×10,61 | 0,13 % | 1,26 % |

**Consecuencia sobre el factor de Cauchy (doc 97 §4.4, [P-rev] "su valor en los parámetros del doc 96").**
- En los parámetros de Lumina la mejora de amplitud es **×3,75 (κ=2,7) y ×4,59 (κ=2,8)**. El "4–5×" del doc 96 es aproximado: en κ=2,7 sale por debajo de 4.
- El "~16×" del doc 97 viene de su salida H4 (`salida_doc97_local.txt` l.56 y l.61-62). Usa la **lectura "h ≤ R0"** sobre sus propios modelos (f = 0,36–0,41 %, 1/√f = 15,7–16,7 en amplitud).
- En el propio proxy de Lumina, esa misma lectura ya da ×9,71 y ×14,10. **La mayor parte de la diferencia entre 4–5× y 16× la explica la regla de atribución**; el resto viene de los modelos.
- Esto corrige la reconciliación del doc 98 l.188. Ahí se sugiere que el 16 "puede referirse al factor de energía". No es así: el doc 97 lo da explícitamente como 1/√f (amplitud), y su factor de energía es 246–278×.

### 97-§8.5 Forma exacta de la desigualdad de Cauchy del doc 94 y del segundo factor

**Confirmada.** Doc 94 §6, l.279-287:

|R^prim_{>R0}|² ≤ ( Σ_{a|A, h∈H_A, ah>R0} φ(ah)|K_A(h)|² ) · ( Σ_{a|A, h∈H_A, ah>R0} |S_{a,h}(N)|²/φ(ah) ),

con el primer factor ≤ A·G_A y el segundo = 𝓑_A(N). El condicional está en l.288-294 (G_A ≪ 1/(φ(A)² log N), 𝓑_A ≪ N²/(A log³N)).

Comprobaciones mías [T]:
1. Es Cauchy–Schwarz con pesos φ(ah) sobre R = Σ K_A(h)S_{a,h}.
2. La cota del primer factor por A·G_A usa que (a,h) = 1 y Σ_{a|A} φ(a) = A.
3. La aritmética del condicional da R ≪ N/(φ(A) log²N) ≍_k N/(A log²N), porque A/φ(A) = O_k(1).

El doc 96 §4 (l.105-115) la refina: se usa G_{>R0} en lugar de A·G_A para aprovechar el corte. Tras la errata, el segundo factor debe formarse con S_{a,h} = M_{ah} + Δ_{ah} (doc 98 §4, l.134-151). **Lo que falta sigue siendo la cota de 𝓑** (doc 94 l.304, [N]).

### 97-§8.6 Los seis [P-rev] del doc 95, contra el doc 94 íntegro

El doc 95 no los numera como lista de seis. Los reconstruyo a partir de sus marcas [P-rev] (§2.1, §5 (d2)) y de la lista del §7.4. **Muse: confirma que esta correspondencia es la tuya.**

| # | [P-rev] del doc 95 | Veredicto | Evidencia en el doc 94 |
|---|---|---|---|
| 1 | Declarar supp b_A ⊆ {(m,N)=1}, régimen de la expansión por caracteres | **Confirmado parcialmente.** (m,N)=1 sí está declarado. (A,N)=1, que también hace falta porque la expansión va al módulo Am, **no** se declara en §0–§1; solo aparece en el experimento 94b (l.381). [P]: añadirlo al enunciado | l.18: "(sd,AN)=1" |
| 2 | Error imprimitivo (carácter inducido frente a primitivo), escrito y acotado | **Confirmado**: está, con O(D_A log²N) | §1, l.85-102 |
| 3 | Pieza del carácter principal | **Confirmado**: P_q = −(1/φ(q))Σ_{p|q}Σ_{p^j≤N} log p, acotada junto con la anterior. Coincide con el −Θ(d;N)/φ(d) de (98.1) | l.46-50 |
| 4 | Versión de límite fijo o variable de la segunda reducción | **Confirmado**: el doc 94 usa el **límite fijo** z_A. El contraejemplo del doc 95 afecta solo a la versión variable, que el doc 94 no usa (doc 96 l.10) | l.15 |
| 5 | ¿Contiene el doc 94 el puente G_A → R_A? (doc 95 §4: "no está en la derivación") | **Refutado en su letra, confirmado en su fondo.** La desigualdad formal sí está (97-§8.5). Lo que falta es la cota del segundo factor, que el doc 94 marca [N] | §6, l.279-296, l.304 |
| 6 | (d2) Restringir la estimación de bloque a ah > R0 y tratar el complemento | **Confirmado**: el bloque r ≤ R0 se trata con Siegel–Walfisz y el resto con Cauchy. El doc 94 l.144 avisa de que un h pequeño con a grande cae en el bloque grande. *Matiz:* el "~94 % de Q_s" del doc 95 (masa de h₀ ≤ 10–30) y el "92,9 % de A·G_A" del 94e son **magnitudes distintas**. El doc 96 l.126 las relaciona ("hace plausible el 94 %"), pero no son comparables directamente | §3, l.117-144; §6 |

### 97-§8.7 Enunciado exacto de sup Q_s en el doc 96

**Confirmado.** Doc 96 §2, l.53-63: Q_s ≪ P_{A,N}(z_A) = ∏_{p<z_A, p∤AN}(1 − 2/(p−1) + p/(p−1)²), con constante uniforme en u_s ∈ [5/3−o(1), 4−o(1)] (l.162).
- Es la versión de **límite fijo**; la variable solo se admite con la cola completa s' > s (l.10 y l.161).
- La constante **ya está normalizada por el producto local** (l.64). Va más allá: para N primoriales, el producto de referencia puede crecer como (log log N)/log z_A (l.72).
- Comprobé la identidad f(p)⁻¹ = (p−1)/(p−2)·(1 − 1/(p²−3p+3)) = (p−1)²/(p²−3p+3) de l.66-70. [T]

**Consecuencia:** la "corrección obligatoria 2" del doc 97 §7 (l.485-487, normalizar sup Q_s por el producto local) **ya estaba en el doc 96**. Sigue valiendo contra el doc 94 §8 (l.373-376), que la enuncia como sup Q_s ≪ 1/log z_A sin normalizar.

### 97-§8.8 Enviar el veredicto del doc 97 a Lumina Chat

**No verificable**: es una acción, no una afirmación, y la comunicación con Lumina Chat la llevas tú.

Hay evidencia indirecta de que el contenido llegó: el doc 96 §7 y el doc 98 l.4 citan el doc 97 y adoptan Δ_r. *Nota de procedencia:* los docs 96 y 98 atribuyen el doc 97 a "Dal", aunque su autora declarada es Muse.

---

## §3. Doc 100 §7, punto por punto

### 100-§7.1 Parámetros exactos del modelo del doc 99 y cifra B_M = 3,9222·10⁷

**Resuelto. Hay dos corridas del mismo script; el doc 100 auditó una y el disco contiene la otra.**

**La versión en disco** (doc 99 §2, l.38-53; script 99 l.10):
- Parámetros: N = 10⁶ (par), A = 3, X = 333333, D = 577, z = 3, y = 69. **No hay primos pequeños** (2 | N y 3 | A), así que d = 1.
- Soporte: 16 pares (s,d), |H| = 17, 33 conductores r = ah > 1. Peso 1/φ(r), suma sobre r > R0 = (ln N)^κ, r = 1 excluido (script l.194-196, l.235-240).
- κ = 1: **29 conductores altos, B_M = 1,981992·10⁷**, B_cross = 1,513556·10⁴, B_Δ = 278,437, cambio relativo 7,777·10⁻⁴.
- Tabla completa reproducida por la ruta sin caracteres, con desviación relativa ≤ 3,3·10⁻¹⁰ (chequeo D1). [M]

**Las cifras del chat que audita el doc 100** salen del **mismo script con N = 999983** (primo impar). Lo comprobé de dos formas: el script 99 de Lumina con solo esa línea cambiada, y mi reimplementación (chequeo D2). [M]
- Con N impar entra el primo pequeño 2: 33 pares, |H| = 35, **69 conductores en total, de los que 61 están sobre R0**.
- κ = 1: B_M = **3,922249·10⁷**, B_cross = **5,402292·10³**, B_Δ = **467,664**, cambio relativo **1,496579·10⁻⁴**, B_cross/B_Δ = 11,55.

**Consecuencias:**
- **"La hipótesis z=3 NO se confirma" (doc 100 §3.3.6, l.227-234): refutada.** El modelo del doc 99 **es** z = 3 (doc 99 l.43). El R3 del doc 100 no lo reproducía porque usaba y = 349 y D = 200000; Lumina usa y = 69 y D = 577, y con N impar también d = 2.
- **"69 conductores"** es el total; el segundo factor con κ = 1 se calcula sobre **61**.
- **Erratas de transcripción del doc 100:** B_Δ = 467,66, no 466,7, y por tanto el cociente es 11,55, no 11,57 (l.47, l.175, l.190).
- **Errata interna del doc 99, l.108:** "La cifra B_M ≈ 3,92·10⁷ es el tamaño del momento en este proxy" contradice su propia tabla (1,98·10⁷ con κ = 1). Es un resto de la corrida con N = 999983, que nunca se guardó en el repositorio. El doc 98 §7 (l.212-214) sí es coherente con N = 10⁶.

### 100-§7.2 Lista exacta de los 16 conductores del script 97 de Lumina

**Resuelto.** Script 97 l.10 y doc 98 l.124-126: **3, 5, 7, 11, 15, 21, 35, 55, 105, 385, 9, 25, 45, 75, 225, 315** (N = 999983).
- Comparte 11 conductores con el conjunto del doc 100 ({3, 5, 7, 9, 11, 15, 21, 25, 35, 45, 75}). Solo en Lumina: {55, 105, 225, 315, 385}. Solo en el doc 100: {13, 27, 49, 63, 99}.
- Reproduzco los 16 cocientes |Δ_r|/|S_r| con desviación ≤ 4,9·10⁻⁷ puntos (chequeo C6). [M]

**Sobre la "evidencia fuerte" del doc 100 (l.110-112):** el máximo |Δ|/|S| = 2,588591 % se alcanza en **r = 7** en los dos conjuntos, y el error relativo máximo de Lumina (9,7032·10⁻¹²) también está en **r = 9**, común a ambos. La coincidencia es la consecuencia trivial de calcular los mismos S_r con el mismo N, no una prueba aparte. Errata menor: el doc 100 transcribe 9,71·10⁻¹²; el archivo dice 9,7032·10⁻¹².

### 100-§7.3 La demostración de Δ_r del doc 98 y la condición p∥r

**Confirmada la prueba** [T]. Revisé cada paso:

1. **(98.1)**, l.55-62. Ortogonalidad módulo d sin el carácter principal (φ(d)E = Σ_{χ≠χ₀} χ̄(N)ψ(N,χ_d) − Θ(d;N)), más la relación ψ(N,χ_d) = ψ(N,χ*) − Σ_{p|d, p∤c} Σ_j χ*(p)^j log p. Ambas correctas con (N,d) = 1.
2. **(98.2)**, l.71-80. Inversión de Möbius: Σ_{d: c|d|r} μ(r/d) = [c = r].
3. **(98.3)**, l.86-106:
   - Término de inducción: Σ_{d: cp|d|r} μ(r/d) = [cp = r], y como p ∤ c eso equivale a p∥r con c = r/p.
   - Término principal: Σ_{d: p|d|r} μ(r/d) = [r = p].
4. **Equivalencias.** (98.3) es la misma fórmula que el Teorema 97.1 (p ∤ r/p ⟺ p∥r; el caso r primo es su convención r/p = 1). También coincide con la forma ortogonal del doc 100 cuando (p, r/p) = 1.

La condición p∥r está **explícita** en (98.3) (l.101) y en l.108.

**Numérico [M], en 22 conductores con N = 999983 (la unión de los conjuntos de Lumina, del doc 100 y del doc 97):**
- S_r por una ruta **sin caracteres**, S_r = Σ_{(n,r)=1} Λ(n)·Σ_{d|r} μ(r/d)φ(d)[n ≡ N mod d], coincide con la suma directa sobre caracteres primitivos de construcción propia: error relativo ≤ 5,4·10⁻¹² (C1).
- El número de primitivos dado por mi test directo coincide con la fórmula multiplicativa (C2).
- S_r = M_r + Δ_r(98.3) con error relativo ≤ 4,2·10⁻¹¹ (C4).
- Δ_r(98.3) evaluado con y sin caracteres: diferencia ≤ 2,8·10⁻¹³ (C5).
- r = 6: S_6 = 0 y Δ_6 = log 2 (C8).

**Precisión que refuta el doc 100 §2.2 tal como está planteado.**
- En (98.3), un término con p² | r llevaría χ*(p) con p | r/p, y **χ*(p) = 0** por definición. Esos términos se anulan solos. En la forma con caracteres, "olvidar" p∥r no cambia nada.
- La condición sí es esencial en la *derivación* (el coeficiente de Möbius) y en las formas sin caracteres, cuya identidad exige (p, r/p) = 1.
- La columna "espurio" del doc 100 (l.138-151) se reproduce exactamente **solo si se da a χ*(p) el valor 1** (C9, los 8 valores).
- Aplicando la forma ortogonal del propio doc 100 sin la condición se obtienen otros valores: coincide en r = 9, 25, 49, 99, 121, pero da **0** en r = 45, 63 y 75, donde el doc 100 dice +13,18, +13,18 y +12,88 (C10).
- Por tanto, "ignorar p∥r inventa correcciones de hasta el 3,6 %" no describe un error que (98.3) permita. Describe el uso de una identidad fuera de su hipótesis.

### 100-§7.4 ¿Cambia alguna cifra de 94d/94e tras la corrección Δ_r?

**Confirmado que no** [T por estructura]. Los scripts 94d y 94e solo calculan magnitudes del lado de los coeficientes (pesos V⁺, K_A(h), G_A, Q_s); no calculan ψ, E ni caracteres, y Δ_r vive en S_r, del lado de los caracteres. Lo dice también el doc 96 §7 l.183.

El doc 98 §5 no cambia ninguna cifra: aclara reglas y reinterpreta el "4–5×" como factor de amplitud. Mis partes A y B reproducen todas las cifras de 94d y 94e sin tocar Δ.

### 100-§7.5 y 100-§7.6

Coinciden con 97-§8 (todo el §2 de este informe) y con 97-§8.6 (§2.6).

### 100-§7.7 Experimento E100.1 (W_max(N), S_BV(N), N = 10⁶–10⁸)

**No verificable con las fuentes de este encargo**: no es un contraste contra los docs 94–99, sino un experimento nuevo. Muse ya lo ejecutó (doc 101 en el repositorio, criba segmentada hasta 10⁸, α_W = 0,4289). **No lo he re-verificado.** [P] Si lo quieres, es un buen encargo de cómputo para esta PC.

---

## §4. Erratas y avisos para los autores (sin tocar sus archivos)

Según la regla 3 del README, no modifiqué ningún archivo de otro agente. Propongo:

1. **Doc 96 l.81:** el mínimo de Q_s/P es 1,00186 (e=50, u≈4), no 1,00266.
2. **Doc 99 l.108:** "B_M ≈ 3,92·10⁷" corresponde a la corrida con N = 999983, no a la tabla del documento (1,98·10⁷). Conviene corregirlo o documentar las dos corridas.
3. **Doc 98 l.188:** el "~16×" del doc 97 es un factor de amplitud calculado con la lectura "h ≤ R0", no un factor de energía. La regla de atribución explica la mayor parte de la diferencia.
4. **Doc 97:**
   - §4.2: "r > L_s" debe ser "r > √L_s".
   - §4.3 y §5.3.2: hipótesis del factor local refutada.
   - §4.1: la regla no es "a = A".
   - §7, corrección 2: ya estaba en el doc 96.
   - Errata menor: §1 fila 7 dice "acuerdo 0,004–0,028", y §4.1 dice "0,004–0,041"; el real es 0,003–0,041.
5. **Doc 100:**
   - §3.3.6: el modelo del doc 99 **es** z = 3.
   - §2.2 y K2: el "espurio" supone χ*(p) = 1.
   - Transcripciones: B_Δ 467,66; cociente 11,55; 9,70·10⁻¹².
6. **Doc 94 §0–§1:** declarar (A,N) = 1 en el enunciado, no solo en el experimento 94b.
7. **Perfil 94d en u≈5/3:** V⁺ tiene 1–4 elementos (en e = 40, V⁺ = {1}). No conviene usar el extremo superior del rango como indicio de la constante asintótica. [M]
8. **Publicación:**
   - OneDrive no tiene los docs 95, 97 ni 100–119.
   - GitHub no tiene los docs 94, 96, 98 y 99 ni sus scripts.
   - El clon `~/dev/Lumina-La-Conjetura-De-Goldbach` tiene 7 commits sin publicar.
   - No he tocado nada de esto.

**Qué sigue abierto, sin cambios [P]:** la cota del segundo factor 𝓑_{>R0} restringido a R_A, y la comparación uniforme Q_s ≪ P_{A,N}(z_A) bajo el recorte V⁺. Ninguna de las verificaciones de este informe las acerca: solo confirman que las identidades sobre las que se construyen son exactas.

---

## §5. Cómo reproducir

```
cd "C:\Users\dal_n\OneDrive\conjetura de golbach Lumina\numerico"
python -B verificacion_97_100_claude_code.py
```

- Solo biblioteca estándar; 0,9 s en esta PC. Salida guardada en `resultados/20261128_verificacion_97_100_claude_code.txt`. SHA-256: script `0cdffb8921f4154e…`, salida `fe90f3c534cf2f0d…`.
- **Parte A (94d):** Q_s exacto (`Fraction`) por la diagonal y por la doble suma de Gram, en los 20 escenarios.
- **Parte B (94e):** soporte del doc 94 reconstruido, tabla κ con tres reglas de atribución.
- **Parte C (Δ_r):** 22 conductores con ψ reales a N = 999983, S_r por dos rutas, Δ_r(98.3) por dos rutas, y la comprobación del "espurio".
- **Parte D (doc 99):** segundo factor con N = 10⁶ y N = 999983.

Re-ejecución de los scripts de Lumina sin cambios: `python -B 94d_diagonal_profile.py`, ídem 94e, 97 y 99. Coinciden byte a byte con sus salidas de `resultados/20261127_*.txt`. La variante con N = 999983 es el script 99 con la línea 10 cambiada a `N = 999983`.

*Goldbach fuerte sigue abierta.*
