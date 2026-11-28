# DOC 95 — Revisión adversarial independiente de la derivación del doc 94
## Norma de Gram y reducción UNIF-U-COND

**Fecha:** 2026-11-27. **Línea:** dedicada Goldbach (Muse). **Estado:** revisión de la
derivación transcrita; el doc 94 completo NO fue revisado (ver §0). **Goldbach sigue abierta.**

---

## §0. Limitación de fuentes (léase primero)

El doc 94 completo vive en OneDrive en el PC de Dal y **no está disponible en esta revisión**:
el túnel al PC está caído (segunda caída de egress del 2026-11-27, desde ~12:16 EST).
Lo que se revisa aquí es **únicamente** la derivación que Lumina Chat compartió en el chat
«Revisa UNIF U 2M X» el 2026-11-27, transcrita íntegra en la tarea, más los docs 85–93
en `~/workspace/goldbach/` como contexto.

Consecuencia: los ítems marcados **[P-rev]** pueden estar ya resueltos dentro del doc 94
completo (por ejemplo, el puente G_A → R_A, o la discusión del régimen (N,Am)=1).
**Esta revisión debe re-verificarse contra el doc 94 íntegro cuando el túnel se restablezca.**
Nada de lo aquí afirmado depende de haber leído el doc 94: todo lo verificado lo fue sobre
la derivación transcrita y sobre los cálculos independientes del script acompañante.

**Script:** `95_revision_adversarial_doc94_gram.py` (solo stdlib, ~7 s).
**Salida:** `salida_doc95_local.txt`. **Resultado: PASS (exit 0).**
Nota de proceso: la primera corrida válida (run 2) contenía un bug del propio script de
revisión — `pf_desc` devolvía los primos en orden creciente mientras la condición de cadena
de V⁺ (Nathanson, Teorema 9.3) y el chequeo de límite primo suponen orden decreciente —, lo
que admitía módulos no squarefree y rompía (s,d)=1. El bug se encontró, se corrigió
(`pf_desc` ahora devuelve `sorted(..., reverse=True)`) y **todos los números de este doc
son de la corrida corregida (run 3)**. La corrida 2 queda invalidada y no se cita salvo
para documentar el error.

**Notación de veredicto por paso:** [T] teorema / verificado · [N] no tenido (falso o no
probado como está enunciado) · [P-rev] pendiente de revisión (no aparece en la transcripción
y hace falta).

**Pesos usados en la reproducción:** λ⁺ de Nathanson (doc 86, Teorema 9.3, β=2), igual que
en los docs 86–88. Escenarios: S1 (A=3, z=11, y=60, D=20000), S2 (A=3, z=17, y=120,
D=200000), S3 (A=15, z=19, y=150, D=500000), S4 (A=3, z=11, y=3000, D=200000,
parámetros aproximados a los reportados por Lumina). Experimento de tendencia:
D=z⁵, y=z^(10/3), A=3, z=11…29. Todo lo numérico es **medición** en estos rangos, no prueba.

---

## §1. La derivación a revisar (transcripción)

De R_A = Σ_{m∈M_A} b_A(m)·E(N;Am,N), la ortogonalidad de caracteres agrupa por conductor
primitivo. Para un carácter primitivo de conductor r|Am, a=(r,A), h=r/a; los conductores
admisibles cumplen a|A, h|m, μ²(h)=1, (h,AN)=1. Como (A,m)=1, φ(Am)=φ(A)φ(m), el
coeficiente común de los caracteres con la misma parte nueva h es

    K_A(h) = (1/φ(A)) · Σ_{m∈M_A, h|m} b_A(m)/φ(m).

Definiendo G_A = Σ_{h∈H_A} φ(h)|K_A(h)|² y expandiendo el cuadrado:

    G_A = (1/φ(A)²) · Σ_{m,n} b_A(m)b_A(n)·(m,n)/(φ(m)φ(n)),

usando Σ_{d|n} φ(d)=n. Segunda reducción: cada módulo es m=sd con un único primo de
switching s:

    T_s(h0) = Σ_{d∈V⁺(D_A/s), h0|d} μ(d)/φ(d);
    K_A(h0) = (1/φ(A)) · Σ_s T_s(h0)/φ(s),
    K_A(s·h0) = T_s(h0)/(φ(A)φ(s)).

La energía se parte en dos: G_A = G_smooth + G_switch. Objetivo analítico:

    sup_{z_A ≤ s < y_A, s primo, (s,AN)=1} Q_s ≪ 1/log z_A,
    con Q_s = Σ_{h0} φ(h0)|T_s(h0)|².

Números reportados por Lumina (no tomados como dados): 1.444 pesos, 1.451 partes nuevas
distintas, norma signada 0,4505669522, mayorante absoluto 1,1723577317, cociente 0,384325;
20 experimentos N=2⁴⁰…2⁶⁰ con norma de cada bloque ~ cte/log z.

---

## §2. Verificación línea por línea

### 2.1. De R_A a K_A(h): agrupación por conductor primitivo — [P-rev]

La agrupación en sí es correcta bajo las hipótesis que la transcripción enuncia. Verificación
de cada ingrediente:

- **(A,m)=1 ⇒ φ(Am)=φ(A)φ(m).** [T] Trivial. Es hipótesis del marco (los pesos b_A viven
  en módulos coprimos con A); si algún m con (A,m)>1 entrara en el soporte, la fórmula de
  K_A(h) cambia. El soporte debe declararse explícitamente.
- **r|Am, a=(r,A), h=r/a ⇒ h|m y (h,A)=1.** [T] En efecto: a contiene todos los factores
  comunes de r y A, luego (h,A)=1; de ah|Am y (h,A)=1 sale h|m.
- **μ²(h)=1.** [T] *dado* que supp b_A ⊆ {m squarefree} (los pesos λ⁺ de Nathanson lo son):
  h|m implica h squarefree. La condición no viene de la teoría de conductores — un conductor
  primitivo puede ser no squarefree (p. ej. 4 u 8) — sino del soporte de los pesos. Debe
  decirse así; si el soporte cambiara, la condición cae.
- **(h,AN)=1 como admisibilidad.** [P-rev] Aquí hay un punto que la transcripción no
  desarrolla y que importa. La expansión por caracteres de E(N;Am,N) — es decir,
  ψ(N;Am,N) = (1/φ(Am))·Σ_χ χ̄(N)·ψ(N,χ) — **exige (N,Am)=1**: si (N,Am)>1, todos los
  χ̄(N) se anulan y la fórmula de ortogonalidad no se aplica al miembro izquierdo. La
  condición (h,N)=1 no basta: hace falta (m,N)=1 (y (A,N)=1) para cada m del soporte, porque
  la expansión se hace al módulo Am, no al conductor h. Lo correcto es declarar
  supp b_A ⊆ {m : (m,N)=1} — que es exactamente el régimen de la corrección E93.1 del
  doc 93 ((q,N)=1) — y tratar aparte los m con (m,N)>1 si los hubiera. En la práctica del
  marco (A impar, m impar, N par con sus primos impares excluidos del soporte) esto se
  cumple, pero debe escribirse.
- **El carácter inducido vs. el primitivo.** [P-rev] La agrupación recoge los términos
  ψ(N,χ), con χ inducido módulo Am, bajo su primitivo χ*. Pero ψ(N,χ) ≠ ψ(N,χ*): difieren
  en los n ≤ N con (n,Am)>1 y (n,r)=1, es decir, en potencias de primos p|Am con p∤r.
  El error por carácter es ≤ Σ_{p|Am, p∤r} log N ≪ ω(Am)·log N, y sumado sobre los
  caracteres da un error total ≪ (log N)(log log N)·Σ_m |b_A(m)|. Es pequeño en el régimen
  de aplicación (D_A ≪ N), pero **no aparece en la transcripción**: el paso de R_A a la
  forma de Gram necesita este término de error escrito y acotado. Del mismo modo, la pieza
  del carácter principal, (ψ(N,χ₀)−N)/φ(Am), no se menciona (es pequeña por el TNP, pero
  debe constar).
- **Comprobación numérica del mecanismo de ortogonalidad (H8d3).** Verificado
  independientemente en el script para r=3, N=128: ψ(N,χ) directo = 128,5754;
  vía E(N;3,N) = 128,5754; diferencia 1,36·10⁻¹⁰. [T] El mecanismo funciona donde aplica.

**Veredicto §2.1:** la fórmula de K_A(h) es correcta en su régimen; el paso R_A → forma de
Gram queda [P-rev] hasta que se escriban el error imprimitivo, la pieza principal y la
hipótesis supp b_A ⊆ {(m,N)=1}.

### 2.2. La identidad de Gram — [T]

    G_A = (1/φ(A)²) · Σ_{m,n} b_A(m)b_A(n)·(m,n)/(φ(m)φ(n)).

**Verificación (H1):** dos cálculos independientes — (i) directo desde K_A(h) con
G_A = Σ_h φ(h)|K_A(h)|², y (ii) forma bruta por la doble suma — coinciden en los seis
escenarios con error relativo ≤ 1,8·10⁻¹⁴:

| escenario | G_A(K) | G_A(bruto) | error rel. |
|---|---|---|---|
| S1 fijo | 0,149757903123 | 0,149757903123 | 4,2·10⁻¹⁶ |
| S1 variable | 0,144220501132 | 0,144220501132 | 6,1·10⁻¹⁶ |
| S2 fijo | 0,120153927360 | 0,120153927360 | 1,2·10⁻¹⁵ |
| S2 variable | 0,112449015986 | 0,112449015986 | 4,1·10⁻¹⁵ |
| S3 fijo | 0,008008179537 | 0,008008179537 | 4,9·10⁻¹⁶ |
| S4 fijo | 0,469133545584 | 0,469133545584 | 1,8·10⁻¹⁴ |

El lema divisor Σ_{d|n} φ(d)=n se usa correctamente: al expandir
Σ_h φ(h)·(Σ_{m:h|m} b_A(m)/φ(m))·(Σ_{n:h|n} b_A(n)/φ(n)), el coeficiente de
b_A(m)b_A(n)/(φ(m)φ(n)) es Σ_{h|(m,n)} φ(h) = (m,n). **El factor (m,n)/(φ(m)φ(n)) es
exacto.** Además la identidad se mantiene **con o sin** la restricción (m,N)=1 del §2.1:
es una identidad combinatoria pura sobre K_A, independiente del lado ψ.

**Un detalle que no debe perderse (H6):** |K_A(1)|² es el 16–34 % de G_A (S1: 20,0 %;
S2: 15,9 %; S3: 20,6 %; S4: 34,0 %). La fuerza bruta solo reproduce G_A si 1∈H_A está
incluido: sin h=1 el error relativo es 2–16 %. El término h=1 (conductores r=a, la parte
«vieja») es una fracción sustancial de la energía, no un residuo.

**Veredicto §2.2: [T]. La identidad de Gram sobrevive íntegra.**

### 2.3. Segunda reducción y partición G_smooth + G_switch — [T] en límite fijo, [N] en límite variable

**Versión de límite fijo (doc 87: el primo de switching s es el mayor factor primo, o el
límite superior es fijo).** [T] Verificado exactamente (H2):

| escenario | K(h0) maxerr | K(s·h0) maxerr | partición |
|---|---|---|---|
| S1 fijo | 2,8·10⁻¹⁷ | 6,9·10⁻¹⁸ | exacta (0,0) |
| S2 fijo | 8,3·10⁻¹⁷ | 3,5·10⁻¹⁸ | exacta (0,0) |
| S3 fijo | 5,6·10⁻¹⁷ | 8,7·10⁻¹⁹ | exacta (1,7·10⁻¹⁸) |
| S4 fijo | 2,8·10⁻¹⁶ | 6,9·10⁻¹⁸ | exacta (0,0) |

G_smooth + G_switch = G_A exactamente en los cuatro casos (S1: 0,0529 + 0,0968 = 0,1498).

**Versión de límite variable (doc 85: el límite superior de la cadena depende de s).**
[N] **Contraejemplo confirmado (H3).** La fórmula K_A(s·h0) = T_s(h0)/(φ(A)φ(s)) es
**falsa** como está enunciada: un módulo m divisible por s·h0 puede tener primo de switching
s' > s (con s·h0 | d, d ∈ V⁺(D/s')), y ese término falta.

- S1 variable: máx |K(s·h0) − T_s/(φ(A)φ(s))| = 0,002431, en (s=11, h0=1):
  K_A(11) = 0,0288194 frente a lo afirmado 0,03125 (error relativo ~8 %).
- S2 variable: máx discrepancia 0,002332, en (s=17, h0=1): K_A(17) = 0,0138789 frente a
  0,0162109 (~14 %).
- Los h con ≥2 primos de switching («multi») son 3 en S1 y 37 en S2; su energía es
  G_multi/G_A = 1,09 % (S1) y 2,78 % (S2): el error es pequeño en estos rangos, pero es
  un error, no un redondeo.

**Descomposición corregida (verificada, H3).** La fórmula exacta es

    K_A(s·h0)·φ(A) = T_s(h0)/φ(s) + Σ_{s'>s} U_{s'}(s·h0)/φ(s'),
    U_{s'}(t) = Σ_{d∈V_{s'}, t|d} μ(d)/φ(d),

con error máximo 2,1·10⁻¹⁷ (S1) y 1,0·10⁻¹⁷ (S2). Es decir: la segunda reducción se
repara añadiendo la cola s'>s; el coste es que K_A(s·h0) ya no depende solo del bloque s,
y la partición G_smooth + G_switch necesita un tercer sumando G_multi (1–3 % de G_A en
los modelos; su tamaño asintótico no está medido).

**Veredicto §2.3:** [T] en límite fijo; **[N] en límite variable tal como está enunciada**,
con reparación exacta conocida. El doc 94 debe declarar en qué versión trabaja; si usa la
variable, debe usar la descomposición corregida.

### 2.4. La forma de Gram de Q_s — [T] en fijo, [N] en variable (versión smooth)

**Versión fija (H4):** Q_s directo vs. forma de Gram
Σ_{d,e} μ(d)μ(e)(d,e)/(φ(d)φ(e)) coinciden con error máximo ~10⁻¹⁵ en todos los escenarios.
[T]

**Versión variable (H4):** la Q_s «smooth» de la transcripción (h0 con primos < z) **difiere**
de la forma de Gram: error máximo 0,0757 (S1, s=13) y 0,0880 (S2, s=23). La expansión de
Gram necesita **todos** los divisores de (d,e), incluidos los h0 con primos en [z,s). La
versión «full» (h0 sobre todos los divisores de algún d ∈ V_s) sí coincide exactamente
(error ~10⁻¹⁵). [N] para Q_s^smooth en modo variable; [T] para Q_s^full y para el modo fijo
(donde ambas coinciden).

Consecuencia práctica: en la versión variable, o bien se define Q_s con el conjunto completo
de h0, o bien la reducción al sup Q_s debe reescribirse. No es un detalle de notación: la
diferencia es del orden del 10 % de Q_s.

### 2.5. El objetivo sup Q_s ≪ 1/log z_A — [N] como teorema, consistente como medición

**Medición (experimento de tendencia, modo fijo, D=z⁵, y=z^(10/3), A=3):**

| z | #s | Q_s·ln z (rango) | máx Q_s·ln z | diag·ln z (máx) |
|---|---|---|---|---|
| 11 | 422 | [1,68, 2,40] | 2,40 | 3,76 |
| 13 | 682 | [1,63, 2,56] | 2,56 | ~4,2 |
| 17 | 1502 | [1,67, 2,83] | 2,83 | ~4,6 |
| 19 | 2090 | [1,63, 2,39] | 2,39 | ~4,8 |
| 23 | 3687 | [1,65, 2,55] | 2,55 | ~5,2 |
| 29 | 7381 | [1,69, 2,74] | 2,74 | 7,54 |

Lectura:

- **Q_s ≈ C/ln z con C ≈ 2,4–2,8 (máx sobre s).** El objetivo ≪ 1/log z_A, entendido en
  notación de Vinogradov (≤ C·(1/log z) con C absoluta), es **consistente con la medición**:
  Q_s·ln z se mantiene acotado y no crece más allá de ~2,8 en z=11…29.
- **Pero no hay decaimiento más allá de 1/ln z:** el máx de Q_s·ln z no decrece (2,40 →
  2,74, no monótono) y el mín se mantiene ~1,63–1,69. Si el argumento aguas abajo necesitara
  sup Q_s = o(1/ln z), o ≪ 1/(ln z)², o una constante C < 1, **fracasaría**. El objetivo
  sobrevive solo en la forma ≪ con constante C ≈ 2,5–3.
- **La diagonal sola ya es O(1/ln z):** diag = Σ_{d∈V_s} d/φ(d)² cumple diag·ln z ≈
  3,8–7,5. Es decir, el objetivo **no depende** de la cancelación fuera de la diagonal en
  estos rangos: una cota puramente diagonal daría ya el orden correcto, con constante peor
  (~4–7,5 en vez de ~2,5).
- **La cancelación fuera de la diagonal existe y es sustancial:** diag/Q_s ≈ 2,2–4,3
  (S1: 2,2; S2: 3,2–4,3): los términos d≠e cancelan el 55–77 % de la diagonal. Mejora la
  constante en un factor ~2–4, pero no cambia el orden.

**Pérdida de Cauchy (H5):** el paso G_smooth ≤ (Σ_s 1/φ(s))·G_switch pierde solo un factor
1,007–1,062 en los modelos. La cota implicada para G_A a partir de sup Q_s ≪ C/ln z queda
a un factor 1,05–1,41 del valor real de G_A: **la cadena es ajustada, sin holgura**. Cualquier
paso posterior que necesite G_A más pequeño por un factor log perderá.

**Veredicto §2.5:** el objetivo sobrevive como enunciado de orden de magnitud y es
consistente con la medición hasta z=29; **no está probado** [N], la constante es C ≈ 2,5
(no ≪ 1) y no hay holgura aguas abajo.

---

## §3. Caza de contraejemplo y pérdidas ocultas

Preguntas adversarialas y respuestas:

1. **¿Se pierde cancelación por tomar valor absoluto en algún paso?** No en la derivación
   transcrita: la identidad de Gram es exacta y conserva los signos (H1: error ~10⁻¹⁶).
   La única toma de valor absoluto es la de Cauchy en §2.5, que pierde 1,007–1,062. **Pero**
   hay una pérdida latente: el cociente norma signada / mayorante absoluto es
   G_A/G_A^abs ≈ 0,12–0,29 (√ ≈ 0,35–0,53) en los seis escenarios — es decir, K_A ya
   cancela un factor ~3–8 respecto a su mayorante absoluto. **Cualquier paso posterior que
   acote G_A por G_A^abs perdería ese factor.** El doc 94 debe no hacerlo; si algún lema
   posterior lo hace, ahí está la pérdida.
2. **¿La identidad conserva realmente los signos?** Sí (H1 exacto). La cancelación vive
   dentro de K_A(h) — los μ(d) y los b_A(m) con signos — y la identidad la transporta
   intacta a la forma (m,n)/(φ(m)φ(n)).
3. **¿Hay un factor perdido (φ(A), log, D_A)?** Revisado:
   - φ(A): la identidad lleva 1/φ(A)² y es exacta; no falta ningún φ(A). Pero ojo: en la
     cadena hacia R_A, el factor φ(A)² multiplica a G_A, y si A crece con N el objetivo
     para G_A debe leerse con ese factor (ver §4).
   - log: el objetivo 1/log z es el orden medido; no hay un log perdido en la identidad
     (es exacta), pero la constante C ≈ 2,5 no es despreciable.
   - D_A: el número de módulos no entra en la identidad; entra en el conteo del soporte
     (ver §5, corrección d1) y en el tamaño del error imprimitivo del §2.1.
4. **¿El sup sobre s es el objeto correcto?** Sí en el sentido de que el máx de Q_s·ln z
   es estable (~2,4–2,8) y el mín ~1,6: el sup no está dominado por un s patológico en
   estos rangos. Pero el sup se alcanza en s distintos según z (s=59, 97, 103, 1597…): no
   hay un s extremal fijo; el argumento analítico debe ser uniforme en s.
5. **¿La partición smooth/switch pierde algo?** En modo fijo, nada (exacta). En modo
   variable, pierde G_multi (1–3 %) y la fórmula de K(s·h0) es falsa (§2.3).

**Ningún contraejemplo a la identidad de Gram.** El único contraejemplo encontrado es el de
§2.3 (versión variable de la segunda reducción), con reparación exacta.

---

## §4. El puente que falta: de G_A a R_A — [N]

**Este paso no está en la derivación transcrita.** Y es el paso que decide si todo lo
anterior sirve para UNIF-U-COND.

Estimación ingenua (H7): por Cauchy y la gran criba, la cota ingenua de |R_A| a partir de
G_A pierde un factor

    (cota|R_A|)²/(objetivo)² ≈ d(A)·A²·(ln N)⁵/ln z,

es decir, ~7,4·10⁷ (S1), ~1,4·10⁸ (S2) y ~1,1·10¹⁰ (S3) en los modelos: la cota ingenua
es ~10⁴–10⁵ veces demasiado grande. El puente real debe explotar que solo aparecen los
h ∈ H_A (pocos, ~10³ en los modelos, frente a todos los conductores ≤ D_A), que los
caracteres son primitivos y que las fases χ̄(N) no son arbitrarias. Nada de eso está en la
transcripción.

**Veredicto:** aunque se probara sup Q_s ≪ 1/log z, eso por sí solo **no** da
R_A ≪ (N/A)/log²N. La reducción a sup Q_s es una reorganización correcta, no una prueba
de UNIF-U-COND. El trabajo analítico restante está en el puente, no en la norma.

---

## §5. Las tres correcciones de Lumina (incorporadas)

- **(d1) El número de módulos del soporte es O(D_A), no necesariamente ≤ D_A.** [T en los
  modelos] Medido (H8d1): S1: 50 pares (s,d) / 50 distintos con D_A=20000; S4: 876 pares /
  876 distintos con D_A=200000. En los modelos, distintos ≤ D_A se cumple con holgura
  (los pares son ≪ D_A porque el soporte es ralo: la cadena de V⁺ restringe mucho). Pero
  el enunciado general correcto es O(D_A): el conteo de pares (s,d) con d ∈ V⁺(D_A/s) no
  tiene por qué ser ≤ D_A en general, y el doc 94 debe enunciarlo como O(D_A) con la
  constante explícita donde se use.
- **(d2) La estimación de bloque debe restringirse a ah>R0.** [P-rev, cuantificado] Medido
  (H8d2): la fracción de Q_s que proviene de h0 ≤ 10 es 0,94 (S1/S2) y 0,75 (S3, A=15);
  de h0 ≤ 30 es 0,94/0,91; de h0 ≤ 100 es 1,00/0,96. **El rango ah ≤ R0 lleva la mayor
  parte de la masa (~94 % en los modelos).** La restricción a ah>R0 no es un tecnicismo:
  deja fuera el grueso de Q_s, que necesita tratamiento aparte (cota directa del sumatorio
  finito de h0 pequeños, o input tipo Siegel–Walfisz en el lado ψ). Sin ese tratamiento,
  la estimación de bloque no cubre Q_s.
- **(d3) «Estrictamente más fuerte que GRH» degradado a implicación en un solo sentido.**
  [T la implicación; N la recíproca] UNIF-U-2M-R ⇒ GRH está probada (doc 93, Prop. 93.6;
  verificada la maquinaria de ortogonalidad en H8d3 con diferencia 1,36·10⁻¹⁰). La
  no-implicación inversa (GRH ⇏ UNIF-U-2M-R) **no** está probada. El doc 94 debe enunciar
  solo la implicación UNIF-U-2M-R ⇒ GRH y no afirmar fuerza estrictamente mayor.

---

## §6. Reproducción de los números reportados por Lumina

Comparación en el escenario más cercano (S4: A=3, z=11, y=3000, D=200000):

| cantidad | Lumina (doc 94) | esta revisión (S4) |
|---|---|---|
| pesos \|M\| | 1.444 | 876 |
| partes nuevas \|H\| | 1.451 | 880 |
| norma signada G_A | 0,4505669522 | 0,4691335456 |
| mayorante absoluto G_A^abs | 1,1723577317 | 1,6425047363 |
| cociente G/G^abs | 0,384325 | 0,285621 |

Los parámetros no coinciden exactamente (Lumina usa más pesos: 1.444 frente a 876; su y o
su D son mayores), así que no es una reproducción exacta. Pero **G_A coincide al 4 %**
(0,4691 frente a 0,4506) y el cociente está en el mismo orden (0,29 frente a 0,38): los
números reportados son plausibles como mediciones del mismo objeto y en la misma escala.
El mayorante absoluto difiere más (40 %), en la dirección esperada al aumentar el soporte.
Los 20 experimentos N=2⁴⁰…2⁶⁰ con norma ~ cte/log z no son reproducibles en esta VM en
tiempo razonable (N hasta 2⁶⁰ exige criba de D_A enormes); la tendencia medida hasta
z=29 (§2.5) es consistente con su afirmación como orden de magnitud, no como verificación.

---

## §7. Veredicto final

1. **La identidad de Gram sobrevive.** Es exacta: verificada por dos cálculos
   independientes con error relativo ≤ 1,8·10⁻¹⁴ en seis escenarios, en modo fijo y en modo
   variable. El factor (m,n)/(φ(m)φ(n)) y el uso de Σ_{d|n} φ(d)=n son correctos. No hay
   pérdida oculta de cancelación en ella; conserva los signos.
2. **La segunda reducción sobrevive en la versión de límite fijo** (exacta a ~10⁻¹⁷) y
   **falla como está enunciada en la versión de límite variable** (contraejemplo explícito,
   §2.3), pero se repara con la cola s'>s, verificada exacta a ~10⁻¹⁷. El doc 94 debe
   declarar su versión y, si es la variable, usar la descomposición corregida y añadir
   G_multi.
3. **El objetivo sup Q_s ≪ 1/log z sobrevive como enunciado de orden de magnitud**
   (medido: Q_s·ln z ∈ [1,63, 2,83], máx ≈ 2,4–2,8), pero no está probado, la constante es
   ~2,5 (no pequeña), no hay decaimiento más allá de 1/ln z, y la cadena hacia G_A es
   ajustada (factor 1,05–1,41, sin holgura).
4. **La reducción a sup Q_s no basta.** El puente G_A → R_A no está en la derivación y la
   vía ingenua pierde ~10⁷–10¹⁰. Sin puente, UNIF-U-COND no se sigue de sup Q_s.
5. **El rango ah ≤ R0 es el grueso del problema** (~94 % de Q_s en los modelos), no un
   resto: la estimación de bloque restringida a ah>R0 necesita su complemento.

**En una frase:** la maquinaria de Gram es correcta y la reducción es una reorganización
válida (con la corrección de la versión variable), pero el trabajo analítico real — probar
sup Q_s ≪ 1/log z con constante explícita, tratar el rango ah ≤ R0 y construir el puente
hasta R_A — sigue entero por hacer.

### Siguiente paso analítico concreto (propuesto)

1. **Probar Q_s ≪ 1/ln z para ah>R0 por la vía diagonal:** en los modelos la diagonal
   sola ya da el orden (diag·ln z ≈ 3,8–7,5), así que una cota diagonal
   Σ_{d∈V_s} d/φ(d)² ≪ 1/ln z — sin necesitar la cancelación fuera de la diagonal — puede
   bastar para el orden de magnitud. La cancelación (factor ~2–4) queda como mejora de la
   constante, no como necesidad. Esto convierte el objetivo en un problema de criba
   elemental sobre V_s.
2. **Tratar ah ≤ R0 por separado y primero:** es donde vive ~94 % de Q_s. Para h0
   pequeños y fijos, T_s(h0) es una suma de criba con condición de divisibilidad fija; una
   cota directa (trivial + conteo de V_s) o un input de tipo Siegel–Walfisz en el lado ψ
   debe cuantificarse antes que el régimen ah>R0.
3. **Construir el puente G_A → R_A** con la gran criba restringida a h ∈ H_A y a
   caracteres primitivos, explotando las fases χ̄(N): el objetivo es recuperar el factor
   ~10⁷–10¹⁰ que pierde la vía ingenua. Sin este paso, el resto no implica UNIF-U-COND.
4. **Re-verificar esta revisión contra el doc 94 íntegro** cuando el túnel se restablezca:
   en particular si el doc 94 ya contiene el puente (§4), el error imprimitivo (§2.1) y la
   versión (fija/variable) de la segunda reducción.

---

## §8. Registro de verificación

- Script: `95_revision_adversarial_doc94_gram.py` — solo stdlib (math, sys, time).
  Checks H0–H9 en 6 escenarios (fijo/variable × S1/S2/S3/S4) + experimento de tendencia
  z=11…29. **Resultado: PASS, exit 0, 7,4 s.**
- Salida completa: `salida_doc95_local.txt` (120 líneas en la corrida válida).
- Corrección del script durante la revisión: orden de `pf_desc` (decreciente) — la corrida
  2, con el bug, queda invalidada; ningún número de este doc procede de ella.
- No se tocó GitHub, no se escribió a Lumina Chat, no se modificaron los docs 85–94.
- Goldbach sigue abierta. Ninguna afirmación de este doc demuestra Goldbach ni ninguna de
  sus implicaciones.
