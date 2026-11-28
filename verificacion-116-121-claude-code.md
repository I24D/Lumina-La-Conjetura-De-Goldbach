# Cola [P-rev] de los docs 116–121: prioridades (peso de D-NP, auditoría DIST-D)

**Agente:** Claude Code (Anthropic, modelo `claude-opus-5-5`), brazo verificador de la línea Goldbach.
**Fecha:** 2026-11-28.
**Encargo:** de Muse, por orden de Dal (encargo 4 de 4). Extraer los [P-rev] de las secciones de cierre de los docs 116–121 y contrastarlos. Prioridades fijadas:
1. La corrección de peso de D-NP del doc 121 (W_d = Ω_d^C·φ(d)), revisando el álgebra.
2. Los veredictos de la auditoría DIST-D del doc 120: (A2) y (A6) refutadas, con sus formas corregidas.

**Método:** reimplementación propia desde las definiciones de los docs 101, 106, 115 §7, 120b §1 y 121b, sin importar código de nadie. Los caracteres se construyen con CRT y raíces primitivas propias.
Script: `numerico/verificacion_120_121_claude_code.py`. Salida: `numerico/resultados/20261128_verificacion_120_121_claude_code.txt`. **12/12 PASS**.
**Estado:** PROPUESTA hasta el contraste de Muse.
**Etiquetas:** [T] teorema, [M] medición, [P] pendiente.

**Goldbach sigue abierta.**

---

## §0. Alcance real de la cola

En los docs 116–121 **casi no hay [P-rev] en sentido estricto**, es decir, puntos pendientes de re-verificar contra un texto. Lo que listan sus cierres son **[P] abiertos**: los residuos de (H2) y (DIST-𝒟):
- (FLOW-COMP), (FC-2)_A;
- (BAND-δ), (BAND-LOW), (REG-Π);
- (BAND-MOD-HI), (SLAB-V), (SLAB-V′), (BAND-LOW-NONDEG);
- (BAND-δ-HI), (BAND-δ-HI′);
- ℳ_Ω^{liso} ≪ T(N).

Son **enunciados por demostrar**, no textos por contrastar: por construcción no admiten veredicto confirmado/refutado hasta que alguien los pruebe o los refute.

Lo verificable de esos documentos son sus afirmaciones [T] y [M]. En esta tanda hago las dos prioridades completas. El resto queda listado en el §4 como cola pendiente, **sin verificar todavía**.

---

## §1. Prioridad 1: el peso de D-NP (docs 121 y 121b)

### §1.1 El álgebra [T], revisada paso a paso

Definiciones (doc 120b §6, doc 115 §7):
- D_NP = Σ_{q,m} Θ̃_m(q)·φ(m)²·Σ_{χ≠χ₀ mod q} |S_{m,q}(χ)|².
- Ω_{qm}^C = Θ̃_m(q)/(q−1).

Pasos:
1. **Ortogonalidad.** Para (N, m) = 1 vale S_{m,q}(χ) = (1/φ(m))·Σ_{ψ mod m} ψ̄(N)·Ψ_N(ψχ). La restricción q ∤ n la impone el propio χ, porque χ(n) = 0 si q | n.
2. **Cauchy–Schwarz en ψ** (con |ψ̄(N)| = 1): |S_{m,q}(χ)|² ≤ (1/φ(m))·Σ_ψ |Ψ_N(ψχ)|².
3. **CRT** ((m, q) = 1): los pares (ψ, χ ≠ χ₀) son exactamente los η mod qm con η_q ≠ χ₀. En ambos lados hay φ(m)(q−2).
4. **Por tanto,** el término (q, m) de D_NP es ≤ Θ̃_m(q)·φ(m)²·(1/φ(m))·Σ_η |Ψ_N(η)|² = Θ̃_m(q)·φ(m)·Σ_η |Ψ_N(η)|².
5. **El peso.** Θ̃_m(q)·φ(m) = (q−1)·Ω_{qm}^C·φ(m) = **φ(qm)·Ω_{qm}^C**, porque φ(qm) = (q−1)φ(m). Así queda

   **D_NP ≤ Σ_{d=qm∈𝒟_S} Ω_d^C·φ(d)·Σ_{η mod d, η_q≠χ₀} |Ψ_N(η)|²,  W_d = Ω_d^C·φ(d).**

**Confirmado** [T]. La forma con Ω_d^C/φ(d) no se sigue de la derivación: es φ(d)² veces menor.

### §1.2 Numérico [M]

| Chequeo | Afirmación (doc 121b) | Mi cálculo | Veredicto |
|---|---|---|---|
| K1 | Fórmula con ψ̄(N) exacta en 5 casos; la variante sin conjugar falla en (200,21,13) | error ≤ 2,5·10⁻¹⁴ con ψ̄(N); sin conjugar, 17,39 en (200,21,13) | **Confirmado** |
| K2 | Pérdida de Cauchy–Schwarz λ = 2,40 (m=3), 4,06 (m=7), 11,64 (m=21), 14,52 (m=13) | 2,40 / 4,06 / 11,64 / 14,52. El caso m = 11 de (100,11,7), que el doc no lista, da 7,21 | **Confirmado** |
| — | Conteo CRT φ(m)(q−2) | coincide en los 5 casos | **Confirmado** [T] |
| K3 | En z = 17, N = 200, las 9 parejas (q, m): término real / forma Ω/φ = 1622–4768 = φ(d)²/λ | 1622, 2346, 2548, 2105, 2988, 3214, 3288, 4426, 4768, iguales a φ(d)²/λ en cada pareja; el término real cumple la cota con Ω·φ en las 9 | **Confirmado** |
| K4 | (200,35,11), con (m,N) = 5: S = log 5·χ(3) ≠ 0, mientras la fórmula da 0 | S = log 5·χ(3), con \|S\| = 1,61 | **Confirmado**: la hipótesis (m,N) = 1 es necesaria |
| doc 115 §7 | Ω_{qm}^C = Θ̃_m(q)/(q−1), identidad puntual | error relativo máximo 3,9·10⁻¹⁶ en todas las parejas (q, m) de S11, S13 y S17, contra la definición literal (106.2) | **Confirmado** |

---

## §2. Prioridad 2: la auditoría DIST-D (doc 120b) y la recuperación de (A6) (doc 121a)

### §2.1 (A2)

**Confirmado el veredicto (refutada tal como se escribió) y confirmada la forma corregida.**
- **La forma corregida es teorema** [T]. Por CRT, ψ(N; qm, N) = B_{m,q}(N mod q), con φ(qm) = (q−1)φ(m) y Δ_{q,m} = B_{m,q}(N) − (A_m − L_{m,q})/(q−1). Sustituyendo,

  Δ + (E_m − L)/(q−1) = B − ψ(N)/((q−1)φ(m)) = E(N; qm, N).
- **La versión escrita difiere exactamente en ((q−2)/(q−1))·E_m**, que no es 0 en general.
- **Medido** [M] en los 8 casos del doc 120b, incluidos m pares:
  - la corregida es exacta (< 10⁻⁹);
  - el error de la escrita coincide con ((q−2)/(q−1))E_m a 6 decimales: −1,313565; −36,297433; −38,923459; −90,270306; −3,035147; 1,664822; −221,376625; −0,025823.
- **(A1)** Σ_{a≠0} B(a) = A_m − L y **(A7)** Parseval también se confirman en los 8 casos.

### §2.2 (A6)

**Confirmado el veredicto del doc 120b (refutada como estaba escrita) y confirmada la recuperación del doc 121a.**

**Lo refutado:**
- El α_m auditado lleva el factor invertido: α_m = Σ_q Θ̃_m(q)(q−1) = Σ_q (q−1)²·Ω_{qm}^C.
- Medido [M]: **viola** α_m ≤ Ω_m^C en todos los m lisos de S11 (3 de 3), S13 (7 de 7) y S17 (15 de 15), con cocientes de hasta 3·10⁶.

**Lo recuperado** (doc 121a) [T]:
- α̃_m = Σ_q Ω_{qm}^C = Σ_{r∈ℋ_C, m|r} (τ(r)−1)/φ(r)·#{q ∈ 𝒮 : q | r}.
- Por la **unicidad de rama**, todo r ∈ ℋ tiene a lo sumo un primo ≥ z (Lema 2 del doc 101). Luego #{q : q | r} ≤ 1 y α̃_m ≤ Ω_m^C.
- Medido [M]: la unicidad se cumple en todo r de S11–S17, y **0 violaciones** de α̃_m ≤ Ω_m^C.

**El contraejemplo del doc 120b** (ℛ = {3·5·7}, m = 3) se reproduce: 0,2917 > 0,1458. Necesita dos primos de conmutación en un mismo r, lo que es imposible en ℋ, así que queda fuera del soporte real, como dice el doc 121.

**El peso del intercambio** (doc 121a §(d)), en S11 con m = 5:
- Σ_q Ω_{5q}^C = **1,091744**.
- Con el peso τ(r)⁻¹/φ(r) se obtiene **0,072937**, un factor **15,0**.
- La identidad escrita con τ⁻¹/φ es falsa. La desigualdad con el peso correcto (τ−1)/φ sobrevive. **Confirmado.**

---

## §3. Resumen de veredictos

| Punto | Veredicto |
|---|---|
| Doc 121/121b: W_d = Ω_d^C·φ(d) (la "cuarta inversión de φ") | **Confirmado** [T + M] |
| Doc 121b K1–K4 (ortogonalidad con ψ̄(N), Cauchy–Schwarz con λ, conteo CRT, necesidad de (m,N) = 1) | **Confirmados** [T + M] |
| Doc 115 §7: Ω_{qm}^C = Θ̃_m(q)/(q−1) puntual | **Confirmado** [M, error 4·10⁻¹⁶] |
| Doc 120b (A2): refutada tal como se escribió; forma corregida E = Δ + (E_m − L)/(q−1) | **Confirmados ambos** [T + M] |
| Doc 120b (A6): refutada tal como se escribió (factor (q−1) invertido) | **Confirmado** [T + M] |
| Doc 121a: (A6) recuperada como α̃_m ≤ Ω_m^C por unicidad de rama | **Confirmado** [T + M] |
| Doc 121a §(d): la identidad con τ⁻¹/φ es falsa (factor ~15) | **Confirmado** [M] |
| Doc 120b (A1), (A7) | **Confirmados** [M] |

---

## §4. Lo que queda en la cola (NO verificado en esta tanda) [P]

Afirmaciones [T]/[M] de los docs 116–120 que admiten contraste y que **no he verificado todavía**:
- **Doc 116:** (REM-PERS-small) cerrado [T]; (R-FLOW-71+) falso puntualmente, con contraejemplos en el marco [T]; (BAND) cruda falsa asintóticamente por conductores primos [T].
- **Doc 117:** (FLOW-COMP) reducido a (FC-2) y probado en el marco de su contraejemplo, con margen ×28; (BAND-δ) reducido a (DESP-SIGN) + (BAND-MOD) [T]; (REG-Π) medido en 0,17–0,46·Q_c.
- **Doc 118:** (REG-Π) y (REG-Π-LOW) falsas, porque el mayorante diverge como log log z [T]; estratos q₀ ≥ 13 cerrados [T]; Teorema A (L ≥ q₀⁹ ⟹ ningún degenerado visible con ω ≤ 2) [T].
- **Doc 119:** (BAND-MOD-HI) falsa [T módulo FL]; (SLAB-V) reducido a (SLAB-V′) por reindexación [T].
- **Doc 120a:** (BAND-δ-HI): la cancelación puntual no ocurre. Estadísticas a contrastar: E(r) = 0 en el 58,2 % de los casos con exceso, mediana de afeitado 62,7 %.

Son el siguiente lote natural de la cola. Necesitan leer íntegros los docs 116–120 (unos 150 KB) y reimplementar sus predicados (flujos Add/Rem, bandas, Π(r), Ŝ). Propongo atacarlos en la siguiente tanda, por orden de dependencia: 116 → 117 → 118 → 119 → 120a.

*Goldbach fuerte sigue abierta.*
