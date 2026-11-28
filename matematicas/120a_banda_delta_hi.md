# Doc 120a — (BAND-δ-HI): ¿cancela E al exceso del modelo en la capa alta?

**Fecha:** 2026-11-28. **Línea:** UNIF-U / (H2), continuación del doc 119 (§§1–2).
**Documento hermano:** ninguno; esta sección es autocontenida.
**Etiquetas:** [T] = teorema, [M] = medición en ventana, [P] = pendiente.

---

## 0. El enunciado y la pregunta

Con δ(r) = S(r) − Π(r) = (Ŝ(r) − Π(r)) + E(r) [T, doc 113], el residuo
que el doc 119 dejó en la capa alta de la banda es:

**(BAND-δ-HI)** [P]: Σ_{r ∈ banda, u_r ≥ u₀*} δ(r)⁺/φ(r) ≤ C·Q_c,

con banda = {r visible: u_r < u₀ = 3/2, Π(r) < 1/2} y u₀* el umbral de la
FL del modelo (u₀* < 3/2; en ventana se mide la capa u_r ∈ [1, 3/2) y la
sub-rebanada [1,4, 1,5) donde se concentra la masa de ℱ).

El doc 119a probó que (Ŝ−Π)⁺ ≍ Π sobre las familias obstructoras:
(BAND-MOD-HI) es falsa y ninguna ruta por el modelo solo cierra la capa
alta. Por tanto (BAND-δ-HI) vive o muere con la **cancelación**: como
E(r) ≤ 0 (si (DESP-SIGN)), vale puntualmente

δ(r)⁺ = ((Ŝ−Π)(r) − |E(r)|)⁺  [T, dado E(r) ≤ 0],

y (BAND-δ-HI) equivale a que el **exceso no cubierto**
Σ ((Ŝ−Π)⁺ − |E|)⁺/φ sea ≤ C·Q_c sobre la capa. La pregunta de esta
sección es empírica antes que teórica: **¿cancela E de verdad, y dónde?**

Tres sub-preguntas, en orden:
(1) [M] En ventana (marcos A/B y estrés), ¿cuánto vale Σ δ⁺/φ sobre la
capa alta, contra Q_c? ¿La cancelación es r a r o solo en suma? ¿Cuál es
el cociente |E|/(Ŝ−Π)⁺ típico y mínimo sobre la capa?
(2) [T/M] ¿Qué parte de E hace el trabajo — E_Add o E_Rem (doc 114)? ¿Qué
enunciado mínimo exacto apoya la medición?
(3) ¿Dice algo el margen ×28 de (FLOW-COMP) (doc 117a) sobre esta
cancelación, o son demandas independientes sobre E?

## 1. Medición decisiva de la capa alta [M]

Script `120a_banda_delta_hi.py` (PASS, K0–K3, 2,2 s). Marcos A
(z = 17, 29, 37) y B (z = 29), 60 niveles por marco, y estrés z = 53
(4 niveles con u_L ≈ 2). Capa medida: u_r ∈ [1, 3/2); sub-rebanada
[1,4, 1,5). K0 verifica las identidades δ = (Ŝ−Π) + E y
E = Σ_Add a − Σ_Rem a a precisión de máquina, y E ≤ 0 en la capa.

**La suma de capa es pequeña en ventana.** Máximos sobre niveles, en
unidades de Q_c:

| marco | Σ δ⁺/φ | Σ (Ŝ−Π)⁺/φ | Σ \|E\|/φ | δ⁺ en [1,4, 1,5) |
|---|---|---|---|---|
| A, z=29 | 0,0198 | 0,1028 | 0,0830 | 0,0198 |
| A, z=37 | 0,0933 | 0,2103 | 0,1244 | 0,0933 |
| B, z=29 | 0,0128 | 0,0128 | 0,0000 | 0,0073 |
| estrés, z=53 | 0,0867 | 0,2695 | 0,1828 | 0,0491 |

(z = 17 no tiene instancias de banda en la capa.) La sub-rebanada
[1,4, 1,5) por sí sola alcanza el mismo orden máximo que la capa entera
en los marcos A: la parte alta de la capa no es despreciable, pero
tampoco crece en ventana. Como en los docs 118–119, esto **no decide**
el enunciado asintótico: la familia ℱ sigue ausente de la banda en
ventana (0 miembros; Π(pq) ≥ 1/2 los excluye), y son exactamente sus
conductores los que decidirían (BAND-δ-HI) cuando z → ∞.

**La cancelación puntual no ocurre.** Sobre las 67 instancias no
degeneradas de la capa con (Ŝ−Π) > 0:

- cobertura completa |E| ≥ (Ŝ−Π)⁺: **0 de 67** (máx del cociente
  |E|/(Ŝ−Π)⁺ = 0,9930 < 1). Con E ≤ 0, δ(r) > 0 en **todas** ellas.
- E = 0 **exacto** en 39 de 67 (58,2 %): el desplazamiento ni siquiera
  actúa sobre esos conductores, y el exceso del modelo pasa íntegro a δ.
- Condicionado a E ≠ 0 (28 instancias): el cociente tiene mediana
  0,6271 — E recorta típicamente un ~63 % del exceso, nunca todo.

Es decir: en ventana, E **afeita** el exceso del modelo; no lo cancela.
Si (BAND-δ-HI) es cierta, no será por cancelación r a r, sino porque el
exceso no cubierto ((Ŝ−Π)⁺ − |E|)⁺, positivo punto a punto, tenga masa
sumable sobre la capa — una afirmación sobre la *distribución* del
exceso no cubierto entre conductores, no sobre su anulación.

**Parte degenerada.** Los degenerados de la capa (V̂ = {1}, E = 0,
δ = 1 − Π) no pertenecen a la pregunta de cancelación: su masa es la
contabilidad de (SLAB-V′). En los marcos A el máximo de δ⁺ de la capa
coincide con el de su parte no degenerada; en el marco B la parte no
degenerada solo llega a 0,0037 frente a 0,0128 del total, así que ahí
los degenerados cargan la mayor parte de la masa de capa. Estrictamente, (BAND-δ-HI) = parte degenerada
(territorio SLAB) + el enunciado de la §2.

## 2. Estructura: quién cancela, y el enunciado mínimo

**E es 100 % lado Add en la capa [M].** Descomposición de E por el
clasificador del doc 114 (`lados`) en 43 instancias no degeneradas de
la capa (identidad E = Σ_Add a − Σ_Rem a verificada a 2,2·10⁻¹⁶):

- masa |E_Rem| = 0,0000: el lado Rem no aporta **nada** en la muestra
  (Σ_Rem a = 0 en las 43; ningún signo anómalo en ningún lado).
- toda la masa de |E| es |E_Add| (100,0 %); cociente mediano
  |E_Add|/(Ŝ−Π)⁺ = 0,5581 en esta muestra; cobertura puntual de Add
  solo: 0/43; de E total: 0/43 (consistente con K1).

Esto es coherente con el doc 114 (masa global de Rem⁻ = 0,0288 contra
18,41 de Add⁺), pero en la capa alta es extremo: Rem está vacío como
fuente de cancelación. La demanda de (BAND-δ-HI) recae entera sobre el
flujo Add — el mismo objeto de (FLOW-COMP)/(FC-2), como anunció el
doc 119, ahora con el lado Rem eliminado de la discusión.

**Enunciado mínimo que la medición apoya.** Con la identidad de flujo
del doc 114, |E_Add(r)| = Σ_{y∈Add⁻(r)} |a(y)|(1 − R_A(y)) [T como
identidad], se formula:

**(BAND-δ-HI′)** [P]: Σ_{r∈banda, u_r ≥ u₀*} ((Ŝ−Π)(r)⁺ −
|E_Add(r)|)⁺/φ(r) ≤ C·Q_c,

declarado módulo dos provisiones, ambas [M] en ventana y ninguna [T]:
(i) E_Rem(r) = 0 sobre los conductores de banda de la capa;
(ii) E_Add(r) ≤ 0 ahí (la restricción de (DESP-SIGN) a la capa).
Bajo (i)+(ii), (BAND-δ-HI) ⟺ (BAND-δ-HI′) [T, álgebra puntual].
(BAND-δ-HI′) ya no menciona cancelación: pide que las masas de déficit
de flujo (1 − R_A(y)) dominen en suma, sobre la capa, al exceso del
modelo (Ŝ−Π)⁺(r) — dos objetos de naturaleza distinta (flujo de
apareamiento contra término principal W(u) − 1) comparados conductor a
conductor. La medición dice que hoy, en ventana, el flujo pierde
punto a punto (cociente ≤ 0,993, mediana ~0,6 cuando actúa, 0 en el
58 % de los casos) y gana solo en que la masa total resultante es
≤ 0,1 Q_c — la misma invisibilidad log log de los docs 118–119.

## 3. Conexión con (DESP-SIGN) y el margen de (FLOW-COMP)

El margen ×28 del doc 117a es una medición [M] en **una** instancia
(z = 331, r = 5, L = 79051268): GOOD = 0,030399 contra BAD = 0,001086
en la desigualdad (FLOW-COMP) para ese conductor. Certifica que
E_Add(5) es negativo con holgura **ahí**. Sobre la cancelación de la
capa alta no dice nada, por tres razones independientes:

1. **Distinto conjunto de conductores.** (FLOW-COMP) se evalúa en los
   conductores del argumento de (DESP-SIGN); (BAND-δ-HI) suma sobre
   conductores de banda (Π(r) < 1/2) con u_r ∈ [1, 3/2) en niveles
   u_L → 2. El margen es local a (r, marco); no es una cota uniforme.
2. **Distinta comparación.** (FLOW-COMP) compara dos partes de la masa
   Add entre sí (buena contra mala) y concluye un **signo** con margen.
   El exceso del modelo (Ŝ−Π)⁺ no aparece en (FLOW-COMP). (BAND-δ-HI)
   pide una **magnitud**: |E_Add(r)| contra (Ŝ−Π)⁺(r) ≍ Π(r) sobre las
   familias obstructoras. Signo con margen ≠ tamaño suficiente.
3. **La ventana separa ambas demandas empíricamente.** El mismo
   E_Add que tiene margen ×28 en r = 5 es exactamente 0 en el 58 % de
   los conductores de banda de la capa con exceso positivo (§1). Las
   dos afirmaciones coexisten sin tocarse.

Conclusión precisa: si (DESP-SIGN) se demostrara mañana con ese
margen, (BAND-δ-HI) ganaría solo la provisión (ii) de la §2 — el signo
que permite escribir δ⁺ = ((Ŝ−Π) − |E|)⁺ — y quedaría **exactamente
igual de abierta** como (BAND-δ-HI′): la comparación de tamaños entre
el déficit de flujo y el término principal W(u) − 1 es una demanda
independiente sobre E, no una consecuencia de su signo.

## 4. Veredicto

**(BAND-δ-HI): ni refutada ni apoyada como cancelación; reformulada.**

- No refutada [M]: Σ δ⁺/φ sobre la capa es ≤ 0,093 Q_c en todos los
  marcos medidos, y la familia decisiva ℱ está ausente de la banda en
  ventana — el enunciado asintótico queda intacto y abierto.
- Muerta como mecanismo [M]: la cancelación puntual de E sobre
  (Ŝ−Π)⁺ no ocurre en ventana (0/67; E = 0 en el 58 % de los casos;
  cobertura máxima 99,3 %). Toda futura prueba deberá acotar el
  exceso no cubierto directamente, no invocar cancelación.
- Reducida [T, módulo (i)+(ii) de la §2] a **(BAND-δ-HI′)**: el
  déficit del flujo Add contra el exceso del modelo, suma sobre la
  capa — con el lado Rem eliminado [M] y la parte degenerada
  derivada a la contabilidad (SLAB-V′).
- El margen ×28 de (FLOW-COMP) no transfiere: demandas independientes
  sobre E (signo agregado contra magnitud por conductor).

Inventario de (H2) sin cambios de cardinalidad: (FC-2),
(BAND-δ-HI) — ahora leída como (BAND-δ-HI′) —, (BAND-LOW-NONDEG),
(SLAB-V′). **Goldbach permanece abierto.**

## 5. Reproducibilidad

Script `120a_banda_delta_hi.py` (solo stdlib; importa la maquinaria de
los docs 110/114/117b sin tocarla), salida
`salida_doc120a_local.txt`.
