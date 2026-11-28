# Doc 116 — Los tres residuos del doc 114, resueltos: uno cerrado, dos refutados en su forma puntual

**Fecha:** 2026-11-28. **Línea:** UNIF-U / (H2), continuación de los docs 113–115.
**Documentos de prueba (este doc es la síntesis; las pruebas completas viven en ellos):**
`116a_rem_pers.md` (script PASS, J0–J5, 2,4 s), `116b_flujo_afinado.md` (script PASS, J0–J4, 1,1 s), `116c_banda_conteo.md` (script PASS, J0–J5, 0,7 s).
**Etiquetas:** [T] = teorema, [M] = medición en ventana, [P] = pendiente.

---

## 0. Veredicto en breve

El doc 114 redujo (DESP-SIGN) a **(R-FLOW-71+)** + **(REM-PERS-small)** y la banda a **(BAND)**. Este doc ataca los tres, y el resultado es asimétrico:

1. **(REM-PERS-small) CERRADO [T]** (§1, doc 116a): la familia residual es vacía siempre, no solo en ventana. (DESP-SIGN) pierde uno de sus dos residuos por completo.
2. **(R-FLOW-71+) es FALSO como enunciado puntual [T]** (§2, doc 116b): hay contraejemplos dentro del marco, verificados con los predicados de la serie. La ruta del doc 114 a (DESP-SIGN) vía flujo puntual está muerta. **(DESP-SIGN) en sí NO queda refutado**: hereda el residuo como compensación global, **(FLOW-COMP)** [P].
3. **(BAND) con la cota cruda es FALSA en el marco asintótico [T]** (§3, doc 116c): obstrucción exacta por conductores primos. El residuo mínimo que sobrevive es **(BAND-δ)** [P], la forma en media L¹.

**(H2) sigue abierta.** Lo que queda entre (H2) y la cadena de los docs 109–113 son dos enunciados **promediados**, no puntuales: (FLOW-COMP) y (BAND-δ). Goldbach no se declara demostrado.

---

## 1. (REM-PERS-small): cerrado [T] (doc 116a)

**Teorema (116a).** Todo fallo fusionado de todo x ∈ Rem⁻ cae en un primo > q_min(x). En consecuencia la persistencia del doc 114 §4 aplica sin excepciones, ψ(x) = x/q_min(x) ∈ Rem⁺ siempre, y la contribución de la familia residual a E_Rem es exactamente 0.

*Prueba en dos lemas.* (i) q_min(x) no puede fallar: su posición propia es ω(x), y el teorema R1 del doc 114 exige que un primo propio que falla ocupe posición par propia; el análisis de paridad de la posición de q_min lo excluye. (ii) Ningún primo p | r con p < q_min(x) puede fallar por tamaño: la condición del modelo en la última posición propia da r·x·q_min(x)² < L, y un fallo fusionado en p obligaría q_min(x) < p, contradicción. **Corolario [T]:** Rem no contiene primos (el caso ω(x) = 1 queda excluido por (i)–(ii)).

*Verificación:* las 4298 instancias del doc 114 — Rem⁻ = 39 + 36, aterrizajes de ψ 39/39 y 36/36, 0 fallos en primos ≤ q_min, razón máxima del lema de tamaño 0,7581 < 1.

## 2. (R-FLOW-71+): falso puntualmente [T]; (DESP-SIGN) hereda (FLOW-COMP) (doc 116b)

**Contraejemplos dentro del marco [T].** R_A(y) = 1,00220 con min(y) = 251 (r = 5, y primo, z = 331, L ∈ ℒ_A, r visible): es el primer cruce de una familia con R_A(y) ~ log log Y → ∞, luego el enunciado puntual es falso para min(y) grande, no solo no probado. Y R_R(y) = Σ_{5≤p≤251} 1/(p−1) = 1,25620 con min(y) = 4099 (z = 6052, s = 5).

**Las cotas afiladas existen y caracterizan la obstrucción [T].** La cota cúbica (★) por la posición fusionada de q, el volteo (★★) y su caso (★★★) q < L/(ρ³y) valen sin excepción (0 violaciones en los 1640 pares del doc 114). R_A(y) > 1 solo puede ocurrir en la **configuración libre**: y primo, r primo pequeño, q en posición par y atado de volteo vacío en rango. Fuera de ella, el flujo puntual sí cierra.

**Por qué (DESP-SIGN) no cae con el enunciado puntual.** E_Add(r) = −Σ_y |a(y)|(1 − R_A(y)): lo que (DESP-SIGN) necesita es el signo de la SUMA, no de cada término. En la configuración libre, |a(y)| = 1/(y−1) con y ≥ 251 es minúsculo, y el doc 114 midió holgura estructural global (masa peligrosa 18,41 contra masa socia 93,29). El residuo exacto que hereda (DESP-SIGN) es:

**(FLOW-COMP) [P]:** para todo r visible, Σ_{y∈Add⁻} |a(y)|(1 − R_A(y)) ≥ 0 y el análogo en Rem, usando que los términos con R > 1 solo viven en la configuración libre caracterizada por (★)–(★★★).

*Medición [M]:* a z = 101, máx R_A crece con min(y): 0,8982 (min 71) → 0,9514 (min 89/97), siempre < 1 en ventana; a z = 71, máx 0,8493. La falsedad puntual es un fenómeno de rango mayor que el medido, detectado por construcción teórica y verificado con los predicados exactos.

## 3. (BAND) cruda: falsa asintóticamente [T]; residuo (BAND-δ) (doc 116c)

**Obstrucción [T].** La suma de la cota cruda sobre la banda, Σ_{banda} (1−Π(r)²)/φ(r), no es O(Q_completo): los conductores **primos** del slice superior aportan liminf ≥ B₁(u*) = log(min(1, u*/3)/(u* − u₀)) > 0 — B₁(5/3) = log(10/3) = 1,2040 y B₁(2) = log(4/3) = 0,2877, verificados por sumas de Mertens hasta Z = 10⁸/10¹² — mientras Q_completo ~ K·M_∞/ln z → 0. Con cofactores fijos y cajas explícitas, la obstrucción cubre todo ω ≥ 1 y todo u* ∈ [5/3, 4]. La cota cruda pierde porque sustituye δ(r) por 1 − Π(r): en los conductores primos de la banda, δ(r) real es mucho menor que la cota, pero el conteo no lo ve.

**Caracterización [T].** Π(r) = Π_full / Π_{p|r}(1 − 1/(p−1)) ≥ Π_full (la división es la dirección correcta); Π(r) < 1/2 ⟺ Σ_{p|r} 1/(p−1) ≲ log(1/(2Π_full)). La banda es vacía a z = 17 (Π_full = 0,5156 > 1/2) y asintóticamente todo primo entra en ella. **Dicotomía exacta:** el conteo cierra una región ⟺ su región límite de Mertens tiene medida cero; la de la banda tiene medida positiva en todo nivel. Lo único cerrado por conteo sigue siendo {Π ≥ c₀} [T] (doc 114; reverificado ≤ 0,43·Q_c medido). Validación cruzada: el máximo muestreado 0,4408 (z = 37) coincide con la columna Π < 1/2 del doc 114.

**Residuo mínimo exacto. (BAND-δ) [P]:** Σ_{r en banda} δ(r)⁺/φ(r) ≤ C·Q_completo. Implica la forma sumable de la banda con factor 3/2 [T, módulo S(r) ≤ 1 del doc 114 §6]. (BAND-δ) no es conteo: es (FL-DESP) en media L¹ sobre u_r < u₀ — pide que los defectos positivos, promediados con peso 1/φ sobre la banda, sean O(Q_c), aunque puntualmente no lo sean.

## 4. Veredicto de cadena actualizado

- Doc 114: (DESP-SIGN) ⟸ (R-FLOW-71+) + (REM-PERS-small); banda ⟸ (BAND).
- **Doc 116:** (REM-PERS-small) **cerrado [T]**. (R-FLOW-71+) **falso puntualmente [T]** → sustituido por **(FLOW-COMP)** [P]. (BAND) cruda **falsa [T]** → sustituida por **(BAND-δ)** [P].
- **Teorema condicional actualizado:** (FLOW-COMP) + (BAND-δ) ⟹ (DESP-SIGN) + banda ⟹ (FL-DESP⁺) ⟹ (H2⁗) ⟹ (H2), con la constante ≈ 6,4 del doc 113 más la constante de (BAND-δ). Y con el doc 115: **(DIST-𝒟) + (H2) ⇒ UNIF-U-COND**, con (H2) dependiendo ahora exactamente de esos dos enunciados promediados.
- Lectura honesta: las formas puntuales/de conteo de ambos residuos eran falsas; las formas promediadas son estrictamente más débiles y compatibles con toda la medición (4298 instancias sin contraejemplo de (DESP-SIGN); masa real de banda ≤ 0,3806·Q_c). El patrón se repite en toda la línea (H2): cada vez que un enunciado puntual muere, el enunciado promedio correcto es el que la medición ya veía. (FLOW-COMP) y (BAND-δ) son el frente exacto que queda en la línea de coeficientes.

Siguen fuera, sin cambios: (DIST-𝒟) en la línea de Lumina Chat, y la separación de capas G0. **Goldbach permanece abierto.**

## 5. Reproducibilidad

Las pruebas y verificaciones viven en `116a_rem_pers.py`, `116b_flujo_afinado.py` y `116c_banda_conteo.py` (stdlib; importan la maquinaria de los docs 109–114), con salidas `salida_doc116a_local.txt`, `salida_doc116b_local.txt` y `salida_doc116c_local.txt`, las tres en PASS. Este doc no introduce cómputo nuevo: es la síntesis y el veredicto de cadena de las tres secciones.
