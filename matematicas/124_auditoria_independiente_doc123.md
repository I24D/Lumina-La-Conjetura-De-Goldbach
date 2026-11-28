# Doc 124 — Auditoría independiente de Doc 123: POS, concentración conjunta y alcance de OUT-MASS

**Autora:** Lumina Chat  
**Fecha:** 2026-11-28  
**Línea:** UNIF-U / (H2)  
**Fuentes auditadas:** 123_wit_outmass.md, 123a_wit.md, 123b_outmass.md, 123a_wit.py, 123b_outmass.py.  
**Verificación adicional:** 124_auditoria_doc123.py.  
**Etiquetas:** [T] = resultado analítico; [M] = medición reproducida; [P] = pendiente.

---

## 0. Veredicto

La auditoría independiente confirma los dos resultados computacionales centrales de Doc 123: (WIT) tiene contraejemplos reales en los estratos 7 y 11; y las identidades estructurales y la fórmula de la ventana prima de (OUT-MASS) se reproducen en los marcos documentados.

Se encontraron tres precisiones que deben incorporarse al registro:

1. **El caso del estrato 5 de (WIT) sí queda cerrado por (POS).** La supuesta posibilidad residual \(c=11\) de 123a_wit.md proviene de usar dos primos consecutivos iguales en la cota penúltima. Como \(p_{k-1}>p_k\), ese caso es imposible. Para \(q_0=5\), (POS) implica \(c\in\{1,5,7\}\), luego WIT sí vale en ese estrato. WIT global continúa refutado por los estratos 7 y 11.
2. **El 99,9 % corresponde a la condición marginal de conductor primo, no a su intersección con la banda de \(u_r\).** La intersección reproducida es 92,6537 % de la masa pooled.
3. **La implicación de (OUT-MASS) hacia la cota de \(\delta^+\) sigue dependiendo de \(E_{\rm Rem}=0\) en la capa.** Esa provisión está medida en la ventana anterior, pero no está demostrada uniformemente. No debe elevarse de [M] a [T].

Goldbach permanece abierto.

---

## 1. Ejecución independiente

Los archivos y las dependencias de numerico/ se copiaron a un directorio temporal aislado de la PC LUMINA para respetar las rutas relativas de importación. No se modificó ningún archivo fuente del repositorio durante la ejecución.

**Entorno:** Windows, Python 3.14.10, dispositivo Desktop Commander Lumina.

| Programa | Chequeos | Resultado | Tiempo | Código de salida |
|---|---|---:|---:|---:|
| 123a_wit.py | K0–K3 | PASS | 4,77 s | 0 |
| 123b_outmass.py | K0–K4 | PASS | 90,88 s | 0 |
| 124_auditoria_doc123.py | POS exacto, concentración conjunta y cinco casos K4 | PASS | 89,16 s | 0 |

Los scripts de Doc 123 se ejecutaron desde el directorio numerico con sus dependencias originales: 109, 110, 117b, 118b1, 120a y 122b. El script 124 vuelve a calcular el conjunto pooled de K2 y añade una prueba explícita de la masa cercana compuesta en los cinco casos K4.

### 1.1 Resultados reproducidos de 123a

- K0: los conjuntos de cofactores de los marcos originales se reproducen exactamente.
- K1: el evaluador ligero coincide con la maquinaria real a \(z=53\), con error máximo de \(\Pi\) de \(4,8\cdot10^{-16}\).
- Estrato 7, marco A, \(z=61\): \(L=39179\), \(r=115=5\cdot23\), testigos \(\{1,7,11,13\}\), \(N(r)=4\).
- Estrato 7, marco A, \(z=67\): reaparece el patrón con \(L=28942\), \(r=85\), \(N(r)=4\).
- Estrato 11, marco A, \(z=127\): \(L=127^4=260144641\), \(r=364105=5\cdot7\cdot101\cdot103\), testigos \(\{1,11,13,17,19,23\}\), \(N(r)=6\).
- K3: la violación del estrato 11 se reproduce en 17 niveles.

Esto confirma la refutación de WIT como lema uniforme. Los contraejemplos no son artefactos de una única elección de nivel.

### 1.2 Resultados reproducidos de 123b

- K0: \(\sum S(r)/\varphi(r)\) y \(\sum\Pi(r)/\varphi(r)\) se imprimen como 1; el residuo numérico de \(\sum\delta(r)/\varphi(r)\) es del orden de \(10^{-16}\). La identidad exacta se apoya en la prueba algebraica de §3 de 123b_outmass.md, no en redondear la salida flotante.
- K0: la identidad \(\delta=E_{\rm Rem}+M-P\), la equivalencia de fibras \(\widehat V^{(r)}=V^+(L/r)\) y la forma cerrada de \(M_{\rm far}\) se reproducen sin discrepancias significativas; no aparece violación de la cota puntual.
- K1: las cuatro anclas anteriores se reproducen: 0,0533; 0,1956; 0,0059; 0,3152.
- K2: se reproducen los máximos \(C(z)\): 0,1956; 0,3152; 0,1955; 0,2089; 0,2598 para \(z=37,53,61,67,71\).
- K4: en las cinco instancias primas, los conjuntos \(\mathcal Q(p_0)\) coinciden con los primos cercanos ausentes. La masa cercana con \(\omega(d')\ge3\) sale cero en los cinco casos.

La identidad global con signo es válida: \(\sum_r S(r)/\varphi(r)=1\) por la suma de Möbius sobre los divisores de cada \(d\), y \(\sum_r\Pi(r)/\varphi(r)=1\) por factorización local de Euler. Por tanto \(\sum_r\delta(r)/\varphi(r)=0\) exactamente.

---

## 2. Corrección de la clasificación POS en el estrato 5

La familia correcta es, para **cada posición impar** \(m\) de la lista fusionada decreciente,

\[
c\,p_m^2<q_0^3\prod_{j>m}p_j.
\]

La palabra «impar» debe conservarse también en el resumen Doc 123. No se obtiene una desigualdad para toda posición.

La cota de última posición

\[
c<\frac{q_0^3}{p_k^2}
\]

solo corresponde a \(k\) impar. Si \(k\) es par, la última posición impar es \(k-1\), y la cota correcta es

\[
c<\frac{q_0^3p_k}{p_{k-1}^2},
\qquad p_{k-1}>p_k.
\]

### 2.1 Cierre de WIT para \(q_0=5\)

Como \(q_0=5\), el primo 5 no divide a \(r\); todos los primos de \(c\) son al menos 5.

**Si \(k\) es impar:**

- Si \(5\mid c\), entonces \(p_k=5\) y \(25c<125\), imposible para \(c\ge5\).
- Si \(5\nmid c\), entonces \(p_k\ge7\) y \(49c<125\), de donde \(c=1\).

**Si \(k\) es par**, se usa la posición \(k-1\):

- Si \(p_k=5\), entonces \(5\mid c\), \(p_{k-1}\ge7\) y
  \[
  c<\frac{625}{49}<13.
  \]
  Como \(c\) es libre de cuadrados y divisible por 5, queda \(c=5\).
- Si \(p_k=7\), entonces \(p_{k-1}\ge11\) y
  \[
  c<\frac{875}{121}<8.
  \]
  Quedan \(c=1\) o \(c=7\). Para que \(c=7\), el 7 debe proceder de \(c\), y la segunda prima debe ser 11; si \(p_{k-1}\ge13\), la cota baja de 8.
- Si \(p_k\ge11\), entonces \(p_{k-1}\ge13\) y
  \[
  c<\frac{125\cdot11}{13^2}
   =\frac{1375}{169}<9.
  \]
  Como todos los primos de \(c\) son al menos 11 en este caso, queda \(c=1\).

Por tanto, en ambas paridades,

\[
\boxed{q_0=5\quad\Longrightarrow\quad c\in\{1,5,7\}.}
\]

Así \(N(r)\le3\) y WIT queda demostrado para el estrato 5. La posibilidad residual \(c=11\) de 123a_wit.md no existe: requeriría \(p_{k-1}=p_k=11\), contrario al orden estricto de la lista. Los límites racionales exactos están reproducidos por 124_auditoria_doc123.py.

### 2.2 Exclusión de cofactores compuestos en los tres estratos

Las cotas de POS también justifican que \(c>1\) sea primo:

- \(q_0=7\): para \(k\) impar, \(c<343/25=13,72\); para \(k\) par, \(c<343\cdot5/49=35\). El menor cofactor compuesto admisible sería \(7\cdot11=77\).
- \(q_0=11\): para \(k\) impar, \(c<1331/25=53,24\); para \(k\) par, \(c<1331\cdot5/49<135,82\). El menor cofactor compuesto admisible sería \(11\cdot13=143\).
- \(q_0=5\): queda la clasificación más fuerte \(c\in\{1,5,7\}\).

El resultado correcto es: \(c=1\) o \(c\) primo. Esto no salva WIT global: los primos 13, 17, 19 y 23 ya producen los contraejemplos documentados.

---

## 3. Auditoría de la concentración de OUT-MASS

La salida K3 de 123b_outmass.py informa dos distribuciones marginales: masa por \(\omega(r)\) y masa por banda de \(u_r\). El resumen Doc 123 y la formulación transmitida mezclan esas dos marginales al atribuir el 99,9 % a su intersección.

El script 124 recalcula la intersección directamente sobre las 2127 instancias pooled de K2, con el mismo peso \(M(r)/(\varphi(r)Q_c)\).

| Condición | Masa pooled | Fracción del total |
|---|---:|---:|
| Total | 4,072038624194 | 100 % |
| \(\omega(r)=1\) | 4,067937062993 | 99,899275 % |
| \(1,35\le u_r<1,5\) | 3,775420378432 | 92,715731 % |
| Ambas condiciones | 3,772894930921 | **92,653712 %** |

**Corrección:** 99,9 % es la masa en conductores primos, sin imponer la banda. La masa que satisface simultáneamente «conductor primo» y \(u_r\in[1,35,3/2)\) es 92,6537 % en esta muestra pooled.

Esto sigue siendo una medición finita de los marcos y niveles de K2; no es una afirmación asintótica uniforme.

---

## 4. Alcance exacto de la fórmula PRIME-WINDOW

Para un conductor primo \(r=p_0\), la fórmula

\[
M_{\rm near}^{(1)}(p_0)
=
\sum_{q\in\mathcal Q(p_0)}\frac1{q-1},
\qquad
\mathcal Q(p_0)=
\{q>p_0:\ q^3\ge L,\ p_0q<L\},
\]

es exacta para la componente de ausentes cercanos con \(d'=q\) primo.

La igualdad con la masa cercana total \(M_{\rm near}\) requiere además que no haya ausentes cercanos con \(\omega(d')\ge3\). 123b_outmass.md lo demuestra por tamaño en las cinco instancias K4: allí \(L/r<455\), mientras que 455 es una cota inferior para el producto de tres primos permitidos coprimos con \(r\).

La ejecución reproduce masa compuesta cercana cero en los cinco casos:
\((37,1286,7)\), \((53,2585,7)\), \((61,3450,11)\), \((67,4024,11)\), \((71,4421,11)\).

**Precisión de alcance:** la fórmula para \(M_{\rm near}^{(1)}\) es general para la componente prima; \(M_{\rm near}=M_{\rm near}^{(1)}\) está probado por tamaño en esas cinco instancias, no para todos los conductores primos de la banda asintótica. La constante \(C(z)\) oscila en la ventana medida \(37\le z\le71\); esto no demuestra que esté uniformemente acotada al crecer \(z\).

**Mejora de reproducibilidad:** K4 en 123b_outmass.py compara el conjunto de primos \(\mathcal Q(p_0)\) con los primos ausentes cercanos, pero no afirma explícitamente \(M_{\rm near}(\omega(d')\ge3)=0\). El script 124 añade esa aserción en los cinco casos K4.

---

## 5. Dependencia de \(E_{\rm Rem}\) en la ruta de signo

La identidad puntual es

\[
\delta(r)=E_{\rm Rem}(r)+M(r)-P(r).
\]

Por tanto,

\[
\sum_{\rm capa}\frac{M-P}{\varphi}
=
\sum_{\rm capa}\frac{\delta}{\varphi}
-
\sum_{\rm capa}\frac{E_{\rm Rem}}{\varphi}.
\]

La identidad global \(\sum_{\rm enum}\delta/\varphi=0\) permite sustituir la suma de la capa por menos la del complemento, pero no elimina el término \(E_{\rm Rem}\).

Asimismo, la deducción puntual \(\delta^+(r)\le M(r)\), usada para que (OUT-MASS) implique la cota de \(\delta^+\), necesita \(E_{\rm Rem}(r)=0\) en la capa (o una cota alternativa que controle su contribución positiva).

El doc 122a mide \(E_{\rm Rem}=0\) en sus 67 instancias de capa; esa provisión sigue siendo [M], no [T] uniforme. K0 de 123b_outmass.py verifica la descomposición \(\delta=E_{\rm Rem}+M-P\), pero no demuestra por sí solo \(E_{\rm Rem}=0\) para todos los niveles.

**Consecuencia para el inventario:** NET-LAYER + P-MASS no debe presentarse como cierre incondicional de la ruta con signo sin incluir el término \(E_{\rm Rem}\). La ruta puntual PRIME-WINDOW + R1–R3 también necesita conservar la provisión \(E_{\rm Rem}=0\) al transferir (OUT-MASS) a (BAND-\(\delta\)-HI′).

---

## 6. Inventario corregido de H2

- **(FC-2):** sin cambio; permanece [P].
- **(OUT-MASS):** estructura e identidad con signo demostradas; reducción a (NET-LAYER) + (P-MASS), o a (PRIME-WINDOW) + R1–R3, sin cierre uniforme. Mantener explícita la provisión \(E_{\rm Rem}=0\) para la transferencia a \(\delta^+\).
- **(BAND-LOW-NONDEG):** sin cambio; [P].
- **(SLAB-V′):** sin cambio de enunciado; [P]. La reducción a tres cofactores queda descartada. (POS) controla los cofactores: \(c=1\) o primo, con límites posicionales. WIT global sigue refutado, pero el estrato 5 queda demostrado por POS.
- **(DIST-\(\mathcal D\)):** sin cambio respecto de Doc 122 y su balance de factibilidad; permanece abierto.

Los resultados de Doc 123 estrechan los residuos y corrigen una parte del análisis del estrato 5, pero no cierran H2 ni la cadena UNIF-U-COND.

**Goldbach fuerte permanece abierto.**

---

## Reproducibilidad

- 123a_wit.py: PASS, K0–K3, ejecutado independientemente.
- 123b_outmass.py: PASS, K0–K4, ejecutado independientemente.
- 124_auditoria_doc123.py: PASS; reproduce las estadísticas conjuntas, las cotas racionales de POS para \(q_0=5\) y la ausencia de masa cercana compuesta en los cinco casos K4.

No se modificaron los documentos fuente 123a/123b; este archivo registra las correcciones y el alcance auditado.
