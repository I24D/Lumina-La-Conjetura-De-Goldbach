# Doc 122 — Balance de factibilidad de D-NP con el conteo asintótico verdadero de \(\mathcal D_S\)

**Autora:** Lumina Chat  
**Fecha:** 2026-11-28  
**Base documental:** docs 108, 115, 120b, 121, 121b y \`verificacion-conteo-ds-claude-code.md\`.  
**Objeto:** cerrar el balance de factibilidad de la reducción CRT de D-NP antes de proponer otra desigualdad.  
**Alcance:** análisis de las cotas estándar y de la reindexación por conductores; no es una demostración de (DIST-\(\mathcal D\)) ni de Goldbach.

## 0. Veredicto ejecutivo

La aplicación de Baier–Zhao (BZ) planteada para D-NP no es viable en su forma actual. El dato correcto es

\[
|\mathcal D_S|=X^{1/2+o(1)},
\]

no el exponente \(0.4211\), que correspondía a \(|\mathcal H|\) en la ventana pequeña y saturada. La medición superior de Claude Code da \(\alpha=0.4659\), intervalo de ajuste \([0.4649,0.4669]\), hasta \(X=3.2\cdot10^{28}\); el documento también da una demostración del exponente asintótico \(1/2\).

La condición \(\alpha<0.420232\) era el requisito de la aplicación concreta de BZ que se venía intentando. El documento de conteo advierte que la derivación original de ese umbral no está archivada allí; por tanto, se conserva como requisito de esa ruta, no como hipótesis literal del Teorema 1 de BZ. La diferencia respecto del exponente asintótico es

\[
\frac12-0.420232=0.079768,
\qquad
X^{0.079768}\big|_{X=3.2\cdot10^{28}}\approx187.845.
\]

La ruta no debe reintentarse con el mismo balance de conteo.

Hay además una segunda obstrucción: el peso inducido por la reducción CRT no desaparece al agrupar caracteres por conductor primitivo. La suma de pesos sobre los módulos que inducen un mismo conductor tiene una fórmula exacta y crece polinómicamente en una subfamilia. La unicidad de rama evita duplicar el primo de conmutación dentro de un conductor; no elimina los múltiples levantamientos \(m\) del mismo carácter primitivo.

**Conclusión:** la criba grande clásica y el Lema 8 de BZ, aplicados mediante el mayorante uniforme del peso inducido, no alcanzan \(T(X)\). Esto descarta esa aplicación genérica; no demuestra que el momento real D-NP sea grande ni que (DIST-\(\mathcal D\)) sea falso.

---

## 1. Definiciones y objetivo exacto

Usamos las definiciones fijadas por los docs 115, 120b y 121b:

\[
\Omega_{qm}^{C}=\frac{\widetilde\Theta_m(q)}{q-1},
\qquad
\widetilde\Theta_m(q)=
\sum_{\substack{h\in V^+(D_A/q)\setminus\{1\}\\m\mid h}}
\frac{2\tau(h)-1}{\varphi(h)}.
\]

Para \(d=qm\in\mathcal D_S\), la reducción CRT corregida de D-NP tiene el peso

\[
W_d=\Omega_d^C\varphi(d),
\]

y el objetivo es acotar

\[
\mathcal M_W(X)=
\sum_{d=qm\in\mathcal D_S}
W_d
\sum_{\substack{\eta\bmod d\\\eta_q\ne\chi_0}}
|\Psi_X(\eta)|^2
\]

por la escala

\[
T(X)=\frac{X^{2-2c_*}}{0.0225(\log X)^4},
\qquad c_*=0.039884,
\]

es decir,

\[
T(X)=\frac{X^{1.920232}}{0.0225(\log X)^4}.
\]

La forma \(W_d=\Omega_d^C\varphi(d)\), y no \(\Omega_d^C/\varphi(d)\), es obligatoria: se sigue de \(\widetilde\Theta_m(q)=(q-1)\Omega_{qm}^C\) y \(\varphi(qm)=(q-1)\varphi(m)\). Véase doc 121b, §§2–3.

La pérdida de Cauchy–Schwarz medida en doc 121b es \(\lambda=2.40,4.06,11.64,14.52\) en los casos de prueba. No hay cota uniforme superior demostrada. Si \(M_W\) denota el momento después de Cauchy, entonces \(M_W=\lambda D_{\rm NP}\) en la normalización de esa auditoría. Para demostrar \(D_{\rm NP}\le T(X)\), basta demostrar \(M_W\le T(X)\); **no se multiplica otra vez por \(\lambda\)**. La medición de \(\lambda\) cuantifica el sobreconteo introducido por la reducción, pero no autoriza dividir por 2.4 ni por 14.5 en una prueba uniforme. En las cotas siguientes se usa únicamente \(\lambda\ge1\).

## 2. Conteo verdadero de \(\mathcal D_S\)

El informe de Claude Code establece:

\[
|\mathcal D_S|=X^{1/2+o(1)}
\]

y también

\[
\sum_{d\in\mathcal D_S}W_d=X^{1/2+o(1)}.
\]

En la ventana \(z=307,\ldots,709\), \(X=z^{10}\) llega a \(3.2\cdot10^{28}\), con \(|\mathcal D_S|\) hasta \(8.61\cdot10^{10}\). El ajuste de la ventana superior da \(\alpha=0.4659\), con intervalo de regresión \([0.4649,0.4669]\). Los puntos son deterministas: el intervalo describe la dispersión del ajuste, no incertidumbre estadística sobre el límite.

El valor \(0.4211\) era el exponente de \(|\mathcal H|\) en \(z=11,\ldots,37\), no el de \(|\mathcal D_S|\). En esa misma ventana, el exponente medido de \(|\mathcal D_S|\) era \(0.4444\). La prueba del §4 del informe fija el límite asintótico en \(1/2\).

Si la aplicación concreta de BZ exige \(\alpha<0.420232\), el requisito falla por \(0.079768\) en el exponente. A \(X=3.2\cdot10^{28}\), el factor de cardinalidad correspondiente es aproximadamente \(187.845\). Este factor expresa la brecha entre el umbral de esa ruta y el conteo real; no debe confundirse con el cociente final entre una cota de gran criba y \(T(X)\).

**Precisión de alcance:** el propio informe de conteo señala que el desarrollo original del umbral \(0.420232\) no está archivado en el repositorio. Por tanto, este documento acepta ese umbral como requisito de la aplicación concreta previamente propuesta, pero no lo atribuye al enunciado literal del Teorema 1 de BZ.

## 3. Gran criba clásica: escala sin el peso inducido

Para \(a_n=\Lambda(n)\),

\[
\sum_{n\le X}|a_n|^2\asymp X\log X.
\]

La gran criba clásica para caracteres primitivos da

\[
\sum_{r\le Q}\frac r{\varphi(r)}
\sum_{\chi\bmod r}^{*}
\left|\sum_{n\le X}a_n\chi(n)\right|^2
\ll (Q^2+X)\sum_{n\le X}|a_n|^2.
\]

Con \(Q\asymp D_A\asymp X^{1/2}\),

\[
Q^2+X\asymp X,
\qquad
\text{momento primitivo sin pesos}\ll X^2\log X.
\]

El cociente de esta escala respecto de \(T(X)\) es

\[
\frac{X^2\log X}{T(X)}
=
0.0225\,X^{0.079768}(\log X)^5.
\]

En \(X=3.2\cdot10^{28}\), el factor explícito de esta comparación es aproximadamente

\[
5.15\cdot10^9.
\]

Así, incluso antes de introducir \(W_d\), la cota genérica de gran criba no alcanza la escala objetivo.

## 4. Lema 8 de Baier–Zhao

El Lema 8 de BZ (Acta Arith. 125.2 (2006), fórmula (4.2)) afirma, bajo las hipótesis de distribución de su conjunto de módulos, que

\[
\sum_{q\in S_t(Q/t)}\frac q{\varphi(q)}
\sum_{\chi\bmod q}^{*}
\left|\sum_{n=M+1}^{M+N}a_n\chi(n)\right|^2
\ll
\Delta(N,Q,t)\sum_{n=M+1}^{M+N}|a_n|^2,
\]

donde

\[
\Delta(Y,Q,t)=
Y+\frac Qt(QY)^\varepsilon
\left(\sqrt Y+|S_t(Q/t)|\right).
\]

En la sustitución de escala \(t=1\), \(Q\asymp X^{1/2}\), \(|S(Q)|\asymp X^{1/2+o(1)}\), \(Y=X\), el término entre paréntesis es \(X^{1/2+o(1)}\). Por tanto,

\[
\Delta(X,Q,1)=X^{1+o(1)}
\]

y, para \(a_n=\Lambda(n)\), el momento primitivo no ponderado vuelve a ser

\[
\ll X^{2+o(1)}\log X.
\]

No aparece una ganancia polinómica sobre la gran criba clásica.

El Lema 9 es bilineal y controla un primer momento con máximo sobre cortes \(mn\le X\); no sustituye directamente el segundo momento cuadrático de D-NP. Para usarlo aquí haría falta una reducción bilineal adicional, con sus pérdidas explícitas.

Dos cautelas:
1. Los lemas 7–9 presuponen que el conjunto \(S\) es bien distribuido en el sentido definido por BZ. El exponente de conteo no demuestra esa hipótesis.
2. El resultado anterior es la escala de la parte primitiva **sin** los coeficientes inducidos \(A_r\). Esos coeficientes deben incorporarse al reindexar el momento de D-NP.

Referencia primaria: Baier–Zhao, [*Bombieri–Vinogradov type theorems for sparse sets of moduli*](https://www.impan.pl/shop/publication/transaction/download/product/82297?download.pdf=), en particular §4, Lemas 8–9. El Teorema 1 tiene hipótesis propias de buena distribución, regularidad de divisores y cardinalidad; no se debe identificar su condición \(|S(Q)|\gg\sqrt Q\) con el umbral \(0.420232\) de la aplicación concreta.

## 5. Reindexación por conductores primitivos

Un carácter módulo \(d=qm\) con componente no principal en el primo \(q\) tiene conductor \(r=qr_0\), con \(r_0\mid m\). El carácter módulo \(d\) es el levantamiento del carácter primitivo \(\chi^*\bmod r\), restringido a los enteros coprimos con \(d\). Definimos

\[
\Psi_X^{(d)}(\chi^*)=
\sum_{\substack{n\le X\\(n,d)=1}}\Lambda(n)\chi^*(n).
\]

La diferencia entre \(\Psi_X^{(d)}(\chi^*)\) y \(\Psi_X(\chi^*)\) procede de las potencias de primos que dividen \(d/r\); es un término corrector de soporte en esos primos, no una igualdad automática entre sumas inducidas y primitivas.

La unicidad de rama (doc 121a, basada en el Lema 2 del doc 101 §2.1) asegura que el primo de conmutación \(q\) es único para cada conductor del soporte en la rama considerada. No elimina los distintos divisores \(m\) que pueden levantar el mismo carácter primitivo.

Para un conductor \(r=qr_0\), el peso total inducido es

\[
A_{qr_0}=
\sum_{\substack{m:\ r_0\mid m\\qm\in\mathcal D_S}}W_{qm}.
\]

Usando la identidad exacta del doc 115,

\[
\widetilde\Theta_m(q)=
\sum_{\substack{h\in V^+(D_A/q)\setminus\{1\}\\m\mid h}}
\frac{2\tau(h)-1}{\varphi(h)},
\qquad
W_{qm}=\varphi(m)\widetilde\Theta_m(q),
\]

e intercambiando las sumas positivas, resulta

\[
\boxed{
A_{qr_0}=
\frac{\varphi(r_0)}{r_0}
\sum_{\substack{h\in V^+(D_A/q)\setminus\{1\}\\r_0\mid h}}
(2\tau(h)-1)\frac h{\varphi(h)}.
}
\]

La identidad \(\sum_{r_0\mid m\mid h}\varphi(m)=\varphi(r_0)h/r_0\) se usa porque \(h\) es libre de cuadrados y \((r_0,h/r_0)=1\).

En particular, para \(r=q\) (\(r_0=1\)),

\[
A_q=
\sum_{h\in V^+(D_A/q)\setminus\{1\}}
(2\tau(h)-1)\frac h{\varphi(h)}.
\]

Cada sumando es al menos \(3\). La cota inferior asintótica de \(|V^+(L)|\) del informe de conteo, aplicada a primos admisibles \(q\asymp z=X^{1/10}\), con \(D_A\asymp X^{1/2}\), da

\[
A_q\ge X^{2/5-o(1)}
\]

en esa subfamilia. Es una inflación real del coeficiente inducido; no desaparece por la unicidad de rama.

Por otro lado, como \(A_r\) suma un subconjunto de los pesos positivos \(W_d\),

\[
A_r\le \sum_{d\in\mathcal D_S}W_d
=X^{1/2+o(1)}.
\]

Así, para el mayorante por máximo de conductor, el coeficiente

\[
\mathfrak A=\max_r A_r\frac{\varphi(r)}r
\]

satisface, con la información asintótica disponible,

\[
X^{2/5-o(1)}\ \le\ \mathfrak A\ \le\ X^{1/2+o(1)}.
\]

La cota inferior se entiende para la familia \(r=q\) anterior; la superior usa el peso total. No se afirma que el exponente exacto de \(\mathfrak A\) sea uno de los extremos.

## 6. Gran criba después de agregar los pesos

La desigualdad primitiva clásica da

\[
\mathcal M_W(X)
\ll
\mathfrak A\,(Q^2+X)
\sum_{n\le X}\Lambda(n)^2.
\]

Con \(Q\asymp X^{1/2}\),

\[
\mathcal M_W(X)
\ll
\mathfrak A\,X^2\log X.
\]

Usando solo la cota superior demostrada para el peso agregado,

\[
\mathcal M_W(X)\ll X^{5/2+o(1)}\log X.
\]

Frente a \(T(X)\), el cociente de este mayorante es

\[
0.0225\,X^{0.579768+o(1)}(\log X)^5.
\]

En \(X=3.2\cdot10^{28}\), omitiendo el \(o(1)\) y las constantes implícitas, el factor de escala es aproximadamente \(9.21\cdot10^{23}\).

Incluso el extremo inferior \(A_q\ge X^{2/5-o(1)}\) muestra que el método que extrae \(\mathfrak A\) por un máximo no puede tener coeficiente \(X^{o(1)}\). En ese extremo, la escala del mayorante sería como mínimo \(X^{12/5-o(1)}\log X\), cuya comparación con \(T(X)\) contiene el factor

\[
0.0225\,X^{0.479768-o(1)}(\log X)^5.
\]

A \(X=3.2\cdot10^{28}\), la parte explícita sin el \(o(1)\) es aproximadamente \(1.30\cdot10^{21}\).

**Interpretación rigurosa:** estos cocientes comparan los mayorantes producidos por la aplicación estándar de gran criba con \(T(X)\). No son cotas inferiores para el momento verdadero \(\mathcal M_W\). Una desigualdad ponderada que explote la distribución de \(A_r\) sin extraer su máximo podría comportarse de otra manera; esa desigualdad no está disponible en esta etapa.

## 7. Dónde entra \(\lambda\)

La pérdida \(\lambda\) de doc 121b es el cociente entre el momento levantado por Cauchy–Schwarz y el momento original D-NP. La reindexación CRT no pierde caracteres: el conteo de pares es exactamente \(\varphi(m)(q-2)\) en ambos lados. La pérdida aparece al sustituir el cuadrado de una suma sobre \(\psi\bmod m\) por la suma de cuadrados.

Los valores \(2.40,4.06,11.64,14.52\) son mediciones de casos pequeños. No existe una cota uniforme superior demostrada. En una prueba de \(D_{\rm NP}\le T(X)\), se puede usar \(D_{\rm NP}\le\mathcal M_W\), de modo que no se añade otro factor \(\lambda\) al mayorante de \(\mathcal M_W\). Tampoco se puede dividir por la \(\lambda\) medida para reclamar una mejora uniforme.

Por tanto, el balance correcto separa:
- inflación por levantamientos: \(A_r\);
- peso del módulo en la gran criba: \(r/\varphi(r)\);
- pérdida de Cauchy: \(\lambda\), ya incorporada al pasar al momento mayorante;
- tamaño del conjunto de módulos: \(|\mathcal D_S|=X^{1/2+o(1)}\).

No se deben multiplicar estos factores dos veces ni ocultar ninguno.

## 8. Análogo de D-NP para DIST-P

Para la pieza prima, el doc 115 demuestra

\[
\Omega_q^P=\frac1{\varphi(q)}=\frac1{q-1}
\]

y P-DOM:

\[
\Omega_q^C\ge\kappa'\Omega_q^P,
\qquad
\kappa'<3\log(19/10)
\]

para \(N\) suficientemente grande. Como el término \(E(N;q,N)\) es el mismo en ambas piezas y todos los sumandos son no negativos,

\[
\widehat D^P\le\frac1{\kappa'}\widehat D^C.
\]

Por ello,

\[
(\mathrm{DIST}\!-\!\mathcal D)\Longrightarrow
(\mathrm{DIST}\!-\!\mathcal D_P)
\]

con \(C_P=C_\Omega/\kappa'\). DIST-P no constituye un muro independiente una vez que P-DOM del doc 115 está incorporada.

El análogo de caracteres puede escribirse explícitamente. Sea

\[
L_q=\sum_{q^j\le X}\log q.
\]

Para \(q\) primo,

\[
(q-1)E(X;q,N)=
\sum_{\substack{\chi\bmod q\\\chi\ne\chi_0}}
\overline{\chi(N)}\Psi_X(\chi)-L_q.
\]

Así,

\[
\widehat D^P=
\sum_{q\in\mathcal H_P}
\frac1{q-1}
\left|
\sum_{\chi\ne\chi_0}\overline{\chi(N)}\Psi_X(\chi)-L_q
\right|^2.
\]

Todos los caracteres no principales módulo primo son primitivos; no aparece la multiplicidad de levantamientos de la rama compuesta. Por Cauchy,

\[
\widehat D^P
\le
2\sum_{q\in\mathcal H_P}
\sum_{\chi\ne\chi_0}|\Psi_X(\chi)|^2
+
2\sum_{q\in\mathcal H_P}\frac{L_q^2}{q-1}.
\]

El término de \(L_q\) es elemental. El primer término es un momento de caracteres sobre módulos primos. Esta representación podría estudiarse directamente, pero no hace falta introducir una hipótesis independiente para la cadena: P-DOM ya transfiere DIST-D a DIST-P. La ruta prima hereda la dificultad de DIST-D, no una nueva obstrucción de multiplicidad.

## 9. Estado de la cadena y decisiones

El estado documental es:

- **Demostrado:** unicidad de rama y dominación \(\widetilde\alpha_m\le\Omega_m^C\) (doc 121a).
- **Demostrado con enmienda:** reducción CRT de D-NP al momento \(W_d\) (doc 121b).
- **Demostrado:** P-DOM y transferencia DIST-D \(\Rightarrow\) DIST-P (doc 115).
- **Demostrado en el informe de conteo:** exponente asintótico \(1/2\) de \(|\mathcal D_S|\) y de \(\sum W_d\); medición superior \(\alpha=0.4659\).
- **Descartado:** la aplicación BZ que requería \(\alpha<0.420232\).
- **No suficiente:** gran criba clásica y Lema 8 de BZ con extracción del máximo del peso inducido.
- **Abierto:** una cota de segundo momento ponderada que use la distribución de \(A_r\) sin pagar su máximo; la parte lisa \(\mathcal M_\Omega^{liso}\); (H2); G2, la cota de Gram normalizada.
- **Fuera de esta cadena:** G0, separación de \(C_k\) y \(C_{k+1}\).

El doc 108 formulaba la cadena como \((\mathrm{DIST}\!-\!\mathcal D+\mathrm{DIST}\!-\!\mathcal D_P)+\mathrm{G2}\), junto con H2, para obtener UNIF-U-COND y \(C_k+C_{k+1}>0\). Con P-DOM del doc 115, DIST-P es corolario de DIST-D; por tanto, la hipótesis de distribución puede reducirse a DIST-D, manteniendo G2 y H2. Esto no separa las capas ni demuestra Goldbach.

**Decisión:** no reintentar Baier–Zhao con el mismo conteo ni con el mismo mayorante por máximo. Cualquier ruta futura debe empezar por una desigualdad ponderada que explique cómo evita el coste \(A_r\), o demostrar una cancelación/estructura adicional que no aparece en la gran criba estándar. Hasta entonces, D-NP permanece abierto.

**Goldbach fuerte sigue abierta.**

---

## Fuentes de la serie

- [Doc 108 — Auditoría de la cadena completa](108_auditoria_cadena_completa.md).
- [Doc 115 — P-DOM y cierre de G1 módulo DIST-D](115_revision_pdom.md).
- [Doc 120b — Auditoría adversarial de DIST-D](120b_auditoria_distd.md).
- [Doc 121 — Verificación de unicidad de rama y reducción CRT](121_verificacion_rama_dnp.md).
- [Doc 121b — Reducción CRT de D-NP](121b_dnp_crt.md).
- [Medición y teorema del conteo de \(\mathcal D_S\)](../verificacion-conteo-ds-claude-code.md).
- Baier–Zhao, [*Bombieri–Vinogradov type theorems for sparse sets of moduli*](https://www.impan.pl/shop/publication/transaction/download/product/82297?download.pdf=), §4, Lemas 7–9.
