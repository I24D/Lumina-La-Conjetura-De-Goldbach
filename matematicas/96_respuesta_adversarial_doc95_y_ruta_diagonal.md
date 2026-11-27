# 96 — Respuesta LUMINA a la revisión adversarial del doc 94

Fecha: 2026-11-27.
Procedencia: revisión de Dal comunicada en conversación (resumen del doc 95), lectura completa del doc 94 en OneDrive y nuevos controles 94d/94e ejecutados en la PC LUMINA.
Alcance: no se tuvo el archivo íntegro del doc 95 ni sus seis ítems [P-rev]; se evalúan las cinco conclusiones reproducidas por Dal. No se afirma haber revalidado los seis ítems originales.

## 1. Veredicto

1. **Identidad Gram/gcd [T]: conservar.** La identidad del doc 94 sigue siendo exacta. La revisión informa seis escenarios independientes con error relativo <=1.8e-14. No se ha detectado una pérdida de signo en esa identidad.
2. **Límite fijo [T]: conservar. Límite variable: retirar o reparar.** El doc 94 usa el límite fijo z_A. La versión variable no debe mezclarse con él: al permitir dos primos de switching en h, la fórmula de un solo switching omite conductores. Si se conserva la versión variable, se indexa por el menor primo de switching s y se agrega la cola completa s'>s. La corrección numérica indicada por Dal es consistente con esta reparación; falta contrastar sus seis casos con el doc 95 íntegro.
3. **Cota Q_s [N]: no demostrada.** Los datos son compatibles con Q_s de orden 1/log z en la familia N=2^e, A=3. La constante observada, aproximadamente 2.5 en Q_s log z, no proporciona por sí sola margen suficiente para el balance posterior ni prueba uniformidad.
4. **Puente a R_A [N]: el cuello es el segundo factor.** El doc 94 sí contiene la desigualdad formal de Cauchy en la sección 6. Lo que no contiene es una cota demostrada y suficientemente fuerte para el momento de bloques B_A(N). Por tanto, la objeción es correcta en sustancia: una cota de G_A, aun demostrada, no basta.
5. **Corte ah<=R_0 [T/N]: separar exactamente ambos bloques.** El tratamiento de conductores completos pequeños por Siegel–Walfisz y el bloque grande deben aparecer como dos piezas complementarias. La masa Gram pequeña no es automáticamente la masa del resto R_A: es energía de coeficientes, no una estimación de los bloques de caracteres.

## 2. Diagonalización positiva exacta

Sea L_s=D_A/s y
\[
C_s(r)=\sum_{\substack{d\in V^+(L_s)\\r\mid d}}
\frac{\mu(d)}{\varphi(d)},\qquad r\mid P(z_A).
\]
La identidad gcd da
\[
(d,e)=\sum_{r\mid(d,e)}\varphi(r),
\]
y, por tanto, la representación positiva exacta
\[
\boxed{\displaystyle
Q_s=\sum_{r\mid P(z_A)}\varphi(r)\,C_s(r)^2.}
\]
Esta es la vía diagonalizada que debe estudiarse: convierte el núcleo gcd en una suma de cuadrados no negativos.

**Advertencia:** no es válido quedarse solamente con los términos d=e de la suma doble original. Esa diagonal es
\[
D_{\rm diag}=\sum_{d\in V^+(L_s)}\frac1{\varphi(d)}
\]
y puede ser varias veces mayor que Q_s; la cancelación entre d distintos sigue siendo esencial.

Como referencia exacta, si se reemplaza V^+(L_s) por el conjunto de todos los divisores cuadrados-libres de P(z_A), entonces
\[
C_{\rm full}(r)=\frac{\mu(r)}{\varphi(r)}
\prod_{\substack{p<z_A\\p\nmid ANr}}
\left(1-\frac1{p-1}\right),
\]
y
\[
Q_{\rm full}(z_A)=
\prod_{\substack{p<z_A\\p\nmid AN}}
\left(1-\frac2{p-1}+\frac{p}{(p-1)^2}\right).
\]
La identidad del producto fue verificada numéricamente en 94d con error máximo 9.44e-16.

El lema preciso que queda abierto es comparar la restricción Rosser–Iwaniec con esta referencia:
\[
\boxed{\displaystyle
Q_s\ll
P_{A,N}(z_A):=
\prod_{\substack{p<z_A\\p\nmid AN}}
\left(1-\frac2{p-1}+\frac{p}{(p-1)^2}\right),
\quad
u_s=\frac{\log(D_A/s)}{\log z_A}\in[5/3-o(1),4-o(1)].
}
\]
La forma normalizada P_{A,N}(z_A) es preferible a escribir sin más 1/log z_A. Para N variable, los factores locales excluidos por p|N introducen el factor singular correspondiente. En efecto, para p>=3,
\[
f(p)^{-1}=
\frac{p-1}{p-2}
\left(1-\frac1{p^2-3p+3}\right),
\quad
f(p)=1-\frac2{p-1}+\frac{p}{(p-1)^2}.
\]
El producto de los segundos factores converge a una constante positiva; el producto de (p-1)/(p-2) sobre p|N es comparable al factor local de la serie singular de Goldbach. En la sucesión de N primoriales, con A coprimo elegido fuera de sus factores, el producto de referencia puede crecer como (log log N)/log z_A. Esto no demuestra que Q_s tenga ese crecimiento, pero significa que Q_s<<1/log z_A no se deduce de la referencia completa y requeriría una cancelación adicional de V+. Los controles N=2^e, A=3 mantienen fijo ese factor local y no prueban uniformidad para todos los N pares.

## 3. Controles nuevos 94d: perfil de la diagonalización

Script: numerico/94d_diagonal_profile.py.
Familia: N=2^e, e=40,45,50,55,60; A=3; cuatro valores objetivo de u_s entre 5/3 y 4; 20 escenarios.

Resultados:
- Q_s log z entre 1.55603 y 2.82873.
- Q_s/P_{A,N}(z) entre 1.00266 y 1.70054.
- La identidad Q_full=P_{A,N}(z) por el producto local tiene error máximo 9.44e-16.
- D_diag/Q_s entre 1.00 y 5.82: la diagonal d=e sola no prueba la cota.
- En la representación positiva, la porción r>L_s/z_A llega hasta 9.84% de Q_s; la porción r>sqrt(L_s), hasta 18.15%. No se puede borrar la cola por una afirmación de pequeñez sin demostrarla.
- En los r presentes, el máximo observado de |C_s(r)/C_full(r)| fue 2.66. Es evidencia finita, no una cota uniforme.

Subproblema recomendado: separar r<=L_s/z_A y r>L_s/z_A. En el primer rango, L_s/r>=z_A y podría ser posible reutilizar el límite fijo de criba, tras tratar la condición V^+ al factorizar d=r e. En el segundo, e<L_s/r<z_A; se necesita una cota de cola específica. La condición V^+ depende del orden de los factores primos de r e, por lo que no se puede aplicar el teorema de criba a e sin justificar esa transformación.

## 4. Corte por conductor completo ah<=R_0

Defínanse
\[
\mathcal G_{\le R_0}=
\sum_{\substack{a\mid A,\ h\in\mathcal H_A\\ah\le R_0}}
\varphi(ah)|K_A(h)|^2,
\quad
\mathcal G_{>R_0}=
\sum_{\substack{a\mid A,\ h\in\mathcal H_A\\ah>R_0}}
\varphi(ah)|K_A(h)|^2.
\]
Como (a,h)=1 y \(\sum_{a\mid A}\varphi(a)=A\),
\[
\mathcal G_{\le R_0}+\mathcal G_{>R_0}=A\mathcal G_A.
\]
El bloque pequeño se controla por separado con Siegel–Walfisz, como en la sección 3 del doc 94. Para el bloque grande, Cauchy debe usar la energía restringida:
\[
|R_{>R_0}^{\rm prim}|^2
\le \mathcal G_{>R_0}\,\mathcal B_{>R_0}(N),
\]
\[
\mathcal B_{>R_0}(N)=
\sum_{\substack{a\mid A,\ h\in\mathcal H_A\\ah>R_0}}
\frac{|S_{a,h}(N)|^2}{\varphi(ah)}.
\]
No reemplazar \(\mathcal G_{>R_0}\) por A G_A si se busca aprovechar el corte.

### Perfil numérico del corte 94e

En el proxy de doc 94, N=10^12 y A=3, se calculó la fracción de la energía total A G_A con ah<=R_0 y R_0=(log N)^kappa:
- kappa=1: 49.17%.
- kappa=2: 73.28%.
- kappa=2.7: 92.88%.
- kappa=2.8: 95.25%.
- kappa=3: 99.11%.

Esto hace plausible el 94% reportado, pero muestra que depende sensiblemente de kappa. El documento debe fijar kappa y la familia de parámetros antes de usar ese porcentaje. Además, esa fracción es masa de coeficientes Gram; no mide directamente la contribución de R_A.

## 5. Puente exacto hacia los bloques de caracteres

Para un conductor r>1 coprimo con N, sea
\[
S_r(N)=\sum_{\chi^*\bmod r}^{*}
\overline{\chi^*(N)}\,\psi(N,\chi^*).
\]
La ortogonalidad de caracteres primitivos, por inversión de Möbius sobre los conductores, da la identidad exacta
\[
\boxed{\displaystyle
S_r(N)=
\sum_{d\mid r}\mu(r/d)\varphi(d)\,E(N;d,N).
}
\]
La suma principal se cancela porque \(\sum_{d\mid r}\mu(r/d)=0\) para r>1. Aquí d|r implica (d,N)=1, como requiere E(N;d,N).

Por tanto, el momento grande puede escribirse sobre el conjunto disperso de conductores
\[
\mathcal R_A=\{r=ah:\ a\mid A,\ h\in\mathcal H_A,\ (a,h)=1\}
\]
como
\[
\mathcal B_{>R_0}(N)=
\sum_{\substack{r\in\mathcal R_A\\r>R_0}}
\frac1{\varphi(r)}
\left|
\sum_{d\mid r}\mu(r/d)\varphi(d)E(N;d,N)
\right|^2.
\]
Esta transformación ofrece el puente para una gran criba restringida a H_A. No la demuestra: una aplicación de la gran criba general con Cauchy puede reintroducir el momento amplio de a=N que el doc 92 identifica como insuficiente o GRH-implicante. El ahorro debe venir de la estructura de R_A y de la cancelación de Möbius entre divisores.

## 6. Siguiente secuencia de investigación

1. Mantener el límite fijo z_A como línea principal; no reabrir la versión variable salvo con la cola completa s'>s.
2. Probar el lema diagonalizado Q_s<=C P_{A,N}(z_A), con C uniforme en u_s y en los parámetros permitidos. Atacar primero la suma C_s(r) para r<=L_s/z_A y después la cola r>L_s/z_A.
3. Fijar kappa y medir \(\mathcal G_{\le R_0}\), \(\mathcal G_{>R_0}\) por separado en una familia más amplia de N y A.
4. Desarrollar la gran criba restringida usando la identidad de Möbius para S_r. El objetivo sigue siendo una cota útil para \(\mathcal B_{>R_0}(N)\), no solo una cota absoluta.
5. Combinar ambas piezas y volver a auditar las constantes contra el balance de Chen.

Estado final: identidad Gram [T]; diagonalización gcd [T]; comparación uniforme V^+ con el producto completo [N]; momento restringido de bloques [N]; puente cuantitativo suficiente hacia R_A [N]. Goldbach fuerte sigue abierta.
