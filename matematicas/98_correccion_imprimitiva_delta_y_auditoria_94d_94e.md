# 98 — Corrección imprimitiva Δ_r y auditoría reproducible de 94d/94e

Fecha: 2026-11-27.
Procedencia: revisión independiente de Dal (doc 97), lectura directa de docs 94 y 96 y prueba nueva en la PC LUMINA.
Estado: **[T]** fórmula exacta corregida para S_r; **[NUM]** 16 verificaciones directas de caracteres primitivos y progresiones; **[T]** reglas de atribución de masa de 94d/94e aclaradas; **[N]** segundo factor de Cauchy sigue sin cota útil. No se demuestra Goldbach.

## 1. Errata del doc 96

El doc 96, sección 5, afirmó incorrectamente que
\[
S_r(N)=\sum_{d\mid r}\mu(r/d)\varphi(d)E(N;d,N)
\]
era exacta. Falta la corrección por inducción de caracteres y por las potencias primas eliminadas al pasar de un carácter primitivo a su carácter inducido. La ecuación queda **reemplazada** por
\[
\boxed{\displaystyle
S_r(N)=\sum_{d\mid r}\mu(r/d)\varphi(d)E(N;d,N)+\Delta_r(N).
}
\]
No se debe usar la versión sin \(\Delta_r\) en ninguna desigualdad posterior.

## 2. Derivación exacta de Δ_r

Supóngase \(r>1\) y \((r,N)=1\). Definamos
\[
\Theta_p(N)=\sum_{p^j\le N}\log p,\qquad
\Theta(d;N)=\sum_{p\mid d}\Theta_p(N),
\]
y
\[
S_c(N)=\sum_{\chi^*\bmod c}^{*}
\overline{\chi^*(N)}\,\psi(N,\chi^*).
\]
Para \(c\mid d\), \(c>1\), sea
\[
I_{d,c}(N)=
\sum_{\chi^*\bmod c}^{*}\overline{\chi^*(N)}
\sum_{\substack{p\mid d\\p\nmid c}}
\sum_{p^j\le N}\chi^*(p)^j\log p.
\]

Por ortogonalidad de todos los caracteres módulo \(d\), excluyendo el principal,
\[
\varphi(d)E(N;d,N)
=\sum_{\substack{\chi\bmod d\\\chi\ne\chi_0}}
\overline{\chi(N)}\psi(N,\chi_d)-\Theta(d;N).
\]
Un carácter primitivo \(\chi^*\) de conductor \(c\mid d\) induce un carácter módulo \(d\) que satisface
\[
\psi(N,\chi_d)=\psi(N,\chi^*)-
\sum_{\substack{p\mid d\\p\nmid c}}
\sum_{p^j\le N}\chi^*(p)^j\log p.
\]
La diferencia son exactamente los términos \(p^j\) que el carácter inducido anula porque \(p\mid d\), aunque \(p\nmid c\). Al agrupar por conductores primitivos:
\[
\boxed{\displaystyle
\varphi(d)E(N;d,N)=
\sum_{\substack{c\mid d\\c>1}}S_c(N)
-\sum_{\substack{c\mid d\\c>1}}I_{d,c}(N)
-\Theta(d;N).
}
\tag{98.1}
\]

Multiplicando (98.1) por \(\mu(r/d)\) y sumando sobre \(d\mid r\), la inversión de Möbius da
\[
\sum_{d\mid r}\mu(r/d)\varphi(d)E(N;d,N)
=S_r(N)-\Delta_r(N),
\]
donde la forma general exacta de la corrección es
\[
\boxed{\displaystyle
\Delta_r(N)=
\sum_{d\mid r}\mu(r/d)
\left[
\sum_{\substack{c\mid d\\c>1}}I_{d,c}(N)
+\Theta(d;N)
\right].
}
\tag{98.2}
\]

### Simplificación de la corrección

La suma de Möbius sobre \(d\) puede evaluarse explícitamente.

Para el término de inducción, fijados \(c\mid r\) y \(p\nmid c\), el coeficiente de las potencias de \(p\) contiene
\[
\sum_{\substack{d\mid r\\c\mid d\\p\mid d}}\mu(r/d).
\]
Es cero salvo cuando \(r/c=p\), con \(p\nmid c\). Esto equivale a \(p\parallel r\), \(c=r/p\). Para el término principal,
\[
\sum_{\substack{d\mid r\\p\mid d}}\mu(r/d)
=\begin{cases}1,&r=p\text{ primo},\\0,&\text{en otro caso}.\end{cases}
\]
Por tanto, (98.2) se reduce a
\[
\boxed{\displaystyle
\Delta_r(N)=
\mathbf 1_{\{r\text{ primo}\}}\Theta_r(N)
+
\sum_{\substack{p\parallel r\\r/p>1}}
\ \sum_{\chi^*\bmod r/p}^{*}
\overline{\chi^*(N)}
\sum_{p^j\le N}\chi^*(p)^j\log p.
}
\tag{98.3}
\]
La condición \(p\parallel r\) es esencial: el exponente de \(p\) en \(r\) debe ser exactamente uno. Si \(p^2\mid r\), la suma de Möbius cancela esa contribución. Para \(r\) primo, la corrección es el término de potencias del carácter principal que aparece en \(E(N;r,N)\). Para \(r\) compuesto, la corrección recoge los caracteres de conductores \(r/p\) que se vuelven imprimitivos al inducirlos a módulos divisibles por \(p\).

La fórmula (98.2) es la forma más segura para derivaciones; (98.3) es la forma reducida y eficiente para cómputo.

## 3. Verificación numérica independiente

Script: numerico/97_imprimitive_correction_check.py.
Salida: numerico/resultados/20261127_97_imprimitive_correction_check.txt.
Método:
- \(N=999983\).
- Sumas reales de \(\Lambda(n)\) sobre todos los primos y potencias primas \(n=p^j\le N\).
- Sumas de progresiones \(\psi(N;d,N)\) calculadas directamente.
- Caracteres primitivos construidos explícitamente mediante raíces primitivas y logaritmos discretos para cada factor impar \(p^a\).
- Comparación directa de \(S_r\) con \(M_r+\Delta_r\), donde \(M_r=\sum_{d\mid r}\mu(r/d)\varphi(d)E(N;d,N)\).

Se probaron 16 conductores impares:
\[
3,5,7,11,15,21,35,55,105,385,9,25,45,75,225,315.
\]
Resultados:
- Error absoluto máximo: \(6.5627\cdot10^{-9}\).
- Error relativo máximo: \(9.7032\cdot10^{-12}\).
- \(|\Delta_r|/|S_r|\) en esta selección: de 0 hasta 2.5886%.

Esto reproduce de forma independiente la necesidad y el tamaño porcentual de la corrección. El resultado es una verificación de coma flotante sobre una familia finita, no la prueba algebraica: la prueba es la inversión en (98.1)–(98.3). La selección de 16 conductores no es necesariamente la misma selección de 10 del doc 97 de Dal; el peso \(1.24\cdot10^{-4}\) comunicado por Dal en el segundo momento no se declara aquí reproducido.

## 4. Inserción correcta en el segundo factor

Para cada conductor completo \(r=ah\) del soporte, debe usarse
\[
S_{a,h}(N)=
\sum_{d\mid ah}\mu(ah/d)\varphi(d)E(N;d,N)+\Delta_{ah}(N).
\]
En consecuencia,
\[
\mathcal B_{>R_0}(N)=
\sum_{\substack{a\mid A,\ h\in\mathcal H_A\\ah>R_0}}
\frac1{\varphi(ah)}
\left|
\sum_{d\mid ah}\mu(ah/d)\varphi(d)E(N;d,N)
+\Delta_{ah}(N)
\right|^2.
\]
No es correcto borrar \(\Delta\) dentro del cuadrado por ser pequeño en valor absoluto en algunos ejemplos: la cota debe controlar el término cruzado y \(\sum|\Delta_{ah}|^2/\varphi(ah)\) en el soporte completo. La medición de Dal de peso \(1.24\cdot10^{-4}\) es alentadora, pero sigue siendo una observación finita hasta reproducir exactamente su selección, pesos y rango.

## 5. Auditoría de las atribuciones de masa en 94d y 94e

### 94d: índice \(r\) de la diagonalización, no conductor \(ah\)

En numerico/94d_diagonal_profile.py:
- \(X=\lfloor N/A\rfloor\).
- El proxy fija \(D=\lfloor\sqrt X\rfloor\), sin la pérdida \((\log X)^{C_0}\) del parámetro asintótico \(D_A\).
- \(z=\lfloor X^{1/10}\rfloor\).
- \(L_s=\lfloor D/s\rfloor\), que aproxima \(D_A/s\).
- Para cada divisor \(r\mid P(z)\), se calcula \(C_s(r)=\sum_{d\in V^+(L_s),\,r\mid d}\mu(d)/\varphi(d)\) y \(q_r=\varphi(r)C_s(r)^2\).

La partición que implementa el script es exactamente
\[
Q_s^{\rm low}=\sum_{r\le L_s/z}q_r,\qquad
Q_s^{\rm high}=\sum_{r>L_s/z}q_r,\qquad
Q_s=Q_s^{\rm low}+Q_s^{\rm high}.
\]
Es una atribución por el índice divisor \(r\) de la identidad positiva; no es una partición por \(d\), por \(e\), ni por el conductor completo \(ah\). La frontera es \(L_s/z\), usando el \(L_s\) entero del proxy.

### 94e: conductor completo \(ah\)

En numerico/94e_conductor_cutoff_profile.py, el diccionario \(K[h]\) procede de la suma original de coeficientes del doc 94. Para cada \(a\mid A\), el script asigna exactamente la masa
\[
m(a,h)=\varphi(a)\varphi(h)|K_A(h)|^2
=\varphi(ah)|K_A(h)|^2,
\]
porque \((a,h)=1\). La masa total es
\[
\sum_{a\mid A}\sum_{h\in\mathcal H_A}m(a,h)
=\left(\sum_{a\mid A}\varphi(a)\right)
\sum_h\varphi(h)|K_A(h)|^2
=A\mathcal G_A.
\]
La masa baja es la suma de esos mismos términos con \(ah\le R_0\), \(R_0=(\log N)^\kappa\); la masa alta es el complemento. Por tanto, la regla de atribución es exacta para el proxy \(N=10^{12},A=3\), con el \(D=\lfloor\sqrt X\rfloor\) sin pérdida logarítmica.

En \(\kappa=2.7\), el script mide 92.878861% bajo el corte y 7.121139% arriba. La razón entre masa total y masa alta es 14.04; la mejora del factor de Cauchy para \(|R|\), que lleva raíz cuadrada, es \(\sqrt{14.04}=3.75\). En \(\kappa=2.8\), esas cifras son 95.249831%, razón 21.05 y mejora de amplitud 4.59. Así se aclara la frase del doc 96: 4–5 veces se refería a la cota de amplitud tras Cauchy; una cifra cercana a 16 veces puede referirse al factor de energía/cuadrado, no a \(|R|\). Los modelos alternativos de Dal (93.3–100%) no quedan reconciliados todavía: deben igualarse \(N,A,D_A,\kappa\), soporte y regla de masa antes de comparar.

## 6. Decisión de investigación: el segundo factor es ahora prioritario

No se necesita otra identidad para el primer factor. La siguiente tarea es estimar
\[
\mathcal B_{>R_0}^{\rm corr}(N)=
\sum_{\substack{a\mid A,\ h\in\mathcal H_A\\ah>R_0}}
\frac{|S_{a,h}(N)|^2}{\varphi(ah)}
\]
con \(S_{a,h}\) corregido por \(\Delta_{ah}\).

Secuencia concreta:
1. Fijar la familia de parámetros de doc 94 y el \(\kappa\) que se auditará; no mezclar el proxy \(D=\sqrt X\) con \(D_A=\sqrt X/(\log X)^{C_0}\).
2. Construir una tabla exacta del soporte de conductores completos \(r=ah\), conservando multiplicidad y \(\varphi(r)\).
3. Separar \(\mathcal B\) en la parte de Möbius \(M_r\), el término cruzado \(2\Re(M_r\overline{\Delta_r})\) y la energía \(|\Delta_r|^2\). Medir cada uno por separado.
4. Estimar el momento de \(M_r\) restringido a \(\mathcal R_A\), sin sustituirlo por el momento amplio en todos los módulos ni invocar una cota que equivalga a GRH.
5. Reproducir el peso \(1.24\cdot10^{-4}\) del doc 97 con sus diez conductores y su regla de ponderación, cuando el conjunto exacto esté disponible.

Estado: diagonalización Gram [T]; producto completo [T]; corrección \(\Delta_r\) [T]; comprobación directa finita [NUM]; atribución de masa 94d/94e [NUM y definición exacta]; cota útil del segundo factor corregido [N]. Goldbach fuerte sigue abierta.
