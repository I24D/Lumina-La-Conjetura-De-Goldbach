# 99 — Primer diagnóstico numérico del segundo factor corregido

Fecha: 2026-11-27.
Procedencia: LUMINA, GPT-6.2; PC Lumina. Sigue a docs 94, 96 y 98 y a la revisión independiente de Dal (doc 97).
Estado: **[NUM]** cálculo del segundo factor en un soporte finito pequeño, con caracteres primitivos enumerados directamente y comparación contra la identidad de Möbius corregida. **[N]** cota analítica del segundo factor. No se demuestra Goldbach.

## 1. Objetivo

La prioridad ya no es producir otra identidad para el primer factor Gram. Se calcula en un proxy finito
\[
\mathcal B_{>R_0}^{\rm corr}(N)=
\sum_{\substack{a\mid A,\ h\in\mathcal H_A\\ah>R_0}}
\frac{|S_{ah}(N)|^2}{\varphi(ah)}
\]
y se separa, conductor por conductor,
\[
S_r=M_r+\Delta_r,\qquad
M_r=\sum_{d\mid r}\mu(r/d)\varphi(d)E(N;d,N),
\]
de modo que
\[
\mathcal B^{\rm corr}
=B_M+B_{\rm cross}+B_\Delta,
\]
con
\[
B_M=\sum_{r>R_0}\frac{|M_r|^2}{\varphi(r)},\quad
B_{\rm cross}=2\Re\sum_{r>R_0}\frac{M_r\overline{\Delta_r}}{\varphi(r)},\quad
B_\Delta=\sum_{r>R_0}\frac{|\Delta_r|^2}{\varphi(r)}.
\]
La comparación decisiva es con la suma directa sobre caracteres primitivos, no solo con la identidad algebraica.

## 2. Construcción del soporte

Script: numerico/99_second_factor_support_proxy.py.
Salida: numerico/resultados/20261127_99_second_factor_support_proxy.txt.

Se usó
\[
N=1000000\ (\text{par}),\quad A=3,\quad X=\lfloor N/A\rfloor=333333,
\]
\[
D=\lfloor\sqrt X\rfloor=577,\quad z=\lfloor X^{1/10}\rfloor=3,\quad
y=\lfloor X^{1/3}\rfloor=69.
\]
Es un proxy pequeño, con \(D=\sqrt X\) sin pérdida logarítmica. No reproduce el régimen asintótico de doc 94: en particular, \(z=3\) es demasiado pequeño para estudiar el límite \(z\to\infty\).

Se enumeraron directamente todos los pares \((s,d)\) del soporte fijo:
\[
z\le s<y,\ s\text{ primo},\ (s,AN)=1,\quad
d\in V^+(D/s),\ d\mid P(z),\ d<D/s.
\]
Resultado: 16 pares \((s,d)\), 17 partes nuevas \(h\), y 33 conductores completos distintos \(r=ah\), con \(a\mid3\), \((a,h)=1\), \(r>1\).

## 3. Verificación independiente de S_r

Para cada \(r\) del soporte:
1. Se construyeron las sumas de \(\Lambda(n)\) sobre \(n=p^j\le N\).
2. Se calculó directamente \(\psi(N;d,N)\) para cada divisor \(d\mid r\).
3. Como \(N\) es par, el soporte \(h\) excluye al primo 2; por tanto, los 33 conductores \(r=ah\) son impares. Se enumeraron explícitamente todos sus caracteres primitivos mediante raíces primitivas y logaritmos discretos.
4. Se calculó \(S_r\) directamente como suma de caracteres y, por separado, \(M_r+\Delta_r\) mediante la fórmula de doc 98.

Resultados:
- 33 conductores comprobados; los 33 tienen suma primitiva no nula.
- Error absoluto máximo entre la suma directa y \(M_r+\Delta_r\): \(1.374088\cdot10^{-8}\).
- Error relativo máximo, normalizado por \(\max(1,|S_r|)\): \(4.186868\cdot10^{-11}\).

Esto valida en el soporte de este proxy tanto la fórmula corregida como la enumeración de caracteres. Sigue siendo una prueba finita con aritmética de coma flotante.

## 4. Descomposición numérica del segundo factor

Se tomó \(R_0=(\log N)^\kappa\). Para cada \(\kappa\), la tabla incluye todos los conductores completos \(r>R_0\).

| \(\kappa\) | Conductores altos | \(B_M\) | \(B_{\rm cross}\) | \(B_\Delta\) | Cambio relativo total |
|---:|---:|---:|---:|---:|---:|
| 1.0 | 29 | \(1.9819922\cdot10^7\) | \(1.5136\cdot10^4\) | \(2.7844\cdot10^2\) | \(7.7770\cdot10^{-4}\) |
| 1.5 | 16 | \(1.0732769\cdot10^7\) | \(-1.3930\cdot10^4\) | \(9.8086\cdot10^1\) | \(-1.2888\cdot10^{-3}\) |
| 2.0 | 1 | \(1.2953852\cdot10^5\) | \(-3.5666\cdot10^1\) | \(2.4551\cdot10^{-3}\) | \(-2.7532\cdot10^{-4}\) |

Para \(\kappa=1\):
\[
\frac{B_{\rm cross}}{B_M}=7.636537\cdot10^{-4},\qquad
\frac{B_\Delta}{B_M}=1.404833\cdot10^{-5},
\]
\[
\frac{B_{\rm corr}-B_M}{B_M}=7.777020\cdot10^{-4}.
\]
La energía \(B_\Delta\) es solo \(1.40\cdot10^{-5}\) de \(B_M\), pero el término cruzado domina el efecto neto. La suma corregida coincide con la suma directa de \(|S_r|^2/\varphi(r)\); la diferencia absoluta de los dos totales fue \(9.999\cdot10^{-6}\).

**Comparación con doc 97:** no se reproduce exactamente el peso \(1.24\cdot10^{-4}\) comunicado por Dal. Aquí la energía aislada de \(\Delta\) pesa \(1.40\cdot10^{-5}\), mientras que el cambio total, incluido el término cruzado, es \(7.78\cdot10^{-4}\). El soporte y la ponderación no son los diez conductores de Dal; hay que alinear la definición exacta de “peso” antes de comparar.

La variación con \(\kappa\) incluye cambio de signo: \(-1.29\cdot10^{-3}\) para \(\kappa=1.5\). Para \(\kappa=2\) queda un solo conductor alto, así que ese dato no es estadísticamente informativo. El término cruzado puede reforzar o cancelar; no se debe presentar el \(1.24\cdot10^{-4}\) como constante universal.

## 5. Qué demuestra y qué no demuestra

**Sí demuestra numéricamente:**
- La corrección \(\Delta_r\) puede incorporarse conductor por conductor.
- La fórmula \(S_r=M_r+\Delta_r\) coincide con la suma explícita de caracteres en este soporte.
- El término cruzado \(M_r\overline{\Delta_r}\) debe conservarse; no basta medir la energía \(B_\Delta\).
- En este proxy, la corrección total al segundo momento es pequeña en términos relativos, del orden \(10^{-4}\) para \(\kappa=1\).

**No demuestra:**
- Que \(B_M\) sea pequeño respecto al objetivo de Chen.
- Que el peso de \(\Delta_r\) siga siendo pequeño uniformemente cuando \(N\to\infty\).
- Que el soporte real con \(z=X^{1/10}\) grande conserve esta distribución.
- Que la gran criba restringida a \(\mathcal H_A\) proporcione el ahorro necesario.

La cifra \(B_M\approx3.92\cdot10^7\) es el tamaño del momento en este proxy; el pequeño cambio relativo por \(\Delta\) no constituye por sí mismo una cota útil del momento.

## 6. Siguiente tarea

El siguiente objetivo analítico es una cota para \(B_M\) restringido a \(\mathcal R_A\), no otra identidad:
1. Escalar el cálculo a familias donde \(z\) crezca y conservar el factor \((\log X)^{C_0}\) en \(D_A\).
2. Medir la distribución de \(M_r\) por conductor, no solo la suma total; identificar si la masa se concentra en una subfamilia de \(\mathcal R_A\).
3. Buscar una estimación de segundo momento para los errores \(E(N;d,N)\) con la estructura de divisores \(d\mid r\) y la restricción \(r\in\mathcal R_A\), sin reemplazarla por el momento amplio que doc 92 muestra que no está disponible.
4. Mantener por separado \(B_M\), \(B_{\rm cross}\) y \(B_\Delta\) en cada escala.
5. Reproducir el conjunto exacto de diez conductores y la ponderación del doc 97 cuando se disponga de ellos.

Estado: segundo factor calculado en un proxy pequeño [NUM]; segundo factor asintóticamente acotado [N]. Goldbach fuerte sigue abierta.