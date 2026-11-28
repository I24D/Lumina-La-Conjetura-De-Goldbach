# 94 — UNIF-U-COND: compresión por conductor y norma de Gram

**Procedencia:** LUMINA, GPT-6.2; 2026-11-27. Método: derivación algebraica, auditoría adversarial de índices y cómputos reproducibles en Python estándar. Fuentes de contexto: docs 81–93; Nathanson, *Additive Number Theory: The Classical Bases*, caps. 9–10; identidad \(\sum_{d\mid n}\varphi(d)=n\).

Fecha: 2026-11-27.
Estado: **[T]** descomposición exacta por conductor admisible, soporte de K_A, identidad Gram/gcd y bloque de conductores completos pequeños; **[NUM]** verificación finita reproducible en PC; **[N]** cota asintótica de Gram y momento de bloques primitivos. No se demuestra UNIF-U ni Goldbach.

## 0. Objetivo y correcciones de alcance

Se estudia el resto de switching para A=Qu:
\[
R_A=\sum_{m\in\mathcal M_A}b_A(m)E(N;Am,N),\qquad
E(N;q,N)=\psi(N;q,N)-\frac{\psi(N)}{\varphi(q)}.
\]
En la aplicación, A=Qu con rad(A)|Q_N y Q_N tiene k-1 factores primos distintos; por tanto A/\varphi(A)=O_k(1). Póngase X_A=N/A, z_A=X_A^{1/10}, y_A=X_A^{1/3}, y D_A=X_A^{1/2}/(\log X_A)^{C_0}, con C_0 suficientemente grande para el ahorro logarítmico de BV. Se usa el límite de criba fijo z_A del doc. 87, y P(z_A)=\prod_{p<z_A,\ p\nmid AN}p. El soporte es
\[
m=sd,\quad z_A\le s<y_A,\quad s\ {\rm primo},\quad
d\in V^+(D_A/s),\quad d\mid P(z_A),\quad (sd,AN)=1,
\]
y b_A(sd)=lambda_d^+(s)=mu(d) en el soporte. Para beta=2, si d=p_1\cdots p_j con p_1>\cdots>p_j, la condición V^+ es
\[
p_\ell^2(p_1\cdots p_\ell)<D_A/s
\quad\text{para todo }\ell\text{ impar};
\]
además d es cuadrado-libre y todos sus primos son menores que z_A. La descomposición m=sd es única: s es el único factor primo de m que puede ser ≥z_A. Para cada s hay a lo sumo D_A/s enteros d<D_A/s, luego
\[
\#\mathcal M_A\le \sum_{\substack{z_A\le s<y_A\\s\ {\rm primo}}}\frac{D_A}{s}
\ll D_A,
\]
pues la suma de recíprocos de primos en [z_A,y_A) está acotada (y tiende a log(10/3)).

Dos correcciones respecto del primer intento:
1. r|Am si y solo si h|m, con h=r/(r,A), solo vale para conductores **admisibles**: a=(r,A)|A, h cuadrado-libre y (h,A)=1. Si r tiene una potencia de un primo de A mayor que la permitida por A, o una potencia cuadrada fuera de A, no aparece en ningún módulo Am.
2. Siegel–Walfisz controla conductores completos r≤(log N)^kappa. No controla automáticamente todos los conductores con parte nueva h≤(log N)^kappa, porque r=ah puede tener a|A grande.

## 1. Descomposición exacta por conductor primitivo [T]

Para un módulo q=Am, la ortogonalidad de caracteres da
\[
E(N;q,N)=\frac1{\varphi(q)}
 \sum_{\substack{\chi\bmod q\\\chi\ne\chi_0}}
 \overline{\chi(N)}\,\psi(N,\chi_q)
 +P_q,
\]
donde
\[
P_q=-\frac1{\varphi(q)}
 \sum_{p\mid q}\ \sum_{p^j\le N}\log p
\]
es la corrección del carácter principal.

Sea chi* primitivo de conductor r, inducido a q=Am. Si la contribución existe, escribimos
\[
a=(r,A),\qquad h=r/a.
\]
Entonces
\[
a\mid A,\quad (h,AN)=1,\quad \mu^2(h)=1,\quad h\mid m.
\]
Recíprocamente, estas condiciones implican r=ah|Am. Como (A,m)=1,
\[
\varphi(Am)=\varphi(A)\varphi(m).
\]

Por tanto el coeficiente de cualquier carácter primitivo de conductor r=ah depende solo de su parte nueva h:
\[
\boxed{
K_A(h)=\frac1{\varphi(A)}
 \sum_{\substack{m\in\mathcal M_A\\h\mid m}}
 \frac{b_A(m)}{\varphi(m)}.
}
\]
La parte primitiva del resto se agrupa exactamente como
\[
R_A^{\rm prim}=
 \sum_{\substack{a\mid A,\ h\in\mathcal H_A\\ah>1}}
 K_A(h)
 \sum_{\chi^*\bmod ah}^{*}
 \overline{\chi^*(N)}\,\psi(N,\chi^*),
\]
más las correcciones por inducción y por el carácter principal.

Aquí H_A={h: existe m en M_A con h|m}. El coeficiente es idéntico para todos los conductores ah que comparten la misma parte nueva h; no significa que K_A(h) sea cero.

### Correcciones por inducción [T]

Para el carácter inducido módulo q,
\[
\psi(N,\chi_q)=\psi(N,\chi^*)-
 \sum_{\substack{p\mid q\\p\nmid r}}
 \sum_{p^j\le N}\chi^*(p)^j\log p.
\]
Cada suma de potencias para un primo es O(log N), y omega(q)=O(log N). Tras promediar sobre los phi(q) caracteres de cada módulo, la corrección es O(log^2 N) por m. El término principal P_q tiene la misma cota. Como hay a lo sumo D_A módulos en el soporte,
\[
R_{\rm inducción}+R_{\rm principal}
 =O(D_A\log^2N).
\]
Con D_A menor que constante por (N/A)^(1/2)(log N)^(-C_0), esto es
\[
o\!\left(\frac{N}{A\log^2N}\right)
\]
para A≤N^(1-epsilon), tomando C_0 fijo suficientemente grande. No es el obstáculo de UNIF-U.

## 2. Soporte exacto de los coeficientes [T]

Como todo m en M_A tiene la forma m=sd, con s en [z_A,y_A) primo y todos los factores de d menores que z_A, se tiene:
\[
K_A(h)=0
\]
salvo que:
- h sea cuadrado-libre y (h,AN)=1;
- ningún primo de h sea ≥y_A;
- h tenga a lo sumo un factor primo en [z_A,y_A).

Así, la parte nueva del conductor no recorre todos los enteros hasta D_A: está restringida a números z_A-suaves o a un primo de switching multiplicado por una parte z_A-suave. Es una restricción exacta del soporte, no una estimación asintótica.

## 3. Conductores completos pequeños: sí se controlan [T]

Para r≤R_0=(log N)^kappa, y solo para estos conductores completos, la estimación elemental
\[
|K_A(h)|\le
 \frac1{\varphi(A)}
 \sum_{\substack{m\le D_A\\h\mid m}}\frac1{\varphi(m)}
 \ll \frac{\log(2D_A)}{\varphi(A)\varphi(h)}
 \ll_k\frac{\log N}{A}
\]
se sigue de m=hj, (h,j)=1, y sum_{j≤x}1/phi(j) menor que constante por log(2x).

Siegel–Walfisz, para cualquier B>0, da
\[
|\psi(N,\chi^*)|\ll_{\kappa,B}
 \frac{N}{(\log N)^B}\qquad(r\le R_0).
\]
Usando sum_{r≤R_0}phi(r) menor que constante por R_0^2,
\[
\boxed{
R_{r\le R_0}
 \ll_{k,\kappa,B}
 \frac{N}{A}(\log N)^{1+2\kappa-B}.
}
\]
Con B>2kappa+3, es o(N/(A log^2 N)).

**Límite de esta conclusión:** no elimina todos los h pequeños. Si h es pequeño pero a|A es grande, el conductor completo r=ah puede exceder el rango de Siegel–Walfisz. Esa parte queda en el bloque grande.

## 4. Nueva identidad exacta: norma de Gram como suma gcd [T]

Definamos la energía de los coeficientes:
\[
\mathcal G_A=\sum_{h\in\mathcal H_A}
 \varphi(h)|K_A(h)|^2.
\]
Al expandir el cuadrado e intercambiar sumas,
\[
\begin{aligned}
\mathcal G_A
&=\frac1{\varphi(A)^2}
 \sum_{m,n\in\mathcal M_A}
 \frac{b_A(m)b_A(n)}{\varphi(m)\varphi(n)}
 \sum_{h\mid(m,n)}\varphi(h)\\
&=\boxed{
 \frac1{\varphi(A)^2}
 \sum_{m,n\in\mathcal M_A}
 \frac{b_A(m)b_A(n)(m,n)}
      {\varphi(m)\varphi(n)}.
}
\end{aligned}
\]
Se usó la identidad exacta sum_{h|g}phi(h)=g.

Además, para m=sd, n=s'd', como los factores de d,d' son menores que z_A y s,s'≥z_A,
\[
(m,n)=
\begin{cases}
(d,d'),&s\ne s',\\
s(d,d'),&s=s'.
\end{cases}
\]
Por tanto la energía es una forma cuadrática explícita en los pesos mu(d), con núcleo gcd. Esta identidad conserva los signos que se pierden al aplicar valor absoluto término a término.

**Nuevo subproblema [N].** Probar uniformemente una cota de la forma
\[
\mathcal G_A\ll_k
 \frac1{\varphi(A)^2\log z_A}.
\]
La identidad no demuestra esa cota. La prueba tendría que explotar la estructura de los pesos superiores de Rosser–Iwaniec en la forma cuadrática gcd, no solo |lambda_d^+(s)|≤1.

## 5. Verificación finita reproducible [NUM]

Script: numerico/94_conductores_grama.py. Ejecutado en la PC LUMINA, salida:

    PARAMS N=1000000000000 A=3 X=333333333333 z=14 y=6933 D=577350 small_primes=[7, 11, 13]
    SUPPORT rows=1444 distinct_h=1451 bad_support=0
    K_A(r) examples:
     K(3)=0.427962281933
     K(7)=-0.067848416650
     K(11)=-0.028330369695
     K(13)=-0.020169817129
     K(21)=-0.067848416650
     K(17)=0.021527777778
    GRAM norm=0.450566952208
    direct_gram=0.450566952208
    abs_majorant=1.172357731746
    signed/abs=0.384325
    PASS

La igualdad GRAM norm = direct_gram verifica la identidad por dos cálculos independientes. El soporte pasó en los 1451 valores de h. La norma con signos es el 38.4% de la mayorante obtenida sustituyendo cada b(m) por |b(m)| en este ejemplo.

**Lectura honesta:** esto verifica una identidad y detecta cancelación finita en los coeficientes. No demuestra que el cociente 0.384 se mantenga, ni que G_A tenga el orden 1/log z_A.

### Verificación separada de la partición Gram [NUM]

Script adicional: numerico/94c_gram_split_check.py. Reconstruye K_A(h) desde la suma original sobre m y, por una segunda ruta, reconstruye las dos clases de conductores h=h_0 y h=s h_0 a partir de T_s(h_0). Resultado en el mismo proxy:
\[
G_{\rm smooth}=0.228419312130,\qquad
G_{\rm switch}=0.222147640078,
\]
\[
G_{\rm total}=0.450566952208,\qquad
\max|\text{diferencia de las tres identidades}|=3.1\cdot10^{-16}.
\]
Los 1444 pares (s,d), los 1451 conductores h y el soporte (sin factores >=y ni dos primos de switching) pasan. Esto verifica la partición finita, no la cota uniforme asintótica de Q_s.

Las ejecuciones en la PC LUMINA son controles de reproducibilidad, no una revisión independiente: el contrato del repositorio excluye como independientes a agentes que comparten el mismo usuario de Windows. La revisión independiente de Dal queda pendiente.

## 6. Corrección de Cauchy y bifurcación condicional

La cantidad de módulos del soporte es #M_A\ll D_A, no AD_A. Por ello, si se tuviera el momento **solo sobre el soporte**
\[
\sum_{m\in\mathcal M_A}|E(N;Am,N)|^2
 \ll N(\log N)^C,
\]
Cauchy daría
\[
|R_A|\ll \sqrt{\#\mathcal M_A\,N}\,(\log N)^{C/2}
 \ll \sqrt{D_A N}\,(\log N)^{C/2}
 \ll N^{3/4}A^{-1/4}(\log N)^{C/2-C_0/2}.
\]
Frente al objetivo N/(A log^2 N), el exponente de potencia del cociente es
\[
\boxed{\frac{3c-1}{4}\quad(A=N^c),}
\]
por lo que el umbral es c<1/3. Esto cubre c_*=0.039884 y 2/17. Si se reemplaza el soporte por todos los módulos q≤AD_A, el conteo vuelve a ser AD_A y se recupera el umbral c<1/5 del doc. 92. No confundir ambos enunciados.

El momento amplio sobre todos los módulos coprimos con N no se persigue. La forma corregida UNIF-U-2M-R es
\[
\sum_{\substack{q\le N^{(1+c)/2}\\(q,N)=1}}
 |E(N;q,N)|^2\ll N(\log N)^C
\]
para N par suficientemente grande. Implica GRH para todo carácter primitivo no principal de conductor impar: fije r impar y a reducido módulo r; para cada x elija un N par con N\equiv a\pmod r y |N-x|\le 2r. Para N grande, r\le N^{(1+c)/2} y (r,N)=1, de modo que el término q=r da
\[
|E(N;r,a)|\ll_r N^{1/2}(\log N)^{C/2}.
\]
Cambiar N por x altera las sumas de Chebyshev solo en O_r(log x), así que la misma cota vale para E(x;r,a). Por ortogonalidad,
\[
|\psi(x,\chi)|=\left|\sum_{a\bmod r}^{*}\chi(a)E(x;r,a)\right|
 \ll_r x^{1/2}(\log x)^{C/2}.
\]
La integral de Mellin de -L'/L converge entonces en Re(s)>1/2, excluye ceros en ese semiplano y, por la ecuación funcional, da GRH. La conclusión solo cubre conductores impares porque N se restringe a pares; no se puede escoger N par y coprimo con un r par.

La condición (q,N)=1 es esencial, no cosmética. Para N=10^6 y q=2, la clase N mod 2 no es reducida: el cálculo exacto da E=-999573.427699..., y |E|^2/[N(log N)^5]=1.98516. Sin excluir esos módulos, (92.1) falla ya en q=2.

Precisión lógica: UNIF-U-2M-R implica GRH y exige una estimación de segundo momento mucho más fuerte en su formulación. No hemos demostrado que GRH, por sí sola, no implique esa estimación; por eso no afirmamos una separación lógica estricta. Esta hipótesis amplia tampoco es el momento restringido al soporte M_A que se propone en este documento.

### Candidato más estrecho para UNIF-U-COND [N]

Sea
\[
S_{a,h}(N)=
 \sum_{\chi^*\bmod ah}^{*}
 \overline{\chi^*(N)}\,\psi(N,\chi^*).
\]
El bloque de conductores grandes es
\[
R_{>R_0}^{\rm prim}
 =\sum_{\substack{a\mid A,\ h\in\mathcal H_A\\ah>R_0}}
 K_A(h)S_{a,h}(N).
\]
Cauchy por bloques de conductor da
\[
|R_{>R_0}^{\rm prim}|^2
\le
\underbrace{\left(\sum_{\substack{a\mid A,\ h\in\mathcal H_A\\ah>R_0}}
 \varphi(ah)|K_A(h)|^2\right)}_{\le A\mathcal G_A}
\underbrace{\left(\sum_{\substack{a\mid A,\ h\in\mathcal H_A\\ah>R_0}}
 \frac{|S_{a,h}(N)|^2}{\varphi(ah)}\right)}_{\mathcal B_A(N)}.
\]
Si se demostraran conjuntamente
\[
\mathcal G_A\ll_k\frac1{\varphi(A)^2\log N},
\qquad
\mathcal B_A(N)\ll_k\frac{N^2}{A(\log N)^3},
\]
se obtendría R_{>R_0}^{prim} menor que constante por N/(A log^2 N): en efecto, la primera suma de Cauchy es a lo sumo A G_A, y A/phi(A)=O_k(1) porque rad(A)|Q_N con k-1 factores primos.

Este es un **condicional nuevo y específico**: B_A agrupa primero todos los caracteres primitivos de un mismo conductor y conserva sus fases conj(N); no es el segundo momento individual sum_chi |psi(N,chi)|^2 ni el momento amplio de errores de progresiones. No se ha probado que estas cotas sean ciertas ni que el condicional sea estrictamente más débil que GRH; eso debe auditarse, no suponerse.

## 7. Veredicto y siguiente experimento

- **[T]** Compresión exacta por parte nueva h, con la condición de admisibilidad corregida.
- **[T]** Conductores completos r≤(log N)^kappa controlados por Siegel–Walfisz.
- **[T]** Identidad Gram/gcd exacta; soporte de h restringido.
- **[NUM]** Identidad verificada en PC: 1444 pesos, 1451 partes nuevas, error de identidad inferior a 10^-10.
- **[N]** Cota uniforme de Gram y momento de bloques primitivos.
- **[N]** No hay cancelación exacta universal: K_A(1) no es cero en el ejemplo; el avance es una estructura de norma, no una anulación.

Siguiente paso: estudiar la forma cuadrática gcd de los pesos V^+(D_A/s), empezando por el caso de un único s, y comparar su norma con la cota 1/log z_A predicha por la criba lineal. En paralelo, buscar una reducción de B_A a una correlación de pares de caracteres/ceros restringida a H_A, sin reintroducir el momento amplio que fuerza GRH.

Goldbach sigue abierta.


## 8. Suplemento: reducción a bloques de un único primo de switching [T]

Hay dos clases de partes nuevas h: las que solo tienen factores menores que z_A y las que contienen exactamente un primo s de switching. Para h_0|P(z_A), definamos
\[
T_s(h_0)=\sum_{\substack{d\in V^+(D_A/s)\\h_0\mid d}}
 \frac{\mu(d)}{\varphi(d)}.
\]
Entonces las dos fórmulas exactas son
\[
K_A(h_0)=\frac1{\varphi(A)}
 \sum_{\substack{z_A\le s<y_A\\s\ {\rm primo}\\(s,AN)=1}}
 \frac{T_s(h_0)}{\varphi(s)}
\]
y, para la parte nueva s h_0 con s primo y (s,AN)=1,
\[
K_A(s h_0)=\frac{T_s(h_0)}{\varphi(A)\varphi(s)}.
\]
La segunda clase no puede omitirse: h puede contener el propio primo de switching.

Para cada primo s con z_A\le s<y_A y (s,AN)=1, definamos
\[
Q_s=\sum_{h_0\mid P(z_A)}\varphi(h_0)|T_s(h_0)|^2
=\sum_{d,e\in V^+(D_A/s)}
 \frac{\mu(d)\mu(e)(d,e)}{\varphi(d)\varphi(e)}.
\]
Separando las dos clases de h en la energía Gram,
\[
\begin{aligned}
\mathcal G_A
&=\frac1{\varphi(A)^2}
 \sum_{h_0}\varphi(h_0)
 \left(\sum_{\substack{z_A\le s<y_A\\s\ {\rm primo}\\(s,AN)=1}}
 \frac{T_s(h_0)}{\varphi(s)}\right)^2\\
&\quad+\frac1{\varphi(A)^2}
 \sum_{\substack{z_A\le s<y_A\\s\ {\rm primo}\\(s,AN)=1}}
 \frac1{\varphi(s)}
 \sum_{h_0}\varphi(h_0)|T_s(h_0)|^2.
\end{aligned}
\]
Por Cauchy ponderado, el primer término está acotado por
\[
\frac1{\varphi(A)^2}
 \left(\sum_{\substack{z_A\le s<y_A\\s\ {\rm primo}\\(s,AN)=1}}
 \frac1{\varphi(s)}\right)^2
 \sup_{\substack{z_A\le s<y_A\\s\ {\rm primo}\\(s,AN)=1}}Q_s,
\]
y el segundo por
\[
\frac1{\varphi(A)^2}
 \left(\sum_{\substack{z_A\le s<y_A\\s\ {\rm primo}\\(s,AN)=1}}
 \frac1{\varphi(s)}\right)
 \sup_{\substack{z_A\le s<y_A\\s\ {\rm primo}\\(s,AN)=1}}Q_s.
\]
Como z_A=X_A^{1/10} e y_A=X_A^{1/3},
\[
\sum_{\substack{z_A\le s<y_A\\s\ {\rm primo}\\(s,AN)=1}}
\frac1{\varphi(s)}
=\log(10/3)+o(1).
\]
En efecto, \(\sum_{z_A\le p<y_A}1/(p-1)=\log(10/3)+o(1)\). En el rango A\le N^{1-\epsilon} de trabajo, z_A\ge N^{\epsilon/10}; los primos excluidos por A están por debajo de z_A para N grande porque rad(A)|Q_N y los factores de Q_N son O_k(\log N). Además, N tiene a lo sumo \(\log N/\log z_A=O_\epsilon(1)\) factores primos distintos >=z_A, cuya suma recíproca es O_\epsilon(1/z_A)=o(1).
Por tanto, basta demostrar la cota uniforme de un solo bloque
\[
\boxed{\sup_{\substack{s\in[z_A,y_A)\\s\ {\rm primo}\\(s,AN)=1}}Q_s
\ll \frac1{\log z_A}.}
\]
Esta reducción incluye tanto las partes nuevas z_A-suaves como las que contienen un primo de switching. Como D_A=X_A^{1/2}/(log X_A)^{C_0}, el parámetro u_s=log(D_A/s)/log z_A recorre [5/3-o(1),4-o(1)], contenido en [3/2,4] para N grande. En el proxy numérico sin pérdida logarítmica, los extremos son 5/3 y 4.

### Experimento 94b [NUM]

Script: numerico/94b_norma_un_s.py. Se tomó A=3 y N=2^e, de modo que A es el selector canónico para k=2 y (A,N)=1. Para cada z se eligió un primo s próximo a D/z^u, con u≈5/3, 2, 3, 4. Se comparó Q_s con
\[
P(z)=\prod_{\substack{p<z\\p\nmid AN}}
\left(1-\frac2{p-1}+\frac{p}{(p-1)^2}\right),
\]
el producto que da la norma si se conservan todos los divisores cuadrados-libres y los signos mu(d), sin el recorte V^+.

Resultados:
- z=14: 1.0027 ≤ Q_s/P(z) ≤ 1.7006.
- z=20: 1.0029 ≤ Q_s/P(z) ≤ 1.5490.
- z=28: 1.0019 ≤ Q_s/P(z) ≤ 1.6193.
- z=40: 1.0046 ≤ Q_s/P(z) ≤ 1.5336.
- z=57: 1.0071 ≤ Q_s/P(z) ≤ 1.6759.

En los 20 casos, Q_s log z quedó entre 1.556 y 2.829. El cociente Q_s/P(z) se acerca a 1 cuando u≈4; el máximo observado aparece cerca de u=5/3. Esto apoya numéricamente, pero no demuestra, Q_s≪1/log z con constante uniforme en u.

**Siguiente sublema analítico:** comparar el recorte V^+(D_A/s) con el producto completo de Möbius en la norma gcd. No basta invocar la fórmula de término principal de la criba, porque esta controla sum_d lambda_d/phi(d), no el segundo momento Q_s.