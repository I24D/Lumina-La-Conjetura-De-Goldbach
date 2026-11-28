# 88 - Lema SW-Q apartado (c): la constante c de (85.4) y la cota de la cola R_y
## Procedencia y estado
- **Agente:** LUMINA/Muse. **Fecha:** 2026-11-26 16:40 ET (reclamo) / entrega.
- **Metodo:** lectura de la fuente primaria (Nathanson, *Additive Number Theory: The Classical Bases*, GTM 164 (1996), cap. 10, Teorema 10.2; Bordignon-Johnston-Starichkova, arXiv:2207.09452, Lema 37, Teorema 51, Lema 52) y calculo directo (integracion numerica propia de F/f y de las integrales de la cola; fuerza bruta finita con pesos Lambda). Version de limite FIJO z certificada en el doc 87 (la variable del doc 85 quedo abandonada). Fuentes de control: Halberstam-Richert, *Sieve Methods* (1974); Baxter, notas de Columbia sobre Nathanson cap. 10 (control secundario); Kuberwastaken, notas sobre Chen (control secundario, peso W(n)).
- **Estado:** **[T]** apartado (c) CERTIFICADO con las correcciones de las secs. 1-4. **Lema SW-Q COMPLETO (a)+(b)+(c)** (salvedad honesta en sec. 6: el resto bilineal de la cola se transfiere del Lema 32 de Bordignon con Q insertado; su re-escritura completa queda [P-rev]). El balance es POSITIVO a z = X^{1/10} para todo Q, luego J_k(N;Q_N) > 0 para N grande: C_k + C_{k+1} > 0 para la subincidencia u = 1. Es el teorema de Chen en el marco Q-torcido (del tipo Lu-Cai, Chen en progresiones aritmeticas), verificado/transferido aqui, no un teorema nuevo; no separa P_1 de P_2 (paridad), no toca la I_k completa (u > 1) y no es Goldbach.

## 1. La constante c = 1/2 (y donde multiplica)
La fuente (Nathanson Teorema 10.2 = Bordignon Lema 37, con z = N^{1/8} ahi):
```text
pi_2(N) > S(A,P(z)) - (1/2) sum_{z<=q<y} S(A_q,P(z)) - (1/2) S(B,P(y)) - 2N^{7/8} - N^{1/3},
```
A = {N-p : p <= N, p no | N} (conteos planos), B = {N-p1 p2 p3 : z<=p1<y<=p2<=p3, p1 p2 p3 < N, (p1 p2 p3, N) = 1}. El 1/2 aparece DOS veces: ante la suma de switching y ante la cola. El mecanismo es el peso
```text
W(v) = 1 - (1/2) #{p|v : z<=p<y} - (1/2) #{(p1,p2,p3) : v = p1 p2 p3, z<=p1<y<=p2<=p3},
```
donde #{p|v} cuenta primos DISTINTOS que dividen a v, y las descomposiciones P3 se cuentan como multiconjuntos (una por v si Omega(v) = 3). Si (v,P_Q(z)) = 1 y W(v) > 0 entonces Omega(v) <= 2: un P3 con p1 en [z,y) y p2,p3 >= y tiene W = 1 - 1/2 - 1/2 = 0; un P3 con p1^2 p3 (p1 en [z,y), p3 >= y) tiene W = 1/2 pero es una fuga contada en sec. 2, no un P3 conmutado.

**Prop. 88.1 [T] (identidad exacta del peso).** Con a_Q(v) = Lambda(N-Qv), S_Q(A;z) el conteo cribado (limite fijo z), B(Q) = sum_{z<=s<y} S_Q(A_s;z) (limite fijo z, doc 87) y Omega_Q = sum_{v = p1 p2 p3, z<=p1<y<=p2<=p3} a_Q(v):
```text
sum_{v} a_Q(v) W(v) = S_Q(A;z) - (1/2) B(Q) - (1/2) Omega_Q.     (88.1)
```
EXACTA (sin resto): el primer termino es S_Q(A;z) por definicion; el segundo, porque cada par (v,s) con s|v, s en [z,y), se cuenta una vez en S_Q(A_s;z) — y el limite FIJO z es lo que lo hace exacto: si (v,P_Q(z)) = 1 y s|v, entonces (v/s,P_Q(z)) = 1 automaticamente (con limite variable s el cofactor tendria que ser cribado por P_Q(s), distinto). El tercer termino es Omega_Q por definicion. (H1: exacta a 1e-15 relativo en 4 casos.)

**Prop. 88.2 [T] (c = 1/2, y NO multiplica a S).** La forma (85.4) del doc 85, "c * [S_Q(A;z) - B(Q)]", es INCORRECTA en la forma: el 1/2 multiplica SOLO a B(Q) y a la cola, no a S_Q(A;z). La correccion es (88.4) en sec. 4. El 1/2 viene del peso: cada P2 con p2 < y se cuenta DOS veces en B(Q) (s = p1 y s = p2, limite fijo; H2: 34/115 casos "twice" a N = 10^6/10^7, Q = 3), y cada P3 conmutado se resta UNA vez con peso 1/2.

**Prop. 88.3 [T] (el objetivo es J_k^{P<=2} = J_k, no J_k^{P2} solo).** El peso vale 1 en primos y 1/2 en P2's con p1 en [z,y); la cota inferior que da el switching es para la suma P1+P2, y NO se puede separar: eso es la barrera de paridad, no un defecto de la prueba. Para el objetivo del bucle (C_k + C_{k+1} > 0, doc 83) eso es exactamente lo necesario, porque J_k = J_k^{P1} + J_k^{P2} es la union omega(v) en {1,2} de (85.1). La (85.4) pedia J_k^{P2} solo; lo que la fuente da (y lo que se certifica) es J_k.

## 2. La cola R_y: lo que realmente se resta
**Prop. 88.4 [T] (la R_y del doc 85 NO se resta).** La cola R_y del doc 85 (P2's con p1 >= y) no necesita ninguna resta: esos P2's no tienen divisor en [z,y), luego aparecen CERO veces en B(Q) y entran con peso +1 en S_Q(A;z) - (1/2) B(Q) (W(v) = 1). Lo que la fuente resta (con 1/2) es la cola conmutada
```text
(1/2) S_Q^Lambda(B, P(y)),  B_Q = {N - Q p1 p2 p3 : z<=p1<y<=p2<=p3, p1 p2 p3 < X},
```
es decir, los P3's de A contados cribando B_Q por primos (switching / role reversal).

**Prop. 88.5 [T] (switching de la cola).** Omega_Q = sum_{b in B_Q} Lambda(b) EXACTO (H2b: acuerdo 1e-11): la biyeccion v <-> b = N - Qv envia cada P3 v con a_Q(v) != 0 a un b en B_Q con Lambda(b) = Lambda(N-Qv), y viceversa. Los b primos son la parte principal (H2b: 12657.3 de 12736.1 a N = 10^6, Q = 3); la parte de potencias de primos b = p^k es pequena (78.8). Luego Omega_Q <= S_Q^Lambda(B,P(y)) := sum_{b in B_Q, (b,P(y)) = 1} Lambda(b), porque los b primos >= y son coprimos con P(y).

**Prop. 88.6 [T] (errores del peso E1'').** Los v con W(v) > 0 y Omega(v) >= 3 (fugas, fuera de Omega_Q) son solo tres familias: (i) cubos de primos v = p^3 (W = 1/2 o 1), a lo sumo X^{1/3} valores (H6: 5-7 a N <= 10^7); (ii) P3's con factor cuadrado v = p1^2 p3, p1 en [z,y), p3 >= y (W = 1/2, fuera de Omega_Q porque su multiconjunto {p1,p1,p3} no cumple p2 >= y), de tamano X/(z log z) * log N por conteo directo; (iii) v pares con N - Qv = 2^k (Q impar, N par: N-Qv par solo si v par), O(log N) valores, total O(log^2 N). Envolvente total O(N^{9/10} * polylog) = O(N/z * polylog) << N/log^2 N. Medido (H6b): E1'' = 4661.2 (1.40% de S) a N = 10^6, 40785.2 (1.22%) a N = 10^7, Q = 3; E1''/N^{9/10} ~ 2% en ambos casos.

## 3. La cota de la cola (el c de Nathanson, explicito)
**Prop. 88.7 [T] (tamano de B_Q).** |B_Q| <= (X/log X) * [c_{1/10} + o(1)] con
```text
c_{1/10} = int_{1/10}^{1/3} log(2-3b)/(b(1-b)) db = 0.4909952010.     (88.2)
```
Derivacion (Lema 52 de Bordignon, con z = X^{1/10} en vez de N^{1/8}): |B_Q| se escribe con I(u) = int_y^{sqrt(N/u)} h_u(t) d log log t, h_u(t) = 1/log(N/ut); el cambio b = log p1/log X, tau = log t/log X da la integral doble que colapsa a log(2-3b)/(b(1-b)) (la integral interior es [1/(1-b)] log(2-3b)). Control publicado reproducido (H3): cbar = int_{1/8}^{1/3} log(2-3b)/(b(1-b)) db = 0.3630837292 < 0.363084 (Lema 52 de Bordignon).

**Prop. 88.8 [T] (criba superior de la cola).** S_Q^Lambda(B,P(y)) <= |B_Q| * V_Q(y) * F(s_b) + R, con s_b = log D^{(3)}/log y, D^{(3)} = N^{1/2-alpha_3}, s_b -> 3/2, F(s_b) = 2e^gamma/s_b -> (4/3) e^gamma. El resto R es la forma bilineal del Lema 32 de Bordignon (valor medio sobre la convolucion {p2 p3} con primos p1, nivel D* = sqrt(XY)/log^{10} Y) con Q insertado: los modulos pasan a ser Qd <= Q * N^{1/2-alpha_3} <= N^{1/2} para N grande, porque Q << log^{k-1} N y N^{alpha_3} domina cualquier polilog; los 10 logs de holgura de D* absorben Q. **[T]** por transferencia del mismo argumento; la re-escritura completa del Lema 32 con Q queda [P-rev] (ver sec. 6).

**Prop. 88.9 [T] (coeficiente de la cola).** En unidades de |A| V_Q(z):
```text
S_Q^Lambda(B,P(y)) <= 0.4 * e^gamma * c_{1/10} * (phi(Q)/Q) * |A| V_Q(z).   (88.3)
```
Derivacion: |B_Q| ~ (X/log X) c; V_Q(y) = 3 U_N prod_{l|Q}(l-1)/(l-2) / log N; F = (4/3) e^gamma; |A| V_Q(z) = (N/phi(Q)) * 10 U_N prod/log N. El log N del peso Lambda cancela el 1/log N de V_Q(y); el factor phi(Q)/Q viene de |B_Q| ~ N/(Q log N) frente a |A| ~ N/phi(Q). A z = N^{1/8} el mismo calculo da 0.5 * e^gamma * cbar (V(y)/V(z) = 3/8 en vez de 3/10), y el balance resultante reproduce la constante publicada de Chen (H5_chen, sec. 5).

## 4. La (85.4) corregida y el balance
**Prop. 88.10 [T] ((85.4) corregida).**
```text
R_Q^{P<=2} >= S_Q(A;z) - (1/2) B(Q) - (1/2) S_Q^Lambda(B,P(y)) - O(N^{9/10} polylog),
J_k >= (log Q + log z) * R_Q^{P<=2} - log z * log N.                  (88.4)
```
R_Q^{P<=2} = sum_{v : Omega(v)<=2} a_Q(v). La primera linea es el peso (88.1) + switching (Prop. 88.5) + errores (Prop. 88.6); la segunda es el peso minimo (H7): log Q + log rad(v) >= log Q + log z para todo v >= z con (v,P_Q(z)) = 1 (v primo >= z, o P2 con rad(v) >= z^2), y el termino v = 1 cuesta a lo sumo log z * log N. El desdoblamiento fino del peso log rad(v) del doc 85 sec. 2 [P-rev] NO hace falta para positividad: el peso minimo basta.

**Prop. 88.11 [T] (balance positivo a z = X^{1/10}).** En unidades de |A| V_Q(z): S_Q(A;z) >= f(5) con f(5) = 0.9982417256 (H4; la criba inferior con s = log D/log z = 5); (1/2) B(Q) <= (1/2) J con J = int_{1/10}^{1/3} F(5-10b) db/b = 1.5238762503 (H4; doc 87, M_B(Q) de (87.7)); (1/2) cola <= (1/2) * 0.4 e^gamma c_{1/10} (phi(Q)/Q) (Prop. 88.9). Balance:
```text
Q = 1: 0.061404;  Q = 2: 0.148854;  Q = 6: 0.178004;  Q = 30: 0.189664.   (88.5)
```
Positivo en todos los Q; el peor caso es Q = 2 (k = 2, el selector mas pequeno). Luego R_Q^{P<=2} >= (balance) * |A| V_Q(z) - O(N/log^A N) > 0 para N grande, y por (88.4) J_k(N;Q_N) > 0: existe v con omega(v) en {1,2}, (v,Q) = 1, N - Qv primo (coprimo con N salvo O(log N) excepciones). Es C_k + C_{k+1} > 0 para la subincidencia u = 1 (doc 83). Control: a z = N^{1/8}, Q = 1, el mismo calculo da balance 0.0188710759, que en unidades de la literatura (x 4/e^gamma) es 0.0423813781 = 2 log 3 - log 6 - cbar (acuerdo 7.6e-10, H5_chen): la maquinaria reproduce la constante de Chen publicada.

## 5. Verificacion numerica (control finito, no evidencia asintotica)
Script `numerico/88_lema_swq_parte_c.py` (solo stdlib, como docs 85-87; ~6 s local):
```text
PASS H4_f5 | f(5)=0.998241725611 grid vs 0.998241725611 via int_3^4 F; F(4)=1.021641549
PASS H4_J | J=1.523876250329 simp vs 1.523876250382 trap; J(1/8)=1.595626688122 vs (e^g/2)log6=1.595626685157
PASS H3_c | c_{1/10}=0.490995201004 (simp 0.490995201668 trap); cbar=0.363083729248 < 0.363084
PASS H5_balance | balance z=X^{1/10}: Q=1:0.061404 Q=2:0.148854 Q=6:0.178004 Q=30:0.189664
PASS H5_chen | z=N^{1/8} balance=0.018871075894; x4/e^g=0.042381378103 vs 2log3-log6-cbar=0.042381378860
PASS H1_identidad_N1000000_Q3 | S=333577.6472 B=250517.6617 Om=12736.1051 SW=201950.7639 rel=2.36e-15
PASS H2_multiplicidad_N1000000_Q3 | P2: twice=34 once=7881 zero=4149 (B_s=17)
PASS H2b_switching_N1000000_Q3 | Om=12736.1051 sum_Lambda(B_Q)=12736.1051 prime=12657.3102 pp=78.7948
PASS H6_errores_N1000000_Q3 | cubos=7 (<=X^{1/3}=69.3); fugas=740 E1''=4661.2036 (0.0140 de S; <=N^0.9=251188.6)
PASS H7_peso_minimo_N1000000_Q3 | Jk=3244344.0677 Rle2=254226.0042 bound=558576.4466
PASS H1/H2/H2b/H6/H7 | N=2*10^6 Q=3: S=666580.2114 B=535344.7118 Om=26833.4688 SW=385491.1211 rel=1.57e-15; twice=53 once=14997; E1''=9437.2941 (0.0142 de S); Jk=6629889.0545 bound=1082835.5419
PASS H1/H2/H2b/H6/H7 | N=10^7 Q=3: S=3332306.6485 B=3003875.6405 Om=169974.2418 SW=1745381.7073 rel=2.31e-15; twice=115 once=65628; E1''=40785.1501 (0.0122 de S); Jk=34785744.4702 bound=5734314.9084
PASS H1/H2/H2b/H6/H7 | N=10^6 Q=77: S=12810.9137 B=9123.4825 Om=710.6580 SW=7893.8435 rel=5.68e-16; twice=3 once=282; E1''=894.0009 (0.0698 de S); Jk=112818.8844 bound=44668.2652
PASS
```
Lectura honesta: H1/H2/H2b/H6/H7 son exactos por fuerza bruta (Q IMPAR: con Q par, N-Qv es par y a_Q es casi identicamente cero; el selector Q_N del doc 83 es coprimo con N). H3/H4/H5 son calculos numericos de las constantes analiticas con dos metodos independientes (Simpson por tramos en los quiebres C^2 de F en u = 2,3; trapecio fino; malla de retardo). F(4) = 1.021641549 (malla M = 3000; doc 87: 1.021641553, acuerdo 4e-9). Los numeros deben ser IDENTICOS local y en el PC.

## 6. Estado del Lema SW-Q y siguiente paso
- **(a) CERTIFICADO** (doc 86). **(b) CERTIFICADO**, limite fijo z (doc 87). **(c) CERTIFICADO** con las correcciones de las secs. 1-4: c = 1/2 (multiplica B(Q) y la cola, no S); la cola que se resta es (1/2) S_Q^Lambda(B,P(y)) con c_{1/10} = 0.4909952010 (no la R_y del doc 85); el objetivo certificado es J_k (P1+P2), no J_k^{P2} solo; balance positivo (88.5). **Lema SW-Q COMPLETO**, con una salvedad: el resto bilineal de la cola se transfiere del Lema 32 de Bordignon con Q insertado por el mismo argumento (Prop. 88.8); su re-escritura completa queda **[P-rev]** — es escritura, no un hueco matematico (los modulos Qd caben en el nivel con holgura polilogaritmica).
- Queda abierto: **[N]** el desdoblamiento fino del peso log rad(v) para una cota CUANTITATIVA de J_k (para positividad basta el peso minimo, H7); **[N]** la I_k completa (u > 1, uniformidad en A = Qu, doc 85 sec. 3); **[N]** separar P1 de P2 (paridad: el switching no lo hace). Goldbach sigue abierta.
- **Siguiente paso natural:** la uniformidad de la incidencia completa (u > 1) o la instanciacion del bucle (b_{P,l}, Motohashi doc 57 sec. 4). No se empieza aqui.
