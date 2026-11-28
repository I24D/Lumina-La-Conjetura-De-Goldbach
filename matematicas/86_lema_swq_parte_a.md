# 86 - Lema SW-Q apartado (a): pesos de Chen con Q insertado
## Procedencia y estado
- **Agente:** LUMINA/Muse. **Fecha:** 2026-11-26 14:21 ET (reclamo) / 14:40 ET (entrega).
- **Metodo:** lectura peso a peso de Nathanson, *Additive Number Theory: The Classical Bases* (GTM 164, 1996), cap. 9 (secs. 9.1-9.4) y cap. 10 (secs. 10.2-10.5), desde el texto. Fuentes secundarias de control: Yamada, arXiv:1511.03409 (criba lineal de Rosser-Iwaniec explicita, sec. 3); Bordignon, Bull. Aust. Math. Soc. 105 (2022), 344-346.
- **Estado:** **[T]** apartado (a) CERTIFICADO, con dos refinamientos (coeficiente afilado 1; forma exacta del resto) y una correccion de enunciado (C0 depende de A). **[N]** (b) y (c) no se intentan aqui. No se demuestra Goldbach; no hay victoria parcial.

## 1. Los pesos: Teorema 9.3 (Nathanson cap. 9, sec. 9.2, beta=2)
Con beta > 1 y nivel D > 0, sean V+ (resp. V-) el conjunto de 1 y los
d = p1...pk libres de cuadrados con p1 > ... > pk, pk < ... < p1 < D,
y p_m < (D/(p1...pm))^{1/beta} para todo m impar (resp. par). Los pesos
de la criba lineal superior (resp. inferior) son
```text
lambda^+(d) = mu(d) si d en V+ y d | P(D);  0 si no.
lambda^-(d) = mu(d) si d en V- y d | P(D);  0 si no.      (86.1)
```
Para el resto del cap. 9, Nathanson toma beta = 2 (sec. 9.3).

**Prop. 86.1 [T].** `|lambda^+(d)| <= 1` para todo d; soporte d < D,
d libre de cuadrados, d | P(D); lambda^+(1) = 1. La cota es inmediata
de (86.1): los pesos son mu(d) o 0. Ademas es hipotesis EXPLICITA del
Teorema 9.2 (p. 237-238): "Let |lambda^+(d)| <= 1 for all d >= 1", y
la prueba del Teorema 9.7 (Jurkat-Richert, p. 258) lo repite: "there
exist upper and lower bound sieves lambda^+(d) ... with
|lambda^+(d)| <= 1 for all d >= 1". La extension del Lema 9.1 (rango
P1 u Q, nivel DQ) preserva la cota: la descomposicion d = d1*d2 es
unica y lambda~(d) = lambda(d1)*mu(d2), luego |lambda~(d)| <= 1.

**Conclusion.** Los pesos lambda_d(s) de (85.6) — criba lineal superior
con nivel D/s y limite de criba s (o z fijo, ver sec. 6) — cumplen
`|lambda_d(s)| <= 1` [T]. La primera mitad de (a) esta certificada.

## 2. La forma exacta del resto (correccion de (85.6))
En el Teorema 9.7 el resto es suma de valores absolutos,
R = sum_{d|P(z), d<DQ} |r(d)|, con r(d) = |A_d| - |A|*g(d). Para A_s,
con g(d) = 1/phi(d) (Nathanson sec. 10.3) y (84.2):
|A_{s,d}| = psi(N;Qsd,N) - Lambda(N), |A_s| = psi(N;Qs,N) - Lambda(N),
y phi(Qsd) = phi(Qs)*phi(d) (Q, s, d coprimos dos a dos).

**Prop. 86.2 [T] (resto exacto).**
```text
r_s(d) = E(N;Qsd,N) - E(N;Qs,N)/phi(d) - Lambda(N)*(1 - 1/phi(d)).
                                                              (86.2)
```
La forma (85.6) del doc 85 es solo el PRIMER termino. Los otros dos
son correcciones, y ambas son despreciables por BV:
- **C2** = E(N;Qs,N) * sum_d lambda_d(s)/phi(d): |sum_d lambda_d(s)/phi(d)|
  <= sum_{d<D/s} 1/phi(d) << log N (Teorema A.17, citado en la prueba
  del Teorema 10.5). Y sum_s |E(N;Qs,N)| <= sum_{m<=Qy} max_a |E| << N/log^{A+1} N
  por BV, pues Qy << N^{1/3}*polylog. Coste: un logaritmo, absorbido
  tomando A+1 en BV.
- **C3** = Lambda(N) * sum_s sum_d lambda_d(s)(1-1/phi(d)): |.| <=
  Lambda(N) * sum_s (D/s) << D * log log y * log N << N^{1/2}*polylog
  << N/log^A N.
Nathanson hace exactamente esto en la prueba del Teorema 10.5
(p. 282-283): r_q(d) = r(qd) - r(q)/phi(d), y estima los dos
sumandos por separado. La forma firmada (85.6) esta dominada por la
forma en valor absoluto, asi que basta acotar sum_s sum_d |r_s(d)|.

## 3. El coeficiente afilado: cada modulo aparece una sola vez
**Prop. 86.3 [T].** En sum_s sum_d lambda_d(s) E(N;Qsd,N), cada modulo
m = Qsd aparece a lo sumo UNA vez. Prueba: si s1 < s2 en [z,y)
dividen m/Q, entonces s2 | d1 = m/(Qs1); pero d1 | P_Q(s1) exige
primos < s1 < s2, contradiccion. Luego s es el mayor factor primo de
m/Q en [z,y). Por tanto
```text
|sum_s sum_d lambda_d(s) E(N;Qsd,N)| <= sum_{m<=QD} |E(N;m,N)|
    <= sum_{m<=QD} max_a |E(N;m,a)|.                         (86.3)
```
El <= 5 del doc 85 sec. 3(b) es valido pero no afilado; el afilado
es 1. Esto es lo que hace Nathanson en el Teorema 10.5 (p. 283):
sum_q sum_{d<QD_q} |r(qd)| <= sum_{d'<QD} |r(d')|.

## 4. Nivel, BV y la correccion de C0
El nivel del resto en el Teorema 9.7 es DQ, con Q el producto del
conjunto excepcional para (9.34): en Nathanson sec. 10.3 son los
primos < u0(epsilon), con Q < log N (10.7). Luego m = Qsd < QD*log N
<= N^{1/2}/log^{C-(k-1)/2-1} N. La condicion de densidad (9.34) con
g(p) = 1/(p-1) se tiene con K = 1+epsilon (Nathanson sec. 10.3, via
Teorema 6.9); excluir los primos que dividen QN solo disminuye el
producto. La pre-criba (sec. 5) anade un factor Q^2, absorbido igual.

**CORRECCION DE ENUNCIADO [T].** El Lema SW-Q del doc 85 sec. 5 dice
"D = X^{1/2}/log^{C0} X con C0 absoluto" y concluye "<< N/log^A N
para todo A". Eso es inconsistente: BV exige nivel N^{1/2}/log^{B(A)} N
con B(A) creciente en A; con C0 absoluto solo se obtiene la cota
hasta un A0(C0,k). Correccion: **C0 = C0(A,k)** (por ejemplo
C0 = B(A)+k+2 cubre el log N de Q_Nath, el Q^{1/2} y el Q^2 de la
sec. 5). Con esa correccion, (a2) se tiene por BV (84.3).

**Veredicto (a): CERTIFICADO** con C0 = C0(A,k): |lambda_d(s)| <= 1
(Prop. 86.1) y sum_s sum_d lambda_d(s) E(N;Qsd,N) <<_{A,k} N/log^A N
(Props. 86.2-86.3 + BV + correcciones C2, C3).

## 5. La pre-criba (w,Q)=1 (cierra doc 85 sec. 3(d))
**Prop. 86.4 [T].** Con (w,Q)=1 y (d,Q)=1,
```text
|A^Q_{s,d}| = sum_{e|Q} mu(e) [psi(N;Qsde,N) - Lambda(N)].   (86.4)
```
Es Mobius sobre e|Q: son 2^{k-1} terminos de la misma forma, con
modulos <= Q^2*D*log N, todos BV-acotados. La forma (85.6) no cambia;
solo se multiplica el numero de terminos por 2^{k-1}. El detalle
[P-rev] del doc 85 sec. 3(d) queda cerrado.

## 6. Diferencia de montaje con Nathanson (nota para (b), no se resuelve aqui)
Nathanson criba A_q con limite FIJO z = N^{1/8} (Teorema 10.5:
s_q = log(D/q)/log z en [4/3,3]); el doc 85 criba A_s con limite
VARIABLE s. Para (a) no importa: los pesos existen para cualquier
nivel y limite (Teorema 9.3) y el resto es el de Teorema 9.2, sin
condicion D >= z. Para (b) SI importa: el Teorema 9.7 exige D >= z
(s >= 1) para la cota superior, y en la version del doc 85 el nivel
D/s cae por debajo del limite s cuando s > D^{1/2} ~ X^{1/4}, es
decir, en buena parte de [z,y). Ahi u_s < 1 y F(u_s) no esta
controlada por los Teoremas 9.5-9.7. La version de Nathanson (limite
fijo) evita ese regimen. La transferencia del termino principal
M_B(Q) de (85.5) debe decidir entre ambas versiones; eso es (b).

## 7. Verificacion numerica (control finito, no evidencia asintotica)
Script `numerico/86_lema_swq_parte_a.py`: N par en [100,2000], k=2,3;
H2/H3/H5 en submuestra N%10==0 mas N en {128,256,512,1024};
D_s = floor((X/s)^{1/2}) (proxy finito sin log^C, como doc 85).
```text
H1_Th93: 3/3 (pesos mu o 0, |lam|<=1, (9.8) para m|P(D), m<=200)
H4_densidad_934: max_R=1.022 (ilustracion finita, no prueba)
H2_coef_afilado: 28/28 modulos, max_pares=1 (peor (1024, 2, 15, [(5, 1)]))
H3_resto_exacto: 28/28, max|E1|=146.812 max_corr2=146.812 max_corr3=0.000
H5_precriba_Q: 28/28 (Mobius sobre e|Q exacto)
PASS
```
Lectura honesta: H1/H2/H3/H5 pasan al 100% (construccion de pesos,
multiplicidad 1, identidad exacta (86.2), Mobius de la pre-criba).
H4 mide el cociente de Mertens (max 1.022 a z=100), consistente con
R -> 1; no prueba (9.34). A N<=2000 los errores son de orden del
termino principal (max|E1|=146.8): el script controla la FORMA, no
el tamano asintotico; el tamano es BV (analitico). max_corr3=0.000
en la submuestra porque d=1 domina (phi(1)=1 anula el termino) y
Lambda(N)=0 salvo N potencia de primo; el termino esta en la
identidad y se verifica como 0. Los numeros deben ser IDENTICOS
local y en el PC.

## Siguiente paso
(b) Termino principal M_B(Q) con F(u_s) explicita — con la nota de
la sec. 6: decidir entre limite variable s (doc 85) y limite fijo z
(Nathanson), y tratar el regimen u_s < 1. (c) Constante c de (85.4)
y cota de la cola R_y. Goldbach sigue abierta.
