# Doc 91 — LEMA UNIF-U, etapa 2: las cinco rutas medidas

Fecha: 2026-11-26. Estado: NINGUNA ruta da el lema. Goldbach sigue abierta.
Objetivo (doc 89, lema UNIF-U): R(N;Qsd)=sum_s sum_d lambda_d(s)·E(N;Qsd,N)
menor que el termino principal: R << (N/A)/log^2 N · o(1), A=Qu, u<=N^c.
Docs previos: SW-Q certificado para u=1 (docs 85-88); uniformidad delimitada
(doc 89: geometria hasta c<=2/17, balance hasta c*=0.039884); etapa 1 (doc 90):
faltan NIVEL (theta=1/2+c/2) y PEQUENEZ RELATIVA (A'(N)>=2+c·lnN/lnlnN;
BV da solo A' fijo). Cancelacion de raiz cuadrada bastaria solo hasta c<1/5.

## 0. El criterio maestro: ninguna cota absoluta sirve

Teorema (filtro). Sea B(N) << N/(ln N)^{A'} una cota ABSOLUTA del
resto, con A' fijo. El cociente al objetivo es
  R_abs = N^c / (ln N)^{A'-2} -> infinito para todo c > 0 fijo.
Medido (A'=10): log10 R_abs a N=10^12 es -11.29 (c=0.02) y -10.12
(c=2/17); a N=10^3600 es +40.65 y +392.18. El cociente primero CAE
(doc 90, Prop. 90.1: cruce solo en log10 N* >= 205) y luego crece
sin cota. Consecuencia: UNIF-U exige una cota RELATIVA — cancelacion
de raiz cuadrada en promedio, o una franja de ceros de ancho ~c.
Este filtro se aplica a las cinco rutas. [T]

## 1. BFI (Bombieri-Friedlander-Iwaniec) — NO da el lema

Teorema (BFI 1986, Thm 10; Maynard II, arXiv:2006.07088, Def. 1).
lambda_q es bien factorizable de nivel Q si para TODA factorizacion
Q=Q1Q2 existen gamma^(1), gamma^(2) con |gamma|<=1, soportadas en
[1,Qi], con lambda_q = sum_{q=q1q2} gamma^(1)_{q1} gamma^(2)_{q2}.
Si lambda es bien factorizable de nivel Q <= x^{4/7-eps} y a es
FIJO, entonces sum_{q<=Q} lambda_q E(x;q,a) << x/(log x)^A.
Conclusion ABSOLUTA. [T]

Nivel: cabe. El margen es 4/7 - (1/2+c/2) = 1/14 - c/2:
+0.0614 (c=0.02), +0.0515 (c*=0.0399), +0.0464 (c=0.05),
+0.0214 (c=0.10), +0.0126 (c=2/17). Maynard II Thm 1.1 (triple
bien factorizable) da 3/5 y Thm 1.2 (los pesos lambda^+ del tamiz
lineal SON bien factorizables) da 7/12: margenes mayores, mismo
problema de fondo.

Pesos: nuestros lambda_d(s) NO son bien factorizables en la
particion natural (s,d). Los pesos de Nathanson (Thm 9.3, beta=2)
tienen soporte V^+(D/s) que DEPENDE de s: medido a X=10^12,
|V^+| = 31 (s=17), 31 (s=19), 24 (s=101), 2 (s=9973); testigo
d=15 esta en V^+(D/17) pero no en V^+(D/9973); max d: 5005
(s=17) contra 3 (s=9973). Un lambda que depende de s no se
factoriza como gamma^(1)_s · gamma^(2)_d. La particion trivial
(A, sd) si factoriza (gamma = delta_A · h) pero es inutil para
dispersion: no separa nada. (El doc 58, Prop. 58.1, ya mostro que
lambda_t no es bien factorizable a ningun nivel; aqui se mide el
mecanismo exacto para lambda_d(s).)

Reparacion por bloques: en bloques diadicos [S,2S) la particion
natural se repara — el nivel es D/(2S) y el balance converge al
valor del doc 88: J_bloque = 1.790/1.710/1.586/1.540 a
X=10^6/10^12/10^24/10^100 (J=1.5238762503); Bal_bloque(Q=3) =
-0.013/0.027/0.088/0.112 -> 0.11965. La reparacion es
asintoticamente gratis (a X=10^6 el balance aun es negativo),
pero la conclusion de BFI sigue siendo ABSOLUTA: no da (I2).
Ademas a=N no es fijo. VEREDICTO: BFI NO da UNIF-U.

## 2. Maynard III — NO da el lema

Maynard III, arXiv:2006.08250. Thm 1.1: modulos q1q2 con
Q1 <= x^{1/10-3delta}/(log x)^C, Q2 <= x^{4/10+4delta}(log x)^C
(delta<1/1000), nivel total <= x^{0.501}; error solo O(delta)
sobre la cota trivial ("we are only winning a factor O(delta)
over the trivial bound"). Thm 1.2: Q1Q2Q3 = x^{1/2+delta}, ventana
x^{1/10+12delta}/Q2 < Q3 < x^{1/10-4delta}/Q2^{3/5}, error bueno
x/(log x)^A. Thm 1.3: minorante, no Lambda. Umbral 1/10 agudo
(doc 75: viene de Q^{15/2}R^{10} <= x^4, Lema 11.8). [T]

Medido: la fraccion de masa (s,d) con nivel total sobre x^{0.501}
es 0.030 (c=0.02), 0.083 (c=0.05), 0.178 (c=0.10): para todo
c>0.002 hay masa constante fuera del alcance. La ventana del
Thm 1.2 cubre solo el 2.75% de la masa (c=0.05). El error debil
del Thm 1.1 da R_debil = delta·N^c·(ln N)^2 = 3 a N=10^12,
48 a 10^24, 5.3e6 a 10^100 (c=0.05, delta=1/1000): >1 y
creciendo. El Thm 1.3 es el objeto equivocado (doc 64: 4.78% de
la masa del slice). VEREDICTO: Maynard III NO da UNIF-U.

## 3. BV con pesos bien factorizables — NO da el lema

Bombieri-Vinogradov vale para coeficientes acotados ARBITRARIOS:
la buena factorizabilidad no anade nada a su nivel (1/2). Con C0
creciente el exponente efectivo cae a 1/2-c (doc 90) y D_A solo
enciende desde N0 ~ 10^41-10^55. El nivel que ganaria la buena
factorizabilidad es el de BFI (4/7 con la particion reparada,
sec. 1), pero la conclusion de BV es absoluta con A' fijo: nunca
creciente. El cociente R_abs es el de la sec. 0 y diverge para
todo c>0. VEREDICTO: BV NO da UNIF-U, con o sin buena
factorizabilidad.

## 4. Dispersion Tipo II — NO da el lema (condicional hasta c<1/5)

El termino Q^2 de la criba grande a Q=N^{1/2+c/2} es Q^2/N = N^c:
medido a N=10^12: 1.74 (c=0.02), 3.98 (c=0.05), 15.85 (c=0.10),
25.81 (c=2/17) — exactamente el factor que falta. La celda
balanceada M~K~x^{1/2} esta fuera del alcance uniforme
(Maynard III sec. 5; doc 64). [T]

Segundo momento fijo-en-a: sum_{q<=Q} |E(N;q,N)|^2 << N (ln N)^C
NO es teorema — Barban-Davenport-Halberstam promedia sobre a; la
version a fija es el modelo aleatorio. SI valiera, Cauchy daria
cociente R_aleat = N^{5c/4-1/4} (ln N)^{C/2+2} -> 0 sii c<1/5:
medido (C=5), log10 R_aleat a c=0.05: 4.24/3.34/-8.12/-25.52
(N=10^12/10^24/10^100/10^200); a c=2/17: 5.25/5.37/0.34/-8.60.
O sea: aun condicional, solo hasta c<1/5, y a escala computable
el cociente aun es >1. GRH por q tampoco rescata: ratio
N^{3c/2}(ln N)^4, log10 = 6.67/8.77/16.95 a N=10^12/10^24/10^100
(c=0.05): >>1 y creciendo. VEREDICTO: Tipo II NO da UNIF-U;
deja el sub-lema exacto (sec. 6).

## 5. Densidad de ceros — NO da el lema

Lo que haria falta: una franja libre de ceros de ancho ~c fijo (sec. 0),
o una estimacion de densidad que diera una cota relativa. Lo que se tiene:
region de Vinogradov-Korobov (Ford 2002, 1/57.54; asintotico 1/49.13;
Mossinghoff et al. 2024: 1/48.1588): ancho w(N)=1/(57.54 (ln N)^{2/3}
(ln ln N)^{1/3}) = 0.00127 (N=10^12), 0.00075 (10^24), 0.00026 (10^100).
El cociente c/w es 15.7/26.5/76.0 (c=0.02) y 39.2/66.3/190.1 (c=0.05):
la franja necesaria es 16 a 190 veces mas ancha que la conocida, y el
ancho conocido se ENCOGE como (log N)^{-2/3} mientras c es fijo.
Para L(s,chi) hay ademas un termino 1/log q (Kahale): a N=10^12,
w_L = 0.00102, aun mas estrecha. Guth-Maynard 2024 (Annals 2026, vol. 203):
N(sigma,T) << T^{30(1-sigma)/13+o(1)}, mejora de Huxley 12/5=2.4 a
30/13=2.3077 — pero la densidad alimenta cotas ABSOLUTAS (Gallagher),
no relativas: no toca el filtro de la sec. 0. VEREDICTO: NO da UNIF-U.

## 6. Veredicto y sub-lema que falta

Las cinco rutas medidas: BFI no (conclusion absoluta + a=N no
fijo, aunque el nivel cabe y los bloques reparan la
factorizabilidad); Maynard III no (nivel 0.501 + error debil +
ventana 2.75%); BV no (nivel 1/2, conclusion absoluta); Tipo II
no (Q^2/N = N^c; segundo momento a fija no es teorema; GRH no
rescata); densidad de ceros no (franja 16-190x mas estrecha de
lo necesario; GM2024 alimenta cotas absolutas). El cuello de
botella es UNICO y es el filtro de la sec. 0: toda cota absoluta
falla para c>0; solo una cota relativa puede dar UNIF-U.

LEMA UNIF-U-2M: sum_{q<=Q} |E(N;q,N)|^2 << N (log N)^C para
Q = N^{1/2+c/2} a a=N fijo. De el seguiria, por Cauchy, la
cancelacion de raiz cuadrada hasta c<1/5. Es el sub-lema exacto
que falta; no es teorema conocido.

Script: numerico/91_lema_unif_u_etapa2.py — PASS 9/9 local y en
el PC, numeros identicos. Goldbach sigue abierta.

Fuentes: Nathanson caps. 9-10; Chen 1973; Bombieri-Friedlander-
Iwaniec 1986 Thm 10; Maynard II arXiv:2006.07088; Maynard III
arXiv:2006.08250; Fouvry-Radziwill arXiv:1811.08672; Lichtman
arXiv:2309.08522; Ford 2002; Guth-Maynard arXiv:2405.20552;
Montgomery; Gallagher; docs 58, 64, 65, 75, 87-90.
