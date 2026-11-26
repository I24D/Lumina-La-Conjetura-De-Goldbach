# 69 - Lema del momento mixto, via Vaughan/Heath-Brown de S_rho (Etapa 1)

Fecha: 2026-11-26. Etapa 1 de 2 (solo esta etapa): el lema del momento
mixto del doc 67 sec. 3, atacado por UNA sola via: la descomposicion de
Vaughan/Heath-Brown de S_rho. Script: `numerico/69_momento_mixto_vaughan.py`
— PASS local (X = 2x10^6). No hay prueba del lema, no hay cota inferior
de capas, no hay demostracion de Goldbach.

## 1. El lema, con total precision

Sea N par, S_1(alfa) = sum_{p<N} e(p*alfa), S_rho(alfa) = sum_n rho(n)
e(n*alfa) con rho el minorante de Harman (doc 66, diadico a escala N).
Arcos mayores M: union sobre q <= Q_0 = (log N)^B (B fijo grande) de
[a/q - 1/(qQ_1), a/q + 1/(qQ_1)], Q_1 = N/(log N)^B; arcos menores
m = [0,1] \ M. El termino principal del conteo unilateral es
c_rho * S(N) * N/log^2 N (doc 67, Prop. 67.1).

**Lema MM-69.1 (momento mixto).** int_m |S_1(alfa) S_rho(alfa)| d alfa
= o(N/log^2 N).

Es la condicion suficiente estandar para que el error de arcos menores
del conteo unilateral quede bajo el termino principal. Nota de alcance:
el lema es suficiente, no necesario: el error con fase
int_m S_1 S_rho e(-N alfa) podria ser pequeno aunque la integral
absoluta no lo sea. Lo que se prueba aqui es que la via de
descomposicion NO puede establecer el lema, no que el lema sea falso.

## 2. La via: Vaughan/Heath-Brown de S_rho

Para Lambda, la identidad de Vaughan/Heath-Brown descompone S_Lambda en
piezas Tipo I (convolucion con factor suave largo) y Tipo II (bilineal
con ambos factores en rango medio); asi se prueba la cota de Vinogradov
sup_m |S_Lambda| << N/(log N)^A. Para el minorante rho, el analogo es
su descomposicion de Buchstab (Maynard III sec. 9, doc 66): rho =
S^triv + S^G - C*rho_B, y en el proxy sin parte B (doc 67)

  rho_tilde(n) = 1_p(n) - S_{2,3}(n) - S_{2,5}(n),

luego como sumas exponenciales S_rho_tilde = S_1 - S_{2,3} - S_{2,5}.
Cada pieza es una convolucion Tipo I o Tipo II: S_{2,3}, S_{2,5} son
bilineales (dos factores primos en [z_1, z_2] = [N^{1/7}, N^{3/7}] con
cofactor rugoso). La pregunta de la etapa 1: estimando cada pieza por
separado — por sup, por Cauchy, o por cuarto momento — se alcanza el
o(N/log^2 N)? Respuesta: no, por tres razones independientes (Props.
69.1-69.3). Las mediciones son a X = 2x10^6 (misma escala del doc 67).

## 3. Prop. 69.1: la barrera de Cauchy por piezas [T]

Sea S_rho = sum_{j=1}^J T_j cualquier descomposicion finita en piezas
(en particular la de Buchstab). Por Minkowski en L^2(d alfa) y Parseval,

  sum_j ||T_j||_2 >= ||S_rho||_2,

luego la estimacion de Cauchy por piezas cumple

  sum_j int_m |S_1||T_j| <= sum_j ||S_1||_2 ||T_j||_2
     = ||S_1||_2 sum_j ||T_j||_2 >= ||S_1||_2 ||S_rho||_2.

Es decir: descomponer solo REDISTRIBUYE la masa L^2, nunca la reduce
por debajo de ||S_rho||_2. La estimacion por piezas es siempre >= la
de Cauchy sin descomponer (doc 67, ruta (ii)), que ya es N/log N, un
log por encima del objetivo.

Medicion (proxy rho_tilde, script): sum S_{2,3}^2 = 10,476 = 0.070*pi;
sum S_{2,5}^2 = 58,091 = 0.390*pi; sum rho_tilde^2 = 218,454 =
1.467*pi (doc 67). Las piezas son DELGADAS (no gruesas: son conteos
de n con dos factores primos en rangos estrechos). Normas:
||S_1|| = 385.9, ||S_{2,3}|| = 102.4, ||S_{2,5}|| = 241.0, suma =
729.3 >= ||S_rho_tilde|| = 467.4 (factor 1.56). Cauchy por piezas:
int|S_1|^2 <= pi = 148,933 + C_{2,3} = 39,500 + C_{2,5} = 93,014, total
281,447 = 29.6x el objetivo N/log^2 N = 9,501; peor que el Cauchy del
doc 67 (180,375 = 19.0x). Y el termino int_m |S_1|^2 solo ya es
148,933 = 15.7x el objetivo: la descomposicion introduce el termino
|S_1|^2, que por si solo mata la ruta (su masa L^2 es pi(N) ~ N/log N;
los arcos mayores capturan solo o(1) de ella, por el calculo L^2
estandar de arcos mayores).

## 4. Prop. 69.2: la barrera sup x L^1 [T]

Por Cauchy-Schwarz y Parseval, int_m |S_1| <= (int|S_1|^2)^{1/2} =
pi(N)^{1/2} ~ (N/log N)^{1/2}. Luego la ruta sup x L^1 sobre una
pieza T_j exigiria sup_m |T_j| = o(N^{1/2}/(log N)^{3/2}).

(a) Lo disponible es la cota de Vinogradov sup_m |S| <<_A N/(log N)^A
(A fijo). El cociente con lo requerido es N^{1/2}/(log N)^{A-3/2} ->
infinito cuando N -> infinito: asintoticamente, la cota disponible
esta un factor ~N^{1/2} por encima de lo requerido. Con A = 3
(realista a N accesible): cociente 6.9x10^3 a N = 10^12 y 2.4x10^9 a
N = 10^24 (script). (Nota: con A = 10 la forma cruda N/(log N)^A
parece menor que lo requerido a N = 10^12-10^24, pero eso es un
artefacto: la forma asintotica solo vale para N >> 10^39 con ese A;
la constante implicita C(A) lo anula. Se usa A = 3.)

(b) Peor aun: ni la cancelacion de raiz cuadrada IDEAL basta. Si
sup_m |T_j| ~ N^{1/2} (lo mejor concebible para una suma exponencial
generica), sup x L^1 ~ N^{1/2} * (N/log N)^{1/2} = N/(log N)^{1/2},
que excede el objetivo N/log^2 N por (log N)^{3/2}: 145x a N = 10^12
y 411x a N = 10^24 (script). La ruta sup x L^1 falla DOS veces: las
cotas disponibles estan N^{1/2} por encima, y aun la cota ideal esta
(log N)^{3/2} por encima. El problema es int|S_1| ~ (N/log N)^{1/2}:
es demasiado grande para cualquier cota sup sobre la otra pieza.

## 5. Prop. 69.3: la diagonal del cuarto momento mixto [T]

La ruta (iii) del doc 67 necesitaria int_m |S_1|^2 |S_rho|^2 =
o(N^2/log^4 N). Para las piezas NO NEGATIVAS S_{2,3}, S_{2,5} (son
conteos >= 0), la expansion de Parseval tiene todos los terminos >= 0
y la diagonal p_1 = p_2, t_1 = t_2 es una cota INFERIOR rigurosa:

  int_0^1 |S_1|^2 |S_{2,3}|^2 >= pi(N) * sum_n S_{2,3}(n)^2.

Medicion: diagonal S_{2,3} = 1.56x10^9 y diagonal S_{2,5} = 8.65x10^9
frente a lo necesario o(N^2/log^4 N) = 9.03x10^7 a X = 2x10^6: 17.3x
y 95.8x por encima (script). Asintoticamente, la diagonal es
pi(N) * (c*N/log N) = c*N^2/log^2 N, un factor c*(log N)^2 por encima
de lo necesario, con c = 0.070 (S_{2,3}) y 0.390 (S_{2,5}). La ruta
del cuarto momento muere en la diagonal, antes de cualquier
estimacion no trivial. Para las piezas firmadas (S^triv, S^G) la
diagonal no es cota inferior (cancelacion entre terminos), pero sus
normas L^2 son >= las de estas piezas, asi que no estan mejor.

## 6. Veredicto y que sigue (etapa 2)

**Veredicto etapa 1: la descomposicion de Vaughan/Heath-Brown de
S_rho NO da la cancelacion de raiz cuadrada que falta, y no puede
darla por ninguna de las tres estimaciones por piezas.** Cauchy por
piezas no mejora el Cauchy sin descomponer (Prop. 69.1, medido:
29.6x el objetivo, peor que 19.0x); sup x L^1 falla dos veces
(Prop. 69.2: cotas disponibles N^{1/2} por encima; ideal (log N)^{3/2}
por encima); cuarto momento muere en la diagonal (Prop. 69.3: 17x y
96x). La razon estructural: todas las estimaciones por piezas
factorizan por normas SEPARADAS de S_1 y de las piezas, y esas normas
estan fijadas por la diagonal (el conteo de terminos). Descomponer
no crea cancelacion conjunta: S_1 y S_rho estan correlacionadas
(rho_tilde = 1_p - D, con D >= 0 delgado), y ninguna estimacion por
piezas usa esa correlacion.

**Etapa 2 (no hecha aqui):** una estimacion CONJUNTA que no factorice
por normas separadas: o bien una cota del momento mixto que use la
correlacion rho_tilde = 1_p - D (el termino cruzado int_m
S_1(alfa) * conj(S_D(alfa)) * e(-N alfa) con la fase conservada, no
en valor absoluto), o bien una estimacion de gran criba / valor medio
para el PRODUCTO S_1*S_rho directamente, o bien input de densidad
de ceros / GRH. Si ninguna de esas vias existe, el lema MM-69.1 queda
como cuello de botella abierto de la linea rho (docs 66-67). Goldbach
sigue abierta.

Referencias (leidas, no de memoria):
- Docs 66 (rho, c_rho, regiones S_{2,3}/S_{2,5}), 67 (conteo
  unilateral, Prop. 67.1/67.3, mediciones L^2), 64-65 (obstrucciones
  de dispersion y lema DUA-esquina).
- E. C. Titchmarsh / R. C. Vaughan, "The Hardy-Littlewood Method",
  2nd ed., caps. 2-3: identidad de Vaughan, cota de Vinogradov
  sup_m |S_Lambda| << N/(log N)^A, metodo del circulo binario.
- J. Maynard, arXiv:2006.08250v1, sec. 9 (descomposicion de Buchstab
  de rho) y Lema 6.1 (identidad de Heath-Brown para Lambda).
