# 67. Conteo unilateral con rho de Harman (SL-65.2, etapa 2)

Fecha: 2026-11-26. Etapa 2 del proyecto SL-65.2 (doc 66): rehacer criba /
metodo del circulo con el minorante rho de Harman para el conteo unilateral
U(N) = sum_{p<N} rho(N-p). Solo etapa 2: el conteo calculado, el termino
principal y el error, y el veredicto. Script:
`numerico/67_conteo_unilateral_rho.py` — PASS 46/46, verificado local y en
el PC con numeros identicos. No se afirma ninguna prueba de Goldbach.

## 1. El conteo y el termino principal

Se trabaja con el proxy sin parte B,
rho_tilde(n) = 1_p(n) - S_{2,3}(n) - S_{2,5}(n), con las regiones del doc 66
(Maynard III sec. 9). El rho verdadero difiere de rho_tilde solo por la
parte B, de masa relativa O(delta) (doc 66: fraccion 0.0431 de n <= X con
factor en B); todo lo medido aqui se transfiere a rho con un factor
(1 + O(delta)). rho es diadico (n ~ x); el conteo a escala N usa rho a
escala N (diseccion diadica, como en el paper).

Prop. 67.1 (termino principal): el conteo unilateral tiene termino
principal c_rho * S(N) * N/log^2 N, con c_rho en [0.43, 0.47] (doc 66) y
S(N) = 2*C_2*prod_{p|N, p>2} (p-1)/(p-2) la serie singular usual de
Goldbach. Es decir: la serie singular es LA MISMA que para R(N); solo
cambia la constante (c_rho en vez de 1). Prueba: el lema de soporte
(abajo) implica que el defecto D(n) = S_{2,3}(n) + S_{2,5}(n) se anula
cuando P^-(n) < z1, luego la masa de rho en la clase 0 mod p es 0 (para
rho_tilde, exactamente 1: solo n = p) y en cada clase reducida es
c_rho/(p-1) de la masa total; los factores locales son por tanto los de
1_p y el producto da S(N). El error de arcos mayores es << N/(log N)^A:
el lado primo usa Bombieri-Vinogradov (nivel 1/2) y el lado rho usa el
nivel 1/2+delta del Teor. 1.3(3) en la ventana Q1 en [x^{2/5+5delta},
x^{3/7}], con las componentes S^triv, S^G, rho_B equidistribuyendo en
modulos pequenos (sec. 9, Prop. 5.1, Lema 6.8). [T-lit para el error de
arcos mayores; la parte de la serie singular esta verificada
numericamente para rho_tilde en el script.]

Lema de soporte (exacto, script): si P^-(n) < z1 = x^{1/7}, entonces
S_{2,3}(n) = S_{2,5}(n) = 0. Prueba: p1, p2 >= z1 no pueden contener el
factor primo pequeno de n, que sobrevive en el cofactor n/(p1*p2), y
rho(cofactor, p2) = 0 cuando P^-(cofactor) < p2. Verificado en los
1,542,857 n <= 2x10^6 con P^-(n) < z1. Corolario: rho_tilde(n) = 1_p(n)
en todos esos n.

Verificaciones (script, X = 2x10^6, pi(X) = 148933):
- rho_tilde(p) = 1 en todo primo; rho_tilde(n) <= 0 en todo compuesto.
  max D(n) = 3; 37,832 n con D >= 1.
- L_S = (S_a + S_b)/pi = 0.3249 (doc 66: 0.325); S_a/pi = 0.0663,
  S_b/pi = 0.2586.
- Densidad en clases reducidas: p = 3, 5, 7 (p < z1): clase 0 = 1
  exacto; desviacion maxima en clases reducidas 0.001 (p=3), 0.000
  (p=5, 7). Para p = 11, 13 (p > z1, sin lema): clase 0 = -5339, -6081
  (el defecto SI vive ahi, como debe ser); desviacion reducida 0.006,
  0.005. La serie singular es la estandar, verificada en las clases.
- sum rho_tilde = 100,544 = (1 - L_S)*pi(X), exacto.

## 2. Medicion finita del conteo unilateral

U_tilde(N) = sum_{p<N} rho_tilde(N-p), R(N) = #{p<N: N-p primo},
r(N) = U_tilde(N)/(S(N)*N/log^2 N):

| N       | U_tilde | R(N)  | S(N)  | U/R    | 1-L_S(N) | 1-L_S*kappa | kappa  | r/(1-L_S) |
|---------|---------|-------|-------|--------|----------|-------------|--------|-----------|
| 200000  | 2293    | 2834  | 1.760 | 0.8091 | 0.8187   | 0.8058      | 1.0714 | 1.1851    |
| 500000  | 4595    | 6104  | 1.760 | 0.7528 | 0.7770   | 0.7614      | 1.0696 | 1.1570    |
| 1000000 | 7548    | 10804 | 1.760 | 0.6986 | 0.7267   | 0.7092      | 1.0640 | 1.1262    |
| 1500000 | 20285   | 30328 | 3.521 | 0.6689 | 0.6961   | 0.6776      | 1.0609 | 1.1158    |
| 2000000 | 12596   | 19440 | 1.760 | 0.6479 | 0.6751   | 0.6559      | 1.0590 | 1.1155    |

Lectura: U_tilde(N) <= R(N) siempre (termino a termino, rho_tilde <=
1_p) y U_tilde(N) > 0 en todos los N medidos. U/R sigue a 1 - L_S(N)
con un deficit sistematico del ~4%, explicado casi todo por
kappa(N) = <pr/phi(pr)> pesado por los pares (1.0714 -> 1.0590,
decreciendo hacia 1: los primos de los pares crecen como N^{1/7});
con kappa, U/R = 1 - L_S*kappa dentro del 1%. El cociente
r/(1-L_S) decrece hacia 1 (1.185 -> 1.116): r(N) tiende a 1 - L_S(N),
que a su vez tiende a c_rho ~ 0.465. El termino principal de la
Prop. 67.1 es correcto; la discrepancia restante es efecto de escala
finita (R(N) mismo esta ~16% por encima de S(N)*N/log^2 N a estas
escalas: convergencia lenta O(1/log N)).

## 3. El error: los arcos menores NO estan controlados

Los arcos mayores funcionan (sec. 1). Los arcos menores son la
obstruccion. Sea S_1(alfa) = sum_{p<N} e(p*alfa),
S_rho(alfa) = sum_n rho(n) e(n*alfa). Hay que probar
int_m |S_1 * S_rho| = o(N/log^2 N). Las tres rutas disponibles fallan,
cada una por una razon estructural distinta:

(i) sup_m |S_1| * int_m |S_rho|: sup_m |S_1| << N/(log N)^A
(Vinogradov) pero int |S_rho| <= (sum rho^2)^{1/2} ~ (N/log N)^{1/2},
luego el producto es ~ N^{3/2}/(log N)^{A+1/2}: sobra un factor
N^{1/2}/(log N)^{A-3/2} -> infinito. La norma L^1 de S_rho es
demasiado grande (~N^{1/2}); el sup no la compensa.

(ii) Cauchy: int_m |S_1*S_rho| <= (int|S_1|^2)^{1/2} (int|S_rho|^2)^{1/2}
= (pi(N) * sum rho^2)^{1/2}. Medido: sum rho_tilde^2 = 218,454 =
1.467*pi(X); en N = 2x10^6, Cauchy = 180,375 frente a termino
principal 11,292: cociente 16.0, y asintoticamente ~2*log N -> infinito.
La cota de Cauchy excede el termino principal por un logaritmo: las
normas L^2 estan fijadas por la diagonal (el conteo de terminos, ~N/log N)
y su producto es N/log N, un log por encima de N/log^2 N. Ninguna cota
que use solo normas L^2 puede bajar de ahi.

(iii) Cuarto momento / momento mixto: haria falta
int |S_1|^2|S_rho|^2 = o(N^2/log^4 N). Lo disponible es mucho peor:
Li-Liu (doc 16) da int_m |S_w|^4 << X^{5/2} (log X)^{100} bajo
cuasi-RH Theta < 11/12, frente al X^2/log^4 X necesario: sobra X^{1/2}.
No existe ninguna estimacion de momento mixto para el par (1_p, rho).

Prop. 67.3: con los metodos conocidos, el error de arcos menores del
conteo unilateral NO se puede acotar por debajo del termino principal.
El minorante rho elimina la obstruccion de equidistribucion (arcos
mayores / nivel 1/2+delta, docs 64-66) pero NO toca la obstruccion de
arcos menores: esta es una obstruccion de sumas exponenciales
oscilantes, no de equidistribucion. Es la misma razon estructural por
la que Goldbach binario esta abierto: el metodo del circulo con dos
factores necesita cancelacion de raiz cuadrada o un cuarto momento, y
no hay ninguno disponible.

Lema que falta (delimitado, no probado): una cota
int_m |S_1(alfa) S_rho(alfa)| d alfa = o(N/log^2 N) que NO sea solo
L^2 (o cancelacion de raiz cuadrada para S_rho en arcos menores).
Probarlo requeriria una idea nueva para el par (1_p, rho) o input de
nivel GRH; refutarlo tampoco es posible con lo conocido. La ruta de
criba para U(N) es circular: rho YA ES la salida de la criba; una cota
inferior por criba necesitaria una criba mas alla de la barrera de
paridad.

## 4. Veredicto y que dice sobre P_2

Veredicto: el conteo unilateral con rho tiene el termino principal
correcto (c_rho * S(N) * N/log^2 N, Prop. 67.1, verificado
numericamente) y los arcos mayores controlados (nivel 1/2+delta),
pero los arcos menores NO estan controlados con ningun metodo
conocido (Prop. 67.3). NO da una cota inferior positiva para
Goldbach. La obstruccion esta delimitada con precision: falta el
lema del momento mixto / cancelacion en arcos menores. Goldbach
sigue abierta.

P_2: el conteo es todo o nada. rho_tilde(m) > 0 si y solo si m es
primo (en el rango diadico), luego U(N) > 0 implicaria R(N) > 0
directamente (Goldbach fuerte, no solo P_2), y no hay intermedio:
rho se anula (<= 0) en TODOS los compuestos, incluidos los P_2.
La ruta rho no dice nada nuevo sobre P_2 mas alla de lo que diria
una prueba completa de Goldbach. Un minorante para P_2 seria otra
construccion (no intentada aqui).

## 5. Que sigue

- SL-65.1 (bucle): agudeza del umbral 1/10 de la Prop. 5.1 / Teor. 1.1
  de Maynard III.
- El lema que falta (sec. 3): buscar una idea nueva para el momento
  mixto del par (1_p, rho) o para la cancelacion de S_rho en arcos
  menores. Es el cuello de botella preciso de la etapa 2; sin el, la
  linea rho no avanza hacia una cota inferior.
- La linea del bucle (seleccion: switching, C10/W(d); Motohashi doc 57
  sec. 4; b_{P,l} doc 53; uniones C_k+C_{k+1}; variante Lema 4
  Pintz-Althofer) no se toca aqui: esta etapa no la duplica.

Referencias (leidas, no de memoria):
- J. Maynard, "Primes in arithmetic progressions to large moduli III:
  Uniform residue classes", arXiv:2006.08250v1 (2020), sec. 9 (Teor.
  1.3), sec. 6 (Lemas 6.1, 6.5, 6.8), Prop. 5.1. Copia ORA.
- Doc 66 (etapa 1): formulacion de rho, c_rho en [0.43, 0.47],
  correccion del conteo unilateral (constante c_rho, no c_rho^2).
- Doc 64: obstrucciones de dispersion uniforme en a (Teor. 1.3
  bloqueado por paridad para Lambda). Doc 65: lema DUA-esquina
  cerrado. Doc 16: cota de cuarto momento Li-Liu.
