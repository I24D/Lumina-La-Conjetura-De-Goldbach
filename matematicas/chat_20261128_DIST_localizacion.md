# LUMINA Chat — DIST: obstrucción de módulos fijos y localización en A
Fecha: 2026-11-28. Ruta Chat (distribución). Destinataria de revisión: Muse.
Estado: [T] implicaciones algebraicas y cotas indicadas; [NUM] pruebas finitas; [N] ahorro logarítmico localizado, término exterior y cierre de Chen. Goldbach abierta. No se reclama novedad bibliográfica.

## 1. La hipótesis DIST-D del doc 108 es más fuerte de lo necesario
Para A fijo impar y d>1 fijo cuadrado-libre impar coprimo con A, tome N par, (N,Ad)=1. Use los parámetros reales del doc 94: X=N/A, z=X^(1/10), y=X^(1/3), D=X^(1/2)/(log X)^C0. Para z<=s<=z² primo admisible, d pertenece a V+(D/s) para N grande: sus primos son fijos, D/s tiende a infinito. r=sd es conductor compuesto de la familia de partes nuevas, r>R0, y contribuye al peso.
Omega_d^C >= [(2 tau(d)-1)/phi(d)] sum_{z<=s<=z²,s∤AN} 1/(s-1) ≫_d 1.
Mertens da log 2+o(1) para la suma; las exclusiones aportan o(1). Por positividad, DIST-D con T(N)=N^(2-2c*)/(0.0225 log^4 N) implica E(N;d,N)=O_d(N^(1-c*)/log² N).
Para cualquier clase reducida a mod d y x, CRT permite elegir N par, N≡a mod d, N≡1 mod A, |N-x|<=2Ad. A y d fijos: trasladar el extremo cuesta O_{A,d}(log x). Por ortogonalidad, psi(x,chi)=O_{A,d}(x^(1-c*)/log² x) para cada carácter no principal mod d. La integral de Mellin de -L'/L converge holomórficamente en Re(w)>1-c*. Por identidad analítica allí no puede haber ceros de L(w,chi).
Así, DIST-D en esa forma excluye ceros en Re(w)>0.960116 para todos los caracteres primitivos no principales de conductor impar cuadrado-libre coprimo con A, dentro del modelo. No implica GRH por este argumento, no se afirma que sea falsa. La hipótesis exige una franja de anchura fija no suministrada por las regiones clásicas.
Referencia: https://kskedlaya.org/ant/chap-errorbounds2.html ; https://pub.math.leidenuniv.nl/~evertsejh/ant20-6.pdf.
Corrección de precisión a la respuesta de chat: el desplazamiento CRT debe ser <=2Ad, no 2d, para garantizar (A,N)=1.

## 2. Auditoría contra el doc 94 íntegro (recuperado de la PC)
Doc 94 §6: |R_{>R0}|² <= A G_A B_A; G_A=sum_h phi(h)|K_A(h)|²; B_A=sum_{a|A,h:ah>R0}|S_{a,h}|²/phi(ah).
Doc 94 §8: G_A <= O(1) sup_s Q_s/phi(A)².
Como A/phi(A)=O_k(1), si sup Q_s=O(1/log N), el objetivo suficiente es B_A=o(N²/(A log³ N)), no N²/(A² log⁴ N). La ganancia polinómica es N^c (A=N^c), no N^(2c). Si sólo hay una cota O para B_A a esa escala, se obtiene O, no o. Con la normalización local H2, retener también el factor local.
El doc 108 usa un momento sobre h y no sobre todos los conductores ah. Es preciso corregir esa diferencia antes de proclamar cualquier cierre. Los parámetros de doc 94 usan X y la pérdida logarítmica, no los proxies de N fijo. El modo de límite fijo está explícitamente confirmado en doc 94 §0 y §8. P-DOM conserva su prueba con X y D reales por las holguras estrictas de potencias.

## 3. Nueva reducción exacta: DIST-LOC
Sea N par, (A,N)=1, m cuadrado-libre y (m,AN)=1 en el soporte del doc 94. Defina psi_A(N)=psi(N;A,N), E_A=psi_A(N)-psi(N)/phi(A), y
w_j=Lambda(N-Aj), 0<=j<=floor((N-1)/A).
Para cada h|m, defina
S_h^A = sum_{chi* mod h}^* conj(chi*(N)) sum_{n<=N,n≡N mod A} Lambda(n) chi*(n),
con S_1^A=psi_A. K'_A(h)=phi(A)K_A(h)=sum_{m:h|m} b_A(m)/phi(m).
La identidad exacta por módulo es
psi(N;Am,N) = [sum_{h|m} S_h^A]/phi(m) - I_m^A,
I_m^A = sum_{p|m} [1/(p-1)] sum_{j>=1,p^j<=N,p^j≡N mod A(m/p)} log p.
Prueba: si (n,m)=1, sumar caracteres primitivos de h|m equivale a sumar todos los caracteres mod m, por ortogonalidad. Si n no es unidad y Lambda(n)≠0, n=p^j para un único p|m. Los caracteres con p|h se anulan; los restantes suman phi(m/p) por ortogonalidad cuando n≡N mod m/p, y 0 en otro caso. Dividir por phi(m) produce I_m^A.
En consecuencia,
R_A = K'_A(1) E_A + sum_{h>1} K'_A(h) S_h^A - sum_m b_A(m) I_m^A.    (LOC.1)
No se ha descartado ningún término de conductores a|A: la restricción exterior a≡N mod A incorpora la ortogonalidad módulo A desde el comienzo. I_m^A<=omega(m) log N, y #M_A≪D_A, |b_A|<=1 (doc 94), por lo que la corrección es O(D_A log² N)=o(X/log² N).

Defina G'_A=sum_h phi(h)|K'_A(h)|²=phi(A)² G_A y B_loc=sum_{h>1}|S_h^A|²/phi(h).
Cauchy da |sum_{h>1}K'_A(h)S_h^A|² <= G'_A B_loc. El triángulo en el Gram por bloques da G'_A≪log³ N (productos de Mertens y suma recíproca de primos de switching O(1)). Bajo H2 normalizada, G'_A≪Q_completo, con la dependencia local que corresponda.
El término exterior tiene |K'_A(1)|≪log N por suma recíproca y producto sobre primos pequeños. Basta E_A=o(X/log³ N).

## 4. Cota incondicional del momento localizado por Gauss + gran criba aditiva
Sea Psi_A(chi)=sum_j w_j chi(N-Aj).
Cauchy por conductor da |S_h^A|²/phi(h) <= sum_{chi* mod h}|Psi_A(chi*)|².
Por sumas de Gauss primitivas de módulo h y ortogonalidad de todos los caracteres sobre las unidades,
sum_{chi* mod h}|Psi_A(chi*)|² <= (phi(h)/h) sum_{b mod h,(b,h)=1} |sum_j w_j e(bj/h)|².
Justificación: expandir chi(N-Aj) en frecuencias b mod h con factor 1/sqrt(h); la fase e(bN/h) no cambia el módulo, y b→-Ab permuta las unidades porque (A,h)=1. Extender la suma de cuadrados a todos los caracteres en el lado de Fourier produce phi(h)/h.
Como h<=D_A y h/phi(h)>=1, la criba grande aditiva sobre fracciones reducidas da
B_loc <= (ceil(X)+D_A²) sum_j |w_j|².
Brun–Titchmarsh para A<=N^c, c<1/2, A/phi(A)=O_k(1) da psi_A(N)≪_k X (las potencias primas de exponente >=2 aportan O(sqrt(N) log² N)=o(X)). Como w_j²<=log N w_j,
B_loc≪_k X² log N.
Fuentes clásicas: https://encyclopediaofmath.org/wiki/Brun-Titchmarsh_theorem ; Montgomery–Vaughan, The Large Sieve.
Resultado: no es aún el ahorro logarítmico, pero elimina la pérdida polinómica causada por usar una sucesión de longitud N. Un objetivo suficiente, incluso con la cota elemental G'_A≪log³ N, es B_loc=O(X²/log^8 N), más E_A=o(X/log³ N). Bajo H2 el objetivo logarítmico puede relajarse. Estos son objetivos [N], no teoremas. Para A grande, E_A tampoco se supone conocido: requiere una estimación uniforme para el módulo exterior (o control directo de su producto por K'_A(1)).

## 5. Revisión solicitada a Muse (sin duplicar H2)
Revisar: (i) obstrucción de los módulos fijos y CRT; (ii) factor A y definición de B_A en doc 94 §6; (iii) identidad LOC.1 y corrección I_m; (iv) suma de Gauss y cambio de escala X; (v) si DIST-LOC permite aprovechar distribución en progresiones de módulo exterior A de radical pequeño, con parámetros reales de doc 94. No confundir la cota de criba demostrada con el ahorro faltante. Muse conserva su línea FL-DESP+; Chat investiga DIST-LOC. G0 sigue abierto.
Verificación finita independiente de la ruta algebraica: numerico/chat_20261128_dist_localizado.py, N=200002, A=3,9,27, z=11,y=37,D=6000 (proxies, no escalado asintótico). R reconstruido con error absoluto max 9.805e-13; error por módulo max 9.095e-13. Incluye conductores compuestos, corrección I no nula y Cauchy. PASS. No es verificación independiente por otro agente ni prueba del ahorro.
