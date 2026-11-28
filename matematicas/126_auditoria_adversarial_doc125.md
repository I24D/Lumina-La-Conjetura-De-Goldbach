# Doc 126 — Auditoría adversarial independiente del doc 125

**Fecha:** 2026-11-28  
**Agente:** LUMINA Chat  
**Estado:** revisión matemática y de reproducibilidad; no es una demostración de Goldbach.  
**Base auditada:** commit 4f07d202593fd6770e847c7366f3d9fd0c86c890 del repositorio I24D/Lumina-La-Conjetura-De-Goldbach.  
**Fuentes:** matematicas/125_inventario_h2.md, matematicas/125a_prime_window.md, matematicas/125b_erem_pos.md; scripts numerico/125a_prime_window.py y numerico/125b_erem_pos.py.

---

## 0. Veredicto ejecutivo

El núcleo matemático del doc 125 sobrevive con correcciones:

1. **(PRIME-WINDOW) queda refutada en el régimen u_L → 2**. La doble suma tiene límite positivo
   \[
   F(2)=\log(3/2)\log(4/3)=0.1166450426\ldots,
   \]
   mientras Q_c(z) ~ c_Q/log z. Por tanto PW(L)/Q_c(z) → ∞. Como SM(L) ≥ PW(L), también cae (OUT-MASS) en escala absoluta C Q_c.
2. **La refutación debe restringirse a 3/2 < u_L < 9/4**. La extensión impresa a todo u_L en (3/2,3) es falsa: para u_L ≥ 9/4 la ventana exterior es vacía (en 9/4, los extremos coinciden).
3. **Los lemas de paridad y los umbrales 343 y 1331 de 125b son correctos** bajo las definiciones de V+, V^(r), Vhat^(r) y Rem de la serie.
4. **La anulación uniforme E_Rem = 0 es falsa**, y existe un contraejemplo exacto que no requiere confiar en la enumeración exhaustiva del script.
5. **El cierre analítico del estrato 5 vía (POS) es correcto**. Sin embargo, el chequeo K3 del script 125b tiene un defecto que vuelve vacía su enumeración de candidatos; su PASS en esa subprueba es vacuo.

La ruta con signo (NET-LAYER)+(P-MASS) sigue siendo la ruta alta no refutada, pero el doc 125 no demuestra esos dos residuos. Goldbach permanece abierta.

## 1. 125a — identidad de ventana y refutación en u_L → 2

### 1.1. Parte exacta confirmada

Para p0 primo en la capa, p0^3 < L y u_p0 ≥ 1. De ahí:
- p0 < L^(1/3);
- L/p0 ≥ z;
- q^3 ≥ L implica q ≥ L^(1/3) > p0;
- q < z ≤ L/p0 implica p0 q < L.

Por tanto los tres filtros de la definición de Q(p0) se reducen exactamente a L^(1/3) ≤ q < z, con q primo permitido, y
\[
M_{near}^{(1)}(p_0)=S(z)-S(L^{1/3}).
\]
La suma es independiente de p0 en el nivel. Esta identidad y la reducción de los extremos son correctas.

Si u_L = log L/log z permanece en un compacto de (3/2,3), Mertens da
\[
S(z)-S(L^{1/3})
=\log(3/u_L)+O(1/\log z).
\]
El término principal interior es correcto.

### 1.2. Corrección de dominio: el intervalo exterior se vacía en 9/4

La condición u_p0 en [1,3/2), junto con p0^3 < L, da
\[
Lz^{-3/2}<p_0\le L^{1/3}
\]
cuando u_L > 3/2. Pero este intervalo solo es no vacío si
\[
Lz^{-3/2}<L^{1/3}
\iff u_L-3/2<u_L/3
\iff u_L<9/4.
\]
En consecuencia:
- para 3/2 < u_L < 9/4, la suma exterior tiene el término principal positivo
  \[
  \log\!\left(\frac{u_L/3}{u_L-3/2}\right);
  \]
- para u_L = 9/4, el intervalo (x,x] es vacío;
- para 9/4 < u_L < 3, el extremo inferior supera al superior y no hay primos de capa con esas dos condiciones.

Así, la fórmula impresa para PW(L) en todo (3/2,3) no puede ser correcta: para u_L > 9/4 su segundo logaritmo es negativo, aunque PW(L) ≥ 0. El enunciado correcto es, para u_L → u en (3/2,9/4),
\[
PW(L)=
\log(3/u)\,
\log\!\left(\frac{u/3}{u-3/2}\right)+o(1).
\]

**La conclusión central en u=2 no cambia.** En los marcos D=z^5, se puede escoger un primo admisible s en (z^3,2z^3) y L=floor(z^5/s); entonces u_L=2+O(1/log z). Bertrand garantiza tales primos para z grande. En esa sucesión,
\[
PW(L)\longrightarrow \log(3/2)\log(4/3)
=0.1166450426\ldots>0.
\]
Como Q_c(z) ~ c_Q/log z, se obtiene
\[
\frac{PW(L)}{Q_c(z)}
\sim \frac{0.1166450426\ldots}{c_Q}\log z
\longrightarrow\infty.
\]
Con c_Q aproximadamente 1.72, el coeficiente de log z es aproximadamente 0.0678. Puesto que SM ≥ PW término a término, (OUT-MASS) con constante absoluta contra Q_c queda refutada en este régimen. Esta es una refutación asintótica real, no una extrapolación de la ventana numérica.

Los filtros Pi(p0)<1/2 y no degeneración desaparecen para z suficientemente grande en esta sucesión: el extremo inferior hace p0 → ∞, Pi(p0) ~ 1/log z, y L/p0 ≥ L^(2/3) → ∞, de modo que existe un primo permitido q distinto de p0 con q^3 < L/p0.

### 1.3. Error de constante en Mertens

El doc define
\[
S(x)=\sum_{5\le p<x}\frac1{p-1}
\]
y escribe
\[
M'=B_1-\frac32+\sum_{p\ge5}\frac1{p(p-1)},
\]
donde B_1 es la constante de Mertens de sum_{p<x} 1/p.

Pero
\[
\frac1{p-1}=\frac1p+\frac1{p(p-1)}
\]
y quitar p=2,3 de sum 1/p resta 1/2+1/3=5/6, no 1+1/2=3/2. La constante correcta es
\[
\boxed{M'_{corr}=B_1-\frac56+
\sum_{p\ge5}\frac1{p(p-1)}.}
\]
La constante impresa es menor por 2/3. Por tanto, el E(x) definido con la constante impresa no tiende a cero: tiende a 2/3. El error no destruye la fórmula de la ventana porque en E(z)-E(L^(1/3)) la constante aditiva se cancela; con la constante corregida, E(x)=O(1/log x) y la derivación queda limpia. Debe corregirse el normalizador y la afirmación E(x) → 0.

### 1.4. Signo de la corrección por filtros

En la suma exterior, el doc define epsilon_Pi ≥ 0 como la masa de primos del intervalo completo excluidos por Pi(p0) ≥ 1/2 o degeneración, pero escribe
\[
\sum_{p_0\in\mathrm{capa}_1}\frac1{p_0-1}
=S(L^{1/3})-S(Lz^{-3/2})+\varepsilon_\Pi.
\]
Si epsilon_Pi es la masa excluida y no negativa, el signo correcto es **menos**:
\[
\text{suma filtrada}=
\text{suma del intervalo completo}-\varepsilon_\Pi.
\]
El signo impreso es un error algebraico. No afecta el límite en u → 2, porque los filtros se vuelven vacíos para z suficientemente grande y entonces epsilon_Pi=0.

## 2. 125b — auditoría de E_Rem

### 2.1. Lemas de paridad y umbrales: confirmados

La cadena de la prueba es válida:

- Si r es visible, existe d en V+(L) con r|d. Al ser ambos squarefree, los primeros i primos de r son término a término no mayores que los de d; sus productos de prefijo también. Las desigualdades de Nathanson para d implican las de r. Por tanto r pertenece a V+(L).
- Si d' pertenece a Rem(r), sea m una posición impar que falla en la lista fusionada y j el número de primos de d' en ese prefijo. j=0 contradice la visibilidad de r. Si j es impar positivo, la condición de V+(L/r) en la posición j da r q_j^2 Q_j < L, mientras que el fallo fusionado da L ≤ p_m^2 R* Q_j ≤ q_j^2 r Q_j, contradicción. Luego j es par y, al no poder ser cero, j≥2. Así omega(d')≥2.
- Si omega(d')≥2, el mayor primo q1≥7, y la condición en la posición 1 da q1^3<L/r. Por tanto L/r>343.
- Si L/r≤1331, no puede haber omega(d')≥3, pues entonces q1≥11 y q1^3<L/r exigiría L/r>1331. Todo elemento de Rem tiene exactamente dos primos; por tanto mu(d')=+1 y
  \[
  E_{Rem}(r)=-\sum_{d'\in Rem(r)}\frac1{\varphi(d')}\le0.
  \]

Se confirman, pues:
- Rem(r) vacío y E_Rem=0 para L/r≤343;
- E_Rem≤0 para L/r≤1331;
- el residuo de signo para L/r>1331 sigue abierto.

### 2.2. Contraejemplo exacto reforzado: E_Rem=0 uniforme es falso

Se puede refutar sin depender de la afirmación computacional más fuerte Rem(r)={35}.

Tómese z=61, D=61^5=844596301, s=215077 (primo admisible) y
\[
L=\lfloor D/s\rfloor=3926,\qquad r=11.
\]
Entonces r es visible porque 11^3=1331<3926, y
\[
u_r=\frac{\log(3926/11)}{\log61}=1.42974\ldots\in[1,3/2).
\]
Con el producto Pi(r) del marco, Pi(11)=151963425/385875968 aproximadamente 0.393814<1/2; además 7 es un testigo no degenerado porque 11·7^3=3773<3926.

Ahora d'=35 es coprimo con r, r d'=385<3926, y en el modelo V+(L/r) la única condición impar es
\[
r\cdot7^3=3773<3926.
\]
Por tanto 35 pertenece a Vhat^(r). En la lista fusionada de r d'=385, los primos son 11,7,5; la condición en la posición m=3 falla:
\[
5^2\cdot(11\cdot7\cdot5)=9625\ge3926.
\]
Así 35 pertenece a Rem(r). Como L/r=356.909...≤1331, el Lema 4 fuerza que **todos** los elementos de Rem tengan dos primos y contribuyan con el mismo signo. En particular,
\[
E_{Rem}(11)\le-\frac1{\varphi(35)}
=-\frac1{24}<0.
\]
Esto prueba exactamente que la anulación uniforme es falsa. La afirmación más fuerte Rem={35} puede conservarse como resultado de enumeración exhaustiva, pero no es necesaria para el contraejemplo.

## 3. (POS), estrato 5: prueba válida; K3 del script no valida su barrido

La derivación analítica de c en {1,5,7} es correcta. Los casos k impar, k par con 5|c, p_k=7 y p_k≥11 quedan cubiertos por las cotas
\[
c<125/p_k^2\quad(k\text{ impar}),
\qquad
c<625/49,\quad
c<875/121,\quad
c<1375/169,
\]
según la rama correspondiente. La exclusión de c=11 usa p_k=11, p_{k-1}≥13, y la cota 1375/169<11. Para p_k≥13, la desigualdad auxiliar puede justificarse directamente:
\[
11p^2-125p+44=(p-11)(11p-4)>0
\qquad(p\ge13).
\]
No depende de un barrido numérico hasta 10^5.

**Defecto encontrado en numerico/125b_erem_pos.py, función check_K3, subfunción sqfree_cs:**

La lista pr=primes_upto(100000) comienza por 2,3,5,... Dentro de rec(start,cur) aparece:

    if p < start or cur * p >= B:
        break

La llamada inicial es rec(5,1). Al encontrar p=2<5, el bucle ejecuta break y termina inmediatamente. Por tanto sqfree_cs(125) devuelve la lista vacía; el chequeo de candidatos alive queda vacío por construcción y pasa de forma vacua. El PASS global no valida ese barrido.

**Reparación:** sustituir la primera condición por if p < start: continue, conservando break para cur*p >= B, o iterar sobre una lista de primos filtrada desde 5. Con esa corrección, la enumeración genera 37 candidatos squarefree c<125 con factores primos ≥5, y el filtro de las ramas deja vacío el conjunto de supervivientes fuera de {1,5,7}. En todo caso, la prueba analítica anterior ya demuestra el resultado.

## 4. Estado corregido de la línea (H2)

**Confirmado [T]:**
- Identidad exacta de ventana M_near^(1)=S(z)-S(L^(1/3)).
- Refutación asintótica de (PRIME-WINDOW) y (OUT-MASS) absoluta en una sucesión admisible con u_L→2.
- Lema de paridad; E_Rem=0 para L/r≤343; E_Rem≤0 para L/r≤1331.
- (POS) y cierre del estrato 5: c en {1,5,7}.

**Correcciones de escritura necesarias:**
- Restringir la fórmula de PW a u_L<9/4.
- Corregir M' por +2/3.
- Cambiar +epsilon_Pi por -epsilon_Pi si se mantiene la definición de masa excluida no negativa.
- Reparar sqfree_cs en K3 y volver a ejecutar el script.

**Pendiente [P]:**
- (NET-LAYER) y (P-MASS), la cancelación con signo M↔P.
- (EREM-SIGN) para L/r>1331, posible cuando z>121.
- (FC-2), (BAND-LOW-NONDEG) y (SLAB-V') fuera del estrato 5.

La refutación de la ruta absoluta es sólida en u_L→2; no prueba la ruta con signo. Ninguna parte de este doc demuestra la conjetura fuerte de Goldbach.

## 5. Alcance de esta auditoría

Se inspeccionaron los tres documentos 125 y el código de ambos scripts en el commit indicado. Se hicieron comprobaciones algebraicas independientes de las identidades y umbrales, del contraejemplo r=11,L=3926, y del generador de candidatos de K3. **No se afirma haber ejecutado de extremo a extremo las suites históricas con todas sus dependencias importadas**; el defecto de K3 se detectó por inspección directa del código y la prueba analítica de (POS) se verificó por separado.

**Conclusión:** doc 125 contiene un avance matemático real. La refutación de (OUT-MASS) absoluta y el rescate del signo de E_Rem sobreviven. Deben corregirse el dominio de u_L, la constante de Mertens, el signo de epsilon_Pi y el chequeo K3. Goldbach sigue abierta.
