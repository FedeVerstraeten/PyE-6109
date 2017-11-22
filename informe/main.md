\pagebreak

*Objetivo: Lorem ipsum dolor sit amet, consectetur adipiscing elit. Fusce sagittis augue sit amet metus tristique lacinia. Curabitur sit amet tortor accumsan, hendrerit purus sit amet, volutpat dui. In molestie hendrerit elementum. Phasellus aliquam tortor nec risus auctor ullamcorper. Pellentesque et diam at metus elementum sagittis. Ut dapibus vulputate efficitur. Praesent nec tortor et mi gravida commodo a nec metus. Maecenas semper ut leo eu cursus.*

#1. Análisis y métodos de simulación

El método de Box-Müller, permite la obtención de dos variables aleatorias normales independientes $X_1$ y $X_2$ a partir de dos uniformes independientes $U_1$ y $U_2$. En este sentido, se definen las variables aleatorias:

$$  R := \sqrt[]{2 \textrm{ } log (U_1)}  $$
$$  \Theta := 2 \pi U_2 $$

Se puede observar que como $U_1$ y $U_2$ pertenecen a espacios de probabilidad independientes y las transformaciones que se les aplican son regulares, entonces $R$ y $\Theta$ también son independientes. Entonces,
$$  S := R^2  $$
$$  \implies S = 2 \textrm{ } log (U_1) $$
$$  \implies s = 2 \textrm{ } log (u_1) = g(u_1)$$
$$  \textrm{g es biyectiva} \implies g^{-1}(s) = e^{-\frac{s}{2}} $$

Por el método del jacobiano se obtiene que:

$$  \implies f_{S}(s) = \left.\frac{f_{U_1}(u_1)}{| J_g |}\right|_{u_1=g^{-1}(s)} $$
$$  \implies f_{S}(s) = \frac{1}{2} e^{-\frac{1}{2}s} \cdot 1 \{ s>0\} $$
$$  R = \sqrt[]{S} \quad,\quad  h(s)=\sqrt[]{s} \implies h^{-1}(r)= r^2$$
$$  \implies f_{R}(r) = \left.\frac{f_{S}(s)}{| J_h |}\right|_{s=h^{-1}(r)} = \cdots r e^{\frac{r^2}{2}} \cdot 1 \{ r>0\} $$
$$ \textrm{y además } f_\Theta (\theta) = \frac{1}{2 \pi} 1 \{ 0 < \theta < 2 \pi \}$$

$$\implies f_{R,\Theta}(r,\theta)= f_{R}(r) f_\Theta (\theta) = \cdots = \frac{r}{2 \pi} e^{-\frac{r^2}{2}} \cdot 1 \{r>0 \quad,\quad 0 < \theta < 2 \pi \}$$

Si ahora se aplica un nuevo cambio de variables :

$$ (Z_1,Z_2) := (R \cos (\Theta), R \sin (\Theta))$$
$$\implies | J_g | = r $$
$$\implies f_{Z_1,Z_2}(z_1,z_2)= \left.\frac{f_{R,\Theta}(r,\theta)}{| J_g |}\right|_{g^{-1}}$$
$$\implies f_{Z_1,Z_2}(z_1,z_2)= \frac{1}{2\pi} e^{-\frac{1}{2}(z_1^2 + z_2^2)}$$
$$\implies Z_1,Z_2 \sim \mathcal{N}(0 , 1) \textrm{ independientes}$$

Como es de interés simular dos normales estándar correlacionadas, se puede construir una matriz de covarianzas $A$ tal que
$$ A =\begin{bmatrix}
\sigma^2_1 & \frac{\rho \sigma^2_1 \sigma^2_2}{2}\\
\frac{\rho \sigma^2_1 \sigma^2_2}{2} & \sigma^2_2
\end{bmatrix} \quad , \quad \textrm{siendo $\rho$ el coeficiente de correlación}$$

Dado que $A$ es definida positiva y hermítica, se puede aplicar la descomposición de Cholesky para A y así obtener las variables aleatorias $X_1$ y $X_2$ normales estándar correlacionadas.


#2. Simulación de Variables Aleatorias

A partir de las distintas implementaciones provistas por la cátedra, para la simulación de normales estándar independientes mediante el métodode *Box-Muller* y normales estándar correlacionales mediante la descomposición de *Cholesky* de la matriz de covarianzas, se optó utilizar la implementación `simulacion_normal_bivariada_numpy.py`.

Se presentan en la siguiente tabla las primeras 3 (tres) simulaciones, donde cada columna representan:

**$U_{1,i}$,$U_{2,i}$:** par de V.A. uniformes.

**$Z_{1,i}$,$Z_{2,i}$:** par de V.A. normales estándar independientes.

**$X_{1,i}$,$X_{2,i}$:** par de V.A. normales correlacionadas.


| i |$U_{1,i}$|$U_{2,i}$|$Z_{1,i}$|$Z_{2,i}$|$X_{1,i}$|$X_{2,i}$|
|:-:|--------:|--------:|--------:|--------:|--------:|--------:|
| 1 | 0.065405| 0.996567|  2.33490| -0.05038|  2.33490|  2.23997|
| 2 | 0.426408| 0.986825|  1.30118| -0.10796|  1.30118|  1.22732|
| 3 | 0.216612| 0.722268| -0.30323| -1.72260| -0.30323|-0.744368|

Para las simulaciones se consideró un coeficiente de correlación $\rho = 0{.}965$, obteniendo así una matriz de Covarianza:

$$
Cov =
 \begin{pmatrix}
  \Var[Z_1] & \rho\sqrt{\Var[Z_1] \cdot Var(Z_2)}\\
  \rho\sqrt{\Var[Z_1] \cdot Var(Z_2)} & Var(Z_2) 
 \end{pmatrix}
=
 \begin{pmatrix}
  1 & 0{.}965 \\
  0{.}965 & 1
 \end{pmatrix}
$$

#3. Análisis de los gráficos

Por medio de la simulación se obtienen el gráfico de la distribución empírica de la variable $X_2$ para $n_{sim} = 10^4$.

![Función distribución empírica de $X_2$](img/fde.png){width=80%}


Función histograma de $X_2$ con los valores límite para $n_{sim} = 10^4$: 

$$[-a_0, -3, -2, -1, -0{.}5, 0, 0{.}5, 1, 2, 3, a_{10}]$$

Donde sean los coeficinetes,

$$a_0 = min(-4{.}5, min(X_1), min(X_2))$$
$$a_{10} = max(-4{.}5, max(X_1) + 10^{-8}, max(X_2) + 10^{-8})$$

![Función histograma de $X_2$](img/hist.png){width=80%}

#4. Tendencia de valores

En el caso evaluado, se considera el coeficiente de correlación como $\rho = \frac{100 + a}{200}$, siendo $a$ el valor de las dos últimas cifras del padrón de cada integrante. 
Para la tendencia, podemos observar en particular el coeficiente correspondiente a Lucía ($\rho \sim 0.5$) respecto al de Pablo ($\rho \sim 1$), que cuanto se acerca $\rho$ de 1, más se acercan los puntos del gráfico a formar una recta de pendiente positiva, coincidiendo con la recta de regresión.

Se observa la tendencia esperada por la recta de regresión para una distribución normal estándar bivariada.

En los siguientes gráficos, generados a partir de las simulaciones, podemos observar el comportamiento de la normal bivariada obtenida respecto a diferentes coeficientes de correlación, marcando la distribución de puntos $(X_{1,i} , X_{2,i})$ en el plano.

\pagebreak

**Estudiante:** Pablo González\
**Padrón:** 96993

![Recta regresión con a = 12](img/foo5_Pablo.png){width=80%}

**Estudiante:** Lucía Kasman\
**Padrón:** 97112

![Recta regresión con a = 12](img/foo5_Lucia.png){width=80%}

**Estudiante:** Leonardo Taffarel\
**Padrón:** 97066

![Recta regresión de $X_2$](img/foo5_Leo.png){width=80%}

**Estudiante:** Federico Verstraeten\
**Padrón:** 92892

![Recta regresión de $X_2$](img/foo5_Fede.png){width=80%}

#5. Cálculo de probabilidades. Estimación y comparación de resultados

Como se puede observar en las figuras de la sección anterior, las normales simuladas $(X_1,X_2)$ tienen una distribución conjunta cuyo soporte no es de forma rectangular. En este sentido, se puede descartar la hipótesis de que dichas variables sean independientes. Además, 

$$X_1 , X_2 \textrm{ son independientes} \iff \rho = 0$$

Esto se puede observar en la figura \ref{biv_rho0}, cuyo soporte es "más rectangular". Sin embargo, en las simulaciones solicitadas el coeficiente de correlación $\rho$ toma valores del intervalo $[0{.}5 ; 1]$ y más específicamente tendiendo $\rho \rightarrow 1$ en los casos evaluados, por lo que se puede considerar que $X_1$ y $X_2$ no son independientes. 

![Normal bivariada con $\rho = 0$\label{biv_rho0}](img/foo5_rho0.png){width=80%}

Por lo tanto, se puede pensar que:

$$\mathbb{P}(X_1 \leq 1 \cap X_1 \leq 1) \neq \mathbb{P}(X_1 \leq 1) \cdot \mathbb{P}(X_2 \leq 1)$$

$$\implies \mathbb{P}(X_1 \leq 1 \cap X_1 \leq 1) \neq (\Phi(1))^2$$

En este sentido, a partir del código provisto por la cátedra, se implementó un ciclo que permitiera contar los casos favorables de la simulación (es decir, que cumplieran la condición de la intersección) y dividiendolo por los casos totales (es decir, el tamaño de la simulación). A partir de esto se obtuvieron los siguientes resultados:

\hfill

|Padrón|$\rho$|Favorables|Totales|$\mathbb{P}(X_1 \leq 1 \cap X_1 \leq 1)$|$(\Phi(1))^2$|
|:------:|:-------:|:----------:|:-------:|:-----------------------:|:--------:|
|  96993 | 0.965   |        8165|    10000|                                  0.8165|      0.70785|
|  97112 | 0.560   |        7546|    10000|                                  0.7546|      0.70785|
|  97066 | 0.830   |        7834|    10000|                                  0.7834|      0.70785|
|  92892 | 0.960   |        8067|    10000|                                  0.8067|      0.70785|

En base al ejercicio **5.4)** de la guía de ejercicios, dada una normal bivariada $(X_1,X_2)$, se puede obtener que:

$$X_1|X_2=x_2 \sim \mathcal{N}\left(\mu=\mu_{x_1}+\rho \frac{\sigma_{x_1}}{\sigma_{x_2}} (x_2 - \mu_{x_2}) , \sigma^2 = \sigma_{x_1}^2 (1- \rho ^2)\right)$$

En este caso, se usó $\rho=0{.}965$  y que $X_1$ y $X_2$ son normales estándar. 

$$\implies X_1|X_2=1 \sim \mathcal{N}(\mu=0{.}965 , \sigma^2 = 0{.}068)$$
$$\implies \mathbb{P}(X_1 \leq 1|X_2=1) \approx 0{.}5531$$

Dado que la simulación tiene un rango finito de puntos, se buscó aumentar la cantidad de simulaciones,para que al marginar, la probabilidad no varíe  demasiado y se parezca al valor esperado. Por esto mismo, se aumentó la cantidad de número aleatorios (o pseudo-aleatorios) generados de 10 000 a 1 000 000. Además, se reemplazó la condición de $X_2 = 1$ por una condición más abarcativa para que se aproximaran mejor los puntos de la simulación: se utilizó la condición $0{.}95 \leq X_2 < 1{.}05$. Los resultados obtenidos son:

(PRESENTAR RESULTADOS EN TABLA CON CASOS FAVORABLES, CASOS TOTALES , PROB CALCULADA Y PROB ESPERADA. SI QUERES PODEMOS PONER LO QUE MANDE DE LAS SIMULACIONES DE 10 000 MUESTRAS Y ALGUNA DIFERENCIA O ERROR DE SIMULACION)

#6. Conclusiones      
