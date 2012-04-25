// ccomp: https://github.com/quijot/ccomp -------------------------------------
//
//   por
//
//     Dra. Laura Cornaglia <cornagli@fceia.unr.edu.ar>
//     Santiago Pestarini   <santip@fceia.unr.edu.ar>
//
//   se encuentra bajo una
//
//   Licencia Creative Commons Atribución-CompartirDerivadasIgual 2.5 Argentina
//   (ver más información en http://creativecommons.org/licenses/by-sa/2.5/ar/)
//
//   Usted es libre de:
//       · Compartir - copiar, distribuir, ejecutar y comunicar públicamente la
//         obra
//       · hacer obras derivadas
//       · hacer un uso comercial de esta obra
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Resolución por Método de Medidas Indirectas
//    (como un caso particular de Medidas Indirectas con pesos 1)
//-----------------------------------------------------------------------------
function [Q, X, sigma] = mind(A, U)
  // Resolución por Mét de Med Ind con pesos 1
  I = eye(size(A,1),size(A,1)) // matriz identidad con mismas dimesiones que A
  [Q, X, sigma] = mindp(A,U,I)
endfunction

function [diagQ] = mind_apriori(A)
  // Devuelve sólo la diagonal principal de Q = M^-1
  I = eye(size(A,1),size(A,1)) // matriz identidad con mismas dimesiones que A
  U = ones(size(A,1),1) // vector columna de unos
  diagQ = diag(mindp(A,U,I))
endfunction

function [X, sigma] = mind_aposteriori(A, U)
  // Devuelve la solución con su calidad y la sumatoria de los v²
  I = eye(size(A,1),size(A,1)) // matriz identidad con mismas dimesiones que A
  [Q, X, sigma] = mindp(A,U,I)
endfunction

//-----------------------------------------------------------------------------
// Resolución por Método de Medidas Indirectas Ponderadas
//-----------------------------------------------------------------------------
// Entrada ----------------------------------------------
//    A: matriz de coeficientes
//    U: vector columna de términos independientes
//    P: matriz de pesos
// Salida -----------------------------------------------
//    Q: matriz estructural o cofactor
//    X: vector de incógnitas/soluciones
//    sigma: vector de precisiones de cada Xi
//    sigmao: precisión de la solución de peso 1
//-----------------------------------------------------------------------------
function [Q, X, sigma, sigmao] = mindp(A, U, P)
  M = A'*P*A
  Q = inv(M)
  X = Q*A'*P*U
  V = A*X - U
  sumpvv = V'*P*V
  // (n-m) es la sobreabundancia
  n = size(A,1) // nro de filas
  m = size(A,2) // nro de columnas
  sigmao = sqrt(sumpvv/(n-m))
  px = diag(Q)^-1              // Pesos de Xi
  sigma = sigmao./sqrt(px)     // Vector de precisiones de X
endfunction

//-----------------------------------------------------------------------------
// Resolución por Método de Observaciones Condicionadas
//    (como un caso particular de Observ. Condicionadas Ponderadas de peso 1)
//-----------------------------------------------------------------------------
function [X, V, sigma] = moc(X, B, W)
  // Resolución por Mét de Obs Cond con pesos 1
  I = eye(size(B,1),size(B,1)) // matriz identidad con mismas dimesiones que B
  [X, V, sigma] = mocp(X, B, W, I)
endfunction

//-----------------------------------------------------------------------------
// Resolución por Método de Observaciones Condicionadas Ponderadas
//-----------------------------------------------------------------------------
// Entrada ----------------------------------------------
//    X: vector de soluciones "aproximadas"
//    B: matriz de coeficientes
//    W: vector columna de términos independientes
//    P: matriz de pesos
// Salida -----------------------------------------------
//    X: vector de incógnitas/soluciones
//    V: vector de residuos
//    sigmao: precisión de la solución de peso 1
//-----------------------------------------------------------------------------
function [X, V, sigmao] = mocp(X, B, W, P)
  // Primero obtiene V
  N = B*inv(P)*B'
  K = inv(N)*(-W)
  V = inv(P)*B'*K
  // Luego la sumatoria [pvv] para calcular sigma0
  sumpvv = V'*P*V
  // I-(I-r) = r es la sobreabundancia (con I: cantidad de incógnitas)
  r = size(B,1)             // cantidad de filas
  sigmao = sqrt(sumpvv/r)   // precisión de la solución de peso 1
  // Por último devuelve el vector X ajustado 
  X = X + V
endfunction

//-----------------------------------------------------------------------------
// Cálculo de Matriz varianza-covarianza
//-----------------------------------------------------------------------------
function [Sigma] = matrizvc(Q, sigmao)
  Sigma = sigmao^2*Q
endfunction

//-----------------------------------------------------------------------------
// Cálculo de Elipse de Error
//-----------------------------------------------------------------------------
function [alfa, sigmaY, sigmaX] = eliperr(Q, sigmao)
  // ángulo de giro para anular covarianzas
  alfa = atan(2*Q(2,1)/(Q(1,1) - Q(2,2)))/2
  K = sqrt((2*Q(1,2))^2 + (Q(1,1)-Q(2,2))^2)
  Qyy = (Q(1,1) + Q(2,2) - K)/2
  Qxx = (Q(1,1) + Q(2,2) + K)/2
   // semiejes de la elipse
  sigmaY = sigmao*sqrt(Qyy)
  sigmaX = sigmao*sqrt(Qxx)
endfunction

//-----------------------------------------------------------------------------
// Función de dibujado de elipses rotada con centro en C = (x,y)
//-----------------------------------------------------------------------------
// Entrada ----------------------------------------------
//    a, b: longitudes de los semiejes de la elipse
//    C: coordenadas del centro de la elipse, por ej.: [0,0]
//    alfa: ángulo de rotación (en rad)
//    titulo: título de la gráfica, por ej.: "Elipse rotada"
//-----------------------------------------------------------------------------
function plotEllipse(a, b, C, alfa, titulo)
  // Rango para dibujo --------------------------------------------------------
  N = 25 // cantidad de segmentos para discretizar la elipse
  theta = 0:%pi/N:2*%pi
  // Calcular los pares (x,y) de la elipse ------------------------------------
  x = a*cos(theta)
  y = b*sin(theta)
  // Rotar la elipse con el angulo "alfa" -------------------------------------
  R = [cos(alfa) -sin(alfa); sin(alfa) cos(alfa)] // matriz de rotación
  M = R * [x;y] // matriz auxiliar
  // Centrar la elipse en el punto C ------------------------------------------
  x = M(1,:) + C(1)
  y = M(2,:) + C(2)
  // Ejes en escala -----------------------------------------------------------
  // gca me sirve para manejar propiedades de los ejes/gráfica
  e = gca()         // es lo mismo que e = get("current_axes")
  e.isoview = 'on'  // usar misma escala en ambos ejes
  // agrandar tamaño de los textos de la gráfica
  e.title.font_size = 5
  e.x_label.font_size = 3
  e.y_label.font_size = 3
  e.font_size = 2
  // establecer textos de título y ejes
  xtitle(titulo, "X [metros]", "Y [metros]")
  // Dibujar ------------------------------------------------------------------
  plot(x,y)            // elipse
  plot(C(1),C(2),'.r') // centro
endfunction

//-----------------------------------------------------------------------------
// Inserción de un punto en una red de coordenadas planas (2D)
// (1) Parámetros de entrada pedidos:
//       archivo: nombre de archivo con datos sistema de ecuaciones (Ver (2))
//       S: solución "aproximada"
//       e: error de la red (precisión de un punto cualquiera de la red)
//       ps: diferencia requerida entre soluciones de 2 iteraciones sucesivas
// (2) El archivo de entrada debe tener las siguientes columnas:
//       1 y 2: coordenadas X e Y (en metros) de puntos de la red
//       3: distancias horizontales (en metros) medidas a los puntos de la red
//       4: errores (en metros) de las distancias horizontales medidas
// (3) SALIDA:
//       S: coordenadas del punto que se incorpora a la red
//       sigma: precisión de las coordenadas de S
//       sigmao: precisión de una observación de peso 1
//       MVC: Matriz Varianza-Covarianza
//       alfa: ángulo de giro de la elipse para anular covarianzas
//       sY, sX: semiejes de la elipse
//-----------------------------------------------------------------------------
function [S, sigma, sigmao, MVC, alfa, sY, sX] = insertar_punto_en_red2D()
  // Configura formato de salida para ver más decimales
  format("v", 15)
  
  // ENTRADA DE DATOS por teclado \\
  I = 2 // cantidad de incógnitas (para el plano)
  archivo = uigetfile("", ".", "Seleccione el archivo de datos (valores en metros)")
  [ok,S,e,ps] = getvalue("Ingrese los siguientes datos del problema", ["Coordenadas [X,Y] de solución ""aproximada"" (m)   S:"; "Precisión de la Red (m)   e:"; "Diferencia requerida entre soluciones iteraciones sucesivas (m)   ps:"], list("row", I, "vec", 1, "vec", 1), ["[0,0]"; "0.10"; "0.10"])
  
  // CARGA DE DATOS desde archivo \\
  data = read(archivo,-1,I+2)
  n = size(data,1)  // cantidad de observaciones detectadas
  C = data(1:n,1:I) // coordenadas de puntos conocidos
  D = data(1:n,I+1) // distancias horizontales medidas
  E = data(1:n,I+2) // errores de dichas medidas
  
  // CÁLCULOS por Taylor \\
  i = 0         // contador de iteraciones
  X = [2*ps;0]  // inicializa X de manera que entre al bucle
  k = 1         // factor de escala
  // repetir hasta que todos los valores de X sean menores que "ps"
  while or(abs(X)>ps) do
    i = i+1 // incrementa contador de iteraciones
    dX = S(1) - C(1:n,1)
    dY = S(2) - C(1:n,2)
    d = sqrt(dX^2+dY^2) // vector de distancias a S
    // Prepara matrices para Método de Medidas Indirectas Ponderadas
    A = [dX./d, dY./d]
    U = D-d
    // propagar pesos
    eu = k * sqrt(E^2+e^2)
    colP = eu^-2
    P = diag(colP)
    // aplicar Mét. de Med. Ind. Ponderadas
    [Q, X, sigma, sigmao] = mindp(A,U,P)
    // Solución de la i-ésima iteración
    S = S + X'
  end
  S = S' // pasa de ser vector fila a vector columna

  // Matriz Varianza-Covarianza \\
  MVC = matrizvc(Q, sigmao)

  // Elipse \\
  [alfa, sY, sX] = eliperr(Q, sigmao)
  plotEllipse(sY, sX, [0,0], alfa, "Elipse de errores")

  // Imprimir salida \\
  printf("Leidas de archivo %s: %d observaciones.\nSolución encontrada en %d iteraciones.\n", archivo, n, i)
  disp(X, "Diferencia entre soluciones de 2 últimas iteraciones:")
endfunction
