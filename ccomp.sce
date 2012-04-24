//-----------------------------------------------------------------------------
// Resolución por Método de Medidas Indirectas
//-----------------------------------------------------------------------------
function [Q, X, sigma, sumvv] = mind(A, U)
  // Resolución por Mét de Med Ind con pesos 1
  I = eye(size(A,1),size(A,1)) // matriz identidad con mismas dimesiones que A
  [Q, X, sigma, sumvv] = mindp(A,U,I)
endfunction

function [diagQ] = mind_apriori(A)
  // Devuelve sólo la diagonal principal de Q = M^-1
  I = eye(size(A,1),size(A,1)) // matriz identidad con mismas dimesiones que A
  U = ones(size(A,1),1) // vector columna de unos
  diagQ = diag(mindp(A,U,I))
endfunction

function [X, sigma, sumvv] = mind_aposteriori(A, U)
  // Devuelve la solución con su calidad y la sumatoria de los v²
  I = eye(size(A,1),size(A,1)) // matriz identidad con mismas dimesiones que A
  [Q, X, sigma, sumvv] = mindp(A,U,I)
endfunction

//-----------------------------------------------------------------------------
// Resolución por Método de Medidas Indirectas Ponderadas
//-----------------------------------------------------------------------------
function [Q, X, sigma, sumpvv, sigmao] = mindp(A, U, P)
  M = A'*P*A
  Q = inv(M)    // Matriz Estructural
  X = Q*A'*P*U  // Valores de las incógnitas
  V = A*X - U   // Valores de sus residuos
  sumpvv = V'*P*V
  // (n-m) es la sobreabundancia
  n = size(A,1) // nro de filas
  m = size(A,2) // nro de columnas
  sigmao = sqrt(sumpvv/(n-m))  // EMC de la solución de peso 1
  dQ = diag(Q)
  px = dQ^-1,  // calcula el peso de xi
  for i=1:m do 
    sigma(i)=sigmao/sqrt(px(i)),  // calcula el sigma de xi
  end
endfunction

//-----------------------------------------------------------------------------
// Resolución por Método de Observaciones Condicionadas
//-----------------------------------------------------------------------------
function [X, V, sigma, sumvv] = moc(X, B, W)
  // Resolución por Mét de Obs Cond con pesos 1
  I = eye(size(B,1),size(B,1)) // matriz identidad con mismas dimesiones que B
  [X, V, sigma, sumvv] = mocp(X, B, W, I)
endfunction

function [X, V, sigmao, sumpvv] = mocp(X, B, W, P)
  // Primero obtiene V
  N = B*inv(P)*B'
  K = inv(N)*(-W)
  V = inv(P)*B'*K
  // Luego la sumatoria [pvv] para calcular sigma0
  sumpvv = V'*P*V
  // I-(I-r) = r es la sobreabundancia (con I cant de incóg)
  r = size(B,1) // nro de filas
  sigmao = sqrt(sumpvv/r)  // EMC de la solución de peso 1
  // Por último devuelve el vector X ajustado 
  X = X + V
endfunction

//-----------------------------------------------------------------------------
// Inserción de un punto en una red de coordenadas planas
// (1) Parámetros de entrada pedidos:
//       archivo: nombre de archivo con datos sistema de ecuaciones (Ver (2))
//       S: solución aproximada
//       e: error de la red (precisión de un punto cualquiera de la red)
//       ps: diferencia requerida entre soluciones de 2 iteraciones contíguas
// (2) El archivo de entrada debe tener las siguientes columnas:
//       1 y 2: coordenadas X e Y (en metros) de puntos de la red
//       3: distancias  D (en metros) medida a los puntos de la red
//       4: precisiones/errores E (en metros) de las D
//-----------------------------------------------------------------------------
function [S, sigma, pvv, sigmao] = insertar_punto_en_red2D()
  // Configura formato de salida para ver más decimales
  format("v", 15)
  
  // ENTRADA DE DATOS \\
  I = 2 // cantidad de incógnitas (para el plano)
  archivo = uigetfile("", ".", "Seleccione el archivo de datos (valores en metros)")
  [ok,S,e,ps] = getvalue("Ingrese los siguientes datos del problema", ["Coordenadas [X,Y] de solución estimada (m)   S:"; "Precisión de la Red (m)   e:"; "Diferencia requerida entre soluciones iteraciones contíguas (m)   ps:"], list("row", I, "vec", 1, "vec", 1), ["[0,0]"; "0.10"; "0.10"])
  
  // CARGA DE DATOS \\
  data = read(archivo,-1,I+2)
  n = size(data,1)  // cantidad de observaciones detectadas
  C = data(1:n,1:I) // coordenadas de puntos conocidos
  D = data(1:n,I+1) // distancias medidas al punto que se quiere incorporar a la red
  E = data(1:n,I+2) // precisiones de dichas medidas
  
  // CÁLCULOS \\
  i = 0 // contador de iteraciones
  X = [2*ps;0] // inicializa X de manera que entre al bucle
  while or(abs(X)>ps) do // hasta que todos los valores de X sean menores que la tolerancia
    i = i+1 // incrementa contador de iteraciones
    dX = S(1) - C(1:n,1)
    dY = S(2) - C(1:n,2)
    d = sqrt(dX^2+dY^2) // vector de distancias al punto solución estimado
    // Prepara matrices para Método de Medidas Indirectas Ponderadas
    A = [dX./d, dY./d]
    U = D-d
    eu = sqrt(E^2+e^2)
    colP = eu^-2
    P = diag(colP)
    [Q, X, sigma, pvv, sigmao] = mindp(A,U,P)
    // Solución de la i-ésima iteración
    S = S + X'
  end
  S = S' // pasa de ser vector fila a vector columna
  printf("Leidas de archivo %s: %d observaciones.\nSolución encontrada en %d iteraciones.\n", archivo, n, i)
  disp(X, "Diferencia entre soluciones de 2 últimas iteraciones:")
endfunction

//-----------------------------------------------------------------------------
// Errores axiales de un teodolito o Estación Total
//-----------------------------------------------------------------------------
// Cálculo "a priori"
//-----------------------------------------------------------------------------
function [diagQ, sigmac, sigmai] = tp3_apriori()
  // Configura formato de salida para ver más decimales
  format("v", 15)

  // ENTRADA DE DATOS
  I = 2 // cantidad de incógnitas (errores de colimación c y de inclinación i)
  [ok,R,n,sigmau] = getvalue("Ingrese los siguientes datos del problema", ["Rango [d,a] del círuclo vertical para observaciones (rad):"; "Cantidad de observaciones:"; "Precisión del instrumntal (segundos):"], list("row", 2, "vec", 1, "vec", 1), ["[-1,1]"; "20"; "10"])

  // CARGA DE DATOS
  dif = (R(2)-R(1))/(n-1)
  O = [R(1):dif:R(2)]'
  A = [sec(O) tan(O)]
  
  // CÁLCULOS
  diagQ = mind_apriori(A)
  sigmac = sqrt(diagQ(1)*sigmau^2)
  sigmai = sqrt(diagQ(2)*sigmau^2)
endfunction

//-----------------------------------------------------------------------------
// Cálculo "a posteriori" (valores de retorno en segundos sexagesimales, '')
//-----------------------------------------------------------------------------
function [X, sigma, sumvv] = tp3_aposteriori()
  // Configura formato de salida para ver más decimales
  format("v", 15)

  // ENTRADA DE DATOS
  I = 2 // cantidad de incógnitas (errores de colimación c y de inclinación i)
  archivo = uigetfile("", ".", "Seleccione el archivo de observaciones Ld, Li, h")

  // CARGA DE DATOS
  data = read(archivo,-1,I+1) // I+1 porque se carga también el ángulo en altura h
  n = size(data,1)  // cantidad de observaciones detectadas
  Ld = data([1:n],1)
  Lv = (data([1:n],1) + data([1:n],2) - 180) / 2
  h = data([1:n],3) * %pi/180 // pasado a rad

  // CÁLCULOS
  A = [sec(h) tan(h)]
  U = Ld - Lv
  [X, sigma, sumvv] = mind_aposteriori(A,U)
  X = X * 3600
  sigma = sigma * 3600
  sumvv = sumvv * 3600^2
endfunction

//-----------------------------------------------------------------------------
// Cálculo de Matriz varianza-covarianza
//-----------------------------------------------------------------------------
function [sigma, Q] = matrizvcAP(A, P, sigmao)
  Q = inv(A'*P*A)
  sigma = sigmao^2*Q
endfunction

function [sigma] = matrizvcQ(Q, sigmao)
  sigma = sigmao^2*Q
endfunction

//-----------------------------------------------------------------------------
// Cálculo de Elipse de Error
//-----------------------------------------------------------------------------
function [alfa, K, Qyy, Qxx, sigmaY, sigmaX] = eliperr(Q, sigmao)
  alfa = atan(2*Q(2,1)/(Q(1,1) - Q(2,2)))/2 // ángulo de giro para anular covarianzas
  K = sqrt((2*Q(1,2))^2 + (Q(1,1)-Q(2,2))^2)
  Qyy = (Q(1,1) + Q(2,2) - K)/2
  Qxx = (Q(1,1) + Q(2,2) + K)/2
  sigmaY = sigmao*sqrt(Q(2,2)) // radio de la elipse
  sigmaX = sigmao*sqrt(Q(1,1)) // radio de la elipse
endfunction
