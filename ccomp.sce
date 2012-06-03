// ccomp: http://quijot.github.com/ccomp/ -------------------------------------
//
// por
//
//   Santiago Pestarini <santip@fceia.unr.edu.ar>
//   Laura Cornaglia    <cornagli@fceia.unr.edu.ar>
//
// se encuentra bajo una
//
//   Licencia Creative Commons Atribución-CompartirDerivadasIgual 2.5 Argentina
//   (ver más información en http://creativecommons.org/licenses/by-sa/2.5/ar/)
//
//   Usted es libre de compartir, copiar, distribuir, ejecutar y comunicar
//   públicamente la obra, hacer obras derivadas, hacer un uso comercial.
//-----------------------------------------------------------------------------


//-----------------------------------------------------------------------------
// Método de Medidas Indirectas Ponderadas
//-----------------------------------------------------------------------------
// Entrada
//    A: matriz de coeficientes
//    U: vector de términos independientes
//    P: matriz de pesos
//-----------------------------------------------------------------------------
// Salida
//    Q: matriz estructural o cofactor
//    X: vector de incógnitas/soluciones
//    sigma0p: precisión 'a posteriori' de la observación de peso 1
//    MVC: matriz Varianza-Covarianza
//    sigma: vector de precisiones del vector X
//-----------------------------------------------------------------------------
function [Q, X] = mindp_apriori(A, U, P)
  M = A'*P*A
  Q = inv(M)
  X = Q*A'*P*U
endfunction

function [Q, X, sigma0p, MVC, sigma] = mindp_aposteriori(A, U, P)
  [Q, X] = mindp_apriori(A, U, P)
  V = A*X - U
  sumpvv = V'*P*V
  // (n-m) es la sobreabundancia
  n = size(A,1) // nro de filas
  m = size(A,2) // nro de columnas
  sigma0p = sqrt(sumpvv/(n-m))
  MVC = sigma0p^2*Q
  sigma = sqrt(diag(MVC)) 
endfunction

//-----------------------------------------------------------------------------
// Elipse de Error
//-----------------------------------------------------------------------------
// Entrada
//    MVC: matriz Varianza-Covarianza
//-----------------------------------------------------------------------------
// Salida 
//    alfa: ángulo de giro para anular covarianzas
//    seeX, seeY: Semiejes de la elipse de error
//-----------------------------------------------------------------------------
function [alfa, seeY, seeX] = eliperr(MVC)
  autoval = spec(MVC)
  autovec1 = kernel(MVC-autoval(1) * eye(MVC))
  autovec2 = kernel(MVC-autoval(2) * eye(MVC))
  autovec = [autovec2,autovec1]
  MVCrot = clean(inv(autovec) * MVC * autovec)
  seeX = sqrt(MVCrot(1,1))
  seeY = sqrt(MVCrot(2,2))
  alfa = atan(autovec2(2) / autovec2(1))
endfunction

//-----------------------------------------------------------------------------
// Dibujo de elipse. Ecuación paramétrica. Rotación para ejes Gauss Krüger.
//-----------------------------------------------------------------------------
// Entrada
//    a, b: longitudes de los semiejes de la elipse
//    Ce: coordenadas del centro de la elipse, por ej.: [0,0]
//    alfa: ángulo de rotación (en rad)
//    titulo: título de la gráfica, por ej.: "Elipse rotada"
//-----------------------------------------------------------------------------
// Salida
//    Dibujo de la elipse (plot) en pantalla.
//-----------------------------------------------------------------------------
function plotEllipse(a, b, Ce, alfa, titulo)
  N = 25                // cantidad de segmentos para discretizar la elipse
  theta = 0:%pi/N:2*%pi // ángulo paramétrico
  // Calcular los pares (x,y) de la elipse
  x = a*cos(theta)
  y = b*sin(theta)
  R = [cos(alfa) -sin(alfa); sin(alfa) cos(alfa)] // matriz de rotación
  M = R * [x;y]                                   // elipse rotada
  // Centrar la elipse en el punto Ce
  x = M(1,:) + Ce(1)
  y = M(2,:) + Ce(2)
  // Características de la gráfica ---------------------------
  // gca: administra las propiedades de los ejes de la gráfica
  e = gca()         // equivalente a e = get("current_axes")
  e.isoview = 'on'  // habilita emplear la misma escala en ambos ejes
  // tamaño de los textos de la gráfica
  e.title.font_size = 5
  e.x_label.font_size = 3
  e.y_label.font_size = 3
  e.font_size = 2
  // textos: título y ejes
  xtitle(titulo, "Y [m]", "X [m]")
  // Dibujo de la elipse y del punto Ce
  plot(x,y,Ce(1),Ce(2),'.r')
endfunction

//-----------------------------------------------------------------------------
// Inserción de un punto en una red de coordenadas planas (2D)
//-----------------------------------------------------------------------------
// Entrada
// (pedida por cuadro de diálogo en pantalla, NO como parámetros de la función)
//    archivo: nombre de archivo de entrada
//    S: solución "aproximada"
//    e: estimación del error de posición de los puntos de la red 
//    ps: diferencia requerida entre soluciones de dos iteraciones sucesivas
//-----------------------------------------------------------------------------
// Formato del archivo de entrada
// (en columnas separadas por espacios o TAB)
//    Columnas 1 y 2: coordenadas X e Y de puntos de la red
//    Columna 3: distancias horizontales medidas a los puntos de la red
//    Columna 4: estimación del error de las distancias horizontales medidas
//    (Nota: La unidad de los valores numéricos es el metro)
//-----------------------------------------------------------------------------
// Salida
//    S: coordenadas del punto que se incorpora a la red
//    sigma: precisión de S
//    sigma0p: precisión 'a posteriori' de una observación de peso 1
//    MVC: Matriz Varianza-Covarianza
//    alfa: ángulo de giro de la elipse de error para anular covarianzas
//    seeY, seeX: semiejes de la elipse
//-----------------------------------------------------------------------------
function [seeY, seeX, alfa, MVC, sigma0p, sigma, S] = insertar_punto_en_red2D()
  // Configura el formato de salida para ver más decimales
  format("v", 15)
  
  // Entrada
  //---------------------------------------------------------------------------
  I = 2 // cantidad de incógnitas (para el plano)
  archivo = uigetfile("", ".", "Seleccione el archivo de datos (valores en metros)")
  [ok,S,e,ps] = getvalue("Ingrese los siguientes datos del problema", ["Coordenadas [X,Y] de solución ""aproximada"" (m)   S:"; "Precisión de la Red (m)   e:"; "Diferencia requerida entre soluciones iteraciones sucesivas (m)   ps:"], list("row", I, "vec", 1, "vec", 1), ["[0,0]"; "0.10"; "0.10"])
  
  // Lectura del archivo de entrada
  //---------------------------------------------------------------------------
  data = read(archivo,-1,I+2)
  n = size(data,1)  // cantidad de observaciones detectadas
  C = data(1:n,1:I) // coordenadas conocidas de los puntos de la red
  D = data(1:n,I+1) // distancias horizontales medidas
  E = data(1:n,I+2) // estimación de los errores de las distancias D
  
  // Cálculos
  //---------------------------------------------------------------------------
  // Relación funcional de las incógnitas aproximada por el polinomio de Taylor
  //---------------------------------------------------------------------------
  i = 0         // contador de iteraciones
  X = [2*ps;0]  // inicializa X de manera que entre al bucle
  K = 1         // precisión 'a priori' de la observación de peso 1
  while or(abs(X)>ps) do
    // itera mientras algún valor de X sea mayor que la "precisón requerida" ps
    i = i+1 // incrementa contador de iteraciones
    dX = S(1) - C(1:n,1)
    dY = S(2) - C(1:n,2)
    d = sqrt(dX^2+dY^2) // vector de distancias entre C y S
    // matrices del Sitema Normal para Método de Medidas Indirectas Ponderadas
    A = [dX./d, dY./d]
    U = D-d
    colP = K. / (E^2+e^2) // propagación de pesos
    P = diag(colP)
    // aplicación del Método de Medidas Indirectas Ponderadas a priori
    [Q, X] = mindp_apriori(A,U,P)
    // solución de la i-ésima iteración
    S = S + X'
  end
  S = S' // trasponer S

  // Método de Medidas Indirectas Ponderadas a posteriori
  // para obtener sigma y MVC
  //---------------------------------------------------------------------------
  [Q, X, sigma0p, MVC, sigma] = mindp_aposteriori(A,U,P)

  // Elipse de error
  //---------------------------------------------------------------------------
  [alfa, seeY, seeX] = eliperr(MVC)
  plotEllipse(seeY, seeX, [0,0], alfa, "Elipse de error")

  // Salida por patanlla
  printf("Leidas de archivo %s: %d observaciones.\nSolución encontrada en %d iteraciones.\n", archivo, n, i)
  disp(X, "Diferencia entre soluciones de 2 últimas iteraciones:")
endfunction

