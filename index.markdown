---
title: ccomp
layout: default
---
# ¿Qué es **ccomp**?

**ccomp** es un conjunto de scripts para [Scilab][scilab] que implementan el principio de mínimos cuadrados mediante los métodos de medidas indirectas ponderadas y de observaciones condicionadas ponderadas, el cálculo de matrices varianza-covarianza, elipses de error y algunas otras funcionalidades útiles para el cálculo de compensación de procesos de medición.

[scilab]: http://www.scilab.org/  "Sitio oficial de Scilab"

## Instalación

**ccomp** no se _instala_. Para usarlo se debe tener instalado [Scilab][scilab]. Luego simplemente descargar **ccomp** y descomprimir su contenido donde guste.

## Uso

Para poder utilizar **ccomp** se debe _incorporar_ a Scilab desde el menú Archivo -> Ejecutar... (File -> Execute...) o desde la línea de comandos con

    exec <dir-ccomp>/ccomp.sce;

donde _&lt;dir-ccomp&gt;_ es la ruta del directorio en que se encuentra **ccomp**.

Luego se pueden ejecutar sus funciones normalmente como cualquier otra función en Scilab, por ejemplo:

    [seeY, seeX, alfa, ro, MVC, sigma0, sigma, S] = insert_pto_red2D()

o, para ejecutar el Método de Medidas Indirectas Ponderadas:

    [Q, X, sigma0, MVC, sigma] = mindp(A, U, P)

donde las matrices de entrada _A_, _U_ y _P_ son las matrices de coeficientes, de términos independientes y de pesos respectivamente; y se obtiene como salida _Q_, la matriz modal; _X_, el vector de soluciones; _sigma0_, la precisión de la observación de peso 1; _MVC_ la matriz varianza-covarianza y _sigma_, vector de precisiones de las soluciones en X.

## Licencia

**ccomp** se encuentra bajo una Licencia [Creative Commons][cc]. Usted es libre de compartir, copiar, distribuir, ejecutar y comunicar públicamente la obra, hacer obras derivadas.

[cc]: http://creativecommons.org/licenses/by-nc-sa/2.5/ar/  "Licencia Creative Commons 2.5"

