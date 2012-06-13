# ccomp

Funciones de cálculo de compensación para Scilab.

## Descripción

**ccomp** es un conjunto de scripts para Scilab que intentan implementar el principio de mínimos cuadrados mediante los métodos de medidas indirectas ponderadas y de observaciones condicionadas ponderadas, el cálculo de matrices varianza-covarianza, elipses de error y algunas otras funcionalidades útiles para el cálculo de compensación de procesos de medición.

## Instalación

**ccomp** no se _instala_. Para usarlo se debe tener instalado [Scilab](http://www.scilab.org/). Luego simplemente descargar **ccomp** y descomprimir su contenido donde guste.

## Uso

Para poder utilizar **ccomp** se debe _incorporar_ a Scilab desde el menú Archivo -> Ejecutar... (File -> Execute...) o desde la línea de comandos con

    exec <directorio-ccomp>/ccomp.sce;

donde \<directorio-ccomp\> es el directorio en que se encuentra **ccomp**.

Luego se pueden ejecutar sus funciones normalmente como cualquier otra función en Scilab, por ejemplo:

    [seeY, seeX, alfa, ro, MVC, sigma0, sigma, S] = insert_pto_red2D()

## Licencia

**ccomp** se encuentra bajo una Licencia [Creative Commons](http://creativecommons.org/licenses/by-sa/2.5/ar/). Usted es libre de compartir, copiar, distribuir, ejecutar y comunicar públicamente la obra, hacer obras derivadas, hacer un uso comercial.

