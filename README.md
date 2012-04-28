# ccomp

Funciones de cálculo de compensación para Scilab.

## Descripción

**ccomp** es un conjunto de scripts para Scilab que implementan el principio de mínimos cuadrados mediante los métodos de medidas indirectas ponderadas y de observaciones condicionadas ponderadas, el cálculo de matrices varianza-covarianza, elipses de error y algunas otras funcionalidades útiles para el cálculo de compensación de procesos de medición.

## Instalación

**ccomp** no se _instala_. Para usarlo se debe tener instalado [Scilab](http://www.scilab.org/). Luego simplemente descargar **ccomp** y descomprimir su contenido donde guste.

## Uso

Para poder utilizar **ccomp** se debe _incorporar_ a Scilab desde el menú Archivo -> Ejecutar... (File -> Execute...) o desde la línea de comandos con

    exec <directorio-ccomp>/ccomp.sce;

donde \<directorio-ccomp\> es el directorio en que se encuentra **ccomp**.
    
Luego se pueden ejecutar sus funciones normalmente como cualquier otra función en Scilab, por ejemplo:

    [S, sigma, sigmao, MVC, alfa, sY, sX] = insertar_punto_en_red2D()
