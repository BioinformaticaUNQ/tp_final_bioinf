# TP Final Bioinformática

## Prerequisitos
- [Instalar Conda](https://docs.conda.io/en/latest/)
- Distribución de GNU/Linux.

## Instalación del entorno
- Clonar repositorio
- Ejecutar en la carpeta root del proyecto el comando: 
```
conda env create --file enviroment.yml
```
- Activar el entorno ejecutando el comando: 
```
conda activate bioinf
```
- Ejecutar el script con el comando: 
```
python main.py
```

### [Cheatsheet conda](https://docs.conda.io/projects/conda/en/latest/_downloads/843d9e0198f2a193a3484886fa28163c/conda-cheatsheet.pdf)

### Instrucciones de uso

1. Ingresar el código PDB de la proteína problema en el campo de tipo texto indicado.

2. Seleccionar los valores de e-value y coverage deseados para la búsqueda de proteínas homólogas.

3. Presionar el botón "Procesar" para la validación de la proteína problema ingresada.

4. Seleccionar la cadena a procesar.

5. Esperar hasta que la interfaz indique que se terminó el proceso de búsqueda y alineamiento.

6. Podrá observar los gráficos que indican el grado de conservación de las estructuras primarias y secundarias. Así como el alineamiento de las estructuras terciarias.

7. Dentro de la carpeta "ejecucion-{fecha y hora de la misma}" se encuentran los logs indicando cada proceso realizado, un workspace de pymol para observar con mayor detalle el alineamiento y el alineamiento que produce clustal.

8. Abrir el archivo .pse de pymol con una version igual o arriba de 2.4 (https://pymol.org/2/).

### Notas
- Puede volver a ejecutar una busqueda pulsando nuevamente el boton "Procesar"
- Para salir del modo pantalla completa, presione la tecla "Esc"
- En la primer ejecucion puede demorar la interfaz ya que se descarga la base de datos PDB

### Known Issues
- Scrollbar: Se requiere maximizar/minimizar la ventana para que se muestre.
- Interfaz bloqueada al procesar.
- La ejecucion de pymol puede demorar mas de lo pensado o no responder.
