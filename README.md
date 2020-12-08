# TP Final Bioinformática

## Prerequisitos
- [Instalar Conda](https://docs.conda.io/en/latest/)

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


### Known Issues
- Scrollbar: Se requiere maximizar/minimizar la ventana para que se muestre.
- Interfaz bloqueada al procesar.
- Al realizar un nuevo procesamiento, no se eliminan los botones de las cadenas anteriores.
