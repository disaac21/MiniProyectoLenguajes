import csv
import math
import numpy as np

def leer_archivo_csv(nombre_archivo):
    matriz = []
    with open(nombre_archivo, 'r') as archivo:
        lector_csv = csv.reader(archivo, delimiter=';')
        for fila in lector_csv:
            matriz.append(fila)
    return matriz

def calcular_media(matriz):
    num_filas = len(matriz)
    num_columnas = len(matriz[0])
    medias = []
    for j in range(num_columnas):
        suma_columna = sum(matriz[i][j] for i in range(num_filas))
        media_columna = suma_columna / num_filas
        medias.append(media_columna)
    return medias

def calcular_desviacion_estandar(matriz):
    medias_columnas = calcular_media(matriz)
    num_filas = len(matriz)
    num_columnas = len(matriz[0])
    desviaciones_estandar = []
    for j in range(num_columnas):
        suma_cuadrados = sum((matriz[i][j] - medias_columnas[j]) ** 2 for i in range(num_filas))
        desviacion_estandar_columna = math.sqrt(suma_cuadrados / num_filas)
        desviaciones_estandar.append(desviacion_estandar_columna)
    return desviaciones_estandar

def centrar_y_reducir(matriz):
    datos = []
    for fila in matriz[1:]:
        fila_numerica = [float(valor.replace(',', '.')) for valor in fila[1:]]
        datos.append(fila_numerica)
    matriz_numerica = np.array(datos)
    medias = calcular_media(matriz_numerica)
    desviaciones_estandar = calcular_desviacion_estandar(matriz_numerica)
    matriz_centralizada_reducida = (matriz_numerica - medias) / desviaciones_estandar
    return matriz_centralizada_reducida

nombre_archivo = 'EjemploEstudiantes.csv'
matriz_datos = leer_archivo_csv(nombre_archivo)

matriz_centralizada_reducida = centrar_y_reducir(matriz_datos)

print("Matriz centrada y reducida:")
print(matriz_centralizada_reducida)