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

def calcular_media(matriz, axis=None):
    if axis is None:
        # Calcular la media de todos los elementos de la matriz
        total_elementos = sum(sum(fila) for fila in matriz)
        total_elementos = float(total_elementos)
        return total_elementos / (len(matriz) * len(matriz[0]))
    elif axis == 0:
        # Calcular la media a lo largo de las columnas
        num_filas = len(matriz)
        num_columnas = len(matriz[0])
        medias = []
        for j in range(num_columnas):
            suma_columna = sum(matriz[i][j] for i in range(num_filas))
            media_columna = suma_columna / num_filas
            medias.append(media_columna)
        return medias
    elif axis == 1:
        # Calcular la media a lo largo de las filas
        medias = []
        for fila in matriz:
            media_fila = sum(fila) / len(fila)
            medias.append(media_fila)
        return medias
    else:
        raise ValueError("El valor de 'axis' debe ser None, 0 o 1.")
    
def calcular_desviacion_estandar(matriz, axis=None):
    if axis is None:
        # Calcular la desviación estándar de todos los elementos de la matriz
        media = calcular_media(matriz)
        total_elementos = sum(sum(fila) for fila in matriz)
        total_elementos = float(total_elementos)
        suma_cuadrados = sum((valor - media) ** 2 for fila in matriz for valor in fila)
        desviacion_estandar = math.sqrt(suma_cuadrados / total_elementos)
        return desviacion_estandar
    elif axis == 0:
        # Calcular la desviación estándar a lo largo de las columnas
        medias_columnas = calcular_media(matriz, axis=0)
        num_filas = len(matriz)
        num_columnas = len(matriz[0])
        desviaciones_estandar = []
        for j in range(num_columnas):
            suma_cuadrados = sum((matriz[i][j] - medias_columnas[j]) ** 2 for i in range(num_filas))
            desviacion_estandar_columna = math.sqrt(suma_cuadrados / num_filas)
            desviaciones_estandar.append(desviacion_estandar_columna)
        return desviaciones_estandar
    elif axis == 1:
        # Calcular la desviación estándar a lo largo de las filas
        desviaciones_estandar = []
        for fila in matriz:
            media_fila = calcular_media([fila])
            suma_cuadrados = sum((valor - media_fila) ** 2 for valor in fila)
            desviacion_estandar_fila = math.sqrt(suma_cuadrados / len(fila))
            desviaciones_estandar.append(desviacion_estandar_fila)
        return desviaciones_estandar
    else:
        raise ValueError("El valor de 'axis' debe ser None, 0 o 1.")

def centrar_y_reducir(matriz):
    # Convertir la matriz de texto a una matriz numérica
    datos = []
    for fila in matriz[1:]:
        fila_numerica = [float(valor.replace(',', '.')) for valor in fila[1:]]  # Convertir valores a float, reemplazando comas por puntos
        datos.append(fila_numerica)
    
    matriz_numerica = np.array(datos)
    
    # Calcular la media y la desviación estándar de cada columna
    # medias = np.mean(matriz_numerica, axis=0)
    medias = calcular_media(matriz_numerica, axis=0)
    desviaciones_estandar = calcular_desviacion_estandar(matriz_numerica, axis=0)
    
    # Centrar y reducir cada columna
    matriz_centralizada_reducida = (matriz_numerica - medias) / desviaciones_estandar
    
    return matriz_centralizada_reducida

nombre_archivo = 'EjemploEstudiantes.csv'
matriz_datos = leer_archivo_csv(nombre_archivo)

matriz_centralizada_reducida = centrar_y_reducir(matriz_datos)

print("Matriz centrada y reducida:")
print(matriz_centralizada_reducida)
