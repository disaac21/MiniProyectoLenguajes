import csv
import math
import numpy as np

def cargarArchivo(filename):
    matriz = []
    with open(filename, 'r') as file:
        readCsv = csv.reader(file, delimiter=';')
        for fila in readCsv:
            matriz.append(fila)
    return matriz

def calcMean(matriz):
    num_filas = len(matriz)
    num_columnas = len(matriz[0])
    medias = []
    for j in range(num_columnas):
        suma_columna = sum(matriz[i][j] for i in range(num_filas))
        media_columna = suma_columna / num_filas
        medias.append(media_columna)
    return medias

def calcStdDev(matriz):
    meanColumnas = calcMean(matriz)
    num_filas = len(matriz)
    num_columnas = len(matriz[0])
    stdDeviations = []
    for j in range(num_columnas):
        sumSquares = sum((matriz[i][j] - meanColumnas[j]) ** 2 for i in range(num_filas))
        columns_StdDev = math.sqrt(sumSquares / num_filas)
        stdDeviations.append(columns_StdDev)
    return stdDeviations

def centrarReducir(matriz):
    datos = []
    for fila in matriz[1:]:
        filaDecimales = [float(valor.replace(',', '.')) for valor in fila[1:]]
        datos.append(filaDecimales)
    matrizDecimales = np.array(datos)
    means = calcMean(matrizDecimales)
    stdDeviations = calcStdDev(matrizDecimales)
    matrizCentradaReducida = (matrizDecimales - means) / stdDeviations
    return matrizCentradaReducida

import numpy as np

def calcCorrelaciones(columnOne, columnTwo):
    meansOne = np.mean(columnOne)
    meansTwo = np.mean(columnTwo)
    numerador = np.sum((columnOne - meansOne) * (columnTwo - meansTwo))
    denominador = np.sqrt(np.sum((columnOne - meansOne)**2)) * np.sqrt(np.sum((columnTwo - meansTwo)**2))
    correlacionSend = numerador / denominador
    return correlacionSend

def calcMatrizCorrelaciones(matrizCentradaReducida):
    n = len(matrizCentradaReducida[0])
    matrizCorrelaciones = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            r = calcCorrelaciones(matrizCentradaReducida[:, i], matrizCentradaReducida[:, j])
            matrizCorrelaciones[i, j] = r
    return matrizCorrelaciones

filename = 'EjemploEstudiantes.csv'
matrizRaw = cargarArchivo(filename)

matrizCentradaReducida = centrarReducir(matrizRaw)
matriz_correlaciones = calcMatrizCorrelaciones(matrizCentradaReducida)

print("\nMatriz Centrada y Reducida:")
print(matrizCentradaReducida)
print("\nMatriz de Correlaciones:")
print(matriz_correlaciones)