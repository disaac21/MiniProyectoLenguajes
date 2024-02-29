import csv
import math
import numpy as np

#Cargar Archivo y Crear Matriz
def cargarArchivo(filename):
    matriz = []
    with open(filename, 'r') as file:
        readCsv = csv.reader(file, delimiter=';')
        for fila in readCsv:
            matriz.append(fila)
    return matriz

#Centrar y Reducir
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

#Calculo de Matriz de Correlaciones
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

#Calculo de Valores y Vectores Propios
def calcValoresPropios(matrizCorrelaciones):
    valoresPropios = np.linalg.eigvals(matrizCorrelaciones)
    return valoresPropios

def calcVectoresPropios(matrizCorrelaciones):
    _, vectoresPropios = np.linalg.eig(matrizCorrelaciones)
    return vectoresPropios

#Ordenamiento de Valores y Vectores Propios & Construccion de Matriz V
def orderValoresPropios(valoresPropios):
    orderedValoresPropios = valoresPropios
    for i in range(len(orderedValoresPropios)):
        for j in range(i+1, len(orderedValoresPropios)):
            if orderedValoresPropios[i] < orderedValoresPropios[j]:
                orderedValoresPropios[i], orderedValoresPropios[j] = orderedValoresPropios[j], orderedValoresPropios[i]
    return orderedValoresPropios

def orderVectoresPropios(vectoresPropios, orderedValoresPropios):
    indexes = np.argsort(orderedValoresPropios)[::-1]
    orderedVectoresPropios = vectoresPropios[:, indexes]
    return orderedVectoresPropios

#Calculo de Matriz de Componentes Principales
def calcMatrizComponentesPrincipales(centradaReducida, matrizV):
    filas_CR, columnas_CR = centradaReducida.shape
    filas_V, columnas_V = matrizV.shape

    matrizComponentesPrincipales = np.zeros((filas_CR, columnas_V))

    for i in range(filas_CR):
        matrizComponentesPrincipales[i] = np.dot(centradaReducida[i], matrizV)

    return matrizComponentesPrincipales

#Calculo de Matriz de Calidades de Individuos
def calcMatrizCalidadesIndividuos(matrizComponentesPrincipales, matrizCentradaReducida):
    listaMatrizCalidadesIndividuos = []
    rowCount = len(matrizComponentesPrincipales)
    columnCount = len(matrizComponentesPrincipales[0])
    
    for i in range(rowCount):
        fila = []
        for j in range(columnCount):
            num = matrizComponentesPrincipales[i][j] ** 2
            den = sum(matrizCentradaReducida[i][k] ** 2 for k in range(columnCount))
            fila.append(num / den)
        listaMatrizCalidadesIndividuos.append(fila)
    
    matrizCalidadesIndividuos = np.array(listaMatrizCalidadesIndividuos)
    return matrizCalidadesIndividuos

# Main ----------     ----------     ----------     ----------     ----------

#Cargar Archivo y Crear Matriz
filename = 'EjemploEstudiantes.csv'
matrizRaw = cargarArchivo(filename)

#Centrar y Reducir
matrizCentradaReducida = centrarReducir(matrizRaw) 

#Calculo de Matriz de Correlaciones
matrizCorrelaciones = calcMatrizCorrelaciones(matrizCentradaReducida)

#Calculo de Valores y Vectores Propios
ValoresPropios = calcValoresPropios(matrizCorrelaciones)
VectoresPropios = calcVectoresPropios(matrizCorrelaciones)

#Ordenamiento de Valores y Vectores Propios
OrdenadosValoresPropios = orderValoresPropios(ValoresPropios)
OrdenadosVectoresPropios = orderVectoresPropios(VectoresPropios, OrdenadosValoresPropios)

#Construccion de Matriz V
# matrizV = np.column_stack([columna for columna in OrdenadosVectoresPropios.T])
matrizV = OrdenadosVectoresPropios

#Calculo de Matriz de Componentes Principales
matrizComponentesPrincipales = calcMatrizComponentesPrincipales(matrizCentradaReducida, matrizV)

#Calculo de Matriz de Calidades de Individuos
matrizCalidadesIndividuos = calcMatrizCalidadesIndividuos(matrizComponentesPrincipales, matrizCentradaReducida)

print("\nMatriz Centrada y Reducida:")
print(matrizCentradaReducida)
print("\nMatriz de Correlaciones:")
print(matrizCorrelaciones)
print()

for i in range(len(OrdenadosValoresPropios)):
    print("Valor propio:", OrdenadosValoresPropios[i])
    print("Vector propio:", OrdenadosVectoresPropios[:, i])
    print()

print("Matriz V:")
print(matrizV)
print()

print("Matriz de Componentes Principales:")
print(matrizComponentesPrincipales)
print()

print("Matriz de Calidades de Individuos:")
print(matrizCalidadesIndividuos)
print()