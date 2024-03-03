import csv
import math
import numpy as np
import matplotlib.pyplot as plt

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

    matrizComponentesPrincipalesTemp = np.zeros((filas_CR, columnas_V))

    for i in range(filas_CR):
        matrizComponentesPrincipalesTemp[i] = np.dot(centradaReducida[i], matrizV)

    matrizComponentesPrincipales = matrizComponentesPrincipalesTemp[:, :2]
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

#Calculo de Matriz de Coordenadas de las Variables
def calcMatrizCoordenadasVariables(centradaReducida, vectoresPropios):
    matrizCoordenadasVariables = np.dot(centradaReducida, vectoresPropios.T)
    return matrizCoordenadasVariables

#Calculo de Matriz de Calidades de las Variables
def calcMatrizCalidadesVariables(valoresPropios, entries):
    matrizCalidadesVariables = valoresPropios / entries
    return matrizCalidadesVariables

#Calculo de Vector de Inercias de los Ejes
def calcVectorInerciasEjes(valoresPropios):
    vectorInerciasEjes = np.array([])
    den = len(valoresPropios)
    for i in range(len(valoresPropios)):
        num = valoresPropios[i] * 100
        value = num/den
        vectorInerciasEjes = np.append(vectorInerciasEjes, value)
    return vectorInerciasEjes

#Extrae Headings
def fileHeadings(filename):
    with open(filename, 'r', newline='') as file:
        readCsv = csv.reader(file, delimiter=';')
        linea = next(readCsv)
    linea = linea[1:]
    return linea

#Extrae Nombres
def fileNames(filename):
    names = []
    with open(filename, 'r', newline='') as file:
        readCsv = csv.reader(file, delimiter=';')
        next(readCsv)
        for row in readCsv:
            names.append(row[0])
    return names

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
matrizV = OrdenadosVectoresPropios

#Calculo de Matriz de Componentes Principales
matrizComponentesPrincipales = calcMatrizComponentesPrincipales(matrizCentradaReducida, matrizV)

#Calculo de Matriz de Calidades de Individuos
matrizCalidadesIndividuos = calcMatrizCalidadesIndividuos(matrizComponentesPrincipales, matrizCentradaReducida)

#Calculo de Matriz de Coordenadas de las Variables
matrizCoordenadasVariables = calcMatrizCoordenadasVariables(matrizCentradaReducida, OrdenadosVectoresPropios)

#Calculo de Matriz de Calidades de las Variables
matrizCalidadesVariables = calcMatrizCalidadesVariables(OrdenadosValoresPropios, len(matrizCentradaReducida))

#Calculo de Vector de Inercias de los Ejes
vectorInerciasEjes = calcVectorInerciasEjes(OrdenadosValoresPropios)

print("\nPaso 1 - Matriz Centrada y Reducida:")
print(matrizCentradaReducida)
print("\nPaso 2 - Matriz de Correlaciones:")
print(matrizCorrelaciones)
print()

print("\nPaso 3 - Valores y Vectores Propios:")
for i in range(len(OrdenadosValoresPropios)):
    print("Valor propio:", OrdenadosValoresPropios[i])
    print("Vector propio:", OrdenadosVectoresPropios[:, i])
    print()

print("Matriz de Vectores Propios:")
print(OrdenadosVectoresPropios)
print()

print("Matriz de Valores Propios:")
print(OrdenadosValoresPropios)
print()

print("Paso 4 - Matriz V:")
print(matrizV)
print()

print("Paso 5 - Matriz de Componentes Principales:")
print(matrizComponentesPrincipales)
print()

print("Paso 6 - Matriz de Calidades de Individuos:")
print(matrizCalidadesIndividuos)
print()

print("Paso 7 - Matriz de Coordenadas de las Variables:")
print(matrizCoordenadasVariables)
print()

print("Paso 8 - Matriz de Calidades de las Variables:")
print(matrizCalidadesVariables)
print()

print("Paso 9 - Vector de Inercias de los Ejes:")
print(vectorInerciasEjes)
print()

# Circular Graph ----------     ----------     ----------     ----------     ----------

arrowHeads = fileHeadings(filename)

fig, ax = plt.subplots(figsize=(8, 8))
ax.set_xlim(-1.1, 1.1)
ax.set_ylim(-1.1, 1.1)
ax.axhline(0, color='gray', linewidth=2, linestyle='dotted', dashes=(3, 3))
ax.axvline(0, color='gray', linewidth=2, linestyle='dotted', dashes=(3, 3))

circle = plt.Circle((0, 0), 1, color='orange', fill=False)
ax.add_artist(circle)

for i, var in enumerate(arrowHeads):
    ax.arrow(0, 0, OrdenadosVectoresPropios[i][0], OrdenadosVectoresPropios[i][1], head_width=0.05, head_length=0.05, fc='orange', ec='orange')
    ax.text(OrdenadosVectoresPropios[i][0]*1.15, OrdenadosVectoresPropios[i][1]*1.15, var, color='black', ha='center', va='center')

plt.title("Círculo de Correlación")
plt.xlabel("Componente 1")
plt.ylabel("Componente 2")
plt.grid(True)
plt.show()

# Regular Graph ----------     ----------     ----------     ----------     ----------

nameLabels = fileNames(filename)

x = matrizComponentesPrincipales[:, 0]
y = matrizComponentesPrincipales[:, 1]

fig, ax = plt.subplots(figsize=(10, 10))
ax.scatter(x, y, color='blue')

for i, nombre in enumerate(nameLabels):
    ax.text(x[i], y[i], nombre, fontsize=9)

ax.set_title("Plano Principal")
ax.set_xlabel("Componente 1")
ax.set_ylabel("Componente 2")
ax.grid(True)
plt.show()