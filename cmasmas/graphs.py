import csv
import numpy as np
import matplotlib.pyplot as plt

filename = "EjemploEstudiantes.csv"

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

OrdenadosVectoresPropios = np.array([
    [-0.52664397, -0.2704963,  -0.43820071, -0.62387762, -0.26121779],
    [-0.42493622, -0.50807221, -0.04049491,  0.32538951,  0.67362724],
    [-0.35914704,  0.56208159, -0.56227583,  0.48374732, -0.07008647],
    [-0.35269747,  0.58648985,  0.39418032, -0.42043348,  0.44664495],
    [ 0.53730181,  0.09374599, -0.57862603, -0.30679407,  0.52305619]
])

matrizComponentesPrincipales = np.array([
    [-0.32306263,  1.7725245 ],
    [-0.66544057, -1.63870215],
    [-1.00254705, -0.51569247],
    [ 3.17209481, -0.26278201],
    [ 0.48886797,  1.3654021 ],
    [-1.70863322, -1.02170044],
    [-0.06758577,  1.46233642],
    [-2.01185516, -1.27586457],
    [ 3.04203029, -1.25488069],
    [-0.92386867,  1.3693593 ]
])

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