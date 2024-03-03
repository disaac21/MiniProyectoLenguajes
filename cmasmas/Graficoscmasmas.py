import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# Leer los datos del archivo CSV
datos = pd.read_csv('datos.csv')

# Extraer las coordenadas de las variables y los vectores
variables = datos[datos['Variable'].str.startswith('V')]
vectores = datos[datos['Variable'].str.startswith('Vector')]

# Crear el gráfico del plano principal
plt.figure(figsize=(8, 6))
plt.scatter(variables['X'], variables['Y'], color='blue', label='Variables')
plt.quiver(0, 0, vectores['X'], vectores['Y'], angles='xy', scale_units='xy', scale=1, color='red', label='Vectores Propios')
plt.xlabel('T1')
plt.ylabel('T2')
plt.title('Plano Principal')
plt.axhline(0, color='black',linewidth=0.5)
plt.axvline(0, color='black',linewidth=0.5)
plt.grid(color = 'gray', linestyle = '--', linewidth = 0.5)
plt.legend()
plt.show()

# Crear el gráfico del círculo de correlación
plt.figure(figsize=(6, 6))
plt.scatter(vectores['X'], vectores['Y'], color='red', label='Vectores Propios')
plt.xlabel('C1')
plt.ylabel('C2')
plt.title('Círculo de Correlación')
plt.axhline(0, color='black',linewidth=0.5)
plt.axvline(0, color='black',linewidth=0.5)
plt.grid(color = 'gray', linestyle = '--', linewidth = 0.5)
plt.gca().set_aspect('equal', adjustable='box')
plt.legend()
plt.show()