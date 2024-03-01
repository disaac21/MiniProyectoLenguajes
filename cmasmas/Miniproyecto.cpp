#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include "eigen-3.4.0/Eigen/Dense"
#include "eigen-3.4.0/Eigen/Eigenvalues"

using Eigen::Matrix3d;
using namespace std;

void imprimirMatriz(const double *matriz, int filas, int columnas)
{

    cout << endl;
    for (int i = 0; i < filas; ++i)
    {
        for (int j = 0; j < columnas; ++j)
        {
            cout << matriz[i * columnas + j] << " ";
        }
        cout << endl;
    }
}

int main()
{
    cout << "Miniproyecto de Estadística" << endl;
    // Primera lectura para obtener cantidad de columnas y de filas
    string archivoCSV = "EjemploEstudiantes.csv", linea = "";
    ifstream archivo(archivoCSV);

    // Leemos la primera línea (encabezado)
    getline(archivo, linea);
    stringstream stream(linea);

    string variable = "";
    char delimitador = ';';
    int columnas = 0, filas = 1;

    getline(archivo, linea);

    // Obtener la cantidad de variables (columnas)
    // cout << linea << endl;
    for (char c : linea)
    {
        if (c == ';')
        {
            columnas++;
        }
    }

    // Obtener la cantidad de elementos (filas)
    while (getline(archivo, linea))
        filas++;
    // while (getline(stream, variable, delimitador))
    //     columnas++;

    // Segunda lectura para obtener los datos
    archivo.clear();
    archivo.seekg(0, ios::beg);

    // Omitir la primera línea (encabezado)
    getline(archivo, linea);

    // Llenar la matriz de datos originales
    double matriz[filas][columnas], media[columnas], desvEstandar[columnas];
    for (int i = 0; i < filas; i++)
    {
        getline(archivo, linea);
        stringstream filaStream(linea);
        string valor = "";

        // Omitir la primera columna (nombre del alumno)
        getline(filaStream, valor, delimitador);

        // cout << "columnas: " << columnas << endl;

        for (int j = 0; j < columnas; j++)
        {
            getline(filaStream, valor, delimitador);
            for (char &c : valor)
            {
                // Replace ',' with '.'
                if (c == ',')
                {
                    c = '.';
                }
                matriz[i][j] = stod(valor);
                media[j] += stod(valor);
            }
            cout << matriz[i][j] << " ";
        }
        cout << endl;
    }

    // Calcular media y desv. estandar de cada variable
    // y obtener la transpuesta de la matriz original
    double transpuesta[columnas][filas];
    cout << endl;
    cout << "Transpuesta: " << endl;
    for (int j = 0; j < columnas; j++)
    {
        media[j] /= filas;
        for (int i = 0; i < filas; i++)
        {
            desvEstandar[j] += pow((matriz[j][i] - media[j]), 2);
            transpuesta[i][j] = matriz[i][j];
            cout << transpuesta[i][j] << " ";
        }
        desvEstandar[j] /= filas;
        cout << endl;
    }
    // imprimirMatriz(&transpuesta[0][0], filas, columnas);

    cout << endl
         << "Normalizada: " << endl;
    // Obtener matriz normalizada (centrada y reducida)
    Eigen::MatrixXd normalizada(filas, columnas);
    for (int i = 0; i < filas; i++)
    {
        for (int j = 0; j < columnas; j++)
        {
            normalizada(i, j) = (matriz[i][j] - media[i]) / desvEstandar[i];
            cout << normalizada(i, j) << " ";
        }
        cout << endl;
    }

    // Obtener matriz de correlación
    Eigen::MatrixXd correlacion(filas, filas);
    for (int i = 0; i < filas; ++i)
    {
        for (int j = 0; j < filas; ++j)
        {
            correlacion(i, j) = 0.0;
            for (int k = 0; k < columnas; ++k)
                correlacion(i, j) += matriz[i][k] * transpuesta[j][k];
            correlacion(i, j) /= filas;
        }
    }

    // Realizar cálculos con Eigen
    Eigen::EigenSolver<Eigen::MatrixXd> solver(correlacion);
    Eigen::MatrixXd eigenvectores = solver.eigenvectors().real();
    Eigen::VectorXd eigenvalues = solver.eigenvalues().real();

    // Obtener matriz de componentes principales (normalizada x eigenvetores)
    // Eigen::MatrixXd compPrincipales = normalizada * eigenvectores;

    archivo.close();
    return 0;
}