#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <iomanip>
#include "eigen-3.4.0/Eigen/Dense"
#include "eigen-3.4.0/Eigen/Eigenvalues"
using Eigen::MatrixXd;
using namespace std;

void imprimirMatriz(double *matriz, int filas, int columnas)
{
    for (int i = 0; i < filas; i++)
    {
        for (int j = 0; j < columnas; j++)
        {
            cout << setw(7) << matriz[i * columnas + j];
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
    cout << endl
         << "matriz base " << endl;
    double matriz[filas][columnas], media[columnas], desvEstandar[columnas], medianas[columnas];
    for (int i = 0; i < filas; i++)
    {
        getline(archivo, linea); // agarra toda la linea
        stringstream filaStream(linea);
        string valor = "";

        // Omitir la primera columna (nombre del alumno)
        getline(filaStream, valor, delimitador); // lee la linea que agarramos antes

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
            }
            try
            {
                matriz[i][j] = stod(valor);
                media[j] += stod(valor);
                // cout << endl << "media[" << j << "]: " << media[j] << " ";
                cout << setw(5) << matriz[i][j];
            }
            catch (const std::invalid_argument &e)
            {
                // Handle invalid input
                cout << "Invalid input: " << valor << endl;
            }
        }
        cout << endl;
    }

    cout << endl
         << "imprimiendo con funcion" << endl;
    imprimirMatriz((double *)matriz, filas, columnas);

    // cout << endl << "Media[0]: " << media[0] << endl;

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
    // imprimirMatriz((double *)transpuesta, columnas, filas);

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
    cout << endl
         << "matriz de correlacion" << endl;
    Eigen::MatrixXd correlacion(filas, filas);
    for (int i = 0; i < filas; ++i)
    {
        for (int j = 0; j < filas; ++j)
        {
            correlacion(i, j) = 0.0;
            for (int k = 0; k < columnas; ++k)
                correlacion(i, j) += matriz[i][k] * transpuesta[j][k];
            correlacion(i, j) /= filas;
            cout << setw(9) << correlacion(i, j);
        }
        cout << endl;
    }

    // Realizar cálculos con Eigen
    Eigen::EigenSolver<Eigen::MatrixXd> solver(correlacion);
    Eigen::MatrixXd eigenvectores = solver.eigenvectors().real();
    Eigen::VectorXd eigenvalues = solver.eigenvalues().real();

    for (size_t i = 0; i < eigenvectores.size(); i++)
    {
        // cout << "eigenvectores: " << eigenvectores(i) << endl;
    }

    for (size_t i = 0; i < eigenvalues.size(); i++)
    {
        // cout << "eigenvalues: " << eigenvalues(i) << endl;
    }

    // Obtener matriz de componentes principales (normalizada x eigenvetores)
    // Eigen::MatrixXd compPrincipales = normalizada * eigenvectores;

    //----- relleno ----
    cout << "medias" << endl;
    for (int i = 0; i < columnas; i++)
    {
        cout << media[i] << " ";
    }
    cout << endl;

    cout << endl
         << "centrada" << endl;
    double centrada[filas][columnas];
    for (int i = 0; i < filas; i++)
    {
        for (int j = 0; j < columnas; j++)
        {
            centrada[i][j] = matriz[i][j] - media[j];
            // cout << centrada[i][j] << " ";
            cout << setw(5) << matriz[i][j]; // aca esta el error
        }
        cout << endl;
    }

    cout << "desvEstandar" << endl;
    for (int i = 0; i < columnas; i++)
    {
        cout << desvEstandar[i] << " ";
    }
    cout << endl
         << "reducida con desviacion estandar" << endl;
    double reducida[filas][columnas];
    for (int i = 0; i < filas; i++)
    {
        for (int j = 0; j < columnas; j++)
        {
            reducida[i][j] = centrada[i][j] / desvEstandar[j];
            cout << setw(11) << reducida[i][j] << " ";
            // cout << centrada[i][j] << " ";
            // cout << setw(5) << matriz[i][j];
        }
        cout << endl;
    }

    //-- SACANDO LA VARIANZA --
    int mediana = filas / 2;
    mediana = round(mediana);
    cout << endl
         << "mediana: " << mediana << endl;

    double temp[filas];

    for (size_t i = 0; i < columnas; i++)
    {
        for (size_t j = 0; j < filas; j++)
        {
            temp[j] = matriz[j][i];
        }
        sort(temp, temp + filas);
        medianas[i] = temp[mediana];

        cout << "temp" << endl;
        for (size_t i = 0; i < filas; i++)
        {
            cout << temp[i] << " ";
        }
        cout << endl;
    }

    cout << endl
         << "medianas" << endl;
    for (size_t i = 0; i < columnas; i++)
    {
        cout
            << medianas[i] << endl;
    }

    Eigen::MatrixXd V(columnas, columnas); // Crear matriz V de tamaño m x m

    // Llenar la matriz V con los vectores propios ordenados como columnas
    for (int i = 0; i < columnas; ++i)
    {
        for (int j = 0; j < columnas; ++j)
        {
            V(j, i) = eigenvectores(j, i); // Seleccionar el i-ésimo vector propio y colocarlo como la columna i de V
        }
    }

    // Imprimir la matriz V
    cout << "Paso 4:" << endl;
    cout << V << endl;

    cout << "Paso 5:" << endl;
    Eigen::MatrixXd C(filas, columnas);

    for (int i = 0; i < filas; ++i)
    {
        for (int j = 0; j < columnas; ++j)
        {
            C(i, j) = 0;
            for (int k = 0; k < columnas; ++k)
            {
                C(i, j) += normalizada(i, k) * V(k, j);
            }
        }
    }
    cout << C << endl;

    cout << "Paso 6:" << endl;
    cout << "Matriz de calidades de individuos Q:" << endl;
    Eigen::MatrixXd Q(filas, columnas); // Crear una matriz para almacenar el resultado
    for (int i = 0; i < filas; ++i)
    {
        for (int j = 0; j < columnas; ++j)
        {
            Q(i, j) = 0;
            for (int k = 0; k < columnas; ++k)
            {
                Q(i, j) += normalizada(i, k) * V(k, j); // Multiplicación de matrices manual
            }
        }
    }
    cout << Q << endl;

    // Paso 7: Calcular la matriz de coordenadas de las variables T
    cout << "Paso 7 :" << endl;
    cout << "Matriz de coordenadas de las variables T :" << endl;
    Eigen::MatrixXd T(filas, columnas); // Crear una matriz para almacenar el resultado

    // Multiplicación de la matriz normalizada por los vectores propios
    for (int i = 0; i < filas; ++i)
    {
        for (int j = 0; j < columnas; ++j)
        {
            T(i, j) = 0;
            for (int k = 0; k < columnas; ++k)
            {
                T(i, j) += normalizada(i, k) * V(k, j); // Calcular la matriz de coordenadas de las variables T
            }
        }
    }
    cout << T << endl;
    cout << "Paso 8:" << endl;
    cout << "Matriz de calidades de las variables S:" << endl;
    Eigen::MatrixXd S(columnas, columnas);

    for (int i = 0; i < columnas; ++i)
    {
        for (int j = 0; j < columnas; ++j)
        {
            S(i, j) = 0;
            for (int k = 0; k < filas; ++k)
            {
                S(i, j) += pow(C(k, i), 2);
            }
            S(i, j) /= filas;
            cout << setw(9) << S(i, j);
        }
        cout << endl;
    }
    cout << "Vector de inercias de los ejes I:" << endl;
    Eigen::MatrixXd I(1, columnas);

    for (int j = 0; j < columnas; ++j)
    {
        double suma = 0;
        for (int i = 0; i < columnas; ++i)
        {
            suma += S(i, j);
        }
        I(0, j) = suma;
        cout << setw(9) << I(0, j);
    }
    cout << endl;

    archivo.close();
    return 0;
}