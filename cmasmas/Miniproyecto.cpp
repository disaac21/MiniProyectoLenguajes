#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <iomanip>
#include "eigen-3.4.0/Eigen/Dense"
#include "eigen-3.4.0/Eigen/Eigenvalues"
using Eigen::MatrixXd;
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
        cout << "eigenvectores: " << eigenvectores(i) << endl;
    }

    for (size_t i = 0; i < eigenvalues.size(); i++)
    {
        cout << "eigenvalues: " << eigenvalues(i) << endl;
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
            cout << setw(5) << matriz[i][j];
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
        cout << endl
             << medianas[i] << " ";
    }



    archivo.close();
    return 0;
}