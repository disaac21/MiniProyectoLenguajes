#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <cmath>
#include <algorithm>
#include "eigen-3.4.0/Eigen/Dense"
using namespace std;
vector<vector<double>> cargarArchivo(const string& filename) {
    vector<vector<double>> matriz;
    ifstream file(filename);
    if (file.is_open()) {
        string line;
        getline(file, line);//leer primera linea
        while (getline(file, line)) {
            stringstream ss(line);
            vector<double> fila;
            string valor;
            while (getline(ss, valor, ';')) {
                // Reemplazar las comas por puntos
                replace(valor.begin(), valor.end(), ',', '.');
                try {
                    fila.push_back(stod(valor));
                } catch (const invalid_argument& e) {
                    // Ignorar valores no numéricos
                }
            }
            matriz.push_back(fila);
        }
        file.close();
    } else {
        cerr << "No se pudo abrir el archivo: " << filename << endl;
    }
    return matriz;
}

// Función para calcular la media de cada columna de la matriz
vector<double> calcMean(const vector<vector<double>>& matriz) {
    int num_filas = matriz.size();
    int num_columnas = matriz[0].size();
    vector<double> medias(num_columnas, 0.0);
    for (int j = 0; j < num_columnas; ++j) {
        double suma_columna = 0.0;
        for (int i = 0; i < num_filas; ++i) {
            suma_columna += matriz[i][j];
        }
        double media_columna = suma_columna / num_filas;
        medias[j] = media_columna;
    }
    return medias;
}

// Función para calcular la desviación estándar de cada columna de la matriz
vector<double> calcStdDev(const vector<vector<double>>& matriz) {
    vector<double> meanColumnas = calcMean(matriz);
    int num_filas = matriz.size();
    int num_columnas = matriz[0].size();
    vector<double> stdDeviations(num_columnas, 0.0);
    for (int j = 0; j < num_columnas; ++j) {
        double sumSquares = 0.0;
        for (int i = 0; i < num_filas; ++i) {
            sumSquares += pow((matriz[i][j] - meanColumnas[j]), 2);
        }
        double columns_StdDev = sqrt(sumSquares / num_filas);
        stdDeviations[j] = columns_StdDev;
    }
    return stdDeviations;
}

// Función para centrar y reducir la matriz
vector<vector<double>> centrarReducir(const vector<vector<double>>& matriz) {
    vector<vector<double>> matrizCentradaReducida(matriz.size(), vector<double>(matriz[0].size(), 0.0));    
    
    vector<double> means = calcMean(matriz);
    vector<double> stdDeviations = calcStdDev(matriz);
    for (size_t i = 0; i < matriz.size(); ++i) {
        for (size_t j = 0; j < matriz[i].size(); ++j) {
            matrizCentradaReducida[i][j] = (matriz[i][j] - means[j]) / stdDeviations[j];
        }
    }
    return matrizCentradaReducida;
}
double calcCorrelaciones(const vector<double>& columnOne, const vector<double>& columnTwo) {
    double meansOne = 0.0;
    double meansTwo = 0.0;
    for (double val : columnOne) {
        meansOne += val;
    }
    meansOne /= columnOne.size();
    for (double val : columnTwo) {
        meansTwo += val;
    }
    meansTwo /= columnTwo.size();

    double numerador = 0.0;
    for (size_t i = 0; i < columnOne.size(); ++i) {
        numerador += (columnOne[i] - meansOne) * (columnTwo[i] - meansTwo);
    }

    double denominador = 0.0;
    for (double val : columnOne) {
        denominador += pow((val - meansOne), 2);
    }
    denominador = sqrt(denominador) * sqrt(columnTwo.size());

    double correlacionSend = numerador / denominador;
    return correlacionSend;
}
vector<double> getColumn(const vector<vector<double>>& matriz, int j) {
    vector<double> columna;
    for (size_t i = 0; i < matriz.size(); ++i) {
        columna.push_back(matriz[i][j]);
    }
    return columna;
}
vector<vector<double>> calcMatrizCorrelaciones(const vector<vector<double>>& matrizCentradaReducida) {
    int n = matrizCentradaReducida[0].size();
    vector<vector<double>> matrizCorrelaciones(n, vector<double>(n, 0.0));
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            double r = calcCorrelaciones(getColumn(matrizCentradaReducida, i), getColumn(matrizCentradaReducida, j));
            matrizCorrelaciones[i][j] = r;
        }
    }
    return matrizCorrelaciones;
}
vector<double> calcValoresPropios(const vector<vector<double>>& matrizCorrelaciones) {
    Eigen::MatrixXd matrizCorrelacionesEigen(matrizCorrelaciones.size(), matrizCorrelaciones[0].size());
    for (size_t i = 0; i < matrizCorrelaciones.size(); ++i) {
        for (size_t j = 0; j < matrizCorrelaciones[i].size(); ++j) {
            matrizCorrelacionesEigen(i, j) = matrizCorrelaciones[i][j];
        }
    }
    Eigen::EigenSolver<Eigen::MatrixXd> solver(matrizCorrelacionesEigen);
    Eigen::VectorXd valoresPropios = solver.eigenvalues().real();

    vector<double> valoresPropiosVec(valoresPropios.size());
    for (int i = 0; i < valoresPropios.size(); ++i) {
        valoresPropiosVec[i] = valoresPropios(i);
    }
    return valoresPropiosVec;
}
vector<vector<double>> calcVectoresPropios(const vector<vector<double>>& matrizCorrelaciones) {
    Eigen::MatrixXd matrizCorrelacionesEigen(matrizCorrelaciones.size(), matrizCorrelaciones[0].size());
    for (size_t i = 0; i < matrizCorrelaciones.size(); ++i) {
        for (size_t j = 0; j < matrizCorrelaciones[i].size(); ++j) {
            matrizCorrelacionesEigen(i, j) = matrizCorrelaciones[i][j];
        }
    }
    Eigen::EigenSolver<Eigen::MatrixXd> solver(matrizCorrelacionesEigen);
    Eigen::MatrixXd vectoresPropios = solver.eigenvectors().real();

    vector<vector<double>> vectoresPropiosVec(matrizCorrelaciones.size(), vector<double>(matrizCorrelaciones[0].size()));
    for (size_t i = 0; i < matrizCorrelaciones.size(); ++i) {
        for (size_t j = 0; j < matrizCorrelaciones[i].size(); ++j) {
            vectoresPropiosVec[i][j] = vectoresPropios(i, j);
        }
    }
    return vectoresPropiosVec;
}

vector<double> orderValoresPropios(vector<double> valoresPropios) {
    for (size_t i = 0; i < valoresPropios.size(); ++i) {
        for (size_t j = i + 1; j < valoresPropios.size(); ++j) {
            if (valoresPropios[i] < valoresPropios[j]) {
                swap(valoresPropios[i], valoresPropios[j]);
            }
        }
    }
    return valoresPropios;
}

vector<vector<double>> orderVectoresPropios(const vector<vector<double>>& vectoresPropios, const vector<double>& orderedValoresPropios) {
    vector<int> indexes(orderedValoresPropios.size());
    for (size_t i = 0; i < indexes.size(); ++i) {
        indexes[i] = i;
    }
    sort(indexes.begin(), indexes.end(), [&](int a, int b) {
        return orderedValoresPropios[a] > orderedValoresPropios[b];
    });

    vector<vector<double>> orderedVectoresPropios(vectoresPropios.size(), vector<double>(vectoresPropios[0].size()));
    for (size_t i = 0; i < vectoresPropios.size(); ++i) {
        for (size_t j = 0; j < vectoresPropios[i].size(); ++j) {
            orderedVectoresPropios[i][j] = vectoresPropios[i][indexes[j]];
        }
    }
    return orderedVectoresPropios;
}
vector<vector<double>> calcMatrizComponentesPrincipales(const vector<vector<double>>& centradaReducida, const vector<vector<double>>& matrizV) {
    int filas_CR = centradaReducida.size();
    int columnas_CR = centradaReducida[0].size();
    int filas_V = matrizV.size();
    int columnas_V = matrizV[0].size();

    vector<vector<double>> matrizComponentesPrincipalesTemp(filas_CR, vector<double>(columnas_V, 0.0));

    for (int i = 0; i < filas_CR; ++i) {
        for (int j = 0; j < columnas_V; ++j) {
            for (int k = 0; k < columnas_CR; ++k) {
                matrizComponentesPrincipalesTemp[i][j] += centradaReducida[i][k] * matrizV[k][j];
            }
        }
    }

    vector<vector<double>> matrizComponentesPrincipales(filas_CR, vector<double>(2, 0.0));
    for (int i = 0; i < filas_CR; ++i) {
        for (int j = 0; j < 2; ++j) {
            matrizComponentesPrincipales[i][j] = matrizComponentesPrincipalesTemp[i][j];
        }
    }
    return matrizComponentesPrincipales;
}
vector<vector<double>> calcMatrizCalidadesIndividuos(const vector<vector<double>>& matrizComponentesPrincipales, const vector<vector<double>>& matrizCentradaReducida) {
    vector<vector<double>> listaMatrizCalidadesIndividuos;
    int rowCount = matrizComponentesPrincipales.size();
    int columnCount = matrizComponentesPrincipales[0].size();
    
    for (int i = 0; i < rowCount; ++i) {
        vector<double> fila;
        for (int j = 0; j < columnCount; ++j) {
            double num = pow(matrizComponentesPrincipales[i][j], 2);
            double den = 0.0;
            for (int k = 0; k < columnCount; ++k) {
                den += pow(matrizCentradaReducida[i][k], 2);
            }
            fila.push_back(num / den);
        }
        listaMatrizCalidadesIndividuos.push_back(fila);
    }
    
    return listaMatrizCalidadesIndividuos;
}
vector<vector<double>> calcMatrizCoordenadasVariables(const vector<vector<double>>& centradaReducida, const vector<vector<double>>& vectoresPropios) {
    int rowCount = centradaReducida.size();
    int columnCount = vectoresPropios[0].size();

    Eigen::MatrixXd centradaReducidaEigen(rowCount, centradaReducida[0].size());
    for (int i = 0; i < rowCount; ++i) {
        for (int j = 0; j < centradaReducida[i].size(); ++j) {
            centradaReducidaEigen(i, j) = centradaReducida[i][j];
        }
    }

    Eigen::MatrixXd vectoresPropiosEigen(vectoresPropios.size(), vectoresPropios[0].size());
    for (int i = 0; i < vectoresPropios.size(); ++i) {
        for (int j = 0; j < vectoresPropios[i].size(); ++j) {
            vectoresPropiosEigen(i, j) = vectoresPropios[i][j];
        }
    }

    Eigen::MatrixXd matrizCoordenadasVariablesEigen = centradaReducidaEigen * vectoresPropiosEigen.transpose();

    vector<vector<double>> matrizCoordenadasVariables(rowCount, vector<double>(columnCount));
    for (int i = 0; i < rowCount; ++i) {
        for (int j = 0; j < columnCount; ++j) {
            matrizCoordenadasVariables[i][j] = matrizCoordenadasVariablesEigen(i, j);
        }
    }

    return matrizCoordenadasVariables;
}
vector<double> calcMatrizCalidadesVariables(const vector<double>& valoresPropios, int entries) {
    vector<double> matrizCalidadesVariables(valoresPropios.size());
    for (size_t i = 0; i < valoresPropios.size(); ++i) {
        matrizCalidadesVariables[i] = valoresPropios[i] / entries;
    }
    return matrizCalidadesVariables;
}
vector<double> calcVectorInerciasEjes(const vector<double>& valoresPropios) {
    vector<double> vectorInerciasEjes;
    int den = valoresPropios.size();
    for (size_t i = 0; i < valoresPropios.size(); ++i) {
        double num = valoresPropios[i] * 100.0;
        double value = num / den;
        vectorInerciasEjes.push_back(value);
    }
    return vectorInerciasEjes;
}
// Función para imprimir un vector
void imprimirVector(const vector<double>& vector) {
    cout << "[";
    for (const auto& elemento : vector) {
        cout << "["<< elemento << "] ";
    }
    cout << "]" << endl;
    cout << endl;
}
void imprimirMatriz(const vector<vector<double>>& matriz) {
    cout << "[";
    for (const auto& fila : matriz) {
        for (const auto& elemento : fila) {
            cout << "["<< elemento << "] ";
        }
        cout << endl;
    }
    cout << "]" << endl;
}
void imprimirParesValores(const vector<double>& valoresPropios, const vector<std::vector<double>>& vectoresPropios) {
    for (size_t i = 0; i < valoresPropios.size(); ++i) {
        cout << "Valor propio: " << valoresPropios[i] << std::endl;
        cout << "Vector propio: ";
        for (size_t j = 0; j < vectoresPropios.size(); ++j) {
            cout << vectoresPropios[j][i] << " ";
        }
        cout << endl << endl;
    }
}
int main() {
    // Nombre del archivo
    const string filename = "EjemploEstudiantes.csv";

    // Cargar archivo y crear matriz
    vector<vector<double>> matrizRaw = cargarArchivo(filename);
    
    // Centrar y reducir
    vector<vector<double>> matrizCentradaReducida = centrarReducir(matrizRaw);
    
    // Cálculo de matriz de correlaciones
    vector<vector<double>> matrizCorrelaciones = calcMatrizCorrelaciones(matrizCentradaReducida);
    
    // Cálculo de valores y vectores propios
    vector<double> valoresPropios = calcValoresPropios(matrizCorrelaciones);
    vector<vector<double>> vectoresPropios = calcVectoresPropios(matrizCorrelaciones);
    
    // Ordenamiento de valores y vectores propios
    vector<double> ordenadosValoresPropios = orderValoresPropios(valoresPropios);
    vector<vector<double>> ordenadosVectoresPropios = orderVectoresPropios(vectoresPropios, ordenadosValoresPropios);
    
    // Construcción de matriz V
    vector<vector<double>> matrizV = ordenadosVectoresPropios;
    
    // Cálculo de matriz de componentes principales
    vector<vector<double>> matrizComponentesPrincipales = calcMatrizComponentesPrincipales(matrizCentradaReducida, matrizV);
    
    // Cálculo de matriz de calidades de individuos
    vector<vector<double>> matrizCalidadesIndividuos = calcMatrizCalidadesIndividuos(matrizComponentesPrincipales, matrizCentradaReducida);
    
    // Cálculo de matriz de coordenadas de las variables
    vector<vector<double>> matrizCoordenadasVariables = calcMatrizCoordenadasVariables(matrizCentradaReducida, ordenadosVectoresPropios);
    
    // Cálculo de matriz de calidades de las variables
    vector<double> matrizCalidadesVariables = calcMatrizCalidadesVariables(ordenadosValoresPropios, matrizCentradaReducida.size());
    
    // Cálculo de vector de inercias de los ejes
    vector<double> vectorInerciasEjes = calcVectorInerciasEjes(ordenadosValoresPropios);


    // Salida de resultados
    cout << "\nPaso 1 - Matriz Centrada y Reducida:\n";
    imprimirMatriz(matrizCentradaReducida);
    cout << "\nPaso 2 - Matriz de Correlaciones:\n";
    imprimirMatriz(matrizCorrelaciones);
    cout << "\nPaso 3 - Valores y Vectores Propios:\n";
    imprimirParesValores(ordenadosValoresPropios, ordenadosVectoresPropios);
    cout << "Matriz de Vectores Propios:\n";
    imprimirMatriz(ordenadosVectoresPropios);
    cout << "Matriz de Valores Propios:\n";
    imprimirVector(ordenadosValoresPropios);
    cout << "\nPaso 4 - Matriz V:\n";
    imprimirMatriz(matrizV);
    cout << "\nPaso 5 - Matriz de Componentes Principales:\n";
    imprimirMatriz(matrizComponentesPrincipales);
    cout << "\nPaso 6 - Matriz de Calidades de Individuos:\n";
    imprimirMatriz(matrizCalidadesIndividuos);
    cout << "\nPaso 7 - Matriz de Coordenadas de las Variables:\n";
    imprimirMatriz(matrizCoordenadasVariables);
    cout << "\nPaso 8 - Matriz de Calidades de las Variables:\n";
    imprimirVector(matrizCalidadesVariables);
    cout << "\nPaso 9 - Vector de Inercias de los Ejes:\n";
    imprimirVector(vectorInerciasEjes);

    return 0;
}
