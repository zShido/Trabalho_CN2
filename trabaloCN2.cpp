#include <iostream>
#include <vector>
#include <cmath>

using namespace std;

vector<double> jacobi(const vector<vector<double>>& matrizA, const vector<double>& vetorB, const vector<double>& x0, int maxIteracoes, double tolerancia, int& iteracoes) {
    int n = matrizA.size();
    vector<double> x(x0);
    vector<double> xNovo(n, 0.0);

    for (iteracoes = 0; iteracoes < maxIteracoes; ++iteracoes) {
        for (int i = 0; i < n; ++i) {
            double soma = 0.0;
            for (int j = 0; j < n; ++j) {
                if (i != j) {
                    soma += matrizA[i][j] * x[j];
                }
            }
            xNovo[i] = (vetorB[i] - soma) / matrizA[i][i];
        }

        double erro = 0.0;
        for (int i = 0; i < n; ++i) {
            erro += abs(xNovo[i] - x[i]);
        }

        if (erro < tolerancia) {
            break;
        }

        x = xNovo;
    }

    return x;
}

vector<double> gaussSeidel(const vector<vector<double>>& matrizA, const vector<double>& vetorB, const vector<double>& x0, int maxIteracoes, double tolerancia, int& iteracoes) {
    int n = matrizA.size();
    vector<double> x(x0);

    for (iteracoes = 0; iteracoes < maxIteracoes; ++iteracoes) {
        bool convergiu = true;
        for (int i = 0; i < n; ++i) {
            double soma = 0.0;
            for (int j = 0; j < n; ++j) {
                if (i != j) {
                    soma += matrizA[i][j] * x[j];
                }
            }
            double xNovo = (vetorB[i] - soma) / matrizA[i][i];
            if (abs(xNovo - x[i]) > tolerancia) {
                convergiu = false;
            }
            x[i] = xNovo;
        }

        if (convergiu) {
            break;
        }
    }

    return x;
}

vector<double> gaussEliminacao(vector<vector<double>>& matrizA, vector<double>& vetorB) {
    int n = matrizA.size();
    vector<double> x(n);

    // Eliminação para frente
    for (int i = 0; i < n; ++i) {
        for (int k = i + 1; k < n; ++k) {
            double fator = matrizA[k][i] / matrizA[i][i];
            for (int j = i; j < n; ++j) {
                matrizA[k][j] -= fator * matrizA[i][j];
            }
            vetorB[k] -= fator * vetorB[i];
        }
    }

    // Substituição para trás
    for (int i = n - 1; i >= 0; --i) {
        x[i] = vetorB[i];
        for (int j = i + 1; j < n; ++j) {
            x[i] -= matrizA[i][j] * x[j];
        }
        x[i] /= matrizA[i][i];
    }

    return x;
}

int main() {
    int tamanho;
    cout << "Insira o tamanho da matriz (ex: 3 para uma matriz 3x3): ";
    cin >> tamanho;

    vector<vector<double>> matrizA(tamanho, vector<double>(tamanho, 0.0));
    vector<double> vetorB(tamanho, 0.0);
    vector<double> x(tamanho, 0.0);
    vector<double> x0 = {0, 0, 0}; // Vetor inicial

    int maxIteracoes;
    double tolerancia;

    cout << "Insira os elementos da matriz A (" << tamanho << "x" << tamanho << "):\n";
    for (int i = 0; i < tamanho; ++i) {
        for (int j = 0; j < tamanho; ++j) {
            cin >> matrizA[i][j];
        }
    }

    cout << "Insira os elementos do vetor B (" << tamanho << " elementos):\n";
    for (int i = 0; i < tamanho; ++i) {
        cin >> vetorB[i];
    }

    cout << "Insira o valor de maxIteracoes: ";
    cin >> maxIteracoes;

    cout << "Insira o valor de tolerancia: ";
    cin >> tolerancia;

    int iteracoes;
    cout << "Escolha o método para resolver o sistema de equações:\n";
    cout << "1. Jacobi\n";
    cout << "2. Gauss-Seidel\n";
    cout << "3. Eliminação de Gauss\n";
    cout << "Opção: ";
    int opcao;
    cin >> opcao;

    switch (opcao) {
        case 1:
            x = jacobi(matrizA, vetorB, x0, maxIteracoes, tolerancia, iteracoes);
            break;
        case 2:
            x = gaussSeidel(matrizA, vetorB, x0, maxIteracoes, tolerancia, iteracoes);
            break;
        case 3:
            x = gaussEliminacao(matrizA, vetorB);
            break;
        default:
            cout << "Opção inválida.\n";
            return 1;
    }

    cout << "Solução do sistema de equações:\n";
    for (int i = 0; i < tamanho; ++i) {
        cout << "x[" << i << "] = " << x[i] << "\n";
    }

    cout << "Número de iterações: " << iteracoes << endl;

    return 0;
}
