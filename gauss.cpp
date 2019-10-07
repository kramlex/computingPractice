//
// Created by markd on 10.09.2019.
//

#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>

using namespace std;


void print_matrix(vector<vector<double> > matrix, vector<double> vec) {
    cout << endl;
    for (int i = 0; i < matrix.size(); i++) {
        for (int j = 0; j < matrix.size(); j++) {
            cout << matrix[i][j] << " ";
        }
        cout << " | ";
        cout << fixed << setprecision(8) << vec[i];
        cout << endl;
    }
}

void print_matrix(vector<vector<double> > matrix) {
    cout << endl;
    for(auto i : matrix){
        for(auto j : i)
            cout << j << " ";
        cout << endl;
    }
    cout << endl;
}

int main() {
    // Размерность матрицы
    int n;
    //
    cout << "Enter dim A : ";
    cin >> n;
    // Сама матрица, ее определитель, и обратная.
    vector<vector<double>> matr(n , vector<double>(n)) , inv(n, vector<double>(n));
    double det = 1.;

    // Вектор b (из уравнения Ax=b)
    vector<double> vec(n);

    cout << endl << "Enter matrix A: " << endl;
    // Считываем матрицу
    for(int i = 0; i < n ; i++){
        for(int j = 0; j < n; j++){
            cin >> matr[i][j];
            inv[i][j] = i == j ? 1 : 0;
        }
    }
    cout << endl << "Enter vector b (string) : ";
    // Считываем вектор b
    for(int i = 0 ; i < n ; i++){
        cin >> vec[i];
    }
    cout << endl;

    // Цикл приведения к треугольному виду
    for (int i = 0; i < n; i++) {
        // Фиксируем диагональный элемент
        double elem = matr[i][i];

        if (elem == 0.) {
            int j = 0;
            // Ищем такую строку из последующих
            // i-ый элемент которой не нулевой
            for (j = i + 1; matr[j][i] == 0 && j < n; j++);
            // Поэлементно прибавляем к i-ой строке j-ую
            for (int k = 0; k < n; k++) {
                matr[i][k] += matr[j][k];
                inv[i][k] += inv[j][k];
            }
            vec[i] += vec[j];
            // Фиксируем новый диагональный элемент
            elem = matr[i][i];
        }

        // По всем строкам ниже i-ой
        for (int j = i + 1; j < n; j++) {
            // Коэффициент, при умножении на который i-ой строки
            //  И прибавлении её к j-ой мы занулим i-ый элемент j-ой строки
            double coef = - matr[j][i] / elem;
            for (int l = 0; l < n; l++) {
                matr[j][l] += matr[i][l] * coef;
                inv[j][l] += inv[i][l] * coef;
            }
            vec[j] += vec[i] * coef;
        }

        // Делим i-ую строку на диагональный
        for (int j = 0; j < n; j++) {
            matr[i][j] /= elem;
            inv[i][j] /= elem;
        }
        vec[i] /= elem;
        det *= elem;
    }

    print_matrix(matr, vec);
    // Считаем решенгия
    for (int i = n - 1; i >= 0; i--) {
        for (int j = i + 1; j < n; j++) {
            vec[i] -= matr[i][j] * vec[j];
        }
    }

    // Доводим процесс поиска обратной до конца
    for (int i = n - 1; i >= 0; i--) {
        for (int j = i - 1; j >= 0; j--) {
            for (int k = 0; k < n; k++)
                inv[j][k] -= inv[i][k] * matr[j][k];
        }
    }

    cout << endl;

    for (int i = 0; i < n; i++) {
        cout << "x" << (i + 1) << " = " << setprecision(8) << vec[i] << endl;
    }

    cout << endl << "Det = " << fixed <<setprecision(8) << det;
    cout << endl << endl << "Inverse matrix" << endl;
    print_matrix(inv);

}