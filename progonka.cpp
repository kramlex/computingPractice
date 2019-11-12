//
// Created by markd on 07.10.2019.
//

#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>

using namespace std;

int main(){
    // Размерность матрицы
    int n;
    //
    cout << "Enter dim A: ";
    cin >> n;
    // Сама матрица
    vector<vector<double>> matr(n , vector<double>(n));

    // Вектор b из уравнения Ax = b , где A - трехдиагональная матрица
    vector<double> vec(n);

    // Массив a (главная диагональ) , массив b(над главной диагональю) , массив c(под главной диагональю)
    vector<double> b(n) , a(n) , c(n-1);
    a[0] = 0;
    // Считываем матрицу
    cout << "Enter matrix A(dim=" << n << ") :" << endl;
    for(int i = 0; i < n ; i++){
        for(int j = 0; j < n; j++){
            cin >> matr[i][j];
            if(i == j)
                b[i] = matr[i][j];
            if(j == i + 1)
                c[i] = matr[i][j];
            if(j == i - 1)
                a[i] = matr[i][j];
        }
    }

    // Считываем ветор b
    cout << "Enter b(string): ";
    for(int i = 0; i < n ; i++){
        cin >> vec[i];
    }
    cout << endl;

    // Прямой ход алгоритма
    // Прогоночные коэфициенты и переменная нужная для вычисления прогоночных коэф.
    vector<double> alpha(n) , beta(n);
    double y;
    y = b[0];
    alpha[0] = -c[0] / y;
    beta[0] = vec[0] / y;
    cout << y << " " << alpha[0] << " " << beta[0] << endl;
    for(int i = 1; i < n - 1 ; i++){
        y = b[i] + a[i] * alpha[i-1];
        alpha[i] = -c[i] / y;
        beta[i] = (vec[i] - a[i] * beta[i-1]) / y;
    }
    y = b[n-1] + a[n-1] * alpha[n-2];
    beta[n-1] = (vec[n-1] - a[n-1] * beta[n-2]) / y;


    // Обратный ход алгоритма
    vector<double> res(n);
    res[n-1] = beta[n-1];
    for(int i = n - 2; i >= 0 ; i--){
        res[i] = alpha[i] * res[i+1] + beta[i];
    }

    cout << endl;
    for (int i = 0; i < n; i++) {
        cout << "x" << (i + 1) << " = " << fixed << setprecision(8) << res[i] << endl;
    }

}



