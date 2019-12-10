//
// Created by markd on 10.12.2019.
//

#include <bits/stdc++.h>
#include "matrix.cpp"
using namespace std;

const double EPS = 1e-9;

int main(){
    int n;
    cin >> n; // считываем размерность

    Matrix A(n); // матрица системы
    Vector b(n); // вектор b в уравнении Ax=b

    double l_min, l_max;
    cout.setf(ios::fixed);
    cout.precision(8);

    // Считываем A
    cout << "Enter matrix A: " << endl;
    cin >> A;
    cout << endl;
    // Считываем вектор b
    cout << "Enter vector b: " << endl;
    cin >> b;

    cin >> l_min >> l_max;

    Vector x(n),x_n(n);

    double w0 = (l_max - 1) / (l_max + l_min);
    double w = w0, w_n = - 1/(2 * (l_max + l_min)/(l_max - l_min) + w);

    int k;

    while (((A*~x_n).toVector() - b).length() > EPS) {
        Vector t = x_n + w * w_n * (x_n - x) - 2 / (l_min + l_max) * (1 + w * w_n) * (A*~x_n - ~b).toVector();

        x = x_n;
        x_n = t;

        w = w_n;
        w_n = - 1/(2 * (l_max + l_min)/(l_max - l_min) + w);
        k++;
    }

    cout << "Iteration count : " << k << endl;

    cout << x_n;

}

