//
// Created by markd on 10.12.2019.
//
# include <iostream>
# include <cmath>
# include "matrix.cpp"
//# define EPS 0.00000001
using namespace std;

const double EPS = 1e-8;

double max_above_d(Matrix& a) {
    double s = 0;
    for (int i = 0; i < a.rows_count(); ++i) {
        for (int j = i + 1; j < a.rows_count(); ++j) {
            s += pow(a[j][i],2);
        }
    }
    return sqrt(s);
}

int main() {
    // Размерность матрицы
    int n;
    // Считываем размерность матрицы
    cin >> n;
    // Матрица системы и Q матрица в RQ разложении
    Matrix A(n);

    cout.setf(ios::fixed);
    cout.precision(8);

    // Считываем A
    cout << "Enter matrix A: " << endl;
    cin >> A;
    cout << endl;

    int i;

    for (i = 0; max_above_d(A) > EPS &&  i < 20; i++) {
//        A = qr(A);
        auto qrs = A.qr();
        A = qrs.second * qrs.first;
    }

    cout << A;
    cout << endl << n << endl;
    for (int k = 0; k < n; k++) {
        cout << "lamda" << k + 1 << " = " << A[k][k] << endl;
    }


    // (QR)x = b сводится к задаче Rx = ~Qx = b'
    cout << "Iteration count : " << i << endl;

    return 0;
}
