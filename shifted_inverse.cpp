//
// Created by markd on 10.12.2019.
//
#include <iostream>
#include <fstream>
#include <cmath>
#include "matrix.cpp"
using namespace std;

const double EPS = 1e-15;

int main(){

    // Размерность матрицы
    int n;
    // Считываем размерность матрицы
    cin >> n;

    double a, c = 1, c_prev = 1.5;

    // Матрица системы и Q матрица в QR разложении
    Matrix A(n);

    cout.setf(ios::fixed);
    cout.precision(8);

    // Считываем A
    cin >> A;
    cin >> a;
    cout << endl;

    Vector x(n), y(n);
    Matrix B = (A - a * Matrix::E(n));
    // B = B.inverse();

    for (int i = 0; i < n; i++)
        x[i] = 1;

    double l = 0;
    int counter = 0;
    while (abs(c - c_prev) > EPS) {
        // y = (B * ~x).toVector();
        y = B.solve(x);
        c_prev = c;
        c = (y * x) / (x * x);
        x = y.normalized();
        counter++;
    }

    l = 1/c + a;
    cout << l << endl;
    cout << "Iteration count: " << counter;
    return 0;
}