#include <bits/stdc++.h>
using namespace std;

const double eps = 0.01;
const size_t max_iter = 1e6;
void conjugate_gradient(vector<vector<double>> &matr, vector<double> &b, vector<double> &x){
    int kl = 1;
    size_t n = b.size();

    double alpha, beta;
    double sumOfSquaresB = 0;
    double nominator, denominator;
    for(auto i : b)
        sumOfSquaresB += pow(i,2);

    x.resize(n, 0.2);
    vector<double> s(n), r(n), p(n);

    // r0 и p0
    for(size_t i = 0; i < n; ++i){
        s[i] = 0;
        // s[i] - умножение i-ого вектор-столбца матрицы на скаляр x[j]
        for(size_t j = 0; j < n; ++j)
            s[i] += matr[i][j] * x[j];
        r[i] = b[i] - s[i];
        p[i] = r[i];
    }

    size_t iter = 0;
    do {
        iter++;
        nominator = 0; // (r[i] , r[i]) -
        denominator = 0; // (Ar[i], r[i])

        for(size_t i = 0; i < n; ++i){
            s[i] = 0;
            for(size_t j = 0; j < n; ++j)
                s[i] += matr[i][j] * p[j];
            denominator += s[i] * p[i];
            nominator += r[i] * r[i];
        }
        alpha = nominator / denominator;
        denominator = 0;
        for(size_t i = 0; i < n; ++i){
            x[i] += alpha * p[i];
            r[i] -= alpha * s[i];
            denominator += r[i] * r[i];
        }
        swap(denominator,nominator);
        kl++;
        cout << kl << " ";
        for(auto i : x) cout << i << " ";
        cout << endl;
        beta = nominator / denominator;
        for(size_t i = 0; i < n; ++i)
            p[i] = r[i] + beta * p[i];

    } while(nominator / sumOfSquaresB > pow(eps,2) && iter < max_iter);
}

int main(){
    size_t n;
    cout << "Enter dim : ";
    cin >> n;
    cout << endl;

    cout << "Enter matrix A: " << endl;
    vector<vector<double>> matr(n,vector<double>(n));
    for(size_t i = 0; i < n; ++i)
        for(size_t j = 0; j < n; ++j)
            cin >> matr[i][j];

    cout << "Enter vector b : ";
    vector<double> b(n);
    for(size_t i = 0; i < n; ++i)
        cin >> b[i];
    cout << endl;

    vector<double> x(n);

    conjugate_gradient(matr, b, x);
    for(auto i : x)
        cout << fixed << setprecision(8) << i << " ";
    cout << endl;


}