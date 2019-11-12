#include <iostream>
#include <cassert>
#include <climits>
#include <cmath>
#include <iomanip>

using namespace std;

class matrix {
    class vector {
    public:
        vector() : len(0), value(nullptr){};
        explicit vector(size_t n) : len(n) , value(new double[n]) {}
        double &operator[](size_t index) const {
            return value[index];
        }
        size_t len;
        double *value;
    };

    void init(const size_t n, const size_t m) {
        row = n; col = m;
        arrvec = new vector[row];
        for (size_t i = 0; i < row; i++)
            arrvec[i] = vector(col);
    }

public:
    matrix(size_t n, size_t m) : row(n), col(m), arrvec(nullptr){
        init(n, m);
    }

    matrix(const matrix &that, const size_t i1, const size_t j1, const size_t i2, const size_t j2, char key)
    : row(0), col(0), arrvec(nullptr) {
        if (key == 'a') {
            init(i2 - i1 + 1, j2 - j1 + 1);
            for (size_t i = i1; i <= i2; i++)
                for (size_t j = j1; j <= j2; j++)
                    (*this)[i - i1][j - j1] = that[i][j];
        }
        if (key == 'b') {
            init(i1, j1);
            this->updateValue(0);
            for (size_t i = i2; i < that.row + i2; i++)
                for (size_t j = j2; j < that.col + j2; j++)
                    (*this)[i][j] = that[i - i2][j - j2];
        }
    }

    vector &operator[](size_t index) const {
        assert(index < row);
        return arrvec[index];
    }

    void updateValue(double val){
        for (size_t i = 0; i < row; i++)
            for (size_t j = 0; j < col; j++)
                (*this)[i][j] = val;
    }

    void updateDiagValue(double val){
        for (size_t i = 0; i < row; i++) (*this)[i][i] = val;
    }

    matrix &operator=(const matrix &that) {
        for (size_t i = 0; i < row; i++)
            for (size_t j = 0; j < col; j++)
                (*this)[i][j] = that[i][j];
        return *this;
    }

    matrix operator+(const matrix &that) const {
        matrix result(row, col);
        result.updateValue(0);
        for (size_t i = 0; i < row; i++)
            for (size_t j = 0; j < col; j++)
                result[i][j] = that[i][j] + (*this)[i][j];
        return result;
    }

    matrix operator*(const matrix &that) const {
        matrix result(row, that.col);
        result.updateValue(0);
        for (size_t i = 0; i < row; i++)
            for (size_t j = 0; j < that.col; j++)
                for (size_t k = 0; k < col; k++)
                    result[i][j] += (*this)[i][k] * that[k][j];
        return result;
    }

    matrix operator*(const double a) {
        matrix result(row, col);
        for (size_t i = 0; i < row; i++)
            for (size_t j = 0; j < col; j++)
                result[i][j] = (*this)[i][j] * a;
        return result;
    }

    matrix transpose() {
        matrix result(col, row);
        for (size_t i = 0; i < row; i++)
            for (size_t j = 0; j < col; j++)
                result[j][i] = (*this)[i][j];
        return result;
    }

    double nornEuclidean() {
        double result = 0;
        for (size_t i = 0; i < row; i++)
            result += pow((*this)[i][0], 2);
        return sqrt(result);
    }

    int zeroCheck(){
        for (size_t i = 1; i < row; i++) {
            if ((*this)[i][0] != 0) return 1;
        }
        return 0;
    }


    matrix hBuilder(size_t k) {
        matrix result(row, col);
        result.updateValue(0);
        result.updateDiagValue(1);
        for (size_t i = 0; i < row; i++) result[i][i] = 1;
        matrix a(*this, k - 1, k - 1, row - 1, k - 1, 'a');
        double r = a.nornEuclidean();
        if ( (r != 0) && a.zeroCheck()) {
            matrix re(row - k + 1, 1);
            re.updateValue(0);
            re[0][0] = -1;
            re = re * r;
            if (a[0][0] == 0)
                a = a + re;
            else {
                if (a[0][0] > 0)
                    re = re * (-1);
                a = a + re;
            }
            r = a.nornEuclidean();
            a = a * (1.0/r);
            matrix A(a * a.transpose() * (-2), row, col, row - a.row, col - a.row, 'b');
            result = result + A;
        }
        return result;
    }

    void solve(const matrix& x, const matrix& b) {
        for (size_t i = row - 1; i != 0; i--) {
            x[i][0] = b[i][0];
            for(size_t j = i + 1; j < col; j++)
                x[i][0] -= (*this)[i][j]*x[j][0];
            x[i][0] /= (*this)[i][i];
        }
        x[0][0] = b[0][0];
        for (size_t j = 1; j < col; j++)
            x[0][0] -= (*this)[0][j] * x[j][0];
        x[0][0] /= (*this)[0][0];
    }


    void input() {
        for (size_t i = 0; i < row; i++)
            for (size_t j = 0; j < col; j++) cin >> (*this)[i][j];
    }

    void print() {
        for (size_t i = 0; i < row; i++) {
            for (size_t j = 0; j < col; j++) cout << fixed << setprecision(8) << (*this)[i][j] << " ";
            cout << endl;
        }
    }

    ~matrix() {
        for (size_t i = 0; i < row; i++) delete[] arrvec[i].value;
        delete[] arrvec;
    }

private:
    size_t row, col;
    vector *arrvec;
};

int main() {
    size_t n;
    cout << "Enter dim : ";
    cin >> n;
    cout << endl;

    cout << "Enter matrix A :" << endl;
    matrix A(n,n);
    A.input();

    cout << "Enter vector b : ";
    matrix b(n,1);
    b.input();
    cout << endl;

    matrix h(n,n), x(n,1);
    x.updateValue(0);

    cout << "Matrix A: " << endl;
    A.print();
    cout << endl;
    cout << "Vector B: " << endl;
    b.print();
    cout << endl;

    for (size_t i = 1; i <= n; i++) {
        h = A.hBuilder(i);
        A = (h * A);
        b = (h * b);
    }
    cout << "Matrix U: " << endl;
    A.print();
    cout << endl;
    cout << "Vector B (recalculated): " << endl;
    b.print();
    cout << endl;
    A.solve(x,b);

    cout << "Vector of solution (x) : " << endl;
    x.print();
    return 0;
}