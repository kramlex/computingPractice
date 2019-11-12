#include <iostream>
#include <stdlib.h>
#include <assert.h>
#include <time.h>
#include <limits.h>
#include <math.h>

class matrix {
    class vector;

    size_t row, col; //row - строки, col - столбцы
    vector *arrvec;

    class vector {
        size_t len;
    public:
        double *value;

        vector() {};

        vector(size_t n) {
            this->len = n;
            value = new double[this->len];
        }

        double &operator[](size_t index) const {
            assert(index < this->len);
            return value[index];
        }
    };

    void init(const size_t n, const size_t m) {
        this->row = n;
        this->col = m;
        assert((row > 0) || (col > 0));
        arrvec = new vector[row];
        for (size_t i = 0; i < row; i++)
            arrvec[i] = vector(col);
    }

public:
    matrix(size_t n, size_t m) {
        this->init(n, m);
    }

    matrix(const matrix &that, const size_t i1, const size_t j1, const size_t i2, const size_t j2, char key) {
        assert((key == 'c') || (key == 'b'));
        if (key == 'c') {
            assert((i1 <= i2) && (j1 <= j2) && (i2 < that.row) && (j2 < that.col));
            this->init(i2 - i1 + 1, j2 - j1 + 1);
            for (size_t i = i1; i <= i2; i++)
                for (size_t j = j1; j <= j2; j++)
                    (*this)[i - i1][j - j1] = that[i][j];
        }
        if (key == 'b') {
            assert((i2 + that.row - 1 < i1) && (j2 + that.col - 1 < j1));
            this->init(i1, j1);
            this->zeroing();
            for (size_t i = i2; i < that.row + i2; i++)
                for (size_t j = j2; j < that.col + j2; j++)
                    (*this)[i][j] = that[i - i2][j - j2];
        }
    }

    vector &operator[](size_t index) const {
        assert(index < row);
        return arrvec[index];
    }

    void zeroing() {
        for (size_t i = 0; i < this->row; i++)
            for (size_t j = 0; j < this->col; j++) (*this)[i][j] = 0;
    }

    matrix &operator=(const matrix &that) {
        assert((this->row == that.row) || (this->col == that.col));
        for (size_t i = 0; i < this->row; i++)
            for (size_t j = 0; j < this->col; j++) (*this)[i][j] = that[i][j];
        return *this;
    }

    matrix operator+(const matrix &that) const {
        assert((this->row == that.row) || (this->col == that.col));
        matrix result(this->row, this->col);
        result.zeroing();
        for (size_t i = 0; i < this->row; i++)
            for (size_t j = 0; j < this->col; j++) result[i][j] = that[i][j] + (*this)[i][j];
        return result;
    }

    matrix operator*(const matrix &that) const {
        assert(this->col == that.row);
        matrix result(this->row, that.col);
        result.zeroing();
        for (size_t i = 0; i < this->row; i++)
            for (size_t j = 0; j < that.col; j++)
                for (size_t k = 0; k < this->col; k++) result[i][j] += (*this)[i][k] * that[k][j];
        return result;
    }

    matrix operator*(const double a) {
        matrix result(this->row, this->col);
        for (size_t i = 0; i < this->row; i++)
            for (size_t j = 0; j < this->col; j++) result[i][j] = (*this)[i][j] * a;
        return result;
    }

    matrix transpose() {
        matrix result(this->col, this->row);
        for (size_t i = 0; i < this->row; i++)
            for (size_t j = 0; j < this->col; j++)
                result[j][i] = (*this)[i][j];
        return result;
    }

    double euclid_norm() {
        assert(this->col == 1);
        double result = 0;
        for (size_t i = 0; i < this->row; i++) result += pow((*this)[i][0], 2);
        return sqrt(result);
    }

    int zero_check(){
        assert(this->col == 1);
        for (size_t i = 1; i < this->row; i++) if ((*this)[i][0] != 0) return 1;
        return 0;
    }

    matrix H_build(size_t k) {
        assert(this->row == this->col);
        matrix result(this->row, this->col);
        result.zeroing();
        for (size_t i = 0; i < this->row; i++) result[i][i] = 1;
        matrix a(*this, k - 1, k - 1, this->row - 1, k - 1, 'c');
        double r = a.euclid_norm();
        if ((r != 0) && (a.zero_check())) {
            matrix re(this->row - k + 1, 1);
            re.zeroing();
            re[0][0] = -1;
            re = re * r;
            if (a[0][0] == 0) {
                a = a + re;
            } else {
                if (a[0][0] > 0) re = re * (-1);
                a = a + re;
            }
            r = a.euclid_norm();
            a = a * (1.0/r);
            matrix A(a * a.transpose() * (-2), this->row, this->col, this->row - a.row, this->col - a.row, 'b');
            result = result + A;
        }
        return result;
    }

    void solve(const matrix& x, const matrix& b) {
        for (size_t i = this->row - 1; i != 0; i--) {
            x[i][0] = b[i][0];
            for (size_t j = i + 1; j < this->col; j++) x[i][0] -= (*this)[i][j]*x[j][0];
            x[i][0] /= (*this)[i][i];
        }
        x[0][0] = b[0][0];
        for (size_t j = 1; j < this->col; j++) x[0][0] -= (*this)[0][j] * x[j][0];
        x[0][0] /= (*this)[0][0];
    }


    void input() {
        for (size_t i = 0; i < row; i++)
            for (size_t j = 0; j < col; j++) scanf("%lf", &(*this)[i][j]);
    }

    void random_input() {
        for (size_t i = 0; i < row; i++)
            for (size_t j = 0; j < col; j++) (*this)[i][j] = rand() % 100;
    }

    void output() {
        for (size_t i = 0; i < row; i++) {
            for (size_t j = 0; j < col; j++) printf("%.8lf ", (*this)[i][j]);
            printf("\n");
        }
    }

    ~matrix() {
        for (size_t i = 0; i < row; i++)
            if (arrvec[i].value) delete[] arrvec[i].value;
        if (arrvec) delete[] arrvec;
    }

};

int main() {
    srand((unsigned) time(0));
    size_t dim;
    scanf("%lld", &dim);
    matrix A(dim, dim);
    matrix x(dim, 1);
    matrix b(dim, 1);
    matrix H(dim, dim);
    x.zeroing();
    A.input();
    b.input();
    A.output();
    printf("\n");
    for (size_t i = 1; i <= dim; i++) {
        H = A.H_build(i);
        A = (H * A);
        b = (H * b);
    }
    H.output(); printf("\n");
    A.output(); printf("\n"); b.output();
    A.solve(x,b);
    x.output();
    return 0;
}