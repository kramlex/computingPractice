//
// Created by markd on 10.12.2019.
//
#include <vector>
#include <iostream>
#include <utility>
#include <exception>
#include <cmath>
using namespace std;

class Vector;
// Класс-матрица
class Matrix {

    vector<vector<double> > matrix;
    int n; int m;

public:

    Matrix(int); // Создать квадратную матрицу заданного размера
    Matrix(int, int); // Создать прямоугольную матрицу
    Matrix(vector<double> const &v); // Создать матрицу из набора векторов
    Matrix(Matrix const&); // Копировать матрицу

    Matrix& operator=(Matrix const&); // Присвоить матрицу
    vector<double>& operator[](int i); // Получить строку
    Matrix operator~(void); // Транспонировать матрицу
    Vector toVector();

    int rows_count(); // Количество строк
    int colls_count(); // Количество столбцов

    Matrix inverse(); // Получить обратную матрицу
    pair<Matrix, Matrix> qr(); // Получить QR разложение

    friend Matrix operator+(Matrix const&, Matrix const&); // Сложить две матрицы
    friend Matrix operator-(Matrix const&, Matrix const&); // Вычисть две матрицы
    friend Matrix operator*(Matrix const&, Matrix const&); // Перемножить две матрицы

    friend Matrix operator*(double const, Matrix const&); // Умножить число на матрицу
    friend Matrix operator*(Matrix const&, double const); // Умножить матрицу на число

    friend Matrix operator/(Matrix const&, double const); // Разделить матрицу на число
    friend ostream& operator<<(ostream&, Matrix const&); // Распечатать матрицу
    friend istream& operator>>(istream&, Matrix&); // Считать матрицу

    static Matrix E(int n); // Создать еденичную матрицу размера n

    Vector solve(Vector); // Решить систему с заданной правой частью
};

class Vector {
    double *el;
    int n;
public:
    Vector(int); // Создать вектор фиксированного размера
    Vector(Vector const&); // Создать копию ветора
    Vector(vector<double> const&); // Создать вектор из набора чисел
    Vector& operator=(Vector const&); // Присвоить вектор

    ~Vector();

    Matrix operator~(); // Перевернуть вектор

    operator Matrix(); // Получить матрицу из вектора

    double operator[] (int i) const; // Получить i-ый компонент
    double& operator[] (int i); // Присвоить i-ый компонент

    int size(); // Размер вектора

    double length(); // Длина вектора

    Vector normalized(); // Нормализованный вектор

    friend double operator*(Vector const&, Vector const&); // Скалярное произведение
    friend Vector operator+(Vector const&, Vector const&); // Прибавить число к вектору
    friend Vector operator-(Vector const&, Vector const&); // Вычисть число из вектора

    friend Vector operator*(double const, Vector const&); // Умножить число на вектор
    friend Vector operator*(Vector const&, double const); // Умножить вектор на число

    friend Vector operator/(Vector const&, double const); // Разделить вектор на число

    friend ostream& operator<<(ostream&, Vector const&); // Распечатать вектор
    friend istream& operator>>(istream&, Vector&); // Считать вектор

    static Vector e(int, int); // Еденичный вектор
};

Matrix operator+(Matrix const&, Matrix const&); // Сумма матриц
Matrix operator-(Matrix const&, Matrix const&); // Разность матриц
Matrix operator*(Matrix const&, Matrix const&); // Произведение матриц
Matrix operator*(Matrix const&, Matrix const&); // Произведение матриц

Matrix operator*(double const, Matrix const&); // Произведение числа на матрицу
Matrix operator*(Matrix const&, double const); // Произведение матрицу на число

Matrix operator/(Matrix const&, double const); // Разделить матрицу на число

ostream& operator<<(ostream&, Matrix const&); // Вывести матрицу на печать
istream& operator>>(istream&, Matrix&); // Считать матрицу

double operator*(Vector const&, Vector const&); // Произведение векторов
Vector operator+(Vector const&, Vector const&); // Сумма векторов
Vector operator-(Vector const&, Vector const&); // Разность векторов

Vector operator*(double const, Vector const&); // Умножить число на вектор
Vector operator*(Vector const&, double const); // Умножить вектор на число

Vector operator/(Vector const&, double const); // Разделить вектор на число

ostream& operator<<(ostream&, Vector const&); // Распечатать вектор
istream& operator>>(istream&, Vector&); // Считать вектор


// *****************

Matrix::Matrix(int n, int m) {
    if (n <= 0 || m <= 0) throw invalid_argument("Wrong matrix size");

    this->n = n;
    this->m = m;

    matrix = vector<vector<double> >(n);

    for (int i = 0; i < n; i++) {
        matrix[i] = vector<double>(m);
        for (int j = 0; j < m; j++) {
            matrix[i][j] = 0.;
        }
    }
}

Matrix::Matrix(int n) {
    if (n <= 0) throw invalid_argument("Wrong matrix size");

    this->n = n;
    this->m = n;

    matrix = vector<vector<double> >(n);

    for (int i = 0; i < n; i++) {
        matrix[i] = vector<double>(n);
        for (int j = 0; j < n; j++) {
            matrix[i][j] = 0.;
        }
    }
}

Matrix::Matrix(Matrix const &a) {
    m = a.m;
    n = a.n;

    matrix = vector<vector<double> >(n);

    for (int i = 0; i < n; i++) {
        matrix[i] = vector<double>(m);
        for (int j = 0; j < m; j++) {
            matrix[i][j] = a.matrix[i][j];
        }
    }
}

Matrix::Matrix(vector<double> const &v) {
    Matrix(1, v.size());
    matrix[0] = v;
}

int Matrix::rows_count() {
    return n;
}

int Matrix::colls_count() {
    return m;
}

Matrix& Matrix::operator=(Matrix const &a) {
    m = a.m;
    n = a.n;

    matrix = vector<vector<double> >(n);

    for (int i = 0; i < n; i++) {
        matrix[i] = vector<double>(m);
        for (int j = 0; j < m; j++) {
            matrix[i][j] = a.matrix[i][j];
        }
    }

    return *this;
}

vector<double>& Matrix::operator[](int i) {
    if (i < 0 || i > n) throw out_of_range("Illegal index");
    return matrix[i];
}


Matrix Matrix::operator~() {
    Matrix result(m, n);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            result.matrix[j][i] = matrix[i][j];
        }
    }
    return result;
}

Matrix operator+(Matrix const &left, Matrix const &right) {
    if (left.n != right.n || left.m != right.m)
        throw invalid_argument("The matrices of different sizes!");

    Matrix result(left.n, right.m);
    for (int i = 0; i < left.n; i++) {
        for (int j = 0; j < left.m; j++) {
            result.matrix[i][j] = left.matrix[i][j] + right.matrix[i][j];
        }
    }

    return result;
}

Matrix operator-(Matrix const &left, Matrix const &right) {
    if (left.n != right.n || left.m != right.m)
        throw invalid_argument("The matrices of different sizes!");

    Matrix result(left.n, right.m);
    for (int i = 0; i < left.n; i++) {
        for (int j = 0; j < left.m; j++) {
            result.matrix[i][j] = left.matrix[i][j] - right.matrix[i][j];
        }
    }

    return result;
}

Matrix operator*(Matrix const &left, Matrix const &right) {
    if (left.m != right.n)
        throw invalid_argument("Matrix mismatched sizes!");

    Matrix result(left.n, right.m);

    for (int i = 0; i < left.n; i++) {
        for (int j = 0; j < right.m; j++) {
            for (int k = 0; k < left.m; k++) {
                result.matrix[i][j] += left.matrix[i][k] * right.matrix[k][j];
            }
        }
    }

    return result;
}

Matrix operator*(double const a, Matrix const& m) {
    Matrix result(m);
    for (int i = 0; i < m.n; i++) {
        for (int j = 0; j < m.m; j++) {
            result.matrix[i][j] = m.matrix[i][j] * a;
        }
    }

    return result;
}

Matrix operator*(Matrix const& m, double const a) {
    Matrix result(m);
    for (int i = 0; i < m.n; i++) {
        for (int j = 0; j < m.m; j++) {
            result.matrix[i][j] = m.matrix[i][j] * a;
        }
    }

    return result;
}

Matrix operator/(Matrix const& m, double const a) {
    if (a == 0.)
        throw runtime_error("Division by zero!");

    Matrix result(m);
    for (int i = 0; i < m.n; i++) {
        for (int j = 0; j < m.m; j++) {
            result.matrix[i][j] = m.matrix[i][j] / a;
        }
    }

    return result;
}

Matrix Matrix::E(int n) {
    Matrix result(n,n);
    for (int i = 0; i < n; i++) {
        result[i][i] = 1;
    }
    return result;
}

Vector Matrix::toVector() {
    if (n != 1 && m != 1)
        throw invalid_argument("Matrix can't be convert to vector!");
    if (n == 1) {
        Vector v(m);
        for (int  i = 0; i < m; ++i) {
            v[i] = matrix[0][i];
        }
        return v;
    }
    Vector v(n);
    for (int  i = 0; i < n; ++i) {
        v[i] = matrix[i][0];
    }
    return v;
}

ostream& operator<<(ostream& out, Matrix const& m) {
    for (int i = 0; i < m.n; i++) {
        for (int j = 0; j < m.m; j++) {
            out << m.matrix[i][j] << " ";
        }
        out << endl;
    }
    out << endl;

    return out;
}

istream& operator>>(istream& in, Matrix &m) {
    for (int i = 0; i < m.n; i++) {
        for (int j = 0; j < m.m; j++) {
            in >> m.matrix[i][j];
        }
    }

    return in;
}

Vector::Vector(int n) {
    if (n <= 0)
        throw invalid_argument("Wrong vector size");
    this->n = n;
    el = new double[n];
    for (int i = 0; i < n; i++) el[i] = 0;
}

Vector::Vector(Vector const &v) {
    n = v.n;
    el = new double[n];
    for (int i = 0; i < n; i++) el[i] = v.el[i];
}

Vector::Vector(vector<double> const &v) {
    n = v.size();
    el = new double[n];
    for (int i = 0; i < n; i++) el[i] = v[i];
}

Vector& Vector::operator=(Vector const &v) {
    n = v.n;
    el = new double[n];
    for (int i = 0; i < n; i++) el[i] = v.el[i];
    return *this;
}

Vector::~Vector() {
    delete[] el;
}

Matrix Vector::operator~() {
    Matrix result(n, 1);
    for (int i = 0; i < n; i++) {
        result[i][0] = el[i];
    }
    return result;
}

Vector::operator Matrix() {
    Matrix result(1, n);
    for (int i = 0; i < n; i++) {
        result[0][i] = el[i];
    }
    return result;
}

double Vector::operator[] (int i) const {
    if (i < 0 || i > n) throw out_of_range("Illegal index");
    return el[i];
}

double& Vector::operator[] (int i) {
    if (i < 0 || i > n) throw out_of_range("Illegal index");
    return el[i];
}

double operator*(Vector const &l, Vector const &r) {
    if (l.n != r.n)
        throw invalid_argument("Vectors of different sizes");
    double res = 0.;
    for (int i = 0; i < l.n; i++) {
        res += l.el[i] * r.el[i];
    }
    return res;
}

Vector operator+(Vector const &l, Vector const &r) {
    if (l.n != r.n)
        throw invalid_argument("Vectors of different sizes");
    Vector v(l.n);
    for (int i = 0; i < l.n; i++) {
        v.el[i] = l.el[i] + r.el[i];
    }
    return v;
}

Vector operator-(Vector const &l, Vector const &r) {
    if (l.n != r.n)
        throw invalid_argument("Vectors of different sizes");
    Vector v(l.n);
    for (int i = 0; i < l.n; i++) {
        v.el[i] = l.el[i] - r.el[i];
    }
    return v;
}

Vector operator*(double const a, Vector const& v) {
    Vector r(v.n);
    for (int i = 0; i < v.n;  i++) {
        r.el[i] = a * v.el[i];
    }
    return r;
}

Vector operator*(Vector const &v, double const a) {
    Vector r(v.n);
    for (int i = 0; i < v.n;  i++) {
        r.el[i] = a * v.el[i];
    }
    return r;
}

Vector operator/(Vector const &v, double const a) {
    if (a == 0.)
        throw runtime_error("Division by zero!");
    Vector r(v.n);
    for (int i = 0; i < v.n;  i++) {
        r.el[i] = v.el[i] / a;
    }
    return r;
}

int Vector::size() {
    return n;
}

double Vector::length() {
    return sqrt((*this) * (*this));
}

Vector Vector::normalized() {
    return *this / length();
}

ostream& operator<<(ostream& out, Vector const& v) {
    for (int i = 0; i < v.n; i++) {
        out << v.el[i] << " ";
    }
    out << endl;

    return out;
}

istream& operator>>(istream& in, Vector &v) {
    for (int i = 0; i < v.n; i++) {
        in >> v.el[i];
    }

    return in;
}

Matrix Matrix::inverse() {
    Matrix inv = Matrix::E(n);
    Matrix matrix = *this;

    for (int i = 0; i < n; i++) {
        // Фиксируем диагональный элемент
        double el = matrix[i][i];

        if (el == 0.) {
            int j = 0;
            // Ищем такую строку из последующих
            // i-ый элемент которой не нулевой
            for (j = i + 1; matrix[j][i] == 0 && j < n; j++);

            // Поэлементно прибавляем к i-ой строке j-ую
            for (int l = 0; l < n; l++) {
                matrix[i][l] += matrix[j][l];
                inv[i][l] += inv[j][l];
            }
            // Фиксируем новый диагональный элемент
            el = matrix[i][i];
        }

        // По всем строкам ниже i-ой
        for (int j = i + 1; j < n; j++) {
            // Коэффициент, при умножении на который i-ой строки
            //  И прибавлении её к j-ой мы занулим i-ый элемент j-ой строки
            double k = - matrix[j][i] / el;
            for (int l = 0; l < n; l++) {
                matrix[j][l] += matrix[i][l] * k;
                inv[j][l] += inv[i][l] * k;
            }
        }

        // Делим i-ую строку на диагональный
        for (int j = 0; j < n; j++) {
            matrix[i][j] /= el;
            inv[i][j] /= el;
        }

        // Доводим процесс поиска обратной до конца
        for (int i = n - 1; i >= 0; i--) {
            for (int j = i - 1; j >= 0; j--) {
                for (int l = 0; l < n; l++) {
                    inv[j][l] -= inv[i][l] * matrix[j][i];
                }
            }
        }

        return inv;
    }
    return inv;
}

Vector Vector::e(int n, int i) {
    Vector e(n);
    e[i] = 1;
    return e;
}

inline double sgn(double x) {
    return x > 0 ? 1 : -1;
}

pair<Matrix, Matrix> Matrix::qr() {
    Matrix R = *this, Q = Matrix::E(n);
    for (int i = 0; i < n-1; i++) {
        // Получаем i-ый вектор-столбец
        Vector t = (~R)[i];
        // Зануляем, все элементы над i-ым
        for (int j = 0; j < i; j++) {
            t[j] = 0;
        }
        // Находим вектор w необходимый для преобразования Хаусхолдера
        Vector w = (t - sgn(-t[i]) * t.length() * Vector::e(n, i)).normalized();
        // Находим матрицу преобразования
        Matrix H = Matrix::E(n) - 2 * ~w * w;
        // Применяем преобразование Хаусхолдера к A
        // Мы на шаг ближе к R из A = QR
        R = H * R;
        // Применяем к Q
        Q = H * Q;
    }
    return make_pair(~Q, R);
}

Vector Matrix::solve(Vector b) {
    Vector solution(n);
    auto QR  = qr();
    Matrix Q = QR.first;
    Matrix R = QR.second;
    // Применяем преобразование Q к вектору B
    Matrix d = ~Q * ~b;
    // Обратный ход Гаусса
    for (int i = n - 1; i >= 0; i--) {
        for (int j = i + 1; j < n; j++) {
            d[i][0] -= R[i][j] * d[j][0];
        }
        d[i][0] /= R[i][i];
    }
    for (int i = 0; i < n; i++) {
        solution[i] = d[i][0];
    }
    return solution;
}
