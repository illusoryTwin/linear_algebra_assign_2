#include <iostream>
#include <vector>
#include <iomanip>
#include <cmath>
#include <cstdlib>
#include <string>


using namespace std;

class ColumnVector {
private:
    int size;
    vector<double> data;

public:
    ColumnVector(int n) : size(n), data(n) {}

    int getSize() {
        return size;
    }

    vector<double> &getData() {
        return data;
    }

    void read() {
        for (int i = 0; i < size; ++i) {
            cin >> data[i];
        }
    }

    void print() const {
        cout << fixed << setprecision(4);
        for (double x: data) {
            cout << x << endl;
        }
    }

    // Overload the + operator
    ColumnVector operator+(const ColumnVector &other) const {
        if (size != other.size) {
            throw std::invalid_argument("Size mismatch");
        }
        ColumnVector result(size);
        for (int i = 0; i < size; ++i) {
            result.data[i] = data[i] + other.data[i];
        }
        return result;
    }

    // Overload the - operator
    ColumnVector operator-(const ColumnVector &other) const {
        if (size != other.size) {
            throw std::invalid_argument("Size mismatch");
        }
        ColumnVector result(size);
        for (int i = 0; i < size; ++i) {
            result.data[i] = data[i] - other.data[i];
        }
        return result;
    }

    // Overload the * operator for scalar multiplication
    ColumnVector operator*(double scalar) const {
        ColumnVector result(size);
        for (int i = 0; i < size; ++i) {
            result.data[i] = data[i] * scalar;
        }
        return result;
    }

    // Overload the / operator for scalar division
    ColumnVector operator/(double scalar) const {
        if (scalar == 0) {
            throw std::invalid_argument("Division by zero");
        }
        ColumnVector result(size);
        for (int i = 0; i < size; ++i) {
            result.data[i] = data[i] / scalar;
        }
        return result;
    }
};

class Matrix {

private:
    int rows, columns;
    vector<vector<double>> data;


public:

    Matrix() {};

    Matrix(int n, int m) : rows(n), columns(m), data(n, vector<double>(m, 0)) {};

    int getRows() const {
        return rows;
    }

    int getCols() const {
        return columns;
    }

    double &at(int row, int col) {
        return data[row][col];
    }

    const double &at(int row, int col) const {
        return data[row][col];
    }

    void read() {
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < columns; j++) {
                cin >> data[i][j];
            }
        }
    }

    void print() const {
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < columns; j++) {
                if (abs(data[i][j]) < 1e-9) {
                    cout << "0.0000" << " ";
                } else {
                    cout << fixed << setprecision(4) << data[i][j] << " ";
                }
            }
            cout << endl;
        }
    }

    Matrix &operator=(const Matrix &other) {
        rows = other.rows;
        columns = other.columns;
        data = other.data;
        return *this;
    }

    Matrix operator+(Matrix &other) {
        if (rows != other.rows || columns != other.columns) {
            cout << "Error: Matrix dimensions do not match!" << endl;
            return Matrix();
        }
        Matrix res(rows, columns);
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < columns; j++) {
                res.data[i][j] = data[i][j] + other.data[i][j];
            }
        }
        return res;
    }

    Matrix operator-(Matrix &other) {
        if (rows != other.rows || columns != other.columns) {
            cout << "Error: Matrix dimensions do not match!" << endl;
            return Matrix();
        }
        Matrix res(rows, columns);
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < columns; j++) {
                res.data[i][j] = data[i][j] - other.data[i][j];
            }
        }
        return res;
    }

    Matrix operator*(const Matrix &other) const {
        if (columns != other.rows) {
            cerr << "Error: Multiplication : Matrix dimensions do not match!" << endl;
            return Matrix();
        }
        Matrix res(rows, other.columns);
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < other.columns; j++) {
                for (int k = 0; k < columns; k++) {
                    res.data[i][j] += data[i][k] * other.data[k][j];
                }
            }
        }
        return res;
    }

    ColumnVector operator*(ColumnVector other) const {
        if (columns != other.getSize()) {
            cerr << "Error: Matrix-ColumnVector multiplication: dimensions do not match!" << endl;
            return ColumnVector(0);
        }
        ColumnVector res(rows);
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < columns; j++) {
                res.getData()[i] += data[i][j] * other.getData()[j];
            }
        }
        return res;
    }


    Matrix transpose() const {
        Matrix res(columns, rows);
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < columns; j++) {
                res.data[j][i] = data[i][j];
            }
        }
        return res;
    }

    Matrix getDiagonalMatrix() {
        Matrix res(rows, columns);
        for (int i = 0; i < min(rows, columns); i++) {
            res.data[i][i] = data[i][i];
        }
        return res;
    }

    Matrix subtractFromIdentity() {
        if (rows != columns) {
            cerr << "Error: Matrix must be square!" << endl;
            return Matrix();
        }

        Matrix res(rows, columns);
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < columns; j++) {
                if (i == j) {
                    res.data[i][j] = 1 - data[i][j];
                } else {
                    res.data[i][j] = -data[i][j];
                }
            }
        }
        return res;
    }

    double determinant() const {
        if (rows != columns) {
            cerr << "Error: Matrix must be square!" << endl;
            return 0;
        }
        if (rows == 1) {
            return data[0][0];
        } else if (rows == 2) {
            return data[0][0] * data[1][1] - data[0][1] * data[1][0];
        }

        double det = 0;
        for (int i = 0; i < columns; i++) {
            Matrix submatrix(rows - 1, columns - 1);
            for (int j = 1; j < rows; j++) {
                for (int k = 0; k < columns; k++) {
                    if (k < i) {
                        submatrix.data[j - 1][k] = data[j][k];
                    } else if (k > i) {
                        submatrix.data[j - 1][k - 1] = data[j][k];
                    }
                }
            }
            det += (i % 2 == 0 ? 1 : -1) * data[0][i] * submatrix.determinant();
        }
        return det;
    }

    Matrix adjugate() const {
        if (rows != columns) {
            cerr << "Error: Matrix must be square!" << endl;
            return Matrix();
        }

        Matrix adj(rows, columns);
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < columns; j++) {
                Matrix submatrix(rows - 1, columns - 1);
                for (int k = 0; k < rows; k++) {
                    for (int l = 0; l < columns; l++) {
                        if (k != i && l != j) {
                            submatrix.data[k < i ? k : k - 1][l < j ? l : l - 1] = data[k][l];
                        }
                    }
                }
                adj.data[i][j] = ((i + j) % 2 == 0 ? 1 : -1) * submatrix.determinant();
            }
        }
        return adj.transpose();
    }

    Matrix inverse() const {
        if (rows != columns) {
            cerr << "Error: Matrix must be square!" << endl;
            return Matrix();
        }
        double det = determinant();

        if (det == 0) {
            cerr << "Error: Matrix is singular, inverse not possible!" << endl;
            return Matrix();
        }

        Matrix adj = adjugate();
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < columns; j++) {
                adj.data[i][j] /= det;
            }
        }
        return adj;
    }

    Matrix getLowerTriangular() {
        Matrix res(rows, columns);
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j <= i; j++) {
                res.data[i][j] = data[i][j];
            }
        }
        return res;
    }

    Matrix getUpperTriangular() {
        Matrix res(rows, columns);
        for (int i = 0; i < rows; i++) {
            for (int j = i + 1; j < columns; j++) {
                res.data[i][j] = data[i][j];
            }
        }
        return res;
    }
};


istream &operator>>(istream &in, Matrix &mat) {
    mat.read();
    return in;
}

ostream &operator<<(ostream &out, const Matrix &mat) {
    mat.print();
    return out;
}

istream &operator>>(istream &in, ColumnVector &vec) {
    vec.read();
    return in;
}

// Define the output stream operator for ColumnVector outside the class
ostream &operator<<(ostream &out, const ColumnVector &vec) {
    vec.print();
    return out;
}

std::string get_polynom(ColumnVector mat, int n) {
    std::string equation = "";
    for (int i = 0; i <= n; i++) { // polynomial
        equation += std::to_string(mat.getData()[i]) + "*x**" + std::to_string(n - i);
        //equation += std::to_string(mat[i]) + "*x**" + std::to_string(n - i);
        //equation += mat[i][0] + "*pow(x, " + std::to_string(n - i - 1) + ")";

        if (i != n)
            equation += "+";
    }
    return equation;
}

int main() {


    int m;
    std::cin >> m;
    std::vector<double> t(m), b(m);

    for (int i = 0; i < m; i++) {
        std::cin >> t[i] >> b[i];
    }
    int n;
    cin >> n;
    ColumnVector b_col(m);
    for (int i = 0; i < m; i++) {
        b_col.getData()[i] = b[i];
    }
    Matrix A(m, n + 1);
    int k = 0;
    for (int i = 0; i < m; i++) {
        for (int j = 0; j <= n; j++) {
            A.at(i, j) = pow(t[i], j);
            k++;
        }
    }
    cout << "A:" << endl;
    cout << A;
    Matrix AT = A.transpose();
    Matrix AT_A = AT * A;
    cout << "A_T*A:" << endl;
    cout << AT_A;
    Matrix AT_Ainv = AT_A.inverse();
    cout << "(A_T*A)^-1:" << endl;
    cout << AT_Ainv;
    ColumnVector ATb = AT * b_col;
    cout << "A_T*b:" << endl;
    cout << ATb;
    Matrix AT_AinvAT = AT_Ainv * AT;
    ColumnVector X = AT_AinvAT * b_col;
    cout << "x~:" << endl;
    cout << X;

    return 0;
}