//Ekaterina Mozhegova
#include <iostream>
#include <cmath>
#include <iomanip>
#include <vector>
using namespace std;
class ColumnVector {
public:
    ColumnVector(int n) : size(n), data(n) {}

    void read() {
        for (int i = 0; i < size; ++i) {
            cin >> data[i];
        }
    }

    void print() const {
        cout << fixed << setprecision(4);
        for (double x : data) {
            cout << x << endl;
        }
    }

    // Overload the + operator
    ColumnVector operator+(const ColumnVector& other) const {
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
    ColumnVector operator-(const ColumnVector& other) const {
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

    int size;
    vector<double> data;
};

bool isDiagonallyDominant(const vector<vector<double>>& A) {
    int n = A.size();
    for (int i = 0; i < n; i++) {
        double diagonal = abs(A[i][i]);
        double sum = 0.0;
        for (int j = 0; j < n; j++) {
            if (i != j) {
                sum += abs(A[i][j]);
            }
        }
        if (diagonal <= sum) {
            return false;
        }
    }
    return true;
}
class Matrix {

protected:

public:
    int rows, columns;
    vector<vector<double>> data;
    Matrix() {};

    Matrix(int n, int m) : rows(n), columns(m),  data(n, vector<double>(m, 0)) {};
    int getRows() const {
        return rows;
    }
    int getCols() const {
        return columns;
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
                if (abs(data[i][j]) < 1e-9){
//                    cout  << fixed << setprecision(4) << 0 << " ";
                    cout << "0.0000" << " ";
                } else{
                    cout << fixed << setprecision(4) << data[i][j] << " ";
                }
            }
            cout << endl;
        }
    }
    Matrix& operator=(const Matrix &other) {
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
    Matrix operator*(const Matrix& other) const {
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
    ColumnVector operator*(const ColumnVector& other) const {
        if (columns != other.size) {
            cerr << "Error: Matrix-ColumnVector multiplication: dimensions do not match!" << endl;
            return ColumnVector(0);
        }
        ColumnVector res(rows);
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < columns; j++) {
                res.data[i] += data[i][j] * other.data[j];
            }
        }
        return res;
    }


    Matrix transpose() const{
        Matrix res( columns, rows);
        for (int i =0; i < rows; i++){
            for (int j = 0; j< columns; j++){
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


istream& operator>>(istream& in, Matrix& mat) {
    mat.read();
    return in;
}
ostream& operator<<(ostream& out, const Matrix& mat) {
    mat.print();
    return out;
}
istream& operator>>(istream& in, ColumnVector& vec) {
    vec.read();
    return in;
}

// Define the output stream operator for ColumnVector outside the class
ostream& operator<<(ostream& out, const ColumnVector& vec) {
    vec.print();
    return out;
}
class IdentityMatrix : public Matrix {
public:
    IdentityMatrix(int n) : Matrix(n, n) {
        for (int i = 0; i < n; i++) {
            data[i][i] = 1;
        }
    }
};
double get_length_of_vector_subtraction(const vector<double> &a, const vector<double> &b) {
    double sum = 0;
    for (size_t i = 0; i < a.size(); i++) {
        sum += pow(a[i] - b[i], 2);
    }
    return sqrt(sum);
}

void SeidelMethod(const Matrix &IminusBInv, const Matrix &C, const ColumnVector &beta, double epsilon) {
    ColumnVector x0 = beta;
    x0 = beta;
    cout << "x(" << 0 << "):\n";
    cout << x0;
    ColumnVector x1(IminusBInv.getRows());
    int k = 1;
    while (true) {
        x1 = IminusBInv * C * x0 + IminusBInv * beta;
        auto norm = get_length_of_vector_subtraction(x1.data, x0.data);
        cout << "e: " <<  norm << "\n";
        cout << "x(" << k << "):\n";
        cout << x1;
        x0 = x1;
        k++;
        if (norm < epsilon) {
            break;
        }
    }
}


int main() {
    int n;
    cin >> n;
    Matrix A(n, n);
    cin >> A;

    int n2;
    cin >> n2;
    ColumnVector b(n2);
    cin >> b;

    double epsilon;
    cin >> epsilon;

    if (!isDiagonallyDominant(A.data)) {
        cout << "The method is not applicable!";
        return 0;
    }

    Matrix D = A.getDiagonalMatrix();
    Matrix D_1 = D.inverse();

    ColumnVector beta = D_1 *b;
    cout << "beta:" << endl;
    cout << beta;

    Matrix X = D_1*A;
    Matrix alpha = X.subtractFromIdentity();
    cout << "alpha:" << endl;
    cout << alpha;

    Matrix B = alpha.getLowerTriangular();
    cout << "B:" << endl;
    cout << B;
    Matrix Breal = A.getLowerTriangular();

    Matrix C = alpha.getUpperTriangular();
    cout << "C:" << endl;
    cout << C;

    Matrix IminusB = B.subtractFromIdentity();
    cout << "I-B:" << endl;
    cout << IminusB;

    Matrix IminusBInv = IminusB.inverse();
    cout << "(I-B)_-1:" << endl;
    cout << IminusBInv;
    SeidelMethod(IminusBInv, C, beta,epsilon);
}