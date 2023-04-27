#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>

using namespace std;

class Matrix {
public:
    Matrix(int n) : data(n, vector<double>(n)) {}

    int size() const { return data.size(); }

    vector<double>& operator[](int i) { return data[i]; }
    const vector<double>& operator[](int i) const { return data[i]; }

    bool isDiagonallyDominant() const {
        int n = data.size();
        for (int i = 0; i < n; i++) {
            double diagonal = abs(data[i][i]);
            double sum = 0.0;
            for (int j = 0; j < n; j++) {
                if (i != j) {
                    sum += abs(data[i][j]);
                }
            }
            if (diagonal <= sum) {
                return false;
            }
        }
        return true;
    }

private:
    vector<vector<double>> data;
};

class ColumnVector {
public:
    ColumnVector(int n) : data(n) {}

    int size() const { return data.size(); }

    double& operator[](int i) { return data[i]; }
    const double& operator[](int i) const { return data[i]; }

    double getLengthOfVectorSubtraction(const ColumnVector& other) const {
        double sum = 0;
        for (int i = 0; i < data.size(); i++) {
            sum += (data[i] - other.data[i]) * (data[i] - other.data[i]);
        }
        return sqrt(sum);
    }

private:
    vector<double> data;
};

int main() {
    int n;
    cin >> n;
    Matrix A(n);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            cin >> A[i][j];
        }
    }

    int n2;
    cin >> n2;
    ColumnVector b(n);
    for (int i = 0; i < n; i++) {
        cin >> b[i];
    }

    double epsilon;
    cin >> epsilon;

    if (!A.isDiagonallyDominant()) {
        cout << "The method is not applicable!";
        return 0;
    }

    Matrix alpha(n);
    ColumnVector beta(n);

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (i == j) {
                alpha[i][j] = 0.0;
            } else {
                alpha[i][j] = -A[i][j] / A[i][i];
            }
        }
        beta[i] = b[i] / A[i][i];
    }

    cout << "alpha:" << endl;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            cout << fixed << setprecision(4) << alpha[i][j] << " ";
        }
        cout << endl;
    }

    cout << "beta:" << endl;
    for (int i = 0; i < n; i++) {
        cout << fixed << setprecision(4) << beta[i] << endl;
    }

    ColumnVector x(n);
    int iteration = 0;
    double max_error = 0.0;
    while (true) {

        ColumnVector next_x(n);
        for (int i = 0; i < n; i++) {
            next_x[i] = beta[i];
            for (int j = 0; j < n; j++) {
                next_x[i] += alpha[i][j] * x[j];
            }
        }
        max_error = x.getLengthOfVectorSubtraction(next_x);
        if (iteration != 0)
            cout << "e: " << fixed << setprecision(4) << max_error << endl;
        x = next_x;
        cout << "x(" << iteration << "):" << endl;
        bool exit = false;
        if (max_error < epsilon) {
            exit = true;
        }
        for (int i = 0; i < n; i++) {
            cout << fixed << setprecision(4) << x[i];
            if (!(exit && i == n-1))
                cout << endl;
        }
        iteration++;
        if (exit) {
            break;
        }
    }

    return 0;
}

