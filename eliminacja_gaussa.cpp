#include <iostream>
#include <fstream>
#include <string>

using namespace std;

typedef double** matrix;

void wypelnij(matrix& b, const int& m_, const int& n_, const double w = 0);
void usun(matrix& b, const int& n);
double* utworz(const int n, const double war = 0);
int czytajZpliku(string nazwaPliku, matrix& b, const int& m_, const int& n_);
double* gaus(string nazwaPliku, const int row, const int col);

int main()
{
    cout << "\t UKLAD ROWNAN" << endl;
    cout << " 0.5*x1 -2*x2 +0.5*x3 -3*x4 +3*x5=-4 " << endl;
    cout << " 0.5*x1 -2*x2 +2*x3 +0.5*x4 +2*x5=0.5 " << endl;
    cout << " -4*x1 -3*x2 +2*x3 +2*x4 +0.5*x5=-3 " << endl;
    cout << " 4*x1 +0.5*x2 +2*x3 -4*x4 +0.5*x5=4 " << endl;
    cout << " 4*x1 -3*x2 -4*x3 -5*x4 +4*x5=0.5 " << endl;
    int m = 5, n = 5;

    string plik = "A.txt";
 
    gaus(plik, m, n);

    return 0;
}

void wypelnij(matrix& b, const int& m_, const int& n_, const double w)
{
    b = new double* [m_];
    for (int i = 0; i < m_; ++i) {
        b[i] = new double[n_];
        for (int j{ 0 }; j < n_; j++)
        {
            b[i][j] = w;
        }
    }
}

double* utworz(const int n, const double war) {
    double* nowy = new double[n];
    for (int j{ 0 }; j < n; j++) {
        nowy[j] = war;
    }
    return nowy;
}

void usun(matrix& b, const int& n)
{
    for (int i = 0; i < n; i++)
    {
        delete[] b[i];
    }
    delete[]b;
}

int czytajZpliku(string nazwaPliku, matrix& b, const int& m_, const int& n_)
{
    ifstream pA;
    pA.open(nazwaPliku.c_str());
    for (int i = 0; i < m_; ++i)
        for (int j = 0; j < n_; ++j)
            pA >> b[i][j];
    pA.close();
    return 0;
}
double* gaus(string nazwaPliku, const int m, const int n)
{
    matrix A;
    wypelnij(A, m, n + 1);
    czytajZpliku(nazwaPliku, A, m, n + 1);

    for (int i = 0; i < m - 1; i++)
    {
        for (int j = i + 1; j < n; j++)
        {
            for (int s = i + 1; s < n + 1; s++)
            {
                A[j][s] -= (A[j][i] / A[i][i]) * A[i][s];
            }
        }
    }

    double* X = utworz(m);
    
    X[m - 1] = A[m - 1][n] / A[m - 1][m - 1];

    for (int i = m - 1; i >= 0; i--)
    {
        double suma = 0;
        for (int j = i + 1; j < m; j++)
        {
            suma += A[i][j] * X[j];
        }
        X[i] = (A[i][n] - suma) / A[i][i];
    }

    cout << endl;
    cout << "WYNIKI: " << endl;

    for (int i = 0; i < m; i++)
    {
        cout << "x" << i + 1 << " = " << X[i] << endl;
    }

    usun(A, m);

    return 0;
}