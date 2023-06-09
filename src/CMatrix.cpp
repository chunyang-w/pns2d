# include <iostream>
# include <fstream>
# include <CMatrix.h>

using namespace std;

Cmatrix::Cmatrix()
{
  mat_1D = nullptr;
  mat_2D = nullptr;
}

Cmatrix::Cmatrix(int imax, int jmax)
{
  n = imax;
  m = jmax;
  mat_1D = new double[n * m];
  mat_2D = new double* [n];
  for (int i = 0; i < n; i++)
    mat_2D[i] = &mat_1D[i * m];
}

Cmatrix::~Cmatrix()
{
  delete[] mat_1D;
  delete[] mat_2D;
}

void Cmatrix::to_tab_file(string fname)
{
  fstream fout;
  fout.open(fname, ios::out);
  if (fout.fail())
  {
    cout << "Error openin file!" << endl;
    cout.flush();
    exit(0);
  }

  for (int i = 0; i < n; i++)
  {
    for (int j = 0; j < m; j++)
      fout << mat_2D[i][j] << "\t";
    fout << endl;
  }
  fout.close();
}

// copy constructor
Cmatrix::Cmatrix(const Cmatrix& other)
{
    n = other.n;
    m = other.m;
    mat_1D = new double[n * m];
    mat_2D = new double* [n];
    for (int i = 0; i < n; i++)
    {
        mat_2D[i] = &mat_1D[i * m];
        for (int j = 0; j < m; j++)
        {
            mat_2D[i][j] = other.mat_2D[i][j];
        }
    }
}

void Cmatrix::fill(double val)
{
  for (int i = 0; i < n; i++)
    for (int j = 0; j < m; j++)
      mat_2D[i][j] = val;
}

void Cmatrix::swap(Cmatrix &a, Cmatrix &b) {
    std::swap(a.mat_1D, b.mat_1D);
    std::swap(a.mat_2D, b.mat_2D);
    std::swap(a.n, b.n);
    std::swap(a.m, b.m);
}