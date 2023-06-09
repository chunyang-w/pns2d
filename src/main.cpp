# include <iostream>
# include <CMatrix.h>

int main() {
  Cmatrix A(3, 3);
  A.fill(1.0);

  std::cout << "A is: " << std::endl;

  for (int i = 0; i < A.n; i++) {
    for (int j = 0; j < A.m; j++) {
      std::cout << A.mat_2D[i][j] << " ";
    }
    std::cout << std::endl;
  }
  std::cout << std::endl;
}