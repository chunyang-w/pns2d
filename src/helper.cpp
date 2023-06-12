# include <iostream>

using namespace std;

void get_subdomain(int Nx, int Ny, int id, int p, int *row_start, int *col_start, int *row_end, int *col_end, int *num_rows, int *num_cols, int *row_idx, int *col_idx, int *row_count, int *col_count) {
  // compute sutiable number of rows and columns to split the domain
  double row_col_ratio = (double) Ny / Nx;
  double residual = row_col_ratio;
  for (int i = 1 ; i <= p; i++) {
    if (p % i == 0) {
      int temp_num_rows = i;
      int temp_num_cols = p / i;
      double temp_residual = abs(row_col_ratio - (double) temp_num_rows / temp_num_cols);
      if (temp_residual < row_col_ratio) {
        residual = temp_residual;
        *num_rows = temp_num_rows;
        *num_cols = temp_num_cols;
      }
    }
  }

  // compute start/end indices for each subgrid
  int row_counts[*num_rows];
  int col_counts[*num_cols];
  int row_size = Ny / *num_rows;
  int col_size = Nx / *num_cols;

  for (int i = 0; i < *num_rows; i++) row_counts[i] = row_size;
  for (int i = 0; i < *num_cols; i++) col_counts[i] = col_size;

  int row_remainder = Ny % *num_rows;
  int col_remainder = Nx % *num_cols;

  for (int i = 0; row_remainder > 0; i++) {
    row_counts[i]++;
    row_remainder--;
  }
  for (int i = 0; col_remainder > 0; i++) {
    col_counts[i]++;
    col_remainder--;
  }

  // compute start/end indices for each subgrid
  int subdomain_row_idx = id / *num_cols;
  int subdomain_col_idx = id % *num_cols;

  int row_start_indices[*num_rows];
  int row_end_indices[*num_rows];
  int col_start_indices[*num_cols];
  int col_end_indices[*num_cols];

  for (int i = 0; i < *num_rows; i++) {
    row_start_indices[i] = 0;
    for (int j = 0; j < i; j++) row_start_indices[i] += row_counts[j];
  }
  for (int i = 0; i < *num_rows; i++) row_end_indices[i] = row_start_indices[i] + row_counts[i] - 1;
  for (int i = 0; i < *num_cols; i++) {
    col_start_indices[i] = 0;
    for (int j = 0; j < i; j++) col_start_indices[i] += col_counts[j];
  }
  for (int i = 0; i < *num_cols; i++) col_end_indices[i] = col_start_indices[i] + col_counts[i] - 1;

  *row_start = row_start_indices[subdomain_row_idx];
  *row_end = row_end_indices[subdomain_row_idx];
  *col_start = col_start_indices[subdomain_col_idx];
  *col_end = col_end_indices[subdomain_col_idx];

  *col_idx = subdomain_col_idx;
  *row_idx = subdomain_row_idx;

  *row_count = row_counts[subdomain_row_idx];
  *col_count = col_counts[subdomain_col_idx];
}