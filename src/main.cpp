# include <iostream>
# include <mpi.h>
# include <cmath>
# include <sstream>
# include <vector>
# include <chrono>
# include <helper.h>
# include <CMatrix.h>


//# define TEST
# define OUTPUT

using namespace std;

int id, p;

// const int Nx = 401;
// const int Ny = 201;

const int Nx = 201;
const int Ny = 101;

const double Lx = 0.1, Ly = 0.05;
const double rho = 1000, nu = 1e-6;
const double P_max = 0.5;
const double t_end = 50.0;
const double dt_min = 1.e-3;
const double courant = 0.01;
const double dt_out = 1;

int Nx_local, Ny_local;
int Lx_local, Ly_local;

int row_start, col_start, row_end, col_end; // **this is inclusive** start and end indices of the subdomain
int num_rows, num_cols; // interms of the layout of the subdomain, what is the number of rows and columns
int row_idx, col_idx; // the row and column index of the subdomain
int row_count, col_count; // how many rows and columns are there in the subdomain

bool is_top, is_bottom, is_left, is_right;
int i_start, i_end, j_start, j_end; // **start inclusive end exclusive** :start and end indices of the computation grid

Cmatrix * P, * P_old, * u, * u_old, * v, * v_old, * PPrhs;
double dx, dy, dt, t;

int time_it;
int jacobian_it;
long jacobian_time;


MPI_Datatype MPI_matcol;

void calculate_ppm_RHS_central(void) {
	for (int i = i_start; i < i_end; i++)
		for (int j = j_start; j < j_end; j++)
		{
			PPrhs->mat_2D[i][j] = rho / dt * ((u->mat_2D[i + 1][j] - u->mat_2D[i - 1][j]) / (2. * dx) + (v->mat_2D[i][j + 1] - v->mat_2D[i][j - 1]) / (2. * dy));
		}
}

void calculate_intermediate_velocity(void) {
	for (int i = i_start; i < i_end; i++) {
		for (int j = j_start; j < j_end; j++)
		{
			//viscous diffusion
			u->mat_2D[i][j] = u_old->mat_2D[i][j] + dt * nu * ((u_old->mat_2D[i + 1][j] + u_old->mat_2D[i - 1][j] - 2.0 * u_old->mat_2D[i][j]) / (dx * dx) + (u_old->mat_2D[i][j + 1] + u_old->mat_2D[i][j - 1] - 2.0 * u_old->mat_2D[i][j]) / (dy * dy));
			v->mat_2D[i][j] = v_old->mat_2D[i][j] + dt * nu * ((v_old->mat_2D[i + 1][j] + v_old->mat_2D[i - 1][j] - 2.0 * v_old->mat_2D[i][j]) / (dx * dx) + (v_old->mat_2D[i][j + 1] + v_old->mat_2D[i][j - 1] - 2.0 * v_old->mat_2D[i][j]) / (dy * dy));
			//advection - upwinding
			if (u->mat_2D[i][j] > 0.0)
			{
				u->mat_2D[i][j] -= dt * u_old->mat_2D[i][j] * (u_old->mat_2D[i][j] - u_old->mat_2D[i - 1][j]) / dx;
				v->mat_2D[i][j] -= dt * u_old->mat_2D[i][j] * (v_old->mat_2D[i][j] - v_old->mat_2D[i - 1][j]) / dx;
			}
			else
			{
				u->mat_2D[i][j] -= dt * u_old->mat_2D[i][j] * (u_old->mat_2D[i + 1][j] - u_old->mat_2D[i][j]) / dx;
				v->mat_2D[i][j] -= dt * u_old->mat_2D[i][j] * (v_old->mat_2D[i + 1][j] - v_old->mat_2D[i][j]) / dx;
			}
			if (v->mat_2D[i][j] > 0.0)
			{
				u->mat_2D[i][j] -= dt * v_old->mat_2D[i][j] * (u_old->mat_2D[i][j] - u_old->mat_2D[i][j - 1]) / dy;
				v->mat_2D[i][j] -= dt * v_old->mat_2D[i][j] * (v_old->mat_2D[i][j] - v_old->mat_2D[i][j - 1]) / dy;
			}
			else
			{
				u->mat_2D[i][j] -= dt * v_old->mat_2D[i][j] * (u_old->mat_2D[i][j + 1] - u_old->mat_2D[i][j]) / dy;
				v->mat_2D[i][j] -= dt * v_old->mat_2D[i][j] * (v_old->mat_2D[i][j + 1] - v_old->mat_2D[i][j]) / dy;
			}
		}
  }
}

int idx2rank(int row_idx, int col_idx) {
  return row_idx * num_cols + col_idx;
}

void type_setup(Cmatrix &data) {
  int block_lengths[data.n];
  MPI_Aint displacements[data.n], addresses[data.n], add_start;
  MPI_Datatype typelist[data.n];

  for (int i = 0; i < data.n; i++) {
    block_lengths[i] = 1;
    typelist[i] = MPI_DOUBLE;
    MPI_Get_address(&data.mat_2D[i][0], &addresses[i]);
  }

  MPI_Get_address(&data.mat_2D[0][0], &add_start);
  for (int i = 0; i < data.n; i++) displacements[i] = addresses[i] - add_start;

  MPI_Type_create_struct(data.n, block_lengths, displacements, typelist, &MPI_matcol);
  MPI_Type_commit(&MPI_matcol);
}

void sync_boundary(Cmatrix &data, int num_tag) {
  vector<MPI_Request> send_requests;
  vector<MPI_Request> recv_requests;
  if (col_idx > 0) { // has left
    MPI_Request send_request;
    MPI_Request recv_request;
    int left_rank = idx2rank(row_idx, col_idx - 1);
    // cout << "process " << id << " left rank: " << left_rank << endl;
    MPI_Isend(&data.mat_2D[1][0], data.m, MPI_DOUBLE, left_rank, num_tag, MPI_COMM_WORLD, &send_request);
    MPI_Irecv(&data.mat_2D[0][0], data.m, MPI_DOUBLE, left_rank, num_tag, MPI_COMM_WORLD, &recv_request);
    send_requests.push_back(send_request);
    recv_requests.push_back(recv_request);
  }
  if (col_idx < num_cols - 1) { // has right
    MPI_Request send_request;
    MPI_Request recv_request;
    int right_rank = idx2rank(row_idx, col_idx + 1);
    // cout << "porcess " << id << " right rank: " << right_rank << endl;
    MPI_Isend(&data.mat_2D[data.n - 2][0], data.m, MPI_DOUBLE, right_rank, num_tag, MPI_COMM_WORLD, &send_request);
    MPI_Irecv(&data.mat_2D[data.n - 1][0], data.m, MPI_DOUBLE, right_rank, num_tag, MPI_COMM_WORLD, &recv_request);
    send_requests.push_back(send_request);
    recv_requests.push_back(recv_request);
  }
  if (row_idx > 0) { // has top
    MPI_Request send_request;
    MPI_Request recv_request;
    int top_rank = idx2rank(row_idx - 1, col_idx);
    // cout << "process " << id << " top rank: " << top_rank << endl;
    // because it has top element, so it should send Ny-2 to top and receive Ny-1 from top
    MPI_Isend(&data.mat_2D[0][data.m - 2], 1, MPI_matcol, top_rank, num_tag, MPI_COMM_WORLD, &send_request);
    MPI_Irecv(&data.mat_2D[0][data.m - 1], 1, MPI_matcol, top_rank, num_tag, MPI_COMM_WORLD, &recv_request);
    send_requests.push_back(send_request);
    recv_requests.push_back(recv_request);
  }
  if (row_idx < num_rows - 1) { // has bottom
    MPI_Request send_request;
    MPI_Request recv_request;
    int bottom_rank = idx2rank(row_idx + 1, col_idx);
    // cout << "process " << id << " bottom rank: " << bottom_rank << endl;
    MPI_Isend(&data.mat_2D[0][1], 1, MPI_matcol, bottom_rank, num_tag, MPI_COMM_WORLD, &send_request);
    MPI_Irecv(&data.mat_2D[0][0], 1, MPI_matcol, bottom_rank, num_tag, MPI_COMM_WORLD, &recv_request);
    send_requests.push_back(send_request);
    recv_requests.push_back(recv_request);
  }
  MPI_Waitall(send_requests.size(), &send_requests[0], MPI_STATUSES_IGNORE);
  MPI_Waitall(recv_requests.size(), &recv_requests[0], MPI_STATUSES_IGNORE);

}

void set_pressure_BCs(void)
{
  // cout << "P domain shape: " << col_count << " x " << row_count << endl;
  // cout << "Nx_local: " << Nx_local << " Ny_local: " << Ny_local << endl;
  // cout << "update y interval: " << 1 + (Ny_local - 2) / 2 << " ~ " << Ny_local - 1 << endl;
  // cout << "update x interval: " << i_start - 1 << " ~ " << i_end + 1 << endl;
  if (is_bottom) {
    for (int i = i_start - 1; i < i_end + 1; i++){
      P->mat_2D[i][1] = P->mat_2D[i][2];
    }
  }
  if (is_top) {
    for (int i = i_start - 1; i < i_end + 1; i++){
      P->mat_2D[i][Ny_local - 2] = P->mat_2D[i][Ny_local - 3];
    }
  }
  if (is_right) {
    if (row_idx < (num_rows / 2)) { // bottom half, should update
      for (int j = j_start - 1; j < j_end + 1; j++){
        P->mat_2D[Nx_local - 2][j] = P->mat_2D[Nx_local - 3][j];
      }   
    }
    if ((num_rows % 2 != 0) && (row_idx == int(num_rows / 2))) { // right in the middle, shold update half
      for (int j = 1 + (Ny_local - 2) / 2; j < Ny_local - 1; j++) {
        P->mat_2D[Nx_local - 2][j] = P->mat_2D[Nx_local - 3][j];
      }
    }
  }
}

void set_velocity_BCs(void)
{
  // cout << "Vel domain shape: " << col_count << " x " << row_count << endl;
  // cout << "Nx_local: " << Nx_local << " Ny_local: " << Ny_local << endl;
  // cout << "update y interval: " << 1  << " ~ " << (Ny_local-2) / 2 + 1 << endl;
  // cout << "update x interval: " << j_start - 1 << " ~ " << j_end + 1 << endl;
  if (is_left) {
    for (int j = j_start-1; j < j_end + 1; j++) {
      u->mat_2D[1][j] = u->mat_2D[2][j]; 
    }
  }
  if (is_right) {
    if (row_idx >= ((num_rows + 1) / 2)) { // bottom half, should update
      for (int j = j_start - 1; j < j_end + 1; j++){
        u->mat_2D[Nx_local - 2][j] = u->mat_2D[Nx_local - 3][j];
      }
    }
    if ((num_rows % 2 != 0) && (row_idx == int(num_rows / 2))) { // right in the middle
      for (int j = 1; j < (Ny_local-2) / 2 + 1; j++) {
        u->mat_2D[Nx_local - 2][j] = u->mat_2D[Nx_local - 3][j];
      }
    }
  }
}

double project_velocity(void) {
	double vmax = 0.0;
	for (int i = i_start; i < i_end; i++)
		for (int j = j_start; j < j_end; j++)
		{
			u->mat_2D[i][j] = u->mat_2D[i][j] - dt * (1. / rho) * (P->mat_2D[i + 1][j] - P->mat_2D[i - 1][j]) / (2. * dx);
			v->mat_2D[i][j] = v->mat_2D[i][j] - dt * (1. / rho) * (P->mat_2D[i][j + 1] - P->mat_2D[i][j - 1]) / (2. * dy);

			double vel = sqrt(u->mat_2D[i][j] * u->mat_2D[i][j] + v->mat_2D[i][j] * v->mat_2D[i][j]);

			vmax = max(vmax, vel);
		}
	set_velocity_BCs();
	return vmax;
}

int pressure_poisson_jacobi(double rtol = 1.e-5)
{
	double tol = 10. * rtol; // init tolerrance
  double tol_local; // local tolerrance for commnication
	int it = 0; // iterations took

	while (tol > rtol)
	{
    Cmatrix::swap(*P, *P_old);
		double sum_val = 0.0;
    double sum_val_local = 0.0;
    tol = 0.0;
		tol_local = 0.0;
		it++;
    jacobian_it++;

    chrono::high_resolution_clock::time_point start = chrono::high_resolution_clock::now(); // start timer

		// Jacobi iteration
		for (int i = i_start; i < i_end; i++) {
			for (int j = j_start; j < j_end; j++) {
				P->mat_2D[i][j] = 1.0 / (2.0 + 2.0 * (dx * dx) / (dy * dy)) * (P_old->mat_2D[i + 1][j] + P_old->mat_2D[i - 1][j] + 
					(P_old->mat_2D[i][j + 1] + P_old->mat_2D[i][j - 1]) * (dx * dx) / (dy * dy)
					- (dx * dx) * PPrhs->mat_2D[i][j]);

				sum_val_local += fabs(P->mat_2D[i][j]);
				tol_local += fabs(P->mat_2D[i][j] - P_old->mat_2D[i][j]);
			}
    }
    
    sync_boundary(*P, it); // sync boundary

    set_pressure_BCs();
  
    MPI_Allreduce(&tol_local, &tol, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD); // use MPI_Allreduce to get the maximum tol_local
    MPI_Allreduce(&sum_val_local, &sum_val, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD); // use MPI_Allreduce to get the maximum tol_local
    
    tol = tol / max(1.e-10, sum_val);

    chrono::high_resolution_clock::time_point end = chrono::high_resolution_clock::now();
    chrono::microseconds duration = chrono::duration_cast<chrono::microseconds>(end - start);
    jacobian_time += duration.count();
	}

	return it;
}

void solve_NS() {

  time_it = 0;
  jacobian_it = 0;
  jacobian_time = 0.0;
	double vel_max = 0.0; // the local maximum velocity
	int its; // the number of iterations tooked to get get the altimate solution on this process
	int out_it = 0; // the number of output iterations
	double t_out = dt_out;
  // solver loop, iterate through time to get solution at every time step
  while (t < t_end) {

    // find current dt by communicating with all other process, find the minimum dt above all processors
    double dt_local;
		if (vel_max > 0.0) { // find minima for local dt
			dt_local = min(courant * min(dx, dy) / vel_max, dt_min);
		}
		else {
      dt_local = dt_min;
    }

    MPI_Allreduce(&dt_local, &dt, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD); // collective communication to get the minimum dt

    // cout << "current time step: " << dt << endl;

		t += dt; // total elapsed time increment
		time_it++; // iteration count increment

    Cmatrix::swap(*u, *u_old); // swap vector matries
    Cmatrix::swap(*v, *v_old);

    calculate_intermediate_velocity(); // do intermediate velocity calculation

    sync_boundary(*v, time_it); // sync boundary
    sync_boundary(*u, time_it); // sync boundary
    
    calculate_ppm_RHS_central(); // calculate RHS for pressure poisson equation
    its = pressure_poisson_jacobi(1.e-5); // solve pressure poisson equation

    vel_max = project_velocity(); // project velocity

    sync_boundary(*v, time_it); // sync boundary
    sync_boundary(*u, time_it); // sync boundary

		if (t >= t_out)
		{
      double vel_max_global;
      MPI_Reduce(&vel_max, &vel_max_global, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
			out_it++;
			t_out += dt_out;
      // only in process zero do the priting to make log cleaner
      if (id == 0) {
			  cout << time_it << ": \t" << t << " Jacobi iterations: " << its << " vel_max: " << vel_max_global << " with dt " << dt << endl;
      }
#ifdef OUTPUT
      // file writing
      stringstream Pfname, ufname, vfname;
      Pfname << "../out/" << id << "_" << row_idx << "_" << col_idx << "_P_" << out_it << ".dat";
      ufname << "../out/" << id << "_" << row_idx << "_" << col_idx << "_u_" << out_it << ".dat";
      vfname << "../out/" << id << "_" << row_idx << "_" << col_idx << "_v_" << out_it << ".dat";
      P->to_tab_file(Pfname.str());
      u->to_tab_file(ufname.str());
      v->to_tab_file(vfname.str());
#endif
		}
    // break;
  }
};

void setup() {
  // setup local grids

  Ny_local = row_count + 2;
  Nx_local = col_count + 2;

  P = new Cmatrix(Nx_local, Ny_local); P->fill(0.0);
  P_old = new Cmatrix(Nx_local, Ny_local); P_old->fill(0.0);
  u = new Cmatrix(Nx_local, Ny_local); u->fill(0.0);
  u_old = new Cmatrix(Nx_local, Ny_local); u_old->fill(0.0);
  v = new Cmatrix(Nx_local, Ny_local); v->fill(0.0);
  v_old = new Cmatrix(Nx_local, Ny_local); v_old->fill(0.0);
  PPrhs = new Cmatrix(Nx_local, Ny_local); PPrhs->fill(0.0);

  is_top = (row_idx == 0);
  is_bottom = (row_idx == num_rows - 1);
  is_left = (col_idx == 0);
  is_right = (col_idx == num_cols - 1);

  i_start = 1; i_end = Nx_local - 1;
  j_start = 1; j_end = Ny_local - 1;

  if (is_bottom) j_start = 2;
  if (is_top) j_end = Ny_local - 2;
  if (is_left) i_start = 2;
  if (is_right) i_end = Nx_local - 2;

  type_setup(*P);

	dx = Lx / (Nx - 1);
	dy = Ly / (Ny - 1);

  // init P
  if (col_idx == 0) {
    for (int i = 1; i < Ny_local - 1; i++) {
      P->mat_2D[1][i] = P_max;
    }
  }
  // init P_old to be the same as P
  if (col_idx == 0) {
    for (int i = 1; i < Ny_local - 1; i++) {
      P_old->mat_2D[1][i] = P_max;
    }
  }

	t = 0.0;
};

void cleanup() {
  delete P;
  delete P_old;
  delete u;
  delete u_old;
  delete v;
  delete v_old;
  delete PPrhs;
};

int main(int argc, char* argv[]) {
  chrono::high_resolution_clock::time_point start = chrono::high_resolution_clock::now(); // start timer

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &p);
  MPI_Comm_rank(MPI_COMM_WORLD, &id);

  get_subdomain(Nx, Ny, id, p, &row_start, &col_start, &row_end, &col_end, &num_rows, &num_cols, &row_idx, &col_idx, &row_count, &col_count);
  
  setup();

  if (id == 0) {
    cout << endl << "========================================" << endl;
    cout << "Running with " << p << " processors" << endl;
    cout << "Geometry lay out " << endl;
    cout << "Nx: " << Nx << " Ny: " << Ny << endl;
    cout << "NxLocal: " << Nx_local << " NyLocal: " << Ny_local << endl;
    cout << "num_rows: " << num_rows << " num_cols: " << num_cols << endl;
    cout << "========================================" << endl << endl;
  }

  // cout << "Process " << id << " has " << row_count << " rows and " << col_count << " cols" << endl;

  solve_NS();

  MPI_Finalize();

  cleanup();


  if (id == 0) cout << endl << endl << ";-) programme finished" << endl << endl;

  chrono::high_resolution_clock::time_point end = chrono::high_resolution_clock::now();
  chrono::seconds duration = chrono::duration_cast<chrono::seconds>(end - start);
  chrono::microseconds m_duration = chrono::duration_cast<chrono::microseconds>(end - start);

  // Calculate minutes and seconds
  long minutes = duration.count() / 60;
  long seconds = duration.count() % 60;
  
  if (id == 0) {
    cout << endl << "========================================" << endl;
    cout << "Total time taken: " << minutes << " min : " << seconds << "s" << endl;
    cout << "Total seconds took: " << duration.count() << endl;
    cout << "Time per time step: " << m_duration.count() / (double)time_it << " ms" << endl;
    cout << "Time per jacobian iteration: " << jacobian_time / (double)jacobian_it << " ms" << endl;
    cout << "========================================" << endl << endl;
  }


  return 0;
}