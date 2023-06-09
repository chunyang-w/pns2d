# include <string>

class Cmatrix
{
public:
	double* mat_1D;
	double** mat_2D;
	int n, m;

	Cmatrix();
	Cmatrix(const Cmatrix& other);
	Cmatrix(int imax, int jmax);
	~Cmatrix();
	void to_tab_file(std::string fname);
	void fill(double val);
	static void swap(Cmatrix &a, Cmatrix &b);
};