#ifndef _BMPOISSON_HXX
#define _BMPOISS0N_HXX

// MultiGridFAS constants
const double minerr=5.0e-4; // Relative error in Gauss-Seidel
const int maxit = 1000; // Max iter in Gauss-Seidel
const int NGMAX=10; // Max number og grids
const int cn = 1; // Number of coarser grids

int inklud(int dim, int til);
void GausS_DvN(double **a, double **b, double **rh,
               double (*Te)(int i, int j, int imax, int jmax), double ll[], int nl[]);
void s_GausS_DvN(double **a, double **b, double **rh,
                 double (*Te)(int i, int j, int imax, int jmax), double ll[], int nl[]);
void s1_GausS_DvN(double **a, double **b, double **rh,
                  double (*Te)(int i, int j, int imax, int jmax), double ll[], int nl[]);
void mgfas(double **u, double **r, double **b,
           double (*T)(int i, int j, int imax, int jmax), double lx[], int n[]);
void s_mgfas(double **u, double **r, double **b,
             double (*T)(int i, int j, int imax, int jmax), double lx[], int n[]);
void s1_mgfas(double **u, double **r, double **b,
              double (*T)(int i, int j, int imax, int jmax), double lx[], int n[]);
void interp(double **uf, double **uc, int nf[]);
void rstrct(double **uc, double **uf, int nc[]);
void s_rstrct(double **uc, double **uf, int nc[]);
void s1_rstrct(double **uc, double **uf, int nc[]);
void smooth(double **uc, double **uf, int nc);

#endif
