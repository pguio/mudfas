#ifdef THREE_D
#ifndef _BMPOISSON3D_HXX
#define _BMPOISS0N3D_HXX


const double minerr= 5.0e-4;
const int NGMAX=10, maxnwt=7, maxcycle = 1, maxit = 3000, cn = 1;
//unsigned int it;
//  double ***irho[NGMAX+1],***irhs[NGMAX+1],***iu[NGMAX+1],***ru[NGMAX+1],
//             **irnd[NGMAX+1],**ibnd[NGMAX+1];


void interp(double ***uf, double ***uc, int nf[]);
void rstrct(double ***uc, double ***uf, int nc[]);
void smooth(double **uc, double **uf, int nc[]);
int inklud(int dim, int til);
void GausS_DvN(double ***a, double **b, double ***rh,  const double Te,
               double ll[], int nl[]);
void mgfas(double ***u, double ***r, double **b, double ***rho0, const double T, double lx[],
           int n[]);

#endif
#endif









