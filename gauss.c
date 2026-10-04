#include "gauss.h"
#include <math.h>

int Gauss(int n, double** A, double** B) {
int i, j, i0=0, k; double M=0, eps=1e-8, bufer;

for (i=0; i<n; i++) { for (j=0; j<n; j++) if (fabs(A[i][j])>M) M=fabs(A[i][j]); }
eps*=M; M=0;                         //choose eps

for (j=0; j<n; j++) {                  //direct moving
M=0;
for (i=j; i<n; i++) { if (fabs(A[i][j])>fabs(M)) {i0=i; M=A[i][j];} }      //max element
if (fabs(M)<eps) return 0;

for (k=j; k<n; k++) { bufer=A[j][k]; A[j][k]=A[i0][k]; A[i0][k]=bufer; }      //swap strings
for (k=0; k<n; k++) { bufer=B[j][k]; B[j][k]=B[i0][k]; B[i0][k]=bufer; }

for (k=j; k<n; k++) A[j][k]/=M;
for (k=0; k<n; k++) B[j][k]/=M;         //arithmetic

for (i=j+1; i<n; i++) { bufer=A[i][j];
for (k=j; k<n; k++) A[i][k]-=bufer*A[j][k];
for (k=0; k<n; k++) B[i][k]-=bufer*B[j][k];
}

}

for (j=n-1; j>0; j--) {
for (i=0; i<j; i++) { bufer=A[i][j];         //inverse moving
for (k=0; k<n; k++) B[i][k]-=bufer*B[j][k];
}
}

return 1;
}
