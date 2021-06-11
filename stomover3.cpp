//SURFACE TOMOGRAPHY VERSION 3.0
//
//Created and implemented By Tim Petersen
//the SVD and dynamic memory allocation routines were adapted from
//"Numerical Recipes in C Numerical Recipes in C: The Art of Scientific Computing",
// Second Edition, William H. Press, Saul A. Teukolsky, William T. Vetterling
// and Brian P. Flannery, Cambridge University Press, 2002.
//
//Old C-Style memory and array managment based upon that contained in
//"Advanced computing in electron microscopy".
//By Earl J. Kirkland, New York: Plenum Press, 1998.

//Version 1.0 code demonstration published in CPC:
/*
@article{stomo2010,
  title = {An electron tomography algorithm for reconstructing 3D morphology using surface tangents of projected scattering interfaces},
  journal = {Comput. Phys. Commun.},
  volume = {181},
  number = {3},
  pages = {676-682},
  year = {2010},
  author = {T.C. Petersen and S.P. Ringer}
}
*/

//Version 2.0 udpated and implemented in 2010 by Tim Petersen to cater for
//rectangular image sets, ensure smoother image rotations, provide ridge
//detection (suitable for sensing phase-contrast Fresnel fringes), compute
//faster/larger kernal edge detection and also greatly reduce RAM usage.
//Published in CPC (2012) after lengthy review:
/*@article{stomoCodev2,
  title = {Local electron tomography using angular variations of surface tangents: Stomo version 2},
  journal = {Comput. Phys. Commun.},
  volume = {183},
  number = {3},
  pages = {698-704},
  year = {2012},
  author = {T.C. Petersen and S.P. Ringer}
}
*/

//Version 3.0 further developed at the Centre for Electron Microscopy,
//Monash University, Clayton Campus.
//StomoVer3 of this repository upgrades the ridge and valley detection
//to correctly implement hysteresis threholding and computation of
//the normals at the reconstructed surface-tangent intersection points.
//Perhaps the most useful update is the proper user-defined selection
//of ridges OR valleys, rather than a Hessian that is sensitive to both.
//This is particularly important for phase contrast images, which most
//often have both of these features in equal amounts.  By topoligically
//distinguishing these two, sparsity increases and thus improves tracking.
//Other important fixes include new feature-thinning to remove unuseful
//adjacent points that align orthogonal to the tilt axis and so provide
//no deoth information and otherwise spoil tracking.  This is implemented
//using suppress_edges_horizontal() with calculates a 2nd non-maximum
//suppressiong in the direction of the tilt axis.

//Last updated: June 2021.

#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>
#include <assert.h>

#define PI 3.141592654
#define negflag -1.0e5 //Global garbage value for maps.

//Single Value Decomposition (SVD) declarations from Numerical Recipes:

# define true  ((int)1)
# define false ((int)0)
# define nmax  ((int)1000)
# define mmax  ((int)50)
# define tol   ((float)1.0e-5)

using namespace std;

void funcs(float x, float* afunc, unsigned int ma);
void svdfit(float* X, float* Y, float* Sig, unsigned int NData,
    float* A, unsigned int MA,
    float** U, float** V, float* W, unsigned int MP, unsigned int NP,
    float* ChiSq, void funcs(float x, float* afunc, unsigned int ma));
void svdcmp(float** A, unsigned int M, unsigned int N,
    unsigned int MP, unsigned int NP, float* W, float** V);
void svbksb(float** U, float* W, float** V, unsigned int M,
    unsigned int N, unsigned int MP, unsigned int NP,
    float* B, float* X);
void svdvar(float** V, unsigned int MA, float* W, float** CVM);
void fpoly(float x, float p[], unsigned int np);
void deriv_poly(float*** images, float*** derivs, float*** errors, float* delzs,
    float* sig, int nimages, int N, int M, unsigned int ma);
//.....end declarations for SVD................................

//.....dynamic memory allocation declarations..................
float** float2D(int nx, int ny, const char* message);
int** int2D(int nx, int ny, const char* message);
void free_int2D(int** a, int nx);
void free_float2D(float** a, int nx);
float*** float3D(int na, int nx, int ny, const char* message);
void free_float3D(float*** a, int na, int nx);
//.....end dynamic memory allocation declarations..............

//.....all other declarations.......................................
void read_binary_2d_flarray(int N, int M, const char* file, float** array);
void write_binary_1d_flarray(int start, int end, const char* file, float* array);
void write_binary_2d_flarray(int N, int M, const char* file, float** array);
void write_binary_3d_flarray(int nstart, int nimgs, int N, int M, const char* file, const char* fn, float*** array);
void read_binary_3d_flarray(int nstart, int nend, int N, int M, const char* file, float*** array);
void make_tip_img_II(float** img, float** typos, float** tyneg, float*** sypos,
    float*** syneg, int N, int M, float ang, float var,
    float width, float yoffset, float sabsorb, float tabsorb,
    float* zheights, float* ycenters, float* xcenters,
    float* radii, int nspheres);
void make_sphere_tilt_series_II(float*** imgs, int nimgs, int N, int M,
    int nspheres, float absorb, float* radii,
    float* zheights);
void make_tip_tilt_series_II(float*** imgs, int nimgs, int N, int M,
    float sabsorb, float tabsorb, float width,
    float var, float yoffset, float* zheights,
    float* ycenters, float* xcenters, float* radii,
    int nspheres);
void make_tip_3d(float*** imgs, int nimgs, int N, int M, float width, float var,
    float yoffset, float start_angle, float sabsorb, float tabsorb,
    float* zheights, float* ycenters, float* xcenters,
    float* radii, int nspheres);
void Sobel(float** img, float** Gx, float** Gy, int N, int M);
void Hessian(float** img, float** Gx, float** Gy, float** Gxx,
    float** Gyy, float** Gxy, int N, int M);
void conv_filters(float** img, int flag, int N, int M);
void convolve_Gaussian(float** img, float fw, int N, int M);
void edge_detect_Canny(float** img, float** edgemap, float** gdir, float fw,
    float T1, float T2, int edgeflag, int valridgeflag, int N, int M);
void link_edges(int & rec, int i, int j, float** emap, float T2, float eflag, int N, int M);
void fix_wireframe_edges(float** wireframe, float res, float clip, int N, int M);
void suppress_edges_horizontal(float** edgemap, int N, int M);
void rotate_imgs(float** img, float** slices, int N, int M, float theta);
void find_z_intersection(float** img, float** imgbel, float** imgabv,
    float** zmap, int nslices, int nchan,
    float tiltx, int res, float angabv, float angbel);
void find_z_intersection_II(float*** imgs, float** zmap, float** emap, float* angs,
    float stdrr, float cut, int nslices, int nchan,
    int npnts, int minpnts, int midpnt, int nimgs,
    int norder, int res, float tiltx);
float search_tangent_x(float** img, int res, int iedge,
    int jedge, int nchan, float tiltx);
void fit_poly(float* zs, float* angs, float* derivs,
    float* errors, float* sig, int npnts, unsigned int ma);
void convert_point_cloud_to_tomogram(float* xyz, int nimgs,
    int N, int M, int ndata);
void shell_sort(unsigned long ndata, float* x, float* y, float* z);
void compute_xyz(float** zmap, float** emap, float** gdir, int N, int M,
    int& numpnts, float theta, float tiltx);
void reconstruct_local(int& numpnts, int imgnum, int res, int N, int M,
    int tiltx, int minpnts, int ndata, int norder, int eflag, int vrflag,
    int nimgs, float tiltang, float clip, float fw, float T1,
    float T2, float stdrr, float cut, float colmin,
    float colmax, float rowmin, float rowmax, float* xyz,
    float* angles);
//.....end all other declarations.......................................

//....all function definitions..........................................
void svdfit(float* X, float* Y, float* Sig, unsigned int NData,
    float* A, unsigned int MA, float** U, float** V, float* W,
    unsigned int MP, unsigned int NP, float* ChiSq,
    void funcs(float x, float* afunc, unsigned int ma))
{
    /*
      Given a set of NData points X[], Y[] with individual standard
      deviations of Sig[], use chi-square minimization to determine the
      MA coefficients, A[], of the fitting function
      y = sum over i Ai * funcsi(x).
      Here we solve the fitting equation using singular value decomposition
      of the NData by MA matrix. The arrays U, V and W provide workspace
      on input. On output they define the singular value decomposition and
      can be used to obtaint the covariance matrix. MP and NP are the
      physical dimensions of the matrices U, V, and W as indicated below.
      It is necessary that MP be greater than or equal to NData and that
      NP be greather than or equal to MP. The program returns values for
      the MA fit parameters A[] and the chi-square, ChiSq. The user
      supplies a subroutine, funcs(), that returns the MA basis functions
      evaluated at x in the array afunc[].
    */

    int i, j;
    float sum, thresh, tmp, wmax, wmin;
    float beta[nmax], afunc[mmax];

    // Accumulate coefficients of the fitting matrix.
    for (i = 0; i < NData; ++i) {
        funcs(X[i], afunc, MA);
        tmp = 1.0 / Sig[i];
        for (j = 0; j < MA; ++j) {
            U[i][j] = afunc[j] * tmp;
        }
        beta[i] = Y[i] * tmp;
    }

    /* Singular value decomposition. */
    svdcmp(U, NData, MA, MP, NP, W, V);
    /* Edit the singular values, given tol from the parameter statement,
       between here ... */
    wmax = 0.0e00;
    wmin = 1.0e30;
    for (j = 0; j < MA; ++j) {
        if (W[j] > wmax)
            wmax = W[j];
        if (W[j] < wmin)
            wmin = W[j];
    }

    thresh = tol * wmax;
    for (j = 0; j < MA; ++j) {
        if (W[j] < thresh) {
            W[j] = 0.0;
        }
    }
    /* ... and here. */

    svbksb(U, W, V, NData, MA, MP, NP, beta, A);

    // Evaluate chi-square.
    *ChiSq = 0.0;
    for (i = 0; i < NData; ++i) {
        funcs(X[i], afunc, MA);
        sum = 0.0;
        for (j = 0; j < MA; ++j) {
            sum = sum + A[j] * afunc[j];
        }
        tmp = (Y[i] - sum) / Sig[i];
        *ChiSq = *ChiSq + tmp * tmp;
    }

    return;
}

void svdvar(float** V, unsigned int MA, float* W, float** CVM)
{
    /*
    To evaluate the covariance matrix CVM of the fit for MA paramaters
    obtained by svdfit, call this routine with matrix V and W as returned
    from svdfit.
    */

    int i, j, k;
    float sum, * wti;

    wti = new float[MA];

    for (i = 0; i < MA; ++i) {
        wti[i] = 0.0E00;
        if (W[i] != 0.0E00)
            wti[i] = 1.0E00 / (W[i] * W[i]);
    }

    for (i = 0; i < MA; ++i) {
        for (j = 0; j <= i; ++j) {
            sum = 0.0E00;
            for (k = 0; k < MA; ++k) {
                sum = sum + V[i][k] * V[j][k] * wti[k];
            }
            CVM[i][j] = sum;
            CVM[j][i] = sum;
        }
    }

    delete[] wti;

    return;
}

void svbksb(float** U, float* W, float** V, unsigned int M,
    unsigned int N, unsigned int MP, unsigned int NP,
    float* B, float* X)
{
    /*
      Solves A * X = B for a vector X where A is specified by the arrays
      U, W and V as returned by svdcmp. M and N are the logical dimensions
      of A and will be equal for a square matrices. MP and NP are the
      physical dimensions of A. B is the input right-hand side. X is the
      output solution vector. No input quantities are destroyed, so the
      routine may be called sequentially with different B's. M must be
      greater to N (see svdcmp).
    */

    int i, j;
    float S, tmp[nmax];

    /* Calculate (transpose U ) times B */
    for (j = 0; j < N; ++j) {
        S = 0.0;
        /* Nonzero result only if W[j] is nonzero. */
        if (W[j] != 0.0) {
            for (i = 0; i < M; ++i) {
                S = S + U[i][j] * B[i];
            }
            S = S / W[j];
        }
        tmp[j] = S;
    }

    /* Multiply by V to get answer. */
    for (j = 0; j < N; ++j) {
        S = 0.0E00;
        for (i = 0; i < N; ++i) {
            S = S + V[j][i] * tmp[i];
        }
        X[j] = S;
    }

    return;
}

void svdcmp(float** A, unsigned int M, unsigned int N,
    unsigned int MP, unsigned int NP, float* W, float** V)
{
    /*
    Give a matrix A, with logical dimensions M by N and physical
    dimensions MP by NP, this routine computes its singular value
    decomposition, A = U * W * transpose V. The matrix U replaces
    A on output. The diagonal matrix of singular values, W, is output
    as a vector W. The matrix V (not the transpose of V) is output as
    V. M must be greater or equal to N. If it is smaller then A should
    be filled up to square with zero rows.
    */

    /* Householder reduction to bidiagonal form. */
    int NM, flag, i, its, j, jj, k, l;
    float C, F, G = 0.0, H, S, X, Y, Z, Scale = 0.0, ANorm = 0.0, tmp;
    float rv1[nmax];

    if (M < N) {
        fprintf(stderr, "You must augment A with extra zero rows.\n");
        return;
    }

    for (i = 0; i < N; ++i) {
        l = i + 1;
        rv1[i] = Scale * G;
        G = 0.0;
        S = 0.0;
        Scale = 0.0;
        if (i < M) {
            for (k = i; k < M; ++k) {
                Scale = Scale + fabs(A[k][i]);
            }
            if (Scale != 0.0) {
                for (k = i; k < M; ++k) {
                    A[k][i] = A[k][i] / Scale;
                    S = S + A[k][i] * A[k][i];
                }
                F = A[i][i];
                G = sqrt(S);
                if (F > 0.0) {
                    G = -G;
                }
                H = F * G - S;
                A[i][i] = F - G;
                if (i != (N - 1)) {
                    for (j = l; j < N; ++j) {
                        S = 0.0;
                        for (k = i; k < M; ++k) {
                            S = S + A[k][i] * A[k][j];
                        }
                        F = S / H;
                        for (k = i; k < M; ++k) {
                            A[k][j] = A[k][j] + F * A[k][i];
                        }
                    }
                }
                for (k = i; k < M; ++k) {
                    A[k][i] = Scale * A[k][i];
                }
            }
        }

        W[i] = Scale * G;
        G = 0.0;
        S = 0.0;
        Scale = 0.0;
        if ((i < M) && (i != (N - 1))) {
            for (k = l; k < N; ++k) {
                Scale = Scale + fabs(A[i][k]);
            }
            if (Scale != 0.0) {
                for (k = l; k < N; ++k) {
                    A[i][k] = A[i][k] / Scale;
                    S = S + A[i][k] * A[i][k];
                }
                F = A[i][l];
                G = sqrt(S);
                if (F > 0.0) {
                    G = -G;
                }
                H = F * G - S;
                A[i][l] = F - G;
                for (k = l; k < N; ++k) {
                    rv1[k] = A[i][k] / H;
                }
                if (i != (M - 1)) {
                    for (j = l; j < M; ++j) {
                        S = 0.0;
                        for (k = l; k < N; ++k) {
                            S = S + A[j][k] * A[i][k];
                        }
                        for (k = l; k < N; ++k) {
                            A[j][k] = A[j][k] + S * rv1[k];
                        }
                    }
                }
                for (k = l; k < N; ++k) {
                    A[i][k] = Scale * A[i][k];
                }
            }
        }
        tmp = fabs(W[i]) + fabs(rv1[i]);
        if (tmp > ANorm)
            ANorm = tmp;
    }

    /* Accumulation of right-hand transformations. */
    for (i = N - 1; i >= 0; --i) {
        if (i < (N - 1)) {
            if (G != 0.0) {
                for (j = l; j < N; ++j) {
                    V[j][i] = (A[i][j] / A[i][l]) / G;
                }
                for (j = l; j < N; ++j) {
                    S = 0.0;
                    for (k = l; k < N; ++k) {
                        S = S + A[i][k] * V[k][j];
                    }
                    for (k = l; k < N; ++k) {
                        V[k][j] = V[k][j] + S * V[k][i];
                    }
                }
            }
            for (j = l; j < N; ++j) {
                V[i][j] = 0.0;
                V[j][i] = 0.0;
            }
        }
        V[i][i] = 1.0;
        G = rv1[i];
        l = i;
    }

    /* Accumulation of left-hand transformations. */
    for (i = N - 1; i >= 0; --i) {
        l = i + 1;
        G = W[i];
        if (i < (N - 1)) {
            for (j = l; j < N; ++j) {
                A[i][j] = 0.0;
            }
        }
        if (G != 0.0) {
            G = 1.0 / G;
            if (i != (N - 1)) {
                for (j = l; j < N; ++j) {
                    S = 0.0;
                    for (k = l; k < M; ++k) {
                        S = S + A[k][i] * A[k][j];
                    }
                    F = (S / A[i][i]) * G;
                    for (k = i; k < M; ++k) {
                        A[k][j] = A[k][j] + F * A[k][i];
                    }
                }
            }
            for (j = i; j < M; ++j) {
                A[j][i] = A[j][i] * G;
            }
        }
        else {
            for (j = i; j < M; ++j) {
                A[j][i] = 0.0;
            }
        }
        A[i][i] = A[i][i] + 1.0;
    }

    /* Diagonalization of the bidiagonal form.
       Loop over singular values. */
    for (k = (N - 1); k >= 0; --k) {
        /* Loop over allowed iterations. */
        for (its = 1; its <= 30; ++its) {
            /* Test for splitting.
               Note that rv1[0] is always zero. */
            flag = true;
            for (l = k; l >= 0; --l) {
                NM = l - 1;
                if ((fabs(rv1[l]) + ANorm) == ANorm) {
                    flag = false;
                    break;
                }
                else if ((fabs(W[NM]) + ANorm) == ANorm) {
                    break;
                }
            }

            /* Cancellation of rv1[l], if l > 0; */
            if (flag) {
                C = 0.0;
                S = 1.0;
                for (i = l; i <= k; ++i) {
                    F = S * rv1[i];
                    if ((fabs(F) + ANorm) != ANorm) {
                        G = W[i];
                        H = sqrt(F * F + G * G);
                        W[i] = H;
                        H = 1.0 / H;
                        C = (G * H);
                        S = -(F * H);
                        for (j = 0; j < M; ++j) {
                            Y = A[j][NM];
                            Z = A[j][i];
                            A[j][NM] = (Y * C) + (Z * S);
                            A[j][i] = -(Y * S) + (Z * C);
                        }
                    }
                }
            }
            Z = W[k];
            /* Convergence. */
            if (l == k) {
                /* Singular value is made nonnegative. */
                if (Z < 0.0) {
                    W[k] = -Z;
                    for (j = 0; j < N; ++j) {
                        V[j][k] = -V[j][k];
                    }
                }
                break;
            }
            if (its >= 30) {
                fprintf(stderr, "No convergence in 30 iterations.\n");
                return;
            }

            X = W[l];
            NM = k - 1;
            Y = W[NM];
            G = rv1[NM];
            H = rv1[k];
            F = ((Y - Z) * (Y + Z) + (G - H) * (G + H)) / (2.0 * H * Y);
            G = sqrt(F * F + 1.0);
            tmp = G;
            if (F < 0.0)
                tmp = -tmp;
            F = ((X - Z) * (X + Z) + H * ((Y / (F + tmp)) - H)) / X;
            /* Next QR transformation. */
            C = 1.0;
            S = 1.0;
            for (j = l; j <= NM; ++j) {
                i = j + 1;
                G = rv1[i];
                Y = W[i];
                H = S * G;
                G = C * G;
                Z = sqrt(F * F + H * H);
                rv1[j] = Z;
                C = F / Z;
                S = H / Z;
                F = (X * C) + (G * S);
                G = -(X * S) + (G * C);
                H = Y * S;
                Y = Y * C;
                for (jj = 0; jj < N; ++jj) {
                    X = V[jj][j];
                    Z = V[jj][i];
                    V[jj][j] = (X * C) + (Z * S);
                    V[jj][i] = -(X * S) + (Z * C);
                }
                Z = sqrt(F * F + H * H);
                W[j] = Z;

                /* Rotation can be arbitrary if Z = 0. */
                if (Z != 0.0) {
                    Z = 1.0 / Z;
                    C = F * Z;
                    S = H * Z;
                }
                F = (C * G) + (S * Y);
                X = -(S * G) + (C * Y);
                for (jj = 0; jj < M; ++jj) {
                    Y = A[jj][j];
                    Z = A[jj][i];
                    A[jj][j] = (Y * C) + (Z * S);
                    A[jj][i] = -(Y * S) + (Z * C);
                }
            }
            rv1[l] = 0.0;
            rv1[k] = F;
            W[k] = X;
        }
    }

    return;
}

void fpoly(float x, float p[], unsigned int np)
{
    //Function supplied to NR SVD fitting code for polynomial fitting.
    //The polynomial coefficients are for an assumed Talyor series variation.

    int i, j;
    float nfact;

    p[0] = 1.0E00;

    for (j = 1; j < np; j++) {
        nfact = 1.0E00;
        for (i = j + 1; i > 1; i--) {
            nfact *= (i - 1.0E00);
        }
        p[j] = p[j - 1] * x / nfact;
    }

}

int** int2D(int nx, int ny, const char* message)
{
    //	2D array allocator for type int
    //	make space for m[0...(nx-1)][0..(ny-1)]

    int** m;
    int i;

    m = (int**)malloc(nx * sizeof(int*));
    if (m == NULL) {
        printf("int2D cannot allocate pointers, size=%d: %s\n", nx, message);
        exit(0);
    }

    for (i = 0; i < nx; i++) {
        m[i] = (int*)malloc(ny * sizeof(int));
        if (m[i] == NULL) {
            printf("int2D cannot allocate arrays, size=%d: %s\n", ny, message);
            exit(0);
        }
    }

    return m;

}  /* end int2D() */

float** float2D(int nx, int ny, const char* message)
{
    //	2D array allocator for type float
    //	make space for m[0...(nx-1)][0..(ny-1)]

    float** m;
    int i;

    m = (float**)malloc(nx * sizeof(float*));
    if (m == NULL) {
        printf("float2D cannot allocate pointers, size=%d: %s\n", nx, message);
        exit(0);
    }

    for (i = 0; i < nx; i++) {
        m[i] = (float*)malloc(ny * sizeof(float));
        if (m[i] == NULL) {
            printf("float2D cannot allocate arrays, size=%d: %s\n", ny, message);
            exit(0);
        }
    }

    return m;

}  /* end float2D() */

float*** float3D(int na, int nx, int ny, const char* message)
{
    //	3D array allocator for type float
    //	make space for m[0..(na-1)][0...(nx-1)][0..(ny-1)]

    float*** a;
    int i, j;

    a = (float***)malloc(na * sizeof(float***));
    if (a == NULL) {
        printf("float3D cannot allocate pointers, size=%d: %s\n",
            na, message);
        exit(0);
    }

    for (i = 0; i < na; i++) {
        a[i] = (float**)malloc(nx * sizeof(float*));
        if (a[i] == NULL) {
            printf("float3D cannot allocate arrays, size=%d: %s\n",
                nx, message);
            exit(0);
        }
        for (j = 0; j < nx; j++) {
            a[i][j] = (float*)malloc(ny * sizeof(float));
            if (a[i][j] == NULL) {
                printf("float3D cannot allocate arrays, size=%d: %s\n",
                    ny, message);
                exit(0);
            }
        }
    }

    return a;
}  /* end float3D() */

void free_float2D(float** a, int nx)
{
    int i;

    for (i = nx - 1; i >= 0; i--) {
        free(a[i]);
    }
    free(a);

}  /* end free_float2D() */

void free_int2D(int** a, int nx)
{
    int i;

    for (i = nx - 1; i >= 0; i--) {
        free(a[i]);
    }
    free(a);

}  /* end free_int2D() */


void free_float3D(float*** a, int na, int nx)
{
    int m, i;

    for (m = na - 1; m >= 0; m--) {
        for (i = nx - 1; i >= 0; i--) {
            free(a[m][i]);
        }
        free(a[m]);
    }

    free(a);

}  /* end free_float3D() */

void read_binary_2d_flarray(int N, int M, const char* file, float** array)
{
    // Read floating point file into a 2D array of type float
    int i, j, count, size;
    float* buffer;

    size = N * M;
    buffer = new float[size];

    //Open the binary file with appropriate mode bits (ios flags)
    ifstream fin;
    fin.open(file, ios::in | ios::binary);
    //...........................................................

    //Read the entire file into the buffer array
    fin.read((char*)buffer, size * sizeof(float));
    fin.close();
    //...........................................

    //Convert buffer array to 2D array
    count = 0;
    for (i = 0; i < N; i++) {
        for (j = 0; j < M; j++) {
            array[i][j] = buffer[count];
            count++;
        }
    }
    //................................
    delete[] buffer;
}

// Write a 1D array of type float
void write_binary_1d_flarray(int start, int end, const char* file, float* array)
{
    int size, i;
    size = end - start;

    //Open the binary file with appropriate mode bits (ios flags)
    ofstream fout;
    if (start == 0)
        fout.open(file, ios::out | ios::binary);
    else
        fout.open(file, ios::app | ios::binary);
    //...........................................................

    fout.write((char*)array, size * sizeof(float));
    ////  fout.close();
    //NOTE: bug here.  Must for some reason NOT use ofstream close() but let the
    //implicit destructor close the file once the function is out of scope.
    //Perhaps this is due to I/O Windows scheduling issues when writing to
    //the disk too often.
  //...........................................
}

void write_binary_2d_flarray(int N, int M, const char* file, float** array)
{
    // Write a *.dat floating point file of a 2D array of type float
    int i, j;
    unsigned long count, size;
    float* buffer;

    size = N * M;
    buffer = new float[size];

    count = 0;
    for (i = 0; i < N; i++) {
        for (j = 0; j < M; j++) {
            buffer[count] = array[i][j];
            count++;
        }
    }

    //Open the binary file with appropriate mode bits (ios flags)
    ofstream fout;
    fout.open(file, ios::out | ios::binary);
    //...........................................................

    //Read the entire file into the buffer array
    fout.write((char*)buffer, size * sizeof(float));
    fout.close();
    //...........................................
    delete[] buffer;
}

void write_binary_3d_flarray(int nstart, int nend, int N, int M,
    const char* file, const char* fn, float*** array)
{
    // Write to a *.dat floating point file, a 3D image stack of type float
    int i, j, k, nimgs;
    unsigned long count, size;
    float* buffer;

    ////cout<<"\n\nBeginning file write for "<<fn<<" ..."<<endl;

    nimgs = nend - nstart;
    size = N * M * nimgs;
    buffer = new float[size];

    //Read the data into the 1D buffer array
    count = 0;
    for (k = 0; k < nimgs; k++) {
        for (i = 0; i < N; i++) {
            for (j = 0; j < M; j++) {
                buffer[count] = array[k][i][j];
                count++;
            }
        }
    }

    //Open the binary file with appropriate mode bits (ios flags)
    ofstream fout;
    if (nstart == 0)
        fout.open(file, ios::out | ios::trunc | ios::binary);
    else
        fout.open(file, ios::app | ios::out | ios::binary);
    //...........................................................

    //write the buffer array to the file
    fout.write((char*)buffer, size * sizeof(float));
    ////fout.close();
   //NOTE: bug here.  Must for some reason NOT use ofstream close() but let the
   //implicit destructor close the file once the function is out of scope.
   //Perhaps this is due to I/O Windows scheduling issues when writing to
   //the disk too often.
 //...........................................
    delete[] buffer;
}

void read_binary_3d_flarray(int nstart, int nend, int N, int M, const char* file, float*** array)
{
    // Read a *.dat floating point file, containg a 3D image stack of type float
    int i, j, m, count, size;
    float* buffer;

    size = N * M * (nend - nstart) * sizeof(float);
    buffer = new float[size];
    //Open the binary file with appropriate mode bits (ios flags)
    ifstream fin;
    fin.open(file, ios::in | ios::binary);
    // get length of file:
    fin.seekg(nstart * N * M * sizeof(float), ios::beg);

    //...........................................................

    //Read the entire file into the buffer array
    fin.read((char*)buffer, size * sizeof(float));
    fin.close();
    //...........................................

    //Convert buffer array to 3D array
    count = 0;
    for (m = 0; m < nend - nstart; m++) {
        for (i = 0; i < N; i++) {
            for (j = 0; j < M; j++) {
                array[m][i][j] = buffer[count];
                count++;
            }
        }
    }
    //................................
    delete[] buffer;
}

void Sobel(float** img, float** Gx, float** Gy, int N, int M)
{
    //Function computes arrays required to the discrete Sobel transform
    //of img and stores the results in the (x,y) gradient maps Gx, Gy.

    int i, j;

    for (i = 1; i < N - 1; i++) {
        for (j = 1; j < M - 1; j++) {
            //Sobel transform calculated using discrete convolution
            Gy[i][j] += 1 * img[i - 1][j - 1]; Gy[i][j] += 2 * img[i][j - 1]; Gy[i][j] += 1 * img[i + 1][j - 1];
            Gy[i][j] += 0;
            Gy[i][j] -= 1 * img[i - 1][j + 1]; Gy[i][j] -= 2 * img[i][j + 1]; Gy[i][j] -= 1 * img[i + 1][j + 1];

            Gx[i][j] += 1 * img[i - 1][j - 1]; Gx[i][j] += 0; Gx[i][j] -= 1 * img[i + 1][j - 1];
            Gx[i][j] += 2 * img[i - 1][j];   Gx[i][j] += 0; Gx[i][j] -= 2 * img[i + 1][j];
            Gx[i][j] += 1 * img[i - 1][j + 1]; Gx[i][j] += 0; Gx[i][j] -= 1 * img[i + 1][j + 1];
        }
    }

}

void Hessian(float** img, float** Gx, float** Gy, float** Gxx, float** Gyy, float** Gxy, int N, int M)
{
    //Function computes the Hessian matrix transform of img and stores
    //the results in the (x,y) 2nd order gradient maps Gxx, Gyy and Gxy.

    Sobel(Gx, Gxx, Gxy, N, M);
    Sobel(Gy, Gxy, Gyy, N, M);
}

void Hessian_finite(float** img, float** Gx, float** Gy, float** Gxx, float** Gyy, float** Gxy, int N, int M)
{
    //Function computes the Hessian matrix transform of img and stores
    //the results in the (x,y) 2nd order gradient maps Gxx, Gyy and Gxy.

    int i, j;

    for (i = 1; i < N - 1; i++) {
        for (j = 1; j < M - 1; j++) {
            //Use a simple centred finite-difference for gradients
            Gx[i][j] = (img[i][j+1] - img[i][j-1]) / 2.0;
            Gy[i][j] = (img[i + 1][j] - img[i - 1][j]) / 2.0;
        }
    }

    for (i = 1; i < N - 1; i++) {
        for (j = 1; j < M - 1; j++) {
            //Use a simple centred finite-difference again for 2nd order gradients
            Gxx[i][j] = (Gx[i][j + 1] - Gx[i][j - 1]) / 2;
            Gyy[i][j] = (Gy[i + 1][j] - Gy[i - 1][j]) / 2;
            Gxy[i][j] = (Gx[i+1][j] - Gx[i-1][j])/2;
        }
    }

}


void conv_filters(float** img, int flag, int N, int M)
{
    //Function computes a collection of different convolution filters
    //selected by flag and applies them to the image img.

    int i, j;
    float** timg, ** timg2;

    timg = float2D(N, M, "timg");
    timg2 = float2D(N, M, "timg2");

    for (i = 0; i < N; i++) {
        for (j = 0; j < M; j++) {
            timg[i][j] = 0.0E00;
            timg2[i][j] = 0.0E00;
        }
    }

    if (flag == 1) {
        //Gaussian smoothing with 3x3 Kernal
        for (i = 1; i < N - 1; i++) {
            for (j = 1; j < M - 1; j++) {
                timg[i][j] += 1 * img[i - 1][j - 1]; timg[i][j] += 2 * img[i][j - 1]; timg[i][j] += 1 * img[i + 1][j - 1];
                timg[i][j] += 2 * img[i - 1][j];   timg[i][j] += 4 * img[i][j];   timg[i][j] += 2 * img[i + 1][j];
                timg[i][j] += 1 * img[i - 1][j + 1]; timg[i][j] += 2 * img[i][j + 1]; timg[i][j] += 1 * img[i + 1][j + 1];
                timg[i][j] /= 16.0E00;
            }
        }
    }

    if (flag == 2) {
        for (i = 1; i < N - 1; i++) {
            for (j = 1; j < M - 1; j++) {
                //Discrete Laplacian calculated using discrete convolution
                timg[i][j] += -1 * img[i - 1][j - 1]; timg[i][j] += -1 * img[i][j - 1]; timg[i][j] += -1 * img[i + 1][j - 1];
                timg[i][j] += -1 * img[i - 1][j];   timg[i][j] += 8 * img[i][j];   timg[i][j] += -1 * img[i + 1][j];
                timg[i][j] += -1 * img[i - 1][j + 1]; timg[i][j] += -1 * img[i][j + 1]; timg[i][j] += -1 * img[i + 1][j + 1];
            }
        }
    }

    if (flag == 3) {
        //Sobel transform calculated using discrete convolution
        Sobel(img, timg, timg2, N, M);
        for (i = 0; i < N; i++) {
            for (j = 0; j < M; j++) {
                timg[i][j] = sqrt(pow(timg[i][j], 2) + pow(timg2[i][j], 2));
            }
        }
    }

    if (flag == 4) {
        //x-dir 1D Gaussian smoothing with 1x3 Kernal
        for (i = 1; i < N - 1; i++) {
            for (j = 1; j < M - 1; j++) {
                timg[i][j] += 2 * img[i][j - 1];
                timg[i][j] += 4 * img[i][j];
                timg[i][j] += 2 * img[i][j + 1];
                timg[i][j] /= 8.0E00;
            }
        }
    }

    if (flag == 5) {
        //y-dir 1D Gaussian smoothing with 1x3 Kernal
        for (i = 1; i < N - 1; i++) {
            for (j = 1; j < M - 1; j++) {
                timg[i][j] += 2 * img[i - 1][j];
                timg[i][j] += 4 * img[i][j];
                timg[i][j] += 2 * img[i + 1][j];
                timg[i][j] /= 8.0E00;
            }
        }
    }

    if (flag == 1 || flag == 4 || flag == 5) {
        for (i = 0; i < N; i++) {
            for (j = 0; j < 1; j++) {
                timg[i][j] = img[i][j];
            }
        }

        for (i = 0; i < N; i++) {
            for (j = M - 1; j < M; j++) {
                timg[i][j] = img[i][j];
            }
        }

        for (i = 0; i < 1; i++) {
            for (j = 0; j < M; j++) {
                timg[i][j] = img[i][j];
            }
        }

        for (i = N - 1; i < N; i++) {
            for (j = 0; j < M; j++) {
                timg[i][j] = img[i][j];
            }
        }
    }

    for (i = 0; i < N; i++) {
        for (j = 0; j < M; j++) {
            img[i][j] = timg[i][j];
        }
    }

    free_float2D(timg, N);
    free_float2D(timg2, N);

}

void convolve_Gaussian(float** img, float fw, int N, int M)
{
    //Function convloves an image using a 2D discrete Gaussian
    //filter of del+1 * del+1 (rectangular) dimensions.
    //The FWHM of the Gaussian is used to set the kernal size
    //so that the Gaussian is sufficiently sampled without too much
    //clipping.  Recursive filters could be used to improve computational
    //speed here if desired.

    float** timg, ** timg2;
    float** gauss;
    float area, temp;
    int i, j, m, n, del;

    //Kernal size set to 3 times the gaussian fwhm + 1
    del = 3 * int(fw);

    timg = float2D(N, M, "timg");
    timg2 = float2D(N, M, "timg2");
    gauss = float2D(del + 1, del + 1, "gauss");

    area = 0;
    for (n = -del; n <= del; n++) {
        for (m = -del; m <= del; m++) {
            gauss[abs(n)][abs(m)] = exp(-4.0E00 * log(2.0E00) / fw / fw * (n * n + m * m));//area;
            area += gauss[abs(n)][abs(m)];
            //Symmetry of the Gaussian exploited here with abs()
            //Could half this double sum with (0==>del) but it looks nicer this way
        }
    }

    for (i = 0; i < N; i++) {
        for (j = 0; j < M; j++) {
            timg[i][j] = 0.0E00;
            timg2[i][j] = 0.0E00;
        }
    }

    n = 0;
    for (i = del; i < N - del; i++) {
        for (j = del; j < M - del; j++) {
            temp = 0.0E00;
            for (m = -del; m <= del; m++) {
                temp += gauss[abs(n)][abs(m)] * img[i + n][j + m];
                //area norm assures that the average timg intensity is unchanged
            }
            timg[i][j] = temp;
        }
    }

    m = 0;
    for (i = del; i < N - del; i++) {
        for (j = del; j < M - del; j++) {
            temp = 0.0E00;
            for (n = -del; n <= del; n++) {
                temp += gauss[abs(n)][abs(m)] * timg[i + n][j + m];
                //area norm assures that the average timg intensity is unchanged
            }
            timg2[i][j] = temp;
        }
    }

    for (i = del; i < N - del; i++) {
        for (j = del; j < M - del; j++) {
            img[i][j] = timg2[i][j] / area;
        }
    }

    free_float2D(timg, N);
    free_float2D(timg2, N);
    free_float2D(gauss, del + 1);
}

void edge_detect_Canny(float** img, float** edgemap, float** gdir, float fw,
    float T1, float T2, int edgeflag, int valridgeflag, int N, int M)
{
    //Function creates a wireframe using edge detection.  Parameter thresh sets
    //the threshold for detecting an edge from a laplacian convolution filter.
    //Parameter del rejects points near the edges of images, where the laplacian
    //cannot be correctly estimated due to unknown boundary conditions.
    //Hessian function selected by edgeflag also allows 'ridge' detection as
    //opposed to conventional edge detection.  The ridges are detected using
    //2nd partial image derivatives in the Hessian matrix, the eigenvalues of which
    //also provide ridge direction for Canny-type hysteresis thresholding.

    int i, j, m, n, loc, rec;
    int** gdiri;
    float reg = 1e-10, dir, eflag = 10000, temp = 0, Lpp, Lqq, cosbeta, sinbeta, eig1, eig2;
    float** tedge, ** Gx, ** Gxx, ** Gyy, ** Gxy, ** Gy, ** grad, ** timg;

    Gx = float2D(N, M, "Gx");
    Gy = float2D(N, M, "Gy");
    Gxx = float2D(N, M, "Gxx");
    Gxy = float2D(N, M, "Gxy");
    Gyy = float2D(N, M, "Gyy");
    grad = float2D(N, M, "grad");
    gdiri = int2D(N, M, "gdiri");
    timg = float2D(N, M, "timg");

    //Copy image to temp timg[] and zero other arrays
    for (i = 0; i < N; i++) {
        for (j = 0; j < M; j++) {
            timg[i][j] = img[i][j];
            edgemap[i][j] = 0.0E00;
            Gx[i][j] = 0.0E00;
            Gy[i][j] = 0.0E00;
            Gxx[i][j] = 0.0E00;
            Gyy[i][j] = 0.0E00;
            Gxy[i][j] = 0.0E00;
            grad[i][j] = 0.0E00;
            gdir[i][j] = 0.0E00;
            gdiri[i][j] = 0;
        }
    }

    convolve_Gaussian(timg, fw, N, M);

    //Calculate x&y gradients of smoothed image using 1st order finite difference
    Sobel(timg, Gx, Gy, N, M);

    if (edgeflag == 1)
      Hessian(timg, Gx, Gy, Gxx, Gyy, Gxy, N, M);

    //Calculate approx gradient magnitude and direction (Sobel operator)
    if (edgeflag == 0) {
        for (i = 0; i < N; i++) {
            for (j = 0; j < M; j++) {
                grad[i][j] = sqrt(pow(Gx[i][j], 2) + pow(Gy[i][j], 2));
                gdir[i][j] = atan(Gy[i][j] * Gx[i][j] / (Gx[i][j] * Gx[i][j] + reg));
                gdiri[i][j] = int((gdir[i][j] / PI * 180 + 90 + 45.0E00 / 2) / 45) * 45;
                //Add PI to gdir for computation of normals outside of edge_detect
                //gdiri must not be similarly updated here because there will then
                //exist angles outside of the range 0 == > PI and a access violation
                //will occur for the non maximum supression which assumes this range.
                if (Gx[i][j] < 0) gdir[i][j] += PI;
            }
        }
    }
    else {
        for (i = 0; i < N; i++) {
            for (j = 0; j < M; j++) {
                //Gradient ridge detector here chosen as sqrt of equation (50) in
                //T. Lindeberg International Journal of Computer Vision 30(2),
                //117–154 (1998) (see figure 11 in that paper).  These are in fact the
                //eigenvalues of the Hessian, on the next two lines
                Lpp = Gxx[i][j] + Gyy[i][j] - sqrt(pow((Gxx[i][j] - Gyy[i][j]), 2) + 4 * Gxy[i][j] * Gxy[i][j]);
                Lqq = Gxx[i][j] + Gyy[i][j] + sqrt(pow((Gxx[i][j] - Gyy[i][j]), 2) + 4 * Gxy[i][j] * Gxy[i][j]);

                if (valridgeflag&&((abs(Lpp) > abs(Lqq)&&Lpp<0)||(abs(Lqq) > abs(Lpp) && Lqq < 0))||(!valridgeflag && !((abs(Lpp) > abs(Lqq) && Lpp < 0) || (abs(Lqq) > abs(Lpp) && Lqq < 0)))) {
                        grad[i][j] = fw * fw * fw * (Gxx[i][j] + Gyy[i][j]) * (Gxx[i][j] + Gyy[i][j])
                        * ((Gxx[i][j] - Gyy[i][j]) * (Gxx[i][j] - Gyy[i][j])
                            + 4 * Gxy[i][j] * Gxy[i][j]);
                    grad[i][j] = sqrt(grad[i][j]);

                    /* TP 2020: fix ridge linking by using correct Eq. 37 and Eq.38 from Lindberg
                    */
                    sinbeta = Gxy[i][j] / (fabs(Gxy[i][j]) + reg)
                        *sqrt(1.0 / 2.0 * (1.0 - (Gxx[i][j] - Gyy[i][j])
                            / sqrt(reg + (Gxx[i][j] - Gyy[i][j]) * (Gxx[i][j] - Gyy[i][j])
                                + 4 * Gxy[i][j] * Gxy[i][j])));
                    cosbeta =
                        sqrt(1.0 / 2.0 * (1.0 + (Gxx[i][j] - Gyy[i][j])
                            / sqrt(reg + (Gxx[i][j] - Gyy[i][j]) * (Gxx[i][j] - Gyy[i][j])
                                + 4 * Gxy[i][j] * Gxy[i][j])));

                    //TP July 2020: define Hessian direction through the arc tangent
                    //switch direction sign for Ridges or Valleys using valridgeflag.
                    //The switching below shouldn't be necessary, so there must be something wrong
                    //with my understanding of Lindberg's definition of the local Hessian
                    //eigenvalue direction that is parameterised by the angle beta.

                    if(valridgeflag&&(atan2(sinbeta, cosbeta)<0))
                      gdir[i][j] = atan2(sinbeta,cosbeta) + PI/2;

                    if (valridgeflag && (atan2(sinbeta, cosbeta)>= 0))
                        gdir[i][j] = atan2(sinbeta, cosbeta) - PI / 2;


                    if(!valridgeflag)
                      gdir[i][j] = atan2(sinbeta, cosbeta);

                    gdiri[i][j] = int((gdir[i][j] / PI * 180 + 90 + 45.0E00 / 2) / 45) * 45;
                    //TP 2021: this line below is a fudge to fix the sense of normals.
                    //I currently can't work out how to do this properly from the
                    //Hessian eigenvalues, even though there ought to be a genuine 180 degree
                    //ambiguity.  The extra issue stems from the unorthodox signs of angles
                    //that I used to code up the original hysteresis thresholding.
                    //Clearly, this is sub-optimal for ridges and valleys, where the local
                    //gradient ought to be close to zero.
                    if (Gx[i][j] < 0) gdir[i][j] += PI;
                    if (!valridgeflag&&(Gy[i][j] < 0)) gdir[i][j] += PI;

                }
            }
        }
    }

    loc = 2;

    for (i = loc; i < N - loc; i++) {
        for (j = loc; j < M - loc; j++) {
            if (gdiri[i][j] == 0 || gdiri[i][j] == 180) {
                m = i + 0;
                n = j + 1;
            }
            if (gdiri[i][j] == 45) {
                m = i - 1;
                n = j + 1;
            }
            if (gdiri[i][j] == 90) {
                m = i + 1;
                n = j + 0;
            }
            if (gdiri[i][j] == 135) {
                m = i + 1;
                n = j + 1;
            }
            //apply non-maximum supression by searching orthogonal to the edge

            if ((grad[i][j] > grad[m][n]) && (grad[i][j] > grad[2 * i - m][2 * j - n])) {
                edgemap[i][j] = grad[i][j];
            }
            else {
                edgemap[i][j] = 0;
            }

        }
    }

    //Use timg to create proxy for edgemap, for use in hysteresis thresholding.
    //This way, seeds points can be marked in timg whilst preserving the edge
    //strengths in edgemap, which will be output to the user after using timg
    //as a binary filter on edgemap.
    for (i = 0; i < N; i++) {
        for (j = 0; j < M; j++) {
            timg[i][j] = edgemap[i][j];
        }
    }

    //eflag is used to mark seeds in the edgemap from which to search for
    //linked edges with strengths between T1 and T2.  eflag is much larger
    //than either T1 and T2 and so is used as a (floating point) binary flag.

    //Hysteresis thresholding applied here using explicit code from "Feature
    //Extraction and Image Processing", pp 117-119,
    //M. S. Nixon and A. S. Aguado, Newnes (Oxford), 1st edition, 2002.
    //ISBN 0 750650788.

    for (i = loc; i < N - loc; i++) {
        for (j = loc; j < M - loc; j++) {
            rec = 0;
            if (timg[i][j] > T1 && fabs(timg[i][j] - eflag) > 1.0) {
                timg[i][j] = eflag;
                link_edges(rec, i, j, timg, T2, eflag, N, M);
            }
        }
    }

    for (i = 0; i < N; i++) {
        for (j = 0; j < M; j++) {
            if (fabs(timg[i][j] - eflag) > 1 || edgemap[i][j] <= T2)
                edgemap[i][j] = 0.0E00;
        }
    }

    free_float2D(Gx, N);
    free_float2D(Gy, N);
    free_float2D(Gxx, N);
    free_float2D(Gxy, N);
    free_float2D(Gyy, N);
    free_float2D(grad, N);
    free_int2D(gdiri, N);
    free_float2D(timg, N);

}

void link_edges(int & rec, int i, int j, float** emap, float T2, float eflag, int N, int M)
{
    //function links edges above threshold T2 with adjacnent edges larger than T1.
    //rec limits the number of recursions to avoid stack overflow.
    int n, m;

    for (n = i - 1; n <= i + 1; n++) {
        for (m = j - 1; m <= j + 1; m++) {
            if (n >= 0 && n < N && m >= 0 && m < M) {
                if (fabs(emap[n][m] - eflag) > 1.0 && emap[n][m] > T2&&rec<50) {
                    rec = rec + 1;
                    emap[n][m] = eflag;
                    link_edges(rec, n, m, emap, T2, eflag, N, M);
                }
            }
        }
    }
}

void fix_wireframe_edges(float** wireframe, float res, float clip, int N, int M)
{
    //Function removes contiguous strips of data (along rows othog. to tilt axis)
    //that could otherwise cause errors with z intersection reconstruction

    int i, j, m;
    float av, ** fixedwire;

    fixedwire = float2D(N, M, "fixedwire");

    for (i = 0; i < N; i++) {
        for (j = 0; j < M; j++) {
            fixedwire[i][j] = 0.0E00;
        }
    }

    //remove contiguous data along rows othogonal to projected tilt axis
    for (i = 0; i < N; i++) {
        for (j = 0; j < M; j++) {
            av = 0.0E00;
            for (m = j - int(res / 2); m < j + int(res / 2); m++) {
                if (m >= 0 && m < M && wireframe[i][j]>1e-2) {
                    if (wireframe[i][m] > 1e-2) {
                        av += 1.0E00;
                    }
                }
            }
            if (av <= res * fabs(1.0E00 - clip)) {
                fixedwire[i][j] = wireframe[i][j];
            }
            else {
                fixedwire[i][j] = 0.0E00;
            }
        }
    }

    for (i = 0; i < N; i++) {
        for (j = 0; j < M; j++) {
            wireframe[i][j] = fixedwire[i][j];
        }
    }

    free_float2D(fixedwire, N);

}

void suppress_edges_horizontal(float** edgemap, int N, int M)
{
    //Function removes applies non-maxium suppression along the horizontal direcition,
    //orthogonal to the tilt axis, to amelerioate tracking errors.  The input is
    //assumed to be the egde maps which have been linked with hysteresis thresholding,
    //with greyscales given by the edge strength.

    int i, j;
    float ** edge_suppress;

    edge_suppress = float2D(N, M, "edge_suppress");

    for (i = 0; i < N; i++) {
        for (j = 0; j < M; j++) {
            edge_suppress[i][j] = 0.0E00;
        }
    }

    for (i = 0; i < N ; i++) {
        for (j = 1; j < M-1; j++) {
            if ((edgemap[i][j] > edgemap[i][j+1]) && (edgemap[i][j] > edgemap[i][j-1])) {
                edge_suppress[i][j] = edgemap[i][j];
            }
        }
    }

    for (i = 0; i < N; i++) {
        for (j = 0; j < M; j++) {
            edgemap[i][j] = edge_suppress[i][j];
        }
    }


    free_float2D(edge_suppress, N);

}


void rotate_imgs(float** img, float** slices, int N, int M, float theta)
{
    //Function rotates a rectangular image and copies data to an enlarged
    //image with the same pixel sampling.  Bi-linear interpolation is employed
    //to reduce digital artefacts in the rotation transformation, resulting
    //in smoother edges after rotation of img.

    int i, j, Nl, Ml, il, ir, jt, jb;
    float** timg, ish, jsh, temp, irot, jrot;

    theta = theta / 360.0E00 * 2.0E00 * PI;

    //Calculate the new size of the post-rotated rectangular matrix
    Nl = int(N * cos(fabs(theta)) + M * sin(fabs(theta)));
    Ml = int(M * cos(fabs(theta)) + N * sin(fabs(theta)));

    timg = float2D(Nl, Ml, "timg");

    //Center the smaller rectangle inside the bigger one
    //in which it will be inscribed after rotation.
    for (i = 0; i < Nl; i++) {
        for (j = 0; j < Ml; j++) {
            ish = i - (Nl / 2.0 - N / 2.0);
            jsh = j - (Ml / 2.0 - M / 2.0);
            if (int(ish) > 0 && int(ish) < N && int(jsh) > 0 && int(jsh) < M)
                timg[i][j] = img[int(ish)][int(jsh)];
            else
                timg[i][j] = 0.0E00;
        }
    }

    //Rotate the big rectangle containing the centered small rectangle
    //Perform bilinear interpolation to reduce blockiness in rotated image.
    //ir = iright, jt = jtop, jb = jbottom .... etc.
    for (i = 0; i < Nl; i++) {
        for (j = 0; j < Ml; j++) {
            ish = i - Nl / 2.0;
            jsh = j - Ml / 2.0;
            irot = cos(theta) * (ish)+sin(theta) * jsh + Nl / 2.0 + 0.5;
            jrot = -sin(theta) * (ish)+cos(theta) * jsh + Ml / 2.0 + 0.5;
            if (int(irot - 0.5) > 0 && int(jrot - 0.5) > 0 && int(irot + 0.5) < Nl && int(jrot + 0.5) < Ml) {
                il = int(irot + 0.5);
                ir = int(irot - 0.5);
                jt = int(jrot + 0.5);
                jb = int(jrot - 0.5);
                slices[i][j] = timg[int(il)][int(jt)] * (ir - irot) * (jb - jrot);
                slices[i][j] += timg[int(ir)][int(jt)] * (irot - il) * (jb - jrot);
                slices[i][j] += timg[int(il)][int(jb)] * (ir - irot) * (jrot - jt);
                slices[i][j] += timg[int(ir)][int(jb)] * (irot - il) * (jrot - jt);
            }
            else
                slices[i][j] = 0.0E00;
        }
    }

    free_float2D(timg, N);

}

void find_z_intersection(float** img, float** imgbel, float** imgabv,
    float** zmap, int nslices, int nchan,
    float tiltx, int res, float angabv, float angbel)
{
    //Function grabs three tangent lines from three tilt images to estimate
    //the z-height of the middle tangent line surface intersection using
    //the above and below tangent line intersections. This function is
    //purely geometric and does not implement polynomial fitting.
    //This function is no longer used but was integral to the 1st implementation
    //of stomo, which successfully computed tangent intersections geometrically
    //using only single pairs of local tangent lines.  This function is left
    //in the source code for prosperity as it was implemented stomo_version_0.ccp
    //using (successfully) only pairs of images, without SVD fitting.
    //This function implements finite differences with geometry alone, whereas
    //find_z_intersection_II() uses differential geometry.

    int i, j, k, abvflag, belflag;
    float z, ci, cin1, cip1;

    for (i = 0; i < nslices; i++) {
        for (j = 0; j < nchan; j++) {
            zmap[i][j] = negflag;
        }
    }

    for (i = 0; i < nslices; i++) {
        for (j = 0; j < nchan; j++) {
            if (img[i][j] > 1e-5) {
                ci = j - tiltx;
                abvflag = 0;
                belflag = 0;
                for (k = j - res; k < j + res; k++) {
                    if (k >= 0 && k < nchan) {
                        if (imgbel[i][k] > 1e-5) {
                            cin1 = k - tiltx;
                            belflag = 1;
                            k = -1;
                            break;
                        }
                    }
                }
                for (k = j - res; k < j + res; k++) {
                    if (k >= 0 && k < nchan) {
                        if (imgabv[i][k] > 1e-5) {
                            cip1 = k - tiltx;
                            abvflag = 1;
                            k = -1;
                            break;
                        }
                    }
                }
                if (belflag == 1 && abvflag == 1) {
                    zmap[i][j] = (ci - cip1 / cos(angabv)) / tan(angabv) / 2.0E00;
                    zmap[i][j] -= (ci - cin1 / cos(angbel)) / tan(angbel) / 2.0E00;
                    zmap[i][j] += (ci - cin1 / cos(angbel)) / tan(angbel);
                }
                else
                    zmap[i][j] = negflag;
            }
        }
    }

}

float search_tangent_x(float** img, int res, int iedge,
    int jedge, int nchan, float tiltx)
{
    //Function analyses an image img of binary thinned edges.  An edge point iedge,
    //jedge from an image at a neighbouring tilt angle is supplied and img is
    //searched within the ith slice either side of iedge.  Parameter res sets the
    //breadth of the search along the j direction.

    int j, flag = 0, tj, jdel;
    float x;

    x = negflag;

    //parameter res sets the maximum search length
    for (j = jedge - res; j < jedge + res; j++) {
        if (j >= 0 && j < nchan) {
            jdel = int((res + (j - jedge)) / 2);
            if (j % 2 == 0) //search either side of jedge in a hopscotch pattern
                tj = jedge - jdel;
            else
                tj = jedge + jdel;
            //note: there is an inefficiency in that jdel is zero twice
            if (tj < nchan && tj >= 0) {
                if (img[iedge][tj] > 1e-5) {
                    x = float(tj) - tiltx;
                    break;
                }
            }
        }
    }

    return x;
}

void fit_poly(float* zs, float* angs, float* derivs,
    float* errors, float* sig, int npnts, unsigned int ma)
{
    //Implements polynomial fitting using SVD Numerical Recipes routines
    //The 1st derivative derivs[1] gives the variation of tangent abscissa zs[]
    //wrt to angle angs[].  The standard errors are used as user weighted
    //filters to crop poor fits that result from misidentified or broken edges.

    int i, j, m;
    unsigned int ndata, Mp, Np;
    float* x, * y, ** u, ** v, * a, * w, ** CVM, chisq, x2;

    ndata = npnts;

    Np = ndata;
    Mp = ndata;

    x = new float[ndata];
    y = new float[ndata];
    w = new float[ma];
    a = new float[ma];

    u = float2D(ndata, ma, "u");
    v = float2D(ma, ma, "v");
    CVM = float2D(ma, ma, "CVM");

    for (m = 0; m < npnts; m++) {
        x[m] = angs[m];
        y[m] = zs[m];
    }

    svdfit(x, y, sig, ndata, a, ma, u, v, w, Mp, Np, &chisq, fpoly);
    svdvar(v, ma, w, CVM);

    for (m = 0; m < ma; m++) {
        derivs[m] = a[m];
        errors[m] = sqrt(CVM[m][m] + 1e-20);
    }

    free_float2D(u, ndata);
    free_float2D(v, ma);
    free_float2D(CVM, ma);
    delete[] x;
    delete[] y;
    delete[] a;
    delete[] w;

}

void find_z_intersection_II(float*** imgs, float** zmap, float** emap, float* angs,
    float stdrr, float cut, int nslices, int nchan,
    int npnts, int minpnts, int midpnt, int nimgs,
    int norder, int res, float tiltx)
{
    //Function grabs a maximum of npts tangent-lines from npts tilt images to
    //estimate the z-height of the tangent-line to surface intersection using
    //a Talyor series fit (SVD least squares) to infer tangent abscissa angular
    //variations.

    int i, j, m, k, c, mpnts, tpnts, tj;
    float z, x, xt, avg, var;
    float* xs, * txs, * zs, * tangs, * dangs, * derivs, * errors, * sig;

    ofstream flog;
    flog.open("logfile.txt", ios::app);

    xs = new float[npnts];
    zs = new float[npnts];
    txs = new float[npnts];
    sig = new float[npnts];
    tangs = new float[npnts];
    dangs = new float[npnts];
    derivs = new float[norder];
    errors = new float[norder];
    mpnts = 0;

    for (k = 0; k < npnts; k++) {
        xs[k] = negflag;
        zs[k] = 0.0E00;
        txs[k] = 0.0E00;
        tangs[k] = 0.0E00;
        dangs[k] = 0.0E00;
        sig[k] = stdrr;
    }

    for (k = 0; k < norder; k++) {
        derivs[k] = 0.0E00;
        errors[k] = 0.0E00;
    }

    k = midpnt;

    for (i = 0; i < nslices; i++) {
        for (j = 0; j < nchan; j++) {
            zmap[i][j] = negflag;
        }
    }

    for (i = 0; i < nslices; i++) {
        for (j = 0; j < nchan; j++) {
            //calculate derivatives for all imgs, even if near angular limits
            if (imgs[k][i][j] > 1e-5) {
                x = j - tiltx;
                mpnts = 0;
                tj = j;
                //search near kth image in +ve direction first
                c = 0;
                for (m = k + 1; m < k + 1 + int(npnts / 2); m++) {
                    if (m < nimgs) {
                        //Note: could try a 1d cross-correlation to find these edges
                        //      or some other improved optical-flow tracking algorithm.
                        xs[c + int((npnts + 1) / 2)] = search_tangent_x(imgs[m], res, i, tj, nchan, tiltx);
                        tangs[c + int((npnts + 1) / 2)] = angs[m];
                        //set the next edge search channel to the previous found edge
                        if (xs[c + int((npnts + 1) / 2)] > (negflag + 1))
                            tj = int(xs[c + int((npnts + 1) / 2)] + tiltx);
                        mpnts++;
                    }
                    c++;
                }

                //for odd npnts, add in middle edge point (case m == k)
                //(replaces previous assignment for even npnts)
                xs[npnts / 2] = x;
                tangs[npnts / 2] = angs[k];
                if (!(npnts % 2 == 0))
                    mpnts++;

                tj = j;
                //now search near kth image in -ve direction
                c = 0;
                for (m = k - 1; m > k - 1 - int(npnts / 2); m--) {
                    if (m >= 0) {
                        xs[int(npnts / 2) - 1 - c] = search_tangent_x(imgs[m], res, i, tj, nchan, tiltx);
                        tangs[int(npnts / 2) - 1 - c] = angs[m];
                        //set the next edge search channel to the previous found edge
                        if (xs[int(npnts / 2) - 1 - c] > (negflag + 1))
                            tj = int(xs[int(npnts / 2) - 1 - c] + tiltx);
                        mpnts++;
                    }
                    c++;
                }

                //mpnts ordered from lowest m to highest m with k in the middle
                tpnts = mpnts;
                mpnts = 0;
                for (m = 0; m < tpnts; m++) {
                    if (xs[m] > (negflag + 1)) {
                        txs[mpnts] = xs[m];
                        dangs[mpnts] = tangs[m] - angs[k];
                        mpnts++;
                    }
                }
                if (mpnts > 1) { //there may also be no located points
                    for (m = 0; m < mpnts; m++)
                        zs[m] = txs[m] - x;

                    if (mpnts >= norder && mpnts > minpnts) {
                        fit_poly(zs, dangs, derivs, errors, sig, mpnts, (unsigned int)(norder));
                        //Division below is regularised to robustly handle small derivs
                        //Note: if errors is set by user to be large, then slow varying
                        //derivs are cut - this is can remove good data such as a
                        //stationary sphere, so some caution with cut is required.
                        if (fabs(errors[1] * derivs[1]) / (pow(derivs[1], 2) + 1e-10) < cut) {
                            zmap[i][j] = derivs[1];
                            //TP 2021: don't assign the edge map to the gradient error
                            //otherwise there will be outliers.
                            //emap[i][j] = errors[1];

                            emap[i][j] = imgs[k][i][j];
                        }
                        else
                            zmap[i][j] = negflag;
                    }
                }
                else
                    zmap[i][j] = negflag;
            }
        }
    }

    delete[] xs;
    delete[] zs;
    delete[] txs;
    delete[] sig;
    delete[] tangs;
    delete[] dangs;
    delete[] derivs;
}

void make_tip_img_II(float** img, float** typos, float** tyneg, float*** sypos,
    float*** syneg, int N, int M, float ang, float var,
    float width, float yoffset, float sabsorb, float tabsorb,
    float* zheights, float* ycenters, float* xcenters,
    float* radii, int nspheres)
{
    //Function simulates a TEM image of an atom probe tip using a quantitative
    //thickness map based on an elliptical cross section with parabolic tapering.
    //nsphere spheres with offset heights and positions along the tip axis are also
    //imbedded inside the tip.  Beer's law is used to find the total absorption
    //for the two different scattering densities specified by tabsorb and sabsorb.

    int i, j, m, ti, tj;
    float a, b, c, yp, yn, x, A, B, thick, sthick;
    float tang, toffj, tzcent, tjcent, zcent, jcent, av;
    float** sphimg;

    sphimg = float2D(N, M, "sphereimg");

    //tip thickness expression comes from equation of an ellipse with the x & y
    //coordinates rotated about the z-axis.  Calculating the new y gives a
    //quadratic, which is solved algebraically.  The profile of the tip (set by
    //ellipse parameters A and B) is just modeled by a parabola of varying
    //size, set by var.

    for (i = 0; i < N; i++) {
        for (j = 0; j < M; j++) {
            sphimg[i][j] = 0.0E00; //(sphimg[][] not used here - zeroed for later)
            x = j - M / 2;
            A = width * sqrt((i - N / 2 + yoffset + fabs(i - N / 2 + yoffset)) / 2);
            B = width * sqrt((i - N / 2 + yoffset + fabs(i - N / 2 + yoffset)) / 2) * var;
            a = pow(sin(ang) / A, 2) + pow(cos(ang) / B, 2);
            b = -2 * x * cos(ang) * sin(ang) * (B * B - A * A) / B / B / A / A;
            c = x * x * (pow(cos(ang) / A, 2) + pow(sin(ang) / B, 2)) - 1;

            if (4 * a * c < b * b && A>1e-5 && B > 1e-5) {
                yp = (-b + sqrt(b * b - 4 * a * c)) / 2.0E00 / a; //the two solutions define multiple
                yn = (-b - sqrt(b * b - 4 * a * c)) / 2.0E00 / a; //heights along beam axis
            }
            else
                yp = yn = 0;
            thick = fabs(yp - yn);
            img[i][j] = thick;
            typos[i][j] = yp;
            tyneg[i][j] = yn;
        }
    }
    //sphere thicknesses calculated here.  First, the spheres are subtracted
    //from the original tip thickness to create 'swiss cheese' of a single
    //tabsorb scattering density.  Then the spheres are added back in with their
    //different scattering density sabsorb.  sphimg[][] provides temporary storage
    //for the sphere thickness maps and accounts for sphere overlap in projection.

    for (m = 0; m < nspheres; m++) {
        tjcent = xcenters[m];
        tzcent = zheights[m];
        jcent = cos(ang) * tjcent + sin(ang) * tzcent; //Rotation about tip axis
        zcent = -sin(ang) * tjcent + cos(ang) * tzcent;

        for (i = 0; i < N; i++) {
            for (j = 0; j < M; j++) {
                ti = i + int(ycenters[m]);
                tj = j - int(jcent); //Note: zheights denote preset j-axis
                //displacements at zero angle, which swing with tilt
                sthick = radii[m] * radii[m] - (ti - N / 2) * (ti - N / 2) - (tj - M / 2) * (tj - M / 2);
                sthick = 2.0E00 * sqrt((sthick + fabs(sthick)) / 2 + 1e-20); //negatives removed
                if (sthick > 1e-5) {
                    sypos[m][i][j] += sthick / 2.0E00 + zcent;
                    syneg[m][i][j] += -sthick / 2.0E00 + zcent;
                }
                sphimg[i][j] += sthick;
            }
        }
    }

    av = 0.0E00;

    for (i = 0; i < N; i++) {
        for (j = 0; j < M; j++) {
            if (img[i][j] > 1e-5) {
                if (img[i][j] < sphimg[i][j] && sphimg[i][j]>1e-5)
                    img[i][j] -= sphimg[i][j]; //convert tip to Swiss cheese.
                else
                    if (sphimg[i][j] > (img[i][j] + 1e-5))
                        img[i][j] = 0.0E00; //in case sphere is fatter than the tip
                img[i][j] = exp(-tabsorb * img[i][j] - sabsorb * sphimg[i][j]); //Beer's law
            }
            else //allow spheres outside the tip aswell
                img[i][j] = exp(-sabsorb * sphimg[i][j]);
            av += img[i][j] / (N * M);
        }
    }

    for (i = 0; i < N; i++) {
        for (j = 0; j < M; j++) {
            img[i][j] /= av;
        }
    }

    free_float2D(sphimg, N);

}

void make_sphere_tilt_series_II(float*** imgs, int nimgs, int N, int M,
    int nspheres, float absorb, float* radii,
    float* zheights)
{
    //Function makes an image of off-centre spheres using
    //Beer's law and an analytical thickness map. zheight denotes the vertical
    //off-centre distances above the tilt axis.

    int i, j, k, m, jshift;
    float*** timgs, ang;

    timgs = float3D(nspheres, N, M, "imgs");

    for (m = 0; m < nspheres; m++) {
        for (i = 0; i < N; i++) {
            for (j = 0; j < M; j++) {
                timgs[m][i][j] = 0.0E00;
            }
        }
    }

    for (m = 0; m < nspheres; m++) {
        for (i = 0; i < N; i++) {
            for (j = 0; j < M; j++) {
                timgs[m][i][j] = radii[m] * radii[m] - (i - N / 2) * (i - N / 2) - (j - M / 2) * (j - M / 2);
                timgs[m][i][j] = sqrt((timgs[m][i][j] + fabs(timgs[m][i][j])) / 2.0E00 + 1e-20);
            }
        }

        for (k = 0; k < nimgs; k++) {
            ang = PI / nimgs * k;
            for (i = 0; i < N; i++) {
                for (j = 0; j < M; j++) {
                    jshift = j + int(zheights[m] * cos(ang) + 0.5);
                    if (jshift < M && jshift >= 0) {
                        imgs[k][i][j] += timgs[m][i][jshift];
                    }
                }
            }
        }
    }

    for (k = 0; k < nimgs; k++) {
        for (i = 0; i < N; i++) {
            for (j = 0; j < M; j++) {
                imgs[k][i][j] = exp(-1.0E00 * imgs[k][i][j] * absorb);
            }
        }
    }

    free_float3D(timgs, nspheres, N);

}

void make_tip_tilt_series_II(float*** imgs, int nimgs, int N, int M,
    float sabsorb, float tabsorb, float width,
    float var, float yoffset, float* zheights,
    float* ycenters, float* xcenters, float* radii,
    int nspheres)
{
    //Function makes an tilt series image stack of a tapered tip using Beer's law
    //and an empirical thickness map.

    int i, j, k;
    float** img, ** typ, ** tyn, *** syn, *** syp, angle;

    img = float2D(N, M, "img");
    typ = float2D(N, M, "typ");
    tyn = float2D(N, M, "tyn");
    syp = float3D(nspheres, N, M, "syp");
    syn = float3D(nspheres, N, M, "syn");

    ofstream fout;
    fout.open("logfile.txt", ios::app);
    cout << "\n";
    fout << "\n";
    for (k = 0; k < nimgs; k++) {
        cout << "\nImage " << k + 1 << " created";
        fout << "\nImage " << k + 1 << " created";
        angle = PI / nimgs * k + PI / 4; //Pi/4 start-angle accentuates the blade appearance
        make_tip_img_II(img, typ, tyn, syp, syn, N, M, angle, var, width, yoffset, sabsorb,
            tabsorb, zheights, ycenters, xcenters, radii, nspheres);
        for (i = 0; i < N; i++) {
            for (j = 0; j < M; j++) {
                imgs[k][i][j] = img[i][j];
            }
        }
    }

    fout.close();

    free_float2D(img, N);
    free_float2D(typ, N);
    free_float2D(tyn, N);
    free_float3D(syp, nspheres, N);
    free_float3D(syn, nspheres, N);

}

void make_tip_3d(float*** imgs, int nimgs, int N, int M, float width, float var,
    float yoffset, float start_angle, float sabsorb, float tabsorb,
    float* zheights, float* ycenters, float* xcenters,
    float* radii, int nspheres)
{
    //function makes a 3d image of a tip using an analytic model

    unsigned long count, n;
    int i, j, k, kth;
    float z, * xyz, ** img, *** syn, *** syp, ** typ, ** tyn, angle, thresh = 10;

    img = float2D(N, M, "img");
    typ = float2D(N, M, "typ");
    tyn = float2D(N, M, "tyn");
    syp = float3D(nspheres, N, M, "syp");
    syn = float3D(nspheres, N, M, "syn");

    for (i = 0; i < N; i++) {
        for (j = 0; j < M; j++) {
            img[i][j] = 0.0E00;
            typ[i][j] = 0.0E00;
            tyn[i][j] = 0.0E00;
            for (k = 0; k < nspheres; k++) {
                syp[k][i][j] = 0.0E00;
                syn[k][i][j] = 0.0E00;
            }
        }
    }

    angle = start_angle;
    make_tip_img_II(img, typ, tyn, syp, syn, N, M, angle, var, width, yoffset, sabsorb,
        tabsorb, zheights, ycenters, xcenters, radii, nspheres);

    for (k = 0; k < nimgs; k++) {
        for (i = 0; i < N; i++) {
            for (j = 0; j < M; j++) {
                imgs[k][i][j] = 0.0E00;
            }
        }
    }

    count = 0;
    for (i = 0; i < N; i++) {
        for (j = 0; j < M; j++) {
            if (fabs(tyn[i][j]) > 1e-5 && fabs(typ[i][j]) > 1e-5) {
                kth = int((tyn[i][j] / N * nimgs + nimgs / 2.0E00));
                imgs[kth][i][j] = kth; count += 3;
                kth = int((typ[i][j] / N * nimgs + nimgs / 2.0E00));
                imgs[kth][i][j] = kth; count += 3;
            }
            for (k = 0; k < nspheres; k++) {
                if (fabs(syn[k][i][j]) > 1e-5) {
                    kth = int((syn[k][i][j] / N * nimgs + nimgs / 2.0E00));
                    imgs[kth][i][j] = kth; count += 3;
                }
                if (fabs(syp[k][i][j]) > 1e-5) {
                    kth = int((syp[k][i][j] / N * nimgs + nimgs / 2.0E00));
                    imgs[kth][i][j] = kth; count += 3;
                }
            }
        }
    }

    xyz = new float[3 * count];
    n = 0;

    for (i = 0; i < N; i++) {
        for (j = 0; j < M; j++) {
            if (fabs(tyn[i][j]) > 1e-5 && fabs(typ[i][j]) > 1e-5) {
                z = tyn[i][j];
                xyz[n] = i; n++; xyz[n] = j; n++; xyz[n] = z; n++;
                z = typ[i][j];
                xyz[n] = i; n++; xyz[n] = j; n++; xyz[n] = z; n++;
            }
            for (k = 0; k < nspheres; k++) {
                if (fabs(syn[k][i][j]) > 1e-5) {
                    z = syn[k][i][j];
                    xyz[n] = i; n++; xyz[n] = j; n++; xyz[n] = z; n++;
                }
                if (fabs(syp[k][i][j]) > 1e-5) {
                    z = syp[k][i][j];
                    xyz[n] = i; n++; xyz[n] = j; n++; xyz[n] = z; n++;
                }
            }
        }
    }

    write_binary_1d_flarray(0, n, "xyz_tip_3d.dat", xyz);
    delete[] xyz;

    free_float2D(img, N);
    free_float2D(typ, N);
    free_float2D(tyn, N);
    free_float3D(syp, nspheres, N);
    free_float3D(syn, nspheres, N);

}

void shell_sort(unsigned long ndata, float* x, float* y, float* z)
{
    //Modified "Shell's method" for sorting arrays, x, y and z by the
    //entries in array z; adapted from Numerical Recipes in C, page 332.
    //quicksort() could be used here but adds unnecessary complexity.

    unsigned long i, j, inc;
    float xt, yt, zt;
    inc = 1;
    do {
        inc *= 3;
        inc++;
    } while (inc < ndata);
    do {
        inc /= 3;
        for (i = inc; i < ndata; i++) {
            zt = z[i];
            xt = x[i];
            yt = y[i];
            j = i;
            while (z[j - inc] > zt) {
                z[j] = z[j - inc];
                x[j] = x[j - inc];
                y[j] = y[j - inc];
                j -= inc;
                if (j < inc) break;
            }
            z[j] = zt;
            x[j] = xt;
            y[j] = yt;
        }
    } while (inc > 0);
}

void convert_point_cloud_to_tomogram(float* xyz, int nimgs,
    int N, int M, int ndata)
{
    //Function first sorts xyz[] by z to create a 3D voxelised tomogram
    //from point cloud xyz[] for writing to file 3d_SI.dat.  The extra complexity
    //of this algorithm as compared to stomo_version_1.cpp arises from the
    //desire to minimise the use of RAM.  It would be algorithmically simpler
    //to just declare the entire tomogram in a 3D float array and then map xyz[]
    //components to corresponding voxels.  Here an alternative method is used
    //whereby z components of xyz[] are first sorted and then consecutively dumped
    //into binned (using integer casting) tomogram slices, which are in inturn
    //appended to 3d_SI.dat.

    int i, j, k, zint, size, zmax, wflag, count;
    float x, y, z, * xa, * ya, * za, ** tomoslice;

    size = 1;
    zmax = 0;

    tomoslice = float2D(N, M, "tomogram slice");
    xa = new float[int(ndata / 3)];
    ya = new float[int(ndata / 3)];
    za = new float[int(ndata / 3)];

    for (i = 0; i < ndata; i += 3) {
        xa[int(i / 3)] = xyz[i + 0];
        ya[int(i / 3)] = xyz[i + 1];
        za[int(i / 3)] = xyz[i + 2];
    }
    shell_sort(int(ndata / 3), xa, ya, za);
    for (i = 0; i < ndata; i += 3) {
        xyz[i + 0] = xa[int(i / 3)];
        xyz[i + 1] = ya[int(i / 3)];
        xyz[i + 2] = za[int(i / 3)];
    }
    delete[] xa;
    delete[] ya;
    delete[] za;
    count = 0;

    for (j = 0; j < N; j++)
        for (k = 0; k < M; k++)
            tomoslice[j][k] = 0.0E00;

    for (i = 0; i < ndata; i += 3) {
        x = xyz[i + 0];
        y = xyz[i + 1];
        z = xyz[i + 2];
        zint = int(nimgs * z / N); //Bin the z-values into slices with scaling nimgs/N.
        if (int(x) >= 0 && int(x) < N && int(y) >= 0 && int(y) < M && zint == zmax) {
            tomoslice[int(x)][int(y)] += zint;
        }
        //Nested logic here ensures tomoslice is filled for all equivalent zint
        //before moving on to the next slice in the tomogram.  Though monotonic,
        //zint may increase by more than one unit between successives larger values.
        //In such cases the tomogram slice number is advanced and zero-values are
        //written to 3d_SI.dat for these slices. wflag dictates whether to write
        //to disc or not and when to advance to the next slice.  The variable
        //count tracks the slice number.
        if (zmax < zint) {
            zmax = zint;
            wflag = 1;
            count++;
        }
        else
            if (zmax > count) {
                wflag = 0;
                for (j = 0; j < N; j++)
                    for (k = 0; k < M; k++)
                        tomoslice[j][k] = 0.0E00;
                //fill in tomogram with gaps when zint jumps more than 1 slice.
                for (j = 0; j < zmax - count; j++) {
                    count++;
                    write_binary_3d_flarray(count - 1 + j, count + j, N, M, "3d_SI.dat",
                        "3d_SI.dat", &tomoslice);
                }
            }
            else
                wflag = 0;

        if (wflag && count <= nimgs) {
            write_binary_3d_flarray(count - 1, count, N, M, "3d_SI.dat", "3d_SI.dat", &tomoslice);
            //Flush the tomogram slice to prepare for the next accumulation.
            for (j = 0; j < N; j++)
                for (k = 0; k < M; k++)
                    tomoslice[j][k] = 0.0E00;
        }
    }
    //Pad 3d_SI.dat with zero-value slices if count does not reach nimgs.
    if (count < nimgs) {
        for (i = 0; i < nimgs - count; i++) {
            write_binary_3d_flarray(count + i, count + i + 1, N, M, "3d_SI.dat", "3d_SI.dat", &tomoslice);
        }
    }

}

void project_point_cloud_to_tilt_series(float* xyz, float* angles, float tiltx, int nimgs,
    int N, int M, int ndata)
{
    //Function rotates xyz coordinates back to recreate the tilt series
    //and then projects points onto each tilt image, for comparison
    //with edge maps, as a consistency check of the reconstruction.

    int a, i, j, k, jrot;
    float x, y, z, theta, ** tiltimage;

    tiltimage = float2D(N, M, "tiltimage");

    for (a = 0; a < nimgs; a++) {
        theta = angles[a];
        cout << "projecting angle number: " << a << " at " << angles[a] / PI * 180 << endl;

        for (j = 0; j < N; j++)
            for (k = 0; k < M; k++)
                tiltimage[j][k] = 0.0E00;

        for (i = 0; i < ndata; i += 3) {
            x = xyz[i + 0];
            y = xyz[i + 1];
            z = xyz[i + 2];

            //rotate the surface points to their respective tilt planes (x doesn't change)
            jrot = cos(theta) * (y - tiltx) + sin(theta) * (z - N / 2) + tiltx;

            if (int(x) >= 0 && int(x) < N && int(jrot) >= 0 && int(jrot) < M) {
                tiltimage[int(x)][int(jrot)] += 1;
            }
        }
        write_binary_3d_flarray(a, a + 1, N, M, "projections.dat", "projections.dat", &tiltimage);
    }

}


void compute_xyz(float** zmap, float** emap, float** gdir, int N, int M,
    int& numpnts, float theta, float tiltx)
{
    //Function calculates a partial xyz point cloud from an (x,y) map of relative
    //z-heights, as computed by the differential tangent intersection algorithms.

    int i, j, k, count, n, start;
    float* xyz, * nxyz, * zerrs, krot, jrot;
    float nx, ny, nz, tnx, tny, tnz, z, zerr;

    count = 0;

    for (i = 0; i < N; i++) {
        for (j = 0; j < M; j++) {
            if (zmap[i][j] > (negflag + 1)) {
                count++;
            }
        }
    }

    cout << "\n\n" << count << endl;
    xyz = new float[3 * count];
    nxyz = new float[3 * count];
    zerrs = new float[count];
    n = 0;

    start = numpnts;
    //reverse angles to rotate back (the derivatives only compute
    //angular variations relative to the tilt plane of interest)
    theta = -theta;

    for (i = 0; i < N; i++) {
        for (j = 0; j < M; j++) {
            if (zmap[i][j] > (negflag + 1)) {
                z = zmap[i][j];
                nx = sin(gdir[i][j] + PI / 2); //component of gdir from cross product
                tny = -cos(gdir[i][j] + PI / 2); //component of gdir from cross product
                tnz = 0.0E00;
                zerr = emap[i][j];
            }
            else {
                z = negflag;
                nx = 0;
                tny = 0;
                tnz = 0;
                zerr = 0;
            }

            //rotate the surface points to their respective tilt planes
            jrot = cos(theta) * (j - tiltx) + sin(theta) * z + tiltx;
            krot = -sin(theta) * (j - tiltx) + cos(theta) * z + N / 2;
            //rotate the normals
            ny = cos(theta) * tny + sin(theta) * tnz;
            nz = -sin(theta) * tny + cos(theta) * tnz;


            //ensure unit vectors are normalised (should be regardless)
            tnx = nx;
            tny = ny;
            tnz = nz;
            nx /= sqrt(tnx * tnx + tny * tny + tnz * tnz) + 1e-10;
            ny /= sqrt(tnx * tnx + tny * tny + tnz * tnz) + 1e-10;
            nz /= sqrt(tnx * tnx + tny * tny + tnz * tnz) + 1e-10;

            if (!(jrot < 0 || jrot >= N || krot < 0 || krot >= N || z < (negflag + 1))) {
                zerrs[n / 3] = zerr;
                xyz[n] = float(i);
                nxyz[n] = nx;
                n++;
                xyz[n] = jrot;
                nxyz[n] = ny;
                n++;
                xyz[n] = krot;
                nxyz[n] = nz;
                n++;
            }
        }
    }

    write_binary_1d_flarray(start, start + n, "xyz.dat", xyz);
    write_binary_1d_flarray(start, start + n, "norms.dat", nxyz);
    write_binary_1d_flarray(start, start + n / 3, "zerrors.dat", zerrs);
    numpnts += n;
    delete[] xyz;
    delete[] nxyz;
    delete[] zerrs;

}


void reconstruct_local(int& numpnts, int imgnum, int res, int N, int M,
    int tiltx, int minpnts, int ndata, int norder, int eflag, int vrflag,
    int nimgs, float tiltang, float clip, float fw, float T1,
    float T2, float stdrr, float cut, float colmin,
    float colmax, float rowmin, float rowmax, float* xyz,
    float* angles)
{
    //reconstruct_local() is the heart of stomo_version_2.cpp. All iterative
    //reconstruction algorithms in stomo_version_1.cpp have been reorgnised
    //here to process exp images and edge maps one at a time, using a buffered
    //3d array pertaining to only the number tilt images ndata required for the SVD
    //fit of the angular variations.  Version 1 performed the same computation
    //but effectively used a buffered 3d array with nimgs rather than ndata.  These
    //changes (along with similar iterative file writing) have been made to vastly
    //reduce memory usage and hence allow higher spatial and angular resolution
    //data sets to be analysed on standard PCs.

    int i, j, k, n, m, Nt, Mt;
    float kmax, kmin, height, theta, thetat;
    float z, irot, jrot, krot, ith, angle;
    float* radii, * zheights, * ycenters, * xcenters, * tangs;
    float** zmap, ** slices, ** mask, ** mask2, ** emap, ** gdir, ** timg, ** texp;
    float*** imgs, *** timgs;
    char filename[100];
    unsigned long count;

    numpnts = 0;

    norder += 1; //this variable is actually the number of polynomial basis funcs.
    texp = float2D(N, M, "texp");
    mask = float2D(N, M, "mask");
    tangs = new float[ndata];

    //mask used to crop edges detected at image boundaries
    for (i = 0; i < N; i++) {
        for (j = 0; j < M; j++) {
            mask[i][j] = 0.0E00;
        }
    }

    //rescale rectangular image size for rotation about tilt angle
    theta = tiltang / 180 * PI;
    Nt = N;
    Mt = M;
    N = int(Nt * cos(fabs(theta)) + Mt * sin(fabs(theta)));
    M = int(Mt * cos(fabs(theta)) + Nt * sin(fabs(theta)));

    imgs = float3D(ndata, N, M, "imgs");
    timgs = float3D(ndata, N, M, "timgs");
    timg = float2D(N, M, "timg");
    gdir = float2D(N, M, "gdir");

    ofstream flog;
    flog.open("logfile.txt", ios::app);

    if (fabs(tiltang) > 1e-5) {
        cout << "\n\nRotating every image by tiltang. ";
        cout << "\n\nNew image length is: " << N;
        cout << "\nNew image width is: " << M;
        flog << "\nNew image length is: " << N;
        flog << "\nNew image width is: " << M << endl;
        //.................................................
    }

    cout << "\n\nStarting calculation";

    zmap = float2D(N, M, "zmap");
    emap = float2D(N, M, "emap");
    mask2 = float2D(N, M, "mask2");
    emap = float2D(N, M, "emap");
    slices = float2D(N, M, "slices");

    for (i = 0; i < N; i++) {
        for (j = 0; j < M; j++) {
            mask2[i][j] = 0.0E00;
            emap[i][j] = 0.0E00;
            slices[i][j] = 0.0E00;
            zmap[i][j] = 0.0E00;
            for (k = 0; k < ndata; k++) {
                imgs[k][i][j] = 0.0E00;

                timgs[k][i][j] = 0.0E00;
            }
        }
    }

    //NOTE: the colmin rowmax etc create rectangular mask for clipping boundaries
    for (i = 0; i < Nt; i++) {
        for (j = 0; j < Mt; j++) {
            if (i > rowmax || i<rowmin || j>colmax || j < colmin)
                mask[i][j] = 0.0E00;
            else
                mask[i][j] = 1.0E00;
        }
    }
    //Rotate the mask with the tilted images and copy to bigger mask, mask2
    //in which mask is inscribed.
    rotate_imgs(mask, mask2, Nt, Mt, tiltang);

    //Alter mask2 to also crop edges outside of convolution windows
    for (i = 0; i < N; i++) {
        for (j = 0; j < M; j++) {
            if (i <= 3 * fw + 2 || i >= N - 3 * fw - 2 || j >= M - 3 * fw - 2 || j <= 3 * fw + 2)
                mask2[i][j] = 0.0E00;
        }
    }

    //Load up edge maps into buffer timgs to be searched forwards and backwards in
    //the tilt series for computation of local tangent intersertions.

    flog << "\n\nProcessing first image ";
    cout << "\n\nProcessing first image ";

    for (k = 0; k < ndata; k++) {
        read_binary_3d_flarray(k, k + 1, Nt, Mt, "exp_data.dat", &texp);
        rotate_imgs(texp, timg, Nt, Mt, tiltang);

        if (fabs(tiltang) > 1e-5 && k<int(ndata / 2))
            write_binary_3d_flarray(k, k + 1, N, M, "imgs_SI.dat", "imgs_SI.dat", &timg);

        //Detect edges with the Canny edge detector.  fix_wireframe_edges() clips
        //edges parallel to the tilt direction which contain no angular displacment
        //information. gdir[][] is passed to reconstruct_local() so that the normals
        //can be computed externally elsewhere in the code.
        edge_detect_Canny(timg, emap, gdir, fw, T1, T2, eflag, vrflag, N, M);

        //This function clips edges orthogonal to the tilt direction. The amount of
        //clipping is set by the user's "clip" value in stomo_parameters.txt.
        //When clip is set to 0.0E00, this function leaves the edge maps unaltered.

        //TP 2020: make this redundant, in favour of horizontal non-maximum supression
        //(next few lines below)
        /*fix_wireframe_edges(emap, res, clip, N, M);*/

        //Alternative method to "fix" edges by applying horizontal non-maximum suppression
        suppress_edges_horizontal(emap, N, M);

        for (i = 0; i < N; i++) {
            for (j = 0; j < M; j++) {
                imgs[k][i][j] = emap[i][j] * mask[i][j] * mask2[i][j];
            }
        }
    }

    for (m = 0; m < ndata; m++) {
        for (i = 0; i < N; i++) {
            for (j = 0; j < M; j++) {
                timgs[m][i][j] = imgs[m][i][j];
            }
        }
    }

    for (m = 0; m < ndata; m++) {
        for (i = 0; i < N; i++) {
            for (j = 0; j < M; j++) {
                if ((m - int(ndata / 2)) >= 0)
                    imgs[m][i][j] = timgs[m - int(ndata / 2)][i][j];
                else
                    imgs[m][i][j] = 0;
            }
        }
    }

    for (i = 0; i < ndata; i++)
        if ((i - int(ndata / 2)) >= 0 && (i - int(ndata / 2)) < nimgs)
            tangs[i] = angles[i -int(ndata / 2)];
        else
            tangs[i] = 0;

    //Find the z height of intersection tangents inscribed in the edge maps
    //by fitting a Taylor series to estimate local angular variations of edges.
    find_z_intersection_II(imgs, zmap, emap, tangs, stdrr, cut, N, M, ndata, minpnts,
        int(ndata / 2), ndata, norder, res, tiltx);

    //Rotate the z heights to their respective tilt planes and write the
    //resulting point cloud to disc.
    compute_xyz(zmap, emap, gdir, N, M, numpnts, angles[0], tiltx);

    //Write the edge maps to disc so that the user can inspect edge detection
    //quality and possibly change threshold parameters in a 2nd run.  In this
    //2nd version, the actual edge strengths are written, rather than a binary
    //true value; in differnce to version 1.  Inspection of the edge strengths
    //allows the user to accurately set desired threshold values.
    write_binary_3d_flarray(0, 1, N, M, "edge_maps.dat",
        "edge_maps.dat", &imgs[int(ndata / 2)]);

    //Repeat all of the above steps (edge detection, tangent intersection and
    //point cloud computation) hereafter in an iterative fashion, appending
    //binary output files as each tilt image is processed in turn, using a
    //moving batch of tilt images of size ndata to provide SVD fits about each
    //central angle indexed by k.
    for (k = 1; k < nimgs; k++) {
        cout << "\nProcessing image number " << k + 1;
        flog << "\nProcessing image number " << k + 1;
        if ((k + int(ndata / 2) - 1) < nimgs)
            read_binary_3d_flarray(k + int(ndata / 2) - 1, k + 1 + int(ndata / 2) - 1, Nt, Mt, "exp_data.dat", &texp);
        else
            for (i = 0; i < Nt; i++) {
                for (j = 0; j < Mt; j++) {
                    texp[i][j] = 0;
                }
            }

        rotate_imgs(texp, timg, Nt, Mt, tiltang);

        if ((k + int(ndata / 2 - 1)) < nimgs && fabs(tiltang) > 1e-5)
            write_binary_3d_flarray(k + int(ndata / 2) - 1, k + 1 + int(ndata / 2) - 1, N, M, "imgs_SI.dat", "imgs_SI.dat", &timg);

        edge_detect_Canny(timg, emap, gdir, fw, T1, T2, eflag, vrflag, N, M);

        //TP 2020: this feature has been made redundant, see below.
        /*fix_wireframe_edges(emap, res, clip, N, M);*/

        //Alternative method to "fix" edges by applying horizontal non-maximum suppression
        suppress_edges_horizontal(emap, N, M);

        for (i = 0; i < N; i++) {
            for (j = 0; j < M; j++) {
                timg[i][j] = emap[i][j] * mask2[i][j];
            }
        }
        for (m = 0; m < ndata; m++) {
            for (i = 0; i < N; i++) {
                for (j = 0; j < M; j++) {
                    timgs[m][i][j] = imgs[m][i][j];
                }
            }
        }
        for (m = 1; m < ndata; m++) {
            for (i = 0; i < N; i++) {
                for (j = 0; j < M; j++) {
                    imgs[m - 1][i][j] = timgs[m][i][j];
                }
            }
        }
        for (i = 0; i < N; i++) {
            for (j = 0; j < M; j++) {
                imgs[ndata - 1][i][j] = timg[i][j];
            }
        }

        if (k < nimgs) {
            write_binary_3d_flarray(k, k + 1, N, M, "edge_maps.dat",
                "edge_maps.dat", &imgs[int(ndata / 2)]);
        }
        for (i = 0; i < ndata; i++)
            if ((i - int(ndata / 2) + k) >= 0 && (i - int(ndata / 2) + k) < nimgs)
                tangs[i] = angles[i - int(ndata / 2) + k];
            else
                tangs[i] = 0;

        //Load in the next exp image, compute next edge map and append to imgs
        find_z_intersection_II(imgs, zmap, emap, tangs, stdrr, cut, N, M, ndata, minpnts,
            int(ndata / 2), ndata, norder, res, tiltx);
        compute_xyz(zmap, emap, gdir, N, M, numpnts, angles[k], tiltx);
    }

    cout << "\n\nCreated edge maps";
    flog << "\n\nCreated edge maps";
    cout << "\n\nThere are " << numpnts / 3 << " reconstructed points" << endl;
    flog << "\nThere are " << numpnts / 3 << " reconstructed points" << endl;
    flog << "\nReconstructed 3d geometry using the height map";
    flog.close();

    free_float2D(timg, N);
    free_float2D(gdir, N);
    free_float2D(zmap, N);
    free_float2D(emap, N);
    free_float3D(imgs, ndata, N);
    free_float3D(timgs, ndata, N);

    delete[] tangs;


}

//....end function definitions..........................................

int main()
{
    int i, j, k, n, nimgs, res, N, M, Nt, Mt, nstart, nend, npnts;
    int tiltx, ndata, minpnts, norder, err_count, eflag, vrflag, tomoflag, projflag, testflag, nspheres;
    float kmax, kmin, tiltang, width, clip, angbel, angabv;
    float fw, T1, T2, var, theta, thetat, height, z, irot, jrot, krot, ith;
    float tabsorb, sabsorb, angle, stdrr, cut, yoffset;
    float colmin, colmax, rowmin, rowmax;
    float* radii, * angles, * angs_partial, * zheights, * ycenters, * xcenters, * xyz;
    float** zmap, ** slices, ** mask, ** mask2, ** emap, * pcloud;
    float*** zmaps, *** timgs;
    char filename[100];
    unsigned long count;

    ofstream flog;
    flog.open("logfile.txt", ios::app);
    flog << "\n\nRunning surface tomography code" << endl;
    flog.close();

    cout << "\n\nRunning surface tomography code";
    ifstream fin("stomo_parameters.txt");

    fin >> nimgs; fin.ignore(300, '\n');
    fin >> N; fin.ignore(300, '\n');
    fin >> M; fin.ignore(300, '\n');
    fin >> fw; fin.ignore(300, '\n');
    fin >> T1; fin.ignore(300, '\n');
    fin >> T2; fin.ignore(300, '\n');
    fin >> ndata; fin.ignore(300, '\n');
    fin >> minpnts; fin.ignore(300, '\n');
    fin >> norder; fin.ignore(300, '\n');
    fin >> colmin; fin.ignore(300, '\n');
    fin >> rowmin; fin.ignore(300, '\n');
    fin >> colmax; fin.ignore(300, '\n');
    fin >> rowmax; fin.ignore(300, '\n');
    fin >> res; fin.ignore(300, '\n');
    fin >> tiltx; fin.ignore(300, '\n');
    fin >> tiltang; fin.ignore(300, '\n');
    //TP 2020: clip now redundant.
    /*fin >> clip; fin.ignore(300, '\n');*/
    fin >> stdrr; fin.ignore(300, '\n');
    fin >> cut; fin.ignore(300, '\n');
    fin >> eflag; fin.ignore(300, '\n');
    fin >> testflag; fin.ignore(300, '\n');
    fin >> vrflag; fin.ignore(300, '\n');
    fin >> tomoflag; fin.ignore(300, '\n');
    fin >> projflag; fin.ignore(300, '\n');
    //vrflag = 0; //this selects only RIDGES.  1 selects only valley.

    angles = new float[nimgs];

    for (i = 0; i < nimgs; i++) {
        fin >> angles[i];
        fin.ignore(300, '\n');
        angles[i] = angles[i] / 180 * PI;
    }
    fin.close();

    flog.open("logfile.txt", ios::trunc);
    flog << "\n\nSurface tomography running" << endl;

    xyz = new float[100000000];

    clip = 0.0; //TP 2020: this parameter is now redundant, as fix_wireframe_edges is no longer called

    if (testflag == 0) {
        reconstruct_local(npnts, nimgs, res, N, M, tiltx, minpnts, ndata, norder, eflag, vrflag,
            nimgs, tiltang, clip, fw, T1, T2, stdrr, cut, colmin, colmax, rowmin,
            rowmax, xyz, angles);
        pcloud = new float[npnts];
        ifstream fin;
        fin.open("xyz.dat", ios::in | ios::binary);

        //Read the xyz file into the xyz array
        fin.read((char*)pcloud, npnts * sizeof(float));
        fin.close();

        //Tomogram file 3d_imgs.dat must be computed slice by slice to reduce
        //memory usage.  It was not feasible to carry out these calculations
        //inside reconstruct_local().  For ease of coding, the point cloud is
        //simply read from disk here and the tomogram is constructed by
        //projecting and summing shaded slices of the cloud.  The size of the
        //tomogram is set to that of exp_imgs (could be any size).

        if (tomoflag == 1) {
            //bin xyz points onto a (large) structured grid, acting as a binned tomogram.
            convert_point_cloud_to_tomogram(pcloud, nimgs, N, M, npnts);
        }

        if (projflag == 1) {
            //project all points onto each tilt plane, for comparison between reconstruction
            //and original edgemaps or exp_data
            project_point_cloud_to_tilt_series(pcloud, angles, tiltx, nimgs, N, M, npnts);
        }

        cout << "\n\nReconstructed 3d geometry using the height map";
        cout << "\nTomogram written to 3d_SI.dat (same file size as exp_data.dat)";
        cout << "\n\nPress any key and hit enter to finish ";
        cin >> i;

        flog.close();
        exit(0);
    }

    if (testflag == 1) {
        //Unlike version 1, the test here must first be created and copied to
        //exp_data.dat before a reconstruction is computed.  The reconstruction
        //then only uses a small amount of RAM.
        flog << "\n\nCalculating test image stack" << endl;

        //Note: memory usage here identical to that in stomo_version_1.cpp
        timgs = float3D(nimgs, N, M, "timgs");

        for (i = 0; i < N; i++) {
            for (j = 0; j < M; j++) {
                for (k = 0; k < nimgs; k++) {
                    timgs[k][i][j] = 0.0E00;
                }
            }
        }

        //hard-code the parameters for the simulated image stack
        nspheres = 5;
        for (i = 0; i < nimgs; i++) {
            angles[i] = i * 1.0E00 / nimgs * PI;
        }
        ycenters = new float[nspheres]; //displacements along tip axis
        xcenters = new float[nspheres]; //othogonal to tip axis
        zheights = new float[nspheres]; //othogonal to tip axis
        radii = new float[nspheres];

        radii[0] = 25; radii[1] = 42; radii[2] = 30; radii[3] = 50; radii[4] = 20;
        ycenters[0] = -100; ycenters[1] = -180; ycenters[2] = 10; ycenters[3] = -40;
        ycenters[4] = 120; xcenters[0] = -40; xcenters[1] = 0; xcenters[2] = 30;
        xcenters[3] = 20; xcenters[4] = 25;
        zheights[0] = zheights[1] = zheights[2] = zheights[3] = zheights[4] = 0;
        sabsorb = 0.002; tabsorb = 0.001; //sphere and tip absorption params
        yoffset = 230; width = 5; var = 2;

        //scale invariance for arbitrary sized images
        yoffset /= 512.0; yoffset *= N;
        width /= sqrt(512.0); width *= sqrt(N);

        for (i = 0; i < nspheres; i++) {
            radii[i] /= 512.0; radii[i] *= N;
            xcenters[i] /= 512.0; xcenters[i] *= N;
            ycenters[i] /= 512.0; ycenters[i] *= N;
        }

        make_tip_tilt_series_II(timgs, nimgs, N, M, sabsorb, tabsorb, width, var, yoffset,
            zheights, ycenters, xcenters, radii, nspheres);
        angle = PI / 4; // Note: the angles in make_tip_images must match this
                      //starting angle otherwise tip_3d.dat won't match 3d_SI.dat
        write_binary_3d_flarray(0, nimgs, N, M, "tip_imgs.dat",
            " image stack of slices of exact tip model", timgs);
        make_tip_3d(timgs, nimgs, N, M, width, var, yoffset, angle, tabsorb, sabsorb,
            zheights, ycenters, xcenters, radii, nspheres);
        write_binary_3d_flarray(0, nimgs, N, M, "tip_3d.dat",
            " image stack of slices of exact tip model", timgs);
        flog.close();
    }

    return 0;
}
