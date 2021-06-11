#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>
#include <assert.h>

#include "memory.h"

#define PI 3.141592654

using namespace std;

/*---------------------------- float1D() -------------------------------*/
/*
    1D array allocator for type float
    make space for m[0...(n-1)]
    printf error message and exit if not successful

    this save checking for a NULL return etc every time

*/
float* float1D(int n, const char* message)
{
    float* m;

    int i;

    m = (float*)malloc(n * sizeof(float));
    if (m == NULL) {
        printf("float1D() cannot allocate memory size=%d: %s\n",
            n, message);
        std::cin >> i;
        exit(0);
    }
    return(m);

} /* end float1D() */

float** float2D(int nx, int ny, const char* message)
{
    float** m;
    int i;

    m = (float**)malloc(nx * sizeof(float*));
    if (m == NULL) {
        printf("float2D cannot allocate pointers, size=%d: %s\n",
            nx, message);
        std::cin >> i;
        exit(0);
    }

    for (i = 0; i < nx; i++) {
        m[i] = (float*)malloc(ny * sizeof(float));
        if (m[i] == NULL) {
            printf("float2D cannot allocate arrays, size=%d: %s\n",
                ny, message);
            std::cin >> i;
            exit(0);
        }
    }

    return m;

}  /* end float2D() */

float*** float3D(int na, int nx, int ny, const char* message);
void free_float3D(float*** a, int na, int nx);

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


void write_binary_3d_flarray(int nimgs, int N, int M, const char* file, const char* fn, float*** array);

void write_binary_3d_flarray(int nimgs, int N, int M,
    const char* file, const char* fn, float*** array)
{
    // Write to a *.dat floating point file, a 3D image stack of type float
    int i, j, k;
    unsigned long count, size;
    float* buffer;

    ////cout<<"\n\nBeginning file write for "<<fn<<" ..."<<endl;

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
    fout.open(file, ios::app | ios::out | ios::binary);
    //...........................................................

    //write the buffer array to the file
    fout.write((char*)buffer, size * sizeof(float));
    fout.close();
 //...........................................
    delete[] buffer;
}



void free_float2D(float** a, int nx)
{
    int i;

    for (i = nx - 1; i >= 0; i--) {
        free(a[i]);
    }
    free(a);

}  /* end free_float2D() */

void write_binary_2d_flarray(int N, int M, char* file, float** array)
{
    // Write to a Digital Micrograph *.dat floating point file
    // a 2D array of type float

    int i, j, count, size;
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
    std::ofstream fout;

    fout.open(file, std::ios::out | std::ios::binary);
    //...........................................................

    //Read the entire file into the buffer array
    fout.write((char*)buffer, size * sizeof(float));
    fout.close();
    //...........................................

}

void matrix_multiply(float** A, float** B, float** C)
{

    C[0][0] = A[0][0] * B[0][0] + A[0][1] * B[1][0] + A[0][2] * B[2][0];
    C[0][1] = A[0][0] * B[0][1] + A[0][1] * B[1][1] + A[0][2] * B[2][1];
    C[0][2] = A[0][0] * B[0][2] + A[0][1] * B[1][2] + A[0][2] * B[2][2];

    C[1][0] = A[1][0] * B[0][0] + A[1][1] * B[1][0] + A[1][2] * B[2][0];
    C[1][1] = A[1][0] * B[0][1] + A[1][1] * B[1][1] + A[1][2] * B[2][1];
    C[1][2] = A[1][0] * B[0][2] + A[1][1] * B[1][2] + A[1][2] * B[2][2];

    C[2][0] = A[2][0] * B[0][0] + A[2][1] * B[1][0] + A[2][2] * B[2][0];
    C[2][1] = A[2][0] * B[0][1] + A[2][1] * B[1][1] + A[2][2] * B[2][1];
    C[2][2] = A[2][0] * B[0][2] + A[2][1] * B[1][2] + A[2][2] * B[2][2];

}


void generate_box_points(float* x, float* y, float* z,
    float xcell, float ycell, float zcell, float masksize, int ngrid, int &count)
{
    //Function generates a set of points for the surfaces of the box
    //with the box symmetry axes alligned along the axes of a Cartesian
    //coordinate system.
    //cellsize denotes the xyz lengths of the box

    int i, j;

    for (i = 0; i < ngrid; i++) {
        for (j = 0; j < ngrid; j++) {
            x[count] = float(i) / ngrid * xcell;
            y[count] = float(j) / ngrid * ycell;
            z[count] = 0;
            count++;
            x[count] = float(i) / ngrid * xcell;
            y[count] = 0;
            z[count] = float(j) / ngrid * zcell;
            count++;
            x[count] = 0;
            y[count] = float(i) / ngrid * ycell;
            z[count] = float(j) / ngrid * zcell;
            count++;
            x[count] = float(i) / ngrid * xcell;
            y[count] = float(j) / ngrid * ycell;
            z[count] = zcell;
            count++;
            x[count] = float(i) / ngrid * xcell;
            y[count] = ycell;
            z[count] = float(j) / ngrid * zcell;
            count++;
            x[count] = xcell;
            y[count] = float(i) / ngrid * ycell;
            z[count] = float(j) / ngrid * zcell;
            count++;
        }
    }

}

void generate_box_shadow_mask(float* x, float* y, float xcell, float ycell, float** mask,
    float maxsize, int ngrid, int smpoints)
{
    //Function generates a set of points for the box
    //in any orientation, over which the thickness map
    //is to be calculated.  i.e. the projected shadow
    //of the box in the xy plane is calculated and then
    //treated as a mask, with all points in the shadow
    //equal to 1.0 and those outside equal to 0.
    //It is assumed that the box centre is (xcell/2,ycell/2,zcell/2)

    //maxsize is the length of a larger square within which the shadow
    //should be contained.  i.e. maxsize defines the 0 boundaries of the mask.

    int i, j, binx, biny, points;
    int ncubepoints = 6 * ngrid * ngrid;
    float xt, yt;

    for (i = 0; i < smpoints; i++) {
        for (j = 0; j < smpoints; j++) {
            mask[i][j] = 0.0E00;
        }
    }

    for (i = 0; i < ncubepoints; i++) {
        xt = x[i] + maxsize / 2 - xcell / 2;
        yt = y[i] + maxsize / 2 - ycell / 2;
        binx = int(xt / maxsize * smpoints);
        biny = int(yt / maxsize * smpoints);
        mask[binx][biny] = 1.0E00;
    }

}

void determinant_three_by_three(float** a, float* D)
{
    //Function evaluates the determinant of a 3 x 3 matrix
    //and returns the result by reference to *D
    *D = a[0][0] * (a[1][1] * a[2][2] - a[1][2] * a[2][1]);
    *D -= a[0][1] * (a[1][0] * a[2][2] - a[1][2] * a[2][0]);
    *D += a[0][2] * (a[1][0] * a[2][1] - a[1][1] * a[2][0]);

}

void minor_ij(float** a, float** Mij, int i, int j, int size)
{
    //Function contracts the matrix of size (size x size)
    //into a minor defined by i and j to give Mij

    int m, n, mcount = 0, ncount = 0;

    for (m = 0; m < size; m++) {
        ncount = 0;
        for (n = 0; n < size; n++) {
            if (m != i && n != j) {
                Mij[mcount][ncount] = a[m][n];
                ncount++;
            }
        }
        if (m != i)
            mcount++;
    }

}

void determinant_four_by_four(float** a, float* D)
{
    //Function evaluates the determinant of a 4 x 4 matrix
    //and returns the result by reference to *D

    int i, j, size = 4;
    float d;

    float** Mij;

    Mij = float2D(size - 1, size - 1, "Mij");

    for (i = 0; i < size - 1; i++) {
        for (j = 0; j < size - 1; j++) {
            Mij[i][j] = 0.0E00;
        }
    }

    minor_ij(a, Mij, 0, 0, size);
    determinant_three_by_three(Mij, &d);
    *D = a[0][0] * (d);
    minor_ij(a, Mij, 0, 1, size);
    determinant_three_by_three(Mij, &d);
    *D -= a[0][1] * (d);
    minor_ij(a, Mij, 0, 2, size);
    determinant_three_by_three(Mij, &d);
    *D += a[0][2] * (d);
    minor_ij(a, Mij, 0, 3, size);
    determinant_three_by_three(Mij, &d);
    *D -= a[0][3] * (d);

    free_float2D(Mij, size - 1);

}

void line_plane_intersect(float x1, float y1, float z1,
    float x2, float y2, float z2,
    float x3, float y3, float z3,
    float x4, float y4, float z4,
    float x5, float y5, float z5,
    float* x, float* y, float* z)
{

    //Function calculates the intersection of a line defined by points
    //x4,y4,z4 and x5,y5,z5 with a plane defined by points x1,y1,z1
    //x2,y2,z2 and x3,y3,z3.  The intesection coordinates are returned in x, y and z.

    int size = 4;
    static int sign = 1;
    float Ddenom, Dnumer, t, tol = 1e-20;
    float** Adenom, ** Anumer;

    Adenom = float2D(size, size, "Adenom");
    Anumer = float2D(size, size, "Anumer");

    Anumer[0][0] = 1; Anumer[0][1] = 1; Anumer[0][2] = 1; Anumer[0][3] = 1;
    Anumer[1][0] = x1; Anumer[1][1] = x2; Anumer[1][2] = x3; Anumer[1][3] = x4;
    Anumer[2][0] = y1; Anumer[2][1] = y2; Anumer[2][2] = y3; Anumer[2][3] = y4;
    Anumer[3][0] = z1; Anumer[3][1] = z2; Anumer[3][2] = z3; Anumer[3][3] = z4;

    Adenom[0][0] = 1; Adenom[0][1] = 1; Adenom[0][2] = 1; Adenom[0][3] = 0;
    Adenom[1][0] = x1; Adenom[1][1] = x2; Adenom[1][2] = x3; Adenom[1][3] = (x5 - x4);
    Adenom[2][0] = y1; Adenom[2][1] = y2; Adenom[2][2] = y3; Adenom[2][3] = (y5 - y4);
    Adenom[3][0] = z1; Adenom[3][1] = z2; Adenom[3][2] = z3; Adenom[3][3] = (z5 - z4);

    determinant_four_by_four(Anumer, &Dnumer);
    determinant_four_by_four(Adenom, &Ddenom);

    if (fabs(Ddenom) > tol) {
        t = -Dnumer / Ddenom;
        *x = x4 + t * (x5 - x4);
        *y = y4 + t * (y5 - y4);
        *z = z4 + t * (z5 - z4);
    }
    else {
        *x = 0;
        *y = 0;
        *z = 1e10 * sign;
        sign *= -1;
    }

    free_float2D(Adenom, size);
    free_float2D(Anumer, size);

}

void define_box_planes(float** p1, float** p2, float** p3,
    float** p4, float** p5, float** p6,
    float lxbox, float lybox, float lzbox, float lcell)
{
    // lcell is a box within which the box resides
    // The six box planes are defined by three points on the vertices
    // of the box (x1,y1,z1; x2,y2,z2; and x3,y3,z3)

    int i, j;

    p1[0][0] = lxbox / 2; p1[0][1] = lxbox / 2; p1[0][2] = -lxbox / 2;
    p1[1][0] = lybox / 2; p1[1][1] = lybox / 2; p1[1][2] = lybox / 2;
    p1[2][0] = -lzbox / 2; p1[2][1] = lzbox / 2; p1[2][2] = -lzbox / 2;

    p2[0][0] = lxbox / 2; p2[0][1] = lxbox / 2; p2[0][2] = lxbox / 2;
    p2[1][0] = lybox / 2; p2[1][1] = lybox / 2; p2[1][2] = -lybox / 2;
    p2[2][0] = -lzbox / 2; p2[2][1] = lzbox / 2; p2[2][2] = -lzbox / 2;

    p3[0][0] = -lxbox / 2; p3[0][1] = lxbox / 2; p3[0][2] = -lxbox / 2;
    p3[1][0] = lybox / 2; p3[1][1] = -lybox / 2; p3[1][2] = -lybox / 2;
    p3[2][0] = -lzbox / 2; p3[2][1] = -lzbox / 2; p3[2][2] = -lzbox / 2;

    p4[0][0] = lxbox / 2; p4[0][1] = -lxbox / 2; p4[0][2] = lxbox / 2;
    p4[1][0] = -lybox / 2; p4[1][1] = lybox / 2; p4[1][2] = lybox / 2;
    p4[2][0] = lzbox / 2; p4[2][1] = lzbox / 2; p4[2][2] = lzbox / 2;

    p5[0][0] = -lxbox / 2; p5[0][1] = -lxbox / 2; p5[0][2] = -lxbox / 2;
    p5[1][0] = lybox / 2; p5[1][1] = lybox / 2; p5[1][2] = -lybox / 2;
    p5[2][0] = -lzbox / 2; p5[2][1] = lzbox / 2; p5[2][2] = -lzbox / 2;

    p6[0][0] = -lxbox / 2; p6[0][1] = lxbox / 2; p6[0][2] = -lxbox / 2;
    p6[1][0] = -lybox / 2; p6[1][1] = -lybox / 2; p6[1][2] = -lybox / 2;
    p6[2][0] = lzbox / 2; p6[2][1] = -lzbox / 2; p6[2][2] = -lzbox / 2;

    for (i = 0; i < 3; i++) {
        for (j = 0; j < 3; j++) {
            if (i < 2) {
                p1[i][j] += lcell / 2;
                p2[i][j] += lcell / 2;
                p3[i][j] += lcell / 2;
                p4[i][j] += lcell / 2;
                p5[i][j] += lcell / 2;
                p6[i][j] += lcell / 2;
            }
        }
    }

}

void specify_Rxyz(float** Rx, float** Ry, float** Rz,
    float alpha, float beta, float gamma)
{

    int i, j;

    for (i = 0; i < 3; i++) {
        for (j = 0; j < 3; j++) {
            Rx[i][j] = 0.0E00;
            Ry[i][j] = 0.0E00;
            Rz[i][j] = 0.0E00;
        }
    }

    Rx[0][0] = 1;
    Rx[1][1] = cos(alpha);
    Rx[2][2] = cos(alpha);
    Rx[1][2] = sin(alpha);
    Rx[2][1] = -sin(alpha);

    Ry[0][0] = cos(beta);
    Ry[2][2] = cos(beta);
    Ry[2][0] = sin(beta);
    Ry[0][2] = -sin(beta);
    Ry[1][1] = 1;

    Rz[0][0] = cos(gamma);
    Rz[1][1] = cos(gamma);
    Rz[0][1] = sin(gamma);
    Rz[1][0] = -sin(gamma);
    Rz[2][2] = 1;
}

void rotate_xyz(float* x, float* y, float* z, float* nx, float* ny,
    float* nz, float** R,
    float xcell, float ycell, float zcell, int npoints)
{
    //Function has been modified to rotate the coordinate system using
    //just one matrix R[] containing the direction cosines

    float xt, yt, zt;
    int m;

    for (m = 0; m < npoints; m++) {
        xt = x[m] - xcell / 2.0E00;
        yt = y[m] - ycell / 2.0E00;
        zt = z[m] - zcell / 2.0E00;

        nx[m] = xt * R[0][0] + yt * R[0][1] + zt * R[0][2];
        ny[m] = xt * R[1][0] + yt * R[1][1] + zt * R[1][2];
        nz[m] = xt * R[2][0] + yt * R[2][1] + zt * R[2][2];

        nx[m] += xcell / 2.0E00;
        ny[m] += ycell / 2.0E00;
        nz[m] += zcell / 2.0E00;
    }
}


void rotate_box_planes(float** R, float** p1, float** p2, float** p3,
    float** p4, float** p5, float** p6,
    float lcell)
{
    //Function rotates all of the points contained in the plane
    //arrays p1 ==> p6, where the 1st indice denotes the x coordinate
    //and the 2nd indice denotes the plane position number
    //(there are three points to each plane)

    int j, npoints = 1;
    float* nx, * ny, * nz;
    float a11, a12, a13, a21, a22, a23, a31, a32, a33;

    nx = float1D(1, "nx");
    ny = float1D(1, "ny");
    nz = float1D(1, "nz");

    for (j = 0; j < 3; j++) {
        rotate_xyz(&p1[0][j], &p1[1][j], &p1[2][j], nx, ny, nz, R, lcell, lcell, 0/*lcell*/, npoints);
        p1[0][j] = nx[0]; p1[1][j] = ny[0]; p1[2][j] = nz[0];
        rotate_xyz(&p2[0][j], &p2[1][j], &p2[2][j], nx, ny, nz, R, lcell, lcell, 0/*lcell*/, npoints);
        p2[0][j] = nx[0]; p2[1][j] = ny[0]; p2[2][j] = nz[0];
        rotate_xyz(&p3[0][j], &p3[1][j], &p3[2][j], nx, ny, nz, R, lcell, lcell, 0/*lcell*/, npoints);
        p3[0][j] = nx[0]; p3[1][j] = ny[0]; p3[2][j] = nz[0];
        rotate_xyz(&p4[0][j], &p4[1][j], &p4[2][j], nx, ny, nz, R, lcell, lcell, 0/*lcell*/, npoints);
        p4[0][j] = nx[0]; p4[1][j] = ny[0]; p4[2][j] = nz[0];
        rotate_xyz(&p5[0][j], &p5[1][j], &p5[2][j], nx, ny, nz, R, lcell, lcell, 0/*lcell*/, npoints);
        p5[0][j] = nx[0]; p5[1][j] = ny[0]; p5[2][j] = nz[0];
        rotate_xyz(&p6[0][j], &p6[1][j], &p6[2][j], nx, ny, nz, R, lcell, lcell, 0/*lcell*/, npoints);
        p6[0][j] = nx[0]; p6[1][j] = ny[0]; p6[2][j] = nz[0];
    }

    free(nx);
    free(ny);
    free(nz);

}

void find_intersections(float** p1, float** p2, float** p3,
    float** p4, float** p5, float** p6,
    float lx1, float ly1, float lz1,
    float lx2, float ly2, float lz2,
    float* r1int, float* r2int, float* r3int,
    float* r4int, float* r5int, float* r6int)
{
    //Function calculates six different intersection vectors with the line
    //defined by the coordinates lx1==>lz1 and lx2==>lz2
    //for all of the planes of the box and stores them as rxint coordinates.

    float x1, x2, x3, y1, y2, y3, z1, z2, z3;
    float x, y, z;

    x1 = p1[0][0]; y1 = p1[1][0]; z1 = p1[2][0];
    x2 = p1[0][1]; y2 = p1[1][1]; z2 = p1[2][1];
    x3 = p1[0][2]; y3 = p1[1][2]; z3 = p1[2][2];
    line_plane_intersect(x1, y1, z1, x2, y2, z2, x3, y3, z3, lx1, ly1, lz1, lx2, ly2, lz2, &x, &y, &z);
    r1int[0] = x; r1int[1] = y; r1int[2] = z;
    x1 = p2[0][0]; y1 = p2[1][0]; z1 = p2[2][0];
    x2 = p2[0][1]; y2 = p2[1][1]; z2 = p2[2][1];
    x3 = p2[0][2]; y3 = p2[1][2]; z3 = p2[2][2];
    line_plane_intersect(x1, y1, z1, x2, y2, z2, x3, y3, z3, lx1, ly1, lz1, lx2, ly2, lz2, &x, &y, &z);
    r2int[0] = x; r2int[1] = y; r2int[2] = z;
    x1 = p3[0][0]; y1 = p3[1][0]; z1 = p3[2][0];
    x2 = p3[0][1]; y2 = p3[1][1]; z2 = p3[2][1];
    x3 = p3[0][2]; y3 = p3[1][2]; z3 = p3[2][2];
    line_plane_intersect(x1, y1, z1, x2, y2, z2, x3, y3, z3, lx1, ly1, lz1, lx2, ly2, lz2, &x, &y, &z);
    r3int[0] = x; r3int[1] = y; r3int[2] = z;
    x1 = p4[0][0]; y1 = p4[1][0]; z1 = p4[2][0];
    x2 = p4[0][1]; y2 = p4[1][1]; z2 = p4[2][1];
    x3 = p4[0][2]; y3 = p4[1][2]; z3 = p4[2][2];
    line_plane_intersect(x1, y1, z1, x2, y2, z2, x3, y3, z3, lx1, ly1, lz1, lx2, ly2, lz2, &x, &y, &z);
    r4int[0] = x; r4int[1] = y; r4int[2] = z;
    x1 = p5[0][0]; y1 = p5[1][0]; z1 = p5[2][0];
    x2 = p5[0][1]; y2 = p5[1][1]; z2 = p5[2][1];
    x3 = p5[0][2]; y3 = p5[1][2]; z3 = p5[2][2];
    line_plane_intersect(x1, y1, z1, x2, y2, z2, x3, y3, z3, lx1, ly1, lz1, lx2, ly2, lz2, &x, &y, &z);
    r5int[0] = x; r5int[1] = y; r5int[2] = z;
    x1 = p6[0][0]; y1 = p6[1][0]; z1 = p6[2][0];
    x2 = p6[0][1]; y2 = p6[1][1]; z2 = p6[2][1];
    x3 = p6[0][2]; y3 = p6[1][2]; z3 = p6[2][2];
    line_plane_intersect(x1, y1, z1, x2, y2, z2, x3, y3, z3, lx1, ly1, lz1, lx2, ly2, lz2, &x, &y, &z);
    r6int[0] = x; r6int[1] = y; r6int[2] = z;

}
void find_thickness(float* r1int, float* r2int, float* r3int,
    float* r4int, float* r5int, float* r6int, float* t, float *zmom)
{
    //Function sorts the intersection coordinates r in terms of z
    //height and then finds the difference b/n the two innermost
    //z values, which is the thickness.
    //The x and y values are not actually used.

    int i, imin, imax;
    float* zsort, max, temp;
    zsort = new float[6];

    zsort[0] = r1int[2]; zsort[1] = r2int[2]; zsort[2] = r3int[2];
    zsort[3] = r4int[2]; zsort[4] = r5int[2]; zsort[5] = r6int[2];

    //Perform very dumb sort of zsort array
    imin = 0;
    do {
        max = -1e20;
        for (i = imin; i < 6; i++) {
            if (zsort[i] > max) {
                max = zsort[i];
                imax = i;
            }
        }
        temp = zsort[imin];
        zsort[imin] = zsort[imax];
        zsort[imax] = temp;

        imin++;

    } while (imin < 6);
    //.....................
    *t = fabs(zsort[2] - zsort[3]);
    *zmom = (zsort[2] - zsort[3])/2-zsort[2];

    delete[] zsort;
}

void generate_thickness_map(float** mask, float** R, float lxbox, float lybox, float lzbox, float lcell, int ncell)
{
    int i, j;
    char outfile[100] = "tmap.dat";
    char zfile[100] = "zmap.dat";
    float t, zmom, lx1, lx2, ly1, ly2, lz1, lz2;
    float** p1, ** p2, ** p3, ** p4, ** p5, ** p6, ** tmap, **zmap;
    float* r1, * r2, * r3, * r4, * r5, * r6;
    int junk;

    p1 = float2D(3, 3, "p1");
    p2 = float2D(3, 3, "p2");
    p3 = float2D(3, 3, "p3");
    p4 = float2D(3, 3, "p4");
    p5 = float2D(3, 3, "p5");
    p6 = float2D(3, 3, "p6");
    tmap = float2D(ncell, ncell, "tmap");
    zmap = float2D(ncell, ncell, "zmap");
    r1 = new float[3]; r2 = new float[3]; r3 = new float[3]; r4 = new float[3];
    r5 = new float[3]; r6 = new float[3];

    define_box_planes(p1, p2, p3, p4, p5, p6, lxbox, lybox, lzbox, lcell);
    rotate_box_planes(R, p1, p2, p3, p4, p5, p6, lcell);

    for (i = 0; i < ncell; i++) {
        for (j = 0; j < ncell; j++) {
            lx1 = lcell / ncell * i;
            ly1 = lcell / ncell * j;
            lz1 = 0;
            lx2 = lx1; ly2 = ly1;
            lz2 = lcell;
            find_intersections(p1, p2, p3, p4, p5, p6, lx1, ly1, lz1, lx2, ly2, lz2, r1, r2, r3, r4, r5, r6);
            find_thickness(r1, r2, r3, r4, r5, r6, &t, &zmom);
            //TP June 2017: I swapped the indices to respect DM's +ve y is down system
            tmap[j][i] = t;
            zmap[j][i] = zmom;
        }
    }

    for (i = 0; i < ncell; i++) {
        for (j = 0; j < ncell; j++) {
            tmap[j][i] *= mask[i][j];
            zmap[j][i] *= mask[i][j];
        }
    }

    write_binary_2d_flarray(ncell, ncell, outfile, tmap);
    write_binary_2d_flarray(ncell, ncell, zfile, zmap);

    free_float2D(p1, 3);
    free_float2D(p2, 3);
    free_float2D(p3, 3);
    free_float2D(p4, 3);
    free_float2D(p5, 3);
    free_float2D(p6, 3);
    delete[] r1;
    delete[] r2;
    delete[] r3;
    delete[] r4;
    delete[] r5;
    delete[] r6;
    free_float2D(tmap, ncell);
    free_float2D(zmap, ncell);

}

void generate_thickness_maps(float** mask, float **tmap, float **zmap, float** R, float lxbox, float lybox, float lzbox, float lcell, int ncell)
{
    int i, j;
    char outfile[100];
    float t, zmom, lx1, lx2, ly1, ly2, lz1, lz2;
    float** p1, ** p2, ** p3, ** p4, ** p5, ** p6;
    float* r1, * r2, * r3, * r4, * r5, * r6;
    int junk;

    p1 = float2D(3, 3, "p1");
    p2 = float2D(3, 3, "p2");
    p3 = float2D(3, 3, "p3");
    p4 = float2D(3, 3, "p4");
    p5 = float2D(3, 3, "p5");
    p6 = float2D(3, 3, "p6");
    r1 = new float[3]; r2 = new float[3]; r3 = new float[3]; r4 = new float[3];
    r5 = new float[3]; r6 = new float[3];

    define_box_planes(p1, p2, p3, p4, p5, p6, lxbox, lybox, lzbox, lcell);
    rotate_box_planes(R, p1, p2, p3, p4, p5, p6, lcell);

    for (i = 0; i < ncell; i++) {
        for (j = 0; j < ncell; j++) {
            lx1 = lcell / ncell * i;
            ly1 = lcell / ncell * j;
            lz1 = 0;
            lx2 = lx1; ly2 = ly1;
            lz2 = lcell;
            find_intersections(p1, p2, p3, p4, p5, p6, lx1, ly1, lz1, lx2, ly2, lz2, r1, r2, r3, r4, r5, r6);
            find_thickness(r1, r2, r3, r4, r5, r6, &t, &zmom);
            //TP June 2017: I swapped the indices to respect DM's +ve y is down system
            tmap[j][i] = t;
            //2021: normalised centre of mass along the z-axis as the 1st moment or centroid
            zmap[j][i] = zmom;
        }
    }

    for (i = 0; i < ncell; i++) {
        for (j = 0; j < ncell; j++) {
            tmap[j][i] *= mask[i][j];
            zmap[j][i] *= mask[i][j];
        }
    }

    free_float2D(p1, 3);
    free_float2D(p2, 3);
    free_float2D(p3, 3);
    free_float2D(p4, 3);
    free_float2D(p5, 3);
    free_float2D(p6, 3);
    delete[] r1;
    delete[] r2;
    delete[] r3;
    delete[] r4;
    delete[] r5;
    delete[] r6;
}

int main()
{
    int i, j, npoints, ngrid, nbox, count, nang = 50;
    float* x, * y, * z, *tx, *ty, *tz;
    float* nx, * ny, * nz, ** R, **Rnew, **Rcombined, ** mask;
    float xcell, ycell, zcell, xmask;
    float ang, angrange;
    float **a, **tmap, **zmap, D;
    char infile[100];

    a = float2D(4, 4, "a");

    std::cout << "\n\nInput the number of points along one side of the output image:  ";
    std::cin >> npoints;

    tmap = float2D(npoints, npoints, "tmap");
    zmap = float2D(npoints, npoints, "zmap");

    /*
      std::cout<<"\n\nInput the number of grid points on each surface:  ";
      std::cin>>ngrid;
    */
    ngrid = npoints / 2;

    nbox = 6 * ngrid * ngrid;

    x = float1D(nbox, "x");
    y = float1D(nbox, "y");
    z = float1D(nbox, "z");
    tx = float1D(nbox, "x");
    ty = float1D(nbox, "y");
    tz = float1D(nbox, "z");
    R = float2D(3, 3, "R");
    Rnew = float2D(3, 3, "Rnew");
    Rcombined = float2D(3, 3, "Rnew");
    mask = float2D(npoints, npoints, "mask");

    std::cout << "\n\nInput the 1st box side length:  ";
    std::cin >> xcell;
    std::cout << "\n\nInput the 2nd box side length:  ";
    std::cin >> ycell;
    std::cout << "\n\nInput the 3rd box side length:  ";
    std::cin >> zcell;

    xmask = 0;

    /*
     do{
       std::cout<<"\n\nInput the mask side length (must be greater than that of cube):  ";
       std::cin>>xmask;
       if(xmask<=xcell)
         std::cout<<"\n\nMask side length must be greater than that of cube!  ";
     }
     while(xmask<=xcell);
     */
    xmask = npoints;

    nx = float1D(nbox, "nx");
    ny = float1D(nbox, "ny");
    nz = float1D(nbox, "nz");

    float a11, a12, a13;
    float a21, a22, a23;

    std::cout << "\n\nThis simulates the thickness of a rotating cube." << std::endl;
    std::cout << "\nStarting orientiation is in box_vectors.txt " << std::endl;
    std::cout << "\nDelete zmaps.dat and tmaps.dat before running, as these will be appended" << std::endl;
    
    float a31, a32, a33;

    std::ifstream fin;
    fin.open("box_vectors.txt");
    if (!fin.good()) {
        std::cout << "Can't find file box_vectors.txt " << std::endl;
        std::cout << "Enter a number to finish: " << std::endl;
        std::cin >> a11;
        exit(0);
    }

    fin.ignore(300, '\n');
    fin >> a11;
    fin >> a12;
    fin >> a13;
    fin.ignore(300, '\n');
    fin.ignore(300, '\n');
    fin >> a21;
    fin >> a22;
    fin >> a23;
    fin.ignore(300, '\n');
    fin.ignore(300, '\n');
    fin >> a31;
    fin >> a32;
    fin >> a33;

    fin.close();

    //Note: R is a generic rotation matrix, not just rotation about a
    //single axis.  The a11, a12 components specify rotation about x, y and z axes
    //as direction cosines.
    //Hence all of these components will generally be non-zero.

    R[0][0] = a11; R[0][1] = a21; R[0][2] = a31;
    R[1][0] = a12; R[1][1] = a22; R[1][2] = a32;
    R[2][0] = a13; R[2][1] = a23; R[2][2] = a33;

    count = 0;

    generate_box_points(x, y, z, xcell, ycell, zcell, xmask, ngrid,count);
    rotate_xyz(x, y, z, nx, ny, nz, R, xcell, ycell, zcell, nbox);
    generate_box_shadow_mask(nx, ny, xcell, ycell, mask, xmask, ngrid, npoints);
    ////  write_binary_2d_flarray(npoints,npoints,"mask.dat",mask);
    generate_thickness_map(mask, R, xcell, ycell, zcell, xmask, npoints);

    for (j = 0; j < count; j++) {
        tx[j] = x[j];
        ty[j] = y[j];
        tz[j] = z[j];
    }

    //angrange = 0.6; // +/-5 degrees of rotation

    std::cout << "Enter the angular range in degrees: " << std::endl;
    std::cin >> angrange;

    for (i = 0; i < nang; i++) {
        //reset box points
        ang = angrange * 1.0*(i - 1.0*nang / 2) / (nang/2);
        ang /= 180;
        ang *= PI;
        cout << "\n\nang: " << ang << endl;
        for(j=0;j<count;j++){
            x[j] = tx[j];
            y[j] = ty[j];
            z[j] = tz[j];
        }

        //rotate about the y-axis in a clockwise direction when looking towards to origin
 
        Rnew[0][0] = cos(ang); Rnew[0][1] = 0.0; Rnew[0][2] = -sin(ang);
        Rnew[1][0] = 0.0;      Rnew[1][1] = 1.0; Rnew[1][2] = 0.0;
        Rnew[2][0] = sin(ang); Rnew[2][1] = 0.0; Rnew[2][2] = cos(ang);

        matrix_multiply(Rnew, R, Rcombined);

        rotate_xyz(x, y, z, nx, ny, nz, Rcombined, xcell, ycell, zcell, nbox);
        generate_box_shadow_mask(nx, ny, xcell, ycell, mask, xmask, ngrid, npoints);
        ////  write_binary_2d_flarray(npoints,npoints,"mask.dat",mask);
        generate_thickness_maps(mask, tmap, zmap, Rcombined, xcell, ycell, zcell, xmask, npoints);
        //generate_thickness_map(mask, Rcombined, xcell, ycell, zcell, xmask, npoints);
        write_binary_3d_flarray(1, npoints, npoints, "tmaps.dat", "tmaps.dat", &tmap);
        write_binary_3d_flarray(1, npoints, npoints, "zmaps.dat", "zmaps.dat", &zmap);
    }

    exit(0);

    return 0;
}
