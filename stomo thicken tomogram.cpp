#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>
#include <assert.h>

#define PI 3.141592654

using namespace std;

void WRITE_MERGED_PDB(int npoints1, int npoints2, float* x1, float* y1, float* z1, float* x2, float* y2, float* z2)
{
    int i;

    ofstream fout("merged_atoms.pdb");

    fout.setf(ios::fixed, ios::floatfield);

    for (i = 0; i < npoints1; i++) {
        fout << "ATOM" << setw(7) << i + 1 << setw(3) << "B";
        fout << setprecision(3) << setw(24) << x1[i];
        fout << setprecision(3) << setw(8) << y1[i];
        fout << setprecision(3) << setw(8) << z1[i] << endl;
    }

    for (i = 0; i < npoints2; i++) {
        fout << "ATOM" << setw(7) << i + 1 + npoints1 << setw(3) << "N";
        fout << setprecision(3) << setw(24) << x2[i];
        fout << setprecision(3) << setw(8) << y2[i];
        fout << setprecision(3) << setw(8) << z2[i] << endl;
    }

    fout.close();
}

void WRITE_PDB(int npoints1, float* x1, float* y1, float* z1)
{
    int i;

    ofstream fout("junk_atoms.pdb");

    for (i = 0; i < npoints1; i++) {
        fout.setf(ios::fixed, ios::floatfield);
        fout << "ATOM" << setw(7) << i + 1 << setw(3) << "B";
        fout << setprecision(3) << setw(24) << x1[i];
        fout << setprecision(3) << setw(8) << y1[i];
        fout << setprecision(3) << setw(8) << z1[i] << endl;
    }

    fout.close();
}


void read_x_y_z(int* natoms, float* x, float* y, float* z, char* namex, char* namey, char* namez)
{
    int i, count = 0;

    ifstream finx(namex);
    ifstream finy(namey);
    ifstream finz(namez);

    count = 0;

    for (i = 0; i < *natoms; i++) {
        if (finx >> x[i])
            count++;
        else
            break;
        finy >> y[i];
        finz >> z[i];
    }

    *natoms = count;

    finx.close();
    finy.close();
    finz.close();
}

void read_zerrs(int* natoms, float* zerrs, const char* name)
{
    int i, count, length;
    float* buffer;

    ifstream fin;
    fin.open(name, ios::binary);

    // get length of file:
    fin.seekg(0, ios::end);
    length = fin.tellg();
    fin.seekg(0, ios::beg);
    // allocate memory:
    buffer = new float[length];
    // read data as a block:
    fin.read((char*)buffer, length);
    fin.close();
    *natoms = length / 4;
    count = 0;

    for (i = 0; i < (*natoms) ; i ++) {
        zerrs[i] = buffer[i];
    }

    delete[] buffer;
}


void read_xyz(int* natoms, float* x, float* y, float* z, char* name)
{
    int i, count, length;
    float* buffer;

    ifstream fin;
    fin.open(name, ios::binary);

    // get length of file:
    fin.seekg(0, ios::end);
    length = fin.tellg();
    fin.seekg(0, ios::beg);
    // allocate memory:
    buffer = new float[length];
    // read data as a block:
    fin.read((char*)buffer, length);
    fin.close();
    *natoms = length / 12;
    count = 0;

    for (i = 0; i < (*natoms) * 3; i += 3) {
        x[count] = buffer[i + 0];
        y[count] = buffer[i + 1];
        z[count] = buffer[i + 2];
        count++;
    }

    delete[] buffer;
}

// Write a 1D array of type float
void write_binary_1d_flarray(unsigned long size, char* file, float* array)
{
    //Open the binary file with appropriate mode bits (ios flags)
    ofstream fout, flog;
    fout.open(file, ios::out | ios::binary);
    flog.open("logfile.txt", ios::app);

    //...........................................................

    cout << "\n\n.... starting xyz file write." << endl;
    flog << "\n.... starting xyz file write.";
    fout.write((char*)array, size * sizeof(float));
    fout.close();
    flog.close();
    //...........................................
    cout << ".... finished file write." << endl;
}

void x_y_z_to_xyz(int natoms, float* x, float* y, float* z, float* xyz)
{
    int i, count;

    count = 0;

    for (i = 0; i < natoms * 3; i += 3) {
        xyz[count + 0] = x[i / 3];
        xyz[count + 1] = y[i / 3];
        xyz[count + 2] = z[i / 3];
        count += 3;
    }
}

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

void write_binary_3d_flarray(int nstart, int nend, int N, int M,
    const char* file, const char* fn, float*** array)
{
    // Write to a *.dat floating point file, a 3D image stack of type float
    int i, j, k, nimgs;
    unsigned long count, size;
    float* buffer;

    //cout<<"\n\nBeginning file write for "<<fn<<" ..."<<endl;

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

void project_point_cloud_to_tilt_series(float* xyz, float* zerrs, float* angles, float tiltx, int accumflag, int nimgs,
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
                if (accumflag) {
                    //accumulate points by summing into x, jrot bin
                    tiltimage[int(x)][int(jrot)] += zerrs[int(i / 3)];
                }
                else {
                    //replace points each time, don't accumulate
                    tiltimage[int(x)][int(jrot)] = zerrs[int(i / 3)];
                }
            }
        }
        write_binary_3d_flarray(a, a + 1, N, M, "thick_projections.dat", "thick_projections.dat", &tiltimage);
    }

}

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

void point_cloud_to_tomogram(float* xyz, float* zerrs, int accumflag, int N, int M, int ndata)
{
    //Function quantises xyz coordinates in the point cloud
    //and outputs them on a 3D grid.

    int i, j, k, K, count;
    float x, y, z;
    float*** tomogram;
    
    ///Note: due to memory limitations, I will halve the resolution
    //assign 3rd tomogram dimenion K to the x axis field of view
    K = int(M/2);
    N = int(N / 2);
    M = int(M / 2);

    tomogram = float3D(K, N, M, "tomogram");
 
    for (i = 0; i < N; i++) {
        for (j = 0; j < M; j++) {
            for (k = 0; k < K; k++) {
                tomogram[k][i][j] = 0.0E00;
            }
        }
    }

    count = 0;
    for (i = 0; i < ndata; i += 3) {
       x = xyz[i + 0]/2;
       y = xyz[i + 1]/2;
       z = xyz[i + 2]/2;
       if (int(x) >= 0 && int(x) < N && int(y) < M && int(y) >= 0 && int(z) < K && int(z) >= 0) {
           if (accumflag) {
               //accumulate points in x, y, z bin
               tomogram[int(z)][int(x)][int(y)] += zerrs[int(i / 3)];
           }
           else {
               //just replace points, don't accumulate
               tomogram[int(z)][int(x)][int(y)] = zerrs[int(i / 3)];
           }             
           count++;
       }
    }

    cout << "There are " << count << " points in the tomogram " << endl;
    
    for (i = 0; i < K; i++) {
      write_binary_3d_flarray(i, i + 1, N, M, "tomogram.dat", "tomogram.dat", &tomogram[i]);
    }
}

int main()
{
    int m, i, j, natoms1, natoms2, natomsz, natomst, disp, count, nimgs, N, M, zflag, accumflag;
    float tiltx, zav, * x1, * y1, * z1, * x2, * y2, * z2, * array;
    float* x3, * y3, * z3, * x4, * y4, * z4;
    float* xt, * yt, * zt, * nx, * ny, * nz, * xyzt, *angles, *zerrs, *zerrs_thick;
    char* f1, * f2, * f3, * f4;

    natoms1 = natoms2 = natomst = 20000000;

    x1 = new float[natoms1];
    y1 = new float[natoms1];
    z1 = new float[natoms1];
    x2 = new float[natoms2];
    y2 = new float[natoms2];
    z2 = new float[natoms2];
    xt = new float[natomst];
    yt = new float[natomst];
    zt = new float[natomst];
    nx = new float[natomst];
    ny = new float[natomst];
    nz = new float[natomst];
    zerrs = new float[natomst];
    zerrs_thick = new float[natomst];

    for (i = 0; i < natomst; i++) {
        x1[i] = 0.0E00;
        y1[i] = 0.0E00;
        z1[i] = 0.0E00;
        x2[i] = 0.0E00;
        y2[i] = 0.0E00;
        z2[i] = 0.0E00;
        xt[i] = 0.0E00;
        yt[i] = 0.0E00;
        zt[i] = 0.0E00;
        nx[i] = 0.0E00;
        ny[i] = 0.0E00;
        nz[i] = 0.0E00;
        zerrs[i] = 0.0E00;
        zerrs_thick[i] = 0.0E00;

    }

    xyzt = new float[3 * natomst];

    f1 = new char[100];
    f2 = new char[100];

    cout << "Input the name of the xyz file: ";
    cin >> f1;

    cout << "Input the name of the normals file: ";
    cin >> f2;

    cout << "Input 1 for zerrs.dat coloring of points, 0 for no: ";
    cin >> zflag;

    cout << "Input 1 for accumulating points, 0 for no: ";
    cin >> accumflag;

    read_xyz(&natoms1, x1, y1, z1, f1);
    read_xyz(&natoms2, x2, y2, z2, f2);

    if (natoms1 != natoms2) {
        cout << "Unequal numbers of points between xyz and normals." << endl;
        cout << "Enter anything to quit:  " << endl;
        cin >> i;
        exit(0);
    }

    //read natoms1 numbers of points for zerrors.dat but this
    //is then overwritten to check that natoms1 = natomsz.
    natomsz = natoms1;
    if (zflag) {
        read_zerrs(&natomsz, zerrs, "zerrors.dat");
        cout << "\n\n zerrs has " << natomsz << " points." << endl;
    }
    else { //fill zerrs with ones if choosing not to colour the point cloud
        for (i = 0; i < natoms1; i++) {
            zerrs[i] = 1.0E00;
        }
    }

    ifstream fin("stomo_thicken_parameters.txt");

    fin >> nimgs; fin.ignore(300, '\n');
    fin >> N; fin.ignore(300, '\n');
    fin >> M; fin.ignore(300, '\n');
    fin >> tiltx; fin.ignore(300, '\n');

    angles = new float[nimgs];

    for (i = 0; i < nimgs; i++) {
        fin >> angles[i];
        fin.ignore(300, '\n');
        angles[i] = angles[i] / 180 * PI;
    }
    fin.close();


    cout << "\n\n # of xyz points is: " << natoms1;
    cout << "\n\n # of norms is: " << natoms2;

    cout << "\n\nInput the maximum normal displacement for parallel surface (pixels): ";
    cin >> disp;

    zav = 0.0E00;
    for (i = 0; i < natoms1; i++) { //natoms1 should equal natoms2!
        xt[i] = x1[i];
        nx[i] = x2[i];
        yt[i] = y1[i];
        ny[i] = y2[i];
        zt[i] = z1[i];
        nz[i] = z2[i];
        zav += zt[i];
    }
    zav /= natoms1;

    //x_y_z_to_xyz(3 * natomst, xt, yt, zt, xyzt);
    //x_y_z_to_xyz(3 * natomst, nxt, nyt, nzt, nxyzt);

    count = 0;

    for(j=-disp; j<=disp;j++){
        cout << count << endl;
        for (i = 0; i < natoms1; i++) { //natoms1 must be the same as natoms2
            x2[count] = xt[i] + j * nx[i]; //use vacant x2 array for expanding along normals                         
            y2[count] = yt[i] + j * ny[i];
            z2[count] = zt[i] + j * nz[i];
            zerrs_thick[count] = zerrs[i];
            count++;
        }
    }

    cout << "writing " << count << " points as tilt projections" << endl;

    x_y_z_to_xyz(count, x2, y2, z2, xyzt);
    project_point_cloud_to_tilt_series(xyzt, zerrs_thick, angles, tiltx, accumflag, nimgs, N, M, count);

    for (i = 0; i < count; i++) { //natoms1 must be the same as natoms2
      z2[i] = z2[i] -zav+M/2; //recentre z cooridnate system (stomo outputs it off centre)
    }
    //convert again to xyz contiguous array, after having adjusted the z height
    x_y_z_to_xyz(count, x2, y2, z2, xyzt);

    cout << "writing tomogram " << endl;

    point_cloud_to_tomogram(xyzt, zerrs_thick, accumflag, N, M, count);

    return 0;
}
