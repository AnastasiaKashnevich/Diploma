#include <iostream>
#include "iomanip"
#include "math.h"
#include <cstdio>
#include <string>
#include "fstream"
#include <vector>
#include <stdio.h>


const double PI(3.1415926535897932384626433832795);
using namespace std;

double randDouble() { return (double)rand() / RAND_MAX; }
//определяем ф
double value_fi(double tmp) { return 2 * PI * tmp; }
//определяем косинус т
double value_t(double k, double b) {
    if (k > 0) {
        return  acos((1 / (2 * k)) * (1 + pow(k, 2) - pow(((1 - pow(k, 2)) / (1 - k + 2 * k * b)), 2)));
    }
    else
    {
        return  acos((2 * b - 1));
    }
}

double check(double cx, double cy, double cz, double x, double y, double z) {

    double x1 = pow((x - cx), 2);
    double y1 = pow((y - cy), 2);
    double z1 = pow((z - cz), 2);

    // расстояние между центром
    // и заданная точка
    return (x1 + y1 + z1);
}

int sign(double a) {
    if (a > 0) return 1;
    if (a < 0) return -1;
    return 0;
}

void readFileA(std::vector<double>& A, string filename) {
    std::ifstream in(filename, std::ios::in);
    std::string line;
    int countNotZeroA = 0;

    if (in.is_open())
    {
        while (getline(in, line))
        {
            double currentA = ceil(stod(line)); 
            A.push_back(currentA);
            if (currentA != 0) {
                countNotZeroA++;
            }
        }
    }

    in.close();
    in.clear();

    // debug
    std::cout << "countNotZeroA = " << countNotZeroA << std::endl;
}

int main()
{
    setlocale(LC_ALL, "Russian");
    double Ms, Ma, g, lc, l, P, R=0.5, Rc=1;
    double Xsqr, Ysqr, ksi, ksi1, fi, cos_t, sin_t;
    double Xmax = 20, Ymax = 20, Zmax = 10, N = 200, N_photons = 1000000;
    int    Nx = 200, Ny = 200, Nz = 200, iX, iY, iZ, iA, ph;

    cout << "0 - засветка" << endl;
    cout << "1 - флуоресценция" << endl;
    cout << "Введите 0 или 1:";
    cin >> ph;
    

    if (ph == 0)
    {
        //параметры поглощения
        //для 405 нм
        //Ms = 37, Ma = 0.89, g = 0.8;
        //для 660 нм
        //Ms = 14, Ma = 0.15, g = 0.8;
        // начальные 
        Ms = 10, Ma = 0.01, g = 0.9;

        vector <double> A(Nx*Ny*Nz);

        //средняя длина свободного пробега
        lc = 1./ (Ms + Ma);


        ofstream outXX;
        ofstream outYY;
        ofstream outZZ;
    
        outYY.open("valueYY.txt");
        outZZ.open("nXnYnZ.txt");

        for (int i = 0; i < N_photons; i++)
        {
            double lastIA = 0;
            double lastX = 0, lastY = 0, lastZ = 0;

            if (i % 10000 == 0)
            {
                cout << 100 * (i / N_photons) << "% complete" << endl;
            }
        
            P = 1;

            // Quantity change direction
            int inm = 0;
            double Mold[3] = { 0, 0, 1 }, M[3] = { 0, 0, 1 }, X[3] = { 0, 0, 0 }, S[3] = {0, 0, 1};
            bool exitP = false;

            //Масштаб засветки 
            /*
            Xsqr = (randDouble() - 0.5) * 5;
            X[0] = Xsqr;
            Ysqr = (randDouble() - 0.5) * 5;
            X[1] = Ysqr;*/

            //(вариант для луча) 
            while ((-Xmax < X[0]) && (X[0] < Xmax) && (-Ymax < X[1]) && (X[1] < Ymax) && (0 <= X[2]) && (X[2] < Zmax))
            {
                ksi = randDouble();
                ksi1 = randDouble();
                //сл длина свободного пробега
                l = -log(1 - ksi) * lc;
                //угол ф
                fi = value_fi(randDouble());
                //косинус и синус угла т
                sin_t = sin(value_t(g, ksi1));
                cos_t = cos(value_t(g, ksi1));
                //направляющие косинусы вектора скорости
                if (fabs(M[2]) <= 0.99999)
                {
                    Mold[0] = M[0];
                    Mold[1] = M[1];
                    Mold[2] = M[2];

                    M[0] = (sin_t / sqrt(1 - pow(Mold[2], 2))) * (Mold[0] * Mold[2] * cos(fi) - Mold[1] * sin(fi)) + Mold[0] * cos_t;
                    M[1] = (sin_t / sqrt(1 - pow(Mold[2], 2))) * (Mold[1] * Mold[2] * cos(fi) + Mold[0] * sin(fi)) + Mold[1] * cos_t;
                    M[2] = (-sin_t) * cos(fi) * sqrt(1 - pow(Mold[2], 2)) + Mold[2] * cos_t;
                }
                else
                {
                    M[0] = cos(fi) * sin_t;
                    M[1] = sin(fi) * sin_t;
                    M[2] = sign(M[2]) * cos_t;
                }
                //новые координаты фотона
                for (int j = 0; j < 3; j++)
                {
                    X[j] = X[j] + l * M[j];
                }
                //проверка находится ли точка в сфере
                double ras = check(S[0], S[1], S[2], X[0], X[1], X[2]);
                if( ras <= (R*R))
                { 
                    iX = floor(((X[0] + (Xmax / 2)) / Xmax) * (N));
                    iY = floor(((X[1] + (Ymax / 2)) / Ymax) * (N));
                    iZ = floor((X[2] / Zmax) * N);
                    iA = iX + Nx * iY + Nx * Ny * iZ;
             
                    if ((0 < iX && iX < N) && (0 < iY && iY < N) && (0 < iZ && iZ < N))
                    { 
                        A[iA] += P * (Ma / (Ms + Ma));
                        lastIA = iA;  
                    }
                } 
                lastX = X[0]; lastY = X[1]; lastZ = X[2];
                //уменьшение веса фотона
                P = P * (Ms / (Ms + Ma));
                inm++;

                if (P <= 0.0001) {
                    exitP = true;
                    break;
                }
                
            }
        
        }
        outXX.open("valueXX.txt");
        double lastA = 0;
        for (int k=0; k< Nx*Ny*Nz; k++)
        {
            lastA = A[k];
            outXX << lastA << "\n";
        }
        outXX.close();
        A.clear();
        return 0;
    }
    else
    {
        // начальные 
        Ms = 10, Ma = 0.01, g = 0.9;
        //параметры Флуоресценции (всегда считаем для 760 нм)
        //Ms = 12, Ma = 0.13, g = 0.8;

        vector <double> A(Nx* Ny* Nz);
        vector <double> F(Nx* Ny);
        double prevX[3] = { 0,0,0 }, xp, yp;
        double lastIA = 0;
        double Xnew[3] = { 0,0,0 };
        double Q = 0;
        int c=100, ixp, iyp, ifp;
        std::vector<double> B;
        readFileA(B, "valueXX.txt");
        ofstream outFileLastA;
        outFileLastA.open("outFileLastA.txt");
        ofstream Number;
        Number.open("Number.txt");
        
        for (int i = 0; i < B.size(); i++) {

            
            double currentB = B[i]*c;
            double j = 1.0;

            while (currentB > 0) {
                if (currentB < 1) {
                    j = currentB;
                }

                int inm = 0;
                double Mold[3] = { 0, 0, 1 }, M[3] = { 0, 0, 1 };
                bool exitP = false;
                double Ps;

                lc = 1. / (Ms + Ma);

                Xnew[0] = (2 * Xmax) / Nx * ((i) % Nx) + (Xmax) / Nx - Xmax;

                Xnew[1] = (2 * Ymax) / Ny * (int((i - Xnew[0]) / Nx) % Ny) + Ymax / Ny - Ymax;

                Xnew[2] = (2 * Zmax) / Nz * int((i - Xnew[0] - Nx * Xnew[1]) / (Nx * Ny)) + Zmax / (2 * Nz);
                //cout << "x = " << Xnew[0] << "; y = " << Xnew[1] << "; z = " << Xnew[2] << endl;

                Ps = j;
                while ((-Xmax < Xnew[0]) && (Xnew[0] < Xmax) && (-Ymax < Xnew[1]) && (Xnew[1] < Ymax) && (0 <= Xnew[2]) && (Xnew[2] < Zmax))
                {
                    ksi = randDouble();
                    //сл длина свободного пробега
                    l = -log(1 - ksi) * lc;

                    //угол ф
                    fi = PI * (2.0 * rand() / RAND_MAX - 1.0);

                    // угл т
                    double theta = asin(2.0 * rand() / RAND_MAX - 1.0);
                    if (theta < 0)  theta += PI ;

                    //направляющие косинусы вектора скорости
                    M[0] = sin(theta) * cos(fi);
                    M[1] = sin(theta) * sin(fi);
                    M[2] = cos(theta);

                    //новые координаты фотона
                    for (int j = 0; j < 3; j++)
                    {
                        //сохраняем предыдущие коорлинаты
                        prevX[j] = Xnew[j];

                        Xnew[j] = Xnew[j] + l * M[j];
                    }

                    //коорлинаты на поверхности
                    xp = ((Xnew[0] - prevX[0]) * (Xnew[2] / (prevX[2] - Xnew[2]))) + Xnew[0];
                    yp = ((Xnew[1] - prevX[1]) * (Xnew[2] / (prevX[2] - Xnew[2]))) + Xnew[1];


                    //проверка на поверхности
                    if (Xnew[2] < 0) 
                    {
                        ixp = floor(((xp + (Xmax / 2)) / Xmax) * (N));
                        iyp = floor(((yp + (Ymax / 2)) / Ymax) * (N));
                        ifp = ixp + Nx * iyp;

                        double K = (xp * xp) + (yp * yp);
                        //считаем вес внутри круга
                        if (K < (Rc * Rc))
                        {
                            Q += Ps;
                        }
                    }

                   
                    //построение карты фл-ии
                    if ((0 < ixp && ixp < N) && (0 < iyp && iyp < N))
                    {
                        F[ifp] += Ps;    
                    }
                    
                   /**/
                    //общий случай внутри ткани
                    iX = floor(((Xnew[0] + (Xmax / 2)) / Xmax) * (N));
                    iY = floor(((Xnew[1] + (Ymax / 2)) / Ymax) * (N));
                    iZ = floor((Xnew[2] / Zmax) * N);
                    iA = iX + Nx * iY + Nx * Ny * iZ;

                    if ((0 < iX && iX < N) && (0 < iY && iY < N) && (0 < iZ && iZ < N))
                    {
                        A[iA] += Ps * (Ma / (Ms + Ma));
                        lastIA = iA;
                    }
                    

                    Ps = Ps * (Ms / (Ms + Ma));
                    inm++;

                    if (Ps <= 0.0001) {
                        exitP = true;
                        break;
                    }

                    
                }
                currentB--;
                //cout << Q << endl;
            }
            
        }
        Number << Q;
        Number.close();
        cout<< "final Q=" << Q << endl;
        /**/
        ofstream outOne;
        outOne.open("1_out.txt");
        for (int k = 0; k < Nx * Ny * Nz; k++)
        {
            outOne << A[k] << "\n";
        }
        outOne.close();
        for (int l = 0; l < Nx * Ny; l++)
        {
            outFileLastA << F[l] << "\n";
        }
        outFileLastA.close();
        A.clear();
        // return 0;
        cout << "пока ничего";
    }
    return 0;
}

