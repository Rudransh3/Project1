
#define _USE_MATH_DEFINES

#include <stdio.h>
#include <stdlib.h>
#include "math.h"
#include "rir_generator_core.h"

#define ROUND(x) ((x)>=0?(long)((x)+0.5):(long)((x)-0.5))

#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif

double sinc(double x)
{
    if (x == 0)
        return(1.);
    else
        return(sin(x)/x);
}

double sim_microphone(double x, double y, double z, double* microphone_angle, char mtype)
{
    if (mtype=='b' || mtype=='c' || mtype=='s' || mtype=='h')
    {
        double gain, vartheta, varphi, rho;

        // Polar Pattern         rho
        // ---------------------------
        // Bidirectional         0
        // Hypercardioid         0.25
        // Cardioid              0.5
        // Subcardioid           0.75
        // Omnidirectional       1

        switch(mtype)
        {
        case 'b':
            rho = 0;
            break;
        case 'h':
            rho = 0.25;
            break;
        case 'c':
            rho = 0.5;
            break;
        case 's':
            rho = 0.75;
            break;
        default:
            rho = 1;
        };

        vartheta = acos(z/sqrt(pow(x,2)+pow(y,2)+pow(z,2)));
        varphi = atan2(y,x);

        gain = sin(M_PI/2-microphone_angle[1]) * sin(vartheta) * cos(microphone_angle[0]-varphi) + cos(M_PI/2-microphone_angle[1]) * cos(vartheta);
        gain = rho + (1-rho) * gain;

        return gain;
    }
    else
    {
        return 1;
    }
}

void computeRIR(double* imp, double c, double fs, double* rr, int nMicrophones, int nSamples, double* ss, double* LL, double* beta, char microphone_type, int nOrder, double* microphone_angle, int isHighPassFilter){

    // Temporary variables and constants (high-pass filter)
    const double W = 2*M_PI*100/fs; // The cut-off frequency equals 100 Hz
    const double R1 = exp(-W);
    const double B1 = 2*R1*cos(W);
    const double B2 = -R1 * R1;
    const double A1 = -(1+R1);
    double       X0;
    double       Y[3];

    // Temporary variables and constants (image-method)
    const double Fc = 0.5; // The normalized cut-off frequency equals (fs/2) / fs = 0.5
    const int    Tw = 2 * ROUND(0.004*fs); // The width of the low-pass FIR equals 8 ms
    const double cTs = c/fs;
    double*      LPI = new double[Tw];
    double       r[3];
    double       s[3];
    double       L[3];
    double       Rm[3];
    double       Rp_plus_Rm[3];
    double       refl[3];
    double       fdist,dist;
    double       gain;
    int          startPosition;
    int          n1, n2, n3;
    int          q, j, k;
    int          mx, my, mz;
    int          n;

    s[0] = ss[0]/cTs; s[1] = ss[1]/cTs; s[2] = ss[2]/cTs;
    L[0] = LL[0]/cTs; L[1] = LL[1]/cTs; L[2] = LL[2]/cTs;

    for (int idxMicrophone = 0; idxMicrophone < nMicrophones ; idxMicrophone++)
    {
        // [x_1 x_2 ... x_N y_1 y_2 ... y_N z_1 z_2 ... z_N]
        r[0] = rr[idxMicrophone + 0*nMicrophones] / cTs;
        r[1] = rr[idxMicrophone + 1*nMicrophones] / cTs;
        r[2] = rr[idxMicrophone + 2*nMicrophones] / cTs;

        n1 = (int) ceil(nSamples/(2*L[0]));
        n2 = (int) ceil(nSamples/(2*L[1]));
        n3 = (int) ceil(nSamples/(2*L[2]));

        // Generate room impulse response
        for (mx = -n1 ; mx <= n1 ; mx++)
        {
            Rm[0] = 2*mx*L[0];

            for (my = -n2 ; my <= n2 ; my++)
            {
                Rm[1] = 2*my*L[1];

                for (mz = -n3 ; mz <= n3 ; mz++)
                {
                    Rm[2] = 2*mz*L[2];

                    for (q = 0 ; q <= 1 ; q++)
                    {
                        Rp_plus_Rm[0] = (1-2*q)*s[0] - r[0] + Rm[0];
                        refl[0] = pow(beta[0], abs(mx-q)) * pow(beta[1], abs(mx));

                        for (j = 0 ; j <= 1 ; j++)
                        {
                            Rp_plus_Rm[1] = (1-2*j)*s[1] - r[1] + Rm[1];
                            refl[1] = pow(beta[2], abs(my-j)) * pow(beta[3], abs(my));

                            for (k = 0 ; k <= 1 ; k++)
                            {
                                Rp_plus_Rm[2] = (1-2*k)*s[2] - r[2] + Rm[2];
                                refl[2] = pow(beta[4],abs(mz-k)) * pow(beta[5], abs(mz));

                                dist = sqrt(pow(Rp_plus_Rm[0], 2) + pow(Rp_plus_Rm[1], 2) + pow(Rp_plus_Rm[2], 2));

                                if (abs(2*mx-q)+abs(2*my-j)+abs(2*mz-k) <= nOrder || nOrder == -1)
                                {
                                    fdist = floor(dist);
                                    if (fdist < nSamples)
                                    {
                                        gain = sim_microphone(Rp_plus_Rm[0], Rp_plus_Rm[1], Rp_plus_Rm[2], microphone_angle, microphone_type)
                                            * refl[0]*refl[1]*refl[2]/(4*M_PI*dist*cTs);

                                        for (n = 0 ; n < Tw ; n++)
                                        {
                                            const double t = (n-0.5*Tw+1) - (dist-fdist);
                                            LPI[n] = 0.5 * (1.0 + cos(2.0*M_PI*t/Tw)) * 2.0*Fc * sinc(M_PI*2.0*Fc*t);                                            
                                        }
                                        startPosition = (int) fdist-(Tw/2)+1;
                                        for (n = 0 ; n < Tw; n++)
                                            if (startPosition+n >= 0 && startPosition+n < nSamples)
                                                imp[idxMicrophone + nMicrophones*(startPosition+n)] += gain * LPI[n];
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        // 'Original' high-pass filter as proposed by Allen and Berkley.
        if (isHighPassFilter == 1)
        {
            for (int idx = 0 ; idx < 3 ; idx++) {Y[idx] = 0;}
            for (int idx = 0 ; idx < nSamples ; idx++)
            {
                X0 = imp[idxMicrophone+nMicrophones*idx];
                Y[2] = Y[1];
                Y[1] = Y[0];
                Y[0] = B1*Y[1] + B2*Y[2] + X0;
                imp[idxMicrophone+nMicrophones*idx] = Y[0] + A1*Y[1] + R1*Y[2];
            }
        }
    }

    delete[] LPI;
}
