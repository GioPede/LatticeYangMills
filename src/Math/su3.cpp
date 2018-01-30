#include "Math/su3.h"
#include "Math/complex.h"
#include <ctime>
#include <cmath>
#include <cstdio>
#include <utility>
#include <random>

SU3 I, Q2, Q3;
SU3 F0(0.0), F1(0.0), F2(0.0);


SU3::SU3() noexcept{
}

SU3::SU3(double value) noexcept{
    for(auto& i : mat)
        i = value;
}

SU3::SU3(SU3 const& source) noexcept : mat(source.mat){
}

SU3::SU3(SU3&& source) noexcept : mat(source.mat){
}

// Hermitean Conjugate
SU3 operator~(const SU3& a){
    SU3 c;
    for(int i = 0; i < 3; i++){
        for(int j = 0; j < 3; j++){
            c.mat[2 * (3*i + j)] = a.mat[2 * (3*j + i)];
            c.mat[2 * (3*i + j) + 1] = -a.mat[2 * (3*j + i) + 1];
        }
    }
    return std::move(c);
}

SU3 operator~(SU3&& a){
    SU3 c;
    for(int i = 0; i < 3; i++){
        for(int j = 0; j < 3; j++){
            c.mat[2 * (3*i + j)] = a.mat[2 * (3*j + i)];
            c.mat[2 * (3*i + j) + 1] = -a.mat[2 * (3*j + i) + 1];
        }
    }
    return std::move(c);
}

void SU3::printSU3(){
    int idx = 0;
    for(int i = 0; i < 3; i++){
        for(int j = 0; j < 3; j++){
            printf("( %2.5g, %2.5g )\t", mat[idx], mat[idx+1]);
            idx += 2;
        }
        printf("\n");
    }
    printf("\n");
}

void SU3::setSU3Identity(){
    for(int i = 0; i < 18; i++)
        mat[i] = 0;
    mat[0] = 1;
    mat[8] = 1;
    mat[16] = 1;
}

void SU3::setSU3Zero(){
    for(int i = 0; i < 18; i++)
        mat[i] = 0;
}

double SU3::realTrace(){
    return mat[0] + mat[8] + mat[16];
}

double SU3::imagTrace(){
    return mat[1] + mat[9] + mat[17];
}
double getMultSumTrace(SU3 &a, SU3 &b){
    return a.mat[0] * b.mat[0] - a.mat[1]*b.mat[1] + a.mat[2] * b.mat[6] - a.mat[3]*b.mat[7] + a.mat[4] * b.mat[12] - a.mat[5]*b.mat[13]
         + a.mat[6] * b.mat[2] - a.mat[7]*b.mat[3] + a.mat[8] * b.mat[8] - a.mat[9]*b.mat[9] + a.mat[10] * b.mat[14] - a.mat[11]*b.mat[15]
         + a.mat[12] * b.mat[4] - a.mat[13]*b.mat[5] +  a.mat[14] * b.mat[10] - a.mat[15]*b.mat[11] + a.mat[16] * b.mat[16] - a.mat[17]*b.mat[17];

}

double getMultSumTrace(SU3&& a, SU3 &b){
    return a.mat[0] * b.mat[0] - a.mat[1]*b.mat[1] + a.mat[2] * b.mat[6] - a.mat[3]*b.mat[7] + a.mat[4] * b.mat[12] - a.mat[5]*b.mat[13]
         + a.mat[6] * b.mat[2] - a.mat[7]*b.mat[3] + a.mat[8] * b.mat[8] - a.mat[9]*b.mat[9] + a.mat[10] * b.mat[14] - a.mat[11]*b.mat[15]
         + a.mat[12] * b.mat[4] - a.mat[13]*b.mat[5] +  a.mat[14] * b.mat[10] - a.mat[15]*b.mat[11] + a.mat[16] * b.mat[16] - a.mat[17]*b.mat[17];

}

SU3 exp(const SU3& Q){
    Q2 = Q*Q;
        Q3 = Q2*Q;

        double c0 = Q3.realTrace() / 3.0;
        double c1 = Q2.realTrace() / 2.0;

        bool c0negative = false;
        if(c0 < 0){
            c0negative = true;
            c0 = -c0;
        }

        double c1sqrt = sqrt(c1);
        double oneoversqrt3 = 1.0/sqrt(3);
        double c0max = c1*c1sqrt*oneoversqrt3*2.0 / 3.0;
        double theta = acos(c0/c0max);
        double w = c1sqrt*sin(theta/3.0);
        double u = c1sqrt*oneoversqrt3*cos(theta/3.0);
        double xi;
        double u2 = u*u;
        double w2 = w*w;

        if(fabs(w) < 0.05){
            xi = 1.0 - w2/6.0 * (1 - w2/20.0 * (1 - w2/42.0));
        }
        else
            xi = sin(w)/w;


        double cmu = cos(-u);
        double smu = sin(-u);
        double c2u = cmu*cmu-smu*smu;
        double s2u = -2.0*cmu*smu;
        double cw  = cos(w);
        double cwcmu = cw*cmu;
        double cwsmu = cw*smu;
        double xicmu = xi*cmu;
        double xismu = xi*smu;

        double h0r, h0i, h1r, h1i, h2r, h2i;
            h0r = (u2-w2)*c2u + 8.0*u2*cwcmu - 2.0*u*(3.0*u2+w2)*xismu;
            h0i = (u2-w2)*s2u + 8.0*u2*cwsmu + 2.0*u*(3.0*u2+w2)*xicmu;

            h1r = 2*u*c2u - 2*u*cwcmu - (3.0*u2-w2)*xismu;
            h1i = 2*u*s2u - 2*u*cwsmu + (3.0*u2-w2)*xicmu;

            h2r = c2u - cwcmu + 3*u*xismu;
            h2i = s2u - cwsmu - 3*u*xicmu;

            if (c0negative){
                h0i *= -1;
                h1r *= -1;
                h2i *= -1;
            }
            double norm = 1.0/(9.0*u2-w2);
        h0r*=norm;
        h0i*=norm;
        h1r*=norm;
        h1i*=norm;
        h2r*=norm;
        h2i*=norm;



        SU3 result;
    for(int i = 0; i < 18; i+=2){

    //F1
    result.mat[i] = h1r*Q.mat[i]-h1i*Q.mat[i+1];
    //F2
    result.mat[i] += h2r*Q2.mat[i]-h2i*Q2.mat[i+1];

    if(i%8 == 0)
        //F0
        result.mat[i] += h0r;

    //F1
    result.mat[i+1] = h1r*Q.mat[i+1]+h1i*Q.mat[i];
    //F2
    result.mat[i+1] += h2r*Q2.mat[i+1]+h2i*Q2.mat[i];
    if(i%8 == 0)
        //F0
        result.mat[i+1] += h0i;

    }
    return result;
}

SU3 exp(SU3&& Q){
    Q2 = Q*Q;
        Q3 = Q2*Q;

        double c0 = Q3.realTrace() / 3.0;
        double c1 = Q2.realTrace() / 2.0;

        bool c0negative = false;
        if(c0 < 0){
            c0negative = true;
            c0 = -c0;
        }

        double c1sqrt = sqrt(c1);
        double oneoversqrt3 = 1.0/sqrt(3);
        double c0max = c1*c1sqrt*oneoversqrt3*2.0 / 3.0;
        double theta = acos(c0/c0max);
        double w = c1sqrt*sin(theta/3.0);
        double u = c1sqrt*oneoversqrt3*cos(theta/3.0);
        double xi;
        double u2 = u*u;
        double w2 = w*w;

        if(fabs(w) < 0.05){
            xi = 1.0 - w2/6.0 * (1 - w2/20.0 * (1 - w2/42.0));
        }
        else
            xi = sin(w)/w;


        double cmu = cos(-u);
        double smu = sin(-u);
        double c2u = cmu*cmu-smu*smu;
        double s2u = -2.0*cmu*smu;
        double cw  = cos(w);
        double cwcmu = cw*cmu;
        double cwsmu = cw*smu;
        double xicmu = xi*cmu;
        double xismu = xi*smu;

        double h0r, h0i, h1r, h1i, h2r, h2i;
            h0r = (u2-w2)*c2u + 8.0*u2*cwcmu - 2.0*u*(3.0*u2+w2)*xismu;
            h0i = (u2-w2)*s2u + 8.0*u2*cwsmu + 2.0*u*(3.0*u2+w2)*xicmu;

            h1r = 2*u*c2u - 2*u*cwcmu - (3.0*u2-w2)*xismu;
            h1i = 2*u*s2u - 2*u*cwsmu + (3.0*u2-w2)*xicmu;

            h2r = c2u - cwcmu + 3*u*xismu;
            h2i = s2u - cwsmu - 3*u*xicmu;

            if (c0negative){
                h0i *= -1;
                h1r *= -1;
                h2i *= -1;
            }
            double norm = 1.0/(9.0*u2-w2);
        h0r*=norm;
        h0i*=norm;
        h1r*=norm;
        h1i*=norm;
        h2r*=norm;
        h2i*=norm;



        SU3 result;
    for(int i = 0; i < 18; i+=2){

    //F1
    result.mat[i] = h1r*Q.mat[i]-h1i*Q.mat[i+1];
    //F2
    result.mat[i] += h2r*Q2.mat[i]-h2i*Q2.mat[i+1];

    if(i%8 == 0)
        //F0
        result.mat[i] += h0r;

    //F1
    result.mat[i+1] = h1r*Q.mat[i+1]+h1i*Q.mat[i];
    //F2
    result.mat[i+1] += h2r*Q2.mat[i+1]+h2i*Q2.mat[i];
    if(i%8 == 0)
        //F0
        result.mat[i+1] += h0i;

    }
    return result;
}
