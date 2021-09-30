//
// Created by li on 2021/9/9.
//

#ifndef BUTTERFILTER_SIGNAL_FILTER_H
#define BUTTERFILTER_SIGNAL_FILTER_H

#include <iostream>
#include <eigen3/Eigen/Dense>
#include <vector>

using namespace std;
using namespace Eigen;

class Signal_Filter {

    std::vector<double*> VEC_A_B;
    int len_a_b,pad_len,pad_x_len;
    double *zi;

    double* poly(const std::complex<double> data[],int N,double k);
    std::vector<complex<double>> convolve(std::vector<complex<double>> &data,const std::complex<double> kernel[]);

    // 开始进行滤波
    void validate_pad(float *x, int N,int pad_x_len,double *pad_x);
    void lfilter_zi(const std::vector<double *> vec_a_b);
    void lfilter(const std::vector<double *> &vec_a_b,double *x,double coeff,double *y);
public:
    void init(int order, int N,float low,float high,float fps);
    void cal_a_b(int ORDER, double L, double H, double FPS);
    void filtfilt(float* data,int N,float *y);
};


#endif //BUTTERFILTER_SIGNAL_FILTER_H
