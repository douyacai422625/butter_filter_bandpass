//
// Created by li on 2021/9/9.
//

#include "Signal_Filter.h"

#define PI 3.14159

void Signal_Filter::init(int order, int N) {
    len_a_b = 2 * order + 1;
    pad_len = 3 * len_a_b;
    pad_x_len = N + 2 * pad_len;

    zi = new double[len_a_b - 1];
}

void Signal_Filter::cal_a_b(int order, double low, double high, double fps, std::vector<double *> &VEC_A_B) {
    double Wn[2] = {double(2.0 * low / fps), double(2.0 * high / fps)};

    ///////////////////////////////////// 实现python 函数的 buttap
    std::complex<double> p[2 * order], z[2 * order];
    std::complex<double> *p_ptr = p, *z_ptr = z;
    int p_index = 0;

    ///////////////////////////////////// 实现python 函数的 lp2bp_zpk
    double bw = 4 * (tan(PI * Wn[1] / 2.0) - tan(PI * Wn[0] / 2.0));
    double wo = sqrt(16 * (tan(PI * Wn[1] / 2.0) * tan(PI * Wn[0] / 2.0)));
    double k_bp = pow(bw, order);

    ///////////////////////////////////// 实现python 函数的 bilinear_zpk
    double fs2 = 4.0;
    std::complex<double> prod_fs2_p, prod_fs2_z;

    for (int i = -order + 1; i < order; i += 2) {

        ///////////////////////////////////// 实现python 函数的 buttap
        double data = (PI * double(i)) / (2 * double(order));

        std::complex<double> p_lp = -complex<double>(cos(data), sin(data));

        ///////////////////////////////////// 实现python 函数的 lp2bp_zpk
        p_lp *= bw / 2;
        std::complex<double> pow_p_lp = pow(p_lp, 2);
        std::complex<double> public_data = sqrt(std::complex<double>(pow_p_lp.real() - pow(wo, 2), pow_p_lp.imag()));
        std::complex<double> p_bp1 = p_lp + public_data;
        std::complex<double> p_bp2 = p_lp - public_data;

        std::complex<double> fenzi_p_z1 = std::complex<double>(fs2 + p_bp1.real(), p_bp1.imag());
        std::complex<double> fenmu_p_z1 = std::complex<double>(fs2 - p_bp1.real(), -p_bp1.imag());
        std::complex<double> fenzi_p_z2 = std::complex<double>(fs2 + p_bp2.real(), p_bp2.imag());
        std::complex<double> fenmu_p_z2 = std::complex<double>(fs2 - p_bp2.real(), -p_bp2.imag());

        std::complex<double> p_z1 = fenzi_p_z1 / fenmu_p_z1;
        std::complex<double> p_z2 = fenzi_p_z2 / fenmu_p_z2;

        std::complex<double> prod_fenmu = fenmu_p_z1 * fenmu_p_z2;
        std::complex<double> prod_fenzi = std::complex<double>(fs2, 0);
        if (i == -order + 1) {
            prod_fs2_p = prod_fenmu;
            prod_fs2_z = prod_fenzi;
        } else {
            prod_fs2_p *= prod_fenmu;
            prod_fs2_z *= prod_fenzi;
        }

        p_ptr[p_index] = p_z1;
        p_ptr[p_index + order] = p_z2;
        z_ptr[p_index] = std::complex<double>(1, 0);
        z_ptr[p_index + order] = std::complex<double>(-1, 0);
        p_index++;
    }

    double k_z = k_bp * (prod_fs2_z / prod_fs2_p).real();

    double *b = poly(z, 2 * order, k_z);
    double *a = poly(p, 2 * order, 1);

    VEC_A_B.push_back(a);
    VEC_A_B.push_back(b);

    lfilter_zi(VEC_A_B);
}

double *Signal_Filter::poly(const std::complex<double> data[], int N, double k) {
    std::vector<complex<double>> a;
    a.resize(1);
    a[0] = std::complex<double>(1, 0);

    for (int i = 0; i < N; ++i) {
        std::complex<double> kernel[2] = {std::complex<double>(1, 0), -data[i]};
        a = convolve(a, kernel);
    }
    double *return_data = new double[a.size()];
    double *return_data_ptr = return_data;
    for (int i = 0; i < a.size(); ++i) {
        return_data_ptr[i] = a[i].real() * k;
    }
    return return_data;
}

std::vector<complex<double>>
Signal_Filter::convolve(std::vector<complex<double>> &data, const std::complex<double> *kernel) {
    int data_len = data.size() + 1;

    std::vector<complex<double>> return_data;
    return_data.resize(data_len);

    return_data[0] = data[0] * kernel[0];
    return_data[data_len - 1] = data[data_len - 2] * kernel[1];

    if (data_len > 2) {
        for (int i = 1; i < data_len - 1; ++i) {
            return_data[i] = data[i] * kernel[0] + data[i - 1] * kernel[1];
        }
    }
    return return_data;
}

////////////////////////// 开始对数据进行滤波/////////////////////////////////
void Signal_Filter::filtfilt(const std::vector<double *> vec_a_b, double *X, int N,double *y) {
    // 对数据行扩充
    double pad_x[pad_x_len];
    validate_pad(X, N, pad_x_len, pad_x);

    // Forward filter.
    double Forward_y[pad_x_len];
    lfilter(vec_a_b, pad_x, pad_x[0], Forward_y);

    // Backward filter.
    double y0 = Forward_y[pad_x_len - 1];
    reverse(Forward_y, Forward_y + pad_x_len);
    double Backward_y[pad_x_len];
    lfilter(vec_a_b, Forward_y, y0, Backward_y);
    reverse(Backward_y, Backward_y + pad_x_len);
    double *y_ptr = y;
    if (pad_len > 0) {
        for (int i = pad_len; i < pad_x_len - pad_len; ++i) {
            y_ptr[i - pad_len] = Backward_y[i];
        }
    }
}


void Signal_Filter::lfilter(const std::vector<double *> &vec_a_b, double *x, double coeff, double *y) {
    double *x_ptr = x, *y_ptr = y, *zi_ptr = zi;
    double *a_ptr = vec_a_b[0], *b_ptr = vec_a_b[1];
    for (int i = 0; i < pad_x_len; ++i) {
        double z = 0;
        if (i > 0){
            int inter = min(i+1,len_a_b);
            for (int j = 1; j < inter; ++j) {
                z += b_ptr[j] * x_ptr[i-j] - a_ptr[j] * y_ptr[i-j];
            }
            if (i < len_a_b-1){
                z += (zi_ptr[i] * coeff);
            }
        } else {
            z = zi_ptr[0] * coeff;
        }
        y_ptr[i] = b_ptr[0] * x_ptr[i] + z;
    }
}

void Signal_Filter::lfilter_zi(const std::vector<double *> vec_a_b) {
    double sum_B = 0.0, sum_IminusA = 1.0;
    double asum = 1.0;
    double csum = 0.0;
    for (int i = 1; i < len_a_b; ++i) {
        sum_B += vec_a_b[1][i] - vec_a_b[0][i] * vec_a_b[1][0];   // B = b[1:] - a[1:]*b[0]
        sum_IminusA += vec_a_b[0][i];
    }

    zi[0] = sum_B / sum_IminusA;
//    zi[0] = -0.0058862161552206685;
    for (int i = 1; i < len_a_b - 1; ++i) {
        asum += vec_a_b[0][i];
        csum += vec_a_b[1][i] - vec_a_b[0][i] * vec_a_b[1][0];
        zi[i] = asum * zi[0] - csum;
    }
}

void Signal_Filter::validate_pad(double *x, int N, int pad_x_len, double *pad_x) {
    double *x_ptr = x;
    double *pad_x_ptr = pad_x;

    for (int i = 0; i < pad_x_len; ++i) {
        if (i < pad_len) {
            pad_x_ptr[i] = 2 * x_ptr[0] - x_ptr[pad_len - i];
        } else if (i < pad_len + N) {
            pad_x_ptr[i] = x_ptr[i - pad_len];
        } else {
            pad_x_ptr[i] = 2 * x_ptr[N - 1] - x_ptr[2 * N - i + pad_len - 2];
        }
    }
}






















