#include "userthread.hpp"
#include <chprintf.h>
#include <ctime>
#include <hal.h>
#include <complex>
#include <math.h>
#include <vector>
#include <iostream>
#include "global.hpp"


using namespace amiro;


extern Global global;

// For storing complex values of nth roots
// of unity we use complex<double>
typedef std::complex<float> cd;

UserThread::UserThread() :
    chibios_rt::BaseStaticThread<USER_THREAD_STACK_SIZE>()
{
}

UserThread::~UserThread()
{
}

// Recursive function of FFT
std::vector<cd> fft(std::vector<cd>& a)
{
    int n = a.size();

    // if input contains just one element
    if (n == 1)
        return std::vector<cd>(1, a[0]);

    // For storing n complex nth roots of unity
    std::vector<cd> w(n);
    for (int i = 0; i < n; i++) {
        double alpha = 2 * M_PI * i / n;
        w[i] = cd(cos(alpha), sin(alpha));
    }

    std::vector<cd> A0(n / 2), A1(n / 2);
    for (int i = 0; i < n / 2; i++) {

        // even indexed coefficients
        A0[i] = a[i * 2];

        // odd indexed coefficients
        A1[i] = a[i * 2 + 1];
    }

    // Recursive call for even indexed coefficients
    std::vector<cd> y0 = fft(A0);

    // Recursive call for odd indexed coefficients
    std::vector<cd> y1 = fft(A1);

    // for storing values of y0, y1, y2, ..., yn-1.
    std::vector<cd> y(n);

    for (int k = 0; k < n / 2; k++) {
        y[k] = y0[k] + w[k] * y1[k];
        y[k + n / 2] = y0[k] - w[k] * y1[k];
    }
    return y;
}


msg_t
UserThread::main()
{
    this->sleep(MS2ST(1000));
    this->sleep(MS2ST(1000));
    this->sleep(MS2ST(1000));
    this->sleep(MS2ST(1000));
    this->sleep(MS2ST(1000));
    this->sleep(MS2ST(1000));
    this->sleep(MS2ST(1000));
    this->sleep(MS2ST(1000));
    this->sleep(MS2ST(1000));
    this->sleep(MS2ST(1000));
    this->sleep(MS2ST(1000));
    this->sleep(MS2ST(1000));

    // FFT YAZ Debug
    std::vector<cd> a{1, 2, 3, 4};
    std::vector<cd> b = fft(a);

    std::vector<cd> fftInput;
    fftInput.resize(I2S_BUF_SIZE/2);
    std::vector<cd> fftResult;
    fftInput.resize(I2S_BUF_SIZE/2);

//    for (int i = 0; i < 4; i++)
//        std::cout << b[i] << std::endl;

    chprintf((BaseSequentialStream*) &global.sercanmux1, "Init/n");
    i2sInit();



    // uint16_t samples[I2S_BUF_SIZE];
    //  float meanval = 0;
    // double PI = 3.141592653589793238460;

    int cycleNumber = 0;

    i2sStart(&I2SD2, &global.i2scfg);

    while (!this->shouldTerminate())
    {

        // input start
        chprintf((BaseSequentialStream*)&global.sercanmux1,"%d : Input started\n", cycleNumber++);

        this->sleep(MS2ST(1000));
        i2sStartExchange(&I2SD2);

        int restTime = 6;
        for(int i=0; i < restTime; i++)
        {
            this->sleep(MS2ST(1000));
            chprintf((BaseSequentialStream*) &global.sercanmux1,"s:%d", i);
        }

        i2sStopExchange(&I2SD2);
        this->sleep(MS2ST(1000));

        chprintf((BaseSequentialStream*)&global.sercanmux1,"Input finished\n");
        // input end

        /*
        // int index = 0;
        chprintf((BaseSequentialStream*)&global.sercanmux1,"Printing Started\n");
        for(int i=1, k = 0; i < I2S_BUF_SIZE; i+=1, k++)
        {

            samples[i] = (global.i2s_rx_buf[i] & 0xFFFF);
//            samples[i] = global.i2s_rx_buf[i];


//            chprintf((BaseSequentialStream*) &global.sercanmux1, " %d, %d,%08X\n", k,samples[i], samples[i] );

//            double fft_val_real = 0.0;
//            double fft_val_imag = 0.0;

//            for(int j = 1, n = 0; j < I2S_BUF_SIZE; j += 2, n++){
//                Complex fft_val = std::polar(1.0, -2 * PI * k * n / (I2S_BUF_SIZE - 1)) * (samples[i] * 1.0);
//                fft_val_real += real(fft_val);
//                fft_val_imag += imag(fft_val);
//            }

//            double absolute = sqrt(pow(fft_val_real, 2.0) + pow(fft_val_imag, 2.0));
//            chprintf((BaseSequentialStream*)&global.sercanmux1,"%d,%d, %f,%f\n", k, samples[i],(32000.0 * k / 1024.0), absolute); //(k->value, raw(d1)->data, absolute or i2i2s_fft_buf -> complex)
            chprintf((BaseSequentialStream*)&global.sercanmux1,"%d,%d\n", i, samples[i]); //(k->value, raw(d1)->data, absolute or i2i2s_fft_buf -> complex)
        }
        chprintf((BaseSequentialStream*)&global.sercanmux1,"Printing Finished\n");
        */

        // FFT Call

        this->sleep(MS2ST(1000));
        this->sleep(MS2ST(1000));
        this->sleep(MS2ST(1000));


        float frequency = 0;
        for(int i = 1,k=0; i < I2S_BUF_SIZE; i = i+2,k++)
        {
            uint32_t raw = global.i2s_rx_buf[i];
            int16_t data = (raw & 0xFFFF);
            if(data == 0){
                data = (raw >> 16) & 0xFFFF;
            }

            // "data" is what you want aka the correct data
            // put this value to a vector of complex<float>
            fftInput.at(k) = static_cast<float>(data);
            chprintf((BaseSequentialStream*)&global.sercanmux1,"%d,%d\n", i, data);

        }

        chprintf((BaseSequentialStream*)&global.sercanmux1,"FFT Start\n");
        std::vector<cd> fftOutput = fft(fftInput);
        chprintf((BaseSequentialStream*)&global.sercanmux1,"FFT End\n");

        chprintf((BaseSequentialStream*)&global.sercanmux1,"Index,Data,Frequency,Real,Imaginary,Abs\n");
        int fftOutputSize = fftOutput.size();
        float dataComplexSize = static_cast<float>(fftInput.size());
        for(int i = 0; i < fftOutputSize; i++)
        {
            frequency = 32000.0 / dataComplexSize * static_cast<float>(i);

            chprintf((BaseSequentialStream*)&global.sercanmux1,"%d,%f,%f,%f,%f,%f\n", i, fftInput.at(i).real(), frequency, fftOutput.at(i).real(), fftOutput.at(i).imag(), std::abs(fftOutput.at(i)));
        }

        this->sleep(MS2ST(1000));
        this->sleep(MS2ST(1000));
        this->sleep(MS2ST(1000));

    }

    i2sStop(&I2SD2);

    return RDY_OK;
}

/*
void fourier_analysis_two_channel() {

    chprintf((BaseSequentialStream*)&global.sercanmux1,"k, d1, absolute, absolute_32, i2s_fft_buf[k]\n");
    double PI = 3.141592653589793238460;

    for (size_t i = 1, k = 0; i < I2S_BUF_SIZE; i += 2) {
        uint32_t raw1 = global.i2s_rx_buf[i] & 0x0000FFFF;
        int16_t d1 = raw1 & 0x0000FFFF;
        int32_t d_1_32 = ((raw1 & 0x0000FFFF) << 16) | ((raw1 & 0xFFFF0000) >> 16);

        double fft_val_real = 0.0;
        double fft_val_imag = 0.0;
        double fft_val_real_32 = 0.0;
        double fft_val_imag_32 = 0.0;

        for (size_t j = 1, n = 0; j < I2S_BUF_SIZE; j += 2) {

            uint32_t raw = global.i2s_rx_buf[j];
            int16_t d = raw & 0x0000FFFF;
            int32_t d_32 = ((raw & 0x0000FFFF) << 16) | ((raw & 0xFFFF0000) >> 16);

            Complex fft_val = std::polar(1.0, -2 * PI * k * n / (999)) * (d*1.0);
            Complex fft_val_32 = std::polar(1.0, -2 * PI * k * n / (999)) * (d_32*1.0);

            fft_val_real += real(fft_val);
            fft_val_imag += imag(fft_val);
            fft_val_real_32 += real(fft_val_32);
            fft_val_imag_32 += imag(fft_val_32);

            n++;
        }


        double absolute = sqrt(pow(fft_val_real, 2.0) + pow(fft_val_imag, 2.0));
        double absolute_32 = sqrt(pow(fft_val_real_32, 2.0) + pow(fft_val_imag_32, 2.0));

        i2s_fft_buf[k] = absolute;
        k++;

        chprintf((BaseSequentialStream*)&global.sercanmux1,"%d,%d,%f,%f,%f\n", k, d1, absolute, absolute_32, i2s_fft_buf[k]); //(k->value, raw(d1)->data, absolute or i2i2s_fft_buf -> complex)
    }
}
*/
