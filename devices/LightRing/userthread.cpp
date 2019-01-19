#include "userthread.hpp"
#include <chprintf.h>
#include <ctime>
#include <hal.h>
#include <complex>
#include <math.h>
#include "global.hpp"


using namespace amiro;


extern Global global;

typedef std::complex<double> Complex;

UserThread::UserThread() :
    chibios_rt::BaseStaticThread<USER_THREAD_STACK_SIZE>()
{
}

UserThread::~UserThread()
{
}



msg_t
UserThread::main()
{
    chprintf((BaseSequentialStream*) &global.sercanmux1, "Init/n");
    i2sInit();



    uint16_t samples[I2S_BUF_SIZE];
    //  float meanval = 0;
    double PI = 3.141592653589793238460;

    int cycleNumber = 0;

    i2sStart(&I2SD2, &global.i2scfg);

    while (!this->shouldTerminate())
    {

        // input start
        chprintf((BaseSequentialStream*)&global.sercanmux1,"%d : Input started\n", cycleNumber++);

        this->sleep(MS2ST(1000));
        i2sStartExchange(&I2SD2);

        for(int i=0; i < 16; i++)
        {
            this->sleep(MS2ST(1000));
            chprintf((BaseSequentialStream*) &global.sercanmux1,"s:%d", i);
        }

        i2sStopExchange(&I2SD2);
        this->sleep(MS2ST(1000));

        chprintf((BaseSequentialStream*)&global.sercanmux1,"Input finished\n");
        // input end


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
