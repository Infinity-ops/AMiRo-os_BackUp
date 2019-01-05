#include "userthread.hpp"
#include <chprintf.h>
#include <ctime>
#include <hal.h>
#include "global.hpp"


using namespace amiro;


extern Global global;

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
    i2sStart(&I2SD2, &global.i2scfg);
    i2sStartExchange(&I2SD2);
    for(int i=0; i< 72; i++)
      {
         this->sleep(MS2ST(1000));
         chprintf((BaseSequentialStream*) &global.sercanmux1,"%d s", i);
      }

      i2sStopExchange(&I2SD2);


//  chprintf((BaseSequentialStream*) &global.sercanmux1, "Init/n");
//  i2sInit();
//  this->sleep(MS2ST(1000));
//  i2sStart(&I2SD2, &global.i2scfg);
//  this->sleep(MS2ST(1000));
//  chprintf((BaseSequentialStream*) &global.sercanmux1, "Main/n");
//  i2sStartExchange(&I2SD2);
//  this->sleep(MS2ST(1000));
//  for(int i=0; i<25 ; i++)
//  {
//     this->sleep(MS2ST(1000));
//     chprintf((BaseSequentialStream*) &global.sercanmux1,"%d s", i);
//  }

//  i2sStopExchange(&I2SD2);
//  this->sleep(MS2ST(1000));
//  i2sStop(&I2SD2);
//  this->sleep(MS2ST(1000));
//  chprintf((BaseSequentialStream*) &global.sercanmux1,"starts exchange \n");

  uint16_t samples[I2S_BUF_SIZE];
//  float meanval = 0;
//  double PI = 3.14159265358;
  while (!this->shouldTerminate())
  {
      // int index = 0;
      for(int i=1, j = 0; i < I2S_BUF_SIZE; i+=2, j++)
      {
          // if ((global.i2s_rx_buf[i] != 0) && (global.i2s_rx_buf[i] != -1) )
//          {
//            samples[i] = (global.i2s_rx_buf[i] & 0x0000FFFF)<< 16;
//            samples[i] = ((global.i2s_rx_buf[i]<<16) | (global.i2s_rx_buf[i]>>16));
//            samples[i]  >>= 14;
              samples[i] = (global.i2s_rx_buf[i] & 0xFFFF);
              chprintf((BaseSequentialStream*) &global.sercanmux1, " %d, %d,%08X\n", j,samples[i], samples[i] );
//          }
        // index++
//          meanval += samples[j];
      }
//      meanval /= I2S_BUF_SIZE;
//      chprintf((BaseSequentialStream*) &global.sercanmux1,"%f s",meanval );

      //DFT

      for(int k=0; k < I2S_BUF_SIZE; k++){
          double sumReal= 0;
          double sumImg= 0;
          double angle;
          for(int n=0; n < I2S_BUF_SIZE; n++){
              angle = (-1) * 2 * PI * n * k /I2S_BUF_SIZE;
              sumReal = samples[n]* cos(angle);
              sumImg = samples[n]* sin(angle);
          }
      }

//      //FFT


  }

  return RDY_OK;
}

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




