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




