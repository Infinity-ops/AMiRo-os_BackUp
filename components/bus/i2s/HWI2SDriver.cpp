#include <hal.h>
#include <amiro/bus/i2s/HWI2SDriver.hpp>

using namespace chibios_rt;

namespace amiro {

HWI2SDriver::
HWI2SDriver(::I2SDriver *driver, const ::I2SConfig *config) :
  driver(driver),
  config(config) {

}

HWI2SDriver::
~HWI2SDriver() {

}

void
HWI2SDriver::
start_exchange(size_t length) {
  i2sStart(this->driver, this->config);
  i2sStartExchange(this->driver); 
}

void
HWI2SDriver::
stop_exchange(size_t length) {
  i2sStopExchange(this->driver); 
  i2sStop(this->driver, this->config);
}

} /* amiro */
