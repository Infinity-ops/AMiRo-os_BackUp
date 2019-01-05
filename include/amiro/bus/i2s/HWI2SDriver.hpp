#ifndef HWI2SDRIVER_HPP_
#define HWI2SDRIVER_HPP_

#include <ch.hpp>
#include <hal.h>

/**
 * @brief   AMiRo
 */
namespace amiro {


  /**
   * @brief   I²S HW Driver Abstraction
   */

class HWI2SDriver {
public:
  HWSI2SDriver(::I2SDriver *driver);
  virtual ~HWI2SDriver();


private:
  ::I2SDriver *driver;

};



} /* amiro */

#endif /* HWI2SDRIVER_HPP_ */


