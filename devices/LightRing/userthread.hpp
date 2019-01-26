#ifndef AMIRO_USERTHREAD_H_
#define AMIRO_USERTHREAD_H_

#include <ch.hpp>
#include <amiroosconf.h>
#include <vector>
#include <complex>

namespace amiro {

class UserThread : public chibios_rt::BaseStaticThread<USER_THREAD_STACK_SIZE>
{
public:
    explicit UserThread();
    virtual ~UserThread();

    virtual msg_t main();

private:
    int cycleNumber;
    void microphoneInput();
    void adjustData(std::vector<std::complex<float> > &outFftOutput);
    void sleepForSec(int inSeconds);
    std::vector<std::complex<float> > computeDft(const std::vector<std::complex<float> > &input, const int inFtRange);
    void manualDftIncomplete();
    void printFftResult(const std::vector<std::complex<float> > &inFftInput,
                        const std::vector<std::complex<float> > &inFftOutput);
    void lightOffAll();
    float ftThreshold(const std::vector<std::complex<float> > &inFftOutput);
};

} // end of namespace amiro

#endif // AMIRO_USERTHREAD_H_
