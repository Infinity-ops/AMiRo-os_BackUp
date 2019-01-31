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
    void printDftResult(const std::vector<std::complex<float> > &inFftInput,
                        const std::vector<std::complex<float> > &inFftOutput);
    void lightOffAll();
    float ftThreshold(const std::vector<std::complex<float> > &inFftOutput);
    float ftThreshold2(const std::vector<std::complex<float> > &inFftOutput);

    void lightOnlyTheHightestFreq(const std::vector<std::complex<float> > &inFftOutput, int inFfRange);
    void lightTillHighestFrequency(const std::vector<std::complex<float> > &inFftOutput, int inFfRange);

    float maxFtValue;
    int maxFtIndex;

    void mockDataFftInput(std::vector<std::complex<float>> &outFftInput);
    void ftSpecifications(int &outFtRange, int &outAcBufferSize);

    void motorControl(std::vector<std::complex<float>> &outFftInput);

    void middleChunkOfVector(std::vector<std::complex<float>> &inFftInput,
                             int &outStartPoint, int &outEndPoint, int inWindowSize);
};

} // end of namespace amiro

#endif // AMIRO_USERTHREAD_H_
