#include "public.sdk/source/vst/vstaudioeffect.h"
#include "public.sdk/source/vst/vstparameters.h"
#include "pluginterfaces/vst/vsttypes.h"
#include <cmath>
#include <vector>

// Define parameter IDs
enum {
    kGainParam = 0,
    kFrequencyParam,
    kQParam,
    kNumParams
};

// Define parameter ranges
static Steinberg::Vst::ParamValue kGainParamMin = -24.0;
static Steinberg::Vst::ParamValue kGainParamMax = 24.0;
static Steinberg::Vst::ParamValue kFrequencyParamMin = 20.0;
static Steinberg::Vst::ParamValue kFrequencyParamMax = 20000.0;
static Steinberg::Vst::ParamValue kQParamMin = 0.1;
static Steinberg::Vst::ParamValue kQParamMax = 10.0;

// Define the Filter class for the equalizer filters
class Filter {
public:
    Filter() : frequency_(0.0f), q_(1.0f), gain_(0.0f) {}
    
    // Setters for filter parameters
    void setFrequency(float frequency) { frequency_ = frequency; }
    void setQ(float q) { q_ = q; }
    void setGain(float gain) { gain_ = gain; }
    
    // Getters for filter parameters
    float getFrequency() const { return frequency_; }
    float getQ() const { return q_; }
    float getGain() const { return gain_; }
    
    // Process the audio signal with the filter
    void process(float* audio, int numSamples, float frequency, float q) {
        // Calculate filter coefficients
        float omega = 2.0f * M_PI * frequency / 44100.0f;
        float alpha = sin(omega) / (2.0f * q);
        float a0 = 1.0f + alpha;
        float a1 = -2.0f * cos(omega);
        float a2 = 1.0f - alpha;
        float b0 = (1.0f + cos(omega)) / 2.0f;
        float b1 = -(1.0f + cos(omega));
        float b2 = (1.0f + cos(omega)) / 2.0f;
        
        // Apply filter to audio signal
        float x1 = 0.0f, x2 = 0.0f, y1 = 0.0f, y2 = 0.0f;
        for (int i = 0; i < numSamples; i++) {
            float x = audio[i];
            float y = b0 / a0 * x + b1 / a0 * x1 + b2 / a0 * x2 - a1 / a0 * y1 - a2 / a0 * y2;
            audio[i] = y;
            x2 = x1;
            x1 = x;
            y2 = y1;
            y1 = y;
        }
    }
    
private:
    float frequency_; // Center frequency of the filter
    float q_; // Q factor of the filter
    float gain_; // Gain of the filter in decibels
};
class Filter {
public:
    Filter(float frequency, float q, float gain) 
        : frequency_(frequency), q_(q), gain_(gain) {
        updateCoefficients();
    }
    
    void process(float* audio, int numSamples) {
        for (int i = 0; i < numSamples; i++) {
            float input = audio[i];
            float output = b0_*input + b1_*x1_ + b2_*x2_ - a1_*y1_ - a2_*y2_;
            x2_ = x1_;
            x1_ = input;
            y2_ = y1_;
            y1_ = output;
            audio[i] = output;
        }
    }
    
    void setFrequency(float frequency) {
        frequency_ = frequency;
        updateCoefficients();
    }
    
    void setQ(float q) {
        q_ = q;
        updateCoefficients();
    }
    
    void setGain(float gain) {
        gain_ = gain;
        updateCoefficients();
    }
    
private:
    float frequency_; // Center frequency of the filter
    float q_; // Q factor of the filter
    float gain_; // Gain of the filter in decibels
    
    // Filter coefficients
    float b0_, b1_, b2_, a1_, a2_;
    float x1_ = 0, x2_ = 0, y1_ = 0, y2_ = 0;
    
    void updateCoefficients() {
        float w0 = 2 * M_PI * frequency_ / kSampleRate;
        float alpha = sin(w0) / (2 * q_);
        float A = pow(10, gain_ / 40);
        float beta = sqrt(A) + sqrt(A + 1);
        
        b0_ = (1 + alpha * A);
        b1_ = (-2 * cos(w0));
        b2_ = (1 - alpha * A);
        a1_ = (-2 * alpha * cos(w0)) / beta;
        a2_ = (alpha * alpha - sqrt(A)) / beta;
    }
};
class Processing_flow : public Steinberg::Vst::AudioEffect {
public:
    Processing_flow() {
        // Create a default filter with a center frequency of 1000 Hz, a Q factor of 1 and a gain of 0 dB
        filter_ = std::make_unique<Filter>(1000, 1, 0);
    }
    
    void process(Steinberg::Vst::ProcessData& data) override {
        // Get the audio buffer
        Steinberg::Vst::AudioBuffer* audioBuffer = data.inputs[0];
        float* audio = (float*)audioBuffer->channelBuffers32[0];
        int numSamples = data.numSamples;
        
        // Apply the filter to the audio signal
        filter_->process(audio, numSamples);
    }
    
    // Override the setParameter method to update the filter parameters
    Steinberg::tresult PLUGIN_API setParameter(Steinberg::Vst::ParamID id, Steinberg::Vst::ParamValue value, int32 flags) override {
        switch (id) {
            case kFrequencyParam:
                filter_->setFrequency(value);
                break;
            case kQParam:
                filter_->setQ(value);
                break;
            case kGainParam:
                filter_->setGain(value);
                break;
            default:
                break;
        }
        
        return Steinberg::kResultOk;
    }
    
private:
    std::unique_ptr<Filter> filter_;
};