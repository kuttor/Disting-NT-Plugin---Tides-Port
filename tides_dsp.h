// Tides 2 DSP Engine
// Ported from Mutable Instruments Tides 2
// Copyright 2017 Emilie Gillet - MIT License
// Adapted for Disting NT by Claude

#ifndef TIDES_DSP_H_
#define TIDES_DSP_H_

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>

#include "tides_resources.h"

namespace tides {

// ============================================================================
// DSP Helpers (from stmlib)
// ============================================================================

#define CONSTRAIN(var, min, max) \
  if (var < (min)) var = (min); \
  if (var > (max)) var = (max);

#define ONE_POLE(out, in, coefficient) \
  out += (coefficient) * ((in) - out)

#define MAKE_INTEGRAL_FRACTIONAL(x) \
  int32_t x ## _integral = static_cast<int32_t>(x); \
  float x ## _fractional = x - static_cast<float>(x ## _integral);

inline float Interpolate(const float* table, float index, float size) {
    index *= size;
    MAKE_INTEGRAL_FRACTIONAL(index)
    float a = table[index_integral];
    float b = table[index_integral + 1];
    return a + (b - a) * index_fractional;
}

inline float InterpolateWrap(const float* table, float index, float size) {
    index -= static_cast<float>(static_cast<int32_t>(index));
    index *= size;
    MAKE_INTEGRAL_FRACTIONAL(index)
    float a = table[index_integral];
    float b = table[index_integral + 1];
    return a + (b - a) * index_fractional;
}

// PolyBLEP functions for anti-aliasing
inline float ThisBlepSample(float t) {
    return 0.5f * t * t;
}

inline float NextBlepSample(float t) {
    t = 1.0f - t;
    return -0.5f * t * t;
}

inline float ThisIntegratedBlepSample(float t) {
    const float t1 = 0.5f * t;
    const float t2 = t1 * t1;
    const float t4 = t2 * t2;
    return 0.1875f - t1 + 1.5f * t2 - t2 * t1 - t4 + t4 * t1;
}

inline float NextIntegratedBlepSample(float t) {
    t = 1.0f - t;
    const float t1 = 0.5f * t;
    const float t2 = t1 * t1;
    const float t4 = t2 * t2;
    return 0.1875f - t1 + 1.5f * t2 - t2 * t1 - t4 + t4 * t1;
}

// ============================================================================
// Gate Flags
// ============================================================================

enum GateFlags {
    GATE_FLAG_LOW = 0,
    GATE_FLAG_HIGH = 1,
    GATE_FLAG_RISING = 2,
    GATE_FLAG_FALLING = 4
};

// ============================================================================
// Parameter Interpolator
// ============================================================================

class ParameterInterpolator {
public:
    ParameterInterpolator(float* state, float new_value, size_t size) {
        state_ = state;
        value_ = *state;
        increment_ = (new_value - *state) / static_cast<float>(size);
    }
    
    ~ParameterInterpolator() {
        *state_ = value_;
    }
    
    inline float Next() {
        value_ += increment_;
        return value_;
    }

private:
    float* state_;
    float value_;
    float increment_;
};

// ============================================================================
// Hysteresis Quantizer
// ============================================================================

class HysteresisQuantizer2 {
public:
    void Init(int num_steps, float hysteresis, bool symmetric) {
        num_steps_ = num_steps;
        hysteresis_ = hysteresis;
        symmetric_ = symmetric;
        quantized_value_ = 0;
    }
    
    int Process(float value) {
        value *= static_cast<float>(num_steps_);
        
        float hysteresis = quantized_value_ < value ? -hysteresis_ : hysteresis_;
        int q = static_cast<int>(value + hysteresis + 0.5f);
        CONSTRAIN(q, 0, num_steps_ - 1);
        quantized_value_ = q;
        return q;
    }
    
    template<typename T>
    const T& Lookup(const T* table, float value) {
        return table[Process(value)];
    }

private:
    int num_steps_;
    float hysteresis_;
    bool symmetric_;
    int quantized_value_;
};

// ============================================================================
// Enums
// ============================================================================

enum RampMode {
    RAMP_MODE_AD,
    RAMP_MODE_LOOPING,
    RAMP_MODE_AR,
    RAMP_MODE_LAST
};

enum OutputMode {
    OUTPUT_MODE_GATES,
    OUTPUT_MODE_AMPLITUDE,
    OUTPUT_MODE_SLOPE_PHASE,
    OUTPUT_MODE_FREQUENCY,
    OUTPUT_MODE_LAST,
};

enum Range {
    RANGE_CONTROL,
    RANGE_AUDIO,
    RANGE_LAST
};

// ============================================================================
// Ratio (for frequency mode polyrhythms)
// ============================================================================

struct Ratio {
    float ratio;
    int q;
};

// ============================================================================
// Ramp Generator
// ============================================================================

template<size_t num_channels = 4>
class RampGenerator {
public:
    RampGenerator() { }
    ~RampGenerator() { }
    
    inline float phase(size_t index) const { return phase_[index]; }
    inline float frequency(size_t index) const { return frequency_[index]; }
    
    inline void Init() {
        master_phase_ = 0.0f;
        std::fill(&phase_[0], &phase_[num_channels], 0.0f);
        std::fill(&frequency_[0], &frequency_[num_channels], 0.0f);
        std::fill(&wrap_counter_[0], &wrap_counter_[num_channels], 0);
        
        Ratio r;
        r.ratio = 1.0f;
        r.q = 1;
        std::fill(&ratio_[0], &ratio_[num_channels], r);
        next_ratio_ = &ratio_[0];
    }
    
    inline void set_next_ratio(const Ratio* next_ratio) {
        next_ratio_ = next_ratio;
    }
    
    template<RampMode ramp_mode, OutputMode output_mode, Range range, bool use_ramp>
    inline void Step(const float f0, const float* pw, GateFlags gate_flags, float ramp) {
        const size_t n = output_mode == OUTPUT_MODE_FREQUENCY ||
            (output_mode == OUTPUT_MODE_SLOPE_PHASE && ramp_mode == RAMP_MODE_AR)
                ? num_channels : 1;

        if (ramp_mode == RAMP_MODE_AD) {
            if (gate_flags & GATE_FLAG_RISING) {
                std::fill(&phase_[0], &phase_[n], 0.0f);
            }
            
            for (size_t i = 0; i < n; ++i) {
                frequency_[i] = std::min(f0 * next_ratio_[i].ratio, 0.25f);
                if (use_ramp) {
                    phase_[i] = ramp * next_ratio_[i].ratio;
                } else {
                    phase_[i] += frequency_[i];
                }
                phase_[i] = std::min(phase_[i], 1.0f);
            }
        }
        
        if (ramp_mode == RAMP_MODE_AR) {
            if (output_mode == OUTPUT_MODE_SLOPE_PHASE) {
                std::fill(&frequency_[0], &frequency_[n], f0);
            } else {
                for (size_t i = 0; i < n; ++i) {
                    frequency_[i] = std::min(f0 * next_ratio_[i].ratio, 0.25f);
                }
            }
            
            const bool should_ramp_up = use_ramp ? ramp < 0.5f : gate_flags & GATE_FLAG_HIGH;
            
            float clip_at = should_ramp_up ? 0.5f : 1.0f;
            for (size_t i = 0; i < n; ++i) {
                if (phase_[i] < 0.5f && !should_ramp_up) {
                    phase_[i] = 0.5f;
                } else if (phase_[i] > 0.5f && should_ramp_up) {
                    phase_[i] = 0.0f;
                }
                float this_pw = output_mode == OUTPUT_MODE_FREQUENCY ? pw[0] : pw[i];
                float slope = phase_[i] < 0.5f
                    ? 0.5f / (1.0e-6f + this_pw)
                    : 0.5f / (1.0f + 1.0e-6f - this_pw);
                phase_[i] += frequency_[i] * slope;
                phase_[i] = std::min(phase_[i], clip_at);
            }
        }
        
        if (ramp_mode == RAMP_MODE_LOOPING) {
            if (range == RANGE_AUDIO && output_mode == OUTPUT_MODE_FREQUENCY) {
                bool reset = false;
                if (gate_flags & GATE_FLAG_RISING) {
                    std::fill(&phase_[0], &phase_[n], 0.0f);
                    reset = true;
                }
                for (size_t i = 0; i < n; ++i) {
                    frequency_[i] = std::min(f0 * next_ratio_[i].ratio, 0.25f);
                }
                if (!reset) {
                    for (size_t i = 0; i < n; ++i) {
                        phase_[i] += frequency_[i];
                        if (phase_[i] >= 1.0f) {
                            phase_[i] -= 1.0f;
                        }
                    }
                }
            } else {
                if (use_ramp) {
                    for (size_t i = 0; i < n; ++i) {
                        frequency_[i] = std::min(f0 * ratio_[i].ratio, 0.25f);
                    }
                    if (ramp < master_phase_) {
                        for (size_t i = 0; i < n; ++i) {
                            ++wrap_counter_[i];
                            if (wrap_counter_[i] >= ratio_[i].q) {
                                ratio_[i] = next_ratio_[i];
                                wrap_counter_[i] = 0;
                            }
                        }
                    }
                    master_phase_ = ramp;
                } else {
                    bool reset = false;
                    if (gate_flags & GATE_FLAG_RISING) {
                        master_phase_ = 0.0f;
                        std::copy(&next_ratio_[0], &next_ratio_[n], &ratio_[0]);
                        std::fill(&wrap_counter_[0], &wrap_counter_[n], 0);
                        reset = true;
                    }
                    for (size_t i = 0; i < n; ++i) {
                        frequency_[i] = std::min(f0 * ratio_[i].ratio, 0.25f);
                    }
                    if (!reset) {
                        master_phase_ += f0;
                    }
                    if (master_phase_ >= 1.0f) {
                        master_phase_ -= 1.0f;
                        for (size_t i = 0; i < n; ++i) {
                            ++wrap_counter_[i];
                            if (wrap_counter_[i] >= ratio_[i].q) {
                                ratio_[i] = next_ratio_[i];
                                wrap_counter_[i] = 0;
                            }
                        }
                    }
                }
                for (size_t i = 0; i < n; ++i) {
                    float mult_phase = master_phase_ + float(wrap_counter_[i]);
                    mult_phase *= ratio_[i].ratio;
                    MAKE_INTEGRAL_FRACTIONAL(mult_phase);
                    phase_[i] = mult_phase_fractional;
                }
            }
        }
    }

private:
    const Ratio* next_ratio_;
    float master_phase_;
    int wrap_counter_[num_channels];
    float phase_[num_channels];
    float frequency_[num_channels];
    Ratio ratio_[num_channels];
};

// ============================================================================
// Ramp Shaper
// ============================================================================

class RampShaper {
public:
    RampShaper() { }
    ~RampShaper() { }
    
    void Init() {
        next_sample_ = 0.0f;
        previous_phase_shift_ = 0.0f;
    }
    
    inline float BandLimitedPulse(float phase, float frequency, float pw) {
        CONSTRAIN(pw, frequency * 2.0f, 1.0f - 2.0f * frequency);
        
        float this_sample = next_sample_;
        float next_sample = 0.0f;
        
        float wrap_point = pw;
        if (phase < pw * 0.5f) {
            wrap_point = 0.0f;
        } else if (phase > 0.5f + pw * 0.5f) {
            wrap_point = 1.0f;
        }
        
        const float d = phase - wrap_point;
        
        if (d >= 0.0f && d < frequency) {
            const float t = d / frequency;
            float discontinuity = 1.0f;
            if (wrap_point != pw) {
                discontinuity = -discontinuity;
            }
            if (frequency < 0.0f) {
                discontinuity = -discontinuity;
            }
            this_sample += ThisBlepSample(t) * discontinuity;
            next_sample += NextBlepSample(t) * discontinuity;
        }
        
        next_sample += phase < pw ? 0.0f : 1.0f;
        next_sample_ = next_sample;
        
        return this_sample;
    }
    
    template<RampMode ramp_mode, Range range>
    inline float Slope(float phase, float phase_shift, float frequency, float pw) {
        if (ramp_mode == RAMP_MODE_AD) {
            return SkewedRamp(phase, 0.0f, frequency, pw);
        } else if (ramp_mode == RAMP_MODE_AR) {
            return phase;
        } else {
            if (range == RANGE_CONTROL) {
                return SkewedRamp(phase, phase_shift, frequency, pw);
            } else {
                return BandLimitedSlope(phase, phase_shift, frequency, pw);
            }
        }
    }
    
    template<RampMode ramp_mode, Range range>
    inline float EOA(float phase, float frequency, float pw) {
        if (ramp_mode == RAMP_MODE_LOOPING && range == RANGE_AUDIO) {
            return BandLimitedPulse(phase, frequency, pw);
        } else if (ramp_mode == RAMP_MODE_AR) {
            return phase >= 0.5f ? 1.0f : 0.0f;
        } else {
            return phase >= pw ? 1.0f : 0.0f;
        }
    }
    
    template<RampMode ramp_mode, Range range>
    inline float EOR(float phase, float frequency, float pw) {
        if (ramp_mode == RAMP_MODE_LOOPING) {
            const float eor_pw = std::min(0.5f, 96.0f * frequency);
            if (range == RANGE_AUDIO) {
                return 1.0f - BandLimitedPulse(phase, frequency, eor_pw);
            } else {
                return phase < eor_pw ? 1.0f : 0.0f;
            }
        } else {
            return phase >= 1.0f ? 1.0f : 0.0f;
        }
    }

private:
    inline float BandLimitedSlope(float phase, float phase_shift, float frequency, float pw) {
        if (phase_shift != 0.0f) {
            phase += phase_shift;
            frequency += phase_shift - previous_phase_shift_;
            previous_phase_shift_ = phase_shift;
            
            if (phase >= 1.0f) {
                phase -= 1.0f;
            } else if (phase < 0.0f) {
                phase += 1.0f;
            }
        }
        
        float abs_freq = fabsf(frequency);
        CONSTRAIN(pw, abs_freq * 2.0f, 1.0f - 2.0f * abs_freq);
        
        float this_sample = next_sample_;
        float next_sample = 0.0f;
        
        float wrap_point = pw;
        if (phase < pw * 0.5f) {
            wrap_point = 0.0f;
        } else if (phase > 0.5f + pw * 0.5f) {
            wrap_point = 1.0f;
        }
        
        const float slope_up = 1.0f / pw;
        const float slope_down = 1.0f / (1.0f - pw);
        const float d = phase - wrap_point;
        
        if (d >= 0.0f && d < frequency) {
            const float t = d / frequency;
            float discontinuity = -(slope_up + slope_down) * frequency;
            if (wrap_point != pw) {
                discontinuity = -discontinuity;
            }
            if (frequency < 0.0f) {
                discontinuity = -discontinuity;
            }
            this_sample += ThisIntegratedBlepSample(t) * discontinuity;
            next_sample += NextIntegratedBlepSample(t) * discontinuity;
        }
        
        next_sample += phase < pw
            ? phase * slope_up
            : 1.0f - (phase - pw) * slope_down;
        next_sample_ = next_sample;
        
        return this_sample;
    }
    
    inline float SkewedRamp(float phase, float phase_shift, float frequency, float pw) {
        if (phase_shift != 0.0f) {
            phase += phase_shift;
            frequency += phase_shift - previous_phase_shift_;
            previous_phase_shift_ = phase_shift;
            
            if (phase >= 1.0f) {
                phase -= 1.0f;
            } else if (phase < 0.0f) {
                phase += 1.0f;
            }
        }
        
        float abs_freq = fabsf(frequency);
        CONSTRAIN(pw, abs_freq * 2.0f, 1.0f - 2.0f * abs_freq);
        const float slope_up = 0.5f / pw;
        const float slope_down = 0.5f / (1.0f - pw);
        return phase < pw ? phase * slope_up : (phase - pw) * slope_down + 0.5f;
    }
    
    float next_sample_;
    float previous_phase_shift_;
};

// ============================================================================
// Ramp Waveshaper
// ============================================================================

class RampWaveshaper {
public:
    RampWaveshaper() { }
    ~RampWaveshaper() { }
    
    void Init() {
        previous_input_ = 0.0f;
        previous_output_ = 0.0f;
        breakpoint_ = 0.0f;
    }
    
    template<RampMode ramp_mode>
    inline float Shape(float input, const int16_t* shape, float shape_fractional) {
        float ws_index = 1024.0f * input;
        MAKE_INTEGRAL_FRACTIONAL(ws_index)
        ws_index_integral &= 1023;
        float x0 = static_cast<float>(shape[ws_index_integral]) / 32768.0f;
        float x1 = static_cast<float>(shape[ws_index_integral + 1]) / 32768.0f;
        float y0 = static_cast<float>(shape[ws_index_integral + 1025]) / 32768.0f;
        float y1 = static_cast<float>(shape[ws_index_integral + 1026]) / 32768.0f;
        float x = x0 + (x1 - x0) * ws_index_fractional;
        float y = y0 + (y1 - y0) * ws_index_fractional;
        float output = x + (y - x) * shape_fractional;
        
        if (ramp_mode != RAMP_MODE_AR) {
            return output;
        } else {
            if (previous_input_ <= 0.5f && input > 0.5f) {
                breakpoint_ = previous_output_;
            } else if (previous_input_ > 0.5f && input < 0.5f) {
                breakpoint_ = previous_output_;
            } else if (input == 1.0f) {
                breakpoint_ = 1.0f;
            } else if (input == 0.5f) {
                breakpoint_ = 0.0f;
            }
            if (input <= 0.5f) {
                output = breakpoint_ + (1.0f - breakpoint_) * output;
            } else {
                output = breakpoint_ * output;
            }
            previous_input_ = input;
            previous_output_ = output;
            return output;
        }
    }

private:
    float previous_input_;
    float previous_output_;
    float breakpoint_;
};

// ============================================================================
// Filter
// ============================================================================

template<size_t num_channels>
class Filter {
public:
    Filter() { }
    ~Filter() { }
    
    void Init() {
        std::fill(&lp_1_[0], &lp_1_[num_channels], 0.0f);
        std::fill(&lp_2_[0], &lp_2_[num_channels], 0.0f);
    }
    
    template<size_t num_effective_channels>
    inline void Process(float* f, float* in_out, size_t size) {
        while (size--) {
            for (size_t i = 0; i < num_effective_channels; ++i) {
                ONE_POLE(lp_1_[i], *in_out, f[i]);
                ONE_POLE(lp_2_[i], lp_1_[i], f[i]);
                *in_out++ = lp_2_[i];
            }
            in_out += num_channels - num_effective_channels;
        }
    }

private:
    float lp_1_[num_channels];
    float lp_2_[num_channels];
};

// ============================================================================
// Ratio Tables
// ============================================================================

static Ratio audio_ratio_table[21][4] = {
    { { 1.0f, 1 }, { 0.5f, 2 }, { 0.25f, 4 }, { 0.125f, 8 } },
    { { 1.0f, 1 }, { 0.5f, 2 }, { 0.33333333f, 3 }, { 0.2f, 5 } },
    { { 1.0f, 1 }, { 0.5f, 2 }, { 0.33333333f, 3 }, { 0.25f, 4 } },
    { { 1.0f, 1 }, { 0.66666666f, 3 }, { 0.44444444f, 9 }, { 0.296296297f, 27 } },
    { { 1.0f, 1 }, { 0.66666666f, 3 }, { 0.5f, 2 }, { 0.33333333f, 3 } },
    { { 1.0f, 1 }, { 0.75f, 4 }, { 0.66666666f, 3 }, { 0.5f, 2 } },
    { { 1.0f, 1 }, { 0.790123456f, 81 }, { 0.66666666f, 3 }, { 0.5f, 2 } },
    { { 1.0f, 1 }, { 0.790123456f, 81 }, { 0.75f, 4 }, { 0.66666666f, 3 } },
    { { 1.0f, 1 }, { 0.88888888f, 9 }, { 0.790123456f, 81 }, { 0.66666666f, 3 } },
    { { 1.0f, 1 }, { 0.99090909091f, 109 }, { 0.987341772f, 79 }, { 0.9811320755f, 53 } },
    { { 1.0f, 1 }, { 1.0f, 1 }, { 1.0f, 1 }, { 1.0f, 1 } },
    { { 1.0f, 1 }, { 1.009174312f, 109 }, { 1.01265823f, 79 }, { 1.0188679245f, 53 } },
    { { 1.0f, 1 }, { 1.125f, 8 }, { 1.265625f, 64 }, { 1.5f, 2 } },
    { { 1.0f, 1 }, { 1.265625f, 64 }, { 1.3333333f, 3 }, { 1.5f, 2 } },
    { { 1.0f, 1 }, { 1.265625f, 64 }, { 1.5f, 2 }, { 2.0f, 1 } },
    { { 1.0f, 1 }, { 1.33333333f, 3 }, { 1.5f, 2 }, { 2.0f, 1 } },
    { { 1.0f, 1 }, { 1.5f, 2 }, { 2.0f, 1 }, { 3.0f, 1 } },
    { { 1.0f, 1 }, { 1.5f, 2 }, { 2.25f, 4 }, { 3.375f, 8 } },
    { { 1.0f, 1 }, { 2.0f, 1 }, { 3.0f, 1 }, { 4.0f, 1 } },
    { { 1.0f, 1 }, { 2.0f, 1 }, { 3.0f, 1 }, { 5.0f, 1 } },
    { { 1.0f, 1 }, { 2.0f, 1 }, { 4.0f, 1 }, { 8.0f, 1 } },
};

static Ratio control_ratio_table[21][4] = {
    { { 1.0f, 1 }, { 0.5f, 2 }, { 0.25f, 4 }, { 0.125f, 8 } },
    { { 1.0f, 1 }, { 0.5f, 2 }, { 0.33333333f, 3 }, { 0.2f, 5 } },
    { { 1.0f, 1 }, { 0.5f, 2 }, { 0.33333333f, 3 }, { 0.25f, 4 } },
    { { 1.0f, 1 }, { 0.66666666f, 3 }, { 0.5f, 2 }, { 0.25f, 4 } },
    { { 1.0f, 1 }, { 0.66666666f, 3 }, { 0.5f, 2 }, { 0.33333333f, 3 } },
    { { 1.0f, 1 }, { 0.75f, 4 }, { 0.66666666f, 3 }, { 0.5f, 2 } },
    { { 1.0f, 1 }, { 0.8f, 5 }, { 0.66666666f, 3 }, { 0.5f, 2 } },
    { { 1.0f, 1 }, { 0.8f, 5 }, { 0.75f, 3 }, { 0.5f, 2 } },
    { { 1.0f, 1 }, { 0.8f, 5 }, { 0.75f, 4 }, { 0.66666666f, 3 } },
    { { 1.0f, 1 }, { 0.909090909091f, 11 }, { 0.857142857143f, 7 }, { 0.8f, 5 } },
    { { 1.0f, 1 }, { 1.0f, 1 }, { 1.0f, 1 }, { 1.0f, 1 } },
    { { 1.0f, 1 }, { 1.09090909091f, 11 }, { 1.142857143f, 7 }, { 1.2f, 5 } },
    { { 1.0f, 1 }, { 1.25f, 4 }, { 1.33333333f, 3 }, { 1.5f, 2 } },
    { { 1.0f, 1 }, { 1.25f, 4 }, { 1.33333333f, 3 }, { 2.0f, 2 } },
    { { 1.0f, 1 }, { 1.25f, 4 }, { 1.5f, 3 }, { 2.0f, 2 } },
    { { 1.0f, 1 }, { 1.33333333f, 3 }, { 1.5f, 2 }, { 2.0f, 1 } },
    { { 1.0f, 1 }, { 1.5f, 2 }, { 2.0f, 1 }, { 3.0f, 1 } },
    { { 1.0f, 1 }, { 1.5f, 2 }, { 2.0f, 1 }, { 4.0f, 1 } },
    { { 1.0f, 1 }, { 2.0f, 1 }, { 3.0f, 1 }, { 4.0f, 1 } },
    { { 1.0f, 1 }, { 2.0f, 1 }, { 3.0f, 1 }, { 5.0f, 1 } },
    { { 1.0f, 1 }, { 2.0f, 1 }, { 4.0f, 1 }, { 8.0f, 1 } },
};

// ============================================================================
// Poly Slope Generator (main DSP engine)
// ============================================================================

class PolySlopeGenerator {
public:
    static constexpr size_t num_channels = 4;
    
    struct OutputSample {
        float channel[num_channels];
    };
    
    PolySlopeGenerator() { }
    ~PolySlopeGenerator() { }
    
    void Reset() {
        filter_.Init();
    }
    
    void Init() {
        frequency_ = 0.01f;
        pw_ = 0.0f;
        shift_ = 0.0f;
        shape_ = 0.0f;
        fold_ = 0.0f;
        
        ramp_generator_.Init();
        for (size_t i = 0; i < num_channels; ++i) {
            ramp_shaper_[i].Init();
            ramp_waveshaper_[i].Init();
        }
        filter_.Init();
        ratio_index_quantizer_.Init(21, 0.05f, false);
    }
    
    void Render(
            RampMode ramp_mode,
            OutputMode output_mode,
            Range range,
            float frequency,
            float pw,
            float shape,
            float smoothness,
            float shift,
            const GateFlags* gate_flags,
            const float* ramp,
            OutputSample* out,
            size_t size) {
        
        const float max_ratio = 1.0f;
        frequency = std::min(frequency, 0.25f * max_ratio);
        
        if (range == RANGE_CONTROL && pw < 0.5f) {
            pw = 0.5f + 0.6f * (pw - 0.5f) / (fabsf(pw - 0.5f) + 0.1f);
        }
        
        if (ramp && ramp_mode == RAMP_MODE_AR) {
            frequency *= 1.0f + 2.0f * fabsf(pw - 0.5f);
        }
        
        const float slope = 3.0f + fabsf(pw - 0.5f) * 5.0f;
        const float shape_amount = fabsf(shape - 0.5f) * 2.0f;
        const float shape_amount_attenuation = Tame(frequency, slope, 16.0f);
        shape = 0.5f + (shape - 0.5f) * shape_amount_attenuation;
        
        if (smoothness > 0.5f) {
            smoothness = 0.5f + (smoothness - 0.5f) * Tame(
                frequency,
                slope * (3.0f + shape_amount * shape_amount_attenuation * 5.0f),
                12.0f);
        }
        
        // Dispatch to the correct render function
        RenderDispatch(ramp_mode, output_mode, range,
            frequency, pw, shape, smoothness, shift, gate_flags, ramp, out, size);
        
        if (smoothness < 0.5f) {
            float ratio = smoothness * 2.0f;
            ratio *= ratio;
            ratio *= ratio;
            
            float f[4];
            size_t last_channel = output_mode == OUTPUT_MODE_GATES ? 1 : num_channels;
            for (size_t i = 0; i < last_channel; ++i) {
                size_t source = output_mode == OUTPUT_MODE_FREQUENCY ? i : 0;
                f[i] = ramp_generator_.frequency(source) * 0.5f;
                f[i] += (1.0f - f[i]) * ratio;
            }
            if (output_mode == OUTPUT_MODE_GATES) {
                filter_.Process<1>(f, &out[0].channel[0], size);
            } else {
                filter_.Process<num_channels>(f, &out[0].channel[0], size);
            }
        }
    }

private:
    void RenderDispatch(
            RampMode ramp_mode,
            OutputMode output_mode,
            Range range,
            float frequency, float pw, float shape, float smoothness, float shift,
            const GateFlags* gate_flags, const float* ramp,
            OutputSample* out, size_t size) {
        
        // Dispatch based on mode combination
        if (ramp_mode == RAMP_MODE_AD) {
            if (range == RANGE_CONTROL) {
                RenderInternal<RAMP_MODE_AD, RANGE_CONTROL>(
                    output_mode, frequency, pw, shape, smoothness, shift, gate_flags, ramp, out, size);
            } else {
                RenderInternal<RAMP_MODE_AD, RANGE_AUDIO>(
                    output_mode, frequency, pw, shape, smoothness, shift, gate_flags, ramp, out, size);
            }
        } else if (ramp_mode == RAMP_MODE_AR) {
            if (range == RANGE_CONTROL) {
                RenderInternal<RAMP_MODE_AR, RANGE_CONTROL>(
                    output_mode, frequency, pw, shape, smoothness, shift, gate_flags, ramp, out, size);
            } else {
                RenderInternal<RAMP_MODE_AR, RANGE_AUDIO>(
                    output_mode, frequency, pw, shape, smoothness, shift, gate_flags, ramp, out, size);
            }
        } else {
            if (range == RANGE_CONTROL) {
                RenderInternal<RAMP_MODE_LOOPING, RANGE_CONTROL>(
                    output_mode, frequency, pw, shape, smoothness, shift, gate_flags, ramp, out, size);
            } else {
                RenderInternal<RAMP_MODE_LOOPING, RANGE_AUDIO>(
                    output_mode, frequency, pw, shape, smoothness, shift, gate_flags, ramp, out, size);
            }
        }
    }
    
    template<RampMode ramp_mode, Range range>
    void RenderInternal(
            OutputMode output_mode,
            float frequency, float pw, float shape, float smoothness, float shift,
            const GateFlags* gate_flags, const float* ramp,
            OutputSample* out, size_t size) {
        
        const bool is_phasor = !(range == RANGE_AUDIO && ramp_mode == RAMP_MODE_LOOPING);
        
        ParameterInterpolator fm(&frequency_, frequency, size);
        ParameterInterpolator pwm(&pw_, pw, size);
        ParameterInterpolator shift_modulation(&shift_, 2.0f * shift - 1.0f, size);
        ParameterInterpolator shape_modulation(
            &shape_, is_phasor ? shape * 5.9999f + 5.0f : shape * 3.9999f, size);
        ParameterInterpolator fold_modulation(
            &fold_, std::max(2.0f * (smoothness - 0.5f), 0.0f), size);
        
        if (output_mode == OUTPUT_MODE_FREQUENCY) {
            const int ratio_index = ratio_index_quantizer_.Process(shift);
            if (range == RANGE_CONTROL) {
                ramp_generator_.set_next_ratio(control_ratio_table[ratio_index]);
            } else {
                ramp_generator_.set_next_ratio(audio_ratio_table[ratio_index]);
            }
        }
        
        for (size_t i = 0; i < size; ++i) {
            const float f0 = fm.Next();
            const float this_pw = pwm.Next();
            const float this_shift = shift_modulation.Next();
            const float step = this_shift * (1.0f / (num_channels - 1));
            const float partial_step = this_shift * (1.0f / num_channels);
            const float fold = fold_modulation.Next();
            
            float per_channel_pw[num_channels];
            const float pw_increment = (this_shift > 0.0f ? (1.0f - this_pw) : this_pw) * step;
            for (size_t j = 0; j < num_channels; ++j) {
                per_channel_pw[j] = this_pw + pw_increment * float(j);
            }
            
            // Increment ramps
            GateFlags gf = gate_flags ? gate_flags[i] : GATE_FLAG_LOW;
            float ramp_val = ramp ? ramp[i] : 0.0f;
            
            if (output_mode == OUTPUT_MODE_SLOPE_PHASE && ramp_mode == RAMP_MODE_AR) {
                if (ramp) {
                    ramp_generator_.Step<ramp_mode, OUTPUT_MODE_SLOPE_PHASE, range, true>(
                        f0, per_channel_pw, GATE_FLAG_LOW, ramp_val);
                } else {
                    ramp_generator_.Step<ramp_mode, OUTPUT_MODE_SLOPE_PHASE, range, false>(
                        f0, per_channel_pw, gf, 0.0f);
                }
            } else {
                if (ramp) {
                    if (output_mode == OUTPUT_MODE_GATES) {
                        ramp_generator_.Step<ramp_mode, OUTPUT_MODE_GATES, range, true>(
                            f0, &this_pw, GATE_FLAG_LOW, ramp_val);
                    } else if (output_mode == OUTPUT_MODE_AMPLITUDE) {
                        ramp_generator_.Step<ramp_mode, OUTPUT_MODE_AMPLITUDE, range, true>(
                            f0, &this_pw, GATE_FLAG_LOW, ramp_val);
                    } else if (output_mode == OUTPUT_MODE_SLOPE_PHASE) {
                        ramp_generator_.Step<ramp_mode, OUTPUT_MODE_SLOPE_PHASE, range, true>(
                            f0, &this_pw, GATE_FLAG_LOW, ramp_val);
                    } else {
                        ramp_generator_.Step<ramp_mode, OUTPUT_MODE_FREQUENCY, range, true>(
                            f0, &this_pw, GATE_FLAG_LOW, ramp_val);
                    }
                } else {
                    if (output_mode == OUTPUT_MODE_GATES) {
                        ramp_generator_.Step<ramp_mode, OUTPUT_MODE_GATES, range, false>(
                            f0, &this_pw, gf, 0.0f);
                    } else if (output_mode == OUTPUT_MODE_AMPLITUDE) {
                        ramp_generator_.Step<ramp_mode, OUTPUT_MODE_AMPLITUDE, range, false>(
                            f0, &this_pw, gf, 0.0f);
                    } else if (output_mode == OUTPUT_MODE_SLOPE_PHASE) {
                        ramp_generator_.Step<ramp_mode, OUTPUT_MODE_SLOPE_PHASE, range, false>(
                            f0, &this_pw, gf, 0.0f);
                    } else {
                        ramp_generator_.Step<ramp_mode, OUTPUT_MODE_FREQUENCY, range, false>(
                            f0, &this_pw, gf, 0.0f);
                    }
                }
            }
            
            // Compute shape
            const float shape_val = shape_modulation.Next();
            MAKE_INTEGRAL_FRACTIONAL(shape_val);
            const int16_t* shape_table = &lut_wavetable[shape_val_integral * 1025];
            
            if (output_mode == OUTPUT_MODE_GATES) {
                const float phase = ramp_generator_.phase(0);
                const float freq = ramp_generator_.frequency(0);
                const float raw = ramp_shaper_[0].Slope<ramp_mode, range>(phase, 0.0f, freq, this_pw);
                const float slope = ramp_waveshaper_[0].Shape<ramp_mode>(raw, shape_table, shape_val_fractional);
                
                out[i].channel[0] = Fold<ramp_mode>(slope, fold) * this_shift;
                out[i].channel[1] = Scale<ramp_mode>(is_phasor
                    ? ramp_waveshaper_[1].Shape<ramp_mode>(raw, &lut_wavetable[8200], 0.0f)
                    : raw);
                out[i].channel[2] = ramp_shaper_[2].EOA<ramp_mode, range>(phase, freq, this_pw) * 8.0f;
                out[i].channel[3] = ramp_shaper_[3].EOR<ramp_mode, range>(phase, freq, this_pw) * 8.0f;
            } else if (output_mode == OUTPUT_MODE_AMPLITUDE) {
                const float phase = ramp_generator_.phase(0);
                const float freq = ramp_generator_.frequency(0);
                const float raw = ramp_shaper_[0].Slope<ramp_mode, range>(phase, 0.0f, freq, this_pw);
                const float shaped = ramp_waveshaper_[0].Shape<ramp_mode>(raw, shape_table, shape_val_fractional);
                const float slope = Fold<ramp_mode>(shaped, fold) * (this_shift < 0.0f ? -1.0f : 1.0f);
                const float channel_index = fabsf(this_shift * 5.1f);
                for (size_t j = 0; j < num_channels; ++j) {
                    const float channel = static_cast<float>(j + 1);
                    const float gain = std::max(1.0f - fabsf(channel - channel_index), 0.0f);
                    const bool equal_pow = range == RANGE_AUDIO;
                    out[i].channel[j] = slope * gain * (equal_pow ? (2.0f - gain) : 1.0f);
                }
            } else if (output_mode == OUTPUT_MODE_SLOPE_PHASE) {
                float phase_shift = 0.0f;
                for (size_t j = 0; j < num_channels; ++j) {
                    size_t source = ramp_mode == RAMP_MODE_AR ? j : 0;
                    out[i].channel[j] = Fold<ramp_mode>(
                        ramp_waveshaper_[j].Shape<ramp_mode>(
                            ramp_shaper_[j].Slope<ramp_mode, range>(
                                ramp_generator_.phase(source),
                                phase_shift,
                                ramp_generator_.frequency(source),
                                ramp_mode == RAMP_MODE_AD ? per_channel_pw[j] : this_pw),
                            shape_table,
                            shape_val_fractional),
                        fold);
                    phase_shift -= range == RANGE_AUDIO ? step : partial_step;
                }
            } else if (output_mode == OUTPUT_MODE_FREQUENCY) {
                for (size_t j = 0; j < num_channels; ++j) {
                    out[i].channel[j] = Fold<ramp_mode>(
                        ramp_waveshaper_[j].Shape<ramp_mode>(
                            ramp_shaper_[j].Slope<ramp_mode, range>(
                                ramp_generator_.phase(j),
                                0.0f,
                                ramp_generator_.frequency(j),
                                this_pw),
                            shape_table,
                            shape_val_fractional),
                        fold);
                }
            }
        }
    }
    
    template<RampMode ramp_mode>
    inline float Fold(float unipolar, float fold_amount) {
        if (ramp_mode == RAMP_MODE_LOOPING) {
            float bipolar = 2.0f * unipolar - 1.0f;
            float folded = fold_amount > 0.0f ? Interpolate(
                lut_bipolar_fold,
                0.5f + bipolar * (0.03f + 0.46f * fold_amount),
                1024.0f) : 0.0f;
            return 5.0f * (bipolar + (folded - bipolar) * fold_amount);
        } else {
            float folded = fold_amount > 0.0f ? Interpolate(
                lut_unipolar_fold,
                unipolar * fold_amount,
                1024.0f) : 0.0f;
            return 8.0f * (unipolar + (folded - unipolar) * fold_amount);
        }
    }
    
    template<RampMode ramp_mode>
    inline float Scale(float unipolar) {
        if (ramp_mode == RAMP_MODE_LOOPING) {
            return 10.0f * unipolar - 5.0f;
        } else {
            return 8.0f * unipolar;
        }
    }
    
    inline float Tame(float f0, float harmonics, float order) {
        f0 *= harmonics;
        float max_f = 0.5f * (1.0f / order);
        float max_amount = 1.0f - (f0 - max_f) / (0.5f - max_f);
        CONSTRAIN(max_amount, 0.0f, 1.0f);
        return max_amount * max_amount * max_amount;
    }
    
    float frequency_;
    float pw_;
    float shift_;
    float shape_;
    float fold_;
    
    HysteresisQuantizer2 ratio_index_quantizer_;
    RampGenerator<num_channels> ramp_generator_;
    RampShaper ramp_shaper_[num_channels];
    RampWaveshaper ramp_waveshaper_[num_channels];
    Filter<num_channels> filter_;
};

}  // namespace tides

#endif  // TIDES_DSP_H_
