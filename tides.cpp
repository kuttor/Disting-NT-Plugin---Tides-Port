// Tides 2 for Disting NT
// Authentic port of Mutable Instruments Tides 2
// Original DSP: Copyright 2017 Emilie Gillet - MIT License
// Disting NT adaptation: Claude/Anthropic

#include <math.h>
#include <new>
#include <cstring>
#include <algorithm>
#include <distingnt/api.h>

#include "tides_dsp.h"

// ============================================================================
// Algorithm Structure (inherits from _NT_algorithm)
// ============================================================================

struct _TidesAlgorithm : public _NT_algorithm
{
    _TidesAlgorithm() {}
    ~_TidesAlgorithm() {}
    
    tides::PolySlopeGenerator poly_slope_generator;
    tides::PolySlopeGenerator::OutputSample output_buffer[64];
    tides::GateFlags gate_buffer[64];
    float ramp_buffer[64];
    
    // Cached parameters
    float sample_rate;
    float inv_sample_rate;
    
    // Gate detection state
    bool prev_trig_high;
    bool prev_clock_high;
    float clock_phase;
    float clock_period;
    int clock_counter;
    bool clock_synced;
    
    // Previous output mode for reset detection
    int prev_output_mode;
    
    // Cached derived values
    float cached_frequency;
    float cached_shape;
    float cached_slope;
    float cached_smooth;
    float cached_shift;
};

// ============================================================================
// Parameter Definitions
// ============================================================================

enum {
    // Inputs
    kParam_TrigInput,
    kParam_ClockInput,
    kParam_VOctInput,
    kParam_FMInput,
    kParam_ShapeInput,
    kParam_SlopeInput,
    kParam_SmoothInput,
    kParam_ShiftInput,
    
    // Outputs with mode (each WITH_MODE macro creates 2 params)
    kParam_Output1,
    kParam_Output1Mode,
    kParam_Output2,
    kParam_Output2Mode,
    kParam_Output3,
    kParam_Output3Mode,
    kParam_Output4,
    kParam_Output4Mode,
    
    // Mode controls
    kParam_RampMode,      // AD / Cycle / AR
    kParam_Range,         // Low / Medium / High
    kParam_OutputMode,    // Gates / Amplitude / Slope/Phase / Frequency
    
    // Main parameters
    kParam_Frequency,
    kParam_Shape,
    kParam_Slope,
    kParam_Smoothness,
    kParam_Shift,
    
    // Modulation attenuverters
    kParam_FMAmount,
    kParam_ShapeAtten,
    kParam_SlopeAtten,
    kParam_SmoothAtten,
    kParam_ShiftAtten,
    
    kNumParams
};

// Mode enum strings
static const char* rampModeNames[] = { "AD", "Cycle", "AR", NULL };
static const char* rangeNames[] = { "Low", "Medium", "High", NULL };
static const char* outputModeNames[] = { "Gates", "Amplitude", "Slope/Phase", "Frequency", NULL };

// Parameter definitions
static const _NT_parameter parameters[] = {
    // Inputs (page 1)
    // FIXED: Default Trig/Gate to 0 (None) so Cycle mode free-runs
    NT_PARAMETER_CV_INPUT("Trig/Gate In", 0, 0)
    NT_PARAMETER_CV_INPUT("Clock In", 0, 0)
    NT_PARAMETER_CV_INPUT("V/Oct In", 0, 0)
    NT_PARAMETER_CV_INPUT("FM In", 0, 0)
    NT_PARAMETER_CV_INPUT("Shape In", 0, 0)
    NT_PARAMETER_CV_INPUT("Slope In", 0, 0)
    NT_PARAMETER_CV_INPUT("Smooth In", 0, 0)
    NT_PARAMETER_CV_INPUT("Shift In", 0, 0)
    
    // Outputs with mode (page 2)
    NT_PARAMETER_CV_OUTPUT_WITH_MODE("Output 1", 1, 1)
    NT_PARAMETER_CV_OUTPUT_WITH_MODE("Output 2", 1, 2)
    NT_PARAMETER_CV_OUTPUT_WITH_MODE("Output 3", 1, 3)
    NT_PARAMETER_CV_OUTPUT_WITH_MODE("Output 4", 1, 4)
    
    // Mode controls (page 3)
    { .name = "Ramp Mode", .min = 0, .max = 2, .def = 1, .unit = kNT_unitEnum, .scaling = 0, .enumStrings = rampModeNames },
    { .name = "Range", .min = 0, .max = 2, .def = 1, .unit = kNT_unitEnum, .scaling = 0, .enumStrings = rangeNames },
    { .name = "Output Mode", .min = 0, .max = 3, .def = 0, .unit = kNT_unitEnum, .scaling = 0, .enumStrings = outputModeNames },
    
    // Main parameters (page 4) - use 0-1000 with scaling for smooth control
    { .name = "Frequency", .min = -500, .max = 500, .def = 0, .unit = kNT_unitNone, .scaling = kNT_scaling100, .enumStrings = NULL },
    { .name = "Shape", .min = 0, .max = 1000, .def = 500, .unit = kNT_unitPercent, .scaling = kNT_scaling10, .enumStrings = NULL },
    { .name = "Slope", .min = 0, .max = 1000, .def = 500, .unit = kNT_unitPercent, .scaling = kNT_scaling10, .enumStrings = NULL },
    { .name = "Smoothness", .min = 0, .max = 1000, .def = 500, .unit = kNT_unitPercent, .scaling = kNT_scaling10, .enumStrings = NULL },
    { .name = "Shift", .min = 0, .max = 1000, .def = 500, .unit = kNT_unitPercent, .scaling = kNT_scaling10, .enumStrings = NULL },
    
    // Modulation attenuverters (page 5)
    { .name = "FM Amount", .min = -100, .max = 100, .def = 0, .unit = kNT_unitPercent, .scaling = 0, .enumStrings = NULL },
    { .name = "Shape Atten", .min = -100, .max = 100, .def = 50, .unit = kNT_unitPercent, .scaling = 0, .enumStrings = NULL },
    { .name = "Slope Atten", .min = -100, .max = 100, .def = 50, .unit = kNT_unitPercent, .scaling = 0, .enumStrings = NULL },
    { .name = "Smooth Atten", .min = -100, .max = 100, .def = 50, .unit = kNT_unitPercent, .scaling = 0, .enumStrings = NULL },
    { .name = "Shift Atten", .min = -100, .max = 100, .def = 50, .unit = kNT_unitPercent, .scaling = 0, .enumStrings = NULL },
};

// Parameter pages (using arrays of indices)
static const uint8_t pageInputs[] = { 
    kParam_TrigInput, kParam_ClockInput, kParam_VOctInput, kParam_FMInput,
    kParam_ShapeInput, kParam_SlopeInput, kParam_SmoothInput, kParam_ShiftInput
};
static const uint8_t pageOutputs[] = {
    kParam_Output1, kParam_Output1Mode, kParam_Output2, kParam_Output2Mode,
    kParam_Output3, kParam_Output3Mode, kParam_Output4, kParam_Output4Mode
};
static const uint8_t pageMode[] = {
    kParam_RampMode, kParam_Range, kParam_OutputMode
};
static const uint8_t pageMain[] = {
    kParam_Frequency, kParam_Shape, kParam_Slope, kParam_Smoothness, kParam_Shift
};
static const uint8_t pageMod[] = {
    kParam_FMAmount, kParam_ShapeAtten, kParam_SlopeAtten, kParam_SmoothAtten, kParam_ShiftAtten
};

static const _NT_parameterPage pages[] = {
    { .name = "Inputs", .numParams = ARRAY_SIZE(pageInputs), .params = pageInputs },
    { .name = "Outputs", .numParams = ARRAY_SIZE(pageOutputs), .params = pageOutputs },
    { .name = "Mode", .numParams = ARRAY_SIZE(pageMode), .params = pageMode },
    { .name = "Main", .numParams = ARRAY_SIZE(pageMain), .params = pageMain },
    { .name = "Modulation", .numParams = ARRAY_SIZE(pageMod), .params = pageMod },
};

static const _NT_parameterPages parameterPages = {
    .numPages = ARRAY_SIZE(pages),
    .pages = pages,
};

// ============================================================================
// Frequency Constants
// ============================================================================

// Frequency range base frequencies (Hz)
static const float kFrequencyBase[3] = {
    0.125f,      // Low: ~0.125 Hz base
    2.0f,        // Medium: ~2 Hz base
    130.81f      // High: C3 = 130.81 Hz base
};

// ============================================================================
// Utility Functions
// ============================================================================

static inline float SemitonesToRatio(float semitones) {
    float octaves = semitones * (1.0f / 12.0f);
    
    int32_t octave_int = static_cast<int32_t>(octaves);
    float frac = octaves - octave_int;
    if (frac < 0.0f) {
        frac += 1.0f;
        octave_int -= 1;
    }
    
    // 2^frac approximation
    float result = 1.0f + frac * (0.6931472f + frac * (0.2402265f + frac * 0.0558f));
    
    if (octave_int > 0) {
        for (int i = 0; i < octave_int && i < 10; ++i) result *= 2.0f;
    } else if (octave_int < 0) {
        for (int i = 0; i < -octave_int && i < 10; ++i) result *= 0.5f;
    }
    
    return result;
}

// ============================================================================
// Plugin Callbacks
// ============================================================================

void calculateRequirements(_NT_algorithmRequirements& req, const int32_t* specifications)
{
    req.numParameters = ARRAY_SIZE(parameters);
    req.sram = sizeof(_TidesAlgorithm);
    req.dram = 0;
    req.dtc = 0;
    req.itc = 0;
}

_NT_algorithm* construct(const _NT_algorithmMemoryPtrs& ptrs, const _NT_algorithmRequirements& req, const int32_t* specifications)
{
    _TidesAlgorithm* alg = new (ptrs.sram) _TidesAlgorithm();
    alg->parameters = parameters;
    alg->parameterPages = &parameterPages;
    
    // Initialize DSP engine
    alg->poly_slope_generator.Init();
    
    // Initialize state
    alg->sample_rate = NT_globals.sampleRate;
    alg->inv_sample_rate = 1.0f / NT_globals.sampleRate;
    
    alg->prev_trig_high = false;
    alg->prev_clock_high = false;
    alg->clock_phase = 0.0f;
    alg->clock_period = 0.0f;
    alg->clock_counter = 0;
    alg->clock_synced = false;
    alg->prev_output_mode = -1;
    
    alg->cached_frequency = 1.0f;
    alg->cached_shape = 0.5f;
    alg->cached_slope = 0.5f;
    alg->cached_smooth = 0.5f;
    alg->cached_shift = 0.5f;
    
    memset(alg->output_buffer, 0, sizeof(alg->output_buffer));
    memset(alg->gate_buffer, 0, sizeof(alg->gate_buffer));
    memset(alg->ramp_buffer, 0, sizeof(alg->ramp_buffer));
    
    return alg;
}

void parameterChanged(_NT_algorithm* self, int p)
{
    _TidesAlgorithm* pThis = (_TidesAlgorithm*)self;
    
    switch (p) {
        case kParam_Shape:
            pThis->cached_shape = pThis->v[kParam_Shape] / 1000.0f;
            break;
        case kParam_Slope:
            pThis->cached_slope = pThis->v[kParam_Slope] / 1000.0f;
            break;
        case kParam_Smoothness:
            pThis->cached_smooth = pThis->v[kParam_Smoothness] / 1000.0f;
            break;
        case kParam_Shift:
            pThis->cached_shift = pThis->v[kParam_Shift] / 1000.0f;
            break;
        case kParam_OutputMode:
            // Reset DSP when output mode changes
            if (pThis->v[kParam_OutputMode] != pThis->prev_output_mode) {
                pThis->poly_slope_generator.Reset();
                pThis->prev_output_mode = pThis->v[kParam_OutputMode];
            }
            break;
    }
}

void step(_NT_algorithm* self, float* busFrames, int numFramesBy4)
{
    _TidesAlgorithm* pThis = (_TidesAlgorithm*)self;
    int numFrames = numFramesBy4 * 4;
    
    // Get bus assignments
    const int trigBus = pThis->v[kParam_TrigInput];
    const int clockBus = pThis->v[kParam_ClockInput];
    const int voctBus = pThis->v[kParam_VOctInput];
    const int fmBus = pThis->v[kParam_FMInput];
    const int shapeBus = pThis->v[kParam_ShapeInput];
    const int slopeBus = pThis->v[kParam_SlopeInput];
    const int smoothBus = pThis->v[kParam_SmoothInput];
    const int shiftBus = pThis->v[kParam_ShiftInput];
    
    const int out1Bus = pThis->v[kParam_Output1];
    const int out2Bus = pThis->v[kParam_Output2];
    const int out3Bus = pThis->v[kParam_Output3];
    const int out4Bus = pThis->v[kParam_Output4];
    
    const bool out1Replace = pThis->v[kParam_Output1Mode];
    const bool out2Replace = pThis->v[kParam_Output2Mode];
    const bool out3Replace = pThis->v[kParam_Output3Mode];
    const bool out4Replace = pThis->v[kParam_Output4Mode];
    
    // Get mode settings
    const tides::RampMode ramp_mode = static_cast<tides::RampMode>(pThis->v[kParam_RampMode]);
    const int range_idx = pThis->v[kParam_Range];
    const tides::Range range = range_idx < 2 ? tides::RANGE_CONTROL : tides::RANGE_AUDIO;
    const tides::OutputMode output_mode = static_cast<tides::OutputMode>(pThis->v[kParam_OutputMode]);
    
    // Get base parameters
    const float frequency_param = pThis->v[kParam_Frequency] / 100.0f;  // -5 to +5
    const float fm_atten = pThis->v[kParam_FMAmount] / 100.0f;
    const float shape_atten = pThis->v[kParam_ShapeAtten] / 100.0f;
    const float slope_atten = pThis->v[kParam_SlopeAtten] / 100.0f;
    const float smooth_atten = pThis->v[kParam_SmoothAtten] / 100.0f;
    const float shift_atten = pThis->v[kParam_ShiftAtten] / 100.0f;
    
    // Get input pointers - ONLY if bus is assigned (bus > 0)
    const float* trigIn = (trigBus > 0) ? busFrames + (trigBus - 1) * numFrames : nullptr;
    const float* clockIn = (clockBus > 0) ? busFrames + (clockBus - 1) * numFrames : nullptr;
    const float* voctIn = (voctBus > 0) ? busFrames + (voctBus - 1) * numFrames : nullptr;
    const float* fmIn = (fmBus > 0) ? busFrames + (fmBus - 1) * numFrames : nullptr;
    const float* shapeIn = (shapeBus > 0) ? busFrames + (shapeBus - 1) * numFrames : nullptr;
    const float* slopeIn = (slopeBus > 0) ? busFrames + (slopeBus - 1) * numFrames : nullptr;
    const float* smoothIn = (smoothBus > 0) ? busFrames + (smoothBus - 1) * numFrames : nullptr;
    const float* shiftIn = (shiftBus > 0) ? busFrames + (shiftBus - 1) * numFrames : nullptr;
    
    // Process gate/trigger inputs
    for (int i = 0; i < numFrames; ++i) {
        tides::GateFlags flags = tides::GATE_FLAG_LOW;
        
        if (trigIn) {
            // External trigger connected - use it
            bool trig_high = trigIn[i] > 1.0f;
            if (trig_high && !pThis->prev_trig_high) {
                flags = static_cast<tides::GateFlags>(flags | tides::GATE_FLAG_RISING);
            }
            if (trig_high) {
                flags = static_cast<tides::GateFlags>(flags | tides::GATE_FLAG_HIGH);
            }
            pThis->prev_trig_high = trig_high;
        } else if (ramp_mode == tides::RAMP_MODE_LOOPING) {
            // Cycle mode with no external trigger - just keep gate high for free-running
            // No rising edge needed - the DSP increments phase continuously in looping mode
            flags = tides::GATE_FLAG_HIGH;
        }
        
        pThis->gate_buffer[i] = flags;
        
        // Process clock for ramp extraction
        if (clockIn) {
            bool clock_high = clockIn[i] > 1.0f;
            if (clock_high && !pThis->prev_clock_high) {
                if (pThis->clock_counter > 10) {
                    pThis->clock_period = static_cast<float>(pThis->clock_counter);
                    pThis->clock_synced = true;
                }
                pThis->clock_counter = 0;
                pThis->clock_phase = 0.0f;
            }
            pThis->prev_clock_high = clock_high;
            pThis->clock_counter++;
            
            if (pThis->clock_synced && pThis->clock_period > 0.0f) {
                pThis->clock_phase += 1.0f / pThis->clock_period;
                if (pThis->clock_phase >= 1.0f) {
                    pThis->clock_phase -= 1.0f;
                }
                pThis->ramp_buffer[i] = pThis->clock_phase;
            } else {
                pThis->ramp_buffer[i] = 0.0f;
            }
        } else {
            pThis->ramp_buffer[i] = 0.0f;
            pThis->clock_synced = false;
        }
    }
    
    // Calculate frequency (V/Oct: 5V = 1.0 in Disting NT)
    float voct_cv = voctIn ? voctIn[0] * 0.2f * 12.0f : 0.0f;  // Scale to semitones
    float fm_cv = fmIn ? fmIn[0] * 0.2f * 12.0f * fm_atten : 0.0f;
    
    float transposition = frequency_param * 12.0f + voct_cv + fm_cv;
    float base_freq = kFrequencyBase[range_idx] * pThis->inv_sample_rate;
    float frequency = base_freq * SemitonesToRatio(transposition);
    
    // Smooth parameters toward target values (one-pole lowpass)
    // This prevents zipper noise from stepped parameter changes
    const float smoothingCoeff = 0.02f;  // ~20ms smoothing at 48kHz buffer rate
    float target_shape = pThis->v[kParam_Shape] / 1000.0f;
    float target_slope = pThis->v[kParam_Slope] / 1000.0f;
    float target_smooth = pThis->v[kParam_Smoothness] / 1000.0f;
    float target_shift = pThis->v[kParam_Shift] / 1000.0f;
    
    pThis->cached_shape += smoothingCoeff * (target_shape - pThis->cached_shape);
    pThis->cached_slope += smoothingCoeff * (target_slope - pThis->cached_slope);
    pThis->cached_smooth += smoothingCoeff * (target_smooth - pThis->cached_smooth);
    pThis->cached_shift += smoothingCoeff * (target_shift - pThis->cached_shift);
    
    // Apply CV modulation to smoothed parameters
    float shape = pThis->cached_shape;
    float slope = pThis->cached_slope;
    float smoothness = pThis->cached_smooth;
    float shift = pThis->cached_shift;
    
    if (shapeIn) shape += shapeIn[0] * 0.1f * shape_atten;
    if (slopeIn) slope += slopeIn[0] * 0.1f * slope_atten;
    if (smoothIn) smoothness += smoothIn[0] * 0.1f * smooth_atten;
    if (shiftIn) shift += shiftIn[0] * 0.1f * shift_atten;
    
    // Clamp parameters
    shape = std::max(0.0f, std::min(1.0f, shape));
    slope = std::max(0.0f, std::min(1.0f, slope));
    smoothness = std::max(0.0f, std::min(1.0f, smoothness));
    shift = std::max(0.0f, std::min(1.0f, shift));
    
    // Run DSP
    // In Cycle mode, ignore external clock so Frequency knob always controls speed
    const float* ramp_ptr = 
        (clockIn && pThis->clock_synced && ramp_mode != tides::RAMP_MODE_LOOPING) 
            ? pThis->ramp_buffer 
            : nullptr;
    
    pThis->poly_slope_generator.Render(
        ramp_mode,
        output_mode,
        range,
        frequency,
        slope,
        shape,
        smoothness,
        shift,
        pThis->gate_buffer,
        ramp_ptr,
        pThis->output_buffer,
        numFrames
    );
    
    // Write outputs
    // Output scaling: Tides DSP outputs ~0-1 unipolar or ~±1 bipolar
    // Disting NT uses voltage directly (1.0 = 1V), so scale for ~±5V range
    const float outputScale = 5.0f;
    
    if (out1Bus) {
        float* out = busFrames + (out1Bus - 1) * numFrames;
        if (out1Replace) {
            for (int i = 0; i < numFrames; ++i)
                out[i] = pThis->output_buffer[i].channel[0] * outputScale;
        } else {
            for (int i = 0; i < numFrames; ++i)
                out[i] += pThis->output_buffer[i].channel[0] * outputScale;
        }
    }
    if (out2Bus) {
        float* out = busFrames + (out2Bus - 1) * numFrames;
        if (out2Replace) {
            for (int i = 0; i < numFrames; ++i)
                out[i] = pThis->output_buffer[i].channel[1] * outputScale;
        } else {
            for (int i = 0; i < numFrames; ++i)
                out[i] += pThis->output_buffer[i].channel[1] * outputScale;
        }
    }
    if (out3Bus) {
        float* out = busFrames + (out3Bus - 1) * numFrames;
        if (out3Replace) {
            for (int i = 0; i < numFrames; ++i)
                out[i] = pThis->output_buffer[i].channel[2] * outputScale;
        } else {
            for (int i = 0; i < numFrames; ++i)
                out[i] += pThis->output_buffer[i].channel[2] * outputScale;
        }
    }
    if (out4Bus) {
        float* out = busFrames + (out4Bus - 1) * numFrames;
        if (out4Replace) {
            for (int i = 0; i < numFrames; ++i)
                out[i] = pThis->output_buffer[i].channel[3] * outputScale;
        } else {
            for (int i = 0; i < numFrames; ++i)
                out[i] += pThis->output_buffer[i].channel[3] * outputScale;
        }
    }
}

// ============================================================================
// Factory Definition
// ============================================================================

static const _NT_factory factory = {
    .guid = NT_MULTICHAR('T', 'i', 'd', 'S'),
    .name = "Tides",
    .description = "Tidal Modulator - MI Tides 2 port",
    .numSpecifications = 0,
    .specifications = nullptr,
    .calculateStaticRequirements = nullptr,
    .initialise = nullptr,
    .calculateRequirements = calculateRequirements,
    .construct = construct,
    .parameterChanged = parameterChanged,
    .step = step,
    .draw = nullptr,
    .midiRealtime = nullptr,
    .midiMessage = nullptr,
    .tags = kNT_tagUtility,
    .hasCustomUi = nullptr,
    .customUi = nullptr,
    .setupUi = nullptr,
    .serialise = nullptr,
    .deserialise = nullptr,
    .midiSysEx = nullptr,
    .parameterUiPrefix = nullptr,
};

// ============================================================================
// Plugin Entry Point
// ============================================================================

uintptr_t pluginEntry(_NT_selector selector, uint32_t data)
{
    switch (selector) {
        case kNT_selector_version:
            return kNT_apiVersionCurrent;
        case kNT_selector_numFactories:
            return 1;
        case kNT_selector_factoryInfo:
            return (uintptr_t)((data == 0) ? &factory : nullptr);
    }
    return 0;
}
