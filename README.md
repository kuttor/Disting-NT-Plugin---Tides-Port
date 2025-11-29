# Tides 2 for Disting NT

**Authentic port of Mutable Instruments Tides 2 (2018 version)**

This is a complete port of the Mutable Instruments Tides 2 DSP code to the Expert Sleepers Disting NT platform. It uses the original wavetables, band-limited oscillators, ratio tables, and waveshaping algorithms.

## Features

### Ramp Modes
- **AD** - Attack/Decay envelope triggered by gate
- **Cycle** - Free-running LFO/VCO
- **AR** - Attack/Release envelope following gate

### Frequency Ranges
- **Low** - 0.125 Hz base (ultra-slow LFO)
- **Medium** - 2 Hz base (standard LFO)
- **High** - 130.81 Hz base (audio rate, C3)

### Output Modes
- **Gates** - Main shaped output + raw ramp + EOA gate + EOR gate
- **Amplitude** - Panning/crossfading between 4 outputs based on Shift
- **Slope/Phase** - 4 phase-shifted copies of the waveform
- **Frequency** - 4 polyrhythmic outputs at different frequency ratios

### Authentic DSP
- **12 wavetable shapes** - Original 12,300-sample wavetable
- **Band-limited slopes** - PolyBLEP anti-aliasing in audio range
- **Wavefolder** - Original lookup tables for bipolar/unipolar folding
- **Ratio tables** - 21 musical ratios for polyrhythmic mode
- **Parameter smoothing** - Per-sample interpolation

## Building

### Prerequisites

1. **ARM GCC Toolchain**
   ```bash
   # Ubuntu/Debian
   sudo apt install gcc-arm-none-eabi
   
   # macOS
   brew install arm-none-eabi-gcc
   
   # Windows
   # Download from: https://developer.arm.com/tools-and-software/open-source-software/developer-tools/gnu-toolchain/gnu-rm
   ```

2. **Disting NT API**
   ```bash
   git clone https://github.com/expertsleepersltd/distingNT_API
   ```

### Build

```bash
# Check toolchain
make check API_PATH=/path/to/distingNT_API

# Build
make API_PATH=/path/to/distingNT_API

# Syntax check only (no ARM toolchain needed)
make syntax API_PATH=/path/to/distingNT_API
```

### Install

Copy `tides.o` to your Disting NT SD card:
```
SD:/programs/plug-ins/tides.o
```

Or use:
```bash
make install MOUNT_POINT=/media/your-sd-card API_PATH=/path/to/distingNT_API
```

## Parameters

### Page 1: Inputs
| Parameter | Description |
|-----------|-------------|
| Trig/Gate In | Trigger (AD/Cycle) or Gate (AR mode) input |
| Clock In | External clock for tempo sync |
| V/Oct In | 1V/octave pitch CV |
| FM In | Frequency modulation CV |
| Shape In | Shape modulation CV |
| Slope In | Slope modulation CV |
| Smooth In | Smoothness modulation CV |
| Shift In | Shift modulation CV |

### Page 2: Outputs
| Parameter | Description |
|-----------|-------------|
| Output 1-4 | Main waveform outputs (meaning depends on Output Mode) |
| Output 1-4 Mode | Replace or Add to bus |

### Page 3: Mode
| Parameter | Options |
|-----------|---------|
| Ramp Mode | AD, Cycle, AR |
| Range | Low, Medium, High |
| Output Mode | Gates, Amplitude, Slope/Phase, Frequency |

### Page 4: Main Parameters
| Parameter | Range | Description |
|-----------|-------|-------------|
| Frequency | ±5 octaves | Base frequency offset |
| Shape | 0-100% | Waveform shape (exponential ↔ linear ↔ logarithmic) |
| Slope | 0-100% | Rise/fall time ratio (attack/decay asymmetry) |
| Smoothness | 0-100% | 0-50% = lowpass, 50-100% = wavefold |
| Shift | 0-100% | Output mode dependent (phase/pan/ratio) |

### Page 5: Modulation
| Parameter | Range | Description |
|-----------|-------|-------------|
| FM Amount | ±100% | FM input attenuverter |
| Shape Atten | ±100% | Shape CV attenuverter |
| Slope Atten | ±100% | Slope CV attenuverter |
| Smooth Atten | ±100% | Smoothness CV attenuverter |
| Shift Atten | ±100% | Shift CV attenuverter |

## Output Mode Details

### Gates Mode
- **Out 1**: Shaped waveform × Shift amount
- **Out 2**: Raw bipolar ramp (or shaped if audio rate cycling)
- **Out 3**: End-of-Attack gate
- **Out 4**: End-of-Release/End-of-Ramp gate

### Amplitude Mode
- Shift acts as a scanner/panner across 4 outputs
- All outputs receive the same shaped waveform with varying amplitude

### Slope/Phase Mode
- 4 copies of the waveform at different phases
- Shift controls the phase spread between outputs

### Frequency Mode
- 4 polyrhythmic outputs at musical ratios
- Shift selects ratio set (21 preset combinations)
- Includes unison, octaves, fifths, and complex polyrhythms

## Technical Details

- **Plugin GUID**: `TidS`
- **Sample Rate**: Uses NT's sample rate (typically 48kHz)
- **Output Levels**: ±5V bipolar, 0-8V unipolar (mode dependent)
- **Memory**: ~50KB for wavetables + DSP state

## Credits

- **Original Tides 2**: Émilie Gillet / Mutable Instruments
- **Disting NT Port**: Claude / Anthropic
- **Disting NT API**: Expert Sleepers
- **Andrew Kuttor

## License

MIT License

Original Tides 2 code: Copyright 2017 Émilie Gillet
Disting NT adaptation: Copyright 2024

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

---

*"Mutable Instruments" is a registered trademark. This derivative work is not
affiliated with or endorsed by Mutable Instruments.*
