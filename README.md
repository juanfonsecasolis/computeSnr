# ComputeSNR2 
Estimates the signal-to-noise ratio (SNR) given a set of manual annotations (in seconds).
    
This algorithm estimates the signal-to-noise ratio in recordings that contains manually annotations of the desired and noise segments. The assumption that the noise is completely uncorrelated from the desired signal is made (e.g. additive noise). Recordings without information regarding the sound intensity or the voltage-sound-intensity relationship in the recording devices may benefit.

The basis formula (eq. 1) is defined like so:
`SNR = 10 * log10((P(desired)−P(noise))/P(noise))`, (eq. 1)

where:
* P(s) is the power of a finite period signal and 
* P(s) = 1/N * (|s(1)|^2+|s(2)|^2+...+|s(N)|^2)

The desired signal is constructed by appending the entries of the mixed signal specified by the annotated segments (onset[i]-offset[i] pair), and the noise signal is constructed by appending the entries outside the annotations. Good annotations should make the noise signal contain only silence or not-desired segments. Turn flag DEBUG=true to save the segmented signals in two files called "desired.wav" and "noise.wav". This also prints a graph with the annotations and the oscillogram of the saved *.wav files.

## References:
1. Papadopoulos, P., & Tsiartas, A. (2014). A supervised signal-to-noise ratio estimation of desired signals. IEEE International Conference on Acoustics desired and Signal Processing, 8287–8291
2. Beritelli, F., Casale, S., Grasso, R., & Spadaccini, A. (2010). Performance Evaluation of SNR Estimation Methods in Forensic Speaker Recognition. In 2010 Fourth International Conference on Emerging Security Information, Systems and Technologies (pp. 88–92)
3. SNR estimation in segmented speech signals. URL: http://dsp.stackexchange.com/questions/35041/snr-estimation-in-segmented-desired-signals (last consulted on May 12th 2019)
   
## Syntax: 
* ComputeSNR2(signal, fs, annotatedOnsets, annotatedOffsets)

Inputs:
* signal - an array containing the mixed signals: noise and desired
* fs - sampling rate
* annotatedOnsets - array with manually annotated onsets of the desired signal
* annotatedOffsets - array with manually annotated offsets of the desired signal

Outputs:
* SNR - Signal-to-noise ratio estimate in dB

## Notes:
* Other m-files required: none
* Subfunctions: none
* MAT-files required: none

---

Author: Juan Manuel Fonseca-Solís
San José, Costa Rica
Email: juanma2268@gmail.com
Oct 2016; Last revision: 26-Oct-2016

