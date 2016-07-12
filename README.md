# EE3: Real-Time Digital Signal Processing
Real-Time Digital Signal Processing on a Texas Instruments DSP Starter Kit (TMS320C6713 DSK). Imperial College London Electrical &amp; Electronic Engineering 3<sup>rd</sup> year module.

## Lab 5: IIR Filtering

Objectives accomplished:
* Learned to design IIR (Infinite-Impulse Response) digital filters using MATLAB.
* Implemented the IIR filter using the C6713 DSK system in real-time.
* Measured the filter characteristics using a spectrum analyzer.

## Project: Speech Enhancement through noise reduction

> Telephones are increasingly being used in noisy environments such as cars, airports and undergraduate laboratories! The aim of this project is to implement a real-time system that will reduce the background noise in a speech signal while leaving the signal itself intact: this process is called speech enhancement.

~ Paul D. Mitcheson, Imperial College London Dept of EEE, *EE3-19 Real Time Digital Signal Processing*, 2016

Objectives accomplished:
* Implemented triple buffering: *in*, *processing*, and *out* buffers.
* Implemented various noise estimation and noise reduction techniques,
* Compared their performance and chose parameters that produced the best audible speech enhancement ([PDF](/speech-enhancement/EusebiusN_And_PrahnavS_RTDSP_Project.pdf)).
