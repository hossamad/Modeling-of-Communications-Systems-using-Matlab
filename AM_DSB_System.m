%When the program runs, there are four sound functions will run for 4 required cases..
%..and they (and their comments) at line number 283, 335, 385, 437 respectively.
clc
clear all;
%Declaring the signals
[QuranSignal,Quran_Fs]=audioread('Short_QuranPalestine.wav');
[BBCSignal,BBC_Fs]=audioread('Short_BBCArabic2.wav');
[SkyNewsSignal,SkyNews_Fs]=audioread('Short_SkyNewsArabia.wav');
[FM9090Signal,FM9090_Fs]=audioread('Short_FM9090.wav');
%Getting number of samples in each signal and collecting these numbers in an one vector
SamplingVector=[length(QuranSignal),length(BBCSignal),length(SkyNewsSignal),length(FM9090Signal)];
%Getting the largest number of samples
largest_length=max(SamplingVector);
%these lines i used them to compare the signals in different locations in the code with the original ones
% sound(QuranSignal,Quran_Fs);
% pause(18);
% sound(BBCSignal,Quran_Fs);
% pause(18);
% sound(SkyNewsSignal,Quran_Fs);
% pause(17);
% sound(FM9090Signal,Quran_Fs);
% pause(15);
%Forcing the sampling signals to be equal in the length
for i = 1:4
%Initiate a zero vector that has a length equal to the largest length
ExtendedVector=zeros(largest_length,1);
%Initiate an matrix which will contain the signal (that will be extended)
if SamplingVector(i) == length(QuranSignal)
TemporaryVector=QuranSignal;
elseif SamplingVector(i) == length(BBCSignal)
TemporaryVector=BBCSignal;
elseif SamplingVector(i) == length(SkyNewsSignal)
TemporaryVector=SkyNewsSignal;
else
TemporaryVector=FM9090Signal;
end
%extended the signal to the largest length and add the two channel together
for k = 1:SamplingVector(i)
ExtendedVector(k) = TemporaryVector(k,1) ;
end
%make vectors contain the modified signals
if SamplingVector(i) == length(QuranSignal)
Extended_QuranSignal=ExtendedVector;
elseif SamplingVector(i) == length(BBCSignal)
Extended_BBCSignal=ExtendedVector;
elseif SamplingVector(i) == length(SkyNewsSignal)
Extended_SkyNewsSignal=ExtendedVector;
else
Extended_FM9090Signal=ExtendedVector;
end
end
%as all signal have the same sampling frequency (44.1*10^3),
%I will increase that sampling frequency to become larger than the largest
%used carrier frequency over 2.
Largest_Fn=100*10^(3)+3*50*10^(3);
%round the number of times to the nearest integer and multipling it by 2
NumberOfTimes = 2 * round(Largest_Fn/(44.1*10^(3)));
%getting the new sampling freq.
NewFSValue = NumberOfTimes * 44.1*10^(3);
%increasing the sample number for the 4 signal
Interp_QuranSignalsig=interp(Extended_QuranSignal,NumberOfTimes);
Interp_BBCSignal=interp(Extended_BBCSignal,NumberOfTimes);
Interp_SkyNewsSignal=interp(Extended_SkyNewsSignal,NumberOfTimes);
Interp_FM9090Signal=interp(Extended_FM9090Signal,NumberOfTimes);
%Getting number of sample of signals (all sigs. have the same number of the samples).
NumOfSamples = length(Interp_QuranSignalsig);
%creating the frequency axis to map the samples to its frequency.
Freq = (-NumOfSamples/2 : NumOfSamples/2-1) * (NewFSValue/NumOfSamples);
%% Amplitude Modulation part
%initiate a zero matrix which will contain the 4 signals,
%where column 1, 2, 3 and 4 represents QuranSignal, BBCSignal, SkyNewsSignal and FM9090Signal respectively.
FinalShapeToSignals=zeros(NumOfSamples,4);
FinalShapeToSignals(:,1)=Interp_QuranSignalsig;
FinalShapeToSignals(:,2)=Interp_BBCSignal;
FinalShapeToSignals(:,3)=Interp_SkyNewsSignal;
FinalShapeToSignals(:,4)=Interp_FM9090Signal;
% sound(Interp_QuranSignalsig,NewFSValue);
% pause(15);
% sound(Interp_BBCSignal,NewFSValue);
% pause(15);
% sound(Interp_SkyNewsSignal,NewFSValue);
% pause(15);
% sound(Interp_FM9090Signal,NewFSValue);
%Gettting the FFT for the given signals
fft_Interp_QuranSignal = fftshift( abs( fft(FinalShapeToSignals(:,1) ) ) );
fft_Interp_BBCSignal = fftshift( abs( fft(FinalShapeToSignals(:,2) ) ) );
fft_Interp_SkyNewsSignal = fftshift( abs( fft(FinalShapeToSignals(:,3) ) ) );
fft_Interp_FM9090Signal = fftshift( abs( fft(FinalShapeToSignals(:,4) ) ) );
%plotting them versus the frequency
figure(1)
%title('Spectrums Of The Signals Versus Frequency (in Hz)')
subplot(2,2,1)
plot(Freq,fft_Interp_QuranSignal);
xlim ([-50000 50000]);
ylabel('QuranSignalSpectrum');
xlabel('Frequency (in Hz)');
subplot(2,2,2)
plot(Freq,fft_Interp_BBCSignal);
xlim ([-50000 50000]);
ylabel('BBCSignalSpectrum');
xlabel('Frequency (in Hz)');
subplot(2,2,3)
plot(Freq,fft_Interp_SkyNewsSignal);
xlim ([-50000 50000]);
ylabel('SkyNewsSignalSpectrum');
xlabel('Frequency (in Hz)');
subplot(2,2,4)
plot(Freq,fft_Interp_FM9090Signal );
xlim ([-50000 50000]);
ylabel('FM9090SignalSpectrum');
xlabel('Frequency (in Hz)');
%initiate a zero matrix which will contain the 4 modulated signals,
%where column 1, 2, 3 and 4 represents QuranSignal, BBCSignal, SkyNewsSignal and FM9090Signal respectively.
ModulatedSignals=zeros(NumOfSamples,4);
%doing the AM to all signals:Quran, BBC, SkyNews and FM9090 respectively.
Ts=1/NewFSValue;
for i = 1:4
%declare the carrier frequency
Wn(i)=(100*10^(3)+(i-1)*50*10^(3))*2*pi;
%Multiple the signal by its carrier
%(with amplitude equal 2 to maintain the amplitude of the modulated sig. as it is)
for k = 1:NumOfSamples
ModulatedSignals(k,i) = 2*cos(Wn(i)*k*Ts) * FinalShapeToSignals(k,i);
end
end
%Gettting the FFT for the given signals
fft_Modulated_QuranSignal = fftshift( abs( fft(ModulatedSignals(:,1)) ) );
fft_Modulated_BBCSignal = fftshift( abs( fft(ModulatedSignals(:,2)) ) );
fft_Modulated_SkyNewsSignal = fftshift( abs( fft(ModulatedSignals(:,3)) ) );
fft_Modulated_FM9090Signal = fftshift( abs( fft(ModulatedSignals(:,4)) ) );
figure(2)
subplot(2,2,1)
plot(Freq,fft_Modulated_QuranSignal);
%title('Modulated Signals Versus Frequency (in Hz)')
ylabel('ModulatedQuranSignalSpectrum');
xlabel('Frequency (in Hz)');
subplot(2,2,2)
plot(Freq,fft_Modulated_BBCSignal);
ylabel('ModulatedBBCSignalSpectrum');
xlabel('Frequency (in Hz)');
subplot(2,2,3)
plot(Freq,fft_Modulated_SkyNewsSignal);
ylabel('ModulatedSkyNewsSignalSpectrum');
xlabel('Frequency (in Hz)');
subplot(2,2,4)
plot(Freq,fft_Modulated_FM9090Signal );
ylabel('ModulatedFM9090SignalSpectrum');
xlabel('Frequency (in Hz)');
%Add the modulated signal together
SummingOfModulatedSignals=zeros(NumOfSamples,1);
for i = 1:NumOfSamples
for k = 1:4
SummingOfModulatedSignals(i)=SummingOfModulatedSignals(i)+ModulatedSignals(i,k);
end
end
fft_SummingOfModulatedSignals=fftshift( abs( fft(SummingOfModulatedSignals) ) );
figure(3)
plot(Freq,fft_SummingOfModulatedSignals);title('Summation Of Modulated Signals');
xlabel('Frequency (in Hz)');
ylabel('Amplitude Of Modulated Signals');
%% The RF stage:
%from the spectrum of the 4 signals, the BW of the RF BPF is 2*15*10^(3) Hz..
%..where 15*10^(3) is the largest frequency in the 4 spectrums.
%And here is about 8 KHz between every modulated signal, so the transition..
%..of the RF BPF could have about 4 KHz.
%RF BPF:
%declare the tunable centered frequency.
h=0:3;
Fn=100*10^(3)+h*50*10^(3);
BW_RF_BPF=2*15*10^(3); %in Hz.
%detertime the specifications of RF BPF based on the declared tunable ..
%.. centered frequency and the half of BW (which was determined from the signal spectrums).
i=1;
A_stop1 = 60; % Attenuation in the first stopband = 60 dB
F_stop1 = Fn(i) - BW_RF_BPF/2 - 8000 ; % Edge of the stopband in Hz
F_pass1 = Fn(i) - BW_RF_BPF/2; % Edge of the passband in Hz
F_pass2 = Fn(i) + BW_RF_BPF/2; % Closing edge of the passband in Hz
F_stop2 = Fn(i) + BW_RF_BPF/2 + 8000; % Edge of the second stopband in Hz
A_stop2 = 60; % Attenuation in the second stopband = 60 dB
A_pass = 1; % Amount of ripple allowed in the passband = 1 dB
d = fdesign.bandpass;
BandPassSpecObj = ...
fdesign.bandpass('Fst1,Fp1,Fp2,Fst2,Ast1,Ap,Ast2', ...
F_stop1, F_pass1, F_pass2, F_stop2, A_stop1, A_pass, ...
A_stop2, NewFSValue);
RF_BPF = design(BandPassSpecObj, 'equiripple');
%IF BPF:
%There is about 55 KHz between the beginning of the first signal, so the IR_BPF won't be sharpe as RF_BPF.
F_IF=25*10^(3);
A_stop1_IF = 60; % Attenuation in the first stopband = 60 dB
F_stop1_IF = F_IF - BW_RF_BPF/2 - 9000; % Edge of the stopband in Hz
F_pass1_IF = F_IF - BW_RF_BPF/2; % Edge of the passband in Hz
F_pass2_IF = F_IF + BW_RF_BPF/2 ; % Closing edge of the passband in Hz
F_stop2_IF = F_IF + BW_RF_BPF/2 + 9000; % Edge of the second stopband in Hz
A_stop2_IF = 60; % Attenuation in the second stopband = 60 dB
A_pass_IF = 1; % Amount of ripple allowed in the passband = 1 dB
d_IF = fdesign.bandpass;
BandPassSpecObj_IF = ...
fdesign.bandpass('Fst1,Fp1,Fp2,Fst2,Ast1,Ap,Ast2', ...
F_stop1_IF, F_pass1_IF, F_pass2_IF, F_stop2_IF, A_stop1_IF, A_pass_IF, ...
A_stop2_IF, NewFSValue);
IF_BPF = design(BandPassSpecObj_IF, 'equiripple');
%Moving the first signal (Quran signal) to IF:
RF_QuranSignal=filter(RF_BPF,SummingOfModulatedSignals);
%plot the spectrum of the RF_PBF output
fft_RF_QuranSignal=fftshift( abs( fft(RF_QuranSignal) ) );
figure(4)
plot(Freq,fft_RF_QuranSignal);
title('RF BPF Output')
xlabel('Frequency (in Hz)');
ylabel('Magitude');
%declare the carrier frequency of RF stage
Wn=(100*10^(3)+F_IF)*2*pi;
%Multiple the signal by its carrier
%(with amplitude equal 2 to maintain the amplitude of the modulated sig. as it is)
for k = 1:NumOfSamples
Demodulated_QuranSignal(k) = 2*cos(Wn*k*Ts) * RF_QuranSignal(k);
end
%plot the spectrum of the RF_Mixer output
fft_Demodulated_QuranSignal=fftshift( abs( fft(Demodulated_QuranSignal) ) );
figure(5)
plot(Freq,fft_Demodulated_QuranSignal);
title('RF Mixer Output')
xlabel('Frequency (in Hz)');
ylabel('Magitude');
%enter the Demodulated_QuranSignal to IF_PBF
IF_QuranSignal=filter(IF_BPF,Demodulated_QuranSignal);
%plot the spectrum of the IF_BPF output
fft_IF_QuranSignal=fftshift( abs( fft(IF_QuranSignal) ) );
figure(6)
plot(Freq,fft_IF_QuranSignal);title('IF BPF Output')
xlabel('Frequency (in Hz)');
ylabel('Magitude');
%% base band stage(With RF BPF Existing)
Wn=25*10^(3)*2*pi;
for k = 1:NumOfSamples
BaseBand_QuranSignal(k) = 2*cos(Wn*k*Ts) * IF_QuranSignal(k);
end
%plot the spectrum of the BB_STAGE_MIXER output
fft_BaseBand_QuranSignal=fftshift( abs( fft(BaseBand_QuranSignal) ) );
figure(7)
plot(Freq,fft_BaseBand_QuranSignal);
title('The Baseband Stage Mixer output')
xlabel('Frequency (in Hz)');
ylabel('Magitude');
%LPF
N = 120;
Fs = NewFSValue;
Fp = 15e3;
Ap = 0.01;
Ast = 80;
%Obtain the maximum deviation for the passband and stopband ripples in linear units.
Rp = (10^(Ap/20) - 1)/(10^(Ap/20) + 1);
Rst = 10^(-Ast/20);
LPF = firceqrip(N,Fp/(Fs/2),[Rp Rst],'passedge');
QuranSignal_out_from_LPF=conv(LPF,BaseBand_QuranSignal);
%plot the spectrum of the LPF output
fft_QuranSignal_out_from_LPF=fftshift( abs( fft(QuranSignal_out_from_LPF) ) );
N=length(fft_QuranSignal_out_from_LPF);
Freq2=(-N/2 : N/2-1) * (NewFSValue/N);
figure(8)
plot(Freq2,fft_QuranSignal_out_from_LPF);
title('The Baseband Stage LPF output')
xlabel('Frequency (in Hz)');
ylabel('Magitude');
%return the "QuranSignal_out_from_LPF" signal to its original sampling rate by resamply function.
%To use the resample function on uniform samples you must provide both the numerator and..
%..denominator of this rational factor. that is done by rat function.
originalFs = NewFSValue;
desiredFs = Quran_Fs;%the original sample frequency of the 'Short_QuranPalestine.wav' signal.
[p,q] = rat(desiredFs / originalFs);
%resamply the "QuranSignal_out_from_LPF" signal
Y = resample( QuranSignal_out_from_LPF, p, q );
sound(Y,desiredFs);% the demodulated signal in case there is RF BPF
pause(15);
%% demodulate the signal without RF BPF:
%declare the carrier frequency of RF stage
Wn=(100*10^(3)+F_IF)*2*pi;
%Multiple the signal by its carrier
%(with amplitude equal 2 to maintain the amplitude of the modulated sig. as it is)
for k = 1:NumOfSamples
DemodulatedSignals(k) = 2*cos(Wn*k*Ts) * SummingOfModulatedSignals(k);
end
%plot the spectrum of the RF_Mixer output
fft_DemodulatedSignals=fftshift( abs( fft(DemodulatedSignals) ) );
figure(9)
plot(Freq,fft_DemodulatedSignals);
title('RF Mixer Output (Without RF BPF)')
xlabel('Frequency (in Hz)');
ylabel('Magitude');
%enter the Demodulated_QuranSignal to IF_PBF
IF_QuranSignal_without_RF_BPF=filter(IF_BPF,DemodulatedSignals);
%plot the spectrum of the IF_BPF output
fft_IF_QuranSignal_without_RF_BPF =fftshift( abs( fft(IF_QuranSignal_without_RF_BPF) ) );
figure(10)
plot(Freq,fft_IF_QuranSignal_without_RF_BPF);title('IF BPF Output (Without RF BPF)')
xlabel('Frequency (in Hz)');
ylabel('Magitude');
%% base band stage(Without RF BPF Existing)
Wn=25*10^(3)*2*pi;
for k = 1:NumOfSamples
BaseBand_IF_QuranSignal_without_RF_BPF(k) = 2*cos(Wn*k*Ts) * IF_QuranSignal_without_RF_BPF(k);
end
%plot the spectrum of the BB_STAGE_MIXER output
fft_BaseBand_IF_QuranSignal_without_RF_BPF=fftshift( abs( fft(BaseBand_IF_QuranSignal_without_RF_BPF) ) );
figure(11)
plot(Freq,fft_BaseBand_IF_QuranSignal_without_RF_BPF);
title('The Baseband Stage Mixer output(Without RF BPF Existing)')
xlabel('Frequency (in Hz)');
ylabel('Magitude');
QuranSignal_out_from_LPF2=conv(LPF,BaseBand_IF_QuranSignal_without_RF_BPF);
%plot the spectrum of the LPF output
fft_QuranSignal_out_from_LPF2=fftshift( abs( fft(QuranSignal_out_from_LPF2) ) );
figure(12)
N=length(QuranSignal_out_from_LPF2);
Freq2=(-N/2 : N/2-1) * (NewFSValue/N);
plot(Freq2,fft_QuranSignal_out_from_LPF2);title('The Baseband Stage LPF output(Without RF BPF Existing)')
xlabel('Frequency (in Hz)');
ylabel('Magitude');
Y2 = resample( QuranSignal_out_from_LPF2, p, q );
sound(Y2,desiredFs);% the demodulated signal in case there is no RF BPF
pause(15)
%% demodulate the signal with an offset (0.1 kHz) in the oscillator:
%declare the carrier frequency of RF stage
Wn=(100*10^(3)+ F_IF + 0.1*10^(3))*2*pi;
%Multiple the signal by its carrier
for k = 1:NumOfSamples
DemodulatedSignals3(k) = 2*cos(Wn*k*Ts) * RF_QuranSignal(k);
end
fft_RF_QuranSignal3=fftshift( abs( fft(DemodulatedSignals3) ) );
figure(13)
subplot(2,1,1)
plot(Freq,fft_RF_QuranSignal3);
title('Demodulating To The Quran Signal (With offset 0.1 kHz in the oscillator)')
xlabel('Frequency (in Hz)');
ylabel('RF Mixer Output');
%enter the Demodulated_QuranSignal to IF_PBF
IF_QuranSignal3=filter(IF_BPF,DemodulatedSignals3);
%plot the spectrum of the IF_BPF output
fft_IF_QuranSignal3=fftshift( abs( fft(IF_QuranSignal3) ) );
subplot(2,1,2)
plot(Freq,fft_IF_QuranSignal3);
xlabel('Frequency (in Hz)');
ylabel('IF BPF Output');
%% base band stage in case an offset (0.1 kHz) in the oscillator:
Wn_if=25*10^(3)*2*pi;
for k = 1:NumOfSamples
BaseBand_IF_QuranSignal3(k) = 2*cos(Wn_if*k*Ts) * IF_QuranSignal3(k);
end
%plot the spectrum of the BB_STAGE_MIXER output
fft_BaseBand_IF_QuranSignal3=fftshift( abs( fft(BaseBand_IF_QuranSignal3) ) );
figure(14)
subplot(2,1,1)
plot(Freq,fft_BaseBand_IF_QuranSignal3);
title('The Baseband Detection(With offset 0.1 kHz in the oscillator)')
xlabel('Frequency (in Hz)');
ylabel('BB Stage Mixer output');
QuranSignal_out_from_LPF3=conv(LPF,BaseBand_IF_QuranSignal3);
%plot the spectrum of the LPF output
fft_QuranSignal_out_from_LPF2=fftshift( abs( fft(QuranSignal_out_from_LPF3) ) );
subplot(2,1,2)
N=length(QuranSignal_out_from_LPF3);
Freq2=(-N/2 : N/2-1) * (NewFSValue/N);
plot(Freq2,fft_QuranSignal_out_from_LPF2);
xlabel('Frequency (in Hz)');
ylabel('LPF output');
%played the signal
Y3 = resample( QuranSignal_out_from_LPF3, p, q );
% the demodulated signal in case there is offset (0.1 kHz) in the oscillator
sound(Y3,desiredFs);
pause(15);
%% demodulate the signal with an offset (1 kHz) in the oscillator:
%declare the carrier frequency of RF stage
Wn=(100*10^(3)+ F_IF + 1*10^(3))*2*pi;
%Multiple the signal by its carrier
for k = 1:NumOfSamples
DemodulatedSignals3(k) = 2*cos(Wn*k*Ts) * RF_QuranSignal(k);
end
fft_RF_QuranSignal3=fftshift( abs( fft(DemodulatedSignals3) ) );
figure(15)
subplot(2,1,1)
plot(Freq,fft_RF_QuranSignal3);
title('Demodulating To The Quran Signal (With offset 1 kHz in the oscillator)')
xlabel('Frequency (in Hz)');
ylabel('RF Mixer Output');
%enter the Demodulated_QuranSignal to IF_PBF
IF_QuranSignal3=filter(IF_BPF,DemodulatedSignals3);
%plot the spectrum of the IF_BPF output
fft_IF_QuranSignal3=fftshift( abs( fft(IF_QuranSignal3) ) );
subplot(2,1,2)
plot(Freq,fft_IF_QuranSignal3);
xlabel('Frequency (in Hz)');
ylabel('IF BPF Output');
%% base band stage in case an offset (1 kHz) in the RF oscillator:
Wn_if=25*10^(3)*2*pi;
for k = 1:NumOfSamples
BaseBand_IF_QuranSignal3(k) = 2*cos(Wn_if*k*Ts) * IF_QuranSignal3(k);
end
%plot the spectrum of the BB_STAGE_MIXER output
fft_BaseBand_IF_QuranSignal3=fftshift( abs( fft(BaseBand_IF_QuranSignal3) ) );
figure(16)
subplot(2,1,1)
plot(Freq,fft_BaseBand_IF_QuranSignal3);
title('The Baseband Detection(With offset 1 kHz in the oscillator)')
xlabel('Frequency (in Hz)');
ylabel('BB Stage Mixer output');
QuranSignal_out_from_LPF3=conv(LPF,BaseBand_IF_QuranSignal3);
%plot the spectrum of the LPF output
fft_QuranSignal_out_from_LPF2=fftshift( abs( fft(QuranSignal_out_from_LPF3) ) );
subplot(2,1,2)
N=length(QuranSignal_out_from_LPF3);
Freq2=(-N/2 : N/2-1) * (NewFSValue/N);
plot(Freq2,fft_QuranSignal_out_from_LPF2);
xlabel('Frequency (in Hz)');
ylabel('LPF output');
%played the signal
Y4 = resample( QuranSignal_out_from_LPF3, p, q );
% the demodulated signal in case there is offset (0.1 kHz) in the oscillator
sound(Y4,desiredFs);
%ploting the used filters
fvtool(RF_BPF);
fvtool(IF_BPF);
fvtool(LPF);