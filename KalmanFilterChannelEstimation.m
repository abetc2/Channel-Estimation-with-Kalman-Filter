%%%% Script for Channel Estimation with Kalman Filter for OFDM signal %%%%
%%%% written by Ayse Betul Demir
%%%% The source when this code is written is:
%%%% Yuanjin, Zheng. "A novel channel estimation and tracking method for
%%%% wireless OFDM systems based on pilots and Kalman filtering." IEEE
%%%% Transactions on Consumer Electronics 49.2 (2003): 275-283.

%% Initializing

clear all; clc; close all;

nFFT = 256; % fft size
nSym = 1*10^3; % number of symbols
CP = 16; % Cyclic Prefic
M=2; % 2,4 for BPSK and QPSK
EbN0dB =  0:5:30; % bit to noise ratio
pr = 1/4; % pilot ratio
nBitPerSym = (1-pr)*nFFT*log2(M); % number of bits per symbol
pilots = repmat((2.*(rand(1,nFFT*pr) > 0.5)-1), nSym, 1); % pilot bits generating
pInd = 1:(1/pr):nFFT; % index for pilot bits
dataInd = sort([2:(1/pr):nFFT 3:(1/pr):nFFT 4:(1/pr):nFFT]); % index for data bits

for ii = 1:numel(EbN0dB)

    %% Transmitter

    ipBit = 1.*(rand(1,log2(M)*nBitPerSym*nSym) > 0.5); % random 1's and 0's for data
    ipMod = qammod(ipBit.',M,'InputType','bit','UnitAveragePower',true); % BPSK modulation 0 --> -1, 1 --> +1
    ipMod = reshape(ipMod,nBitPerSym,nSym).'; % grouping into multiple symbols
    xF(:,dataInd) = ipMod; % creating transmitted signal in frequency domain
    xF(:,pInd) = pilots;  % creating transmitted signal in frequency domain
    % Taking FFT, the term (nFFT/sqrt(nDSC)) is for normalizing the power of transmit symbol to 1
    xt1 = (nFFT/sqrt(nFFT))*ifft(xF.').'; % creating transmitted signal in time domain
    % Appending cylic prefix
    xt = [xt1(:,[nFFT-CP+1:nFFT]) xt1]; % creating transmitted signal in time domain

    %% Channel

    % creating channel
    nTap = 10; % number of taps
    ht = 1/sqrt(2)*1/sqrt(nTap)*(randn(nSym,nTap) + 1i*randn(nSym,nTap)); % channel in time domain
    hF = (fft(ht,nFFT,2)); % channel in frequency domain
    % creating the signal after going through channel
    xht=zeros(nSym,nTap+nFFT+CP-1);
    for jj = 1:nSym
        xht(jj,:) = conv(ht(jj,:),xt(jj,:));
    end
    xt = xht;
    % Concatenating multiple symbols to form a long vector
    xt = reshape(xt.',1,nSym*((nFFT+CP)+nTap-1)); % total signal in time domain
    % Gaussian noise
    snr = EbN0dB(ii); % determining SNR
    sigma = sqrt(1/(2*(10^(snr/10)))); % Eb/N0
    % total tx power/total # of inf. bits   Eb=numel(xt)/numel(bitler)
    nt = (sigma).*[randn(1,numel(xt)) + 1i*randn(1,numel(xt))];
    % Adding noise,
    yt = xt + nt;

    %% Receiver

    yt = reshape(yt.',((nFFT+CP)+nTap-1),nSym).'; % formatting the received vector into symbols
    yt = yt(:,[CP+1:(nFFT+CP)]); % removing cyclic prefix
    yF = (sqrt(nFFT)/nFFT)*(fft(yt.')).'; % converting to frequency domain
    ipBitR = []; % vector for determined bits
    hF_ofdm = zeros(nSym,nFFT); % estimated channel for OFDM signal
    hF_kalman = zeros(nSym,nFFT);
    hF_est = zeros(1, nFFT) + 1 + 1i;
    P_est = 1;
    Q = 1;
    R = 0.001;
    for kk=1:nSym
        yFR =  yF(kk,:); % received signal for kkth symbol in freq domain
        yF_est_pil = spline(pInd, yFR(pInd)./pilots(1,:), [1:nFFT]); % estimated channel
        ht_est = (nFFT/sqrt(nFFT))*ifft(yF_est_pil.').'; % DFT based channel estimation step 1
        hF_LS = (sqrt(nFFT)/nFFT)*(fft(ht_est(:,1:CP).', nFFT)).'; % DFT based channel estimation step 2
        P_est = P_est + Q;
        Kalman_Gain = P_est / (P_est + R);
        hF_est = hF_est + Kalman_Gain * (hF_LS - hF_est);
        P_est = (1 - Kalman_Gain) * P_est;
        yF_eq = yFR./hF_est; % equalization
        yMod = yF_eq(1,dataInd); % selecting data bits
        ipBitHat = qamdemod(yMod.',M,'OutputType','bit','UnitAveragePower',true).'; % demodulation
        ipBitR = [ipBitR ipBitHat]; % storing data bits for BER calculation
        hF_kalman(kk,:) = hF_est;
    end

    %% Analyzing Part

    mse_ofdm(ii) = mean(mean(abs(hF_ofdm-hF).^2)); % MSE caclulation for channel estimation
    mse_kalman(ii)= mean(mean(abs(hF_kalman-hF).^2)); 
    nErr(ii) = size(find(ipBitR - ipBit),2); % number of errors of estimated bits
end

simBer = nErr/(nSym*nBitPerSym); % Simulated Bit Error Rate 

%% Plotting

figure()
semilogy(EbN0dB,mse_kalman,'o-','LineWidth',2);
axis([0 30 10^-6 1])
grid on
hold on
xlabel('SNR (dB)')
ylabel('MSE')
title('Mean Square Error of channel estimation for BPSK using OFDM')

figure()
semilogy(EbN0dB,simBer,'o-','LineWidth',2);
axis([0 30 10^-4 1])
grid on
hold on
xlabel('SNR (dB)')
ylabel('BER')
title('Bit error probability curve for BPSK using OFDM')











