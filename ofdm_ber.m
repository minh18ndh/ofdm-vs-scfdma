clc;
clear;
close all;

%================ PARAMETERS ================%
Nfft = 64;              % Number of subcarriers
cpLen = 16;             % Cyclic prefix length
M = 16;                 % 16-QAM
k = log2(M);            % Bits per QAM symbol

numOFDMSym = 2000;      % Number of OFDM symbols
SNRdB = 0:2:20;         % SNR range

BER = zeros(size(SNRdB));

%============== TRANSMIT DATA ==============%
numBits = numOFDMSym * Nfft * k;
txBits = randi([0 1], numBits, 1);

% Group bits into symbols
bitMatrix = reshape(txBits, k, []).';
symIndex = bi2de(bitMatrix, 'left-msb');

% 16-QAM modulation
txQAM = qammod(symIndex, M, 'UnitAveragePower', true);

% Map to OFDM grid
txGrid = reshape(txQAM, Nfft, []);

% IFFT
txTime = ifft(txGrid, Nfft, 1);

% Add cyclic prefix
txCP = [txTime(end-cpLen+1:end,:); txTime];

% Serialize
txSignal = txCP(:);

%============== MAIN LOOP ==================%
for ii = 1:length(SNRdB)

    % AWGN channel
    rxSignal = awgn(txSignal, SNRdB(ii), 'measured');

    % Parallelize
    rxMat = reshape(rxSignal, Nfft + cpLen, []);

    % Remove CP
    rxNoCP = rxMat(cpLen+1:end, :);

    % FFT
    rxGrid = fft(rxNoCP, Nfft, 1);

    % Serialize symbols
    rxQAM = rxGrid(:);

    % QAM demod
    rxIndex = qamdemod(rxQAM, M, 'UnitAveragePower', true);

    % Symbols to bits
    rxBitsMat = de2bi(rxIndex, k, 'left-msb');
    rxBits = reshape(rxBitsMat.', [], 1);

    % BER
    BER(ii) = sum(rxBits ~= txBits) / numBits;

end

%================ PLOT =====================%
figure;
semilogy(SNRdB, BER, 'o-','LineWidth',1.5);
grid on;
xlabel('SNR (dB)');
ylabel('Bit Error Rate');
title('OFDM (16-QAM) BER Performance');