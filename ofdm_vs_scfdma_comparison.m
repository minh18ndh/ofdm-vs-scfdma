clc;
clear;
close all;

%===========================================================
% OFDM vs SC-FDMA COMPARISON
% Same FFT size
% Same active subcarriers
% Same modulation
% Same CP
% Same transmitted bits
% Same SNR
%===========================================================

%---------------- SYSTEM PARAMETERS ------------------------%
Nfft   = 64;          % FFT size
Nused  = 16;          % Active subcarriers for BOTH systems
cpLen  = 16;          % Cyclic prefix
M      = 16;          % 16-QAM
k      = log2(M);     % Bits per symbol

numBlk = 4000;        % Number of transmitted blocks
SNRdB  = 0:2:20;

BER_ofdm   = zeros(size(SNRdB));
BER_scfdma = zeros(size(SNRdB));

PAPR_ofdm   = [];
PAPR_scfdma = [];

%===========================================================
% TRANSMIT DATA (same data for both systems)
%===========================================================
numBits = numBlk * Nused * k;

txBits = randi([0 1], numBits, 1);

bitMat = reshape(txBits, k, []).';
symIdx = bi2de(bitMat, 'left-msb');

txQAM = qammod(symIdx, M, 'UnitAveragePower', true);

txBlocks = reshape(txQAM, Nused, numBlk);

%===========================================================
% MAIN LOOP OVER SNR
%===========================================================
for ii = 1:length(SNRdB)

    txSigOFDM = [];
    txSigSC   = [];

    %-------------------------------------------------------
    % TRANSMITTER
    %-------------------------------------------------------
    for n = 1:numBlk

        data = txBlocks(:,n);

        %================ OFDM ==============================
        gridOFDM = zeros(Nfft,1);
        gridOFDM(1:Nused) = data;

        x_ofdm = ifft(gridOFDM, Nfft);

        cp_ofdm = [x_ofdm(end-cpLen+1:end); x_ofdm];

        txSigOFDM = [txSigOFDM; cp_ofdm];

        % PAPR save once
        if ii == 1
            papr_now = max(abs(cp_ofdm).^2) / mean(abs(cp_ofdm).^2);
            PAPR_ofdm = [PAPR_ofdm; 10*log10(papr_now)];
        end

        %================ SC-FDMA ===========================
        spread = fft(data, Nused);

        gridSC = zeros(Nfft,1);
        gridSC(1:Nused) = spread;

        x_sc = ifft(gridSC, Nfft);

        cp_sc = [x_sc(end-cpLen+1:end); x_sc];

        txSigSC = [txSigSC; cp_sc];

        if ii == 1
            papr_now = max(abs(cp_sc).^2) / mean(abs(cp_sc).^2);
            PAPR_scfdma = [PAPR_scfdma; 10*log10(papr_now)];
        end

    end

    %-------------------------------------------------------
    % CHANNEL
    %-------------------------------------------------------
    rxOFDM = awgn(txSigOFDM, SNRdB(ii), 'measured');
    rxSC   = awgn(txSigSC,   SNRdB(ii), 'measured');

    %-------------------------------------------------------
    % RECEIVER OFDM
    %-------------------------------------------------------
    rxMat = reshape(rxOFDM, Nfft+cpLen, []);

    rxBitsAll = [];

    for n = 1:numBlk

        r = rxMat(:,n);
        r = r(cpLen+1:end);

        R = fft(r, Nfft);

        y = R(1:Nused);

        idx = qamdemod(y, M, 'UnitAveragePower', true);

        bits = de2bi(idx, k, 'left-msb');
        bits = reshape(bits.', [], 1);

        rxBitsAll = [rxBitsAll; bits];
    end

    BER_ofdm(ii) = sum(rxBitsAll ~= txBits) / numBits;

    %-------------------------------------------------------
    % RECEIVER SC-FDMA
    %-------------------------------------------------------
    rxMat = reshape(rxSC, Nfft+cpLen, []);

    rxBitsAll = [];

    for n = 1:numBlk

        r = rxMat(:,n);
        r = r(cpLen+1:end);

        R = fft(r, Nfft);

        y = R(1:Nused);

        despread = ifft(y, Nused);

        idx = qamdemod(despread, M, 'UnitAveragePower', true);

        bits = de2bi(idx, k, 'left-msb');
        bits = reshape(bits.', [], 1);

        rxBitsAll = [rxBitsAll; bits];
    end

    BER_scfdma(ii) = sum(rxBitsAll ~= txBits) / numBits;

end

%===========================================================
% BER PLOT
%===========================================================
figure;
semilogy(SNRdB, BER_ofdm, 'o-','LineWidth',1.5); hold on;
semilogy(SNRdB, BER_scfdma, 's-','LineWidth',1.5);
grid on;
xlabel('SNR (dB)');
ylabel('BER');
title('BER Comparison: OFDM vs SC-FDMA');
legend('OFDM','SC-FDMA');

%===========================================================
% PAPR CCDF APPROXIMATION
%===========================================================
paprAxis = 0:0.2:12;

ccdf_ofdm = zeros(size(paprAxis));
ccdf_sc   = zeros(size(paprAxis));

for i = 1:length(paprAxis)
    ccdf_ofdm(i) = mean(PAPR_ofdm > paprAxis(i));
    ccdf_sc(i)   = mean(PAPR_scfdma > paprAxis(i));
end

figure;
semilogy(paprAxis, ccdf_ofdm, 'o-','LineWidth',1.5); hold on;
semilogy(paprAxis, ccdf_sc, 's-','LineWidth',1.5);
grid on;
xlabel('PAPR Threshold (dB)');
ylabel('CCDF');
title('PAPR Comparison: OFDM vs SC-FDMA');
legend('OFDM','SC-FDMA');

%===========================================================
% DISPLAY
%===========================================================
disp('Average PAPR OFDM (dB):');
disp(mean(PAPR_ofdm));

disp('Average PAPR SC-FDMA (dB):');
disp(mean(PAPR_scfdma));
