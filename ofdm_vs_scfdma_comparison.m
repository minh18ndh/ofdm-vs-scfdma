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
% WITH MULTIPATH FADING + ZERO FORCING EQUALIZATION
%===========================================================

%---------------- SYSTEM PARAMETERS ------------------------%
Nfft   = 64;          % FFT size
Nused  = 16;          % Active subcarriers
cpLen  = 16;          % Cyclic prefix length

M      = 16;          % 16-QAM
k      = log2(M);

numBlk = 20000;
SNRdB  = 0:1:20;

BER_ofdm   = zeros(size(SNRdB));
BER_scfdma = zeros(size(SNRdB));

PAPR_ofdm   = zeros(numBlk,1);
PAPR_scfdma = zeros(numBlk,1);

%===========================================================
% MULTIPATH CHANNEL
%===========================================================
% Example 3-tap Rayleigh fading channel

L = 3;

powerProfile = [0 -3 -6];      % Tap powers in dB
powerLinear  = 10.^(powerProfile/10);

%===========================================================
% TRANSMIT DATA
%===========================================================

numBits = numBlk * Nused * k;

txBits = randi([0 1], numBits, 1);

bitMat = reshape(txBits, k, []).';
symIdx = bi2de(bitMat, 'left-msb');

txQAM = qammod(symIdx, M, 'UnitAveragePower', true);

txBlocks = reshape(txQAM, Nused, numBlk);

%=======================================================
% CHANNEL GENERATION
%=======================================================

h = (randn(L,1) + 1j*randn(L,1));

h = h .* sqrt(powerLinear(:)/2);

h = h / norm(h);

H = fft(h, Nfft);

fprintf('\nChannel taps:\n');
disp(h.');

fprintf('Total channel power = %.4f\n\n', ...
    sum(abs(h).^2));

%===========================================================
% MAIN SNR LOOP
%===========================================================

for ii = 1:length(SNRdB)

    %-------------------------------------------------------
    % PREALLOCATE
    %-------------------------------------------------------

    txSigOFDM = zeros(numBlk*(Nfft+cpLen),1);
    txSigSC   = zeros(numBlk*(Nfft+cpLen),1);

    %=======================================================
    % TRANSMITTER
    %=======================================================

    for n = 1:numBlk

        if mod(n,10000)==0
            fprintf('SNR = %2d dB | Processing block %d / %d\n', ...
                SNRdB(ii), n, numBlk);
        end

        idxStart = (n-1)*(Nfft+cpLen) + 1;
        idxEnd   = n*(Nfft+cpLen);

        data = txBlocks(:,n);

        %================ OFDM ==============================

        gridOFDM = zeros(Nfft,1);
        gridOFDM(1:Nused) = data;

        x_ofdm = ifft(gridOFDM, Nfft);

        cp_ofdm = [x_ofdm(end-cpLen+1:end); x_ofdm];

        txSigOFDM(idxStart:idxEnd) = cp_ofdm;

        % PAPR
        if ii == 1
            papr_now = max(abs(cp_ofdm).^2) / ...
                       mean(abs(cp_ofdm).^2);

            PAPR_ofdm(n) = 10*log10(papr_now);
        end

        %================ SC-FDMA ===========================

        spread = fft(data, Nused);

        gridSC = zeros(Nfft,1);
        gridSC(1:Nused) = spread;

        x_sc = ifft(gridSC, Nfft);

        cp_sc = [x_sc(end-cpLen+1:end); x_sc];

        txSigSC(idxStart:idxEnd) = cp_sc;

        if ii == 1
            papr_now = max(abs(cp_sc).^2) / ...
                       mean(abs(cp_sc).^2);

            PAPR_scfdma(n) = 10*log10(papr_now);
        end

    end

    %=======================================================
    % CHANNEL: MULTIPATH CONVOLUTION
    %=======================================================

    rxOFDM = conv(txSigOFDM, h);
    rxSC   = conv(txSigSC,   h);

    % Remove extra tail
    rxOFDM = rxOFDM(1:length(txSigOFDM));
    rxSC   = rxSC(1:length(txSigSC));

    %=======================================================
    % ADD AWGN
    %=======================================================

    rxOFDM = awgn(rxOFDM, SNRdB(ii), 'measured');
    rxSC   = awgn(rxSC,   SNRdB(ii), 'measured');

    %=======================================================
    % RECEIVER OFDM
    %=======================================================

    rxMat = reshape(rxOFDM, Nfft+cpLen, []);

    rxBitsAll = zeros(numBits,1);

    bitPtr = 1;

    for n = 1:numBlk

        r = rxMat(:,n);

        % Remove CP
        r = r(cpLen+1:end);

        % FFT
        R = fft(r, Nfft);

        %===================================================
        % ZERO FORCING EQUALIZATION
        %===================================================

        Yeq = R ./ H;

        % Extract active subcarriers
        y = Yeq(1:Nused);

        % Demodulate
        idx = qamdemod(y, M, ...
              'UnitAveragePower', true);

        bits = de2bi(idx, k, 'left-msb');
        bits = reshape(bits.', [], 1);

        rxBitsAll(bitPtr:bitPtr+length(bits)-1) = bits;

        bitPtr = bitPtr + length(bits);

    end

    BER_ofdm(ii) = mean(rxBitsAll ~= txBits);

    %=======================================================
    % RECEIVER SC-FDMA
    %=======================================================

    rxMat = reshape(rxSC, Nfft+cpLen, []);

    rxBitsAll = zeros(numBits,1);

    bitPtr = 1;

    for n = 1:numBlk

        r = rxMat(:,n);

        % Remove CP
        r = r(cpLen+1:end);

        % FFT
        R = fft(r, Nfft);

        %===================================================
        % ZERO FORCING EQUALIZATION
        %===================================================

        Yeq = R ./ H;

        % Active subcarriers
        y = Yeq(1:Nused);

        % Despread
        despread = ifft(y, Nused);

        % QAM demod
        idx = qamdemod(despread, M, ...
              'UnitAveragePower', true);

        bits = de2bi(idx, k, 'left-msb');
        bits = reshape(bits.', [], 1);

        rxBitsAll(bitPtr:bitPtr+length(bits)-1) = bits;

        bitPtr = bitPtr + length(bits);

    end

    BER_scfdma(ii) = mean(rxBitsAll ~= txBits);

end

%===========================================================
% BER TABLE
%===========================================================

fprintf('\n============== BER TABLE ==============\n');
fprintf(' SNR(dB) |    OFDM BER    |  SC-FDMA BER\n');
fprintf('---------------------------------------\n');

for i = 1:length(SNRdB)

    fprintf(' %7d | %14.6e | %14.6e\n', ...
        SNRdB(i), ...
        BER_ofdm(i), ...
        BER_scfdma(i));

end

fprintf('=======================================\n\n');

%===========================================================
% BER PLOT
%===========================================================

figure;
semilogy(SNRdB, BER_ofdm, 'o-','LineWidth',1.5);
hold on;

semilogy(SNRdB, BER_scfdma, 's-','LineWidth',1.5);

grid on;

xlabel('SNR (dB)');
ylabel('BER');

title('BER Comparison with Multipath + ZF Equalization');

legend('OFDM','SC-FDMA');

%===========================================================
% PAPR CCDF
%===========================================================

paprAxis = 0:0.2:12;

ccdf_ofdm = zeros(size(paprAxis));
ccdf_sc   = zeros(size(paprAxis));

for i = 1:length(paprAxis)

    ccdf_ofdm(i) = mean(PAPR_ofdm > paprAxis(i));

    ccdf_sc(i) = mean(PAPR_scfdma > paprAxis(i));

end

fprintf('\n=========== PAPR CCDF ===========\n');
fprintf('Threshold | OFDM CCDF | SC-FDMA CCDF\n');

for i = 1:length(paprAxis)

    fprintf('  %5.1f dB | %10.6f | %10.6f\n', ...
        paprAxis(i), ...
        ccdf_ofdm(i), ...
        ccdf_sc(i));

end

fprintf('=================================\n\n');

figure;

semilogy(paprAxis, ccdf_ofdm, 'o-','LineWidth',1.5);
hold on;

semilogy(paprAxis, ccdf_sc, 's-','LineWidth',1.5);

grid on;

xlabel('PAPR Threshold (dB)');
ylabel('CCDF');

title('PAPR Comparison');

legend('OFDM','SC-FDMA');

%===========================================================
% DISPLAY RESULTS
%===========================================================

disp('Average PAPR OFDM (dB):');
disp(mean(PAPR_ofdm));

disp('Average PAPR SC-FDMA (dB):');
disp(mean(PAPR_scfdma));
