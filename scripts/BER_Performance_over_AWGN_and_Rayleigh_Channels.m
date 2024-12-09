% Simulation Parameters
SNR_dB = 0:2:20; % SNR from 0dB to 20dB in steps of 2dB
num_bits = 100000; % Number of bits for each transmission

% Pre-allocate BER arrays
BER_BPSK_AWGN = zeros(size(SNR_dB));
BER_QPSK_AWGN = zeros(size(SNR_dB));
BER_16QAM_AWGN = zeros(size(SNR_dB));
BER_2FSK_AWGN = zeros(size(SNR_dB));
BER_4FSK_AWGN = zeros(size(SNR_dB));
BER_BPSK_Rayleigh = zeros(size(SNR_dB));
BER_QPSK_Rayleigh = zeros(size(SNR_dB));
BER_16QAM_Rayleigh = zeros(size(SNR_dB));
BER_2FSK_Rayleigh = zeros(size(SNR_dB));
BER_4FSK_Rayleigh = zeros(size(SNR_dB));

% Generate random bits
% bits = randi([0 1], 1, num_bits);
bits_bpsk = randi([0,1],1,num_bits);
bits_qpsk = randi([0,3],1,num_bits);
bits_16qam = randi([0,15],1,num_bits);
% Modulate the signals
bpsk_sig = pskmod(bits_bpsk,2); % BPSK modulation: {0,1} -> {-1,+1}
qpsk_sig = pskmod(bits_qpsk, 4, pi/4); % QPSK modulation with pi/4 offset
qam16_sig = 1/sqrt(10)*qammod(bits_16qam, 16); % 16QAM modulation
% Generate random symbols for 2FSK and 4FSK
bits_2fsk = randi([0, 1], 1, num_bits); % For 2FSK
bits_4fsk = randi([0, 3], 1, num_bits); % For 4FSK

% Modulate the signals
M = 4;
freqsep = 2;
Fs = 8;
nsamp = 4;
fsk2_sig = fskmod(bits_2fsk, M, freqsep, nsamp, Fs); % 2FSK modulation
fsk4_sig = fskmod(bits_4fsk, M, freqsep, nsamp, Fs); % 4FSK modulation
for i = 1:length(SNR_dB)
    % Calculate noise standard deviation for AWGN
    SNR_linear = 10^(SNR_dB(i) / 10);
    n_power = 10.^(-SNR_dB(i)/10);
    n_real = sqrt(n_power/2) * randn(num_bits,1);
    n_imag = sqrt(n_power/2) * randn(num_bits,1);
    noise_awgn = n_real + 1j * n_imag;
    noise_awgn = noise_awgn.';
    
    n_real_2fsk = sqrt(n_power/2) * randn(num_bits*nsamp,1);
    n_imag_2fsk = sqrt(n_power/2) * randn(num_bits*nsamp,1);
    noise_awgn_2fsk = n_real_2fsk + 1j * n_imag_2fsk;
    noise_awgn_2fsk = noise_awgn_2fsk.';

    n_real_4fsk = sqrt(n_power/2) * randn(num_bits*nsamp,1);
    n_imag_4fsk = sqrt(n_power/2) * randn(num_bits*nsamp,1);
    noise_awgn_4fsk = n_real_4fsk + 1j * n_imag_4fsk;
    noise_awgn_4fsk = noise_awgn_4fsk.';
    

    % BPSK through AWGN
    rx_bpsk_awgn = bpsk_sig + noise_awgn;
    bpsk_demod_awgn = pskdemod(rx_bpsk_awgn, 2);
    BER_BPSK_AWGN(i) = sum(bits_bpsk ~= bpsk_demod_awgn) / num_bits;

    % QPSK through AWGN
    rx_qpsk_awgn = qpsk_sig + noise_awgn;
    qpsk_demod_awgn = pskdemod(rx_qpsk_awgn, 4, pi/4);
    BER_QPSK_AWGN(i) = sum(bits_qpsk ~= qpsk_demod_awgn) / num_bits;

    % 16QAM through AWGN
    rx_qam16_awgn = qam16_sig + noise_awgn;
    qam16_demod_awgn = qamdemod(sqrt(10)*rx_qam16_awgn, 16);
    BER_16QAM_AWGN(i) = sum(bits_16qam ~= qam16_demod_awgn) / num_bits;

    % 2FSK through AWGN
    rx_fsk2_awgn = fsk2_sig + noise_awgn_2fsk; % Add AWGN to 2FSK signal
    demod_bits_2fsk_awgn = fskdemod(rx_fsk2_awgn, M, freqsep, nsamp, Fs); % Demodulate 2FSK signal
    BER_2FSK_AWGN(i) = sum(bits_2fsk ~= demod_bits_2fsk_awgn) / num_bits; % Calculate BER for 2FSK
    
    % 4FSK through AWGN
    rx_fsk4_awgn = fsk4_sig + noise_awgn_4fsk; % Add AWGN to 4FSK signal
    demod_bits_4fsk_awgn = fskdemod(rx_fsk4_awgn, M, freqsep, nsamp, Fs); % Demodulate 4FSK signal
    BER_4FSK_AWGN(i) = sum(bits_4fsk ~= demod_bits_4fsk_awgn) / num_bits; % Calculate BER for 4FSK

    % Rayleigh Fading Channel
    rayleigh_fading = (randn(1, num_bits) + 1i * randn(1, num_bits)) / sqrt(2);
    rayleigh_fading_2fsk = (randn(1, num_bits*nsamp) + 1i * randn(1, num_bits*nsamp)) / sqrt(2);
    rayleigh_fading_4fsk = (randn(1, num_bits*nsamp) + 1i * randn(1, num_bits*nsamp)) / sqrt(2);
    
    % BPSK through Rayleigh channel
    rx_bpsk_rayleigh = rayleigh_fading .* bpsk_sig + noise_awgn;
    bpsk_demod_rayleigh = pskdemod(real(rx_bpsk_rayleigh ./ rayleigh_fading),2);
    BER_BPSK_Rayleigh(i) = sum(bits_bpsk ~= bpsk_demod_rayleigh) / num_bits;

    % QPSK through Rayleigh channel
    rx_qpsk_rayleigh = rayleigh_fading .* qpsk_sig + noise_awgn;
    qpsk_demod_rayleigh = pskdemod(rx_qpsk_rayleigh ./ rayleigh_fading, 4, pi/4);
    BER_QPSK_Rayleigh(i) = sum(bits_qpsk ~= qpsk_demod_rayleigh) / num_bits/2;

    % 16QAM through Rayleigh channel
    rx_qam16_rayleigh = rayleigh_fading .* qam16_sig + noise_awgn;
    qam16_demod_rayleigh = qamdemod(sqrt(10)*(rx_qam16_rayleigh ./ rayleigh_fading), 16);
    BER_16QAM_Rayleigh(i) = sum(bits_16qam ~= qam16_demod_rayleigh) / num_bits/4;
    

    % 2FSK through Rayleigh
    rx_fsk2_rayleigh = rayleigh_fading_2fsk .* fsk2_sig + noise_awgn_2fsk;
    demod_bits_2fsk_rayleigh = fskdemod(rx_fsk2_rayleigh ./ rayleigh_fading_2fsk, M, freqsep, nsamp, Fs);
    BER_2FSK_Rayleigh(i) = sum(bits_2fsk ~= demod_bits_2fsk_rayleigh) / num_bits;

    % 4FSK through Rayleigh
    rx_fsk4_rayleigh = rayleigh_fading_4fsk .* fsk4_sig + noise_awgn_4fsk;
    demod_bits_4fsk_rayleigh = fskdemod(rx_fsk4_rayleigh ./ rayleigh_fading_4fsk, M, freqsep, nsamp, Fs);
    BER_4FSK_Rayleigh(i) = sum(bits_4fsk ~= demod_bits_4fsk_rayleigh) / num_bits;


end

% Plot the BER results
figure;
semilogy(SNR_dB, BER_BPSK_AWGN, '-o', 'LineWidth', 2);
hold on;
semilogy(SNR_dB, BER_QPSK_AWGN, '-s', 'LineWidth', 2);
semilogy(SNR_dB, BER_16QAM_AWGN, '-d', 'LineWidth', 2);
semilogy(SNR_dB, BER_BPSK_Rayleigh, '--o', 'LineWidth', 2);
semilogy(SNR_dB, BER_QPSK_Rayleigh, '--s', 'LineWidth', 2);
semilogy(SNR_dB, BER_16QAM_Rayleigh, '--d', 'LineWidth', 2);

semilogy(SNR_dB, BER_2FSK_Rayleigh, '-d', 'LineWidth', 2);
semilogy(SNR_dB, BER_4FSK_Rayleigh, '-d', 'LineWidth', 2);
semilogy(SNR_dB, BER_2FSK_AWGN, '--d', 'LineWidth', 2);
semilogy(SNR_dB, BER_4FSK_AWGN, '--d', 'LineWidth', 2);
grid on;
xlabel('SNR (dB)');
ylabel('Bit Error Rate (BER)');
legend('BPSK AWGN', 'QPSK AWGN', '16QAM AWGN', ...
       'BPSK Rayleigh', 'QPSK Rayleigh', '16QAM Rayleigh','2FSK Rayleigh','4FSK Rayleigh','2FSK AWGN','4FSK AWGN');
title('BER Performance over AWGN and Rayleigh Channels');
