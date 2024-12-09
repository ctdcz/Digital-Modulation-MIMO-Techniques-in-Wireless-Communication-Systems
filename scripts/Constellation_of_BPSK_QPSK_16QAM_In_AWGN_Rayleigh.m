% Parameters
snr_values = [5, 15]; % SNR values for plotting
num_symbols = 1000; % Number of symbols for constellation visualization

% Generate random symbols
bits_bpsk = randi([0, 1], 1, num_symbols); % BPSK bits
bits_qpsk = randi([0, 3], 1, num_symbols); % QPSK bits
bits_16qam = randi([0, 15], 1, num_symbols); % 16QAM bits

% Modulate signals
bpsk_sig = pskmod(bits_bpsk, 2); % BPSK
qpsk_sig = pskmod(bits_qpsk, 4, pi/4); % QPSK
qam16_sig = 1/sqrt(10) * qammod(bits_16qam, 16); % 16QAM

% Loop through SNR values and plot constellations
for snr_idx = 1:length(snr_values)
    snr = snr_values(snr_idx);
    noise_power = 10^(-snr / 10); % Noise power calculation

    % AWGN noise
    noise_awgn = sqrt(noise_power / 2) * (randn(1, num_symbols) + 1j * randn(1, num_symbols));

    % Rayleigh fading
    rayleigh_fading = (randn(1, num_symbols) + 1j * randn(1, num_symbols)) / sqrt(2);

    % AWGN Channel
    rx_bpsk_awgn = bpsk_sig + noise_awgn;
    rx_qpsk_awgn = qpsk_sig + noise_awgn;
    rx_qam16_awgn = qam16_sig + noise_awgn;

    % Rayleigh Channel
    rx_bpsk_rayleigh = rayleigh_fading .* bpsk_sig + noise_awgn;
    rx_qpsk_rayleigh = rayleigh_fading .* qpsk_sig + noise_awgn;
    rx_qam16_rayleigh = rayleigh_fading .* qam16_sig + noise_awgn;

    % Normalize Rayleigh Channel signals
    rx_bpsk_rayleigh_normalized = rx_bpsk_rayleigh ./ rayleigh_fading;
    rx_qpsk_rayleigh_normalized = rx_qpsk_rayleigh ./ rayleigh_fading;
    rx_qam16_rayleigh_normalized = rx_qam16_rayleigh ./ rayleigh_fading;

    % Plot constellations
    figure('Name', sprintf('Constellation Plots (SNR = %d dB)', snr));

    % BPSK
    subplot(3, 2, 1);
    scatter(real(rx_bpsk_awgn), imag(rx_bpsk_awgn), 'o');
    title(sprintf('BPSK AWGN (SNR = %ddB)', snr));
    xlabel('In-Phase'); ylabel('Quadrature'); grid on;

    subplot(3, 2, 2);
    scatter(real(rx_bpsk_rayleigh_normalized), imag(rx_bpsk_rayleigh_normalized), 'o');
    title(sprintf('BPSK Rayleigh (SNR = %ddB)', snr));
    xlabel('In-Phase'); ylabel('Quadrature'); grid on;

    % QPSK
    subplot(3, 2, 3);
    scatter(real(rx_qpsk_awgn), imag(rx_qpsk_awgn), 'o');
    title(sprintf('QPSK AWGN (SNR = %ddB)', snr));
    xlabel('In-Phase'); ylabel('Quadrature'); grid on;

    subplot(3, 2, 4);
    scatter(real(rx_qpsk_rayleigh_normalized), imag(rx_qpsk_rayleigh_normalized), 'o');
    title(sprintf('QPSK Rayleigh (SNR = %ddB)', snr));
    xlabel('In-Phase'); ylabel('Quadrature'); grid on;

    % 16QAM
    subplot(3, 2, 5);
    scatter(real(rx_qam16_awgn), imag(rx_qam16_awgn), 'o');
    title(sprintf('16QAM AWGN (SNR = %ddB)', snr));
    xlabel('In-Phase'); ylabel('Quadrature'); grid on;

    subplot(3, 2, 6);
    scatter(real(rx_qam16_rayleigh_normalized), imag(rx_qam16_rayleigh_normalized), 'o');
    title(sprintf('16QAM Rayleigh (SNR = %ddB)', snr));
    xlabel('In-Phase'); ylabel('Quadrature'); grid on;
end


