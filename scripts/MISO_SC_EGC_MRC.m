clear; clc; close all;

num_bits = 1000000;                 
bits = rand(1, num_bits) > 0.5;    
Bits1 = bits(1:2:end);
Bits2 = bits(2:2:end);
qpsk_sig = ((Bits1==0).* (Bits2==0)*(exp(1i*pi/4)) + (Bits1==0).*(Bits2==1)*(exp(3*1i*pi/4)) + ...
    (Bits1==1).*(Bits2==1)*(exp(5*1i*pi/4)) + (Bits1==1).*(Bits2==0)*(exp(7*1i*pi/4)));

% SISO BER
h0 = sqrt(0.5)*(randn(1,length(qpsk_sig)) + 1i*randn(1,length(qpsk_sig)));
for SNR_dB = 0:2:20
    N0 = 1/10^(SNR_dB/10);                        
    a = length(qpsk_sig);
    Noise = sqrt(N0/2)*(randn(1,a) + 1i*randn(1,a)); 
    Rx_qpsk = qpsk_sig .* h0 + Noise;
    Rx_qpsk_sig = Rx_qpsk ./ h0;
    Bits4 = (real(Rx_qpsk_sig) < 0);
    Bits3 = (imag(Rx_qpsk_sig) < 0);
    Demod_qpsk_bits = zeros(1, 2*length(Rx_qpsk_sig));
    Demod_qpsk_bits(1:2:end) = Bits3;
    Demod_qpsk_bits(2:2:end) = Bits4;

    Error_bits_qpsk = bits - Demod_qpsk_bits;
    BER_qpsk(SNR_dB/2+1) = sum(abs(Error_bits_qpsk)) / num_bits;
end
semilogy((0:2:20), BER_qpsk, 'o--'); hold on;

% Rayleigh
h1 = sqrt(0.5)*(randn(1,length(qpsk_sig)) + 1i*randn(1,length(qpsk_sig)));
h2 = sqrt(0.5)*(randn(1,length(qpsk_sig)) + 1i*randn(1,length(qpsk_sig)));
h3 = sqrt(0.5)*(randn(1,length(qpsk_sig)) + 1i*randn(1,length(qpsk_sig)));
h4 = sqrt(0.5)*(randn(1,length(qpsk_sig)) + 1i*randn(1,length(qpsk_sig)));
h = {h1, h2, h3, h4};

% SC
for SNR_dB = 0:2:20
    N0 = 1/10^(SNR_dB/10);                
    a = length(qpsk_sig);
    Noise1 = sqrt(N0/2)*(randn(1,a) + 1i*randn(1,a));
    Noise2 = sqrt(N0/2)*(randn(1,a) + 1i*randn(1,a));
    Noise3 = sqrt(N0/2)*(randn(1,a) + 1i*randn(1,a));
    Noise4 = sqrt(N0/2)*(randn(1,a) + 1i*randn(1,a));
    Noise = {Noise1, Noise2, Noise3, Noise4};

    x1 = h1 .* qpsk_sig + Noise1;
    x2 = h2 .* qpsk_sig + Noise2;
    x3 = h3 .* qpsk_sig + Noise3;
    x4 = h4 .* qpsk_sig + Noise4;

    y = zeros(1,length(qpsk_sig));
    equ_y = zeros(1,length(qpsk_sig));
    h_square = zeros(1,4);
    for i = 1:length(qpsk_sig)
        h_square = [abs(h1(i))^2, abs(h2(i))^2, abs(h3(i))^2, abs(h4(i))^2];
        [~, h_i] = max(h_square);
        y(i) = h{h_i}(i) * qpsk_sig(i) + Noise{h_i}(i);
        equ_y(i) = y(i) / h{h_i}(i);
    end
    mean_power_1=mean(abs(equ_y).^2);

    Bits4 = (real(equ_y) < 0);
    Bits3 = (imag(equ_y) < 0);
    Demod_qpsk_bits = zeros(1, 2*length(equ_y));
    Demod_qpsk_bits(1:2:end) = Bits3;
    Demod_qpsk_bits(2:2:end) = Bits4;

    Error_bits_qpsk = bits - Demod_qpsk_bits;
    BER_SC(SNR_dB/2+1) = sum(abs(Error_bits_qpsk)) / num_bits;
end
semilogy((0:2:20), BER_SC, 'o--'); hold on;


%% EGC

for SNR_dB = 0:2:20
    N0 = 1/10^(SNR_dB/10);               
    a = length(qpsk_sig);
    Noise1 = sqrt(N0/2)*(randn(1,a) + 1i*randn(1,a));
    Noise2 = sqrt(N0/2)*(randn(1,a) + 1i*randn(1,a));
    Noise3 = sqrt(N0/2)*(randn(1,a) + 1i*randn(1,a));
    Noise4 = sqrt(N0/2)*(randn(1,a) + 1i*randn(1,a));
    Noise = {Noise1, Noise2, Noise3, Noise4};

    for i=1:length(qpsk_sig)
        x1(i)=conj(h1(i))*qpsk_sig(i)/(abs(h1(i)).^1);
        x1(i) = x1(i)/(abs(h1(i))+abs(h2(i))+abs(h3(i))+abs(h4(i)));
        x2(i)=conj(h2(i))*qpsk_sig(i)/(abs(h2(i)).^1);
        x2(i) = x2(i)/(abs(h1(i))+abs(h2(i))+abs(h3(i))+abs(h4(i)));
        x3(i)=conj(h3(i))*qpsk_sig(i)/(abs(h3(i)).^1);
        x3(i) = x3(i)/(abs(h1(i))+abs(h2(i))+abs(h3(i))+abs(h4(i)));
        x4(i)=conj(h4(i))*qpsk_sig(i)/(abs(h4(i)).^1);
        x4(i) = x4(i)/(abs(h1(i))+abs(h2(i))+abs(h3(i))+abs(h4(i)));
        
    end

    combined_signal = h1.*x1 + h2.*x2 + h3.*x3 + h4.*x4 + Noise1;
    combined_signal_2 = combined_signal;
    equ_y = combined_signal;

    mean_power_2=mean(abs(equ_y).^2);
    mean_noise_2=mean(abs(Noise1).^2);

    Bits4 = (real(equ_y) < 0);
    Bits3 = (imag(equ_y) < 0);
    Demod_qpsk_bits = zeros(1, 2*length(equ_y));
    Demod_qpsk_bits(1:2:end) = Bits3;
    Demod_qpsk_bits(2:2:end) = Bits4;

    Error_bits_qpsk = bits - Demod_qpsk_bits;
    BER_EGC(SNR_dB/2+1) = sum(abs(Error_bits_qpsk)) / num_bits;
end
semilogy((0:2:20), BER_EGC, 's--'); hold on;

% MRC

for SNR_dB = 0:2:20
    N0 = 1/10^(SNR_dB/10);                
    a = length(qpsk_sig);
    Noise1 = sqrt(N0/2)*(randn(1,a) + 1i*randn(1,a));
    Noise2 = sqrt(N0/2)*(randn(1,a) + 1i*randn(1,a));
    Noise3 = sqrt(N0/2)*(randn(1,a) + 1i*randn(1,a));
    Noise4 = sqrt(N0/2)*(randn(1,a) + 1i*randn(1,a));
    Noise = {Noise1, Noise2, Noise3, Noise4};

    x1 = conj(h1) .* qpsk_sig./(abs(h1).^2+abs(h2).^2+abs(h3).^2+abs(h4).^2);
    x2 = conj(h2) .* qpsk_sig./(abs(h1).^2+abs(h2).^2+abs(h3).^2+abs(h4).^2);
    x3 = conj(h3) .* qpsk_sig./(abs(h1).^2+abs(h2).^2+abs(h3).^2+abs(h4).^2);
    x4 = conj(h4) .* qpsk_sig./(abs(h1).^2+abs(h2).^2+abs(h3).^2+abs(h4).^2);

    combined_signal = h1.*x1 + h2.*x2 + h3.*x3 + h4.*x4 + Noise1;

    equ_y = combined_signal;
    mean_power_3=mean(abs(equ_y).^2);
    Bits4 = (real(equ_y) < 0);
    Bits3 = (imag(equ_y) < 0);
    Demod_qpsk_bits = zeros(1, 2*length(equ_y));
    Demod_qpsk_bits(1:2:end) = Bits3;
    Demod_qpsk_bits(2:2:end) = Bits4;

    Error_bits_qpsk = bits - Demod_qpsk_bits;
    BER_MRC(SNR_dB/2+1) = sum(abs(Error_bits_qpsk)) / num_bits;
end
semilogy((0:2:20), BER_MRC, 'x--'); hold off;

legend('SISO', 'SC', 'EGC', 'MRC');
xlabel('SNR (dB)');
ylabel('BER');
title('BER Performance of MISO');
grid on;
