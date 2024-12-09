%% V-BLAST
clear;clc;close all;

%% ML+QPSK
num_bits = 1e6;				  % Number of bits in a frame 
bits = rand(1, num_bits) > 0.5;	  % Generated information bits
 
Bits1 = bits(1:2:end);
Bits2 = bits(2:2:end);

% normalizing factor
a0 = sqrt(1/2);	

% bit mapping
qpsk_sig = a0*((Bits1==0).*(Bits2==0)*(exp(1i*pi/4))...
         +(Bits1==0).*(Bits2==1)*(exp(3*1i*pi/4))...
         +(Bits1==1).*(Bits2==1)*(exp(5*1i*pi/4))...
         +(Bits1==1).*(Bits2==0)*(exp(7*1i*pi/4))); 
x1 = qpsk_sig(1:2:end);
x2 = qpsk_sig(2:2:end);
a = length(x1);

for SNR_dB=0:2:20
N0 = 1/10^(SNR_dB/10);					% Noise variance
n1 = sqrt(N0/2)*(randn(1,a)+1i*randn(1,a));	% AWGN Noise
n2 = sqrt(N0/2)*(randn(1,a)+1i*randn(1,a));	% AWGN Noise
h11 = sqrt(0.5)*(randn(1,a)+1i*randn(1,a)); % Rayleigh channel
h21 = sqrt(0.5)*(randn(1,a)+1i*randn(1,a)); % Rayleigh channel
h12 = sqrt(0.5)*(randn(1,a)+1i*randn(1,a)); % Rayleigh channel
h22 = sqrt(0.5)*(randn(1,a)+1i*randn(1,a)); % Rayleigh channel

%%
X_hat{1} = a0*[(exp(1i*pi/4)) (exp(1i*pi/4))]'; %0000
X_hat{2} = a0*[(exp(1i*pi/4)) (exp(3*1i*pi/4))]'; %0001
X_hat{3} = a0*[(exp(1i*pi/4)) (exp(5*1i*pi/4))]'; %0010
X_hat{4} = a0*[(exp(1i*pi/4)) (exp(7*1i*pi/4))]'; %0011
X_hat{5} = a0*[(exp(3*1i*pi/4)) (exp(1i*pi/4))]'; %0100
X_hat{6} = a0*[(exp(3*1i*pi/4)) (exp(3*1i*pi/4))]'; %0101
X_hat{7} = a0*[(exp(3*1i*pi/4)) (exp(5*1i*pi/4))]'; %0110
X_hat{8} = a0*[(exp(3*1i*pi/4)) (exp(7*1i*pi/4))]'; %0111
X_hat{9} = a0*[(exp(5*1i*pi/4)) (exp(1i*pi/4))]'; %1000
X_hat{10} = a0*[(exp(5*1i*pi/4)) (exp(3*1i*pi/4))]'; %1001
X_hat{11} = a0*[(exp(5*1i*pi/4)) (exp(5*1i*pi/4))]'; %1010
X_hat{12} = a0*[(exp(5*1i*pi/4)) (exp(7*1i*pi/4))]'; %1011
X_hat{13} = a0*[(exp(7*1i*pi/4)) (exp(1i*pi/4))]'; %1100
X_hat{14} = a0*[(exp(7*1i*pi/4)) (exp(3*1i*pi/4))]'; %1101
X_hat{15} = a0*[(exp(7*1i*pi/4)) (exp(5*1i*pi/4))]'; %1110
X_hat{16} = a0*[(exp(7*1i*pi/4)) (exp(7*1i*pi/4))]'; %1111

errors_sum = 0;

for i=1:1:a
    X = [x1(i) x2(i)]';
    N = [n1(i) n2(i)]';
    H = [h11(i),h12(i);h21(i),h22(i)];
    Y = H*X + N;

    % ML
    for j=1:16
        X_test(j) = norm(Y-H*X_hat{j},'fro')^2;
    end

    [X_t,X_hat_I]= min(X_test);
    Rx_qpsk_sig = X_hat{X_hat_I};
    % Demodulation
    Bits4 = (real(Rx_qpsk_sig)<0);
    Bits3 = (imag(Rx_qpsk_sig)>0);
    Demod_qpsk_bits = zeros(1,2*length(Rx_qpsk_sig)); % Demodulated bits
    Demod_qpsk_bits(1:2:end) = Bits3;
    Demod_qpsk_bits(2:2:end) = Bits4;
    TX_qpsk_bits = bits(4*i-3:4*i);
    Error_bits_qpsk = abs(TX_qpsk_bits - Demod_qpsk_bits);
    errors = sum(Error_bits_qpsk);
    errors_sum = errors_sum + errors;
end
BER(SNR_dB/2+1) = errors_sum/num_bits;
end
figure(1);
semilogy((0:2:20), BER);
hold on;

%% ZF+QPSK
num_bits = 1e6;				  % Number of bits in a frame 
bits = rand(1, num_bits) > 0.5;	  % Generated information bits
 
Bits1 = bits(1:2:end);
Bits2 = bits(2:2:end);

% normalizing factor
a0 = sqrt(1/2);	

% bit mapping
qpsk_sig = a0*((Bits1==0).*(Bits2==0)*(exp(1i*pi/4))...
         +(Bits1==0).*(Bits2==1)*(exp(3*1i*pi/4))...
         +(Bits1==1).*(Bits2==1)*(exp(5*1i*pi/4))...
         +(Bits1==1).*(Bits2==0)*(exp(7*1i*pi/4))); 
x1 = qpsk_sig(1:2:end);
x2 = qpsk_sig(2:2:end);
a = length(x1);



for SNR_dB=0:2:20

N0 = 1/10^(SNR_dB/10);					% Noise variance
n1 = sqrt(N0/2)*(randn(1,a)+1i*randn(1,a));	% AWGN Noise
n2 = sqrt(N0/2)*(randn(1,a)+1i*randn(1,a));	% AWGN Noise
h11 = sqrt(0.5)*(randn(1,a)+1i*randn(1,a)); % Rayleigh channel
h21 = sqrt(0.5)*(randn(1,a)+1i*randn(1,a)); % Rayleigh channel
h12 = sqrt(0.5)*(randn(1,a)+1i*randn(1,a)); % Rayleigh channel
h22 = sqrt(0.5)*(randn(1,a)+1i*randn(1,a)); % Rayleigh channel

errors_sum = 0;
for i=1:1:a
    X = [x1(i) x2(i)]';
    N = [n1(i) n2(i)]';
    H = [h11(i),h12(i);h21(i),h22(i)];
    Y = H*X + N;
    H_inverse = conj(H)'/(H*conj(H)');
    X_hat = H_inverse*Y;
    Rx_qpsk_sig = X_hat;
    
    % Demodulation
    Bits4 = (real(Rx_qpsk_sig)<0);
    Bits3 = (imag(Rx_qpsk_sig)>0);
    Demod_qpsk_bits = zeros(1,2*length(Rx_qpsk_sig)); % Demodulated bits
    Demod_qpsk_bits(1:2:end) = Bits3;
    Demod_qpsk_bits(2:2:end) = Bits4;
    TX_qpsk_bits = bits(4*i-3:4*i);
    Error_bits_qpsk = abs(TX_qpsk_bits - Demod_qpsk_bits);
    errors = sum(Error_bits_qpsk);
    errors_sum = errors_sum + errors;
end
BER(SNR_dB/2+1) = errors_sum/num_bits;
end
semilogy(( 0: 2: 20),BER,'o--');
xlabel('SNR (dB)'); ylabel('BER');
title('BER Performance of V-BLAST MIMO System');
legend("ML","ZF");