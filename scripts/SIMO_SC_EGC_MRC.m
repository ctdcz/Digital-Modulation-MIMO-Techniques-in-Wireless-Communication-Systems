clear;clc;close all;
%%
% generate the QPSK signals

num_bits = 1000000;                 % Number of bits in a frame
bits = rand(1, num_bits) > 0.5;   % Generated information bits
Bits1 = bits(1:2:end);
Bits2 = bits(2:2:end);
qpsk_sig =((Bits1==0).* (Bits2==0)*(exp(1i*pi/4))+(Bits1==0).*(Bits2==1)...
    *(exp(3*1i*pi/4))+(Bits1==1).*(Bits2==1)*(exp(5*1i*pi/4))...
    +(Bits1==1).*(Bits2==0)*(exp(7*1i*pi/4)));

%%
% SISO BER
h0 = sqrt(0.5)*(randn(1,length(qpsk_sig))+i*randn(1,length(qpsk_sig)));
for SNR_dB=0:2:20
N0 = 1/10^(SNR_dB/10);                        % Noise variance
a = length(qpsk_sig);
Noise = sqrt(N0/2)*(randn(1,a)+i*randn(1,a)); % AWGN Noise
Rx_qpsk = qpsk_sig.*h0 + Noise; 
Rx_qpsk_sig = Rx_qpsk./h0;
Bits4 = (real(Rx_qpsk_sig)<0);
Bits3 = (imag(Rx_qpsk_sig)<0);
Demod_qpsk_bits = zeros(1,2*length(Rx_qpsk_sig));
Demod_qpsk_bits(1:2:end) = Bits3;
Demod_qpsk_bits(2:2:end) = Bits4;

Error_bits_qpsk = bits - Demod_qpsk_bits;
BER_qpsk(SNR_dB/2+1) = sum(abs(Error_bits_qpsk))/num_bits;
end
semilogy(( 0: 2: 20),BER_qpsk,'o--'); hold on;

%%

% generate four different the rayleigh channel

h1=sqrt(0.5)*(randn(1,length(qpsk_sig))+i*randn(1,length(qpsk_sig)));
h2=sqrt(0.5)*(randn(1,length(qpsk_sig))+i*randn(1,length(qpsk_sig)));
h3=sqrt(0.5)*(randn(1,length(qpsk_sig))+i*randn(1,length(qpsk_sig)));
h4=sqrt(0.5)*(randn(1,length(qpsk_sig))+i*randn(1,length(qpsk_sig)));
h={h1, h2, h3, h4};

%% selective combining

% for different AWGN

for SNR_dB=0:2:20
    N0 = 1/10^(SNR_dB/10);					% Noise variance
    a = length(qpsk_sig);
    Noise1 = sqrt(N0/2)*(randn(1,a)+1i*randn(1,a));
    Noise2 = sqrt(N0/2)*(randn(1,a)+1i*randn(1,a));
    Noise3 = sqrt(N0/2)*(randn(1,a)+1i*randn(1,a));
    Noise4 = sqrt(N0/2)*(randn(1,a)+1i*randn(1,a));
    Noise = {Noise1, Noise2, Noise3, Noise4};

    % four different received signal

    y1=h1.*qpsk_sig + Noise1;
    y2=h2.*qpsk_sig + Noise2;
    y3=h3.*qpsk_sig + Noise3;
    y4=h4.*qpsk_sig + Noise4;

    % select the maximum h
    y=zeros(1,length(qpsk_sig));
    equ_y=zeros(1,length(qpsk_sig));
    for i=1:length(qpsk_sig)
        h_square=[(h1(i))^2, (h2(i))^2, (h3(i))^2, (h4(i))^2];
        [h_s,h_i]=max(h_square);
        y(i) = h{h_i}(i)*qpsk_sig(i) + Noise{h_i}(i);
        equ_y(i) = y(i)./h{h_i}(i);
    end
    % demodulate the signal

    Bits4 = (real(equ_y)<0);
    Bits3 = (imag(equ_y)<0);
    Demod_qpsk_bits = zeros(1,2*length(equ_y));
    Demod_qpsk_bits(1:2:end) = Bits3;
    Demod_qpsk_bits(2:2:end) = Bits4;

    Error_bits_qpsk = bits - Demod_qpsk_bits;
    BER_qpsk(SNR_dB/2+1) = sum(abs(Error_bits_qpsk))/num_bits;
end
semilogy(( 0: 2: 20),BER_qpsk,'o--');
hold on;

%% EGC(equal gain diversity)

% for different AWGN

for SNR_dB=0:2:20
    N0 = 1/10^(SNR_dB/10);					% Noise variance
    a = length(qpsk_sig);
    Noise1 = sqrt(N0/2)*(randn(1,a)+1i*randn(1,a));
    Noise2 = sqrt(N0/2)*(randn(1,a)+1i*randn(1,a));
    Noise3 = sqrt(N0/2)*(randn(1,a)+1i*randn(1,a));
    Noise4 = sqrt(N0/2)*(randn(1,a)+1i*randn(1,a));
    Noise = {Noise1, Noise2, Noise3, Noise4};

    % four different received signal

    y1=h1.*qpsk_sig + Noise1;
    y2=h2.*qpsk_sig + Noise2;
    y3=h3.*qpsk_sig + Noise3;
    y4=h4.*qpsk_sig + Noise4;

    % select the maximum h
    y=zeros(1,length(qpsk_sig));
    equ_y=zeros(1,length(qpsk_sig));
    for i=1:length(qpsk_sig)
        yy1(i)=conj(h1(i))*y1(i)/abs(h1(i));
        yy2(i)=conj(h2(i))*y2(i)/abs(h2(i));
        yy3(i)=conj(h3(i))*y3(i)/abs(h3(i));
        yy4(i)=conj(h4(i))*y4(i)/abs(h4(i));
        equ_y(i) = (yy1(i)+yy2(i)+yy3(i)+yy4(i))/...
            (abs(h1(i))+abs(h2(i))+abs(h3(i))+abs(h4(i)));
    end
    % demodulate the signal

    Bits4 = (real(equ_y)<0);
    Bits3 = (imag(equ_y)<0);
    Demod_qpsk_bits = zeros(1,2*length(equ_y));
    Demod_qpsk_bits(1:2:end) = Bits3;
    Demod_qpsk_bits(2:2:end) = Bits4;

    Error_bits_qpsk = bits - Demod_qpsk_bits;
    BER_qpsk(SNR_dB/2+1) = sum(abs(Error_bits_qpsk))/num_bits;
end
semilogy(( 0: 2: 20),BER_qpsk,'o--');
hold on;

%% MRC(maximum ratio diversity)

% for different AWGN

for SNR_dB=0:2:20
    N0 = 1/10^(SNR_dB/10);					% Noise variance
    a = length(qpsk_sig);
    Noise1 = sqrt(N0/2)*(randn(1,a)+1i*randn(1,a));
    Noise2 = sqrt(N0/2)*(randn(1,a)+1i*randn(1,a));
    Noise3 = sqrt(N0/2)*(randn(1,a)+1i*randn(1,a));
    Noise4 = sqrt(N0/2)*(randn(1,a)+1i*randn(1,a));
    Noise = {Noise1, Noise2, Noise3, Noise4};

    % four different received signal

    y1=h1.*qpsk_sig + Noise1;
    y2=h2.*qpsk_sig + Noise2;
    y3=h3.*qpsk_sig + Noise3;
    y4=h4.*qpsk_sig + Noise4;

    % select the maximum h
    y=zeros(1,length(qpsk_sig));
    equ_y=zeros(1,length(qpsk_sig));
    for i=1:length(qpsk_sig)
        yy1(i)=conj(h1(i))*y1(i);
        yy2(i)=conj(h2(i))*y2(i);
        yy3(i)=conj(h3(i))*y3(i);
        yy4(i)=conj(h4(i))*y4(i);
        equ_y(i) = (yy1(i)+yy2(i)+yy3(i)+yy4(i))/(abs(h1(i))...
            ^2+abs(h2(i))^2+abs(h3(i))^2+abs(h4(i))^2);
    end
    % demodulate the signal

    Bits4 = (real(equ_y)<0);
    Bits3 = (imag(equ_y)<0);
    Demod_qpsk_bits = zeros(1,2*length(equ_y));
    Demod_qpsk_bits(1:2:end) = Bits3;
    Demod_qpsk_bits(2:2:end) = Bits4;

    Error_bits_qpsk = bits - Demod_qpsk_bits;
    BER_qpsk(SNR_dB/2+1) = sum(abs(Error_bits_qpsk))/num_bits;
end
semilogy(( 0: 2: 20),BER_qpsk,'o--');
xlabel("SNR(dB)"); ylabel("BER");
title("Diversity techniques");
legend("SISO","SC","EGC","MRC");
