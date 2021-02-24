clear all;
close all;
clc;

%% Dados
SNR = 0:5:40;
BER_s0 = zeros(size(SNR));
BER_s1 = zeros(size(SNR));

NrSub = 64;
Nofdm = 64;

SNR = 0:5:40;
NrIteracoes = 10^4;

%% Transmissor

% QAM
M = 16;  % Bits por simbolo
%mod = comm.PSKModulator(M,'BitInput',true,'PhaseOffset',pi/M);
%demod = comm.PSKDemodulator(M,'BitOutput',true,'PhaseOffset',pi/M);
mod = comm.RectangularQAMModulator('ModulationOrder',16,'BitInput',true);
demod = comm.RectangularQAMDemodulator('ModulationOrder',16,'BitOutput',true);

% OFDM 48 subportadoras de informacao
Bits_info = log2(M);
s0 = randi([0 1], Bits_info*NrSub, 1);
s1 = randi([0 1], Bits_info*NrSub, 1);

% Modulando
s0_mod = step(mod, s0);
s1_mod = step(mod, s1);

% Todas as subportadoras sao de informacao
Subinfo_s0 = s0_mod;
Subinfo_s1 = s1_mod;

% Adicionando prefixo ciclico (cp)
cp = 16;
ofdm_ifft_s0 = ifft(Subinfo_s0);
ofdm_ifft_s1 = ifft(Subinfo_s1);

ofdm_s0 = zeros(Nofdm+cp, 1);
ofdm_s1 = zeros(Nofdm+cp, 1);

ofdm_s0(1:cp,1) = ofdm_ifft_s0(Nofdm-cp+1:Nofdm,1); 
ofdm_s0(cp+1:Nofdm+cp,1) = ofdm_ifft_s0;

ofdm_s1(1:cp,1) = ofdm_ifft_s1(Nofdm-cp+1:Nofdm,1); 
ofdm_s1(cp+1:Nofdm+cp,1) = ofdm_ifft_s1;

% Variacao da informacao
ofdm_s0_conj = conj(ofdm_s0);
ofdm_s1_conj_neg = conj(ofdm_s1);

for nSNR = 1:length(SNR)
    Error_s0 = 0;
    Error_s1 = 0;
    for nIteracao = 1:NrIteracoes
        % Canal Rayleigh - Tabela II
        h0 = sqrt(0.5)*(randn+1j*randn);                % RX0 - TX0
        h1 = sqrt(0.5)*(randn+1j*randn);                % RX0 - TX1
        h2 = sqrt(0.5)*(randn+1j*randn);                % RX1 - TX0
        h3 = sqrt(0.5)*(randn+1j*randn);                % RX1 - TX1

        % Sinais no receptor - Eq. 14
        R0 = h0*ofdm_s0 + h1*ofdm_s1;                  
        R1 = -h0*conj(ofdm_s1) + h1*conj(ofdm_s0);   
        R2 = h2*ofdm_s0 + h3*ofdm_s1;                                        
        R3 = -h2*conj(ofdm_s1) + h3*conj(ofdm_s0);   
        
        % Adicionando o ruido do canal
        SNR_value = SNR(nSNR);
        RO = awgn(R0,SNR_value);
        R1 = awgn(R1,SNR_value);
        R2 = awgn(R2,SNR_value);
        R3 = awgn(R3,SNR_value);
        
        % Sinais no combinador - Eq. 15
        s0_combinador = conj(h0)*R0 + h1*conj(R1) + conj(h2)*R2 + h3*conj(R3);  
        s1_combinador = conj(h1)*R0 - h0*conj(R1) + conj(h3)*R2 - h2*conj(R3);    
        
        % Retirar o prefixo ciclico
        RXofdm_s0 = s0_combinador(cp+1:Nofdm+cp, 1);
        RXofdm_s1 = s1_combinador(cp+1:Nofdm+cp, 1);

        % FFT
        RX_s0 = fft(RXofdm_s0);
        RX_s1 = fft(RXofdm_s1);
        
        % Demodulando os sinais
        s0_demod = step(demod,RX_s0); 
        s1_demod = step(demod,RX_s1);  

        % Verificando erros
        erro_s0 = length(find(s0_demod ~= s0));
        erro_s1 = length(find(s1_demod ~= s1));

        Error_s0 = Error_s0 + erro_s0;
        Error_s1 = Error_s1 + erro_s1;
    end
    BER_s0(nSNR,1) = Error_s0/NrIteracoes;
    BER_s1(nSNR,1) = Error_s1/NrIteracoes;
end

% Plotando BER
BER_s0 = BER_s0./(log2(M)*10^3);
BER_s1 = BER_s1./(log2(M)*10^3);

figure(1);
semilogy(SNR,BER_s0,'b'); 
xlabel('SNR [dB]');
ylabel('BER');
ylim([0.014 0.024])
xlim([0 20])
title('Alamouti e OFDM - 2x2 - Sinal 0')

figure(2);
semilogy(SNR,BER_s1,'r'); 
xlabel('SNR [dB]');
ylabel('BER');
ylim([0 0.024])
xlim([0 20])
title('Alamouti e OFDM - 2x2 - Sinal 1')







