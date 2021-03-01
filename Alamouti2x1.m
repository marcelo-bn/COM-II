clear all;
close all;
clc;

%% Dados
SNR = 0:5:40;
BER_s0 = zeros(size(SNR));
BER_s1 = zeros(size(SNR));

NrSub = 64;
Nofdm = 64;

NrIteracoes = 10^4;

%% Transmissor

% QAM
M = 16;  % Bits por simbolo
mod = comm.PSKModulator(M,'BitInput',true,'PhaseOffset',pi/M);
demod = comm.PSKDemodulator(M,'BitOutput',true,'PhaseOffset',pi/M);
%mod = comm.RectangularQAMModulator('ModulationOrder',M,'BitInput',true);
%demod = comm.RectangularQAMDemodulator('ModulationOrder',M,'BitOutput',true);

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

for nSNR = 1:length(SNR)
    Error_s0 = 0;
    Error_s1 = 0;
    for nIteracao = 1:NrIteracoes
        % Canal Rayleigh - Tabela I
        h0 = sqrt(0.5)*(randn+1j*randn);                % RX0 - TX0
        h1 = sqrt(0.5)*(randn+1j*randn);                % RX0 - TX1

        % Sinais no receptor - Eq. 11
        R0 = h0*ofdm_s0 + h1*ofdm_s1;                  
        R1 = -h0*conj(ofdm_s1) + h1*conj(ofdm_s0);   
  
        % Adicionando o ruido do canal
        SNR_value = SNR(nSNR);
        RO = awgn(R0,SNR_value);
        R1 = awgn(R1,SNR_value);
   
        % Sinais no combinador - Eq. 12
        s0_combinador = conj(h0)*R0 + h1*conj(R1);
        s1_combinador = conj(h1)*R0 - h0*conj(R1);    
        
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
BER_s0 = BER_s0./(Bits_info*NrSub);
BER_s1 = BER_s1./(Bits_info*NrSub);

%figure;
semilogy(SNR,BER_s0,'LineWidth',2); hold on; grid on;
semilogy(SNR,BER_s1,'LineWidth',2);
legend('s0','s1');
xlabel('SNR [dB]');
ylabel('BER');
title('Transmissao Alamouti 2x1')


