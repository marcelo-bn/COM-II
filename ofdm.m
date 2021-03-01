clear all;
close all;
clc;

%% Dados
SNR = 0:5:40;
NrIteracoes = 10^4;

BER_s0 = zeros(size(SNR));
BER_s1 = zeros(size(SNR));

NrSub = 64;
Nofdm = 64;

%% Transmissor

% QAM
M = 16;  % Bits por simbolo
mod = comm.PSKModulator(M,'BitInput',true,'PhaseOffset',pi/M);
demod = comm.PSKDemodulator(M,'BitOutput',true,'PhaseOffset',pi/M);
%mod = comm.RectangularQAMModulator('ModulationOrder',16,'BitInput',true);
%demod = comm.RectangularQAMDemodulator('ModulationOrder',16,'BitOutput',true);

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

for nSNR=1:length(SNR)
    Error_s0 = 0;
    Error_s1 = 0;
    for nIter=1:NrIteracoes
        %% Canal
        Y_s0 = awgn(ofdm_s0,nSNR);
        Y_s1 = awgn(ofdm_s1,nSNR);

        %% Receptor
        % Retirar o prefixo ciclico
        RXofdm_s0 = Y_s0(cp+1:Nofdm+cp, 1);
        RXofdm_s1 = Y_s1(cp+1:Nofdm+cp, 1);

        % FFT
        rx_info_s0 = fft(RXofdm_s0);
        rx_info_s1 = fft(RXofdm_s1);

        % Demodulador M-PSK
        rx_s0 = step(demod,rx_info_s0);
        rx_s1 = step(demod,rx_info_s1);

        % Verificando erros
        erro_s0 = length(find(rx_s0 ~= s0));
        erro_s1 = length(find(rx_s1 ~= s1));
        Error_s0 = Error_s0 + erro_s0;
        Error_s1 = Error_s1+ erro_s1;
      
    end
    BER_s0(nSNR,1) = Error_s0/NrIteracoes;
    BER_s1(nSNR,1) = Error_s1/NrIteracoes;
end

BER_s0 = BER_s0./(Bits_info*NrSub);
BER_s1 = BER_s1./(Bits_info*NrSub);
figure;
semilogy(SNR,BER_s0,'LineWidth',2); hold on; grid on;
semilogy(SNR,BER_s1,'LineWidth',2);
xlabel('SNR [dB]');
ylabel('BER');
%legend('s0 - Alamouti 2x2','s1 - Alamouti 2x2', 's0 - Alamouti 2x1', 's1 - Alamouti 2x1', 's0 - OFDM', 's1 - OFDM');
legend('s0','s1');
title('Somente OFDM vs. Alamouti 2x1 vs. Alamouti 2x2')





