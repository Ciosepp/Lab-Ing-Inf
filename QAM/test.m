close all
clear
clc
 
M=64;
L=sqrt(M);
 

SNRdB = -10:20;
SNR = 10.^(SNRdB/10);
semilogy(SNRdB,(1/log2(L))*((L-1)/L).*erfc(sqrt((3.*SNR)./(2*(L^2-1)))));ogy(10*log10(EbNo),(1/log2(L))*((L-1)/L).*erfc(sqrt(EbNo*3*log2(L)/(L^2-1))));