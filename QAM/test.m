close all
clear
clc
 

M=16;
L=sqrt(M);
 

SNR = 2*log2(L)*(10.^((0:1:15)/10));
SNRdB = 10*log10(SNR);
semilogy(SNRdB,(1/log2(L))*((L-1)/L).*erfc(sqrt((3.*SNR)./(2*(L^2-1)))));
 

EbNo = (1/(2*log2(L)))*SNR;
 

figure
semilogy(10*log10(EbNo),(1/log2(L))*((L-1)/L).*erfc(sqrt(EbNo*3*log2(L)/(L^2-1))));