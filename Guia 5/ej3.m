addpath('C:\Repositorios\Curso Análisis Tiempo-Frecuencia\Funciones');

Fs = 1000;
Ts = 1/Fs;
t = 0:Ts:1-Ts;
k = 200;
f0 = 100;
N = length(t);

x = exp(1i*2*pi*f0*t);
F = STFT_Gauss(x,t,k);
F_g = STFT_Gauss_diff(x, t, 200);

T = zeros(size(F));

f = 0: Fs/N : Fs/2 - Fs/N;

f_somb =  (F_g/F);

