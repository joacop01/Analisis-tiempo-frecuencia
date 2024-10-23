%Ej9
Fs = 1000; %Frecuencua de muestreo
Ts = 1/Fs; %Período de muestreo
t = 0:Ts:1-Ts;
N = length(t);
f = -Fs/2: Fs/N: (Fs/2) -1;

x1 = exp(-2000*(t-0.5).^2);

subplot(321);
plot(t,x1);
xlabel('Tiempo(s)');
ylabel('Amplitud');
title('x1(t)');

f = -Fs/2: Fs/N: Fs/2 -1;
fft_x1 = fft(x1);
fft_x1 = fftshift(fft_x1);
fft_x1 = abs(fft_x1);

subplot(322);
plot(f, fft_x1);
title('|x1(f)|');
xlabel('Frecuencia (Hz)');
ylabel('Amplitud');

x2 = exp(-2000*(t-0.25).^2);

subplot(323);
plot(t, x2);
title('x2(t)');
xlabel('Tiempo(s)');
ylabel('Amplitud');

fft_x2 = fft(x2);
fft_x2 = fftshift(fft_x2);
fft_x2 = abs(fft_x2);

subplot(324);
plot(f, fft_x2);
title('|x2(f)|');
xlabel('Frecuencia (Hz)');
ylabel('Amplitud');


c = dot(fft_x1, fft_x2)/(norm(fft_x1)* norm(fft_x2));

disp(['El grado de parecido entre las señales es de: ' num2str(c)]);


x3 = exp(-2000*(t-0.75).^2);

subplot(325);
plot(t, x3);
title('x3(t)');
xlabel('Tiempo(s)');
ylabel('Amplitud');


fft_x3 = fft(x3);
fft_x3 = fftshift(fft_x3);
fft_x3 = abs(fft_x3);

subplot(326);
plot(f, fft_x3);
title('|x3(f)|');
xlabel('Frecuencia (Hz)');
ylabel('Amplitud');
figure;

c = dot(fft_x1, fft_x3)/(norm(fft_x1)* norm(fft_x3));
disp(['El grado de parecido entre las señales es de: ' num2str(c)]);
%Ej 10

x4 = x1 .* cos(2*pi*250*t);
subplot(221);
plot(t, x4);
title('x4(t)');
xlabel('Tiempo(s)');
ylabel('Amplitud');

fft_x4 = fft(x4);
fft_x4 = fftshift(fft_x4);
fft_x4 = abs(fft_x4);

subplot(222);
plot(f, fft_x4);
title('|x4(f)|');
xlabel('Frecuencia (Hz)');
ylabel('Amplitud');

x5 = x1 .* exp(1i*2*pi*350*t);
subplot(223);
hold on
plot(t, real(x5), 'r');
title('x5(t)');
plot(t, imag(x5), 'b');
legend('Re(x5(t))', 'Im(x5(t))');
xlabel('Tiempo(s)');
ylabel('Amplitud');
hold off;
fft_x5 = fft(x5);
fft_x5 = fftshift(fft_x5);
fft_x5 = abs(fft_x5);

subplot(224);
plot(f, fft_x5);
title('|x5(f)|');
xlabel('Frecuencia (Hz)');
ylabel('Amplitud');
