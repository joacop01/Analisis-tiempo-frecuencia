function [] = Plot_STFT(F, t, f)
    figure;
    imagesc(t, f, abs(F));
    axis xy;
    xlabel('Tiempo (s)');
    ylabel('Frecuencia (Hz)');
    title('STFT');
    colorbar;
end