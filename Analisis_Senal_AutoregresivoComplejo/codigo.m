%% CSTAR

% señal analitica z(n)
load handel
x = y;
clear y
N = length(x);
xl = [x;zeros(N,1)];
Xl = fft(xl);
Zl = [Xl(1);2*Xl(2:N);Xl(N+1);zeros(N-1,1)];
z = ifft(Zl);
z = [z(1:N);zeros(N,1)];
% z = hilbert(x);

% z(n) segementado por ventanas Lw
L = length(z);
Lw = 244;
m = ceil(L/Lw);
zn = [z(:);zeros((Lw*m)-L,1)];
zm = reshape(zn,Lw,m);
clear zn

% Error de prediccion directa
Lc = 244;
cm = zeros(Lc+1,m);
rc = zeros(Lc,m);
for i=1:m
    [cm(:,i),~,rc(:,i)] = aryule(zm(:,i),Lc);
end
mag = zeros(Lw,m);
for i=1:m
    [mag(:,i), F] = freqz(1,cm(:,i),Lw,Fs);
end
t = (0:L-1)*(1/Fs);
Pmag = abs(mag).^2;
NPmag = 10*log10(Pmag);

% Plot
imagesc(t,F,NPmag)
xlim([0,N*(1/Fs)])
title('CSTAR')
xlabel('Tiempo (s)')
ylabel('Frecuencia (Hz)')
cc = colorbar;
cc.Label.String = 'Nivel de Potencia instantánea (dB)';
colormap('jet')
set(gca,'Fontsize',17)
