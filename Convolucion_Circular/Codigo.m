%% Codigo para el taller de convolucion circular por segmentos

% Este codigo tiene la función de hacer una convolucion circular por
% segmentos, entre una señal de audio de entrada convolucionada con una
% respuesta al impulso de un recinto.

%% Cargar los archivos

% En esta parte del codigo se cargan los archivos de audio, se calcula el
% canal mas energetico de la respuesta al impulso y se aplica un
% dowsampling a la respuesta al impulso para que que su Fs sea igual a la
% de la señal de entrada. Por ultimo Se normalizan las dos señales. Por
% ultimo se le aplica zero padding a h(n) para que tenga una dimencion de
% L+L-1.

% Audio de entrada x(n)
filenamex = 'Vox.wav';
[xn,Fsx] = audioread(filenamex);
xn = xn(:,1);
M = length(xn); % Longitud de xn
clear filenamex

% Audio de la respuesta la impulso h(n)
filenameh = 'middle_tunnel_4way_bformat.wav';
[hn1,Fsh] = audioread(filenameh);
clear filenameh
% se calcula el canal mas energetico de la respuesta al impulso
Ehn = sum(abs(hn1).^2);
[~,ma] = max(Ehn);
hn1 = hn1(:,ma);
clear Ehn ma
% La función resample remuestrea la señal de entrada hn1 en Fsh/Fsx veces
% la frecuencia de muestreo original. El resampleo aplica un filto pasa
% bajo antialiasing y compesa el retraso introducido por el filtro.
hn = resample(hn1,Fsh,Fsx);
L = length(hn); % Longitud de hn
hn = [hn;zeros(L-1,1)];
Fs = Fsx;
clear hn1 Fsh Fsx

% Normalización de las señales
% Pos_xn = xn - min(xn);
% xn = Pos_xn/(max(Pos_xn));
% Pos_hn = hn - min(hn);
% hn = Pos_hn/(max(Pos_hn));
% clear Pos_xn Pos_hn

%% Segmentos de Xn (Xk)

% En esta parte del codigo se parte la señal Xn en k segmentos (Xk), y se
% se crea un vector Xkf que es el segemnto final que contiene las ultimas
% muestras de la señal Xn y por lo general contiene el menor numero de
% muestras, entonces se le agregan ceros al vector para que tenga la misma
% dimension que los demas segmentos. Por ultimo se le hace zero_padding a
% todos los segmentos para que tengan una dimension de L+L-1.

s = ceil(M/L); % ceil function: redondea al valor entero mas alto
xn1 = [xn(:);zeros((L*s)-M,1)];
xk = reshape(xn1,L,s);
clear xn1
% Zero-padding para calcular la convolución circular
[f,c] = size(xk);
xk = [xk; zeros(f-1,c)];
clear f c

% floor function, redondea por defecto 
% Aplicar zero-padding √

%% Convolucion circular por segmentos

% En esta parte del codigo se hace la convolucion circular de cada segemnto
% por la respuesta al impulso, luego se une todas las respuestas de cada
% segmento (Ynk) para obtener la respuesta completa del sistema (Yn).

% yk = 0.039.*ifft(fft(xk).*fft(hn));
yk = ifft(fft(xk).*fft(hn));

% dimension de hn original
hn = hn(1:L);
clear xk

% Se solopan las respuestas yk para obtener la respuesta final del sistema
yni = zeros((L+L-1),1);
ynf = zeros(L,1);
yns = zeros(L,s-2);
for i=1:s
    if i==1
        for j=1:(L+L-1)
            if j<=L
                yni(j,1) = yk(j,i);
            else
                yni(j,1) = (yk(j,i)+yk(j-L,i+1))';
            end
        end
    elseif i==s 
        for j=L:(L+L-1)
                ynf(j+1-L,1) = yk(j,i);
        end
    else
        for j=L:(L+L-1)
            if j<=L
                yns(j+1-L,i-1) = yk(j,i);
            else
                yns(j+1-L,i-1) = yk(j,i)+yk(j-L,i+1);
            end
        end
    end
end
yn = [yni(:);yns(:);ynf(:)];
N = length(yn);
yn = yn(1:N-((L*s)-M));
clear i j s yk yni yns ynf N L M

% error: sum((abs(ync-yn))>1e-4)

% Preguntas:
% 1. Normalización de las señales.