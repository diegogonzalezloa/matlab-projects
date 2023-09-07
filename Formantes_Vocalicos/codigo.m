%% Detección de Formantes en Segmentos de Tiempo
% Este codigo sirve para detectar formantes vocálicos en tiempo corto
% implementando un modelo autoregresivo complejo

close,clear,clc

%------------------------ Cargar speech y(n) ------------------------------
% Variable de entrada
in = input('Si desea cargar el audio digite 1 o si desea grabar el audio digite 2: ');
% Variable de entrada cargar audio
if in==1
    filename = uigetfile('*.wav');   % Cargar audio
    % Si la carga se cancela o no
    if isequal(filename,0)           % Carga cancelada
       disp('Selección cancelada');
       l = 0;
    else
       [x,Fs] = audioread(filename); % Leer archivo de audio
       l = 1;
    end
% Variable de entrada grabar audio
elseif in==2
    Fs = 10000;                      % Frecuencia de muetreo
    x = audiorecorder(Fs,24,1);      % Grabar audio
    disp('Start speaking.')          % -------
    recordblocking(x,2);             % -------
    disp('End of Recording.');       % -------
    x = getaudiodata(x);             % Obtener información del audio
    l = 1;
% Variable de entrada incorrecta
else
    disp('El número digitado es incorrecto, por favor vuelva a intentarlo');
    l = 0;
end

% Si el usuario escogio grabar o cargar sin cancelar
if (in==1 || in==2) && (l == 1)
    clear l in
    
%----------------------- Cuasiperiodo de la señal -------------------------
    % Autocorrelación cruzada para hallar el cuasiperiodo de la señal x(n)
    [corr,lags] = xcorr(x);
    [~,locsh] = findpeaks(corr);     % Picos cortos
    short = mean(diff(locsh));       % Diferencia media de los picos cortos
    [pklg,loclg] = findpeaks(corr,'MinPeakDistance',ceil(short)...
        ,'MinPeakHeight',0.707*max(corr)); % Picos largos
    zers = find(~lags);              % Punto central de ACF
%     TF = isempty(loclg);
    if length(loclg)>1
        for j=1:length(loclg)
            if loclg(j)==zers
                Lw = loclg(j+1)-loclg(j);
                break
            end
        end
    else
        Lw = 80;
    end
    clear i locsh short zers
    
%----------------------- Señal analitica z(n) -----------------------------
    z = hilbert(x);                  % Señal analitica de x(n)
    N = length(z);                   % Longitud de z(n)
    
%---------------------- Ventanas de longitud Lw ---------------------------
    Lw = round(Lw);                  % Lw redondeado
    m = ceil(N/Lw);                  % Segmentos
    zn = [z(:);zeros((Lw*m)-N,1)];   % Dimensión cuadrada de la matriz ... 
                                     % segmentada
    zm = reshape(zn,Lw,m);           % z(n) segmentado
    clear zn
    
%---------- Filtro pasaalto de pre-emphasis a z(n) segmentado -------------
    b = 0.95;                        % Grado del filtro
    zer = [1 -b];                    % Polinomio del cero
    zpre = zeros(Lw,m);              % Vector de ceros
    % Filtro de pre-emphasis tipo FIR de primer orden
    for i=1:m
        zpre(:,i) = filter(zer,1,zm(:,i));
    end
    clear i
    
%------------------------- Ventaneo con Hamming ---------------------------
    zw = hamming(Lw);                % Ventana Hamming y longitud Lw
    % Reducción de fuga espectral a cada ventana
    for i=1:m
        zpre(:,i) = zpre(:,i).*zw;
    end
    clear zw
    
%--------------------- Coeficientes del modelo AR -------------------------       
    Lc = 4;                          % Orden del modelo o suavisado 
    % Si la longitud ventana es mas pequeña que el orden del modelo 
    if Lc>Lw
        Lc = Lw;
    end
    % Calculo de los coeficientes por el metodo de Yule-Walker
    c = zeros(Lc+1,m);               % Vector de ceros
    for i=1:m
        c(:,i) = aryule(zpre(:,i),Lc);
    end
    clear i
    
%-------- Calculo de las raices del polinomio de coeficientes c -----------
    r = zeros(Lc,m);                 % Vector de ceros
    % Raices de los coeficientes
    for i=1:m
        r(:,i) = roots(c(:,i));
    end
    clear i
          
%--------- Argumento de los polos con parte imaginaria positiva -----------
    % tan^-1 (argumentos de los polos)    
    arg = zeros(Lc,m);               % Vector de ceros
    for i=1:m
        arg(:,i) = atan2(imag(r(:,i)),real(r(:,i)));
    end
    clear i
    arg = sort(arg);                 % Argumentos en orden ascendente
    % Argumentos positivos de los polos
    argp = zeros(Lc,m);              % Vector de ceros
    for i=1:m
        for j=1:Lc
            if arg(j,i)>=0
                argp(j,i)= arg(j,i);
            end
        end
    end
    clear i j
    % Invirtiendo los valores nulos
    for i=1:m
        if argp(1,i)==0
            zers = find(~argp(:,i)); % Posición de los ceros en el vector
            zers = zers(end);        % Posición del ultimo cero
            argp(:,i) = [argp(zers+1:end,i);zeros(zers,1)];
        end
    end
    clear zers i
    
%--------------------------- Formantes en Hz ------------------------------
    % Agumentos en Hz
    argHz = (argp*(Fs/2))/pi;
    argHz = argHz(1:3,:)';
    % Normalización de los argumentos en Hz de 0 a 1
    min_argHz = min(argHz);          % Frecuencia minima de cada formantes
    pos_argHz = argHz-min_argHz;     % Resta de la frecuencia minima
    max_argHz = max(pos_argHz);      % Frecuencia maxima de cada formantes
    norm_argHz = pos_argHz./max_argHz; % Frecuencias normalizadas
    % k-means clustering
    [~,C] = kmeans(norm_argHz,5);
    CHz = min_argHz+(C.*max_argHz);
%   Estimación de formantes
    load ('Forman.mat')
    forman = forman';
    D = zeros(m,5);
    for i=1:m
        for j=1:5
            diff1 = pdist2(forman(j,1),argHz(i,1),'euclidean');
            diff2 = pdist2(forman(j,2),argHz(i,2),'euclidean');
            D(i,j) = diff1+diff2;
        end
    end
    clear i j diff1 diff2
    minimos = zeros(1,5);
    k = zeros(1,5);
    for i=1:5
        [minimos(i),k(i)] = min(D(:,i));
    end
    clear i
    k2 = zeros(1,5);
    for i=1:5
        if minimos(i)<=50
            k2(i) = k(i);
        end
    end
    clear i
    zers = find(k2);
    formantes1 = zeros(5,3);
    for i=1:length(zers)
        l = zers(i);
        formantes1(l,:) = argHz(k2(l),:);
    end
    clear i l
    D2 = zeros(5,5);
    for i=1:length(zers)
        l = zers(i);
        for j=1:5
            diff1 = pdist2(formantes1(l,1),CHz(j,1),'euclidean');
            diff2 = pdist2(formantes1(l,2),CHz(j,2),'euclidean');
            D2(l,j) = diff1+diff2;
        end
    end
    clear i j diff1 diff2 l
    minimo2 = zeros(5,1);
    for i=1:5
        minimo2(i) = min(D2(i,:));
    end
    clear i
    formantes2 = zeros(5,3);
    for i=1:length(zers)
        l = zers(i);
        if minimo2(l)<150
            formantes2(l,:) = formantes1(l,:);
        end
    end
    clear i l
    zers2 = find(~formantes2(:,1));
    L_zers2 = length(zers2);
    if L_zers2==5
        [~,k3] = min(minimos);
        formantes2(k3,:) = argHz(k(k3),:);
        zers = k3;
    end
    
    % disp
    format bank
%     disp('Los formantes 1 son:')
%     disp(formantes1')
    disp('Los formantes 2 son:')
    disp(formantes2')
    
% -------------------------------- Plots ----------------------------------

    % Plot ACF
    figure('Name','ACF','NumberTitle','off');
    hold on
    plot(lags(loclg)/Fs,pklg+1.5,'vk')
    plot(lags/Fs,corr)
    hold off
    xlim([(-1)*((max(lags)/Fs)/8),((max(lags)/Fs)/8)])
    title('Picos de la Función de Autocorrelación (ACF) de la Señal')
    xlabel('Tiempo (s)')
    ylabel('Amplitud')
    lgd = legend('Picos mas Altos');
    title(lgd,'Picos')
    set(gca,'Fontsize',15)
    grid on 
    grid minor
    
    % Plot señal analitica
    X = fft(x);
    Z = fft(z);
    freq = (0:N-1)*(Fs/N);
    figure('Name','Señal Analitica','NumberTitle','off');
    subplot(211), plot(freq,abs(X).^2), axis('tight')
    title('Contenido Frecuencial de la Señal Original')
    xlabel('Frecuencia (Hz)')
    ylabel('Nivel Potencia Instantanea (dB)')
    set(gca,'Fontsize',15)
    grid on 
    grid minor
    subplot(212), plot(freq,abs(Z).^2), axis('tight')
    title('Contenido Frecuencial de la Señal Analitica')
    xlabel('Frecuencia (Hz)')
    ylabel('Nivel Potencia Instantanea (dB)')
    set(gca,'Fontsize',15)
    grid on 
    grid minor
    
    % Plot Envolvente media filtro de pre-emphasis
    figure('Name','Filtro de Pre-emphasis','NumberTitle','off');
    % Cofecicientes
    cmean = zeros(Lc+1,1);
    for i=1:Lc+1
        cmean(i) = mean(c(i,:));
    end
    [mag,f] = freqz(1,cmean,Lw,Fs);
    P = abs(mag).^2;
    dB = 10*log10(P);
    % Coeficientes sin filtro pre-emphasis
    zw = hamming(Lw);
    for i=1:m
        zm(:,i) = zm(:,i).*zw;
    end
    clear zw
    c2 = zeros(Lc+1,m);
    for i=1:m
        c2(:,i) = aryule(zm(:,i),Lc);
    end
    clear i
    cmean2 = zeros(Lc+1,1);
    for i=1:Lc+1
        cmean2(i) = mean(c2(i,:));
    end
    [mag2,f2] = freqz(1,cmean2,Lw,Fs);
    P2 = abs(mag2).^2;
    dB2 = 10*log10(P2);
    hold on
    plot(f,dB2,'--','LineWidth',0.5), axis('tight')
    plot(f,dB,'LineWidth',1)
    hold off
    title('Envolvente Promedio de la Señal')
    xlabel('Frecuencia (Hz)')
    ylabel('Nivel de Potencia (dB)')
    lgd = legend('Envolvente sin Filtro Pre-emphasis',...
        ['Envolvente con Filtro Pre-emphasis con grado ',num2str(b)]);
    title(lgd,'Envolventes')
    set(gca,'Fontsize',15)
    grid on
    grid minor
    
    % Plot ordenes
    figure('Name','Ordenes','NumberTitle','off');
    % Señal original por segmentos
    xn = [x(:);zeros((Lw*m)-N,1)];
    xm = reshape(xn,Lw,m);
    clear xn
    % Señal original con filtro de pre-emphasis
    xpre = zeros(Lw,m);
    for i=1:m
        xpre(:,i) = filter(zer,1,xm(:,i));
    end
    clear i
    % Señal original con ventana hamming
    xw = hamming(Lw);
    for i=1:m
        xpre(:,i) = xpre(:,i).*xw;
    end
    clear xw
    % Coeficientes con orden Lc2
    Lc2 = 2*round(((Fs/2)/1000));
%     Lc2 = 4;
    if Lc>Lw
        Lc2 = Lw;
    end
    cx = zeros(Lc2+1,m);
    for i=1:m
        cx(:,i) = aryule(xpre(:,i),Lc2);
    end
    clear i
    % Media de los coeficientes
    cxmean = zeros(Lc2+1,1);
    for i=1:Lc+1
        cxmean(i) = mean(cx(i,:));
    end
    [magx,fx] = freqz(1,cxmean,Lw,Fs);
    Px = abs(magx).^2;
    dBx = 10*log10(Px);
    hold on
    plot(f,dBx,'LineWidth',1), axis('tight')
    plot(f,dB,'LineWidth',1)
    hold off
    title('Envolvente Promedio de la Señal')
    xlabel('Frecuencia (Hz)')
    ylabel('Nivel de Potencia (dB)')
    lgd = legend(['Envolvente de la señal original con orden '...
        ,num2str(Lc2)],['Envolvente de la señal compleja con orden '...
        ,num2str(Lc)]);
    title(lgd,'Envolventes')
    set(gca,'Fontsize',15)
    grid on
    grid minor
    
else
    clear l in 
    disp('No hay variables para analizar')
end