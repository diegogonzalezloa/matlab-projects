classdef ITULoudness
    % ITULoudness calcula los diferentes pasos estrabelcidos en la ITU 1770
    % para medir la sonoridad integrada de un archivo de audio
    
    properties (Access = public)
        signal  % Señal de entrada
        duration % Duración en seg de la señal
        numchannels % Número de canales
        samplerate % Frecuencia de muestreo
        integrateloudness % Sonoridad integrada
        truepeakmax % Pico real maximo de la señal
    end
    
    properties (Access = private)
    end

    methods (Access = public)
        function obj = ITULoudness(signalval, samplerateval)
            % Contructor del objeto ITULoundess
            
            % ####### Inicialización de variables ######
            
            obj.signal = signalval;
            [numsamples, obj.numchannels] = size(signalval);
            obj.duration = round(numsamples*(1/samplerateval),1);
            obj.samplerate = samplerateval;
            
            % #### Resampleo de la señal a 48kHz para que los coeficientes
            % funcionen ########
            
            Fs = 48e3;
            dt = 1/Fs;
            x = resample(obj.signal,Fs,obj.samplerate);
            [N,Ch] = size(x);
            
            % ########## Filtro ponderado K #######
            
            % Coeficientes filtro shelving
            a1 = -1.69065929318241;
            a2 = 0.73248077421585;	
            b0 = 1.53512485958697;
            b1 = -2.6916918940638;	
            b2 = 1.19839281085285;
            aSH = [1, a1, a2];
            bSH = [b0, b1, b2];
            
            % Coeficientes filtro pasa altos
            a1 = -1.99004745483398;
            a2 = 0.99007225036621;
            b0 = 1.0;
            b1 = -2.0;
            b2 = 1.0;
            aHP = [1, a1, a2];
            bHP = [b0, b1, b2];
            
            % k-weighting
            x = filter(bHP, aHP, x);
            y = filter(bSH, aSH, x);
            
            % ##### Cuadrado medio de cada bloque de sincronización ######
            
            N400 = floor(0.4/dt); % Número de muestras para 400ms (Muestra 
            % mas cercana)
            N_sinc = ceil((N-N400)/(N400*0.25)); % Número de bloques de 
            % sincronización

            mz = zeros(N_sinc,Ch); % Cuadrado medio para cada bloque de 
            % sincronización
            for i=1:N_sinc % Bloque de sincronización
                for j=1:Ch % Canales
                    if i < N_sinc
                        mz(i,j) = sum(y((N400*(i-1)*0.25)+1:(N400*...
                            (((i-1)*0.25)+1))+1,j).^2)/N400;
                    else
                        mz(i,j) = sum(y((N400*(i-1)*0.25)+1:end,j))/N400;
                    end
                end
            end
            
            %  ######## Suma ponderada de todos los canales ########
            
            ml = zeros(N_sinc,1); % Sonoridad de cada bloque de 
            % sincronización sin las compuertas
            G = [1,1,1,1.4,1.4]; % Pondercaión de cada canal

            for i=1:N_sinc % Bloque de sincronización
                suma = 0;
                for j=1:Ch % Canales
                    suma = suma+(G(j)*mz(i,j));
                end
                ml(i) = -0.691+(10*log10(suma));
            end
            
            % ####### Compuertas Gating (Absoluta y Relativa) #######
            
            % Umbral Absoluto
            Umb_Abs = any(ml>=-70,2); % Como la dimensión es 2 busca en las
            % filas valores mayor o igual a -70 y retorna valores logicos
            Az = mz(Umb_Abs,:); % Todas la filas de mz en donde los valores
            % de ml son mayor o igual a -70

            % Umbral Relativo
            [N_sinc2, ~] = size(Az);
            suma_ch = 0;
            for i=1:Ch
                suma_bloq = sum(Az(:,i))/N_sinc2;
                suma_ch = suma_ch+(G(i)*suma_bloq);
            end
            r = -0.691+(10*log10(suma_ch))-10; % Umbral relativo
            clear suma_bloq
            
            % ########### Sonoridad Integrada ############
            
            Al = ml(Umb_Abs,:); % Sonoridad de cada bloque de 
            % sincronización con el umbral absoluto
            Umb_Rel = any(Al>=r,2); % Umbral relativo
            Iz = Az(Umb_Rel,:);  % Bloques de sincronización despues del 
            % umbral relativo

            [N_sinc3, ~] = size(Iz);
            suma_ch = 0;
            % Suma media de cada canal
            for i=1:Ch
                suma_bloq = sum(Iz(:,i))/N_sinc3;
                suma_ch = suma_ch+(G(i)*suma_bloq);
            end
            % Sonoridad integrada
            obj.integrateloudness = round(-0.691+(10*log10(suma_ch)),1);
            
            % ########## True peak #############
            
            if (obj.samplerate >= 750) && (obj.samplerate < 1500)
                l = 256; % factor de sobremuestreo
            elseif (obj.samplerate >= 1500) && (obj.samplerate < 3000)
                l = 128; % factor de sobremuestreo
            elseif (obj.samplerate >= 3000) && (obj.samplerate < 6000)
                l = 64; % factor de sobremuestreo
            elseif (obj.samplerate >= 6000) && (obj.samplerate < 12000)
                l = 32; % factor de sobremuestreo
            elseif (obj.samplerate >= 12000) && (obj.samplerate < 24000)
                l = 16; % factor de sobremuestreo
            elseif (obj.samplerate >= 24000) && (obj.samplerate < 48000)
                l = 8; % factor de sobremuestreo
            elseif (obj.samplerate >= 48000) && (obj.samplerate < 96000)
                l = 4; % factor de sobremuestreo
            elseif (obj.samplerate >= 96000) && (obj.samplerate < 192000)
                l = 2; % factor de sobremuestreo
            end
            
            if Fs < 192000
                a = zeros(length(obj.signal)*l,Ch);
                a(1:l:end) = obj.signal;
                coeffInt = designMultirateFIR(l,1); % Coeficientes FIR del filtro de interpolación 
                % factor de inrepolación igual a 4 y un diezmado igual a 1
                b = filter(coeffInt,1,a);
                c = 20*log10(abs(b));
                obj.truepeakmax = max(max(c));
            else
                c = 20*log10(abs(obj.signal));
                obj.truepeakmax = max(max(c));
            end
        end
        
    end
    
end

