 % Universidad de Guadalajara
% Centro Universitario de Ciencias Exactas e Ingenierias (CUCEI)
% Título del trabajo: "PROYECTO FILTRO DE AUDIO"
% Materia: Métodos Matemáticos III (MAPI)
% NRC: 57263, Sección D08
% Mtra. Laura Esther Cortes Navarro
% Alumno: 224697029 - Emiliano Sanchez Godinez

function Proyecto_audio() %FUNCION PRINCIPAL
    clear; close all; clc;
    
    %% -------- SELECCIÓN DE ARCHIVO DE AUDIO --------
    [filename, pathname] = uigetfile({'*.wav','Archivo WAV (*.wav)'});
    if isequal(filename,0)
        disp('Sin archivo seleccionado.');
        return;
    end
    
    [x,fs] = audioread(fullfile(pathname,filename));
    [N, canales] = size(x);
    t = (0:N-1)/fs;
    
    %% -------- Función en consola --------
    disp('+============================================+');
    disp('| CONJUNTO DE VALORES DEL AUDIO x(n)           |');
    disp('+============================================+');
    disp(x);
    
    export_filename = 'datos.xlsx';
    writematrix(x, export_filename);

    
    %% -------- SELECCIÓN DE VENTANA (var: tFiltro) --------
    disp('+-------------------------------------------------------------------------+');
    disp('|               Selecciona el tipo de Ventana:                            |');
    disp('+-------------------------------------------------------------------------+');
    disp('|  1. Bartlett   |  2. Barthannwin |  3. Blackman    |  4. BlackmanHarris |');
    disp('|  5. Bohman     |  6. Chebwin     |  7. Flattopwin  |  8. Gausswin       |');
    disp('|  9. Hamming    | 10. Hann        | 11. Kaiser      | 12. Nuttallwin     |');
    disp('| 13. Parzen     | 14. Rectwin     | 15. Taylor      | 16. Tukeywin       |');
    disp('| 17. Triangular |                 |                 |                    |');
    disp('+-------------------------------------------------------------------------+');   
    tFiltro = input('Filtro: ');
    
    % Parámetros extra
    cheb_atten = 80;
    kaiser_beta = 5;
    gauss_alpha = 2.5;
    tukey_r = 0.5;
    taylor_nbar = 4;

    % -------- FUNCIÓN SELECTORA DE VENTANA --------
    TFiltro = selector_ventana(tFiltro, cheb_atten, kaiser_beta, gauss_alpha, tukey_r, taylor_nbar);


    %% -------- SELECCIÓN DE TIPO DE FILTRO --------
    disp('+-------------------------------------------------------------+');
    disp('| Selecciona el TIPO de filtro a aplicar:                     |');
    disp('+-------------------------------------------------------------+');
    disp('| 1. Pasa bajas                                               |');
    disp('| 2. Pasa altas                                               |');
    disp('| 3. Pasa banda                                               |');
    disp('| 4. Rechaza banda (Notch)                                    |');
    disp('+-------------------------------------------------------------+');
    tipoFiltro = input('Tipo de filtro: ');


    %% -------- GRÁFICA 1: SEÑAL ORIGINAL --------
    figure;
    subplot(3,1,1);
    if canales == 2
        plot(t, x(:,1), t, x(:,2));
        legend('Canal 1','Canal 2');
    else
        plot(t, x);
    end

    grid on; axis tight;
    xlabel('Tiempo (s)'); ylabel('Amplitud');
    title('Señal original');
    
    % -------- FFT SIN SHIFT --------
    Y = fft(x(:,1));
    subplot(3,1,2);
    plot(abs(Y));
    grid on; axis tight;
    xlabel('Índice FFT'); ylabel('Amplitud');
    title('FFT sin shift');
    
    % -------- FFT CON SHIFT + FILTRO --------
    Yshift = fftshift(Y);
    f = linspace(-fs/2, fs/2, N);
    
    subplot(3,1,3);
    
    plot(f, abs(Yshift), 'b'); hold on;

%% Construcción del filtro según selección

% --- Límites de frecuencia ---
f1 = 3000;
f2 = 5000;

%Constructor

Filtro = zeros(N,1);

switch tipoFiltro
    case 1  % PASA BAJAS
        Filtro(abs(f) <= f1) = 1;

    case 2  % PASA ALTAS
        Filtro(abs(f) >= f1) = 1;

    case 3  % PASA BANDA
        Filtro(abs(f) >= f1 & abs(f) <= f2) = 1;

    case 4  % RECHAZA BANDA (NOTCH)
        Filtro(:) = 1;
        Filtro(abs(f) >= f1 & abs(f) <= f2) = 0;

    otherwise
        disp('Selección inválida → DEFAULT: Filtro pasa banda.');
        Filtro(abs(f) >= f1 & abs(f) <= f2) = 1;
end

Lw = 201;
ventana = TFiltro(Lw)';
ventana = ventana / max(ventana);
Filtro = conv(Filtro, ventana, 'same');

plot(f, Filtro * max(abs(Yshift)), 'r','LineWidth',1.4);

grid on; axis tight;
xlabel('Frecuencia (Hz)'); ylabel('Magnitud');
title('FFT con shift + Filtro');
legend('FFT','Filtro');


%% -------- APLICAR FILTRO --------
Yf_shift = Yshift .* Filtro;   % MISMA LONGITUD N
Yf = ifftshift(Yf_shift);      % DESHACER SHIFT
y_filtrado = real(ifft(Yf));   % IFFT → longitud N garantizada

% Normalizar amplitud
if max(abs(y_filtrado)) ~= 0
    y_filtrado = y_filtrado / max(abs(y_filtrado));
end
    
% -------- GRAFICA DE LA SEÑAL FILTRADA --------
figure;
subplot(2,1,1);
plot(t, x(:,1)); grid on; axis tight;
xlabel('Tiempo (s)'); ylabel('Amplitud');
title('Señal original (Canal 1)');

subplot(2,1,2);
plot(t, y_filtrado); grid on; axis tight;
xlabel('Tiempo (s)'); ylabel('Amplitud');
title('Señal filtrada (Canal 1)');


%% SUB-GRÁFICOS (Complemento personal)
% ---------- FIGURE 1: Filtro y Espectros ----------
figure;

subplot(3,1,1);
plot(f, Filtro, 'LineWidth', 1.2);
grid on; axis tight;
xlabel('Frecuencia (Hz)');
ylabel('Magnitud');
title('Respuesta en frecuencia del filtro aplicado');

subplot(3,1,2);
Yf = fftshift(fft(y_filtrado));
plot(f, abs(Yf));
xlabel('Frecuencia (Hz)');
ylabel('Magnitud');
title('Espectro de la señal filtrada');
grid on; axis tight;

subplot(3,1,3);
plot(f, abs(Yshift), 'b'); hold on;
plot(f, abs(Yf), 'r');
xlabel('Frecuencia (Hz)');
ylabel('Magnitud');
legend('Antes del filtro','Después del filtro');
title('Comparación espectral antes vs después de filtrar');
grid on; axis tight;

% ---------- FIGURE 2: Señales en el tiempo ----------
figure;
plot(t, x(:,1)); hold on;
plot(t, y_filtrado, 'r');
xlabel('Tiempo (s)');
ylabel('Amplitud');
legend('Original', 'Filtrada');
title('Comparación temporal: señal original vs filtrada');
grid on;


%% -------- REPRODUCCIÓN --------
p = audioplayer(x,fs); 
playblocking(p);

p = audioplayer(y_filtrado,fs);
playblocking(p);
end % FUNCION PRINCIPAL


%%  FUNCIÓN SELECTORA DE VENTANA 
function TF = selector_ventana(tFiltro, cheb_atten, kaiser_beta, gauss_alpha, tukey_r, taylor_nbar)

switch tFiltro
    case 1;  TF = @(L) bartlett(L);                     % Bartlett window.
    case 2;  TF = @(L) barthannwin(L);                  % Modified Bartlett-Hanning window.
    case 3;  TF = @(L) blackman(L);                     % Blackman window.
    case 4;  TF = @(L) blackmanharris(L);               % Minimum 4-term Blackman-Harris window.
    case 5;  TF = @(L) bohmanwin(L);                    % Bohman window.
    case 6;  TF = @(L) chebwin(L, cheb_atten);          % Chebyshev window.
    case 7;  TF = @(L) flattopwin(L);                   % Flat Top window.
    case 8;  TF = @(L) gausswin(L, gauss_alpha);        % Gaussian window.
    case 9;  TF = @(L) hamming(L);                      % Hamming window.
    case 10; TF = @(L) hann(L);                         % Hann window.
    case 11; TF = @(L) kaiser(L, kaiser_beta);          % Kaiser window.
    case 12; TF = @(L) nuttallwin(L);                   % Nuttall defined minimum 4-term Blackman-Harris window.
    case 13; TF = @(L) parzenwin(L);                    % Parzen (de la Valle-Poussin) window.
    case 14; TF = @(L) rectwin(L);                      % Rectangular window.
    case 15; TF = @(L) taylorwin(L, taylor_nbar);       % Taylor window.
    case 16; TF = @(L) tukeywin(L, tukey_r);            % Tukey window.
    case 17; TF = @(L) triang(L);                       % Triangular window.
    otherwise
        disp('Selección inválida → DEFAULT: Bartlett.');
        TF = @(L) bartlett(L);                          % Bartlett window.
   end
end
