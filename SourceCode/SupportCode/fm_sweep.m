% Function to generate fm_sweep Call
% f_start   [Hz]    Start Frequency
% f_end     [Hz]    End Frequency
% fs        [Hz]    Sample Frequency
% duration  [ms]    Call duration

function sig=fm_sweep(f_start,f_end,fs,duration, amplitude, winprct)

%funci�n para generar la se�al emitida por el murci�lago
%esta se�al es un barido FM que empieza en la frecuencia f_start
%hasta la frecuencia f_end. duration es su duraci�n en milisegundo (!)
%fs es la frecuencia de muestro. 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ms=10^-3;%la escala de tiempo es en milisegundos

duration=duration*ms;%convertir en segundos

t=0:1/fs:duration;%tiempos

ft=1./(1/f_start+t*(1/f_end-1/f_start)/duration);%vector con la modulaci�n de frecuencia (lineal con el periodo)

sig = amplitude * (sin(2*pi*(duration/(1/f_end-1/f_start)) * (log(1/f_start+t*(1/f_end-1/f_start)/duration) - log(1/f_start))));

length_win = round(length(sig)*2*winprct/100);

window_short = hanning(length_win);
window_on = window_short(1:ceil(length_win/2));
window_off = window_short(ceil(length_win/2) + 1 : end);

window = ones(1,length(sig));

window(1:length(window_on)) = window_on;
window(end-length(window_off) + 1 : end) = window_off;

sig = sig .* window;