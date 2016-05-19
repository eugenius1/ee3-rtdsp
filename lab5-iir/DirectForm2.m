%% Eusebius & Prahnav Lab 5 RTDSP - IIR Filter Design - Order 4
clear;

order = 4;
wp = [180 450];     % Filter Passband
Wp = (1/4000)*wp;   % Normalized Frequency

Rs = 23;            % Stopband Ripple 
Rp = 0.4;           % Passband Ripple
Fs = 8000;          % Sampling Frequency

[b,a] = ellip(order/2, Rp ,Rs, Wp);    % Elliptical IIR Filter

% Plotting Filter response and z-plane
figure;
freqz(b,a);
figure;
zplane(b,a);

% Formatting and Writing to Text File
format long e;

formatSpec = '%2.15e';
str1 = num2str(b',formatSpec);
str1 = cellstr(str1);
data1 = strjoin(str1', ', ');
str2 = num2str(a',formatSpec);
str2 = cellstr(str2);
data2 = strjoin(str2', ', ');

fileID = fopen('fir_coef.txt', 'w');
fprintf(fileID, '#define N ');
fprintf(fileID, num2str(length(b)));
fprintf(fileID, '\n double b[]={');
fprintf(fileID, data1);
fprintf(fileID, '};\n');
fprintf(fileID, 'double a[]={');
fprintf(fileID, data2);
fprintf(fileID, '};\n');
fclose(fileID);
