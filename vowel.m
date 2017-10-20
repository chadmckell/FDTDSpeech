%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% Special Project Code - Vowel 
% 
% Author: Chad McKell
% Date: 19 Apr 17
%
% Description: Finite-difference simulation of a vowel. This script models 
% glottal excitation, vocal-tract sound propagation, and radiation from the
% mouth.
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
tic; clear; close all;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% Variable Preamble
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

% Set variables
fs = 44100; % sampling frequency [samples/sec]
f0 = 60; % fundamental frequency [Hz]
c = 340; % speed of sound [meter/sec]
L = 0.1746; % length of vocal tract [meter] (/open_a/)
T = 1; % total duration of simulation [sec]
k = 1/fs; % sampling duration [sec]
NF = floor(T*fs); % total length of simulation [samples]
d = c*k; % distance sound travels between samples ('sample distance') [meter]
h = d/L; % fraction of tube length represented by one 'sample distance'
N = floor(L/d); % number of complete 'sample distances' in tube 
h = 1/N; % 'h' increased so that 'sample distances' fit evenly in tube
lambda = (d/L)*(1/h); % courant no. (ratio of old 'h' to new/bigger 'h')
gamma = c/L; % time required for sound to travel length L [sec]
p = 0.01; % scaling factor for surface areas 

% Set surface areas of tube slices (/open_a/)
S = [0 0.45;1 0.20;2 0.26;3 0.21;4 0.32;5 0.30;6 0.33;...
    7 1.05;8 1.12;9 0.85;10 0.63;11 0.39;12 0.26;13 0.28;...
    14 0.23;15 0.32;16 0.29;17 0.28;18 0.40;19 0.66;20 1.20;...
    21 1.05;22 1.62;23 2.09;24 2.56;25 2.78;26 2.86;27 3.02;...
    28 3.75;29 4.60;30 5.09;31 6.02;32 6.55;33 6.29;34 6.27;35 5.94;...
    36 5.28;37 4.70;38 3.87;39 4.13;40 4.25;41 4.27;42 4.69;43 5.03];

% Normalize first column of S
S(:,1) = S(:,1)/max(S(:,1));

% Compute interpolated values of the function S(:,2)(S(:,1))
S = interp1(S(:,1), S(:,2), 0:h:1)';

% Scale surface areas
S = p*S;

% Set average surface area at glottis, in the vocal tract, and at the lips
Sav = [S(1); 0.25*(S(3:N+1)+2*S(2:N)+S(1:N-1)); S(N+1)]; 

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% Define input impulse train
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
t = 0:k:T; % time bins
uin = sin(2*pi*f0*t); % sine wave with fundamental frequency f0
uin = 0.5*(uin+abs(uin)); % convert negative sinusoidal values to zeros

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% Compute coefficients of finite-difference algorithm
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

% Calculate coefficients of glottal excitation component
g0 = 2*(1-lambda^2);
gr = 2*lambda^2;
gx = (k^2*gamma^2/h/S(1))*(3*S(1)-S(2));

% Calculate coefficients of vocal-tract propagation component
sl = 0.5*lambda^2*(S(2:N)+S(1:N-1))./Sav(2:N);
s0 = g0; 
sr = 0.5*lambda^2*(S(3:N+1)+S(2:N))./Sav(2:N);

% Calculate coefficients for lip radiation component
alf1 = 1/(2*0.8216^2*gamma);
alf2 = L/(0.8216*sqrt(S(1)*S(N+1)/pi)); 
a = 0.5*lambda^2*h*(3*S(N+1)-S(N))/S(N+1);
q1 = alf1*a/k; 
q2 = alf2*a; 
rl = gr/(1+q1+q2); 
r0 = (q1-q2-1)/(1+q1+q2);
r = g0/(1+q1+q2);

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% Compute output signal using finite-difference algorithm
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

% Initialize output sound signal
out = zeros(NF,1); 

% Initialize velocity potential vectors (N+1 is the last spatial slice)
psi = zeros(N+1,1); % current time step
psi1 = zeros(N+1,1); % 1 time step back
psi2 = zeros(N+1,1); % 2 time steps back

for n = 1:NF
    
    % Calculate velocity potential at the glottis
    psi(1) = g0*psi1(1) + gx*uin(n) + gr*psi1(2) - psi2(1);
    
    % Calculate velocity potential in the vocal tract
    psi(2:N) = s0*psi1(2:N) + sl.*psi1(1:N-1) + sr.*psi1(3:N+1)- psi2(2:N);
    
    % Calculate velocity potential at the lips
    psi(N+1) = r0*psi2(N+1) + rl*psi1(N) + r*psi1(N+1);
    
    % Calculate pressure at the lips
    out(n) = (fs/gamma)*(psi(N+1) - psi1(N+1)); 
    
    % Set velocity potential vectors equal to next grid line in time
    psi2 = psi1; 
    psi1 = psi;
end

% Normalize output signal
out = out/max(out);

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% Listen to output signal 
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
soundsc(out,fs); 

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% Apply Hanning window to output signal  
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
n = 0:NF-1; % sample bins
out = 0.5*(1 - cos(2*pi*n'/NF)).*out; % apply window

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% Compute log-magnitude spectra  
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
OUT = abs(fft(out)); % unscaled magnitude spectra
OUTdB = 20*log10(OUT);% log-magnitude spectra

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% Define plotting parameters 
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
t_modx = t(1:NF)'; % time bins modified to match length of each signal
fx = n'*fs/NF; % frequency bins

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% Produce plots 
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

% Input signal
subplot(3,1,1)
plot(t, uin);
xlim([0.5 0.65]);
title('Waveform of Input Signal');
xlabel('Time (sec)'); 
ylabel('|uin|');

% Output signal
subplot(3,1,2)
plot(t_modx, out);
xlim([0.5 0.6]);
title('Waveform of Output Signal');
xlabel('Time (sec)'); 
ylabel('|out|');

% Log-magnitude spectra
subplot(3,1,3)
plot(fx, OUTdB);
title('Magnitude Spectra of Output Signal');
xlabel('Frequency (Hz)'); 
ylabel('|OUT| (dB)');
xlim([0 3000])

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% Check code efficiency
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
toc % print elapsed time
