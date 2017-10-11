%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% Special Project Code - Diphthong with Wall Loss
%
% Author: Chad McKell
% Date: 19 Apr 17
%
% Description: Finite-difference simulation of a diphthong with wall    
% losses. This script models glottal excitation, vocal-tract sound  
% propagation with wall vibration, and radiation from the mouth for two 
% vowels. The vowels are joined using a transition function.
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
tic; clear; close all;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% Variable Preamble
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

% Set variables
fs = 44100; % sampling frequency [samples/sec]
f0 = 100; % fundamental frequency [Hz]
c = 340; % speed of sound [meter/sec]
L = 0.1746; % length of vocal tract [meter] (/open_a/)
T = 0.3; % total duration of simulation [sec]
k = 1/fs; % sampling duration [sec]
NF = floor(T*fs); % total length of simulation [samples]
d = c*k; % distance sound travels between samples ('sample distance') [meter]
h = d/L; % fraction of tube length represented by one 'sample distance'
N = floor(L/d); % number of complete 'sample distances' in tube 
h = 1/N; % 'h' increased so that 'sample distances' fit evenly in tube
lambda = (d/L)*(1/h); % courant no. (ratio of old 'h' to new/bigger 'h')
gamma = c/L; % time required for sound to travel length L [sec]
p = 0.01; % scaling factor for surface areas 

% /open_a/
Sa = [0 0.45;1 0.20;2 0.26;3 0.21;4 0.32;5 0.30;6 0.33;...
    7 1.05;8 1.12;9 0.85;10 0.63;11 0.39;12 0.26;13 0.28;...
    14 0.23;15 0.32;16 0.29;17 0.28;18 0.40;19 0.66;20 1.20;...
    21 1.05;22 1.62;23 2.09;24 2.56;25 2.78;26 2.86;27 3.02;...
    28 3.75;29 4.60;30 5.09;31 6.02;32 6.55;33 6.29;34 6.27;35 5.94;...
    36 5.28;37 4.70;38 3.87;39 4.13;40 4.25;41 4.27;42 4.69;43 5.03];

% Normalize first column of Sa
Sa(:,1) = Sa(:,1)/max(Sa(:,1)); 

% Compute interpolated values of the function Sa(:,2)(Sa(:,1))
Sa = interp1(Sa(:,1), Sa(:,2), 0:h:1)';

% Scale surface areas of Sa
Sa = p*Sa;

% Set surface areas of tube slices 
% /i/
Si = [0 0.33;1 0.30;2 0.36;3 0.34;4 0.68;5 0.50;6 2.43;...
    7 3.15;8 2.66;9 2.49;10 3.39;11 3.80;12 3.78;13 4.35;...
    14 4.50;15 4.43;16 4.68;17 4.52;18 4.15;19 4.09;20 3.51;...
    21 2.95;22 2.03;23 1.66;24 1.38;25 1.05;26 0.60;27 0.35;...
    28 0.32;29 0.12;30 0.10;31 0.16;32 0.25;33 0.24;34 0.38;35 0.28;...
    36 0.36;37 0.65;38 1.58;39 2.05;40 2.01;41 1.58];

% Normalize first column of Si
Si(:,1) = Si(:,1)/max(Si(:,1)); 

% Compute interpolated values of the function Si(:,2)(Si(:,1))
Si = interp1(Si(:,1), Si(:,2), 0:h:1)';

% Scale surface areas of Si
Si = p*Si;

% Initialize time-varying cross-sectional areas
S = Si*0;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% Define damping oscillator parameters
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
omega0 = 500; % fundamental frequency of vocal tract walls
sigma0 = 405000; % damping coefficient
rho = 1.225; % density of air [kg/m^3]
M = 4.76; % mass per unit area of vocal tract walls [kg/m^2] (Titze,1988)

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% Define transition function 'beta'
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

% Set transition parameters
dur = 0.08; % duration of transition [sec]
durint = dur*fs; % duration in samples
startT = floor(NF/2 - durint/2); % start of transition
endT = floor(NF/2 + durint/2); % end of transition

% Initialize transition function
beta = ones(NF, 1);

% Set transition window
beta(startT:endT-1) = 0.5*(1 - cos(pi*(0:durint-1)/durint));

% Set beginning state
beta(1:startT) = 0;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% Define input impulse train
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
t = 0:k:T; % time bins
uin = sin(2*pi*f0*t); % sine wave with fundamental frequency f0
uin = 0.5*(uin+abs(uin)); % convert negative sinusoidal values to zeros

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% Compute output signal
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

% Initialize output sound signal
out = zeros(NF,1);

% Initialize velocity potential vectors. N+1 is last slice
psi = zeros(N+1,1); % velocity potential (n+1 time step)
psi1 = zeros(N+1,1); % velocity potential (n time step)
psi2 = zeros(N+1,1); % velocity potential (n-1 time step)

% Initialize coupling vectors (N-1 is the last slice)
w = zeros(N-1,1); % current time step
w1 = zeros(N-1,1); % 1 time step back
w2 = zeros(N-1,1); % 2 time steps back

for n = 1:NF
    
    % Calculate time-varying area function
    S = (1 - beta(n))*Sa + beta(n)*Si;
    
    % Compute coupling coefficient
    eps = c*sqrt(2*rho/M)*(pi/(S(1)))^(1/4);
   
    % Set average surface area at glottis, in the vocal tract, and at the lips
    Sav = [S(1); 0.25*(S(3:N+1)+2*S(2:N)+S(1:N-1)); S(N+1)];
    
    % Calculate coefficients of glottal excitation component
    g0 = 2*(1-lambda^2);
    gr = 2*lambda^2;
    gx = (k^2*gamma^2/h/S(1))*(3*S(1)-S(2));

    % Calculate coefficients of vocal-tract propagation component
    A = 1/k^2 + sigma0/k;
    B = eps*S(2:N).^(0.25)*k./(2*Sav(2:N));
    C = 2/k^2 - omega0^2;
    D = -2/k^2;
    E = eps*S(2:N).^(0.25)/(2*k);
    F = sigma0/k - 1/k^2;
    sl = 0.5*lambda^2*(S(2:N)+S(1:N-1))./Sav(2:N);
    s01 = g0; 
    s02 = B.*E/A - 1;
    w01 = B*C/A;
    w02 = B*D/A;
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

    % Calculate coefficients for coupling component
    w03 = C/A;
    w04 = F/A;
    wp = E/A;
    
    % Calculate velocity potential at the glottis
    psi(1) = g0*psi1(1) + gx*uin(n) + gr*psi1(2) - psi2(1);
    
    % Calculate velocity potential in the vocal tract
    psi(2:N) = s01*psi1(2:N) + s02.*psi2(2:N) + sl.*psi1(1:N-1) ...
                + sr.*psi1(3:N+1) - w01.*w1 - w02.*w2;
    
    % Calculate velocity potential at the lips
    psi(N+1) = r0*psi2(N+1) + rl*psi1(N) + r*psi1(N+1);
    
    % Calculate pressure at the lips
    out(n) = (fs/gamma)*(psi(N+1) - psi1(N+1)); 
    
    % Compute coupling term
    w = w03*w1 + w04*w2 + wp.*(psi(2:N) - psi2(2:N));
    
    % Set values equal to next grid line in time
    psi2 = psi1; 
    psi1 = psi;
    w2 = w1;
    w1 = w;
end

% Normalize output signal
out = out/max(out);

% Modify time bins to match length of output signal
t_modx = t(1:NF)';

% Define half-Hann window parameters
window_dur = 0.01; % window duration
start = floor(fs*(T-window_dur)); % starting sample index of window
durwin = floor(window_dur*fs); % duration of window in samples

% Apply half-Hann window
out(start:NF)=0.5*(1-cos(pi*(start:NF)/durwin)).*out(start:NF)';

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% Listen to output signal 
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
soundsc(out,fs);

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% Produce plots 
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

% Input signal
subplot(2,1,1)
plot(t, uin);
xlim([0.1 0.15]);
title('Waveform of Input Signal');
xlabel('Time (sec)');
ylabel('|uin|');

% Output signal
subplot(2,1,2)
plot(t_modx, out);
title('Waveform of Output Signal');
xlabel('Time (sec)');
ylabel('|out|');

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% Check code efficiency
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
toc % print elapsed time
