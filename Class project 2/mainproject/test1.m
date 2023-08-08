clear
clc
close all
I = audioread('project\Sample1.wav');%read
%window
Q = I';
N = 256;
Hamm = hamming(N);
frame = 70;
M = Q(((frame - 1) * (N / 2) + 1):((frame - 1) * (N / 2) + N)); 
Frame = M .* Hamm';
P = input('enter predictor para = '); %predictor
ai = lpc(Frame',P);                   % lpc
LP = filter([0 -ai(2:end)],1,Frame); % Frame and filter
E = Frame - LP;
figure(1)
subplot(2,1,1),plot(1:N,Frame,1:N,LP,'-r');grid;
title('origin and predict wave（hanning）')
subplot(2,1,2),plot(E);grid;
title('en predict');

%% graph
ai1 = lpc(I,P); % lpc
LP1 = filter([0 -ai(2:end)],1,I); % filter
figure(2)
subplot(2,1,1),plot(I);axis([0 length(I) -2 2]); grid;
title('origin wave')
subplot(2,1,2),plot(LP1);axis([0 length(LP1) -2 2]); grid;
title('predict wave')
E1 = I - LP1;

%% wave analysis
figure(3)
subplot(2,1,1); 
specgram(I,N,8000,N,N*3/4);
title('origin graph');
subplot(2,1,2);
specgram(LP1,N,8000,N,N*3/4); 
title('predict graph');
figure(4)
[cr,lags] = xcorr(E1,'coeff');
plot(lags,cr), grid
title 'coeff'
xlabel 'Lags', ylabel 'Normalized value'

%%sound
sound(LP1,16000);