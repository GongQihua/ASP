audioPath = 'D:\Audioproject1\source\';
open stft.m
audioDir = dir([audioPath '*.mp3']);
[Y, Fs] = audioread([audioPath audioDir(4).name]);
y = Y(:,1);
ty = (0:length(y)-1)/Fs;
wind = hamming(1024);
olen = 400;
nfft = 2^nextpow2(length(wind));
[s,F,T] = stft(y,Fs,'Window',wind,'OverlapLength',olen,'FFTLength',nfft);
smag = abs(s);
sphs = angle(s);
[x,tx,info] = stftmag2sig(smag,nfft,Fs,'Window',wind,'OverlapLength',olen);
plot(ty,y,tx+500/Fs,x+1);
legend('Original','Reconstructed','Location','best');
%sound(x,Fs);
disp(info);
S1 = spectrogram(y, wind, nooverlap, nfft, Fs);
S2 = spectrogram(x, wind, nooverlap, nfft, Fs);
num = sum(sum(abs(S1).^2))/(2*pi);
den = sum(sum((abs(S1)-abs(S2)).^2))/(2*pi);
ser = 10*log10(num/den);
disp(ser);