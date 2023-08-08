audioPath = 'D:\Audioproject1\source\';
audioDir = dir([audioPath '*.mp3']);
[Y, Fs] = audioread([audioPath audioDir(3).name]);
y = Y(:,1);
ty = (0:length(y)-1)/Fs;
wind = hamming(1024);
overlap = 400;
nfft = 2^nextpow2(length(wind));
[s,F,T] = stft(y,Fs,'Window',wind,'OverlapLength',overlap,'FFTLength',nfft);
smag = abs(s);
sphs = angle(s);
[x,tx,info] = stftmag2sig(smag,nfft,Fs,'Window',wind,'OverlapLength',overlap);
plot(ty,y,tx+500/Fs,x+1);
legend('Original','Reconstructed','Location','best');
sound(x,Fs);
%audiowrite('recaudio12.wav',x,Fs);
disp(info);
figure;
imagesc(T, F, log10(abs(S)))  
colorbar;
set(gca, 'YDir', 'normal')
xlabel('Time (secs)')
ylabel('Freq (Hz)')
title('STFT spectrum')