audioPath = 'D:\Audioproject1\source\';
audioDir = dir([audioPath '*.mp3']);
%for i 1:length(audioDir):
[Y, Fs] = audioread([audioPath audioDir(1).name]);
y = Y(:,1);
figure;
frame = 64;
winsize = 128;
winchoice = hamming(winsize);
nfft = 2^nextpow2(length(winchoice));
overlap = winsize - 1;
[S, F, T, P] = spectrogram(y, winchoice, overlap, nfft, Fs);
STFTM = sum(y * winsize * (nfft-frame*overlap) * exp(-1j*winsize*nfft));
%magnitude = abs(S);
imagesc(T, F, log10(abs(S)))  
colorbar;
set(gca, 'YDir', 'normal')
xlabel('Time (secs)')
ylabel('Freq (Hz)')
title('STFT spectrum')
%saveas(gca,'D:\Audioproject1\compare\audio1hann.jpg');