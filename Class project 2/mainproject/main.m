clear; 
clc;
[x,fs] = audioread('project/Sample1.wav');
fl = 60; %frame length in ms
fo = 30; %frame offset in ms
order = 10; %filter order (1-20)
volume = 0; %noise volume in transmission channel in dB

if max(x) > abs(min(x)) % normalize input
  x = x/max(x);
else
  x = x/abs(min(x));
end

[b_bpf,a_bpf] = butter(2,[40/(fs/2),1000/(fs/2)]); %butter filter
xin = filter(b_bpf,a_bpf,x);

if max(x) > abs(min(x))
  x = x/max(x);
else
  x = x/abs(min(x));
end
L=round(fl*fs/1000);
R=round(fo*fs/1000);
if fl<1 
    fl=1; 
end
if fl>length(xin) % set inputs limit to xin
    fl=length(xin); 
end

if fo<1 
    fo=1; 
end
if fo>length(xin) 
    fo=length(xin); 
end

if order<1 
    order=1; 
end
if order>20 
    order=20; 
end

out=zeros((length(xin)+L),1); % define initial array output, excite, gain and pitch
g=zeros((length(xin)+L),1);
e=zeros((length(xin)+L),1);
pt=zeros((length(xin)+L),1);

n = 1; 
nframes = 0;
while (n+L-1 <= length(xin))

x=xin(n:n+L-1); %input frame                        

win=hamming(L); 
pitchmin = int32(fs/350);
pitchmax = int32(fs/80);
xw=x.*win;

[ak,gk] = lpc(xw,order);%#ok<*ASGLU> perform LPC analysis
p = order; q = 0;
x   = xw(:);
N   = length(x);
X   = convmtx(x,p+1);
Xq  = X(q+1:N+p-1,1:p);
a   = [1;pinv(-Xq)*X(q+2:N+p,1)];
G = x(q+2:N)'*X(q+2:N,1:p+1)*a;
b0 = sqrt(G);
[pitch,v] = pitchdetector(xw,pitchmin,pitchmax); % determine pitch
  
lsfs = poly2lsf(a);% encode LPC coeffs into line spectral pairs and concatenate
prms = [double(b0) lsfs' double(pitch) double(v)];
 
if volume == 0 % transmission channel with noise
  vol = 0;
else
  vol = 10^(volume/20);
end

pmin = int32(fs/350);
pmax = int32(fs/80);

prms(1) = prms(1) + vol*(2*(rand-0.5));

for i=2:order+1
  prms(i) = prms(i) + vol*(pi*rand);
  if prms(i) < 0 
      prms(i) = 0; 
  end
  if prms(i) > pi 
      prms(i) = pi; 
  end
end

pitchidx = order+2;
prms(pitchidx) = prms(pitchidx) + vol*(pmin-pmax)*rand;
if vol ~= 0
  if prms(pitchidx) > pmax 
      prms(pitchidx) = pmax; 
  end
  if prms(pitchidx) < pmin 
      prms(pitchidx) = pmin; 
  end
end

if vol ~= 0
  if rand > 0.5
    if prms(order+3) == 0
      prms(order+3) = 1;
      return;
    else
      prms(order+3) = 0;
    end  
  end 
end  

B = prms(1); % decode received parameters
lsfs = prms(2:order+1);
P = prms(order+2);
V = prms(order+3);
A = lsf2poly(lsfs);

excit = zeros(L,1); % presynthesize parameters
impidx = 1; 
imp = 1;
win=hamming(L);
for j=1:L
  if v==1
    excit(j) = imp;
    if mod(impidx,pitch) == 0
      imp = 1;
    else
      imp = 0;
    end
    impidx = impidx + 1;
  else
    excit(j) = 2*(rand-0.5);
  end
end
sout = filter(b0,a,excit);
sout = sout.*win;

  for b=1:L    % AddOverlap synthesis
    if n==1
        out(b) = tanh(sout(b));
        e(b) = excit(b); 
        g(b) = b0; 
        pt(b) = pitch;
    else
        out(n+b) = tanh(out(n+b) + sout(b));
        e(n+b) = excit(b); 
        g(n+b) = b0; 
        pt(n+b) = pitch;
    end
  end

  n=n+R;                            
  nframes=nframes+1;

end

% display parameters for report
clc;
insize = length(xin)*2; 
outsize = (order+3)*2;

fprintf(1,'input: \n');
fprintf(1,'Sampling rate (Hz): %d\n',fs);
fprintf(1,'Input length (samples): %d\n',length(xin));
fprintf(1,'Input length (seconds): %f\n',length(xin)/fs);
fprintf(1,'Channel noise volume (dB): %.2f\n\n',volume);

fprintf(1,'Frames: \n');
fprintf(1,'Analysis/synthesis frames: %d\n',nframes);
fprintf(1,'Frame length (ms): %d\n',fl);
fprintf(1,'Frame offset (ms): %d\n',fo);
fprintf(1,'Parameter length per frame (bytes): %d\n\n',outsize); % 16-bit

fprintf(1,'compression: \n');
fprintf(1,'Input size (bytes): %d\n',insize); % assumes 16-bit
fprintf(1,'Transmitted size (bytes): %d\n',nframes*outsize);
fprintf(1,'Compression ratio: %.2f\n',insize/(nframes*outsize));
fprintf(1,'Transmission bit rate (kbps): %.2f\n\n',(8/1024)*(nframes*outsize)/(length(xin)/fs));

subplot(4,1,1);
plot(xin)
axis([0 length(xin) -1 1]);
ylabel('volume'); title('Input Audio (filtered)');

subplot(4,1,2)
plot(pt,'m');
axis([0 length(xin) 0 150]);
ylabel('period (samples)'); title('Detected Pitch');
grid on;

subplot(4,1,3);
idx = (1:1:length(out));
zeroCheck = (pt==0);
e_noise = e; e_voice = e;
e_noise(~zeroCheck) = NaN;
e_voice(zeroCheck) = NaN;
plot(idx,e_voice,'g',idx,e_noise,'c',idx,g,'k');
legend('voiced','unvoiced','gain');
axis([0 length(xin) -1.25 1.25]);
ylabel('amplitude');
title('Excitation Signal and Gain');

subplot(4,1,4);
plot(out,'r')
axis([0 length(xin) -1 1]);
ylabel('volume'); title('Synthesized Output');
xlabel('time (samples)'); 
out = out/max(out); 
sound(out,16000);
%audiowrite('mycompress.wav',out,16000);