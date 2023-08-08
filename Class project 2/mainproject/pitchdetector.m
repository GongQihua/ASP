function [pitch,v] = pitchdetector(xw,pitchmin,pitchmax)

kmin = pitchmin;
kmax = pitchmax;
Rn = zeros(length(xw),1); 

max = xw(1); %find max value, the pitch
for i=1:length(xw)
  if xw(i) > max
    max = xw(i);
  end
end
CL = 0.3*max;   

for i=1:length(xw) % clipping
  if xw(i) > CL
    xw(i) = 1;
  else
    if xw(i) < -CL
      xw(i) = -1;
    else
      xw(i) = 0;
    end  
  end
end

for k=kmin:kmax
  for m=1:length(xw)-k
    Rn(k) = Rn(k) + xw(m)*xw(m+k);
  end
end  

R0 = 0; %find Rn(0)
for m=1:length(xw)
  R0 = R0 + xw(m)*xw(m);
end

Rmax = Rn(1); kidx = 1;% find max R and k
for k=1:length(Rn)
  if Rn(k) > Rmax
    Rmax = Rn(k);
    kidx = k + kmin - 1; 
  end
end

if Rmax > 0.1*R0  % voiced or unvoiced detect
  v = 1;
  pitch = kidx;
else
  v = 0;
  pitch = 0;
end

end 

