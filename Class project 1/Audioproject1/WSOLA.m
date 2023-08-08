function [tsmwsola] = WSOLA(Y, Fs, Rate, shift)
    Hs = Fs * shift / 1000; 
    fl = Hs * 2;
    epstep = Hs * Rate;
    win = hanning(fl);
    
    wavelen = length(Y);
    tsmwsola = zeros(1, floor(wavelen / Rate));
    sp = Hs * 2; 
    rp = sp + Hs; 
    ep = sp + epstep; 
    outp = Hs;
    tsmwsola(1:outp) = Y(1:outp);
    spdata = Y(sp + 1 : sp + Hs) .* win(Hs + 1 : end);
    
    while wavelen > ep + fl
        areaa = Y(rp - Hs + 1 : rp + Hs); 
        areab = Y(ep - fl + 1 : ep + fl); 
        delta = smd(areaa, areab, fl, Hs);
        epd = ep + delta;
        spdata = Y(sp + 1 : sp + Hs) .* win(Hs + 1 : end);
        epdata = Y(epd - Hs + 1 : epd) .* win(1 : Hs);
        if length(spdata) == length(tsmwsola(outp : outp + Hs))
            tsmwsola(outp + 1 : outp + Hs) = spdata + epdata;
        else
            tsmwsola_len = length(tsmwsola(outp + 1 : outp + Hs));
            tsmwsola(outp + 1 : outp + Hs) = spdata(1 : tsmwsola_len) + epdata(1 : tsmwsola_len);
            outp = outp + Hs;
        end
        sp = epd;
        rp = sp + Hs;
        ep = ep + epstep;
    end
end
function [delta] = smd(ref, buff, fl, Hs)
    if length(ref) < fl
        ref = cat(1, ref, zeros(fl - length(ref), 1));
    end
    
    corr_max = 0;
    for i = 1 : 1 : length(buff) - fl
        compare = buff(i + 1 : i + fl);
        corr = sum(sum(ref .* compare));
        if corr > corr_max
            corr_max = corr;
            delta = i + Hs - fl;
        end
    end
end