

function[filt_pos] = zb_butterworth(pos,CuttFreq,SamplFreq)

[b,a]=butter(2, 2*(CuttFreq./SamplFreq),'low');
filt_pos = filtfilt(b, a, pos);

end