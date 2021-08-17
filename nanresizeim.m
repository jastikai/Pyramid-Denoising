function resized_im = nanresizeim(input_image,ratio)

clear pic2 pic3
pic2 = input_image;
fac = ratio;
dx = 1./fac;
[r,c] = ndgrid(1:size(pic2,1), 1:size(pic2,2));
[n, ibin] = histc(r(:), 0.5:dx:size(pic2,1)+0.5);
[n, jbin] = histc(c(:), 0.5:dx:size(pic2,2)+0.5);
nr = max(ibin);
nc = max(jbin);
idx = sub2ind([nr nc], ibin, jbin); 
pic3 = accumarray(idx, pic2(:), [nr*nc 1], @nanmean);
pic3 = reshape(pic3, nr, nc);
resized_im = pic3;
end
