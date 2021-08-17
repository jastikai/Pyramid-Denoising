clear t0 t1 dt;
t0=datenum([1996 1 1]);
t1=datenum([2015 12 31]);
dt=1;
% tag=[];
% exposure=[];
% im_time2=[];
% res_c2=[];
% for tm=t0:dt:t1;
%     clear hdr_txt
%     vector_date=datevec(tm);
%     year=vector_date(1);
%     month=vector_date(2);
%     day=vector_date(3);
%     try
%         hdr_txt=fileread(['/media/jaakko/CME Data/CME Data/C2 Images FITS/' num2str(year) '/' sprintf('%.2d',month) '/' sprintf('%.2d',day) '/img_hdr.txt']);
%     catch
%         continue
%     end
%     hdr_txt=strsplit(hdr_txt,'\n');
% %     warning('off')
% %     status = mkdir(['/media/jaakko/CME Data/CME Data/denoised/' num2str(year) '/' sprintf('%.2d',month) '/' sprintf('%.2d',day) '/']);
% %     warning('on')
%     for a=1:length(hdr_txt)-1
%         clear tmp aa bb;
%         tmp=hdr_txt{1,a};
%         tmp=strsplit(tmp,'  ');
% %     hdr_1{end+1,1}=tmp{1,1};
% %     hdr_1{end,2}=tmp{1,5};
%         aa=tmp{1,2};
%         bb=tmp{1,3};
%         im_time2(end+1)=datenum([str2num(aa(1:4)),str2num(aa(6:7)),str2num(aa(9:10)),str2num(bb(1:2)),str2num(bb(4:5)),str2num(bb(7:8))]);
%         exposure(end+1)=str2num(tmp{1,5});
%         res_c2(end+1,1)=str2num(tmp{1,6}); if res_c2(end,1)==0; res_c2(end,:)=[]; im_time2(end)=[]; continue; end
%         res_c2(end,2)=str2num(tmp{1,7});
%         tmp=strsplit(tmp{1,1},'.');
%         tmp=cell2mat(tmp);
%         tmp=tmp(1:end-3);
%         tag(end+1,1)=str2num(tmp);
%     end
% end
% im_time2=datevec(im_time2);
% clear temp_ind;
% temp_ind=find(res_c2(:,1)<1024 | res_c2(:,2)<1024);
% im_time2(temp_ind,:)=[];
% exposure(temp_ind)=[];
% res_c2(temp_ind,:)=[];
% save('im_time2')

for tm=t0:dt:t1;
    vector_date=datevec(tm);
    year=vector_date(1);
    month=vector_date(2);
    day=vector_date(3);
    warning('off')
    status = mkdir(['D:\CME Data\C2 Denoised FITS\' num2str(year) '\' sprintf('%.2d',month) '\' sprintf('%.2d',day) '\']);
    warning('on')
end
load('im_time','im_time')
clear temp_ind tm;
temp_ind=find(res_c2(:,1)<1024 | res_c2(:,1)~=res_c2(:,2));
im_time(temp_ind,:)=[];
exposure(temp_ind)=[];
tm=datenum(im_time);
clear loops;
loops=length(exposure);
clear roll_txt;
roll_txt=urlread('ftp://sohoftp.nascom.nasa.gov/pub/data/ancillary/attitude/roll/nominal_roll_attitude.dat');
roll_txt=strsplit(roll_txt,'\n');
roll_txt=roll_txt(6:end-1);
clear roll_time attitude;
roll_time=zeros(1,length(roll_txt));
for tt=1:length(roll_txt)
    clear temp;
    temp=roll_txt(tt);
    temp=temp{1};
    roll_time(tt)=datenum(temp(1:19));
    attitude(tt)=str2num(temp(20:end));
end
clear im_roll;
im_roll=zeros(1,length(tm));
for tt=1:length(roll_time)-1
    ind=find(tm>=roll_time(tt) & tm<roll_time(tt+1));
    im_roll(ind)=attitude(tt);
end
h=waitbar(0, 'Denoising LASCO C2 images...');
i=1;
circle=zeros(1024,1024);
for x=1:1024
for y=1:1024
if sqrt((x-512).^2 + (y-512).^2)>180
circle(x,y)=1;
end
end
end
clear starttime;
starttime=now;
for kk=279:loops
    clear im1 im2;
    try
        im1=fitsread(['F:\LASCO\C2 FITS\' num2str(im_time(kk,1)) '\' sprintf('%.2d',im_time(kk,2)) '\' sprintf('%.2d',im_time(kk,3)) '\' num2str(im_time(kk,1)) sprintf('%.2d',im_time(kk,2)) sprintf('%.2d',im_time(kk,3)) '_' sprintf('%.2d',im_time(kk,4)) '_' sprintf('%.2d',im_time(kk,5)) '_' sprintf('%.2d',im_time(kk,6)) '.fits']);
    catch
        warning(['F:\LASCO\C2 FITS\' num2str(im_time(kk,1)) '\' sprintf('%.2d',im_time(kk,2)) '\' sprintf('%.2d',im_time(kk,3)) '\' num2str(im_time(kk,1)) sprintf('%.2d',im_time(kk,2)) sprintf('%.2d',im_time(kk,3)) '_' sprintf('%.2d',im_time(kk,4)) '_' sprintf('%.2d',im_time(kk,5)) '_' sprintf('%.2d',im_time(kk,6)) '.fits' ' did not read']);
        im1=zeros(1024,1024);
    end
    if im_roll(kk)==180
        im1=imrotate(im1,180);
    end
    clear bival_mask;
    bival_mask=im1; bival_mask(im1~=0)=1;
    im2=im1;
    im1=im1.*circle;
    clear ex2 ex1;
    ex1=exposure(kk);
    %normalise by exposure time
    clear norm_im1 norm_im2;
    norm_im1=im1./ex1;
    im2=im2./ex1;
    %create the first mask of 0 and 1 pixel values to identify spike-like
    %features such as stars and cosmic rays, and remove them by using the
    %mask created and by a local, pyramidal interpolation using a
    %nonlinear, multiresolution method
    clear mask1 noisy_im1 denoised1;
    mask1=create_mask(norm_im1,0.09);
    noisy_im1=norm_im1.*mask1; noisy_im1(noisy_im1==0)=NaN;
    denoised1=pyramid_denoising2(noisy_im1); denoised1(isnan(denoised1))=0;
    %
    %
    %apply the median filter to the original input images
    %again for the removal of bigger features such as comets and planets,
    %then apply morphological dilation to the mask of invalid and valid
    %pixels
    clear im_medfilt2 mask2 se morph_dil;
    im_medfilt2=medfilt2(im2,[4 4]);
    mask2=create_mask(im_medfilt2,0.25);
    se=strel('square',10);
    morph_dil=imerode(mask2,se);
    %
    %
    %
    denoised1=denoised1.*circle;
    clear bad_pixels;
    bad_pixels=(length(find((bival_mask.*morph_dil)==0))/(1024^2))*100;
    if bad_pixels >= 50
        final_im=zeros(1024,1024);
    else
        final_im=denoised1.*morph_dil;
    end
%     clear X; X=pcolor(final_im); set(X,'EdgeColor','None');
%     title(datestr(im_time(kk,:)));
    fitswrite(final_im,['D:\CME Data\C2 Denoised FITS\' num2str(im_time(kk,1)) '\' sprintf('%.2d',im_time(kk,2)) '\' sprintf('%.2d',im_time(kk,3)) '\' num2str(im_time(kk,1)) sprintf('%.2d',im_time(kk,2)) sprintf('%.2d',im_time(kk,3)) '_' sprintf('%.2d',im_time(kk,4)) '_' sprintf('%.2d',im_time(kk,5)) '_' sprintf('%.2d',im_time(kk,6)) '.fits'],'compression','gzip');
    clear timepassed timetotal timeleft;
    timepassed=now-starttime;
    timetotal=timepassed*loops/i;
    timeleft=timetotal-timepassed;
    timeleft=datevec(timeleft);
    waitbar(kk/loops,h,['Denoising ... ' datestr(im_time(kk,:)) ' / ' num2str(timeleft(3)) 'd ' num2str(timeleft(4)) 'h ' num2str(timeleft(5)) 'min ' sprintf('%.0f',timeleft(6)) 's remaining.']);
    i=i+1;
end

close(h);
% fig=figure;
% movie(fig,F,1,5);
% close(fig);