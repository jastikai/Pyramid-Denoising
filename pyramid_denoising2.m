function denoised_image = pyramid_denoising2(noisy_image)
    clear input_image imsize pixels denoised_image means;
    input_image = noisy_image;
%     denoised_image=input_image;
    imsize=size(input_image);
    circle=zeros(imsize(1),imsize(2));
    for x=1:imsize(1)
        for y=1:imsize(2)
            if sqrt((x-(imsize(1)/2)).^2 + (y-imsize(2)/2).^2) >180
                circle(x,y)=1;
            end
        end
    end
    input_image(isnan(input_image)==1 & circle==0)=0;
    bad_pixels=length(find(isnan(input_image)==1))/length(input_image)^2*100;
    means={};
    means{end+1}=nanresizeim(input_image,0.5);
    while bad_pixels > 0 & length(means{end})>4
        bad_pixels=length(find(isnan(means{end})==1))/length(means{end})^2*100;
        means{end+1}=nanresizeim(means{end},0.5);
    end
    clear synths;
    synths={};
    synths{end+1}=imresize(means{end},2,'bilinear');

    while length(synths{end})<1024
        clear varmean varsynth;
        means(end)=[];
        varmean=means{end};
        varsynth=synths{end};
        varmean(isnan(varmean))=varsynth(isnan(varmean));
        means{end}=varmean;
        synths{end+1}=imresize(means{end},2,'bilinear');
    end
    clear varsynth;
    varsynth=synths{end};
    input_image(isnan(input_image))=varsynth(isnan(input_image));
    denoised_image=input_image;
end
