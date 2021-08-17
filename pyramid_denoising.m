function denoised_image = pyramid_denoising(noisy_image)
    clear input_image imsize pixels denoised_image;
    input_image = noisy_image;
    denoised_image=input_image;
    means1=nanresizeim(input_image,0.5);
    means2=nanresizeim(means1,0.5);
    means3=nanresizeim(means2,0.5);
    synth1=imresize(means3,2,'bilinear');
    means2(isnan(means2))=synth1(isnan(means2));
    synth2=imresize(means2,2,'bilinear');
    means1(isnan(means1))=synth2(isnan(means1));
    synth3=imresize(means1,2,'bilinear');
    denoised_image(isnan(input_image))=synth3(isnan(input_image));
end