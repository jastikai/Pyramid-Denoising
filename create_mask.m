function mask = create_mask(input_image,k_limit)
    clear X Y Z K H k Pmax Pmin;
    imsize=size(input_image);
    X=repmat(1:imsize(2),imsize(2),1);
    Y=repmat(transpose(1:imsize(1)),1,imsize(1));
    Z=log(input_image);
    [K, H, Pmax, Pmin]=surfature(X,Y,Z);
    Pmax=real(Pmax); Pmin=real(Pmin);
    k=sqrt(Pmin.^2 + Pmax.^2);
    %create the first mask of 0 and 1 pixel values to identify spike-like
    %features such as stars and cosmic rays, and remove them by using the
    %mask created and by a local, pyramidal interpolation using a
    %nonlinear, multiresolution method
    mask=k; mask(k>k_limit)=0; mask(k<=k_limit)=1; mask(isnan(mask))=0;
end