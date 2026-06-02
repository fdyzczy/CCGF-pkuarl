function w = weight_alloc(mu0,mu,Sigma0,Sigma,w0)
%WEIGHT_ALLOC Update Gaussian-mixture weights after truncation.
    mixNum = length(w0);
    optKLD = zeros(1, mixNum);
    weight = zeros(1, mixNum);
    for i=1:mixNum
        optKLD(i) = gaussKLD(mu(:,i),mu0(:,i),Sigma(:,:,i),Sigma0(:,:,i));
        weight(i) = w0(i)*exp(-optKLD(i));
    end

    w = weight/sum(weight);
    if any(isnan(w)) || any(~isfinite(w))
        w = w0;
    end

end
