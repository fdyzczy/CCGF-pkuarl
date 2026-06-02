
function KLD = gaussKLD(mu1,mu2,Sigma1,Sigma2)
%GAUSSKLD Gaussian KL divergence KLD(gm1 || gm2).
    mu1 = mu1';
    mu2 = mu2';
    invSigma2 = inv(Sigma2);
    KLD = 1/2*(log(det(Sigma2))-log(det(Sigma1))-length(mu1)+trace(invSigma2*Sigma1)+(mu2-mu1)*invSigma2*(mu2-mu1)');
end
