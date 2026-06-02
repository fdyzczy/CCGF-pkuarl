function [mu,Sigma] = CCTruncation(mu_0,Sigma_0,a,b,delta)
%CCTRUNCATION Single halfspace chance-constrained Gaussian truncation.

    colProb0 = 1/2*(1-erf((a.'*mu_0-b)/sqrt(2*a.'*Sigma_0*a)));
    if colProb0<=delta
        mu = mu_0;
        Sigma = Sigma_0;
        return;
    end

    q0 = a.'*Sigma_0*a;
    bt = (b-a'*mu_0)/sqrt(q0);
    c0=sqrt(2)*erfinv(1-2*delta);

    sigma = sqrt(c0^2*bt^2+4*c0^2+4)-c0*bt;
    sigma = sigma/2/(c0^2+1);
    lambda = -bt-c0*sigma;
    
    mu=mu_0-Sigma_0*a/sqrt(q0)*lambda;
    Sigma=Sigma_0-(1-sigma^2)*Sigma_0*a*a'*Sigma_0/(q0);
end
