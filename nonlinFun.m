function [f,G] = nonlinFun(gamma, B,tau)

Bgamma = B * gamma;

f = sum((1./(Bgamma) - 1./tau).^2);

G = zeros(size(gamma));
for i = 1 : length(gamma)
    
    s = 0;
    for k = 1 : length(tau)
        s = s + (1 / Bgamma(k) - 1./tau(k)) * (-B(k,i) / Bgamma(k)^2);
    end
    
    G(i) = 2 * s;
    
end

end

