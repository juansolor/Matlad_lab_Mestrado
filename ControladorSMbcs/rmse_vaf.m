function [E,VAF] = rmse_vaf(x,X_k,M)%calcula RMSE e VAF de cada saída estimada

VAF = ((1-var(x-X_k)./var(x))*100)';

E =(sqrt(sum((X_k-x).^2)/(M-1)))';%sqrt se quiser rmse
end