function epsilon = POD_getepsilon(s,N)       
    % Computes the amount of energy of the neglected modes (greater than
    % N), given the singular values 's'.

    epsilon = sqrt(sum(s(N+1:end).^2)/sum(s.^2));
%     sigma2_cum = cumsum(s.^2);
%     epsilon = sqrt(1-sigma2_cum(N)./sigma2_cum(end)); 
end