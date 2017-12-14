function [tot] = RevFunc(x, price, eta_charge, eta_discharge)

rev = zeros(1, length(price));
for i = 1:length(price)
    if(x(i)>0)
        rev(i) = x(i)*price(i)/eta_charge;
    else
        rev(i) = x(i)*price(i)*eta_discharge;
    end
end

tot = sum(rev);