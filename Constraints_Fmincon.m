function [c,ceq] = Ineq2(x)

% Cap_max = 100; %kWh
eta_charge = 0.85; eta_discharge = 0.85;
% remember that power is twice the specified value - this value
% corresponds to an energy transfer in a 30 minute period.
Charge_Power = 50; Discharge_Power = 50;
specified_charge = Charge_Power/2; specified_discharge = -Discharge_Power/2; 
charge_max = specified_charge*eta_charge; discharge_max = specified_discharge/eta_discharge;

% no value of energy transfer can be less than discharge_max
% no value of energy transfer can be greater than charge max
% 
% x >= dis_max
% x <= ch_max

c1 = x-charge_max;
ceq1 = zeros(size(x));

c1a = discharge_max-x;
ceq1a = zeros(size(x));

Cap_max = 100; %kWh
SOC = zeros(size(x));
SOC(1) = x(1);
for i = 2:length(x)
    SOC(i) = SOC(i-1)+x(i);
end

% SOC <= Cap_max
% SOC >= 0
    
c2 = SOC - Cap_max;
ceq2 = zeros(size(x));

c3 = -SOC;
ceq3 = zeros(size(x));

c = [c1; c1a; c2; c3];
ceq = [ceq1; ceq1a; ceq2; ceq3];

