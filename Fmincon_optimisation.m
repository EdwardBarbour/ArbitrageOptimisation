% Trying with fmincon
clc; clear all;

% load price
price = textread('Price_timeseries.txt', '%f');

x_0 = zeros(size(price));

% set efficiencies here
eta_charge = 0.85;
eta_discharge = 0.85;
%%

tot = Revenue_Func(x_0, price, eta_charge, eta_discharge);

%%

options = optimset('Disp','iter','Algorithm','interior-point',...
                   'MaxFunEvals',100000,'MaxIter',5000,'TolX',1e-10,'TolFun',1e-6,'TolCon',1e-6);

tic
               
[x, fval] = fmincon(@(x)Revenue_Func(x, price, eta_charge, eta_discharge), x_0, [], [], [], [], [], [], @Constraints_Fmincon, options);

toc
%%

SOC = zeros(size(x));
SOC(1) = x(1);
for i = 2:length(x)
    SOC(i) = SOC(i-1)+x(i);
end

rev = zeros(1, length(price));
for i = 1:length(price)
    if(x(i)>0)
        rev(i) = x(i)*price(i)/eta_charge;
    else
        rev(i) = x(i)*price(i)*eta_discharge;
    end
end

energy_stored = SOC;
energy_transfer = x;
energy_input = zeros(1,length(energy_transfer));
energy_input(energy_transfer<0) = energy_transfer(energy_transfer<0)*eta_discharge;
energy_input(energy_transfer>0) = energy_transfer(energy_transfer>0)/eta_charge;

time = linspace(1,length(energy_stored),length(energy_stored))/2;
figure; 
strRev = ['Revenue = ',num2str(-sum(rev))];
plot(time, energy_stored, time, 2*energy_transfer, time, 2*energy_input, time, price*1000)
text(1,-30,strRev,'HorizontalAlignment','left','fontsize',16);
xlabel('time (hrs)')
ylabel('energy-stored/energy-transfer (kWh/kW)')
legend('energy stored', 'energy transfer', 'energy input', 'price (£/MWh)')
title('Illustrating the optimum schedule of storage and charging/discharging')

