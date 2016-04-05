% optimisation of energy storage schedule to maximise revenue for a device with a fixed round trip efficiency.
% Inspired by the paper "Practical operation strategies for pumped
% hydroelectric energy storage (PHES) utilising electricity price
% arbitrage" by Connelly et al.

% See the paper appendix for algorithm description and flowchart

% Code written by Edward Barbour
% e.r.barbour@bham.ac.uk

clc; clear all;
% clear variables and screen 
price = textread('Price_timeseries.txt', '%f');
% load the data to be used
% In the algorithm I use "price", hence I assign "price" to the relevant
% variable in the loaded data
% load spot_price_2013;
% price = spot_price_2013/1000;

% use prices in £/kWh for direct conversion

% The storage device parameters
Cap_max = 100; %kWh
eta_charge = 0.85; eta_discharge = 0.85;

% the device has a round trip of 72.25%, hence assign both the charging and
% discharging processes an individual efficiency of sqrt(72.25%=0.85).

% remember that power is twice the specified value - this value
% corresponds to an energy transfer in a 30 minute period.
Charge_Power = 50; Discharge_Power = 50;
specified_charge = Charge_Power/2; specified_discharge = -Discharge_Power/2; 
charge_max = specified_charge*eta_charge; discharge_max = specified_discharge/eta_discharge;
clear specified_charge specified_discharge Charge_Power Discharge_Power;

% It should be noted when thinking about charging and discharging powers
% that the specified charge and discharge can each have two meanings:
% either the power that can be taken from the grid and power that can be 
% returned to the grid or the power that can be transferred to the storage
% device and the power that can be drawn from the storage device.

% here the specified_charge and specified_discharge mean power
% taken from grid and power returned to grid respectively

% Counter for the amount of transactions
counter = 0;

% arrays to store energy stored and energy transfer for time period
energy_stored = zeros(1,length(price)); energy_transfer = zeros(1,length(price));

remove = ones(1,length(price));
% this is a marker for each period in the timeseries, remove = 1 means the
% period hasn't yet been removed

while(any(remove)==1)
    
% find maxhour in the price distribution
matrix(1,:) = find(remove==1);
matrix(2,:) = price(remove==1);
index = find(matrix(2,:)==max(matrix(2,:)),1,'first');
maxh = matrix(1,index);
clear matrix index;
% find the last hour before maxh when storage was full
r1 = find(energy_stored(1:maxh)==Cap_max,1,'last')+1;
if(isempty(r1)==1)
    r1 = 1;
end

% find the first hour after maxh when storage is empty
r2 = find(energy_stored(maxh:length(energy_stored))==0,1,'first')+(maxh-2);
if(isempty(r2)==1)
    r2 = length(energy_stored);
end

% find minh in the time range
range_price = price(r1:r2);
range_remove = remove(r1:r2);

% if there is no hour in the range that hasn't been removed then remove
% maxh and skip to the end
if(isempty(range_price(range_remove==1))==1)
    remove(maxh)=0;
else
%     minh = find(range_price==min(range_price(range_remove==1)),1,'first')+(r1-1);
    matrix(1,:) = find(range_remove==1);
    matrix(2,:) = range_price(range_remove==1);
    index = find(matrix(2,:)==min(matrix(2,:)),1,'first');
    minh = matrix(1,index)+r1-1;
    clear matrix index;
   
    % Calculate the marginal operating cost for charging at minh and
    % discharging at maxh
    MoC = price(minh)/(eta_charge*eta_discharge);

    if(MoC<price(maxh) && minh~=maxh)
    
        bottleneck = zeros(1,3);
        % Calculate the storage operation bottlenecks
        bottleneck(1) = energy_transfer(maxh)-discharge_max;
        % bottleneck(1) represents the discharging capability
        bottleneck(2) = charge_max - energy_transfer(minh);
        % bottleneck(2) represents the charging capability
        if maxh>minh
            bottleneck(3) = Cap_max - max(energy_stored(minh:maxh));
        else
            bottleneck(3) = min(energy_stored(maxh:minh));
        end
        % bottleneck(3) represents either the spare storage capacity
        % between minh and maxh or the spare energy stored between maxh and
        % minh
        actual_bottleneck = min(bottleneck);

        % update the storage schedule
        energy_transfer(minh) = energy_transfer(minh) + actual_bottleneck;
        energy_transfer(maxh) = energy_transfer(maxh) - actual_bottleneck;
        if maxh>minh
            energy_stored(minh:maxh-1) = energy_stored(minh:maxh-1) + actual_bottleneck;
        else
            energy_stored(maxh:minh-1) = energy_stored(maxh:minh-1) - actual_bottleneck;
        end

        % check if at the charge or discharge operation is at capacity at either
        % maxh or minh and remove that hour from the price distribution.
    
        if(energy_transfer(maxh)<=discharge_max)
            remove(maxh)=0;
        end
        if(energy_transfer(minh)>=charge_max)
            remove(minh)=0;
        end
    

    else
        % if there is no price incentive then simply remove the hours
        remove(maxh) = 0;
        remove(minh) = 0;
    end
    
end

end

energy_input = zeros(1,length(energy_transfer));
energy_input(energy_transfer<0) = energy_transfer(energy_transfer<0)*eta_discharge;
energy_input(energy_transfer>0) = energy_transfer(energy_transfer>0)/eta_charge;

% Plotting the energy stored, the energy transfer into and out of the
% stores and the energy taken and returned from the network (called 'energy_input')

time = linspace(1,length(energy_stored),length(energy_stored))/2;
figure; 
plot(time, energy_stored, time, 2*energy_transfer, time, 2*energy_input, time, price*1000)
xlabel('time (hrs)')
ylabel('energy-stored/energy-transfer (kWh/kW)')
legend('energy stored', 'energy transfer', 'energy input', 'price (£/MWh)')
title('Illustrating the optimum schedule of storage and charging/discharging')

rev = zeros(1, length(price));
for i = 1:length(price)
    if(energy_transfer(i)>0)
        rev(i) = energy_transfer(i)*price(i)/eta_charge;
    else
        rev(i) = energy_transfer(i)*price(i)*eta_discharge;
    end
end
fprintf('%f \n', -sum(rev));

mycall1 = num2str(-sum(rev));
message = sprintf(['Revenue = £' mycall1]);
uiwait(msgbox(message));