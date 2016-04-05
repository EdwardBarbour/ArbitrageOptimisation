% Code for optimising the schedule of operation of a storage device to
% maximise its revenue via arbitrage. Used in the academic paper "Towards an
% objective method to compare energy storage technologies: development and 
% validation of a model to determine the upper boundary of revenue available 
% from electrical price arbitrage", E. barbour et al. A flow chart and
% description of the algorithms operation can be found in this paper. NOTE
% the algorithm has been extended from the version used in this paper

% Contributors Edward Barbour*, Grant Wilson and Simon Gill
% e.r.barbour@bham.ac.uk
% code written by Edward Barbour

% this program has the capability to model fixed charging and discharging
% efficiencies and time dependent losses.

% It was extended and used in the work "Maximising revenue for non-firm distributed wind
% generation with energy storage in an active management scheme", S. Gill
% et al. A flowchart and description of the algorithms extended operation
% can be found in this paper.

% MODEL INPUTS:
% Price timeseries
% Renewable power timeseries (kW) (optional). 
% Import and Export Constraints (kW) (optional)
% Storage device characteristics


% the Renewable power series 'REC_power' and the constraint files are
% optional. 

%Reset loop, set reset = 1 to continue a previous simulation. Generally use
% 0
reset = 0;

if(reset == 0)
    clear all; clc;
    scenario=0;

% user friendly code for reading in input timeseries using a .mat file    
% [FileName,PathName] = uigetfile('*.mat','Select the MATLAB variable file');
% 
% load(FileName) 

price = textread('Price_timeseries.txt', '%f');
% and then create a time timeseries
% load spot_price_2013.mat;
% price = spot_price_2013';
time = linspace(1,length(price),length(price))/2;

    % price should be £/kWh later
    %
    % REC_power, constraint, constraint_import should all be in kW
    % demand in kW, though constraint needs set up first

    % price denoted 'price'. Hence assign price to the appropriate loaded
    % variable
    
if exist('price', 'var')
    % assign the maximum length that the alogorithm has to work with
    index_ref = length(price);
    prices = price;
else
    error('price data absent')
%   Price info does not exist.
end

if exist('REC_power', 'var')
    if(length(REC_power)~=index_ref)
        warningMessage = sprintf('Warning: REC_power length being corrected');
        uiwait(msgbox(warningMessage));
        REC_power = REC_power(1:index_ref);
    end
    power = REC_power;
else
    power(1:index_ref)=0;
end

if exist('constraint_export', 'var')
    if(length(constraint_export)~=index_ref)
        warningMessage = sprintf('Warning: constraint length being corrected');
        uiwait(msgbox(warningMessage));
        constraint_export = constraint_export(1:index_ref);
    end
else
    constraint_export(1:index_ref) = inf;
end
% If a constraint file is not included then the export constraint is set at
% infinity - i.e. no constraint

if exist('constraint_import', 'var')
    if(length(constraint_import)~=index_ref)
        warningMessage = sprintf('Warning: constraint_import length being corrected');
        uiwait(msgbox(warningMessage));
        constraint_import = constraint_import(1:index_ref);
    end
else
    constraint_import(1:index_ref) = -inf;
end
% If a constraint file is not included then the import constraint is set at
% infinity - i.e. no constraint
        
    % search window characteristics, depend only on length of the price
    % file, no need for user action
    index_1 = 1;
    index_2 = index_ref;
    min_period = index_1; max_period = index_2;
    window = max_period - min_period + 1; period_1 = 0; period_2 = 0;
    max_array_size = max_period - min_period + 1;
    
    %Storage device characteristics

prompt={'Capcity (kWh)', 'Charging Power (kW)', 'Discharge Power (kW)' ...
, 'charging efficiency (%)', 'discharging efficiency (%)', 'self discharge (hrs)'};
defans={'100', '50', '50', '0.85', '0.85', 'inf'};
fields = {'I1', 'I2', 'I3', 'I4', 'I5', 'I6'};
info = inputdlg(prompt, 'Please input storage device characteristics', 1, defans);
if ~isempty(info);            %see if user hit cancel
   info = cell2struct(info,fields);
   
   Cap_max = str2num(info.I1);   %#ok<ST2NM> %convert input strings to numbers, so that
   PLI = str2num(info.I2);   %#ok<ST2NM> %so that they can be used in program
   PLO = str2num(info.I3); %#ok<ST2NM>
   eta_in = str2num(info.I4); %#ok<ST2NM>
   eta_out = str2num(info.I5); %#ok<ST2NM>
   tau = str2num(info.I6); %#ok<ST2NM>
      
   mycall1 = num2str(Cap_max);   mycall2 = num2str(PLI);   %Convert back to strings to tell the
   mycall3 = num2str(PLO);   mycall4 = num2str(eta_in);   %user the inputs
   mycall5 = num2str(eta_out);
   mycall6 = num2str(tau);
   
   
   uiwait(msgbox(sprintf(['Maximum Storage Capcity =' mycall1 ' kWh \n', ...
       'Charge Power = ' mycall2 ' kW \n', ...
       'Discharge Power = ' mycall3 ' kW \n', ...
       'Charging Efficiency  = ' mycall4 '  \n', ...
       'Discharging Efficiency = ' mycall5 '  \n', ...
       'time constant = ' mycall6 ' hours' ]), 'Storage Parameters'))

end


% It should be noted when thinking about charging and discharging powers
% that the specified charge and discharge can each have two meanings:
% either the power that can be taken from the grid and power that can be 
% returned to the grid or the power that can be transferred to the storage
% device and the power that can be drawn from the storage device.

% here the specified charge and discharge ('PLI' and 'PLO') mean power
% taken from grid and power returned to grid respectively

    % convert powers to energy per 30 min period
    ELI = (PLI*eta_in)/2; %kJ 
    ELO = -(PLO/eta_out)/2; %kJ

    % Set curtailed energy equal to zero for all times
    free_energy(1:max_array_size)=0;
    

    %algorithm parameters- govern length and accuracy of the search
prompt={'Number of iterations', 'Nearest Neighbour bias'};
fields = {'I1', 'I2'};
defans={'100000', '100000'};
info = inputdlg(prompt, 'Please input search characteristics', 1, defans);
if ~isempty(info);            %see if user hit cancel
   info = cell2struct(info,fields);
   
   n_max = str2double(info.I1);   %convert input strings to numbers, so that
   nn_bias = str2double(info.I2);   %so that they can be used in program
   
   mycall1 = num2str(n_max);   mycall2 = num2str(nn_bias);
   
   uiwait(msgbox(sprintf(['Maximum iterations =' mycall1 '\n', ...
       'Nearest neighbour bias = ' mycall2 '\n']), 'Seacrh Parameters'))

end  
    %     iteration = zeros(n_max, 6);
    accepted=0;
    step_size = 1;
   
    
    %Set up the power input to the simulation
    old_test_energy = power(index_1:index_2)/2; % convert renewable energy power to energy per 30 min period
    test_price = prices'; %Convert prices to £/kWh, they should already be in this form
%     time = t(index_1:index_2);
    constraint_2 = constraint_import/2; % convert power constraints to energy per 1/2 hour
    constraint = constraint_export/2;
    
    %Initialise the variables we want to have the ability to reset
    E_stored(1:max_array_size)=0; Revenue(1:max_array_size)=0;
    E_to_store(1:max_array_size)=0; L(1:max_array_size)=0; new_test_energy(1:max_array_size)=0;
    f_used(1:max_array_size)=0;
    
    %these are for seeing how revenue changes with iterations
    iterations = linspace(1,n_max,n_max);
    test_Rev = zeros(1, n_max);
     
    %setting up the free (curtailed) energy and revenue with spillage
    for x = min_period:max_period;
        if(old_test_energy(x)>constraint(x));
            free_energy(x) = old_test_energy(x) - constraint(x);
            %the free energy is the difference
        end;
        new_test_energy(x) = old_test_energy(x)- free_energy(x);
        Revenue(x) = new_test_energy(x) * test_price(x);
    end
    
    
% 'free_energy' is the amount of renewables that would be curtailed due to the export constraint
% this is available to energy storage at zero cost
    free_energy_old = free_energy;
    Output_to_grid = new_test_energy;
    old_Rev = sum(new_test_energy*test_price');
    
% old_Rev is the revenue that would be generated by selling the renewable energy at the time that 
% occurs.
    
end %the end for the reset values command

test_factor1 = eta_in*eta_out*exp((min(time)-max(time))/tau);

% % plotting a graph of curtailed energy, energy generated by REC, energy generated
% % with no constraints
% figure; plot(time, old_test_energy, time, free_energy_old)
% xlabel('time (hrs)')
% ylabel('renewable power (kW)')
% legend('Total Renewable Generation', 'Curtailed Renewables')


for n = 1:n_max
    
    %set a variable to accept changes or decline them appropriately
    accept=0; %default is to decline
    %choose amount to shift, either positive or negative
    %+ve and -ve allows the program to correct for mistakes
    dx = step_size - rand(1)*step_size*2;
    %dx = -10000;
    %select a random period
    period_1 = round(1 + (rand(1)*(window-1)));

    %biasing towards near neighbours and making sure period_1 ~= period_2
    while(accept == 0)
        period_2 = round(1 + (rand(1)*(window-1)));
        prob_accept = exp((-(period_1 - period_2)^2)/nn_bias);
        if(rand(1)<prob_accept && period_1 ~= period_2)
            accept = 1;
        end
    end


    
    %set the default back to decline
    accept = 0; %#ok<NASGU>
    %making sure period_1 is smaller
    if(period_2<period_1)
        temp = period_1; period_1 = period_2; period_2 = temp;
    end
% If want to check action at specific periods do so here    
    time_loss = exp((time(period_1)-time(period_2))/tau);
    
    if dx>=0
        dx = ELI;
         if (free_energy(period_1) > 0)
            
            
            %scenario = 9
            if(E_to_store(period_2)<=0)
                scenario = 9;
                dx4= free_energy(period_1) * eta_in;
                dx1 = (ELI-E_to_store(period_1));
                dx2 = (E_to_store(period_2)-ELO)/time_loss;
                dx5=(constraint(period_2)-Output_to_grid(period_2))/(eta_out*time_loss);
                dx = min ( [dx dx1 dx2 dx4 dx5] );
            end
            
            %scenario = 10
            if(E_to_store(period_2)>0)
                scenario = 10;
                dx4= free_energy(period_1) * eta_in;
                dx1= ELI - E_to_store(period_1);
                dx2 = E_to_store(period_2)/time_loss;
                dx5 = (constraint(period_2)-Output_to_grid(period_2))*(eta_in/time_loss);
                dx = min ( [dx dx1 dx2 dx4 dx5] );
            end
            
            i_temp = find(E_stored==max(E_stored(period_1:period_2-1)), 1, 'first');
                    
            time_loss_3 = exp((time(period_1)-time(i_temp))/tau);
            m2 = Cap_max - max(E_stored(period_1:(period_2 - 1)));
            dx3 = (m2/time_loss_3)-(1e-12);
                                
            if(dx3<0); dx3 = 0; end
            
            dx = min( [dx dx3] );
            
            %no price check required as long as prices are always positive
            %(and non zero)
          
                accept=1;
                
            %No checks on in and outflows required - just on interim
            %capacity
            
           
            
        else %no free energy at 1    
            %Reduce size of dx to increase chance of success
            
            %scenario 1
            if(E_to_store(period_1)>=0 && E_to_store(period_2)>0)
                scenario=1;
                dx1 = (ELI-E_to_store(period_1));
                dx2 = E_to_store(period_2)/time_loss;
                dx5 = (constraint(period_2)-Output_to_grid(period_2))*(eta_in/time_loss);
                %new constraint
                dx_test = (Output_to_grid(period_1)-constraint_2(period_1))*eta_in;
                dx = min ( [dx dx1 dx2 dx5 dx_test] );
            end
            
            %scenario 2
            if(E_to_store(period_1)>=0 && E_to_store(period_2)<=0)
                scenario=2;
                dx1 = (ELI-E_to_store(period_1));
                dx2 = (E_to_store(period_2)-ELO)/time_loss;
                dx5 = (constraint(period_2)-Output_to_grid(period_2))/(eta_out*time_loss);
                %new constraint
                dx_test = (Output_to_grid(period_1)-constraint_2(period_1))*eta_in;
                dx = min ( [dx dx1 dx2 dx5 dx_test] );
            end
            
            %scenario 3
            if(E_to_store(period_1)<0 && E_to_store(period_2)>0)
                scenario=3;
                dx1 = -E_to_store(period_1);
                dx2 = E_to_store(period_2)/time_loss;
                dx5 = (constraint(period_2)-Output_to_grid(period_2))*(eta_in/time_loss);
                dx = min ( [dx dx1 dx2 dx5] );
            end
            
            %scenario 4
            if(E_to_store(period_1)<0 && E_to_store(period_2)<=0)
                scenario=4;
                dx1 = -E_to_store(period_1);
                dx2 = (E_to_store(period_2)-ELO)/time_loss;
                dx5 = (constraint(period_2)-Output_to_grid(period_2))/(eta_out*time_loss);
                dx = min ( [dx dx1 dx2 dx5] );
            end
            
            i_temp = find(E_stored==max(E_stored(period_1:period_2-1)), 1, 'first');
                    
            time_loss_3 = exp((time(period_1)-time(i_temp))/tau);
            m2 = Cap_max - max(E_stored(period_1:(period_2 - 1)));
            dx3 = (m2/time_loss_3)-(1e-12);
                                
            if(dx3<0); dx3 = 0; end
            
            dx = min( [dx dx3] );
            
            %The price checkers - accept if revenue increases
            %Works out the values Output_to_grid will take is move is
            %accepted
            
            new_EtS_1 = E_to_store(period_1) + dx; 
            new_EtS_2 = E_to_store(period_2) - dx*time_loss;
            
            if(scenario==1)
                new_OtG_1=Output_to_grid(period_1)-dx/eta_in;
                new_OtG_2=Output_to_grid(period_2)+(dx*time_loss)/eta_in;
            end
            if(scenario==2)
                new_OtG_1=Output_to_grid(period_1)-dx/eta_in;
                new_OtG_2=Output_to_grid(period_2)+(dx*time_loss*eta_out);
            end
            if(scenario==3)
                new_OtG_1=Output_to_grid(period_1)-dx*eta_out;
                new_OtG_2=Output_to_grid(period_2)+(dx*time_loss)/eta_in;
            end            
            if(scenario==4)
                new_OtG_1=Output_to_grid(period_1)-dx*eta_out;
                new_OtG_2=Output_to_grid(period_2)+(dx*time_loss*eta_out);
            end            
            
            Rev1=Output_to_grid(period_1)*test_price(period_1) + Output_to_grid(period_2)*test_price(period_2);
            Rev2=new_OtG_1*test_price(period_1) + new_OtG_2*test_price(period_2);
            
            %Checks if the 'new' revenue is higher than the old and accepts
            %if it is. 
            
            accept = 0;
            if(new_OtG_1<=constraint(period_1) && new_OtG_2<=constraint(period_1))
            if(Rev2 > Rev1 + 1e-14) %temporary correction...
                accept=1;
            end
            end
            
        end
        
    else %if dx < 0 
        dx = ELO;
        if(free_energy(period_2)> 0)

            %scenario = 11
            if(E_to_store(period_1)<=0)
            scenario=11; 
            dx4= -(free_energy(period_2) * eta_in)/time_loss;
            dx2= -(ELI - E_to_store(period_2))/time_loss;
            dx1= (ELO - E_to_store(period_1));
            dx5= -(constraint(period_1)-Output_to_grid(period_1))/(eta_out);
            dx = max ( [dx dx1 dx2 dx4 dx5] );
            end

            %scenario = 12
            if(E_to_store(period_1)>0)
            scenario=12; 
            dx4= -(free_energy(period_2) * eta_in)/time_loss;
            dx2= -(ELI - E_to_store(period_2))/time_loss;
            dx1 = -E_to_store(period_1);
            dx5 = -(constraint(period_1)-Output_to_grid(period_1))*eta_in;
            dx = max ( [dx dx1 dx2 dx4 dx5] );
            end
            
            i_temp = find(E_stored==min(E_stored(period_1:period_2-1)), 1, 'first');
                    
            time_loss_2 = exp((time(period_1)-time(i_temp))/tau);
            m = min(E_stored(period_1:(period_2 - 1)));
            dx3 = (-m/time_loss_2)+(1e-12);
                                
            if(dx3>0); dx3 = 0; end
            
            dx = max( [dx dx3] );

            %no price check required as long as prices are always positive
            %(and non zero)
          
            accept = 1;
            
            
        else % if free_energy(period_2) = 0

            %increase dx to increase acceptance
            
            %scenario 5
            if(E_to_store(period_1)>0 && E_to_store(period_2)>=0)
                scenario=5;
                dx1 = -E_to_store(period_1);
                dx2 = (E_to_store(period_2)-ELI)/time_loss;
                dx5 = -(constraint(period_1)-Output_to_grid(period_1))*eta_in;
                %new constraint
                dx_test = (constraint_2(period_2) - Output_to_grid(period_2))*(eta_in/time_loss);
                dx = max ( [dx dx1 dx2 dx5 dx_test] );
            end
            
            %scenario 6
            if(E_to_store(period_1)>0 && E_to_store(period_2)<0)
                scenario=6;
                dx1 = -E_to_store(period_1);
                dx2 = E_to_store(period_2)/time_loss;
                dx5 = -(constraint(period_1)-Output_to_grid(period_1))*eta_in;
                
                dx = max ( [dx dx1 dx2 dx5] );
            end
            
            %scenario 7
            if(E_to_store(period_1)<=0 && E_to_store(period_2)>=0)
                scenario=7;
                dx1 = (ELO-E_to_store(period_1));
                dx2 = (E_to_store(period_2)-ELI)/time_loss;
                dx5 = -(constraint(period_1)-Output_to_grid(period_1))/eta_out;
                %new constraint
                dx_test = (constraint_2(period_2) - Output_to_grid(period_2))*(eta_in/time_loss);
                dx = max ( [dx dx1 dx2 dx5 dx_test] );
            end
            
            %scenario 8
            if(E_to_store(period_1)<=0 && E_to_store(period_2)<0)
                scenario=8;
                
                dx1 = (ELO-E_to_store(period_1));
                dx2 = E_to_store(period_2)/time_loss;
                dx5 = -(constraint(period_1)-Output_to_grid(period_1))/eta_out;
                
                dx = max ( [dx dx1 dx2 dx5] );
            end
            
            i_temp = find(E_stored==min(E_stored(period_1:period_2-1)), 1, 'first');
                    
            time_loss_2 = exp((time(period_1)-time(i_temp))/tau);
            m = min(E_stored(period_1:(period_2 - 1)));
            dx3 = (-m/time_loss_2)+(1e-12);
                                
            if(dx3>0); dx3 = 0; end
            
            dx = max( [dx dx3] );

            %The price checkers - accept if revenue increases
            %Works out the values Output_to_grid will take if move is
            %accepted
            
            new_EtS_1 = E_to_store(period_1) + dx; 
            new_EtS_2 = E_to_store(period_2) - dx*time_loss;
            
            if(scenario==5)
                new_OtG_1=Output_to_grid(period_1)-dx/eta_in;
                new_OtG_2=Output_to_grid(period_2)+(dx*time_loss)/eta_in;
            end            
            if(scenario==6)
                new_OtG_1=Output_to_grid(period_1)-dx/eta_in;
                new_OtG_2=Output_to_grid(period_2)+(dx*time_loss*eta_out);
            end
            if(scenario==7)
                new_OtG_1=Output_to_grid(period_1)-dx*eta_out;
                new_OtG_2=Output_to_grid(period_2)+(dx*time_loss)/eta_in;
            end
            if(scenario==8)
                new_OtG_1=Output_to_grid(period_1)-dx*eta_out;
                new_OtG_2=Output_to_grid(period_2)+(dx*time_loss*eta_out);
            end
            Rev1=Output_to_grid(period_1)*test_price(period_1) + Output_to_grid(period_2)*test_price(period_2);
            Rev2=new_OtG_1*test_price(period_1) + new_OtG_2*test_price(period_2);
            
            %Checks if the 'new' revenue is higher than the old and accepts
            %if it is. 
            
            accept = 0;
            if(new_OtG_1<=constraint(period_1) && new_OtG_2<=constraint(period_1))
            if(Rev2 > Rev1)
                accept=1;
            end
            end
                        
        end
        
    end
    
    %If acceptable on price, now check energy capacity between period_1 and
    %period_2
        
    ifill = period_1;
    while(accept == 1 && ifill<period_2)
        %check that the move doesn't exceed the max storage in the range period_1 to period_2.
        if(E_stored(ifill) + dx*(exp((time(period_1)-time(ifill))/tau))>Cap_max)
            accept = 0;
        end
        ifill = ifill + 1;
    end

    accept;
    
    ifill = period_1;
    while(accept == 1 && ifill<period_2)
        %check that the move doesn't exceed the max storage in the range period_1 to period_2.
        test = E_stored(ifill) + dx*(exp((time(period_1)-time(ifill))/tau));
        if(test<0)
            accept = 0;
        end
        ifill = ifill + 1;
    end
    
    accept;
      
    if dx == inf 
       error('Fatal error occurred \n'); 
    end
% trying to eliminate rounding errors    
    if(dx<1e-10 && dx>-1e-10)
        dx = 0;
    end
    
    new_ETS_1 = E_to_store(period_1) + dx;
    new_ETS_2 = E_to_store(period_2) - dx*time_loss;

% % Checking for errors, shouldn't flag up hopefully...    
%     if(free_energy_old(period_1)>0)
%         if(new_ETS_1<0)
%             new_ETS_1
%             accept=0
%         end
%     end
%     
%     if(free_energy_old(period_2)>0)
%         if(new_ETS_2<0)
%             new_ETS_2
%             accept=0
%         end
%     end
    
    
    if dx==0; accept = 0; end
    %now if the move is accepted then all the variables can be updated..
    if(accept == 1)
        accepted=accepted+1;
        time_loss = exp((time(period_1)-time(period_2))/tau);
       
        E_to_store(period_1) = E_to_store(period_1) + dx;
        E_to_store(period_2) = E_to_store(period_2) - dx*time_loss;
        
        if(E_to_store(period_1)<1e-10 && E_to_store(period_1)>-1e-10)
            E_to_store(period_1) = 0;
        end
        if(E_to_store(period_2)<1e-10 && E_to_store(period_2)>-1e-10)
            E_to_store(period_2) = 0;
        end
        
            if(scenario == 9)                
                free_energy(period_1) = free_energy(period_1)-dx/eta_in;
                f_used(period_1) = free_energy_old(period_1) - free_energy(period_1);
                Output_to_grid(period_2)= Output_to_grid(period_2)+dx*time_loss*eta_out;                              
            end
            if(scenario == 10)
                free_energy(period_1) = free_energy(period_1)-dx/eta_in;
                f_used(period_1) = free_energy_old(period_1) - free_energy(period_1);
                Output_to_grid(period_2)= Output_to_grid(period_2)+dx*time_loss/eta_in;                              
            end
            if(scenario == 11)
                free_energy(period_2) = free_energy(period_2)+(dx*time_loss)/eta_in;
                f_used(period_2) = free_energy_old(period_2) - free_energy(period_2);
                Output_to_grid(period_1)=Output_to_grid(period_1)-dx*eta_out;
            end
            if(scenario == 12)
                free_energy(period_2) = free_energy(period_2)+(dx*time_loss)/eta_in;
                f_used(period_2) = free_energy_old(period_2) - free_energy(period_2);
                Output_to_grid(period_1)=Output_to_grid(period_1)-dx/eta_in;
            end
            if(scenario==1)
                Output_to_grid(period_1)=Output_to_grid(period_1)-dx/eta_in;
                Output_to_grid(period_2)=Output_to_grid(period_2)+(dx*time_loss)/eta_in;
            end
            if(scenario==2)
                Output_to_grid(period_1)=Output_to_grid(period_1)-dx/eta_in;
                Output_to_grid(period_2)=Output_to_grid(period_2)+(dx*time_loss*eta_out);
            end
            if(scenario==3)
                Output_to_grid(period_1)=Output_to_grid(period_1)-dx*eta_out;
                Output_to_grid(period_2)=Output_to_grid(period_2)+(dx*time_loss)/eta_in;
            end            
            if(scenario==4)
                Output_to_grid(period_1)=Output_to_grid(period_1)-dx*eta_out;
                Output_to_grid(period_2)=Output_to_grid(period_2)+(dx*time_loss*eta_out);
            end            
            if(scenario==5)
                Output_to_grid(period_1)=Output_to_grid(period_1)-dx/eta_in;
                Output_to_grid(period_2)=Output_to_grid(period_2)+(dx*time_loss)/eta_in;
            end            
            if(scenario==6)
                Output_to_grid(period_1)=Output_to_grid(period_1)-dx/eta_in;
                Output_to_grid(period_2)=Output_to_grid(period_2)+(dx*time_loss*eta_out);
            end
            if(scenario==7)
                Output_to_grid(period_1)=Output_to_grid(period_1)-dx*eta_out;
                Output_to_grid(period_2)=Output_to_grid(period_2)+(dx*time_loss)/eta_in;
            end
            if(scenario==8)
                Output_to_grid(period_1)=Output_to_grid(period_1)-dx*eta_out;
                Output_to_grid(period_2)=Output_to_grid(period_2)+(dx*time_loss*eta_out);
            end

        
        %update the energy in the store
        
        if(dx>0)
            for ifill = period_1:(period_2 - 1);
                E_stored(ifill)=E_stored(ifill)+(dx*exp((time(period_1)-time(ifill))/tau));
            end
        end
        if(dx<0)
            for ifill = period_1:(period_2 - 1);
                E_stored(ifill)=E_stored(ifill)+dx*(exp((time(period_1)-time(ifill))/tau));
            end
        end
        
    end
   
    test_Rev(n) = sum(Output_to_grid * test_price');
    
% print the algorithm progress to the screen
    if(mod(n,100000)==0) 
        fprintf('%d \n',n);
    end


for x = 1:index_ref;
    if(E_to_store(x)<0 && Output_to_grid(x)+1e-10<new_test_energy(x))
        fprintf('Problem \n');
        fprintf('scenario = %d \n',scenario);
        fprintf('period_1 = %d \n',period_1);
        ETS1=E_to_store(period_1)-dx;
        fprintf('ETS1 = %f \n',ETS1);
        fprintf('period_2 = %d \n',period_2);
        ETS2=E_to_store(period_2) + dx*time_loss;
        fprintf('ETS2 = %f \n',ETS2);
        fprintf('dx = %f \n',dx);
        OtG2 = Output_to_grid(period_2)-(dx*time_loss*eta_out);
        fprintf('OtG2 = %f \n',OtG2);
        error('Model Run Error: Energy is output from store and output less than initial')
    end
end



end

Revenue = Output_to_grid.*test_price;

new_Rev = sum(Revenue);
Difference = new_Rev - old_Rev;

% % looking at the progress of the algorithm
% figure; plot(iterations, test_Rev, '.')
% xlabel('number of iterations')
% ylabel('Revenue')
% title('Revenue generated by storage schedule with algorithm iterations')

figure; 
plot(time, E_stored, time, 2*E_to_store, time, -2*Output_to_grid, time, price*1000)
xlabel('time (hrs)')
ylabel('energy-stored/energy-transfer (kWh/kW)')
legend('energy stored', 'energy transfer', 'energy input from grid', 'price (£/MWh)')
title('Illustrating the optimum schedule of storage and charging/discharging')

mycall1 = num2str(Difference);
message = sprintf(['Revenue = £' mycall1]);
uiwait(msgbox(message));

% following code used in RPG paper by Gill et al.

% O = 2*Output_to_grid;
% P = 2*new_test_energy;
% P_old = 2*old_test_energy;
% F = 2*free_energy_old;
% F2 = 2*free_energy;
% F_used = F-F2;

% figure;
% plot(time, P_old, 'g', time, P, 'b', time, F, 'r')
% axis([0 max(time) 0 max(P_old)+1]);
% legend('Standalone no Constraint', 'Standalone with Constraint', 'Spilled power')
% xlabel('time (hrs)')
% ylabel('Power (MW)')
% title('Figure Illustrating the action of a 100MW wind farm Feb 2009 subject to a constraint')

% figure;
% plot(time, O, 'b', time, P, 'r', time, price*1000, 'm');
% axis([0 200 Lower_limit max(O)+0.2]);
% legend('With storage', 'Standalone Constrained', 'price')
% xlabel('time (hrs)')
% ylabel('Power (MW)')
% title('Figure illustrating the action of storage')
% 
% figure;
% plot(time, F, 'm', time, F_used, 'c')
% xlabel('time (hrs)')
% ylabel('Power (MW)')
% title('curtailed energy used')