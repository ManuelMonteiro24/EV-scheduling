%% Optimal Scheduling for Charging and Discharging of Electric Vehicles
%
% Authors: Manuel Monteiro, Leonor Fermoselle, Miguel Paulino, In?s
% Peixoto
%
% Date: xx/12/2016
%%
%Variable initialization

load('EV_data.mat');

N = 24;         %Interval set
M = 200;        %EV's set
M_V2G = 200;      %Vehicle-to-grid set
M_CHG = 0;    %Charging-only EV set

%Electricity price model parameters (obtidos na seccao simulation settings)
k0 = 10^-4; % units: [C$/kWh]
k1 = 1.2*(10^-4); %units: [C$/kWh/kW]

%Length of an interval (obtido no sistem model da solucao global)
tau = 1; %units: [hour]

%Battery capacity (obtida na seccao simulation settings )
bat_cap = 16; %units: [kWh] 

%Maximum charging power (obtido na seccao simulation settings)
pmax = 5; %units: [kW]

%Final energy ratio required ???? confirmar????? pelo menos para as
%experiencias eles consideram isso, (obtido na seccao simulation settings)
fe_ratio = 0.9; 

%defining matrix f (charging-interval matrix)
f = zeros(M,N);

for m = 1:M
    for i = (EV_info(m,1)):(EV_info(m,2))
       f(m,i) = 1 ;
    end;
end;   

%% Part I - Global Scheduling solution
%Solve the optimization problem

cvx_begin quiet

variables x(M,N) z(N) y(N) E(M,N);

%Cost function
f_cost = (k0*z(1)+(k1/2)*(z(1)^2)) - (k0*L_b(1)+(k1/2)*(L_b(1)^2));
for i = 2:N
    f_cost = f_cost + (k0*z(i)+(k1/2)*(z(i)^2)) - (k0*L_b(i)+(k1/2)*(L_b(i)^2));
end;

minimize(f_cost);

%subject to   

%Constrain that secures that the total load in an interval
%is equal to the base load in that interval plus the sum of
%the loads of all the cars plugged to a station in that interval
for i = 1:N
    var = 0;
    for m = 1:M
        var = (x(m,i)*f(m,i)) + var;
    end;
    y(i)== var;
    z(i) == L_b(i) + var;
end;

%Constrain that secures that the energy of a vehicle at the end of an
%interval is no less than zero or greater than the battery capacity
%confirmar esta constrain????
for m = 1:M
    for i = 1:N %confirmar a ordem destes 2 primeiros for. e indiferente???
        var = 0;
        for k = EV_info(m,1):i %confirmar isto??
            var = ((x(m,k)*f(m,k))*tau) + var;   
        end;
        E(m,i) == EV_info(m,3)+ var;
        EV_info(m,3)+ var <= bat_cap;
        EV_info(m,3)+ var >= 0;
    end;
end;

%Constrain that secures that the final energy of a vehicle is no less
%than the speficied energy level
for m = 1:M
    var = 0;
    for i = 1:N
        var = ((x(m,i)*f(m,i))*tau) + var;
    end;    
    (EV_info(m,3) + var) >= (fe_ratio*bat_cap);
end; 

%Constrain that secures that the charging power of a charging only vehicle
%isn't lower than zero or higher that the maximum charging power
if M_CHG > 0 
    for m = 1:M_CHG
        for i = 1:N 
            x(m,i) <= pmax;
            x(m,i) >= 0;
        end;
    end;
end;

%Constrain that secures that the charging power of a vehicle to grig type of
%EV isn't lower than the minimun charging power or higher than the maximum 
%charging power
if M_V2G > 0 
    for m = 1:M_V2G
        for i = 1:N  
            x(m,i) <= pmax;
            x(m,i) >= (-pmax);
        end;
    end;
end;

cvx_end;

%% Part I: Simulations

%sum of all x(m,i) in one interval
%for i = 1:N
 %   var = 0;
 %   for m = 1:M
 %   var = x(m,i) + var;    
 %   end;
 %   x_i(i) = var;
%end;

%plot solution
%figure(1); clf; 
% subplot(1,2,1); stem(x_i,'LineWidth',5);
% title('Xi');
% subplot(1,2,2); stem(z,'LineWidth',5);
% title('zi');

 %Charging load and Total load simulation
figure(2);
 plot(y, 'b');      %Globally optimal scheme
 axis([0 24 -200 400]);
 ylabel('Charging load [kW]');
 xlabel('Time [h]');
 legend('Globally optimal scheme');
 
figure(3);
 plot(z, 'b');      %Globally optimal scheme
 hold on;
 plot(L_b, 'r');    %Base load
 axis([0 24 800 1600]);
 ylabel('Total load [kW]');
 xlabel('Time [h]');
 legend('Globally optimal scheme','Base load');
 
 %Energy and Charging power of EV 19
figure(4);
 plot(E(19,:), 'b');    %Globally optimal scheme
 axis([0 24 0 15]);
 ylabel('Energy of EV 19 [kWh]');
 xlabel('Time [h]');
 legend('Globally optimal scheme');
 
figure(5);
 bar(x(19,:));  %Globally optimal scheme
 axis([0 25 -2 6]);
 ylabel('Charging power of EV 19 [kW]');
 xlabel('Time [h]');
 legend('Globally optimal scheme');
 
 %Charging only ratio simulation
%figure(6);
 %axis([0 1 230 270]);
 %ylabel('Total cost [C$]');
 %xlabel('Charging-only ratio');
 
 
 %% Part II - ???
 
 
 
 
 
 