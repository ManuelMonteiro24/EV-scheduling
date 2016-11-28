%% Optimal Scheduling for Charging and Discharging of Electric Vehicles
%
% Authors: Manuel Monteiro, Leonor Fermoselle, Miguel Paulino, In?s
% Peixoto
%
% Date: xx/12/2016
%%
%
tic,

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

[f_cost, y, z, E, x] = global_solution(EV_info, L_b, N, M_V2G, M_CHG, k0, k1, tau, bat_cap, pmax, fe_ratio);

% Comparison of the real base load and the forecasted base load
figure(1);
 plot(L_b, 'b');
 axis([0 24 800 1600]);
 xlabel('Time [h]');
 ylabel('Load [kW]');
 legend('Real base load', 'Location', 'southeast');
 
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
 
 toc,