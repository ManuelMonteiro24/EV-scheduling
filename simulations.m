%% Optimal Scheduling for Charging and Discharging of Electric Vehicles
%
% Authors: Manuel Monteiro, Leonor Fermoselle, Miguel Paulino, Inês
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
 hold off;
 axis([0 24 800 1600]);
 ylabel('Total load [kW]');
 xlabel('Time [h]');
 legend('Globally optimal scheme','Base load', 'Location', 'southeast');
 
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
 
% Variation of total cost with different charging-only ratio simulation
charging_only_ratio = zeros(6,1);
total_cost = zeros(6,1);
charging_only_ratio(1) = M_CHG/M;
total_cost(1) = f_cost;
for i = 2:6
    charging_only_ratio(i) = charging_only_ratio(i-1) + 0.2;
    M_CHG = M*charging_only_ratio(i);
    M_V2G = M - M_CHG;
    [f_cost, ~, ~, ~, ~] = global_solution(EV_info, L_b, N, M_V2G, M_CHG, k0, k1, tau, bat_cap, pmax, fe_ratio);
    total_cost(i) = f_cost;
end

figure(6);
  plot(charging_only_ratio, total_cost, '-o');
  axis([0 1 230 270]);
  ylabel('Total cost [C$]');
  xlabel('Charging-only ratio');
  legend('Globally optimal scheme');
 
 % Performance evaluation under different group sizes (Total cost and Total
 % load)
 % Não faz sentido fazer estes gráficos sem as optimizacoes locais
%  figure(7);
%   axis([1 200 230 270]);
%   xlabel('Average group size [EVs]');
%   ylabel('Total cost [C$]');
%   
%  figure(8);
%   axis([0 24 800 1600]);
%   xlabel('Time [h]');
%   ylabel('Total load [kW]');

% Comparison of total load when considering the cost of battery lifetime
% reduction
M_V2G = 200;
M_CHG = 0;

% beta model parameter
b = 5*(10^-4); %unities: C$/kWh^2

%niu model parameter
niu = (10^-3); %unities: C$/kWh^2

[~, ~, z, ~, ~] = global_solution_batlifetime(EV_info, L_b, N, M_V2G, M_CHG, k0, k1, tau, bat_cap, pmax, fe_ratio, b, niu);

figure(9);
 plot(z, 'b');      %Globally optimal scheme
 hold on;
 plot(L_b, 'r');    %Base load
 hold off;
 axis([0 24 800 1600]);
 ylabel('Total load [kW]');
 xlabel('Time [h]');
 legend('Globally optimal scheme','Base load', 'Location', 'southeast');
 
% Variation of total cost with different load forecasting error
% ?????? figure(10)
 
% Performance evaluation under different number of EVs (Total cost and
% Total load) - So temos dados para 200 EVs, nao para 100, 300 e 400
% figure(12);
%  plot(L_b, 'b');
%  for M_CHG = [100 200 300 400]
%      [f_cost, ~, z, ~, ~] = global_solution(EV_info, L_b, N, M_V2G, M_CHG, k0, k1, tau, bat_cap, pmax, fe_ratio);
%      figure(11);
%      hold on;
%      bar(f_cost, 'b');
%      hold off;
%      figure(12);
%      hold on;
%      plot(z);
%      hold off;
%  end
%  figure(11);
%  axis([0 500 0 600]);
%  xlabel('Number of EVs');
%  ylabel('Total cost [C$]');
%  legend('Globally optimal scheme');
%  figure(12);
%  axis([0 24 800 1600]);
%  xlabel('Time [h]');
%  ylabel('Total load [kW]');
%  legend('Base load', '100 EVs', '200 EVs', '300 EVs', '400 EVs');
 
 toc,