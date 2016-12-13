%% Optimal Scheduling for Charging and Discharging of Electric Vehicles
%
% Authors: Manuel Monteiro, 79138
%          Leonor Fermoselle, 78493
%          Miguel Paulino, 79168
%          Inês Peixoto, 78130
%
% Date: 13/12/2016
%%
%

% Variable initialization
load('EV_data.mat');

N = 24;         %Interval set
M = 200;        %EV's set
M_V2G = 200;    %Vehicle-to-grid set
M_CHG = 0;      %Charging-only EV set

% Electricity price model parameters
k0 = 10^-4; % units: [C$/kWh]
k1 = 1.2*(10^-4); %units: [C$/kWh/kW]

% Length of an interval
tau = 1; %units: [hour]

% Battery capacity
bat_cap = 16; %units: [kWh] 

% Maximum charging power
pmax = 5; %units: [kW]

% Final energy ratio required
fe_ratio = 0.9;

% Battery life flag
bat_life = false;

% beta model parameter
b = 5*(10^-4); %units: C$/kWh^2

% niu model parameter
niu = (10^-3); %units: C$/kWh^2

% Piecewise constant flag
pc = false;

% Constant coeficient to define weight of parameter a in fcost.
cc = 0;

% Flag of the 1st challenge
cf = false;

% Comparison of the real base load and the forecasted base load
figure(1);
 plot(L_b, '-b', 'LineWidth', 2);   % Base load
 hold on;
 plot(F_L_b2, '-r', 'LineWidth', 2);    % Forecasted base load
 axis([0 24 800 1600]);
 xlabel('Time [h]');
 ylabel('Load [kW]');
 legend('Real base load', 'Forecasted base load 2', 'Location',...
        'southeast');

% Runs the global and equal allocation schemes
[f_cost, y, z, E, x] = global_allocation(EV_info, L_b, N, M_V2G, M_CHG, k0, k1, tau, bat_cap, pmax, fe_ratio, bat_life, b, niu, pc, cc, cf);
[f_cost_equal, y_equal, z_equal, E_equal, x_equal] = equal_allocation(EV_info, L_b, N, M_V2G, M_CHG, k0, k1, tau, bat_cap, pmax, fe_ratio, bat_life, b, niu);


% Charging load simulation
figure(2);
 plot(y, '-b', 'LineWidth', 2); %Globally optimal scheme without bat_life
                                %reduction
 hold on;
 plot(y_equal, '-g', 'LineWidth', 2);   %Equal allocation scheme
 hold off;
 axis([0 24 -200 400]);
 ylabel('Charging load [kW]');
 xlabel('Time [h]');

% Total load simulation
figure(3);
 plot(z, '-b', 'LineWidth', 2);  %Globally optimal scheme without bat_life
                                 %reduction
 hold on;
 plot(L_b, '-r', 'LineWidth', 2);    %Base load
 plot(z_equal, '-g', 'LineWidth', 2);   %Equal allocation scheme
 hold off;
 axis([0 24 600 1600]);
 ylabel('Total load [kW]');
 xlabel('Time [h]');
 
% Energy of EV 19 simulation
figure(4);
 plot(E(19,:), '-b', 'LineWidth', 2);  %Globally optimal scheme without
                                       %bat_life reduction
 hold on;
 plot(E_equal(19,:), '-g', 'LineWidth', 2);     %Equal allocation scheme
 hold off;
 axis([0 24 0 15]);
 ylabel('Energy of EV 19 [kWh]');
 xlabel('Time [h]');

% Charging power of EV 19 from globally optimal scheme without bat_life
% reduction and from equal allocation scheme
B = [x(19,:)' x_equal(19,:)' zeros(N,1)];


% Variation of total cost with different charging-only ratio simulation
charging_only_ratio = linspace(0,1,6); % Charging-only-ratio vector
% Initialization of variables
total_cost = zeros(6,1);
total_cost_equal = zeros(6,1);
total_cost(1) = f_cost;
total_cost_equal(1) = f_cost_equal;
% Cycle that runs the global and equal allocation schemes for each
% charging-only ratio
for i = 2:length(charging_only_ratio)
    M_CHG = M*charging_only_ratio(i); %Charging-only EV set    
    M_V2G = M - M_CHG;  %Vehicle-to-grid set
    [total_cost(i), ~, ~, ~, ~] = global_allocation(EV_info, L_b, N, M_V2G, M_CHG, k0, k1, tau, bat_cap, pmax, fe_ratio, bat_life, b, niu, pc, cc, cf);
    [total_cost_equal(i), ~, ~, ~, ~] = equal_allocation(EV_info, L_b, N, M_V2G, M_CHG, k0, k1, tau, bat_cap, pmax, fe_ratio, bat_life, b, niu);
end
figure(6);
 % Total cost of global allocation scheme
 plot(charging_only_ratio, total_cost, '-bo', 'LineWidth', 2);
 hold on;
 % Total cost of equal allocation scheme
 plot(charging_only_ratio, total_cost_equal, '-gv', 'LineWidth', 2);
 hold off;
 ylabel('Total cost [C$]');
 xlabel('Charging-only ratio');
 axis([0 1 225 275]);


M_V2G = 200;    %Vehicle-to-grid set
M_CHG = 0;      %Charging-only EV set


bat_life = true;    % Battery life flag

% Runs the global and equal allocation schemes
[f_cost, y, z, E, x] = global_allocation(EV_info, L_b, N, M_V2G, M_CHG, k0, k1, tau, bat_cap, pmax, fe_ratio, bat_life, b, niu, pc, cc, cf);
[f_cost_equal, y_equal, z_equal, E_equal, x_equal] = equal_allocation(EV_info, L_b, N, M_V2G, M_CHG, k0, k1, tau, bat_cap, pmax, fe_ratio, bat_life, b, niu);

% Charging load simulation
figure(2);
 hold on;
 plot(y, '--b', 'LineWidth', 2); %Globally optimal scheme (battery lifetime
                                 %reduction)
 hold off;
 legend('Globally optimal scheme', 'Equal allocation scheme',...
        'Globally optimal scheme with battery lifetime reduction',...
        'Location', 'north');

% Total load simulation
figure(3);
 hold on;
 plot(z, '--b', 'LineWidth', 2); %Globally optimal scheme (battery lifetime
                                 %reduction)
 hold off;
 legend('Globally allocation scheme', 'Base Load',...
        'Equal allocation scheme',...
        'Globally optimal scheme with battery lifetime reduction',...
        'Location', 'southeast');

% Energy of EV 19 simulation
figure(4);
 hold on;
 plot(E(19,:), '--b', 'LineWidth', 2);   %Globally optimal scheme (battery
                                         %lifetime reduction)
 hold off;
 legend('Globally optimal scheme', 'Equal allocation scheme',...
        'Globally optimal scheme with battery lifetime reduction',...
        'Location', 'northwest');

% Charging power of EV 19
B(:,3) = x(19,:)';
figure(5);
 bar(B); % Histogram with global and equal allocation schemes for EV 19
 axis([0 25 -2 6]);
 ylabel('Charging power of EV 19 [kW]');
 xlabel('Time [h]');
 legend('Globally optimal scheme', 'Equal allocation scheme',...
        'Globally optimal scheme with battery lifetime reduction',...
        'Location', 'northwest');

% Variation of total cost with different charging-only ratio simulation
total_cost(1) = f_cost;
total_cost_equal(1) = f_cost_equal;
for i = 2:length(charging_only_ratio)
    M_CHG = M*charging_only_ratio(i);
    M_V2G = M - M_CHG;
    [total_cost(i), ~, ~, ~, ~] = global_allocation(EV_info, L_b, N, M_V2G, M_CHG, k0, k1, tau, bat_cap, pmax, fe_ratio, bat_life, b, niu, pc, cc, cf);
    [total_cost_equal(i), ~, ~, ~, ~] = equal_allocation(EV_info, L_b, N, M_V2G, M_CHG, k0, k1, tau, bat_cap, pmax, fe_ratio, bat_life, b, niu);
end
figure(6);
 hold on;
 plot(charging_only_ratio, total_cost, '--bo', 'LineWidth', 2);
 plot(charging_only_ratio, total_cost_equal, '--gv', 'LineWidth', 2);
 hold off;
 legend('Globally optimal scheme', 'Equal allocation scheme',...
        'Globally optimal scheme with battery lifetime reduction',...
        'Equal allocation scheme with battery lifetime reduction',...
        'Location', 'southeast');

 bat_life = false; % Battery life flag

% Piecewise Constant Constraint

pc = true;      % Piecewise constant flag

cc = 0.5;       % Initial constant coeficient

% Runs the global allocation scheme
[~, ~, z, ~, ~] = global_allocation(EV_info, L_b, N, M_V2G, M_CHG, k0, k1, tau, bat_cap, pmax, fe_ratio, bat_life, b, niu, pc, cc, cf);

% Total Load
figure(7);
 plot(z, '-b', 'LineWidth', 2);      % Globally optimal scheme
 hold on;
 plot(L_b, '--r', 'LineWidth', 2);    % Base load
 hold off;
 axis([0 24 800 1600]);
 ylabel('Total load [kW]');
 xlabel('Time [h]');
 legend('Globally optimal scheme \lambda=0.5 ', 'Base load',...
        'Location', 'southeast');

% Total load using different constant coeficients
cc = 1;     % Constant coeficient

% Runs the global allocation scheme
[~, ~, z, ~, ~] = global_allocation(EV_info, L_b, N, M_V2G, M_CHG, k0, k1, tau, bat_cap, pmax, fe_ratio, bat_life, b, niu, pc, cc, cf);

figure(8);
 plot(z, '-b', 'LineWidth', 2);      % Globally optimal scheme
 hold on;

cc=10;     % Constant coeficient

% Runs the global allocation scheme
[~, ~, z, ~, ~] = global_allocation(EV_info, L_b, N, M_V2G, M_CHG, k0, k1, tau, bat_cap, pmax, fe_ratio, bat_life, b, niu, pc, cc, cf);

 plot(z, '--g', 'LineWidth', 2); 
 hold on;
 
 cc=1000;   % Constant coeficient
 
 [~, ~, z, ~, ~] = global_allocation(EV_info, L_b, N, M_V2G, M_CHG, k0, k1, tau, bat_cap, pmax, fe_ratio, bat_life, b, niu, pc, cc, cf);
 plot(z, '-r', 'LineWidth', 2); 
 hold on;
 
 plot(L_b, '--b', 'LineWidth', 2);    % Base load
 hold off;
 axis([0 24 800 1600]);
 ylabel('Total load [kW]');
 xlabel('Time [h]');
 legend('Globally optimal scheme \lambda=1 ',...
        'Globally optimal scheme \lambda=10 ',...
        'Globally optimal scheme \lambda=100 ','Base load',...
        'Location', 'southeast');