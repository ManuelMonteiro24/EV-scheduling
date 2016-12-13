%% Optimal Scheduling for Charging and Discharging of Electric Vehicles
%
% Authors: Manuel Monteiro, 79138
%          Leonor Fermoselle, 78493
%          Miguel Paulino, 79168
%          Inês Peixoto, 78130
%
% Date: 13/12/2016
%% Equal allocation scheme
% Solve the equal allocation scheme

function [f_cost, y, z, E, x] = equal_allocation(EV_info, L_b, N, M_V2G, M_CHG, k0, k1, tau, bat_cap, pmax, fe_ratio, bat_life, b, niu)

M = M_V2G + M_CHG;          % Total number of EVs
E_fin = fe_ratio*bat_cap;   % Final energy of any EV

% Variables initialization
x = zeros(M,N);
z = zeros(1,N);
y = zeros(1,N);
E = zeros(M,N);
f = zeros(M,N);

% Calculation of charging power matrix
for m = 1:M
    
    %defining matrix f (charging-interval matrix)
    f(m, EV_info(m,1):EV_info(m,2)) = 1;
    if m <= M_CHG
        % Vehicle only charges. Calculates the charging power of EV m for
        % each interval that belongs to the charging period.
        x(m,EV_info(m,1):EV_info(m,2)) = ((E_fin-EV_info(m,3))/EV_info(m,4)).*ones(1,EV_info(m,4));
    else
        % Vehicle charges and discharges. Calculates the charging power of
        % EV m for each interval that belongs to the charging period. The
        % interval that the vehicle discharges is the interval
        % corresponding to the maximum of the electricity price.
        x(m,EV_info(m,1):EV_info(m,2)) = ((E_fin-EV_info(m,3))/(EV_info(m,4)-2)).*ones(1,EV_info(m,4));
        [~,price_max_ind] = max(k0 + k1.*L_b(EV_info(m,1):EV_info(m,2)));
        price_max_ind = price_max_ind + EV_info(m,1)-1;
        x(m,price_max_ind) = -x(m,price_max_ind);
    end
end

% Instant energy matrix
for i = 1:N
    E(:,i) = EV_info(:,3) + sum(tau.*x(:,1:i).*f(:,1:i),2);
end

y = sum(x.*f,1)';    % Charging load
z = L_b + y;         % Total load

% Total cost
f_cost = sum(k0.*z + (k1/2).*(z.^2) - k0.*L_b - (k1/2).*(L_b.^2));
if bat_life     % Battery lifetime reduction
    f_cost = f_cost + sum(sum(b.*(x.^2), 2)) + sum(sum(niu.*((x(:,2:N)-x(:,1:(N-1))).^2),2));
end

end