%% Optimal Scheduling for Charging and Discharging of Electric Vehicles
%
% Authors: Manuel Monteiro, 79138
%          Leonor Fermoselle, 78493
%          Miguel Paulino, 79168
%          Inês Peixoto, 78130
%
% Date: 13/12/2016
%% Global Scheduling solution
% Solve the optimization problem

function [f_cost, y, z, E, x] = global_allocation(EV_info, L_b, N, M_V2G, M_CHG, k0, k1, tau, bat_cap, pmax, fe_ratio, bat_life, b, niu, pc, cc, cf)

M = M_V2G + M_CHG;          % Total number of EVs
E_fin = fe_ratio*bat_cap;   % Final energy of any EV

%defining matrix f (charging-interval matrix)
f = zeros(M,N);
for m = 1:M
    f(m, EV_info(m,1):EV_info(m,2)) = 1;
end

cvx_begin quiet

variables x(M,N) z(N) y(N) E(M,N);

%Cost function
f_cost_vec = k0.*z + (k1/2).*(z.^2) - k0.*L_b - (k1/2).*(L_b.^2);

if pc    % Second part - piecewise constant
    for i = 2:N
        f_cost_vec(i) = f_cost_vec(i) + cc*norm([z(i)-z(i-1)],1);
    end
    f_cost = sum(f_cost_vec);
elseif cf       % Second part - Cost linear by branches
    f_cost = sum(max(k0.*z,(k0+k1).*z-(2.*10^-4)) - k0.*L_b - (k1/2).*(L_b.^2));
else
    if bat_life     % Battery lifetime reduction
        f_cost = sum(f_cost_vec) + sum(sum(b.*(x.^2), 2)) + sum(sum(niu.*((x(:,2:N)-x(:,1:(N-1))).^2),2));
    else            % Global optimization
        f_cost = sum(f_cost_vec);
    end
end
minimize(f_cost);

%subject to   

%Constrain that secures that the total load in an interval
%is equal to the base load in that interval plus the sum of
%the loads of all the cars plugged to a station in that interval
y == sum(x.*f,1)';
z == L_b + y;

%Constrain that secures that the energy of a vehicle at the end of an
%interval is no less than zero or greater than the battery capacity
for i = 1:N
	E(:,i) == EV_info(:,3) + sum(tau.*x(:,1:i).*f(:,1:i),2);
	E(:,i) <= bat_cap;
	E(:,i) >= 0;
end

%Constrain that secures that the final energy of a vehicle is no less
%than the speficied energy level
(EV_info(:,3) + sum(tau.*x.*f, 2)) >= E_fin.*ones(M,1);

%Constrain that secures that the charging power of a charging only vehicle
%isn't lower than zero or higher that the maximum charging power
if M_CHG > 0
	x(1:M_CHG,:) <= (pmax.*ones(M_CHG,N));
	x(1:M_CHG,:) >= zeros(M_CHG,N);
end;

%Constrain that secures that the charging power of a vehicle to grid type
%of EV isn't lower than the minimun charging power or higher than the
%maximum charging power
if M_V2G > 0
	x(M_CHG+1:M_CHG+M_V2G,:) <= (pmax.*ones(M_V2G,N));
	x(M_CHG+1:M_CHG+M_V2G,:) >= (-pmax.*ones(M_V2G,N));
end;

cvx_end;
end 