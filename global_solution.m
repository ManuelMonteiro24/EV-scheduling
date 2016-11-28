%% Part I - Global Scheduling solution
%Solve the optimization problem

function [f_cost, y, z, E, x, f] = global_solution(EV_info, L_b, N, M_V2G, M_CHG, k0, k1, tau, bat_cap, pmax, fe_ratio)

M = M_V2G + M_CHG;

%defining matrix f (charging-interval matrix)
f = zeros(M,N);
for m = 1:M
    f(m, EV_info(m,1):EV_info(m,2)) = 1;
end

cvx_begin quiet

variables x(M,N) z(N) y(N) E(M,N);

%Cost function
f_cost = sum(k0.*z + (k1/2).*(z.^2) - k0.*L_b - (k1/2).*(L_b.^2));

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
(EV_info(:,3) + sum(tau.*x.*f, 2)) >= (fe_ratio*bat_cap).*ones(M,1);

%Constrain that secures that the charging power of a charging only vehicle
%isn't lower than zero or higher that the maximum charging power
if M_CHG > 0
	x(1:M_CHG,:) <= (pmax.*ones(M_CHG,N));
	x(1:M_CHG,:) >= zeros(M_CHG,N);
end;

%Constrain that secures that the charging power of a vehicle to grig type of
%EV isn't lower than the minimun charging power or higher than the maximum 
%charging power
if M_V2G > 0
	x(M_CHG+1:M_CHG+M_V2G,:) <= (pmax.*ones(M_V2G,N));
	x(M_CHG+1:M_CHG+M_V2G,:) >= (-pmax.*ones(M_V2G,N));
end;

cvx_end;
end 