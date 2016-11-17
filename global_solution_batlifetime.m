function [f_cost, y, z, E, x] = global_solution_batlifetime(EV_info, L_b, N, M_V2G, M_CHG, k0, k1, tau, bat_cap, pmax, fe_ratio, b, niu)

M = M_V2G + M_CHG;

%matriz f
f = zeros(M,N);

for m = 1:M
    for i = (EV_info(m,1)):(EV_info(m,2))
       f(m,i) = 1;
    end;
end;

%Solve the optimization problem

cvx_begin quiet

variables x(M,N) z(N) y(N) E(M,N);

%Cost function
f_cost_1 = (k0*z(1)+(k1/2)*(z(1)^2)) - (k0*L_b(1)+(k1/2)*(L_b(1)^2));
for i = 2:N
    f_cost_1 = f_cost_1 + (k0*z(i)+(k1/2)*(z(i)^2)) - (k0*L_b(i)+(k1/2)*(L_b(i)^2));
end;

f_cost_2 = 0;
for m = 1:M
    f_cost_2 = f_cost_2 + b*((x(m,1))^2);
    for i = 2:N
        f_cost_2 = f_cost_2 + b*((x(m,i))^2) + (niu*(x(m,i) - x(m,(i-1)))^2);
    end;
end;

f_cost = f_cost_1 + f_cost_2;
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

end