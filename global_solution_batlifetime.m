%information variables

%Interval set (24 intervalos)
N = 24;

%EV's set (200 carros)
M = 200;
% others sets ???? 
M_V2G = 0;
M_CHG = 200;

%Electricity price model parameters (obtidos na seccao simulation settings)
k0 = 10^-4; % unidades: C$/kWh
k1 = 1.2*(10^-4); %unities: C$/kWh/kW

%Length of an interval (obtido no sistem model da solucao global)
tau = 1; %unities: 1 hora

%Battery capacity (obtida na seccao simulation settings)
bat_cap = 16; %unities: 16kWh 

%Maximum charging power (obtido na seccao simulation settings)
pmax = 5; %unities: 5kW

%Final energy ratio required ???? confirmar????? pelo menos para as
%experiencias eles consideram isso, (obtido na seccao simulation settings)
fe_ratio = 0.9; 

%matriz f
f = zeros(M,N);

for m = 1:M
    for i = (EV_info(m,1)+1):(EV_info(m,2))
       f(m,i) = 1;
    end;
end;

% beta model parameter
b = 5*(10^-4); %unities: C$/kWh^2

%niu model parameter
niu = (10^-3); %unities: C$/kWh^2

%Solve the optimization problem

cvx_begin quiet

variables x(M,N) z(N);

%Charging load at interval i
for i = 1:N
    var = 0;
    for m = 1:M
        var = (((x(m,i))*(f(m,i))) + var);
    end;
    z(i) == L_b(i) + var; %var = y(i)
end;

%Cost function
f_cost_1 = (k0*z(1)+(k1/2)*(z(1)^2)) - (k0*L_b(1)+(k1/2)*(L_b(1)^2));
for i = 2:N
    f_cost_1 = f_cost_1 + (k0*z(i)+(k1/2)*(z(i)^2)) - (k0*L_b(i)+(k1/2)*(L_b(i)^2));
end;

f_cost_2 = 0;
for m = 1:M
    for i = 1:N
        f_cost_2 = b*((x(m,i))^2) + f_cost_2;
    end;
end;

f_cost_3 = 0;
for m = 1:M
    for i = 2:N
        f_cost_3 = (niu*(x(m,i) - x(m,(i-1)))^2) + f_cost_3;
    end;
end;

minimize(f_cost_1 + f_cost_2 + f_cost_3);

%subject to   

%Constrain that secures that the total load in an interval
%is equal to the base load in that interval plus the sum of
%the loads of all the cars plugged to a station in that interval

for i = 1:N
    var = 0;
    for m = 1:M
        var = (x(m,i)*f(m,i)) + var;
    end;
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
%isnt lower than zero or higher that the maximum charging power
if M_CHG > 0 
    for m = 1:M_CHG
        for i = 1:N 
            x(m,i) <= pmax;
            x(m,i) >= 0;
        end;
    end;
end;
%Constrain that secures that the charging power of a vehicle to grig type of
%EV isnt lower than the minimun charging power or higher than the maximum 
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

%simulation

%sum of all x(m,i) in one interval
for i = 1:N
    var = 0;
    for m = 1:M
    var = x(m,i) + var;    
    end;
    x_i(i) = var;
end;

%plot solution
figure(1); clf; 
 subplot(1,2,1); stem(x_i,'LineWidth',5);
 title('Xi');
 subplot(1,2,2); stem(z,'LineWidth',5);
 title('zi');
