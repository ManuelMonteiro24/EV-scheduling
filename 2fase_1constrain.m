%Constrain that secures that the total load provided by the
%grid in an interval is no more than 1300 kW
for i = 1:N
    z(i) <= 1300;
end;