%This function contains all the non linear constraints that are used to
%optimize with "fmincon_rand".

function [c, ceq] = nlcon_rand(x, measurements, selection, qubits)
% global selection
error = 0.0001*rand;
c16 = 0;, c15 = c16;, c14 = c15;, c13 = c14;, c12 = c13;, c11 = c12;, ...
    c10 = c11;, c9 = c10;, c8 = c9;, c7 = c8;, c6 = c7;, c5 = c6;, ...
    c4 = c5;, c3 = c4;, c2 = c3;, c1 = c2;
c = [c1;c2;c3;c4;c5;c6;c7;c8;c9;c10;c11;c12;c13;c14;c15];

for j = 1:(4^qubits - 1)
    ceq_all(j) = norm(x(j) - measurements(j+1));
end

ceq = ceq_all(selection);
end