%This function contains all the non linear constraints that are used to
%optimize with "fmincon_rand".

function [c, ceq] = nlcon_rand(x, measurements, selection, qubits)
% global selection
error = 0.0001*rand;
c16 = 0;, c15 = c16;, c14 = c15;, c13 = c14;, c12 = c13;, c11 = c12;, ...
    c10 = c11;, c9 = c10;, c8 = c9;, c7 = c8;, c6 = c7;, c5 = c6;, ...
    c4 = c5;, c3 = c4;, c2 = c3;, c1 = c2;
c = [c1;c2;c3;c4;c5;c6;c7;c8;c9;c10;c11;c12;c13;c14;c15];

%These are correlating measurements to Pauli basis by ceq n = Pauli_new{n+1}
% ceq1 = abs(x(1) - measurements(2)); %Pauli_new{2}
% ceq2 = abs(x(2) - measurements(3));
% ceq3 = abs(x(3) - measurements(4));
% ceq4 = abs(x(4) - measurements(5));
% ceq5 = abs(x(5) - measurements(6));
% ceq6 = abs(x(6) - measurements(7));
% ceq7 = abs(x(7) - measurements(8));
% ceq8 = abs(x(8) - measurements(9));
% ceq9 = abs(x(9) - measurements(10));
% ceq10 = abs(x(10) - measurements(11));
% ceq11 = abs(x(11) - measurements(12));
% ceq12 = abs(x(12) - measurements(13));
% ceq13 = abs(x(13) - measurements(14));
% ceq14 = abs(x(14) - measurements(15));
% ceq15 = abs(x(15) - measurements(16)); %Pauli_new{16}

for j = 1:(4^qubits - 1)
    ceq_all(j) = norm(x(j) - measurements(j+1));
end

% ceq1 = norm(x(1) - measurements(2));
% ceq2 = norm(x(2) - measurements(3));
% ceq3 = norm(x(3) - measurements(4));
% ceq4 = norm(x(4) - measurements(5));
% ceq5 = norm(x(5) - measurements(6));
% ceq6 = norm(x(5) - measurements(7));
% ceq7 = norm(x(7) - measurements(8));
% ceq8 = norm(x(8) - measurements(9));
% ceq9 = norm(x(9) - measurements(10));
% ceq10 = norm(x(10) - measurements(11));
% ceq11 = norm(x(11) - measurements(12));
% ceq12 = norm(x(12) - measurements(13));
% ceq13 = norm(x(13) - measurements(14));
% ceq14 = norm(x(14) - measurements(15));
% ceq15 = norm(x(15) - measurements(16));

%valinnat = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15];
% ceq_all = [ceq1;ceq2;ceq3;ceq4;ceq5;ceq6;ceq7;ceq8;ceq9;ceq10;ceq11;ceq12;ceq13;ceq14;ceq15];
ceq = ceq_all(selection);
end