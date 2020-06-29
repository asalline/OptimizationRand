%Random density matrix generator for qubits

%Author: Antti Sällinen
%Last update: 10.6.2020

%First this script generates random state vectors correlating to amount of
%wanted qubits. 
%If rank of the wanted density matrix is 1 < rank <= 2^qubits, this script 
%forms density matrix with that rank. Difference in that case is that then there
%is no state vectors for different qubits, since it forms the density matrix 
%all by using randomness and the properties of linear algebra.
%For rank = 1 situations one can see different state vectors for all the
%qubits from the cell "state{}"

%There is three variables:
%   qubits: amount of qubits to use (natural number, default: 1)
%   ranknum: rank of the density matrix (natural number, default: 1)
%   real: defines if the density matrix is real or complex valued
%         (real valued = 1, complex valued = 0, default: 1)

%Three variables should not be "activated" while running with other
%programs since it can cause problems.

%R = 1;, r = R;
%C = 0;, c = C;
% qubits = 3;
% ranknum = 1;
% real = 0;

% original_rho = randomdensitymatrix(qubits, ranknum, real);

function original_rho = RDM_parempi(qubits, ranknum, real)

%Defining base vectors |0> and |1>.
ket0 = [1, 0]';
ket1 = [0, 1]';

%Amount of qubits and the rank of density matrix is needed to proceed.
%Also the option to choose real or complex valued type of matrix.

clear rank prob trace
%qubits = input("Amount of qubits: ");
%ranknum = input("Choose rank: ");

%Defining holder cell to the state vectors.
state = cell(1, qubits);

%Generating state vectors for all qubits.
%This part generates the vectors for pure states. Also this forms the
%density matrix for pure state.
if ranknum == 1
    if real == 1
        for j = 1:qubits
            alpha = -1 + 2*rand;
            beta = -1 + 2*rand;
            norm = sqrt(alpha^2 + beta^2);
            alpha = alpha/norm;, beta = beta/norm;
            state{j} = alpha * ket0 + beta * ket1;
        end
    else
        for j = 1:qubits
            alpha = (-1 + 2*rand) + 1i*(-1 + 2*rand);
            beta = (-1 + 2*rand) + 1i*(-1 + 2*rand);
            norm = sqrt(alpha^2 + beta^2);
            alpha = alpha/norm;, beta = beta/norm;
            state{j} = alpha * ket0 + beta * ket1;
        end
    end
    final_state = 1;
    for k = 1:qubits
        final_state = kron(final_state, state{k});
    end
    rho = final_state * final_state';
    Trace = trace(rho);
    rho = rho/Trace;
%     disp('Density matrix'), disp(rho)
else
%This part generates the density matrix for rank > 1.    
    if real == 1
        part = -1+2*rand(2^qubits, ranknum);
    else
        part = -1+2*rand(2^qubits, ranknum) + 1i*(-1+2*rand(2^qubits, ranknum));
    end
    [U,S,V] = svd(part);
    part = (U + eye(2^qubits))*part;
    rho = part*part';
    rho = rho / trace(rho);
    disp('Density matrix');, disp(rho);
end
% global original_rho
original_rho = rho;
clear rho
end