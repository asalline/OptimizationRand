%Defining identity matrix "I" and Pauli matrices pauli_x, pauli_y, pauli_z.

global original_rho
% clear measurements

I = [1, 0; 0, 1];

pauli_x = [0, 1; 1, 0];
px = pauli_x;

pauli_y = [0, -1i; 1i, 0];
py = pauli_y;

pauli_z = [1, 0; 0, -1];
pz = pauli_z;

pauli_base = {I, px, py, pz};

finaloperator = 1;

pauli_new = cell(1, length(pauli_base)*length(pauli_base));

n = 1;

for j = 1:length(pauli_base)
    for k = 1:length(pauli_base)
        pauli_new{n} = kron(pauli_base{j}, pauli_base{k});
        measurements(n) = trace(pauli_new{n}*original_rho);
        n = n+1;
    end
end

% concatenate = cat(3, pauli_new{:});
% sum_of_pauli = sum(concatenate, 3);

global pauli_new measurements