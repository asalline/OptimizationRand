function [pauli_new, measurements] = paulimatrices(original_rho, qubits)

%Defining identity matrix "I" and Pauli matrices pauli_x, pauli_y, pauli_z.

clear measurements

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
original_rho;
for kasvu = 1:qubits
    if kasvu == 1
        pauli_new = pauli_base;
    else
        for j = 1:length(pauli_new)
            for l = 1:length(pauli_base)
                pauli_temp{n} = kron(pauli_new{j}, pauli_base{l});
                n = n + 1;
            end
        end
        clear pauli_new
        pauli_new = pauli_temp;
        clear pauli_temp
        n = 1;
    end
end
for n = 1:length(pauli_new)
    measurements(n) = trace(pauli_new{n}*original_rho);
end

% concatenate = cat(3, pauli_new{:});
% sum_of_pauli = sum(concatenate, 3);

end