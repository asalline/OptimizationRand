function [pauli_new, measurements] = paulimatrices(original_rho, qubits)

clear measurements

% Defining identity matrix "I" and Pauli matrices pauli_x, pauli_y, pauli_z.
% Also defining pauli_base, which contains the original Pauli matrices, and
% pauli_new which will contain Pauli basis for wanted dimension.
I = [1, 0; 0, 1];
pauli_x = [0, 1; 1, 0];
pauli_y = [0, -1i; 1i, 0];
pauli_z = [1, 0; 0, -1];

pauli_base = {I, pauli_x, pauli_y, pauli_z};
pauli_new = cell(1, length(pauli_base)*length(pauli_base));

n = 1;
% For-loop that forms the new Pauli basis.
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

% Loop that calculates the expectation values for each Pauli basis -matrix.
% These are named measurements because of other purposes.
for n = 1:length(pauli_new)
    measurements(n) = trace(pauli_new{n}*original_rho);
end
end