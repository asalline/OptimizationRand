function optimized_rho = optimized_rho(x, pauli_new, qubits)
optimized_rho = pauli_new{1};
for j = 1:numel(x)
    optimized_rho = optimized_rho + (x(j) * pauli_new{j+1});
end
optimized_rho = optimized_rho / 2^qubits;