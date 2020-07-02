% HUOMIO VAIHDA NIMI!!!!

function optimized = to_optimize(x, pauli_new, qubits)
x = [1, x];
optimized = 0;
for j = 1:numel(x)
    optimized = optimized + (x(j) * pauli_new{j});
end
optimized = optimized / 2^qubits;
optimized = trace(sqrtm(optimized' * optimized));
end