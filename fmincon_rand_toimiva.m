%This function is used to optimize the density matrix.
%The optimization method that is used is "fmincon" which requires
%Optimization Toolbox.

function [x, fval, history] = fmincon_rand(pauli_new, measurements, selection, qubits)
%Defining some used things.
    
    history = {};
%Below is the function with variables that are being optimized. In this
%case the vector "x" contains the expectation values of Pauli
%measurements, that are needed to obtain the density matrix.

%     optim = testifunktio(x, pauli_new)

%Next up there is the parameters that the "fmincon" could use.
    x0 = zeros(1, 4^qubits - 1);
    A = [];
    b = [];
    Aeq = [];
    beq = [];
    lb = [];
    ub = ones(size(x0));
%"nonlincon" is the main constraint here. It allows multiple non linear
%constraints to be used. "options" is used to obtain certain values
%without need of any loops.
%     nonlincon = @nlcon_rand;
    options = optimset('OutputFcn', @myoutput, 'MaxFunEvals', 10000);
    
%This "fmincon" returns the values of "x".
    [x, fval, ~, output] = fmincon(@(x) testifunktio(x, pauli_new, qubits), x0, A, b, Aeq, beq, lb, ub,...
        @(x) nlcon_rand(x, measurements, selection, qubits), options);
%     disp(x);
%     original_rho;
%     disp(output)

%This function controls and saves the values of each iteration step.
    function stop = myoutput(x, optimvalues, state);
        stop = false;
        if isequal(state,'iter')
            history = [history, optimized_rho_rand(x, pauli_new, qubits)];
        end
    end
end