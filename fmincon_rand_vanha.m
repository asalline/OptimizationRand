%This function is used to optimize the density matrix.
%The optimization method that is used is "fmincon" which requires
%Optimization Toolbox.

function [x, fval, history] = fmincon_rand(x0)
%Defining some used things.
    global pauli_new original_rho
    history = {};
%Below is the function with variables that are being optimized. In this
%case the vector "x" contains the expectation values of Pauli
%measurements, that are needed to obtain the density matrix.
    
%     testi 
    
    f = @(x) trace(sqrtm(1/4 * (pauli_new{1} + x(1)*pauli_new{2} + x(2)*pauli_new{3} + x(3)*pauli_new{4} ...
        + x(4)*pauli_new{5} + x(5)*pauli_new{6} + x(6)*pauli_new{7} + ...
        x(7)*pauli_new{8} + x(8)*pauli_new{9} + x(9)*pauli_new{10} + ...
        x(10)*pauli_new{11} + x(11)*pauli_new{12} + x(12)*pauli_new{13} + ...
        x(13)*pauli_new{14} + x(14)*pauli_new{15} + x(15)*pauli_new{16})' ...
        * (1/4 * (pauli_new{1} + x(1)*pauli_new{2} + x(2)*pauli_new{3} + x(3)*pauli_new{4} ...
        + x(4)*pauli_new{5} + x(5)*pauli_new{6} + x(6)*pauli_new{7} + ...
        x(7)*pauli_new{8} + x(8)*pauli_new{9} + x(9)*pauli_new{10} + ...
        x(10)*pauli_new{11} + x(11)*pauli_new{12} + x(12)*pauli_new{13} + ...
        x(13)*pauli_new{14} + x(14)*pauli_new{15} + x(15)*pauli_new{16}))));
 
%Next up there is the parameters that the "fmincon" could use.
    x0 = zeros(1,15);
    A = [];
    b = [];
    Aeq = [];
    beq = [];
    lb = [];
    ub = [];
%"nonlincon" is the main constraint here. It allows multiple non linear
%constraints to be used. "options" is used to obtain certain values
%without need of any loops.
    nonlincon = @nlcon_rand;
    options = optimset('OutputFcn', @myoutput, 'MaxFunEvals', 1000);
    
%This "fmincon" returns the values of "x".
    [x, fval, ~, output] = fmincon(f, x0, A, b, Aeq, beq, lb, ub,...
        nonlincon, options);
%     disp(x);
%     original_rho;
%     disp(output);

%This function controls and saves the values of each iteration step.
    function stop = myoutput(x, optimvalues, state);
        stop = false;
        if isequal(state,'iter')
            history = [history, 1/4 * (pauli_new{1} + x(1)*pauli_new{2} + x(2)*pauli_new{3} + x(3)*pauli_new{4} ...
        + x(4)*pauli_new{5} + x(5)*pauli_new{6} + x(6)*pauli_new{7} + ...
        x(7)*pauli_new{8} + x(8)*pauli_new{9} + x(9)*pauli_new{10} + ...
        x(10)*pauli_new{11} + x(11)*pauli_new{12} + x(12)*pauli_new{13} + ...
        x(13)*pauli_new{14} + x(14)*pauli_new{15} + x(15)*pauli_new{16})];
        end
    end
end