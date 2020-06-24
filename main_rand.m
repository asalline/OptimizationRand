%This function triggers the solver and then calculates and plots the
%fidelity (closeness of two quantum states) with every iteration step.
%The solver itself optimizes the density matrix with restrictions that are
%randomly obtained.

%Needed scripts to operate: RDM_parempi.m & fmincon_rand.m & nlcon_rand.m &
%                           paulimatrices.m

%There is five variables:
%   repeats: amount of the times this script generates new random density
%            matrix and tries to optimize it (natural number, default: 1)
%   qubits: amount of qubits to use (natural number, default: 1)
%   ranknum: rank of the density matrix (natural number, default: 1)
%   real: defines if the density matrix is real or complex valued
%         (real valued = 1, complex valued = 0, default: 1)
%   amount_of_randoms: amount of randomly selected measurements done in
%                      related Pauli basis. This variable is used in
%                      constraint function "nlcon_rand".
clear all
clc
rounds = 3;
mean_iters = zeros(1,301);

for times = 1:rounds

    amount_of_randoms = 7;
    possibilities = [1:1:15];
    %This loop generates the random set of measurements.
    for k = 1:amount_of_randoms
        selection_txt(k) = possibilities(randi([1,length(possibilities)]));
        possibilities(possibilities == selection_txt(k)) = [];
    end
    selection_txt = sort(selection_txt, 'ascend')
    selection = selection_txt;
    global selection

    clear means fidelity
    clearvars -except collection times means_length keskiarvot selection ...
        last_fidelity mean_iters mean_fidelities amount_of_randoms
    repeats = 15;
    qubits = 2;
    ranknum = 1;
    rank = ranknum;
    real = 0;
    fidelities = cell(1, repeats);
    x0 = zeros(1,15);

    %This loop repeats process for wanted times.
    for jj = 1:repeats
    RDM_parempi;
    paulimatrices;
    [x, fval, history] = fmincon_rand(x0);

    %Variable right below is to save iteration steps from every process.
    %Don't really know if this is important at all.
    mean_iters(length(history)-1) = mean_iters(length(history)-1) + 1;
    
    %For each time processed this loop calculates the fidelity between the
    %original density matrix generated and the optimized density matrix.
    if length(history) <= 65
        for k = 1:length(history)
            fidelity(k) =...
                (trace(sqrtm(sqrtm(history{k})*original_rho*sqrtm(history{k}))))^2;
        end
    else
        fidelity = [];
    end
    fidelity;
    steps = [1:1:length(fidelity)];
    fidelities{jj} = fidelity;

    %Figure that plots fidelities on y-axis with respect to iteration steps on
    %x-axis.
    figure(1);
    grid on, hold on
    plot(steps, 1-fidelity, 'k.', 'markersize', 5);
    xlabel('Iteration steps');
    ylabel('Infidelity (1-fidelity)');

%     rho = history{end};
%     disp('Optimized density matrix');
%     disp(rho);
%     disp('Original density matrix');
%     disp(original_rho);

    clear fidelity
    end
    
    mean_fidels = 0;
    for l = 1:repeats
        mean_fidels = mean_fidels +(1- fidelities{1,l}(end));
    end
    mean_fidelities{times} = mean_fidels / repeats;

    %Next part calculates needed size for the vector that is needed to
    %calculate means between each full optimization.
    for kk = 1:repeats
        step = length(fidelities{kk});
        all_steps(kk) = step;
    end
    max_length = max(all_steps);
    means = zeros(1, max_length);
    amount_of_steps_per_iteration = means;

    for l = 1:repeats
        fidelity =...
            [fidelities{l}, zeros(1, max_length - length(fidelities{l}))];
        means = means + fidelity;
    end

    scaled_means = nan(max_length, length(fidelities));
    for j = 1:size(scaled_means,2)
        scaled_means(1:length(cell2mat(fidelities(j))),j) =...
            cell2mat(fidelities(j));
    end

    scaled_means = 1 - nanmean(scaled_means,2);
    means = means/repeats;
    steps = [1:1:length(means)];

    plot(steps, scaled_means, 'g.', 'markersize', 6);
    if length(steps) < 65
        xlim([steps(1), steps(end)]);
    else
        xlim([steps(1), steps(35)]);
    end
    xlabel('Iteration steps');
    ylabel('Mean fidelity');
%     green = plot(steps, scaled_means, 'g--');
    green(times) = plot(steps, scaled_means, 'g--')


%     legend([green], 'Scaled mean')
    means_length(times) = length(scaled_means);
    collection{times} = scaled_means;

end

% mean_fidels = 0;
% for l = 1:repeats
%     mean_fidels = mean_fidels + fidelities{1,l}(end);
% end
% mean_fidels = mean_fidels / repeats;

final_means = x0; final_vars = x0;
final_means(amount_of_randoms) = mean(cell2mat(mean_fidelities));
final_vars(amount_of_randoms) = var(cell2mat(mean_fidelities));

best_step = find(mean_iters == max(mean_iters))

selection_txt = num2str(selection);
amount_of_measurements = length(selection);
formatSpec1 = "Iterated: %d times, Amount of qubits: %d, Rank of the density matrices: %d";
    str1 = sprintf(formatSpec1, repeats, qubits, rank);
    title({['Data from density matrix optimization']; [str1]; ...
        ['Measurements: ', selection_txt]})

len = min(means_length);
new_steps = (1:1:len);

scaled_means2 = nan(len, length(collection));
for j = 1:size(scaled_means2,2)
    scaled_means2(1:length(cell2mat(collection(j))),j) =...
        cell2mat(collection(j));
end
scaled_means2 = nanmean(scaled_means2,2);
last_fidelity = zeros(1,15);
last_fidelity(amount_of_measurements) = scaled_means2(end);
steps2 = [1:1:length(scaled_means2)];
figure(10)
hold on
grid on
% for k = 1:3
%     collection{k} = collection{k}(1:len);
%     plots{k} = plot(new_steps, collection{k}, '--');
% end
plot(steps2, scaled_means2, '.')
hold off


clear selection