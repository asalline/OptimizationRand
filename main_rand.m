% RANDOM MEASUREMENT DIRECTIONS TO OPTIMIZE RANDOM DENSITY MATRIX

% Author: Antti Sällinen
% Last update: 2.7.2020

% This function triggers the solver and then calculates and plots the
% fidelity (closeness of two quantum states) with every iteration step.
% The solver itself optimizes the density matrix with restrictions that are
% randomly obtained.

% Needed scripts to operate: RDM_parempi.m, fmincon_rand.m, nlcon_rand.m,
%                           paulimatrices.m, optimized_rho_rand.m &
%                           to_optimize.m
% Needed toolboxes to operate: Optimization Toolbox

% There is six variables:
%   wanted_measurements: vector that contains the numbers of measurements 
%                        thet one wants to use for measurement amounts.
%   repeats: amount of the times this script generates new random density
%            matrix and tries to optimize it (natural number, default: 1)
%   qubits: amount of qubits to use (natural number, default: 1)
%   ranknum: rank of the density matrix (natural number, default: 1)
%   real: defines if the density matrix is real or complex valued
%         (real valued = 1, complex valued = 0, default: 1)
%   rounds: number of sets to be generated in different amount of 
%           random measurements

clear all
clc
tic
repeats = 1;
qubits = 5;
ranknum = 1;
real = 0;
final_means = zeros(1, 4^qubits - 1); final_vars = final_means;
possible_measurements = [2:1:4^qubits-1];

% wanted_measurements = [5:5:4^qubits-1];
wanted_measurements = [400];

measurement_ratio = wanted_measurements(end) / 4^qubits;

for choice = wanted_measurements
    clearvars -except collection means_length best_fidelity final_means ...
        final_vars choice repeats qubits ranknum real mean_infidelities ...
        possible_measurements
    amount_of_randoms = choice;
    if amount_of_randoms == possible_measurements(end)
        rounds = 1;
    elseif amount_of_randoms == possible_measurements(end-1)
        rounds = 5;
    else
        rounds = 1;
    end
    for times = 1:rounds
        clear measurement selection selection_txt
        possibilities = [1:1: 4^qubits - 1];
        
        % This loop generates the random set of measurements.
        for k = 1:amount_of_randoms
            selection_txt(k) = possibilities(randi([1,length(possibilities)]));
            possibilities(possibilities == selection_txt(k)) = [];
        end
        selection_txt = sort(selection_txt, 'ascend');
        selection = selection_txt;
        
        rank = ranknum;
        fidelities = cell(1, repeats);
%         x0 = zeros(1, 4^qubits - 1);

        % This loop repeats process for wanted times.
        for jj = 1:repeats
            original_rho = RDM_parempi(qubits, ranknum, real);
            original_rho;
            [pauli_new, measurements] = paulimatrices(original_rho, qubits);
            measurements;
            [x, fval, history] = fmincon_rand(pauli_new, ...
                measurements, selection, qubits);
            
            % For each time processed this loop calculates the fidelity 
            % between the original density matrix generated and the 
            % optimized density matrix.
            % NOTICE: due the rare unability to optimize of fmincon, those
            % values are not taken into data. If one wants to take those,
            % one can delete the next if-loop and only calculate the
%             % fidelities with the for-loop.
            if length(history) <= 100

                for k = 1:length(history)
                    fidelity(k) =...
                        (trace(sqrtm(sqrtm(history{k})*original_rho* ...
                        sqrtm(history{k}))))^2;
                end
            else
                fidelity = []
            end
            fidelity;
            steps = [1:1:length(fidelity)];
            fidelities{jj} = fidelity;

            % Figure that plots fidelities on y-axis with respect to
            % iteration steps on x-axis.
            figure(1);
            grid on, hold on
            plot(steps, 1-fidelity, 'k.', 'markersize', 5);
            xlabel('Iteration steps');
            ylabel('Infidelity (1-fidelity)');
            
            % Under this comment is some things that can be displayed
            % but it is not recommended since this function is not made to
            % do it.
            
        %     rho = history{end};
        %     disp('Optimized density matrix');
        %     disp(rho);
        %     disp('Original density matrix');
        %     disp(original_rho);

            clear fidelity
        end

        % This loop calculates and saves the mean infidelity from all the
        % last fidelities from every repeat within this particular set.
        mean_fidels = 0;
        for l = 1:repeats
            mean_fidels = mean_fidels +(1- fidelities{l}(end));
        end
        mean_infidelities{times} = mean_fidels / repeats;

        % Next part calculates needed size for the vector that is needed to
        % calculate means between each full optimization.
        for kk = 1:repeats
            step = length(fidelities{kk});
            all_steps(kk) = step;
        end
        max_length = max(all_steps);
        means = zeros(1, max_length);

        % Edits the means -vector such that it adds every value from
        % fidelities -cell to the right place in means -vector.
        for l = 1:repeats
            fidelity =...
                [fidelities{l}, zeros(1, max_length - length(fidelities{l}))];
            means = means + fidelity;
        end

        % Creates nan-matrix and edits it the way that the scaled mean can
        % be calculated after the for-loop.
        scaled_means = nan(max_length, length(fidelities));
        for j = 1:size(scaled_means,2)
            scaled_means(1:length(cell2mat(fidelities(j))),j) =...
                cell2mat(fidelities(j));
        end
        scaled_means = 1 - nanmean(scaled_means,2);
        means = means/repeats;
        steps = [1:1:length(means)];

        % Plots the scaled mean values to the plot made above.
        plot(steps, scaled_means, 'g.', 'markersize', 6);
        xlabel('Iteration steps');
        ylabel('Mean fidelity');
        green(times) = plot(steps, scaled_means, 'g--');
        
        % Saves the values below, those are used later.
        means_length(times) = length(scaled_means);
        collection{times} = scaled_means;

        clear fidelity
    end

    % Means and variances from all the data of the different sets is saved
    % into vectors.
    final_means(amount_of_randoms) = mean(cell2mat(mean_infidelities));
    final_vars(amount_of_randoms) = var(cell2mat(mean_infidelities));

    % Makes title to the graph made above.
    selection_txt = num2str(selection);
    amount_of_measurements = length(selection);
    formatSpec1 = "Iterated: %d times, Amount of qubits: %d, Rank of the density matrices: %d";
        str1 = sprintf(formatSpec1, repeats, qubits, rank);
        title({['Data from density matrix optimization']; [str1]})

    % Forms the new scaled mean values the same way as in above, but this
    % time uses the values of collection -cell.
    len = max(means_length);
    scaled_means2 = nan(len, length(collection));
    for j = 1:size(scaled_means2,2)
        scaled_means2(1:length(cell2mat(collection(j))),j) =...
            cell2mat(collection(j));
    end
    scaled_means2 = nanmean(scaled_means2,2);

    % Creates vector that saves the "best" (last) fidelity of each set.
    best_fidelity = zeros(1, 4^qubits - 1);
    best_fidelity(amount_of_measurements) = scaled_means2(end);
    
    % Plots the second scaled mean w.r.t iteration steps.
    steps2 = [1:1:length(scaled_means2)];
    figure(10)
    hold on
    grid on
    plot(steps2, scaled_means2, '.')
    hold off
end
%%
% Plots the complete mean and variance of each iteration of the each set of
% every amount of wanted measurements into nice graph.
meases = (1:1: 4^qubits - 1);
minuserror = final_means - final_vars;
pluserror = final_means + final_vars;
meansandvars = figure('name', 'Means w.r.t amount of measurements');
hold on, grid on
plot(meases, minuserror, 'r.', 'markersize', 3)
plot(meases, pluserror, 'r.', 'markersize', 3)
for k = 1:length(meases)
    plot([meases(k), meases(k)], [pluserror(k), minuserror(k)], 'r');
end
plot(meases, final_means, 'k.', 'markersize', 5)
ylim([0, 0.8]);
xlabel('Amount of random measurement directions');
ylabel('Mean infidelity and variance');
toc
% saveas(meansandvars, 'testi.png');

% NEXT PART IS CURRENTLY WORKING ONLY WITH FINAL VALUES OF EVERYTHING!!!
clear real
measurement_matrix = 0;
for k = selection
    measurement_matrix = measurement_matrix + pauli_new{k};
end
real_measurements = real(measurement_matrix);
imag_measurements = imag(measurement_matrix);
figure('name', 'Real measurements')
colormap gray
image(real_measurements, 'CDataMapping', 'scaled')
colorbar;
figure('name', 'Imaginary measurements')
colormap gray
image(imag_measurements, 'CDataMapping', 'scaled')
colorbar;
figure('name', 'Real part closeness');
% map = [0 0 0.3; 0 0 0.4; 0 0 0.5; 0 0 0.6; 0 0 0.7; 0 0 0.8; 0 0 0.9; 0 0 1.0];
colormap gray
Re_closeness = abs(real(original_rho) - real(history{end}));
image(Re_closeness, 'CDataMapping', 'scaled')
colorbar
figure('name', 'Imaginary part closeness');
colormap gray
Im_closeness = abs(imag(original_rho) - imag(history{end}));
image(Im_closeness, 'CDataMapping', 'scaled')
colorbar;

%That's all Folks!