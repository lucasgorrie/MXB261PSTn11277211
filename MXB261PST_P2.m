%% n11277211, P2 %%

% Poisson and KL anonymous functions %
f = @(lambda, k) ((lambda.^k).*(exp(-lambda)))./factorial(k);
KL = @(P, Q) sum(P.*(log(max((P./(Q + 10^-30)), 10^-10))));  % Addition of 10^-30 necessary to avoid divison by zero. Only way I can think of solving this without making this a non-anonymous function

% Desired Inputs %
Lambda = 4;
k = 0:15;

% Real Distribution %
PMF = f(Lambda, k);

% Scalar Uniform Normalization %
difference = 1 - sum(PMF);
differencePerElement = difference/length(PMF);
PMF = PMF + differencePerElement;  % Element-wise addition

CMF = cumsum(PMF);  % CMF

% Initialization %
sampleSizes = [10, 25, 50, 100, 175, 250];
N = 100;  % How many times to repeat sampling

% 1x6 Cell array of NxSampleNum arrays %
Samples = cell(1, length(sampleSizes));
EmpiricalPMFS = cell(1, length(sampleSizes));
KLDivergences = cell(1, length(sampleSizes));

MeanKLDivergences = zeros(1, length(sampleSizes));

% Fill cells with empty arrays of variable dimensions %
for i = 1:length(sampleSizes)

    % Cell block of sample arrays %
    SampleArray = zeros(N, sampleSizes(i));
    Samples{i} = SampleArray;

    % Cell block of PMFs associated with sample array %
    PMFArray = zeros(N, length(k));  
    EmpiricalPMFS{i} = PMFArray;

    % 1KL per PMF row %
    KLArray = zeros(N, 1);
    KLDivergences{i} = KLArray;

end

% Sample & Interpret %
for i = 1:length(sampleSizes)

    for j = 1:N

        % Inverse Transform Sampling %
        for q = 1:sampleSizes(i)
            r = rand();
            Samples{i}(j, q) = find((r <= CMF), 1);
        end

        % Construct PMF %
        [M, K] = groupcounts(Samples{i}(j, :)', "IncludeEmptyGroups", true);  % M = number of occurences, K = integers which occur

        % PMF
        C = zeros(1, length(k));
        C(K) = C(K) + M';
        EmpiricalPMFS{i}(j, :) = C/sampleSizes(i);  % P = M/N

        % KL %
        KLDivergences{i}(j) = KL(PMF, EmpiricalPMFS{i}(j, :));

    end

    % Mean KL Divergence %
    MeanKLDivergences(i) = sum(KLDivergences{i}(:))./N;

end

% Sanity checks %
sum(EmpiricalPMFS{1}(1, :))   % 10 Sample size,  N = 1,  full PMF
sum(EmpiricalPMFS{1}(80, :))  % 10 Sample Size,  N = 80, full PMF
sum(EmpiricalPMFS{5}(7, :))   % 175 Sample Size, N = 7,  full PMF

% Plot 1 %
figure
hold on
plot(sampleSizes, MeanKLDivergences)
plot(sampleSizes, MeanKLDivergences, 'r.')
errorbar(sampleSizes, MeanKLDivergences, std(MeanKLDivergences)/sqrt(N))
title("KL Divergence of Empirical Poisson Against Empirical Sample Size")
xlabel("Sample size")
ylabel("KL Divergence")
xlim([0 275])

% Plot 2 %
figure
hold on
t = tiledlayout(2, 3);  % For 2x3 facet setup
title(t, "True Poisson Against Empirical Inverse Transform Sampling, Lambda = 4")

nexttile
bar(k, [EmpiricalPMFS{1}(50, :); PMF])
title(sprintf("Sample Size = %d", sampleSizes(1)))
legend("Empirical", "True")
xlabel("k")
ylabel("Probability")

nexttile
bar(k, [EmpiricalPMFS{2}(50, :); PMF])
title(sprintf("Sample Size = %d", sampleSizes(2)))
legend("Empirical", "True")
xlabel("k")
ylabel("Probability")

nexttile
bar(k, [EmpiricalPMFS{3}(50, :); PMF])
title(sprintf("Sample Size = %d", sampleSizes(3)))
legend("Empirical", "True")
xlabel("k")
ylabel("Probability")

nexttile
bar(k, [EmpiricalPMFS{4}(50, :); PMF])
title(sprintf("Sample Size = %d", sampleSizes(4)))
legend("Empirical", "True")
xlabel("k")
ylabel("Probability")

nexttile
bar(k, [EmpiricalPMFS{5}(50, :); PMF])
title(sprintf("Sample Size = %d", sampleSizes(5)))
legend("Empirical", "True")
xlabel("k")
ylabel("Probability")

nexttile
bar(k, [EmpiricalPMFS{6}(50, :); PMF])
title(sprintf("Sample Size = %d", sampleSizes(6)))
legend("Empirical", "True")
xlabel("k")
ylabel("Probability")

%% Animation for Personal Usage %%
L = 0:50;
k = 0:50;

y = zeros(length(L), length(k));
for i = 1:length(L)

    y(i, :) = f(L(i), k);
    
end
for i = 1:length(L)

    plot(k, y(i, :), 'b')
    xlim([0 k(end)])
    ylim([0 0.5])
    ylabel("Probability")
    xlabel("Empirical Number of Occurances K")
    xline(L(i), 'r')
    legend("", "Lambda")
    title(sprintf("Average Number of Occurances Lambda = %f", L(i)))
    drawnow
    pause(0.02)

end
