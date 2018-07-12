clear;clc
load '../data/SF/SFwork.mat';
load '../data/SF/SFshop.mat';

data = SFshop; % SFshop, SFwork
max_iter = 100;
lambda = 1;
epsilon = 0.001;

[estimates_asr, convergence_iter_asr, likelihood_per_iter_asr] = asr_general(data, lambda, epsilon, max_iter);
[estimates_lsr, convergence_iter_lsr, likelihood_per_iter_lsr] = lsr_general(data, lambda, epsilon, max_iter);

iter_to_converge = max([convergence_iter_lsr, convergence_iter_asr]);
h = figure;
plot(1: iter_to_converge, likelihood_per_iter_asr(1: iter_to_converge),'LineWidth', 2);
hold on;
plot(1: iter_to_converge, likelihood_per_iter_lsr(1: iter_to_converge),':','LineWidth', 2);
hold on;
xlim([1 iter_to_converge]);

legend({'ASR','LSR'}, 'FontSize', 16);
xlabel('Iteration number','Interpreter','latex', 'FontSize', 24, 'FontWeight', 'bold');
ylabel('Log-likelihood','Interpreter','latex','FontSize', 24, 'FontWeight', 'bold');

hold off;