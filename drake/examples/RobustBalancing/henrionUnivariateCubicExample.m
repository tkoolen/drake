function henrionUnivariateCubicExample()
model = UnivariateCubic;

T = 100;
R_diag = 1 * ones(1, model.num_states);
goal_radius = 0.01;
target = @(x) goal_radius^2 - x'*x;
n = 0;
options.degree = 32;

nStepCapturabilitySOS(model, T, R_diag, target, n, options)
end
