function henrionUnivariateCubicExample()
model = UnivariateCubic;

T = 100;
R_diag = 1 * ones(1, model.num_states);
goal_radius = 0.01;
target = @(x) goal_radius^2 - x'*x;
n = 0;
options.plotfun = @plotfun;

options.degree = 16;

nStepCapturabilitySOS(model, T, R_diag, target, n, options)
end

function plotfun(n, Vsol, Wsol, h_X, R_diag, t, x)
xs = linspace(-1, 1, 100);
Vs = zeros(length(xs));
for i = 1 : length(xs)
  Vs(i) = subs(Vsol, [x; t], [xs(i); 0]);
end
plot(xs, Vs);
end
