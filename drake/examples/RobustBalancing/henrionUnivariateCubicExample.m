function henrionUnivariateCubicExample()
model = UnivariateCubic;

T = 100;
R_diag = 1 * ones(1, model.num_states);
goal_radius = 0.01;
target = @(x) goal_radius^2 - x'*x;
n = 0;
options.plotfun = @plotfun;
options.degree = 32;

nStepCapturabilitySOS(model, T, R_diag, target, n, options)
end

function plotfun(n, Vsol, Wsol, h_X, R_diag, t, x)
% figure 1 in http://arxiv.org/pdf/1208.1751.pdf
xs = linspace(-1, 1, 100);
Ws = zeros(length(xs));
for i = 1 : length(xs)
  Ws(i) = subs(Wsol, [x; t], [xs(i); 0]);
end
plot(xs, Ws);
end
