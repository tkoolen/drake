function henrionUnivariateCubicExample()
model = UnivariateCubic;

T = 10;
R_diag = 1 * ones(1, model.num_states);
goal_radius = 0.3;
target = @(x) goal_radius^2 - x'*x;
options.degree = 10;

% nStepCapturabilitySOS(model, T, R_diag, target, n, options)

[V_inner, W_inner] = finiteTimeInnerApproximation(model,T,[],R_diag,target,options);
end
