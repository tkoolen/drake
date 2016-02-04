clear all

degree = 6;

R_diag = [1.6 1.6];

prog = spotsosprog;
[prog,t]=prog.newIndeterminate('t',1); % time
[prog,x]=prog.newIndeterminate('x',2); % state
u = msspoly('u',1);

A = diag(1./(R_diag.^2));
h_X = 1 - x'*A*x;

V_vars = [t;x];
W_vars = x;

f = [x(2);u];
T = 1;

[prog,V] = prog.newFreePoly(monomials(V_vars,0:degree));
[prog,W] = prog.newFreePoly(monomials(W_vars,0:degree));

[prog,p] = prog.newFreePoly(monomials(V_vars,0:degree),1);
Vdot = diff(V,x)*f + diff(V,t);
dVdotdu = diff(Vdot,u);
Vdot = subs(Vdot,u,u*0) + sum(p);

prog = prog.withPos(subs(subs(V,t,T),x,[0;0]));

[prog, Vdot_sos] = spotless_add_sprocedure(prog, -Vdot, h_X,V_vars,degree-2);
[prog,Vdot_ind] = prog.withSOS(Vdot_sos);


p_pos_sos = msspoly;
p_neg_sos = msspoly;
p_pos_ind = [];
p_neg_ind = [];
for i=1:1,
  p_pos_sos_i = p(i) - dVdotdu(i);
  [prog, p_pos_sos_i] = spotless_add_sprocedure(prog, p_pos_sos_i, h_X,V_vars,degree-2); 
  [prog, p_pos_sos_i] = spotless_add_sprocedure(prog, p_pos_sos_i,  t * (T - t),V_vars,degree-2);
  [prog,p_pos_ind(i)] = prog.withSOS(p_pos_sos_i);
  
  p_neg_sos_i = p(i) + dVdotdu(i);
  [prog, p_neg_sos_i] = spotless_add_sprocedure(prog, p_neg_sos_i, h_X,V_vars,degree-2);
  [prog, p_neg_sos_i] = spotless_add_sprocedure(prog, p_neg_sos_i,  t * (T - t),V_vars,degree-2);
  [prog,p_neg_ind(i)] = prog.withSOS(p_neg_sos_i);
  
  p_pos_sos = [p_pos_sos;p_pos_sos_i];
  p_neg_sos = [p_neg_sos;p_neg_sos_i];
end

% (3) W(x) >= 0 for x in X
[prog, W_sos] = spotless_add_sprocedure(prog, W, h_X,W_vars,degree-2);
prog = prog.withSOS(W_sos);

% (4) W(x) >= V(0,x) + 1 for x in X
[prog, WminusV_sos] = spotless_add_sprocedure(prog, W - subs(V,t,0) - 1, h_X,W_vars,degree-2);
prog = prog.withSOS(WminusV_sos);

%% Set up cost function -- integration over a sphere
cost = spotlessIntegral(prog,W,x,R_diag,[],[]);

%% Solve
spot_options = spotprog.defaultOptions;
spot_options.verbose = true;
spot_options.do_fr = true;
solver = @spot_mosek;
% solver = @spot_sedumi;
sol = prog.minimize(cost,solver,spot_options);


%% Plotting
Vsol = sol.eval(V);
Wsol = sol.eval(W);
figure(1)
contourSpotless([Wsol;h_X],x(1),x(2),[-R_diag(1) R_diag(1)],[-R_diag(2) R_diag(2)],[],[],[1 0],{'b','r'});
axis([-.7 .7 -1.2 1.2])


%% controller extraction (Majumdar et al. Mark's old code)
  % get sos decomps  
  y_vdot = sol.prog.sosEqsDualVars{Vdot_ind};
  basis_vdot = sol.prog.sosEqsBasis{Vdot_ind};
  mu_vdot = double(sol.dualEval(y_vdot));
  [M_vdot,G] = momentMatrix(mu_vdot,basis_vdot);
  
  M_sig = cell(1,1);
  u_sol = msspoly;
  for i=1:1,
    y_p_pos = sol.prog.sosEqsDualVars{p_pos_ind(i)};
    basis_p_pos = sol.prog.sosEqsBasis{p_pos_ind(i)};
    y_p_neg = sol.prog.sosEqsDualVars{p_neg_ind(i)};
    basis_p_neg = sol.prog.sosEqsBasis{p_neg_ind(i)};
    
    if length(basis_vdot) ~= length(basis_p_pos)
      err1 = 1;
    elseif length(basis_p_neg) ~= length(basis_p_pos)
      err2 = 1;
    else
      [~,err1] = double(basis_vdot-basis_p_pos);
      [~,err2] = double(basis_p_neg-basis_p_pos);
    end
    if err1 || err2
      warning(['Same basis did not come out of SOS ' ...
        'Decomposition.']);
      keyboard
    end
    
    y_sig = double(sol.dualEval(y_p_pos-y_p_neg));
    M_sig{i} = momentMatrix(y_sig,basis_vdot);
    
    u_sol = [u_sol;solveController(M_vdot,M_sig{i},G)];
  end
  
  
  %% Test sim
  x0 = [-.44;.78];
  x0 = [0;.45];
  h = 1e-3;
  t_vec = 0:h:T;
  N = length(t_vec);
  x_vec = zeros(2,N);
  x_vec(:,1) = x0;
  u_vec = zeros(1,N);
  for i=2:N,
    u_vec(i) = max(min(double(subs(u_sol,[t;x],[t_vec(i);x_vec(:,i-1)])),1),-1);
    x_vec(:,i) = x_vec(:,i-1) + h*[x_vec(2,i-1);u_vec(i)];
  end
  figure(2)
  plot(t_vec,x_vec,t_vec,u_vec)