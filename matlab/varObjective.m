function f = varObjective(params,model,Xs)

[void f] = gpPosteriorMeanVar(model,[Xs params]);
% For speed skip scaling
%f = 1./sqrt(2*pi*f(:,1));
f = f(:,1);

return