function [V] = CFF_invpercentile(X,P)

X = X(:);
X = X(~isnan(X));
X = sort(X);
iP = round(P.*numel(X)./100);
V=X(iP);
