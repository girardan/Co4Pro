function ind = Sub2Ind(Nx, subs)
% A custom version of ind2sub function that accepts subscript indices in a
% single input argument. subs is a vector of the same length as Nx.

n_x = length(Nx);
ind = subs(end);
for i = n_x-1:-1:1
    ind = subs(i)+(ind-1)*Nx(i);
end
end