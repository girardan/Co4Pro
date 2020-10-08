function subs = Ind2Sub(Nx, ind)
% A custom version of ind2sub function that calculates subscript indices in a
% single output argument. Accepts only scalar ind.

n_x = length(Nx);
subs = zeros(n_x, 1);
for i = 1:n_x
    subs(i) = mod(ind, Nx(i));
    ind = floor(ind / Nx(i)) + 1;
    if subs(i) == 0
        subs(i) = Nx(i);
        ind = ind - 1;
    end
end
end