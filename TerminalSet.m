function X1 = TerminalSet(Partition, Spec, Condition)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT ARGUMENTS
% Partition is a (n_x by 1) cell array of Nx(s) by 2 matrices that define a set of intervals along dimension s.
% Condition is an anonymous function of (n_x by 1) vector argument which
% returns real values.
% OUTPUT ARGUMENTS
% X1 is a vector of symbolic terminal states. It is assumed that a symbolic
% state is in X1 if all vertices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n_inputs = size(Spec.W, 3); % number of external inputs in the specification

n_x = length(Partition); % number of dimensions of the state space
Nx = zeros(n_x, 1); % number of elements along each dimension
for s = 1:n_x
    Nx(s) = size(Partition{s}, 1);
end
Np = size(Spec.W, 1);
Nx = [Nx; Np; n_inputs];
Nx_a = prod(Nx); % total number of elements in the abstract system

X1 = [];
xmin = zeros(n_x, 1);
xmax = zeros(n_x, 1);
fprintf('Computing terminal set abstraction...\n');
complete = 1;
for linearIndexState = 1:Nx_a
    i = Ind2Sub(Nx, linearIndexState);
    for s = 1:n_x
        xmin(s) = Partition{s}(i(s), 1);
        xmax(s) = Partition{s}(i(s), 2);
    end
    p = i(n_x+1);
    v = i(n_x+2);
    X = [xmin, xmax];
    
    isTerminal = true;
    for k = 1:2^n_x
        vertex = Ind2Sub(2*ones(n_x, 1), k);
        x = X(sub2ind([n_x, 2], (1:n_x)', vertex));
        if Condition(x, p, v) > 0
            isTerminal = false;
        end
    end
    if isTerminal
        X1 = [X1; linearIndexState]; %#ok<AGROW>
    end
    
    if linearIndexState/Nx_a >= complete/100
        fprintf('%d %% complete.\n', complete);
        complete = complete + 1;
    end
end

fprintf('Abstraction computed.\n');
end

