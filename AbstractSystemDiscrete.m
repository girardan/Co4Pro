function [TS, TransitionNumber] = AbstractSystemDiscrete(F, Partition, Controls, Disturbance)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT ARGUMENTS
% F is an anonymous function of 4 vector arguments: states 'x' and 'y' (n_x by 1 each), control 'u' (n_u by 1) and
% disturbance 'w' (n_v by 1).
% Partition is a (n_x by 1) cell array of Nx(s) by 2 matrices that define a set of intervals along dimension s.
% Controls is a (n_u by 1) cell array of row vectors which define possible discretized control values for each component of control vector.
% Disturbance is a (n_v by 2) matrix which defines interval bounds on the disturbance.
% OUTPUT ARGUMENTS
% TS is an (Nx_a by Nu_a by Nv_a) array. Nx_a is a total number of symbolic
% states, Nu_a is the maximal number of control input values available at a
% symbolic state, Nv_a is the maximal number of successors from a symbolic
% state for a given control value. TS(x, u, w) is the successor symbolic
% state from state x when symbolic control u and simbolic disturbance w is
% applied.
% TransitionNumber is an (Nx_a by Nu_a) array. Its elements are the numbers of transitions
% for a particular state x and a control u.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% State space
n_x = length(Partition); % number of dimensions of the state space
Nx = zeros(n_x, 1); % number of elements along each dimension
for s = 1:n_x
    Nx(s) = size(Partition{s}, 1);
end
Nx_a = prod(Nx); % total number of elements in the abstract system

% Control space
n_u = length(Controls); % number of dimensions of the control space
Nu = zeros(n_u, 1); % number of elements along each dimension
for s = 1:n_u
    Nu(s) = length(Controls{s});
end
Nu_a = prod(Nu);

Nv_a = 16; % initial maximal number of successors for a symbolic state

TS = zeros(Nx_a, Nu_a, Nv_a);
TransitionNumber = zeros(Nx_a, Nu_a);

xmin = zeros(n_x, 1); % lower interval bound for the current symbolic state
xmax = zeros(n_x, 1); % upper interval bound for the current symbolic state
minsucc = zeros(n_x, 1); % index vector of the lower interval bound for the reach set over-approximation
maxsucc = zeros(n_x, 1); % index vector of the upper interval bound for the reach set over-approximation

fprintf('Computing system abstraction...\n');
complete = 1;
for linearIndexX = 1:Nx_a
    i = Ind2Sub(Nx, linearIndexX);
    for s = 1:n_x
        xmin(s) = Partition{s}(i(s), 1);
        xmax(s) = Partition{s}(i(s), 2);
    end
    for linearIndexU = 1:Nu_a
        j = Ind2Sub(Nu, linearIndexU);
        u = zeros(n_u, 1);
        for s = 1:n_u
            u(s) = Controls{s}(j(s));
        end
        
        Fmin = F(xmin, xmax, u, Disturbance(:, 1)); % lower interval bound for the reach set over-approximation
        Fmax = F(xmax, xmin, u, Disturbance(:, 2)); % upper interval bound for the reach set over-approximation
        
        isEnabled = true;
        for s = 1:n_x
            tmp = find(Partition{s}(:, 1) <= Fmin(s), 1, 'last');
            if isempty(tmp)
                isEnabled = false;
                break;
            end
            minsucc(s) = tmp;
            
            tmp = max(find(Partition{s}(:, 2) >= Fmax(s), 1, 'first'), tmp);
            if isempty(tmp)
                isEnabled = false;
                break;
            end
            maxsucc(s) = tmp;
        end
        
        if isEnabled
            Nv = maxsucc-minsucc+1;
            Nv_a_max = prod(Nv);
            TransitionNumber(linearIndexX, linearIndexU) = Nv_a_max;
            for linearIndexV = 1:Nv_a_max
                TS(linearIndexX, linearIndexU, linearIndexV) = Sub2Ind(Nx, minsucc+Ind2Sub(Nv, linearIndexV)-1);
            end
            TS(linearIndexX, linearIndexU, Nv_a_max+1:end) = Sub2Ind(Nx, maxsucc);
        else
            TransitionNumber(linearIndexX, linearIndexU) = 0;
            TS(linearIndexX, linearIndexU, :) = 0;
        end
    end
    
    if linearIndexX/Nx_a >= complete/100
        fprintf('%d %% complete.\n', complete);
        complete = complete + 1;
    end
end

% When the size of an array increases, new elements are initialized with
% zeros. We must change them to the appropriate values.
for linearIndexX = 1:Nx_a
    for linearIndexU = 1:Nu_a
        if TransitionNumber(linearIndexX, linearIndexU) > 0
            TS(linearIndexX, linearIndexU, TransitionNumber(linearIndexX, linearIndexU)+1:end) = ...
                TS(linearIndexX, linearIndexU, TransitionNumber(linearIndexX, linearIndexU));
        end
    end
end

fprintf('Abstraction computed.\n');

end