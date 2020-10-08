function TS_Spec = AbstractSpecificationDiscrete(Spec, Partition, TS, TransitionNumber)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT ARGUMENTS
% Spec is a structure with fields W and D (Np by Np by n_inputs cell arrays of functions) that defines a control task. Np is the
% number of discrete modes, n_inputs is the number of external inputs.
% Partition is a (n_x by 1) cell array of Nx(s) by 2 matrices that define a set of intervals along dimension s.
% TS is an (Nx_a by Nu_a by Nv_a) array. Nx_a is a total number of symbolic
% states, Nu_a is the maximal number of control input values available at a
% symbolic state, Nv_a is the maximal number of successors from a symbolic
% state for a given control value. TS(x, u, w) is the successor symbolic
% state from state x when symbolic control u and simbolic disturbance w is
% applied.
% OUTPUT ARGUMENTS
% TS_with_Spec is an (Nx_a*Np*n_inputs by Nu_a by Nv_a*n_inputs by Np) array. TS_with_Spec([x0; p0; v0], u, [w; v1], p1) = [x1; p1; v1] is the successor symbolic
% state from state [x0; p0; v0] when symbolic control u and simbolic
% disturbance [w; v1] are applied while the next mode is chosen to be p1. Every transition [x0; p0; v0] -> [x1; p1;
% v1] of this system has the corresponding transition [x0; p0] -> [x1; p1]
% for external input v0 in the specification system. Every transition [x0;
% p0; v0] -> [x1; p1; v1] for control input u has the corresponding
% transition x0 -> x1 for control input u in the original transition
% system.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n_inputs = size(Spec.W, 3); % number of external inputs in the specification
if isfield(Spec, 'isMonotone') % if the specification is monotone, the condition to check simplifies
    isMonotone = Spec.isMonotone;
else
    isMonotone = false;
end

% State space
n_x = length(Partition); % number of dimensions of the state space
Nx = zeros(n_x, 1); % number of elements along each dimension
for s = 1:n_x
    Nx(s) = size(Partition{s}, 1);
end
Np = size(Spec.W, 1);
Nx = [Nx; Np; n_inputs];
Nx_a = prod(Nx); % total number of elements in the abstract system

Nu_a = size(TS, 2); % maximal number of controls enabled at a symbolic state
Nv_a = size(TS, 3); % maximal number of successors for a symbolic state
Nt = [Nv_a; n_inputs];

TS_Spec = zeros(Nx_a, Nu_a, prod(Nt), Np);

xmin = zeros(n_x, 1); % lower interval bound for the current symbolic state
xmax = zeros(n_x, 1); % upper interval bound for the current symbolic state
xmin_plus = zeros(n_x, 1); % lower interval bound for the successor symbolic state
xmax_plus = zeros(n_x, 1); % upper interval bound for the successor symbolic state

fprintf('Computing specification abstraction...\n');
complete = 0.1;
for linearIndexState = 1:Nx_a
    i = Ind2Sub(Nx, linearIndexState);
    linearIndexX = Sub2Ind(Nx(1:n_x), i(1:n_x));
    for s = 1:n_x
        xmin(s) = Partition{s}(i(s), 1);
        xmax(s) = Partition{s}(i(s), 2);
    end
    p = i(n_x+1);
    v = i(n_x+2);
    
    for linearIndexU = 1:Nu_a
        Nv_a_max = TransitionNumber(linearIndexX, linearIndexU);
        for linearIndexV = 1:Nv_a_max
            i_plus = Ind2Sub(Nx(1:n_x), TS(linearIndexX, linearIndexU, linearIndexV));
            for s = 1:n_x
                xmin_plus(s) = Partition{s}(i_plus(s), 1);
                xmax_plus(s) = Partition{s}(i_plus(s), 2);
            end
            X = [xmin, xmax; xmin_plus, xmax_plus];
            
            empty = true;
            for p_plus = 1:Np
                TS_Spec(linearIndexState, linearIndexU, Sub2Ind(Nt, [linearIndexV; 1]), p_plus) = ...
                    Sub2Ind(Nx, [i_plus; p_plus; 1]);

                if isMonotone
                    if max([Spec.W{p, p_plus, v}(xmin, xmax_plus), Spec.D{p, p_plus, v}(xmin), ...
                            Spec.W{p, p_plus, v}(xmax, xmin_plus), Spec.D{p, p_plus, v}(xmax)]) > 0
                        TS_Spec(linearIndexState, linearIndexU, Sub2Ind(Nt, [linearIndexV; 1]), p_plus) = 0;
                    end
                else
                    for k = 1:2^(2*n_x)
                        vertex = Ind2Sub(2*ones(2*n_x, 1), k);
                        x = X(sub2ind([2*n_x, 2], (1:2*n_x)', vertex));
                        if max(Spec.W{p, p_plus, v}(x(1:n_x), x(n_x+1:end)), Spec.D{p, p_plus, v}(x(1:n_x))) > 0
                            TS_Spec(linearIndexState, linearIndexU, Sub2Ind(Nt, [linearIndexV; 1]), p_plus) = 0;
                            break;
                        end
                    end
                end
                
                if TS_Spec(linearIndexState, linearIndexU, Sub2Ind(Nt, [linearIndexV; 1]), p_plus) ~= 0
                    empty = false;
                    for v_plus = 2:n_inputs
                        TS_Spec(linearIndexState, linearIndexU, Sub2Ind(Nt, [linearIndexV; v_plus]), p_plus) = ...
                            Sub2Ind(Nx, [i_plus; p_plus; v_plus]);
                    end
                end
            end
            if empty
                break;
            end
        end
        
        if Nv_a_max > 0 && ~empty
            for linearIndexV = Nv_a_max+1:Nv_a
                for p_plus = 1:Np
                    for v_plus = 1:n_inputs
                        TS_Spec(linearIndexState, linearIndexU, Sub2Ind(Nt, [linearIndexV; v_plus]), p_plus) = ...
                            TS_Spec(linearIndexState, linearIndexU, Sub2Ind(Nt, [Nv_a_max; v_plus]), p_plus);
                    end
                end
            end
        end
    end
    
    if linearIndexState/Nx_a >= complete/100
        fprintf('%.1f %% complete.\n', complete);
        complete = complete + 0.1;
    end
end

fprintf('Abstraction computed.\n');

end

