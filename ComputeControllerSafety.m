function [V, C, domC] = ComputeControllerSafety(TS, X1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT ARGUMENTS
% TS is an (Nx by Nu by Nv by Np) array. Nx is a total number of symbolic
% states, Nu is the maximal number of control input values available at a
% symbolic state, Nv is the maximal number of successors from a symbolic
% state for a given control value, Np is the number of modes in the
% specification. TS(x, u, w, p) is the successor symbolic 
% state from state x when symbolic control u and simbolic disturbance w are
% applied while the next mode is chosen to be p.
% X1 is a column vector of symbolic terminal states.
% OUTPUT ARGUMENTS
% C is a (Nx by Nu) boolean matrix that represents the safety
% controller. C(x, u) is true if control input u is enabled at a state x of
% the closed-loop system.
% V is a (Nx by 1) vector that represents the value function. V(x) == 0 if the
% all the trajectories of the closed-loop system are infinite, V(x) == 1 if
% all the trajectories of the closed-loop system are either infinite or
% reach the terminal set X1, V(x) = + Inf otherwise.
% domC is a (Nx by 1) boolean vector of controllable states. Terminal
% states from X1 are considered controllable.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    Nx = size(TS, 1);
    Nu = size(TS, 2);
    Nv = size(TS, 3);
    Np = size(TS, 4);
    fixedPointFound = false;
    
    states = (1:Nx)';
    Xu = Nx+1; % Blocking state
    TS(TS == 0) = Xu;
    
    V = zeros(Nx+1, 1);
    V(Xu) = Inf;
    X = (1:Nx)';
    
    succ_number = size(TS(1, :, :, :));
    
    while ~fixedPointFound
        Vtmp = V;
        for x_ind = 1:length(X)
            V(X(x_ind)) = min(max(min(reshape(Vtmp(TS(X(x_ind), :, :, :)), succ_number), [], 4), [], 3), [], 2);
        end
        V(X1) = min(V(X1), 1);
        if Vtmp == V
            fixedPointFound = true;
        end
    end
    
    C = zeros(Nx, Nu, Nv, Np, 'logical');
    C(X1, :, :, :) = false;
    Xsafe = (V < Inf);
    C(Xsafe, :, :, :) = Vtmp(TS(Xsafe, :, :, :)) <= repmat(V(Xsafe), 1, Nu, Nv, Np);
    
    domC = states(V(1:Nx) < Inf);
end


