%% Toolbox test script #2
% This example is described in the following paper:
% V. Sinyakov and A. Girard. Formal Controller Synthesis from
% Specifications Given by Discrete-Time Hybrid Automata. 2019. HAL: https://hal.archives-ouvertes.fr/hal-02361404
%% System abstraction

M = 1370;
f_0 = 51;
f_1 = 1.2567;
f_2 = 0.4342;
r = 2;
T0 = 0.5;

F1 = @(u, v) (M*u - f_0 - f_1*v - f_2*v^2)/M ;
F2 = @(v_1, v_2) real(sqrt((v_1*T0)^2-r^2))/T0 - v_2;

n_x = 4;
F = @(x, y, u, w) [x(1) + (x(2)-y(3))*T0; min(max(x(2) + F1(u(1), y(2))*T0, 10), 30); min(max(x(3) + w*T0, 10), 25); x(4)] * ...
                  (u(2) == x(4)) + ...
                  [x(1) + F2(x(2), y(3))*T0; min(max(x(2) + F1(u(1), y(2))*T0, 10), 30); min(max(x(3) + w*T0, 10), 25); u(2)]*...
                  (u(2) ~= x(4));
        
Partition = cell(n_x, 1);

tmp = [-Inf, linspace(-70, 30, 51)]';
Partition{1} = [tmp(1:end-1)+eps(tmp(1:end-1)), tmp(2:end)];
Partition{1}(1, 1) = -Inf;

tmp = [10, linspace(10, 30, 41)]';
Partition{2} = [tmp(1:end-1)+eps(tmp(1:end-1)), tmp(2:end)];
Partition{2}(1, 1) = 10;

tmp = linspace(10, 25, 31)';
Partition{3} = [tmp(1:end-1)+eps(tmp(1:end-1)), tmp(2:end)];
Partition{3}(1, 1) = 10;

Partition{4} = [0.5, 0.5; 1.5, 1.5];

Controls = cell(2, 1);
Controls{1} = [0, 1, -1, 2, -2, -3, 3];
Controls{2} = [0.5, 1.5];
Disturbance = [-3, 3];

tic;
[TS, TransitionNumber] = AbstractSystemDiscrete(F, Partition, Controls, Disturbance);
toc;

%% Specification abstraction
Spec.isMonotone = true;
Spec.Inputs = 0;
N_inputs = length(Spec.Inputs);
Np = 3;
Spec.W = cell(Np, Np, N_inputs);
Spec.D = cell(Np, Np, N_inputs);
l = 10;
for i = 1:N_inputs
    for p1 = 1:Np
        for p2 = 1:Np
            Spec.W{p1, p2, i} = @(x, x_plus) -x_plus(1)+x(1);
            Spec.D{p1, p2, i} = @(x) 1;
            if p1 == 1 && p2 == 1
                Spec.D{p1, p2, i} = @(x) max(x(1) + l,x(4) - 1);
            elseif p1 == 1 && p2 == 2
                Spec.D{p1, p2, i} = @(x) max(x(1) + l,1-x(4));
            elseif p1 == 2 && p2 == 2
                Spec.D{p1, p2, i} = @(x) 1 - x(4);
            elseif p1 == 2 && p2 == 3
                Spec.D{p1, p2, i} = @(x) max(-x(1) + l, 1 - x(4));
            end
        end
    end
end

tic;
TS_with_Spec = AbstractSpecificationDiscrete(Spec, Partition, TS, TransitionNumber);
toc;

%% Controller computation
Condition = @(x, p, v) max([-x(1) + l, x(4) - 1, abs(p-3)]);
X1 = TerminalSet(Partition, Spec, Condition);

tic;
[V, C, domC] = ComputeControllerReachability(TS_with_Spec, X1);
toc;

%% Controllable sets visualization

[XX, YY, ZZ] = meshgrid(linspace(Partition{2}(2), Partition{2}(end)-1, 40)+1/2, ...
    linspace(Partition{1}(2), Partition{1}(end)-1, 50)+1/2, ...
    linspace(Partition{3}(1), Partition{3}(end)-1, 30)+1/2);

VV = reshape(V(1:(end-1)), [51, 41, 30, 2, 3]);


VV = VV(2:end, 2:end, :, 1, 1);
VV(VV == Inf) = 10000;
figure;
h = patch(isosurface(XX, YY, ZZ, VV, 2000));
set(h, 'EdgeColor', 'None');
set(h, 'FaceColor', [1, 0, 0]);
set(h, 'EdgeLighting', 'phong');
set(h, 'FaceLighting', 'phong');
set(h, 'FaceAlpha', 1);
set(h, 'DiffuseStrength', 0.6);
set(h, 'SpecularStrength', 0.9);
VV(1, :, :) = 2000+1e-4;
VV(:, 1, :) = 2000+1e-4;
VV(:, end, :) = 2000+1e-4;
VV(:, :, 1) = 2000+1e-4;
VV(:, :, end) = 2000+1e-4;
h1 = patch(isosurface(XX, YY, ZZ, VV, 2000));
set(h1, 'EdgeColor', 'None');
set(h1, 'FaceColor', [1, 0, 0]);
set(h1, 'EdgeLighting', 'phong');
set(h1, 'FaceLighting', 'phong');
set(h1, 'FaceAlpha', 0.1);
set(h1, 'DiffuseStrength', 0.6);
set(h1, 'SpecularStrength', 0.9);
xlabel('x^2','FontSize',14);
ylabel('x^1','FontSize',14);
zlabel('x^3','FontSize',14);
grid on;
legend('Controllable set','FontSize',14);
xlim([Partition{2}(1), Partition{2}(end)]);
ylim([Partition{1}(2), Partition{1}(end)]);
zlim([Partition{3}(1), Partition{3}(end)]);
view([0, 42]);
camlight('right');
view([53, 43]);



%% Trajectory simulation
rng(14042020);
x = zeros(4, 101);
x(:, 1) = [-40; 20; 10; 0.5];
z = zeros(5, 101);
z(:, 1) = [1; 1; 1; 1; 1];
u = zeros(2, 100);
for i = 1:100
    i
    for s = 1:n_x
        z(s, i) = find((Partition{s}(:, 1) <= x(s, i)) & (Partition{s}(:, 2) >= x(s, i)), 1, 'first');
    end
    if z(5, i) == 3
        break;
    end
    z_lin = Sub2Ind([51; 41; 31; 2; 3], z(:, i));
    u_lin = find(min(max(C(z_lin, :, :, :), [], 4), [], 3), 1, 'first');
    u_sub = Ind2Sub([7; 2], u_lin);
    u(1, i) = Controls{1}(u_sub(1));
    u(2, i) = Controls{2}(u_sub(2));
    
    w = Disturbance(1)+rand*12;
    if w > Disturbance(2)
        w = w / 3;
    end
    
    x(:, i+1) = F(x(:, i), x(:, i), u(:, i), w);
    
    
    
    if z(5, i) == 1 
        if x(4, i+1) > 1
          z(5, i+1) = 2;
        else
          z(5, i+1) = 1;
        end
    elseif z(5,i) ==2
        if x(4, i+1) < 1
          z(5, i+1) = 3;
        else
          z(5, i+1) = 2;
        end
    end
end
for s = 1:n_x
    z(s, 101) = find((Partition{s}(:, 1) <= x(s, i)) & (Partition{s}(:, 2) >= x(s, i)), 1, 'first');
end

%% Plot trajectory
N_max = i;
tt = (0:N_max)'*T0;

figure;
subplot(4, 1, 1);
n_takeover1 = find(diff(x(4, :)) > 0, 1, 'first')+1;
n_takeover2 = find(diff(x(4, :)) < 0, 1, 'first');
plot(tt(1:n_takeover1), x(1, 1:n_takeover1), 'k');
hold on;
plot(tt(n_takeover1:n_takeover2), x(1, n_takeover1:n_takeover2), 'b', 'LineWidth', 3);
plot(tt(n_takeover2:N_max), x(1, n_takeover2:N_max), 'k');
plot(tt(1:N_max), -l*ones(N_max, 1), '-.r');
plot(tt(1:N_max), l*ones(N_max, 1), '-.r');
hold off;
grid on
ylabel('x^1', 'FontSize', 12)
ylim([-50, 20]);
legend('First lane', 'Second lane', 'Location', 'southeast', 'FontSize', 10);

subplot(4, 1, 2);
plot(tt(1:N_max), x(2, 1:N_max), 'b');
hold on;
plot(tt(1:N_max), x(3, 1:N_max), 'r');
hold off;
grid on
ylabel('velocities', 'FontSize', 12)
legend('x^2', 'x^3', 'Location', 'southeast', 'FontSize', 10)
ylim([0, 30]);

subplot(4, 1, 3);
plot(tt(1:N_max), u(1, 1:N_max), 'r');
grid on
ylabel('u', 'FontSize', 12)
ylim([-1, 4]);

subplot(4, 1, 4);
plot(tt(1:N_max), z(5,1:N_max), 'b');
grid on
xlabel('t','FontSize',12)
ylabel('p','FontSize',12)
ylim([0.5, 3.5]);




