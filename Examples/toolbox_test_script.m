%% Toolbox test script #1
% This example is described in the following paper:
% V. Sinyakov and A. Girard. Formal Controller Synthesis from
% Specifications Given by Discrete-Time Hybrid Automata. 2019. HAL: https://hal.archives-ouvertes.fr/hal-02361404
%% System abstraction
M = 1370;
f_0 = 51.0709;
f_1 = 0.3494;
f_2 = 0.4161;
F1 = @(u, v) (M*u - f_0 - f_1*v - f_2*v^2)*(v > 0)/M + max(M*u - f_0, 0)*(v <= 0)/M;
F2 = @(w, v) (M*w - f_0 - f_1*v - f_2*v^2)*(v > 0)/M + max(M*w - f_0, 0)*(v <= 0)/M;
T0 = 0.25;

n_x = 3;
F = @(x, y, u, w) [x(1) + (x(2)-y(3))*T0; min(max(x(2) + F1(u, y(2))*T0, 0), 30); ...
            min(max(x(3) + F2(w, y(3))*T0, 0), 30)];

Partition = cell(n_x, 1);

tmp = [-Inf, linspace(-100, 0, 51)]';
Partition{1} = [tmp(1:end-1)+eps(tmp(1:end-1)), tmp(2:end)];
Partition{1}(1, 1) = -Inf;

tmp = [0, linspace(0, 30, 31)]';
Partition{2} = [tmp(1:end-1)+eps(tmp(1:end-1)), tmp(2:end)];
Partition{2}(1, 1) = 0;

tmp = linspace(0, 30, 16)';
Partition{3} = [tmp(1:end-1)+eps(tmp(1:end-1)), tmp(2:end)];
Partition{3}(1, 1) = 0;

Controls = cell(1, 1);
Controls{1} = [-20, -15, -10, -5, 10, 5, 0];
Disturbance = [-20, 10];

tic;
[TS, TransitionNumber] = AbstractSystemDiscrete(F, Partition, Controls, Disturbance);
toc;

%% Specification abstraction
delta = 3;
c = 0.1;
l1 = 2;
l2 = 50;

Spec.isMonotone = true;
Spec.Inputs = [5, 10, 15, 20, 25];
N_inputs = length(Spec.Inputs);
Np = 2;
Spec.W = cell(Np, Np, N_inputs);
Spec.D = cell(Np, Np, N_inputs);
for i = 1:N_inputs
    Spec.W{1, 1, i} = @(x, x_plus) (-x_plus(2)+x(2)+c)*(Spec.Inputs(i)-x(2)>=delta) + (x_plus(2)-x(2)+c)*(x(2)-Spec.Inputs(i)>=delta);
    Spec.D{1, 1, i} = @(x) x(1)+l2;
    Spec.W{1, 2, i} = @(x, x_plus) (x_plus(2)-x(2)+c)*(x(2)-Spec.Inputs(i)>=delta);
    Spec.D{1, 2, i} = @(x) max(x(1)+l1, -x(1)-l2);
    Spec.W{2, 1, i} = @(x, x_plus) (-x_plus(2)+x(2)+c)*(Spec.Inputs(i)-x(2)>=delta) + (x_plus(2)-x(2)+c)*(x(2)-Spec.Inputs(i)>=delta);
    Spec.D{2, 1, i} = @(x) x(1)+l2;
    Spec.W{2, 2, i} = @(x, x_plus) (x_plus(2)-x(2)+c)*(x(2)-Spec.Inputs(i)>=delta);
    Spec.D{2, 2, i} = @(x) max(x(1)+l1, -x(1)-l2);
end

tic;
TS_with_Spec = AbstractSpecificationDiscrete(Spec, Partition, TS, TransitionNumber);
toc;

%% Controller computation
tic;
[V, C, domC] = ComputeControllerSafety(TS_with_Spec, []);
toc;

%% Controllable sets visualization
[XX, YY, ZZ] = meshgrid(linspace(Partition{2}(2), Partition{2}(end)-1, 30)+1/2, ...
    linspace(Partition{1}(2), Partition{1}(end)-1, 50)+1/2, ...
    linspace(Partition{3}(1), Partition{3}(end)-1, 15)+1/2);

VV = reshape(V(1:(end-1)), [51, 31, 15, 2, 5]);
VV = max(VV(2:end, 2:end, :, 1, :), [], 5);
VV(VV == Inf) = 1000;
figure;
h = patch(isosurface(XX, YY, ZZ, VV, 200));
set(h, 'EdgeColor', 'None');
set(h, 'FaceColor', [1, 0, 0]);
set(h, 'EdgeLighting', 'phong');
set(h, 'FaceLighting', 'phong');
set(h, 'FaceAlpha', 1);
set(h, 'DiffuseStrength', 0.6);
set(h, 'SpecularStrength', 0.9);
VV(1, :, :) = 200+1e-4;
VV(:, 1, :) = 200+1e-4;
VV(:, end, :) = 200+1e-4;
VV(:, :, 1) = 200+1e-4;
VV(:, :, end) = 200+1e-4;
h1 = patch(isosurface(XX, YY, ZZ, VV, 200));
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
view(32, 32);
camlight('left');
view([-50, 38]);

%% Trajectory simulation
rng(14042020);
x = zeros(3, 201);
x(:, 1) = [-80; 28; 25];
z = zeros(5, 201);
z(:, 1) = [1; 1; 1; 1; 1];
u = zeros(1, 200);
u_lin = 0;
for i = 1:200
    for s = 1:n_x
        z(s, i) = find((Partition{s}(:, 1) <= x(s, i)) & (Partition{s}(:, 2) >= x(s, i)), 1, 'first');
    end
    z_lin = Sub2Ind([51; 31; 15; 2; 5], z(:, i));
    
    u_ind = find(min(max(C(z_lin, :, :, :), [], 4), [], 3), 1, 'last');
    u(i) = Controls{1}(u_ind);
    
    if x(3, i) > 5
        w = Disturbance(1)+rand*diff(Disturbance);
    else
        w = rand*Disturbance(2);
    end
    
    x(:, i+1) = F(x(:, i), x(:, i), u(i), w);
    if i <= 100
        z(5, i+1) = 1;
    else
        z(5, i+1) = 5;
    end
    if x(1, i+1)+l2 <= 0
        z(4, i+1) = 1;
    else
        z(4, i+1) = 2;
    end
end
for s = 1:n_x
    z(s, 201) = find((Partition{s}(:, 1) <= x(s, i)) & (Partition{s}(:, 2) >= x(s, i)), 1, 'first');
end

%% Plot trajectory
tt = (0:200)'*T0;

figure;
subplot(3, 1, 1);
plot(tt, x(1, :), 'k');
hold on
plot(tt, -50*ones(201, 1), '-.r');
hold off
grid on
xlabel('t','FontSize',12)
ylabel('relative distance','FontSize',12)
legend('x^1', 'x^1+l_2=0', 'Location', 'northwest','FontSize',10);
ylim([-150, 0]);

%figure;
subplot(3, 1, 2);
plot(tt, x(2, :), 'b');
hold on;
plot(tt, x(3, :), 'r');
plot(tt, Spec.Inputs(z(5, :)), 'k');
hold off;
grid on
xlabel('t','FontSize',12)
ylabel('velocities','FontSize',12)
ylim([-5, 35]);
legend('x^2', 'x^3', 'v', 'Location', 'northwest','FontSize',10);

%figure;
subplot(3, 1, 3);
plot(tt(1:end-1), u, 'r');
grid on
xlabel('t','FontSize',12)
ylabel('u','FontSize',12)
ylim([-25, 15]);








