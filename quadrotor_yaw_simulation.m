%% Report_4.m  -  Quadrotor Yaw Dynamics: P, PI, PD Controller Design (Q3 - Q6)
clear; clc; close all;

% Plant (from Q1 and Q2):
%   G(s) = (s^2 + 0.89 s - 0.03075) / (s^3 - 7.114 s^2 - 1.309 s - 0.355)

%% Q3) P-controller stability check - Kp sweep (Figure 1)
figure;
s = tf('s');                         % Laplace variable
n = [1 0.89 -0.03075];               % numerator
d = [1 -7.114 -1.309 -0.355];        % denominator

kp = -50:1:50;                       % Kp range
stable = zeros(1, length(kp));       % 1 if stable, 0 if not

for k = 1:length(kp)
    K = kp(k);
    Kn = zeros(1, length(d));        % align K*n with denominator degree
    Kn(end-length(n)+1:end) = K*n;
    char_poly = d + Kn;              % characteristic polynomial
    r = roots(char_poly);            % find roots
    if all(real(r) < 0)              % all roots negative real part -> stable
        stable(k) = 1;
    else
        stable(k) = 0;
    end
end

plot(kp, stable, '--bo', 'LineWidth', 2)
xlabel('Kp')
ylabel('Stability (1=Stable,0=Unstable)')
title('P-Controller Stability')
ylim([-0.2 1.2])
grid on

%% Q3) P-controller step response for positive & negative Kp (Figure 2)
figure('Color', 'white');
hold on;
G_plant = tf(n, d);
Kp_list = [-50, -20, 10, 50];
colors  = {'b', 'c', 'r', 'm'};

for i = 1:length(Kp_list)
    Kp = Kp_list(i);
    Sys_CL = feedback(Kp * G_plant, 1);
    t = 0:0.01:2;
    [y, t_out] = step(Sys_CL, t);
    plot(t_out, y, 'LineWidth', 2, 'Color', colors{i}, ...
        'DisplayName', ['Kp = ' num2str(Kp)]);
end

grid on;
title('Step Response for Positive & Negative Kp', 'FontSize', 12, 'FontWeight', 'bold');
subtitle('Proof: System is Unstable for ALL real values of Kp');
xlabel('Time (s)', 'FontWeight', 'bold');
ylabel('Amplitude (Diverging)', 'FontWeight', 'bold');
legend('Location', 'best');
dim = [0.15 0.6 0.3 0.1];
str = {'RESULT: UNSTABLE', 'Reason: Routh-Hurwitz conflict.', 'No single Kp can satisfy', 'all sign conditions.'};
annotation('textbox', dim, 'String', str, 'FitBoxToText', 'on', ...
    'BackgroundColor', 'white', 'EdgeColor', 'red', 'Color', 'red');

%% Q4) PI-controller stability map - Kp & Ki sweep (Figure 3)
figure;
s = tf('s');                         % Laplace variable
n = [1 0.89 -0.03075];               % numerator
d = [1 -7.114 -1.309 -0.355];        % denominator
G = tf(n, d);

Kp_vec = -100:1:100;
Ki_vec = -100:1:100;
stab = zeros(length(Kp_vec), length(Ki_vec));   % stability map

for i = 1:length(Kp_vec)
    for j = 1:length(Ki_vec)
        Kp = Kp_vec(i);
        Ki = Ki_vec(j);

        % PI(s) = (Kp*s + Ki)/s
        % Closed-loop CE: s*D(s) + (Kp*s + Ki)*N(s)
        t1 = [d 0];                  % s*D(s)
        ctrl = [Kp Ki];              % Kp*s + Ki
        t2 = conv(ctrl, n);          % controller * numerator

        % pad for addition
        if length(t1) > length(t2)
            t2 = [zeros(1, length(t1)-length(t2)) t2];
        end

        ce = t1 + t2;                % characteristic equation
        r = roots(ce);               % find roots

        % stable if all real parts < 0
        if all(real(r) < 0)
            stab(i,j) = 1;
        end
    end
end

% Plot map
imagesc(Ki_vec, Kp_vec, stab);
title('PI Controller Stability Map');
xlabel('Ki');
ylabel('Kp');
colorbar;
caxis([0 1]);                        % 0 = unstable, 1 = stable

%% Q4) PI-controller step response (Figure 4)
figure;
G_plant = tf([1 0.89 -0.03075], [1 -7.114 -1.309 -0.355]);
Kp_list = [0.5 1 2];
Ki_list = [0.2 0.5 1 2];
t = 0:0.01:5;

hold on
grid on
legend_text = {};

for i = 1:length(Kp_list)
    for j = 1:length(Ki_list)
        Kp = Kp_list(i);
        Ki = Ki_list(j);

        C_pi = (Kp*s + Ki) / s;
        Sys_CL = feedback(C_pi * G_plant, 1);

        y = step(Sys_CL, t);
        plot(t, y, 'LineWidth', 3.5)
        legend_text{end+1} = ['Kp=', num2str(Kp), ', Ki=', num2str(Ki)];
    end
end

title('Step Response for Different Kp and Ki')
xlabel('Time (s)')
ylabel('Output')
legend(legend_text, 'Location', 'best')
hold off

%% Q5) PD-controller stability map - Kp & Kd sweep (Figure 5)
figure(1);

% plant numerator
n2 = 1;
n1 = 0.89;
n0 = -0.03075;
% plant denominators
d3 = 1;
d2 = -7.114;
d1 = -1.309;
d0 = -0.355;

% gain ranges
kp_values = -150:1:50;               % Kp
kd_values = -160:1:10;               % Kd

% matrix for stability (1=stable, 0=unstable)
Stability = zeros(length(kp_values), length(kd_values));
disp('Checking stability region (Kp-Kd sweep)...');

% main loop
for i = 1:length(kp_values)
    for j = 1:length(kd_values)
        kp = kp_values(i);
        kd = kd_values(j);

        % new poly coefficients
        a3 = d3 + kd*n2;
        a2 = d2 + kd*n1 + kp*n2;
        a1 = d1 + kd*n0 + kp*n1;
        a0 = d0 + kp*n0;

        r = roots([a3 a2 a1 a0]);

        % check stability
        if all(real(r) < 0)
            Stability(i,j) = 1;      % stable
        else
            Stability(i,j) = 0;      % unstable
        end
    end
end

% plot
imagesc(kd_values, kp_values, Stability);
set(gca, 'YDir', 'normal');
xlabel('Derivative Gain (Kd)');
ylabel('Proportional Gain (Kp)');
title(' PD Controller Stability Map');
subtitle('Yellow = Stable   |   Blue = Unstable');
colormap([0 0 0.5; 1 1 0]);          % blue & yellow
colorbar;
grid on;

%% Q5) PD-controller step response (Figure 6)
figure;
n = [1 0.89 -0.03075];
d = [1 -7.114 -1.309 -0.355];
s = tf('s');
G = tf(n, d);

Kp_list = [-10 -5 5 10];
Kd_list = [-5 -1 1 5];
t = 0:0.01:5;

hold on
grid on
legend_text = {};

for i = 1:length(Kp_list)
    for j = 1:length(Kd_list)
        Kp = Kp_list(i);
        Kd = Kd_list(j);
        C_pd = Kp + Kd*s;
        sys_cl = feedback(C_pd * G, 1);
        try
            y = step(sys_cl, t);
            y(abs(y) > 1e3) = NaN;
            plot(t, y, 'LineWidth', 1.2)
            legend_text{end+1} = ['Kp=', num2str(Kp), ', Kd=', num2str(Kd)];
        catch
            continue
        end
    end
end

title('PD Step Response (Kp & Kd Positive and Negative)')
xlabel('Time (s)')
ylabel('Output')
legend(legend_text, 'Location', 'best')
hold off

%% Q6) Best PD controller - 3D optimization landscape (Figure 7)
n2 = 1; n1 = 0.89; n0 = -0.03075;
d3 = 1; d2 = -7.114; d1 = -1.309; d0 = -0.355;

kp_vals = -150:1:50;
kd_vals = -160:1:10;
MAXOS = 20;                          % max allowed overshoot (%)
Ts_Map = nan(length(kp_vals), length(kd_vals));
disp('Calculating 3D Surface from Analytical Roots...');

for i = 1:length(kp_vals)
    for j = 1:length(kd_vals)
        Kp = kp_vals(i);
        Kd = kd_vals(j);

        a3 = d3 + Kd*n2;
        a2 = d2 + Kd*n1 + Kp*n2;
        a1 = d1 + Kd*n0 + Kp*n1;
        a0 = d0 + Kp*n0;

        r = roots([a3 a2 a1 a0]);

        if all(real(r) < 0)
            sigma = max(real(r));
            Ts_calc = 4 / abs(sigma);                 % settling time
            zetas = abs(real(r)) ./ abs(r);
            minzeta = min(zetas);
            if minzeta >= 1
                OScalc = 0;
            else
                OScalc = 100*exp((-minzeta*pi) / sqrt(1-minzeta^2));
            end
            if (Ts_calc < 50) && (OScalc < MAXOS)
                Ts_Map(i,j) = Ts_calc;
            else
                Ts_Map(i,j) = NaN;
            end
        end
    end
end

[min_Ts, idx] = min(Ts_Map(:));
[best_i, best_j] = ind2sub(size(Ts_Map), idx);
best_Kp = kp_vals(best_i);
best_Kd = kd_vals(best_j);

figure('Color', 'white');
[X_Kd, Y_Kp] = meshgrid(kd_vals, kp_vals);
s = surf(X_Kd, Y_Kp, Ts_Map);
s.EdgeColor = 'none';
s.FaceColor = 'interp';
lighting phong;
camlight left;
xlabel('Derivative Gain (Kd)', 'FontWeight', 'bold');
ylabel('Proportional Gain (Kp)', 'FontWeight', 'bold');
zlabel('Settling Time (sec)', 'FontWeight', 'bold');
title(' 3D Optimization Landscape (Analytical)', 'FontSize', 15);
subtitle(['Optimal Point: Kp = ' num2str(best_Kp) ', Kd = ' num2str(best_Kd) ', Ts = ' num2str(min_Ts, '%.2f') 's']);
colormap(jet);
c = colorbar;
c.Label.String = 'Analytical Settling Time (s)';
hold on;
plot3(best_Kd, best_Kp, min_Ts, 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
line([best_Kd best_Kd], [best_Kp best_Kp], [0 min_Ts], 'Color', 'k', 'LineStyle', '--');
text(best_Kd, best_Kp, min_Ts + 2, ' Best Performance', 'FontWeight', 'bold', 'FontSize', 10);
hold off;
view(-45, 30);
grid on;