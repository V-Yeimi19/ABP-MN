% Proyecto de Curso - Métodos Numéricos (CC2104)
% Autores: Grupo 8
% Fecha: Ciclo 2026-0

clear; clc; close all;

%% 1. PARÁMETROS DEL MODELO

p.VR    = 5.0;
p.O2Gin = 0.21;
p.VL0   = 1.0;
p.O2G0  = 0.21;
p.Fin   = 0.05;    % [L h^-1]
p.Vair  = 1.0;
p.S1in  = 500;     % [mmol L^-1] Xilosa en alimentación
p.pres  = 1.0;     % [atm]
p.HO2   = 0.769;   % [atm L mmol^-1] Constante de Henry para O2 a 30 C
p.mu_max = 0.3048;  % [h^-1]
p.Ks     = 99.997;  % [mmol L^-1]
p.Ys     = 0.0861;  % [g mmol^-1]
p.kPX    = 0.01625;
p.vmax1  = 3.7965;  % [mmol g^-1 h^-1]
p.KM1    = 99.995;  % [mmol L^-1]
p.vmax2  = 51.789;  % [mmol L^-1 h^-1] Reacción abiótica
p.vmax3  = 50.121;  % [mmol g^-1 h^-1] Reacción biótica
p.KM2    = 49.979;  % [mmol L^-1]
p.K2     = 0.9545;  % Constante de equilibrio r2/r3
p.vmaxP  = 1.0186;  % [mmol g^-1 h^-1]
p.KM3    = 12.037;  % [mmol L^-1] Para sustrato S1
p.KM4    = 4.9397;  % [mmol L^-1] Para sustrato S2
p.kla    = 422.46;  % [h^-1] Coeficiente volumétrico de transferencia de masa

% Condiciones iniciales: [VL, X, S1, A1, A2, S2, P, O2L, O2G]
O2L_sat0 = p.pres * p.O2G0 / p.HO2;
y0 = [p.VL0; 0.041; 111.0; 0.0; 0.0; 68.4; 0.0; O2L_sat0; p.O2G0];

T        = 42.5;
h_values = [1.5, 1, 0.5, 0.1];

%% 2. USAMOS ODE45 COMO PUNTO DE REFERENCIA
opts_ref = odeset('RelTol',1e-8,'AbsTol',1e-10,'NonNegative',1:9);
f_sys    = @(t,y) sistema_biorreactor(t, y, p);

results = struct;

for ih = 1:length(h_values)

    h      = h_values(ih);
    t_span = (0:h:T)';
    N      = length(t_span);
    [~, Y45] = ode45(f_sys, t_span, y0, opts_ref);
    Y_ref    = Y45';

    %% 3. METODOS NUMERICOS EXPLICITOS

    % Euler explicito
    Y_Euler      = zeros(9, N);
    Y_Euler(:,1) = y0;
    for k = 1:N-1
        dydt           = sistema_biorreactor(t_span(k), Y_Euler(:,k), p);
        Y_Euler(:,k+1) = Y_Euler(:,k) + h*dydt;
        Y_Euler(:,k+1) = aplicar_restricciones(Y_Euler(:,k+1), p);
    end

    % Runge-Kutta 4
    Y_RK4 = RK4_vectorial(f_sys, t_span, y0, h, p);

    % Adams-Bashforth 2 
    Y_AB2      = zeros(9, N);
    Y_AB2(:,1) = y0;
    Y_AB2(:,2) = Y_RK4(:,2);
    f_prev     = sistema_biorreactor(t_span(1), Y_AB2(:,1), p);
    for k = 2:N-1
        f_curr       = sistema_biorreactor(t_span(k), Y_AB2(:,k), p);
        Y_AB2(:,k+1) = Y_AB2(:,k) + (h/2)*(3*f_curr - f_prev);
        f_prev       = f_curr;
        Y_AB2(:,k+1) = aplicar_restricciones(Y_AB2(:,k+1), p);
    end

    %% 4. RMSE vs ode45 (producto P)
    P_ref  = Y_ref(7,:);
    rmse_E = sqrt(mean((Y_Euler(7,:) - P_ref).^2));
    rmse_R = sqrt(mean((Y_RK4(7,:)   - P_ref).^2));
    rmse_A = sqrt(mean((Y_AB2(7,:)   - P_ref).^2));

    fprintf('\n==== h = %.2f h ====\n', h);
    fprintf('%-8s  RMSE_P\n','Metodo');
    fprintf('%-8s  %.6e\n','Euler', rmse_E);
    fprintf('%-8s  %.6e\n','RK4',   rmse_R);
    fprintf('%-8s  %.6e\n','AB2',   rmse_A);

    fprintf('\n  Valores finales (t = %.1f h):\n', T);
    vars = {'VL [L]','X [g/L]','S1 [mmol/L]','S2 [mmol/L]','P [mmol/L]','O2L [mmol/L]'};
    idx  = [1 2 3 6 7 8];
    fprintf('  %-14s  %10s  %10s  %10s\n','Variable','Euler','RK4','AB2');
    for v = 1:length(idx)
        fprintf('  %-14s  %10.4f  %10.4f  %10.4f\n', vars{v}, ...
            Y_Euler(idx(v),end), Y_RK4(idx(v),end), Y_AB2(idx(v),end));
    end

    results(ih).h     = h;
    results(ih).t     = t_span';
    results(ih).Euler = Y_Euler;
    results(ih).RK4   = Y_RK4;
    results(ih).AB2   = Y_AB2;
    results(ih).ref   = Y_ref;
    results(ih).rmse  = [rmse_E, rmse_R, rmse_A];

    %% 5. GRAFICAS POR PASO h
    figure('Name',sprintf('h = %.2f h',h),'Color','w','Position',[100 100 1100 400]);

    subplot(1,3,1)
    plot(t_span,Y_ref(2,:),   'm-', 'LineWidth',2,  'DisplayName','ode45'); hold on
    plot(t_span,Y_Euler(2,:), 'r--','LineWidth',1.2,'DisplayName','Euler')
    plot(t_span,Y_RK4(2,:),   'b-', 'LineWidth',1.2,'DisplayName','RK4')
    plot(t_span,Y_AB2(2,:),   'g-.','LineWidth',1.2,'DisplayName','AB2')
    xlabel('Tiempo (h)'); ylabel('Biomasa X (g L^{-1})')
    title(sprintf('Biomasa  |  h = %.2f', h))
    legend('Location','best'); grid on

    subplot(1,3,2)
    plot(t_span,Y_ref(7,:),   'm-', 'LineWidth',2,  'DisplayName','ode45'); hold on
    plot(t_span,Y_Euler(7,:), 'r--','LineWidth',1.2,'DisplayName','Euler')
    plot(t_span,Y_RK4(7,:),   'b-', 'LineWidth',1.2,'DisplayName','RK4')
    plot(t_span,Y_AB2(7,:),   'g-.','LineWidth',1.2,'DisplayName','AB2')
    xlabel('Tiempo (h)'); ylabel('Producto P (mmol L^{-1})')
    title(sprintf('Producto  |  h = %.2f', h))
    legend('Location','best'); grid on

    subplot(1,3,3)
    plot(t_span,abs(Y_Euler(7,:)-P_ref),'r--','LineWidth',1.2,'DisplayName','|Euler - ref|'); hold on
    plot(t_span,abs(Y_RK4(7,:) -P_ref),'b-', 'LineWidth',1.2,'DisplayName','|RK4 - ref|')
    plot(t_span,abs(Y_AB2(7,:) -P_ref),'g-.','LineWidth',1.2,'DisplayName','|AB2 - ref|')
    xlabel('Tiempo (h)'); ylabel('Error absoluto (mmol L^{-1})')
    title(sprintf('Error vs ode45  |  h = %.2f', h))
    legend('Location','best'); grid on; set(gca,'YScale','log')

end


%% 6. ANALISIS GLOBAL

% Tabla resumen de RMSE
fprintf('\n======= TABLA RESUMEN RMSE =======\n')
fprintf('%-5s  %-14s  %-14s  %-14s\n','h','RMSE_Euler','RMSE_RK4','RMSE_AB2')
for i = 1:length(results)
    fprintf('%-5.2f  %-14.6e  %-14.6e  %-14.6e\n', ...
        results(i).h, results(i).rmse(1), results(i).rmse(2), results(i).rmse(3))
end

% Convergencia de valores finales
nh    = length(results);
P_fin = zeros(nh,3);
X_fin = zeros(nh,3);
for i = 1:nh
    P_fin(i,:) = [results(i).Euler(7,end), results(i).RK4(7,end), results(i).AB2(7,end)];
    X_fin(i,:) = [results(i).Euler(2,end), results(i).RK4(2,end), results(i).AB2(2,end)];
end

figure('Name','Convergencia','Color','w','Position',[100 100 900 400])
subplot(1,2,1)
plot(h_values,X_fin(:,1),'r--o','LineWidth',1.5,'DisplayName','Euler'); hold on
plot(h_values,X_fin(:,2),'b-s', 'LineWidth',1.5,'DisplayName','RK4')
plot(h_values,X_fin(:,3),'g-d', 'LineWidth',1.5,'DisplayName','AB2')
set(gca,'XDir','reverse'); xlabel('Paso h (h)'); ylabel('Biomasa final (g L^{-1})')
title('Convergencia — Biomasa'); legend('Location','best'); grid on

subplot(1,2,2)
plot(h_values,P_fin(:,1),'r--o','LineWidth',1.5,'DisplayName','Euler'); hold on
plot(h_values,P_fin(:,2),'b-s', 'LineWidth',1.5,'DisplayName','RK4')
plot(h_values,P_fin(:,3),'g-d', 'LineWidth',1.5,'DisplayName','AB2')
set(gca,'XDir','reverse'); xlabel('Paso h (h)'); ylabel('Producto final (mmol L^{-1})')
title('Convergencia — Producto'); legend('Location','best'); grid on

% Orden de convergencia (log-log) con lineas de referencia teoricas
rmse_mat = vertcat(results.rmse);
h_log    = linspace(min(h_values)*0.8, max(h_values)*1.2, 50);
c1 = rmse_mat(end,1)/h_values(end)^1;
c2 = rmse_mat(end,3)/h_values(end)^2;
c4 = rmse_mat(end,2)/h_values(end)^4;

figure('Name','Error vs h','Color','w','Position',[100 100 600 450])
loglog(h_values,rmse_mat(:,1),'r--o','LineWidth',1.5,'DisplayName','Euler'); hold on
loglog(h_values,rmse_mat(:,2),'b-s', 'LineWidth',1.5,'DisplayName','RK4')
loglog(h_values,rmse_mat(:,3),'g-d', 'LineWidth',1.5,'DisplayName','AB2')
loglog(h_log,c1*h_log.^1,'c:',  'LineWidth',1,'DisplayName','O(h)')
loglog(h_log,c2*h_log.^2,'m--', 'LineWidth',1,'DisplayName','O(h^2)')
loglog(h_log,c4*h_log.^4,'-','Color',[0.9 0.6 0],'LineWidth',1,'DisplayName','O(h^4)')
xlabel('Paso h (h)'); ylabel('RMSE producto P (mmol L^{-1})')
title('Orden de convergencia (log-log)'); legend('Location','best'); grid on

% Estabilidad numerica para los 4 pasos
figure('Name','Estabilidad numerica','Color','w','Position',[100 100 1100 500])
for i = 1:nh
    subplot(2,2,i)
    plot(results(i).t,results(i).ref(7,:),   'm-', 'LineWidth',2,  'DisplayName','ode45'); hold on
    plot(results(i).t,results(i).Euler(7,:), 'r--','LineWidth',1.2,'DisplayName','Euler')
    plot(results(i).t,results(i).RK4(7,:),   'b-', 'LineWidth',1.2,'DisplayName','RK4')
    plot(results(i).t,results(i).AB2(7,:),   'g-.','LineWidth',1.2,'DisplayName','AB2')
    xlabel('Tiempo (h)'); ylabel('P (mmol L^{-1})')
    title(sprintf('h = %.2f h', results(i).h))
    legend('Location','best'); grid on
end
sgtitle('Estabilidad numerica — Producto P')


%% 7. SENSIBILIDAD PARAMETRICA (Monte Carlo, N=30)
% Se perturban mu_max, Ks y vmax_P con variacion ±10% normal, usando RK4 con h=0.1

fprintf('\n======= SENSIBILIDAD PARAMETRICA (Monte Carlo) =======\n')
nMC     = 30;
h_MC    = 0.1;
P_final = zeros(nMC,1);
rng(42)

for k = 1:nMC
    p_var        = p;
    p_var.mu_max = p.mu_max * (1 + 0.10*randn);
    p_var.Ks     = p.Ks     * (1 + 0.10*randn);
    p_var.vmaxP  = p.vmaxP  * (1 + 0.10*randn);

    f_var      = @(t,y) sistema_biorreactor(t, y, p_var);
    Y_var      = RK4_vectorial(f_var, (0:h_MC:T)', y0, h_MC, p_var);
    P_final(k) = Y_var(7,end);
end

mu_P  = mean(P_final);
std_P = std(P_final);
fprintf('Media P(T):          %.4f mmol/L\n', mu_P)
fprintf('Desv. estandar P(T): %.4f mmol/L\n', std_P)
fprintf('Coef. variacion:     %.2f %%\n',     100*std_P/mu_P)

figure('Name','Sensibilidad parametrica','Color','w','Position',[100 100 600 430])
histogram(P_final, 10, 'FaceColor',[0.2 0.5 0.8], 'EdgeColor','w')
xline(mu_P,       'r-', sprintf('Media = %.3f', mu_P),'LineWidth',2,'LabelVerticalAlignment','bottom')
xline(mu_P+std_P, 'r--','+1\sigma','LineWidth',1)
xline(mu_P-std_P, 'r--','-1\sigma','LineWidth',1)
xlabel('Producto final P(T) (mmol L^{-1})')
ylabel('Frecuencia')
title({'Sensibilidad parametrica (Monte Carlo, N=30)',...
       'Variacion ±10% normal en \mu_{max}, K_s, v_{max,P}'})
grid on


%% FUNCIONES AUXILIARES

function dydt = sistema_biorreactor(~, y, p)
% Sistema de 9 EDOs del biorreactor fed-batch
    VL  = max(y(1),1e-6);  X   = max(y(2),0);
    S1  = max(y(3),0);     A1  = max(y(4),0);
    A2  = max(y(5),0);     S2  = max(y(6),0);
    P   = max(y(7),0);     O2L = max(y(8),0);

    D  = p.Fin / VL;
    mu = p.mu_max * S1 / (p.Ks + S1 + 1e-12);

    r1      = X * p.vmax1 * S1 / (p.KM1 + S1 + 1e-12);
    eq_term = max(A1 - A2/p.K2, 0);
    r2      = p.vmax2 * eq_term;
    r3      = X * p.vmax3 * eq_term / (p.KM2 + eq_term + 1e-12);
    rP      = X * p.vmaxP * S1*S2 / (p.KM3*S1 + p.KM4*S2 + S1*S2 + 1e-12);

    O2L_sat = p.pres * y(9) / p.HO2;  % Ley de Henry
    OTR     = p.kla * (O2L_sat - max(y(8),0));

    dVL  =  p.Fin;
    dX   = -D*X  + X*mu + rP*p.kPX;
    dS1  =  D*(p.S1in - S1) - X*(mu/p.Ys) - r1 - rP;
    dA1  = -D*A1 + r1 - r2 + r3;
    dA2  = -D*A2 + r2 - r3;
    dS2  = -D*S2 - rP;
    dP   = -D*P  + rP;
    dO2L =  D*(O2L_sat - max(y(8),0)) - 5*X*(mu/p.Ys) - rP + OTR;
    dO2G =  (p.Vair*(p.O2Gin - y(9)) - OTR*VL) / max(p.VR - VL, 1e-4);

    dydt = [dVL; dX; dS1; dA1; dA2; dS2; dP; dO2L; dO2G];
end

function y = aplicar_restricciones(y, p)
    y(1)     = min(y(1), p.VR - 1e-4);
    y(2:end) = max(y(2:end), 0);
end

function Y = RK4_vectorial(f, t_span, y0, h, p)
    N      = length(t_span);
    Y      = zeros(length(y0), N);
    Y(:,1) = y0;
    for i = 1:N-1
        k1 = f(t_span(i),       Y(:,i));
        k2 = f(t_span(i)+0.5*h, Y(:,i)+0.5*h*k1);
        k3 = f(t_span(i)+0.5*h, Y(:,i)+0.5*h*k2);
        k4 = f(t_span(i)+h,     Y(:,i)+h*k3);
        Y(:,i+1) = Y(:,i) + (h/6)*(k1+2*k2+2*k3+k4);
        Y(:,i+1) = aplicar_restricciones(Y(:,i+1), p);
    end
end