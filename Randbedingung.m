%% Skript Randbedingung
%  In dem Skript wird fuer ein Fisch N = 1 das Schwimmverhalten geplottet.
%  Damit wird das Ziel verfolgt, zu zeigen, welchen Einfluss die
%  Randbedingungen haben. Daher wird die Orientierungspräferenz \lambda = 0
%  und die Stoerung eta = 0 gewaehlt.
% clear all

% Parameter definieren:
N        = 1;                     % Anzahl der Fische
L_data   = 9.09;                  % Groesse des Beckens in Versuchen
L        = 1;                     % Vereinfacht
                                  % Einheit in BL (daher nicht: 10)
r        = [0.5,1/2 * L, L];      % Radien der Zonengrenzen

% Geschwindigkeit bestimmen anhand der Daten der Einzeltiere
[x_vec,y_vec,r_vec,phi_vec,v_vec,theta_vec,filenames] = read_tab_single_xy;

% Bestimmen des Radius des Beckens in Pixel 
r_vec_zeile = reshape(r_vec,numel(r_vec),1);

% Beruecksichtigung einer 1% Fehlerrate
r_max = quantile(max(r_vec),0.99);  
   
% Umrechnungsfaktor berechnen
uf = L/r_max;

% Durchschnittliche Geschwindigkeit v bestimmen in der Einheit des Radius
v_vec_neu = uf*v_vec;
v         = mean(mean(v_vec_neu))/L_data;  
psi       = @(N) ones(1,N).* -pi/2;% Bevorzugte Richtung psi (frei gewaehlt)
time_sim  = 0:20*60;               % 20 Minuten Beobachtungszeit
dt        = 1;                     % Zeitabstaende 1 ( = 30 sek)

n = 0;                     % Stoerung
lambda = 1;                % Orientierung an alter Schwimmrichtung

% Bestimmen der Startpositionen fuer den Plot
start = [0,0,0;0.5,0.2,0;0.5,0.2,-pi/4;-0.8,-0.4,-pi/4];

for si = 1 : size(start,1)
    x_vec_start = start(si,1);
    y_vec_start = start(si,2);
    theta_vec_start = start(si,3);
    [x_vec,y_vec,phi_vec,theta_vec,r_vec] = ...
        modell_schwarm(N,time_sim,dt,r,L,n,v,lambda,psi(N),x_vec_start,y_vec_start,theta_vec_start);

    % Plot der Trajektorien
    figure(1)
    subplot(2,2,si)
    pl = polarplot(phi_vec,r_vec,'k-.');
    hold on
    polarplot(phi_vec(1),r_vec(1),'rx','MarkerSize',9,'LineWidth',2)
    hold off
    ax = gca;
    ax.RTick =[];

    % Beschriftung der Achsen
    ax = gca;
    ax.RTick = [];
    ax.ThetaTick = [0,45,90,135,180,225,270,315];
    ax.ThetaTickLabel = {'270 °';'315 °';'0 °';'45 °';'90 °';'135 °';'180 °';'225 °'};
    ax.FontSize = 13;

    title(['\theta(0) = ',num2str(circ_rad2ang(theta_vec_start)),'°'],'FontSize',13)
end
saveas(gcf,'Plot_Randbedingungen.svg')