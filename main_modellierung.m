% Skript main_modellierung
% clear all

%% Festlegen von grundlegenden Paramter
% Farbschema der Arbeit
col1 = '#6495ED';
col2 = '#000080';

% Groesse der Schrift, Nummern und Dicke der Linien
letter_size = 15; 
number_size = 13;
line_size   = 2; 

% Definieren von Parameter fuer die Simulatione eines einzelnen Fisches
% oder einer Fischgruppe
psi      = @(N) ones(1,N).* -pi/2;% Bevorzugte Richtung psi (frei gewaehlt)
time_sim = 0:20*60;               % 20 Minuten Beobachtungszeit
dt       = 1;                     % Zeitabstand (= 30 sek)
iter_end = 100;                   % Anzahl an Wiederholungen (Iterationen)
n_vec    = linspace(0.01,pi,101); % Intervall, in dem eta variiert wird             
l_vec    = 0:0.01:1;              % Intervall, in dem lambda variiert wird

%% main_simulation_einzel
% Simulationen von dem Schwimmverhalten eines einzelnen Fisches (N = 1). 
% Dafuer werden anhand der zugrundeliegenden Daten die Parameter bestimmt.
% Fuer die Parameter eta und lambda werden Werte in einem bestimmten
% Intervall gewaehlt und jeweils simuliert. 

%% Parameter definieren fuer die Simulation eines einzelnen Fisches
%  Einlesen der Daten
[x_vec,y_vec,r_vec,phi_vec,v_vec,theta_vec,filenames] = read_tab_single_xy;

% Bestimmen des Radius des Beckens in Pixel 
r_all= reshape(r_vec,numel(r_vec),1);

% Beruecksichtigung einer 1% Fehlerrate
r_max = quantile(max(r_vec),0.99);           

% Zur Veranschaulichung: Plot der Verteilung der Radien r (Anhang)
figure(1)
histogram(r_all,'FaceColor',col1,'EdgeColor',col2)
hold on 
xline(r_max,'-','Festglegeter Radius r','Color',col1)
hold off
xlabel('Größe der Radien r'); ylabel('Anzahl der Radien r')


%% Definieren  von Parameter fuer die Simulation eines einzelnen Fisches
N        = 1;                     % Anzahl der Fische
L        = 9.09;                  % Groesse des Beckens (Einheit 
                                  % in BL (daher nicht: 10))
r        = [0.5,1/2 * L, L];      % Radien der Zonengrenzen

   
% Umrechnungsfaktor berechnen
uf = L/r_max;

% Durchschnittliche Geschwindigkeit v bestimmen in der Einheit des Radius
v_vec_neu = uf*v_vec;
v         = mean(mean(v_vec_neu));  

%% 1) Variieren von dem Parameter der Stoerung eta
% 1.1) Plot der durchschnittlichen p-Werte fuer das jeweilige eta                    

% Festlegen von dem Parameter lambda
lambda = 0;                          % Fisch orientiert sich an seiner eigenen
                                     % Orientierungspraeferenz

% Erstellen von leeren Vektoren zum Abspeichern der p-Werte und
% Ordnungsparameter fuer jedes eta und jeder Wiederholung
p_n_vec    = NaN(iter_end,length(n_vec));
ord_n_vec  = NaN(iter_end,length(n_vec));
phi_diff_n = NaN(iter_end,length(n_vec));

% Fuer jede Wiederholung und fuer jeden Wert von eta wird der p-Wert
% mithilfe des Rayleigh Tests bestimmt und der Ordnungsparameter
for i_n = 1 : length(n_vec)
    n   = n_vec(i_n); 

    for i_iter = 1 : iter_end
        [x_vec,y_vec,phi_vec,theta_vec] = ...
            modell_schwarm(N,time_sim,dt,r,L,n,v,lambda,psi(N)); 
        p_n_vec(i_iter,i_n)    = circ_rtest(phi_vec);
        ord_n_vec(i_iter,i_n)  = circ_r(phi_vec);
        phi_diff_n(i_iter,i_n) = circ_dist(circ_mean(phi_vec),psi(N));
    end
 
end

% Plot der Graphen fuer die p-Werte und Ordnungsparameter ueber eta
figure(2)
% Plot der p-Werte an der linken y-Achse
yyaxis left
plot(n_vec,mean(p_n_vec,1,'omitnan'),'b-','LineWidth',2)
xlabel('\eta','FontSize',letter_size)
ylabel('{\it p}','FontSize',letter_size)
xlim([0,pi]); ylim([0 1])
ax          = gca;
ax.YColor   = 'b';
ax.FontSize = 13; 

% Plot der Ordnungsparameter an der rechten y-Achse
yyaxis right
plot(n_vec,mean(ord_n_vec,1,'omitnan'),'r-','LineWidth',line_size)
ylabel('$\overline{R}_{\phi}$', 'Interpreter', 'latex','FontSize',letter_size)
ylim([0 1]); xlim([0,pi])
ax            = gca;
ax.YColor     = 'r';
ax.XTick      = [0,1/3*pi,2/3*pi,pi,4/3*pi,5/3*pi,2*pi];
ax.XTickLabel = {'0','1/3 \pi','2/3 \pi','\pi','4/3 \pi','5/3 \pi','2 \pi'};
ax.FontSize   = 13; 
%legend('{\it p}','Ordnungsparameter','Location','southwest','FontSize',12)

% Speichern der Abbildung
saveas(gcf,'P_Ord_n_M_E','svg')

figure(3)
plot(n_vec,circ_rad2ang(circ_mean(phi_diff_n)),'k-','LineWidth',line_size)
xlabel('\eta','FontSize',letter_size)
ylabel('$(\overline{\phi} - \psi)$ in $^{\circ}$','Interpreter', 'latex','FontSize',letter_size)
xlim([0,pi]); ylim([0 360])
ax = gca; 
ax.FontSize = letter_size;

% Speichern der Abbildung
saveas(gcf,'Diff_n_E','svg') 

% 1.2) Plot der Verteilung der Trajektorien
% 4 Werte fuer eta 
n_vec_plot = 0:(1/3)*pi:pi;

% Fuer jeden Werte von eta werden die simulierten Positionen des Fisches
% geplottet
figure(4)
for i_n = 1 : length(n_vec_plot)
    n = n_vec_plot(i_n);
    [x_vec,y_vec,phi_vec,theta_vec,r_vec] = modell_schwarm(N,time_sim,dt,r,L,n,v,lambda,psi(N));

    % Plot der Trajektorien
    subplot(2,2,i_n)
    polarplot(phi_vec(:,:),r_vec(:,:),'-','Color','b')
    hold on
    polarplot(psi(N),L,'r*','LineWidth',line_size,'MarkerSize',10)
%     hold on
%     polarplot([0 (circ_mean(phi_vec))],[0 (circ_r(phi_vec))*L],'k-','LineWidth',1)
%     hold on
    polarplot(0:2*pi,L,'k-')
    hold off

    % Beschriftung der Achsen
    ax                = gca;
    ax.RTick          = [];
    ax.ThetaTick      = [0,45,90,135,180,225,270,315];
    ax.ThetaTickLabel = {'270 °';'315 °';'0 °';'45 °';'90 °';'135 °';'180 °';'225 °'};
    ax.FontSize       = number_size; 

    % Hinzufuegen eines Titels
    if n == n_vec_plot(1) || n == n_vec_plot(end)
        title(['\eta = ',num2str(n/pi,'%.0f'),' \pi'])
    else 
        title(['\eta = ',num2str(n/pi,'%.1f'),' \pi'])
    end
end

% Speichern der Abbildung
saveas(gcf,'Trajekt_n_M_E','svg')
%% 2) Variieren von dem Parameter der Orientierungspraeferenz lambda 
% 2.1) Plot der durchschnittlichen p-Werte fuer lambda                                     
% Festlegen von eta (Stoerung)
n = pi/10;                         

% Erstellen von leeren Vektoren zum Abspeichern der p-Werte und
% Ordnungsparameter fuer jedes lambda und jede Wiederholung
p_l_vec    = NaN(iter_end,length(l_vec));
ord_l_vec  = NaN(iter_end,length(l_vec));
phi_diff_l = NaN(iter_end,length(l_vec));

% Fuer jede Wiederholung und fuer jedes lambda wird der p-Wert
% mithilfe des Rayleigh Tests bestimmt und der Ordnungsparameter
for i_l = 1 : length(l_vec)
    lambda = l_vec(i_l);

    for i_iter = 1 : iter_end
        [~,~,phi_vec,~,~]      = modell_schwarm(N,time_sim,dt,r,L,n,v,lambda,psi(N));
        p_l_vec(i_iter,i_l)    = circ_rtest(phi_vec); 
        ord_l_vec(i_iter,i_l)  = circ_r(phi_vec);
        phi_diff_l(i_iter,i_l) = circ_dist(circ_mean(phi_vec),psi(N));
    end 

end

% Plot der Graphen fuer die p-Werte und Ordnungsparameter ueber eta
figure(5)
% Plot der p-Werte an der linken y-Achse
yyaxis left
plot(l_vec,mean(p_l_vec,1,'omitnan'),'b-','LineWidth',line_size)
xlabel('\lambda','FontSize',letter_size);ylabel('{\it p}','FontSize',letter_size)
ylim([0 1])
ax          = gca;
ax.YColor   = 'b';
ax.FontSize = number_size; 

% Plot des Ordnungsparameters an der rechten y-Achse
yyaxis right
plot(l_vec,mean(ord_l_vec,1,'omitnan'),'r-','LineWidth',line_size)
ylabel('$\overline{R}_{\phi}$', 'Interpreter', 'latex','FontSize',letter_size)
xlim([0,1]); ylim([0 1])
ax          = gca;
ax.YColor   = 'r';
ax.FontSize = number_size;
% legend(' {\it p}','Ordnungsparameter','Location','southwest','FontSize',12)

% Speichern der Abbildung
saveas(gcf,'p_Ord_M_E','svg')

figure(6)
plot(l_vec,circ_rad2ang(circ_mean(phi_diff_l)),'k-','LineWidth',line_size)
xlabel('\lambda','FontSize',letter_size)
ylabel('$(\overline{\phi} - \psi)$ in $^{\circ}$','Interpreter', 'latex','FontSize',letter_size)
xlim([0,1]); ylim([0 360])
ax          = gca;
ax.FontSize = number_size;

% Speichern der Abbildung
saveas(gcf,'Diff_l_E','svg') 

% 2.2) Plot der Verteilung der Trajektorien
% 4 Werte fuer lambda
l_vec_plot = 0:1/3:1;

% Fuer jeden Werte von lambda werden die simulierten Positionen des Fisches
% geplottet
figure(7)
for i_l = 1 : length(l_vec_plot)
        lambda = l_vec_plot(i_l);
        [ ~,~,phi_vec,~,r_vec] = modell_schwarm(N,time_sim,dt,r,L,n,v,lambda,psi(N));

        % Plot der Trajektorien
        subplot(2,2,i_l)
        polarplot(phi_vec(:,:),r_vec(:,:),'-','Color','b')
        hold on
        polarplot(psi(N),L,'r*','LineWidth',line_size,'MarkerSize',10)
%         hold on
%         polarplot([0 (circ_mean(phi_vec))],[0 (circ_r(phi_vec))*L],'k-','LineWidth',1)
%         hold on
        polarplot(0:2*pi,L,'k-')
        hold off

        % Beschriftung der Achsen
        ax                = gca;
        ax.RTick          = [];
        ax.ThetaTick      = [0,45,90,135,180,225,270,315];
        ax.ThetaTickLabel = {'270 °';'315 °';'0 °';'45 °';'90 °';'135 °';'180 °';'225 °'};
        ax.FontSize       = number_size; 

        % Hinzufuegen eines Titels
        title(['\lambda = ',num2str(round(lambda,1))])
        hold on
end

% Speichern der Abbildung
saveas(gcf,'Trajekt_l_M_E','svg')
%% 3) Bestimmung des p-Wertes und Ordnungsparameter fuer variierende Werte 
% fuer lambda und eta
% Vektoren mit Intervallen fuer die Parameter der Stoerung eta und der 
% Ordnungspräferenz lambda

% Vektoren zum Befuellen fuer jede Wiederholung
p_data_iter        = NaN(1,iter_end); 
ord_data_iter      = NaN(1,iter_end); 
phi_diff_data_iter = NaN(1,iter_end);

% Vektoren zum Befuellen der Mittelwerte pro Wiederholung fuer den
% jeweiligen Wert von eta und lambda
p_data_n_l        = NaN(length(n_vec),length(l_vec));
ord_data_n_l      = NaN(length(n_vec),length(l_vec));
phi_diff_data_n_l = NaN(length(n_vec),length(l_vec));

for i_n = 1 : length(n_vec)
    n = n_vec(i_n);

    for i_l = 1 : length(l_vec)
        lambda = l_vec(i_l);

        for i_iter = 1: iter_end
            [~,~,phi_vec,~,~]          = modell_schwarm(N,time_sim,dt,r,L,n,v,lambda,psi(N));
            p_data_iter(i_iter)        = circ_rtest(phi_vec);
            ord_data_iter(i_iter)      = circ_r(phi_vec);
            phi_diff_data_iter(i_iter) = abs(circ_mean(circ_dist(psi(N),phi_vec)));
        end

        p_data_n_l(i_n,i_l)        = mean(p_data_iter);
        ord_data_n_l(i_n,i_l)      = mean(ord_data_iter);
        phi_diff_data_n_l(i_n,i_l) = circ_rad2ang(circ_mean(phi_diff_data_iter));
    end

end

%% Plot der Matrix der p-Werte fuer die variierenden Werte fuer lambda und
% eta
figure(8)
pl = pcolor(l_vec,n_vec,p_data_n_l);
colormap jet
pl.FaceColor     = 'interp';
c                = colorbar; 
c.Title.String   = '{\it p}';
c.Title.FontSize = letter_size; 

xlabel('\lambda','FontSize',letter_size)
ylabel('\eta','FontSize',letter_size)
yticks([0,1/3*pi,2/3*pi,pi,4/3*pi,5/3*pi,2*pi]);
yticklabels({'0','1/3 \pi','2/3 \pi','\pi','4/3 \pi','5/3 \pi','2 \pi'});
ax          = gca;
ax.FontSize = number_size; 

% Speichern der Abbildung
saveas(gcf,'p_n_l_M_E','svg')

% Plot der Matrix der Ordnungsparameter fuer die variierenden Werte fuer lambda und
% eta
figure(9)
colormap jet
oldcmap             = colormap;
colormap( flipud(oldcmap) )
pl                  = pcolor(l_vec,n_vec,ord_data_n_l);
pl.FaceColor        = 'interp';
c                   = colorbar; 
c.Ticks             = [0,0.2,0.4,0.6,0.8];
c.Title.String      = '$\overline{R}_{\phi}$';
c.Title.Interpreter = 'latex';
c.Title.FontSize    = letter_size; 

xlabel('\lambda','FontSize',letter_size)
ylabel('\eta','FontSize',letter_size)
yticks([0,1/3*pi,2/3*pi,pi,4/3*pi,5/3*pi,2*pi]);
yticklabels({'0','1/3 \pi','2/3 \pi','\pi','4/3 \pi','5/3 \pi','2 \pi'});
ax          = gca;
ax.FontSize = number_size; 

% Speichern der Abbildung
saveas(gcf,'Ord_n_l_M_E','svg')

% Plot der Differenzen
figure(10)
colormap jet
pl                  = pcolor(l_vec,n_vec,phi_diff_data_n_l);
pl.FaceColor        = 'interp';
c                   = colorbar; 
c.Title.String      = '$(\overline{\phi} - \psi)$ in $^{\circ}$';
c.Title.Interpreter = 'latex';
c.Title.FontSize    = letter_size;
xlabel('\lambda','FontSize',letter_size)
ylabel('\eta','FontSize',letter_size)
yticks([0,1/3*pi,2/3*pi,pi,4/3*pi,5/3*pi,2*pi]);
yticklabels({'0','1/3 \pi','2/3 \pi','\pi','4/3 \pi','5/3 \pi','2 \pi'});


% Speichern der Abbildung
saveas(gcf,'Diff_l_n_M_E','svg')

%% Skript: main_simulation_vier
%
%  In dem Skript wird das Modell fuer N = 4 simuliert. Dafuer werden die
%  Parameter an den Daten angepasst. Die Parameter eta und lambda werden
%  fuer Werte eines bestimmten Intervalls variiert. Fuer jeweilige Werte
%  von lambda und eta wird der Durchschnitt von den Positionswinkel phi
%  bestimmt, davon die Differenz zu der bevorzugten Richtung psi des
%  informierten Fisches und die durchschnittliche Resultatlänge von phi.
%  Dieses Groessen werden anschliessend ueber die Parameter eta und lambda
%  jeweils geplottet.
% 

%% Parameter fuer die Simulation einer Fischgruppe
% Parameter anhand den gegebenen Informationen bestimmen

N        = 4;                % Anzahl der Fische
L        = 13.18;            % Radius der Becken
r        = [0.5,1/2 * L, L]; % Radien der Zonengrenzen 
psi_vec = psi(N); 

%% Parameter v anhand der Daten bestimmen (nur 80%, um overfitting zu meiden)
% Einlesen der Daten der 4-er Gruppen
cd 'Daten'\Fischgruppen
filenames_all= dir('*4.D.*.csv');
cd ../..

% Bestimmen der maximalsten Laenge eines Radius eines Fisches in allen
% Versuchen 
v_all = NaN(iter_end*N*42,1); r_all = NaN(iter_end*N*42,1);

% Bestimmmen der durchschnittlichen Geschwindigkeit und der bevorzugten 
% Position aller Beobachtungen der 4-er Gruppe
v_mean   = NaN(size(filenames_all));   % Geschwindigkeit
phi_mean = NaN(size(filenames_all));   % Bevorzugte Richtung (Position)

for i_iter = 1 : length(filenames_all)
    % Fuer das jeweilie file die Datein einlesen
    file = filenames_all(i_iter).name;
    [~,~,r_vec,phi_vec,v_vec,~,~,~] = ...
         read_tab_group(file); 

    % Durchschnitt der Geschwindigkeit v und der Winkel der Position phi 
    % bestimmen
    r_all(1+(i_iter-1)*N*42:(i_iter-1)*N*42 + numel(r_vec),1) = reshape(r_vec,numel(r_vec),1);
    v_all(1+(i_iter-1)*N*42:(i_iter-1)*N*42 + numel(v_vec),1) = reshape(v_vec,numel(v_vec),1);
end 

% Bestimmen des "maximalen" Radius
r_max = quantile(r_all,0.99);

% Zur Veranschaulichung: Plot der Verteilung der Radien r (Anhang)
figure(11)
histogram(r_all,'FaceColor',col1,'EdgeColor',col2)
hold on 
xline(r_max,'-','Festglegeter Radius r','Color',col1)
hold off
xlabel('Größe der Radien r'); ylabel('Anzahl der Radien r')

% Umrechnungsfaktor bestimmen
uf = L/r_max;

% Bestimmen des Durchschnitts der Geschwindigkeiten 
v_new = v_all * uf;
v = mean(v_new,'omitnan');

%% Variieren des Parameter eta: Plot Trajektorien
n_vec_plot = 0:(1/3)*pi:pi;
lambda = [0,1,1,1];

figure(12)
tiledlayout flow
for i_n = 1: length(n_vec_plot)
    n = n_vec_plot(i_n);
    [~,~,phi_vec,~,r_vec] = modell_schwarm(N,time_sim,dt,r, L,n,v,lambda,psi(N));

    % Plot des Ergebnisses
    nexttile
       polarplot(phi_vec(:,4),r_vec(:,4),'Color',col1,'LineStyle','-')
       hold on
       polarplot(phi_vec(:,2),r_vec(:,2),'Color',col1,'LineStyle','-')
       hold on
       polarplot(phi_vec(:,3),r_vec(:,3),'Color',col1,'LineStyle','-')
       hold on
       polarplot(phi_vec(:,1),r_vec(:,1),'Color',col2,'LineStyle','-')
       hold on
       polarplot(0:0.01:2*pi,L,'k-')
%        hold on 
%        polarplot([0 (circ_mean(reshape(phi_vec,numel(phi_vec),1)))],...
%            [0 (circ_r(reshape(phi_vec,numel(phi_vec),1)))*L],'k-','LineWidth',1)
%        hold on
       polarplot(psi_vec(1),L,'r*','LineWidth',line_size,'MarkerSize',10)
       hold off

       % Beschriftung der Achsen
       ax                = gca;
       ax.RTick          = [];
       ax.ThetaTick      = [0,45,90,135,180,225,270,315];
       ax.ThetaTickLabel = {'270 °';'315 °';'0 °';'45 °';'90 °';'135 °';'180 °';'225 °'};
       ax.FontSize       = 13; 
       if n == n_vec(1) || n == n_vec(end)
           title(['\eta = ',num2str(n/pi,'%.0f'),' \pi'])
       else
           title(['\eta = ',num2str(n/pi,'%.1f'),' \pi'])
       end
    
end
% Speichern der Abbildung
saveas(gcf,'Trajekt_n_M_4','svg')

%% Variieren des Parameters eta: Mehrere Iterationen
% Festlegen von lambda
lambda           = [0,1,1,1]; 
phi_mean_eta_all = NaN(size(n_vec));
phi_diff_eta_all = NaN(size(n_vec));
ord_eta_all      = NaN(size(n_vec));
p_eta_all        = NaN(size(n_vec));

% Fuer jedes eta wird der Durchschnitt von dem p-Wert und Ordnungsparameter
% über alle Wiederholungen berechnet
for  i_n = 1 : length(n_vec)
     n        = n_vec(i_n);
     phi_mean = NaN(iter_end,1);
     ord      = NaN(iter_end,1);
     phi_diff = NaN(iter_end,1);
     p_wert   = NaN(iter_end,1);

    for i_iter = 1:iter_end
        [~,~,phi_vec,~,~] =  modell_schwarm(N,time_sim,dt,r, L,n,v,lambda,psi(N));
        phi_mean(i_iter)  = circ_mean(reshape(phi_vec,numel(phi_vec),1));
        phi_diff(i_iter)  = abs(circ_dist(psi_vec(1),phi_mean(i_iter)));
        ord(i_iter)       = circ_r(reshape(phi_vec,numel(phi_vec),1));
        p_wert(i_iter)    = circ_rtest(reshape(phi_vec,numel(phi_vec),1));
    end

    phi_mean_eta_all(i_n) = circ_rad2ang(circ_mean(phi_mean));
    phi_diff_eta_all(i_n) = circ_rad2ang(circ_mean(phi_diff));
    ord_eta_all(i_n)      = mean(ord);
    p_eta_all(i_n)        = mean(p_wert);
end

%% Plotten der Ergebnisse
% Plot des p-Wertes ueber eta an der linken y-Achse
figure(13)
yyaxis left
plot(n_vec,p_eta_all,'b-','LineWidth',line_size)
xlabel('\eta','FontSize',letter_size)
ylabel('{\it p}','FontSize',letter_size)
ylim([0 1])
ax        = gca;
ax.YColor = 'b';
ax.FontSize = number_size;

% Plot der Werte fuer den Ordnungsparameter an der rechten y-Achse
yyaxis right
plot(n_vec,ord_eta_all,'r-','LineWidth',line_size)
ylabel('$\overline{R}_{\phi}$', 'Interpreter', 'latex','FontSize',letter_size)
xlim([0,pi]); ylim([0 1])
ax            = gca;
ax.YColor     = 'r';
ax.XTick      = [0,1/3*pi,2/3*pi,pi,4/3*pi,5/3*pi,2*pi];
ax.XTickLabel = {'0','1/3 \pi','2/3 \pi','\pi','4/3 \pi','5/3 \pi','2 \pi'};
ax.FontSize   = number_size;
% legend('{\it p}','Ordnungsparameter')

% Speichern der Abbildung
saveas(gcf,'Ordp_n_M_4','svg')

% Plot der Differenz der Winkel phi zu dem Winkel der
% Orientierungspraeferenz
figure(14)
plot(n_vec,phi_diff_eta_all,'k-','LineWidth',line_size)
xlabel('\eta','FontSize',letter_size)
ylabel('$(\overline{\phi} - \psi)$ in $^{\circ}$','Interpreter', 'latex','FontSize',letter_size)
xlim([0 pi]); ylim([0 360])
ax            = gca; 
ax.XTick      = [0,1/3*pi,2/3*pi,pi,4/3*pi,5/3*pi,2*pi];
ax.XTickLabel = {'0','1/3 \pi','2/3 \pi','\pi','4/3 \pi','5/3 \pi','2 \pi'};
ax.FontSize   = number_size;

% Speichern der Abbildung
saveas(gcf,'Diff_n_M_4','svg')
%% Variieren des Parameter lambda: Plot Trajektorien
l_vec_plot = 0:1/3:1;

% Festlegen von eta (Stoerung)
n = pi/10;

figure(15)
tiledlayout flow
for i_lambda = 1: length(l_vec_plot)
    lambda = [l_vec_plot(i_lambda),1,1,1];
    [~,~,phi_vec,~,r_vec] = modell_schwarm(N,time_sim,dt,r, L,n,v,lambda,psi(N));

    % Plot des Ergebnisses
    figure(15)
    nexttile
       polarplot(phi_vec(:,4),r_vec(:,4),'Color',col1,'LineStyle','-')
       hold on
       polarplot(phi_vec(:,2),r_vec(:,2),'Color',col1,'LineStyle','-')
       hold on
       polarplot(phi_vec(:,3),r_vec(:,3),'Color',col1,'LineStyle','-')
       hold on
       polarplot(phi_vec(:,1),r_vec(:,1),'Color',col2,'LineStyle','-')
       hold on
       polarplot(0:0.01:2*pi,L,'k-')
%        hold on 
%        polarplot([0 circ_mean(circ_mean(phi_vec)')],[0 mean(circ_r(phi_vec))*L],'k-','LineWidth',1)
%        hold on
       polarplot(psi(N),max(max(r_vec)),'r*','LineWidth',line_size,'MarkerSize',10)
       hold off

       % Beschriftung der Achsen
       ax                = gca;
       ax.RTick          = [];
       ax.ThetaTick      = [0,45,90,135,180,225,270,315];
       ax.ThetaTickLabel = {'270 °';'315 °';'0 °';'45 °';'90 °';'135 °';'180 °';'225 °'};
       ax.FontSize       = number_size; 
    title([' \lambda = ', num2str(l_vec_plot(i_lambda),'%.1f')])
end

% Speichern der Abbildung
saveas(gcf,'Trajekt_l_M_4','svg')
%% Variieren des Parameters lambda: Mehrere Iterationen
phi_mean_all_lambda = NaN(size(l_vec));
phi_diff_all_lambda = NaN(size(l_vec));
r_strich_all_lambda = NaN(size(l_vec));
p_wert_all_lambda   = NaN(size(l_vec));

for  i_lambda = 1 : length(l_vec)
    lambda    = [l_vec(i_lambda),1,1,1];
    phi_mean  = NaN(iter_end,1);
    ord       = NaN(iter_end,1);
    phi_diff  = NaN(iter_end,1);
    p_wert    = NaN(iter_end,1);

    for i_iter = 1:iter_end
        [~,~,phi_vec,~]    =  modell_schwarm(N,time_sim,dt,r, L,n,v,lambda,psi(N));
        phi_mean(i_iter)   = circ_mean(reshape(phi_vec,numel(phi_vec),1));
        phi_diff(i_iter)   = abs(circ_mean(circ_dist(psi(1),reshape(phi_vec,numel(phi_vec),1))));
        ord(i_iter)        = circ_r(reshape(phi_vec,numel(phi_vec),1));
        [p_wert(i_iter),~] = circ_rtest(reshape(phi_vec,numel(phi_vec),1));
    end

    phi_mean_all_lambda(i_lambda) = circ_rad2ang(circ_mean(phi_mean));
    phi_diff_all_lambda(i_lambda) = circ_rad2ang(circ_mean(phi_diff));
    r_strich_all_lambda(i_lambda) = mean(ord);
    p_wert_all_lambda(i_lambda)   = mean(p_wert,'all');
end

%% Plotten der Ergebnisse
% Plot der p-Werte ueber lambda an der linken y-Achse
figure(16)
yyaxis left
plot(l_vec,p_wert_all_lambda,'b-','LineWidth',line_size)
xlabel('\lambda','FontSize',letter_size)
ylabel('{\it p}','FontSize',letter_size)
ylim([0 1])
ax          = gca;
ax.YColor   = 'b';
ax.FontSize = number_size; 

% Plot der Werte fuer den Ordnungsparameter an der rechten y-Achse
yyaxis right
plot(l_vec,r_strich_all_lambda,'r-','LineWidth',line_size)
ylabel('$\overline{R}_{\phi}$', 'Interpreter', 'latex','FontSize',letter_size)
xlim([0,1]); ylim([0 1])
ax          = gca;
ax.YColor   = 'r';
ax.FontSize = number_size; 
% legend('{\it p}','Ordnungsparameter')

% Speichern der Abbildung
saveas(gcf,'Ordp_l_M_4','svg')

% Plot der Differenz der Winkel zu dem Winkel der Orientierungspraeferenz
figure(17)
plot(l_vec,phi_diff_all_lambda,'k-','LineWidth',line_size)
xlabel('\lambda','FontSize',letter_size)
ylabel('$(\overline{\phi} - \psi)$ in $^{\circ}$','Interpreter', 'latex','FontSize',letter_size)
xlim([0,1]); ylim([0 360])
ax = gca; 
ax.FontSize = number_size; 

% Speichern der Abbildung
saveas(gcf,'Diff_l_M_4','svg')
%% Bestimmung des p-Wertes und Ordnungsparameter fuer variierendes lambda und 
%  eta

p_data_iter        = NaN(length(n_vec),length(l_vec),iter_end);
r_data_iter        = NaN(length(n_vec),length(l_vec),iter_end);
phi_diff_data_iter = NaN(length(n_vec),length(l_vec),iter_end);
p_data_n_l         = NaN(length(n_vec),length(l_vec));
ord_data_n_l       = NaN(length(n_vec),length(l_vec));
phi_diff_data_n_l  = NaN(length(n_vec),length(l_vec));


for i_n = 1 : length(n_vec)
    n = n_vec(i_n);

    for i_l = 1 : length(l_vec)
        lambda = [l_vec(i_l),1,1,1];

        for i_iter = 1:iter_end
            [x_vec,y_vec,phi_vec,theta_vec,r_vec] = modell_schwarm(N,time_sim,dt,r,L,n,v,lambda,psi(N));
            phi_vec_row = reshape(phi_vec,numel(phi_vec),1);
            p_data_iter(i_n,i_l,i_iter)        = circ_rtest(phi_vec_row);
            r_data_iter(i_n,i_l,i_iter)        = circ_r(phi_vec_row);
            phi_diff_data_iter(i_n,i_l,i_iter) = abs(circ_mean(circ_dist(psi_vec(1),phi_vec_row)));
        end

        p_data_n_l(i_n,i_l)        = mean(p_data_iter(i_n,i_l,:));
        ord_data_n_l(i_n,i_l)      = mean(r_data_iter(i_n,i_l,:));
        phi_diff_data_n_l(i_n,i_l) = circ_rad2ang(circ_mean(reshape(phi_diff_data_iter(i_n,i_l,:),iter_end,1)));
    end
end
       
%% Plot der Ergebnisse fuer variierende Werte von lambda und eta
% Plot der p-Werte 
figure(18)
pl = pcolor(l_vec,n_vec,p_data_n_l);
colormap jet
pl.FaceColor     = 'interp';
c                = colorbar; 
c.Title.String   = 'p';
c.Title.FontSize = letter_size; 
xlabel('\lambda','FontSize',letter_size)
ylabel('\eta','FontSize',letter_size)
yticks([0,1/3*pi,2/3*pi,pi,4/3*pi,5/3*pi,2*pi]);
yticklabels({'0','1/3 \pi','2/3 \pi','\pi','4/3 \pi','5/3 \pi','2 \pi'});

% Speichern der Abbildung
saveas(gcf,'p_l_n_M_4','svg')

% Plot der Werte fuer den Ordnungsparameter
figure(19)
colormap jet
oldcmap             = colormap;
colormap(flipud(oldcmap))
pl                  = pcolor(l_vec,n_vec,ord_data_n_l);
pl.FaceColor        = 'interp';
c                   = colorbar; 
c.Title.String      = '$\overline{R}_{\phi}$';
c.Title.Interpreter = 'latex';
c.Title.FontSize    = letter_size;
xlabel('\lambda','FontSize',letter_size)
ylabel('\eta','FontSize',letter_size)
yticks([0,1/3*pi,2/3*pi,pi,4/3*pi,5/3*pi,2*pi]);
yticklabels({'0','1/3 \pi','2/3 \pi','\pi','4/3 \pi','5/3 \pi','2 \pi'});

% Speichern der Abbildung
saveas(gcf,'Ord_l_n_M_4','svg')

% Plot der Differenzen
figure(20)
colormap jet
pl                  = pcolor(l_vec,n_vec,phi_diff_data_n_l);
pl.FaceColor        = 'interp';
c                   = colorbar; 
c.Title.String      = '$(\overline{\phi} - \psi)$ in $^{\circ}$';
c.Title.Interpreter = 'latex';
c.Title.FontSize    = letter_size;
xlabel('\lambda','FontSize',letter_size)
ylabel('\eta','FontSize',letter_size)
yticks([0,1/3*pi,2/3*pi,pi,4/3*pi,5/3*pi,2*pi]);
yticklabels({'0','1/3 \pi','2/3 \pi','\pi','4/3 \pi','5/3 \pi','2 \pi'});

% Speichern der Abbildung
saveas(gcf,'Diff_l_n_M_4','svg')