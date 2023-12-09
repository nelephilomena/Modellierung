function [x_vec,y_vec,r_vec,phi_vec,v_vec,time_sim,time_duration,phi_kompass] = ...
         read_tab_group(filename)

% Die Funktion liest die Daten mit dem Namen filename ein.
% und erstellt einen Vektor mit den jeweiligen x- und y-Werten. Es werden 
% fuer jeden Beobachtungszeitschritt sowie Fische der jeweilige Radius r,
% Positionswinkel phi und die Schwimmrichtung theta bestimmt. 
% Zusaetzlich ist als Mass der Schwarmbildung der Ordnungsparameter va und 
% die Ausdehnung, also der jeweilige Abstand zum center of mass bestimmt worden
%
% Syntax: 
%         [x_vec,y_vec,r_vec,phi_vec,time_sim,time_duration,phi_kompass] = ...
%          read_tab_group(filename)
%
% Parameter: 
%           x_vec         Vektor mit den x-Werten der Positionen pro dt
%           y_vec         Vektor mit den y-Werten der Positionen pro dt
%           r_vec         Vektor mit den Radien der Positionen pro dt
%           phi_vec       Vektor mit den Winkeln der Positionen pro dt
%           v_vec         Vekotr mit der zue
%           time_sim      Vektor mit den Beobachtungszeitpunkten
%           time_duration Anzahl der Zeitschritte
%           phi_kompass   Matrix mit den Werten der Matrix phi_vec in
%                         Kompassrichtung in Grad umgerechnet
%
%
%           file_name    Dateinname
%
%
% Nele Schuff, 10-07-2023

cd 'Daten'\Fischgruppen
N = 4; 

% Einlesen der Daten
rohdaten = readmatrix(filename);

% Loeschen der ersten Zeile
rohdaten(1,:) = [];

% Mittelpunkt des Kreises
x_zentrum = rohdaten(1,3);
y_zentrum = rohdaten(1,4);

% Loeschen der ersten Zeile mit dem Zentrum
daten = rohdaten(2:end,:);

% Fehlerkontrolle der Daten (Anzahl der Beobachtungen laesst sich nicht
% durch N teilen, weswegen entweder zu viele oder zu wenige Messwerte
% vorhanden sind)
if mod(size(daten,1),N) ~= 0
    c = 1;
    while c+N < size(daten,1)
        % Es gibt mehr Daten pro Zeitschritt als N
        if daten(c,5) == daten(c+N,5)
            daten(c+N,:) = [];
        end
        % Es gibt weniger Daten pro Zeitschritt als N
        if daten(c,5) ~= daten(c+N-1,5)
            daten = [daten(1:c+N-2,:);NaN(1,5);daten(c+N-1:end,:)];
        end
      
        c  = c + N;
    end
end
cd ../../

% Bestimmen der Beobachtungszeitpunkten   
if isnan(size(daten,1)/N)               % Wenn Zeitschritte sich nicht durch 
                                        % die Anzahl teilen laesst, werden
                                        % die Daten nicht ausgewertet
    time_sim = 0;
else
    time_sim = NaN(size(daten,1)/N,1);
    for i = 1 :size(time_sim,1)
        time_sim(i) = daten(i*N,5);
    end
end 

% Anzahl der Zeitschritten 
time_duration = size(time_sim,1);

% Speichern der x- und y-Werte asugehend von einem Zentrum bei (0,0)
x_vec = NaN(N,time_duration); y_vec = NaN(N,time_duration);
c = 1;
for i_t = 1: time_duration
    x_vec(:,i_t) = (x_zentrum -(daten(c:c+N-1,3))');
    y_vec(:,i_t) = (y_zentrum -(daten(c:c+N-1,4))');
    c = c + N;
end 

% Berechnung der Winkel phi der Positionen im Kreis und Radien r der 
% jeweiligen Positionen 
r_vec = NaN(N,time_duration);
phi_vec = NaN(N,time_duration);
for i_t = 1 : time_duration
        for j = 1 : N
            r_vec(j,i_t) = sqrt(x_vec(j,i_t)^2+y_vec(j,i_t)^2);
            phi_vec(j,i_t) = atan2(y_vec(j,i_t),x_vec(j,i_t)); 
        end 
end 

% Bestimmmen der zurueck gelegten Strecke pro Zeiteinheit tau
v_vec = NaN(N,time_duration-1);
for j = 1: N
    for t = 1: time_duration-1
        v_vec(j,t+1) = sqrt((x_vec(j,t+1)-x_vec(j,t))^2+(y_vec(j,t+1)-y_vec(j,t))^2);
    end
end

% Umrechnung des Winkels phi in positive Gradzahlen von 0° bis 360°
phi_ang_vec = circ_rad2ang(phi_vec');
phi_ang_vec = mod(phi_ang_vec',360);

% Umrechnung des Winkels phi in Kompasszahl
phi_kompass = NaN(N,time_duration);
phi_kompass(phi_ang_vec<= 90)=phi_ang_vec(phi_ang_vec <= 90)+270;
phi_kompass(phi_ang_vec>90)=(phi_ang_vec(phi_ang_vec>90)-90);
end