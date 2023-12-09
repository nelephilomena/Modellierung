function [x_vec,y_vec,phi_vec,theta_vec,r_vec] = ...
    modell_gruppe(N,time_sim,dt,r,L,n,v,lambda,psi,varargin)

% Die Funktion modell_gruppe() modelliert das Schwimmverhalten einer
% Fischgruppe (Schwarm). Das Modell basiert auf den Grundlagen von Couzin
% et al. 2002. Modifikationen ist das Hinzufuegen einer bevorzugten
% Richtung. Diese haengt wie bei Couzin et al. 2005 von einem
% weighting parameter bei mir s ab. Anders als bei Couzin et al. wird durch eine zufaellige Auswahl
% einer Von-Mises verteilung eine Stoerung um die bevorzugte Richtung hinzugefuegt.
%
% Syntax: 
%            [x_vec,y_vec,phi_vec,theta_vec,r_vec] = modell_gruppe(N,time_sim,dt,r,L,n,v,s,psi,varargin)
%
% Parameter:
%            N         Anzahl der Fische in der Gruppe
%            time_sim  Vektor mit der Zeitspanne der Simulation 
%            dt        Zeitabstand
%            r         Vektor mit drei Zahlen: 1. Zahl = Grenze zwischen 1. 
%                      und 2. Zone und 2. Zahl = Grenze zwischen 2. und 3.  
%                      Zone und 3. Zahl = Grenze zwischen 3. Zone und (4.) 
%                      keiner Zone
%            L         Groesse des Beckens (Radius)
%            n         Stoerung
%            v         Geschwindigkeit (Strecke pro Einheit)
%            s         Vektor mit fuer jedes N-Individuen eine Zahl zwischen 0
%                      und 1 ist, die das Schwarmverhalten charakterisiert (0 = kein
%                      Schwarmverhalten und 1 = Schwarmverhalten)
%            psi       Vektor mit der bevorzugten Richtung der N-Individiuen
%                      (nur entscheidend, wenn s<1 ist)
%
%            x_vec     x-Werte der simulierten Fische 
%            y_vec     y-Werte der simulierten Fische 
%            phi_vec   Winkel der Position der simulierten Fische 
%            theta_vec Winkel der Schwimmrichtung der simulierten Fische
%            r_vec     Radien der Positionen
% 
%
% Veraenderungen (10-10-23) von Z.112 - 137:
%                            Stoerung: circ_vmrnd(0,n,1) --> rand(1,1)*n
%                            beta: circ_vmrnd(p(j),k,1) --> p(j)
%
%
% Nele Schuff, 19-9-2023

% Definieren von leeren Vektoren
x_vec     = NaN(length(time_sim),N);
y_vec     = NaN(length(time_sim),N); 
phi_vec   = NaN(length(time_sim),N);
theta_vec = NaN(length(time_sim),N);

% Bestimmen der Anfangsbedingung fuer ein Gefaess mit einer Form
% eines Kreises
if nargin < 11
    % N zufaellige Auswahl der x- und y- Werte in dem Ring
    % Definieren der Winkel und Spanne fuer den Radius
    th = 0:pi/50:2*pi;
    radius_span = -L:0.1:L;
    % Erstellen von leeren Vektoren fuer x- und y-Werte, die in dem
    % Ring liegen
    x_ring = []; y_ring = [];
    % Bestimmen dieser x- und y-Werte in dem Ring
    for j = 1: length(radius_span)
        radius = radius_span(j);
        x_r=radius*cos(th);
        x_ring = [x_ring,x_r];
        y_r=radius*sin(th);
        y_ring = [y_ring,y_r];
    end
    % Bestimmung eines zufaelligen Indexes
    ind_ring = randi(length(x_ring),1,N);
    
    % Festlegen der Anfangswerte fuer die x-, y-, theta-, und phi- Werte                            
    x_vec(1,:) = x_ring(ind_ring); y_vec(1,:) = y_ring(ind_ring);
    theta_vec(1,:) = circ_vmrnd(0,0.0001,[1,N]);  
    phi_vec(1,:) = atan2(y_vec(1,:),x_vec(1,:));

% Wenn Daten vorhanden sind, sind die Anfangsbedingungen die gegebenen
% Daten
elseif nargin > 10
    x_vec(1,:) = varargin{1}(1,:); % x_vec_data
    y_vec(1,:) =varargin{2}(1,:);  % y_vec_data
    theta_vec(1,:) = varargin{3}(1,:); % theta_vec_data
    phi_vec(1,:) = atan2(y_vec(1,:),x_vec(1,:));
end

% Fuer jeden Zeitschritt dt ...
for t = 1: length(time_sim)-1

    % .. wird fuer jeden Fisch i ...
    for i = 1: N
        theta_sum = theta_vec(t,i);  % Vektor mit den aufzusummierenden Richtungswinkeln
        pos_x_y_r = [];              % Vektor mit den Positionswinkel der Fisch in der Zone r
        pos_x_y_a = [];              % Vektor mit den Positionswinkel der Fisch in der Zone a
        total_o = 0;                 % Anzahl der Fische in Zone der Orientierung
        total_a = 0;                 % Anzahl der Fische in Zone der Anziehung
        total_r = 0;                 % Anzahl der Fische in Zone der Abstossung ('Repulsion')
        w = [lambda(i);1-lambda(i)]; % Weighting parameter fuer Mittelwertbestimmung

        % ... die Abstaende zu den Nachbarn j berechnet und wenn es im
        % Radius einer der drei Zonen liegen, 
        for j = 1: N
            if j ~= i
                dist = sqrt((x_vec(t,j)-x_vec(t,i))^2+(y_vec(t,j) - y_vec(t,i))^2);

                % wird in der Zone z_r die Positionswinkel der Fische j
                % abgespeichert
                if dist <= r(1)
                    pos_x_y_r(end+1) = atan2(y_vec(t,j),x_vec(t,j)); 
                    total_r = total_r +1; 
                % wird in der Zone z_o die Richtungswinkel der Fische j
                % abgespeichert
                elseif dist > r(1) && dist <= r(2)
                    theta_sum(end+1) =  theta_vec(t,j);
                    total_o = total_o +1;
                % wird in der Zone z_a die Positionswinkel der Fische j
                % abgespeichert
                elseif dist > r(2) && dist <= r(3)
                    pos_x_y_a(end+1) = atan2(y_vec(t,j), x_vec(t,j));
                    total_a = total_a +1;
                end
            end
        end
        
        % Wenn die Nachbarfische j in der Zone der Abstossung vorhanden
        % sind, dann wird der Durchschnitt dieser Positionen bestimmt
        if total_r > 0 
            theta_rep = circ_mean(pos_x_y_r'); 
            theta_vec(t+1,i) = circ_mean([-theta_rep;psi(i)],w) + rand(1,1)*n;

        % Wenn die Nachbarfische j in der Zone der Anziehung und der Zone
        % der Orientierung liegen, dann wird jeweils der Durchschnitt der
        % Positionswinkel der Fische in der Zone der Anziehung und der
        % Durchschnitt der Fische in der Zone der Orientierung berechnet.
        % Von dem Durchschnitt der Positionswinkel und von dem Durchschnitt
        % der Richtungswinkel wird anschließend der Durchschnitt berechnet.
        elseif total_o > 0 && total_a > 0 
            theta_orient = circ_mean(theta_sum');
            theta_attrac = circ_mean(pos_x_y_a');
            theta_att_ori = circ_mean([theta_orient;theta_attrac]);
            theta_vec(t+1,i) = circ_mean([theta_att_ori;psi(i)],w) + rand(1,1)*n; 
        % Wenn nur Fische in der Zone der Anziehung sind, wird der
        % Durchschnitt der Positionswinkel bestimmt
        elseif total_o == 0 && total_a > 0 
            theta_attrac = circ_mean(pos_x_y_a');
            theta_vec(t+1,i) = circ_mean([theta_attrac; psi(i)],w) + rand(1,1)*n; 
        % Wenn nur Fische in der Zone der Orientierung sind, wird der
        % Durchschnitt der Richtungswinkel bestimmt
        elseif total_o > 0 && total_a == 0 
             theta_orient = circ_mean(theta_sum');
             theta_vec(t+1,i) = circ_mean([theta_orient; psi(i)],w) + rand(1,1)*n; 
        else 
        % Wenn sich keine Fische in einer der dreien Zonen befinden
        % oder der Fisch sich nicht an den Schwarmmitglieder orientiert
        % (lambda = 0, ergibt sich der Richtungswinkel fuer den naechsten
        % Zeitschritt durch den Durchschnitt des alten Richtungswinkel
        % und der Richtung der Orientierungspräferenz
            theta_vec(t+1,i) = circ_mean([theta_vec(t,i); psi(i)],w) + rand(1,1)*n;
        end

        % Aktualisieren der Positionen
        x_vec(t+1,i) = x_vec(t,i) + v*cos(theta_vec(t,i))*dt;
        y_vec(t+1,i) = y_vec(t,i) + v*sin(theta_vec(t,i))*dt;
        phi_vec(t+1,i) = atan2(y_vec(t,i),x_vec(t,i));

        % Randbedingung (wenn die neue Position ausserhalb des Beckens
        % liegt)
        radius = sqrt(x_vec(t+1,i)^2+(y_vec(t+1,i)^2));
        if radius >= L
            theta_vec(t,i) = 2*atan2(y_vec(t,i),x_vec(t,i))- theta_vec(t,i) + pi;
            x_vec(t+1,i) = x_vec(t,i) +v*cos(theta_vec(t,i))*dt;
            y_vec(t+1,i) = y_vec(t,i) +v*sin(theta_vec(t,i))*dt;
            theta_vec(t+1,i) = theta_vec(t,i);
        end
    end
end

% Berechnnug des Radius mithilfe des x- und y-Wertes
r_vec = sqrt(x_vec.^2+y_vec.^2);

end


