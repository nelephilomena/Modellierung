function [x_vec,y_vec,r_vec,phi_vec,v_vec,theta_vec,filenames] = read_tab_single_xy

% Die Funktion laedt die Daten fuer die jeweilige 
% Fischart aus der Excel Tabelle ein. Dabei enthaelt die Excel-Datei die x-
% und y-Werte fuer jeden Fisch zu jedem Aufnahmezeitpunkt. Diese x- und
% y-Werte werden jeweils abgespeichert und der Radius r sowie der
% Positionswinkel phi berechnet
% Die Winkel phi aus der Tabell mit allen Beobachtunen und den Winkeln phi
% stimmen ueber ein mit der Berechnung der Winkel phi aus den Daten der x-
% und y-Koordinaten. 
%
% Syntax: 
%        [x_vec,y_vec,r_vec,phi_vec,v_vec,theta_vec,filenames] = read_tab_single_xy
%
% Parameter:
%           art      Bestimmt, ob von der Fischart 'Dodos' oder von der
%                    Fischart 'Cyanos' die Daten eingelesen werden und die
%                    jeweiligen zwei Tabellen erstellt werden 
%
%           x_vec    Vektor mit den x-Werten der Positionen des jeweiligen
%                    Fisches
%           y_vec    Vektor mit den y-Werten der Positionen des jeweiligen
%                    Fisches
%           r_vec    Vektor mit den Radien R der Positionen des jeweiligen
%                    Fisches
%           phi_vec  Vektor mit den Winkel phi der Positionen des jeweiligen
%                    Fisches
%           v_vec    Vektor mit den Strecken zwischen den einzelnen Zeitschritten 
%                    des jeweiligen Fisches
%          theta_vec Vektor mit den Schwimmrichtungen des jeweiligen
%                    Fisches
%
%
% Nele, Schuff, 10-10-2023

% Einlesen um Parameter zu bestimmen
cd 'Daten'\Dodos\
filenames = dir('*Results_Dodo*');
cd ..\..\
  
% Leere Vektoren erstellen
r_vec = zeros(41,length(filenames));
x_vec = zeros(41,length(filenames));
y_vec = zeros(41,length(filenames));
phi_vec = zeros(41,length(filenames));
v_vec = zeros(40,length(filenames));
theta_vec = zeros(40,length(filenames));

cd 'Daten'\Dodos\

for i = 1 : length(filenames)
    file = filenames(i).name;
    % Einlesen der Daten
    rohdaten = readmatrix(file);
    
    % Loeschen der ersten Zeile
    rohdaten(1,:) = [];

    % Mittelpunkt des Kreises
    zentrum_x = rohdaten(1,3);
    zentrum_y = rohdaten(1,4);

    % Loeschen der ersten Zeile mit dem Zentrum
    daten = rohdaten(2:end,:);
    
    % Einlesen der jeweiligen x- und y- Werten
    x_vec(1:size(daten,1),i) = daten(:,3)-zentrum_x;
    y_vec(1:size(daten,1),i) = daten(:,4)-zentrum_y;

    % Berechnung des Radius r
    r_vec(1:size(x_vec,1),i) = sqrt(x_vec(:,i).^2 + y_vec(:,i).^2);
    % Berechnung des Positionswinkel phi
    phi_vec(1:size(x_vec,1),i) = atan2(y_vec(:,i),x_vec(:,i));
    % Berechnung der Geschwindigkeit v
    v_vec(1:size(x_vec,1)-1,i) = abs(r_vec(2:end,i)-r_vec(1:end-1,i));
    % Berechnung des Schwimmrichtungwinkels theta
    theta_vec(1:size(x_vec,1)-1,i) =  ...
        atan2(y_vec(2:end,i)-y_vec(1:end-1,i),x_vec(2:end,i)-x_vec(1:end-1,i));
end
cd ..\..\
end
