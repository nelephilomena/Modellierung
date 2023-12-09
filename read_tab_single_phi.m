function [T_fisch_alle,T_fisch_mean] = read_tab_single_phi()

% Die Funktion laedt die Daten fuer die jeweilige 
% Fischart aus der Excel Tabelle ein und speichert den Mittelwert der
% Winkel pro Zeitschritt in der Tabelle ab und es werden die statistischen
% Groessen wie Varianz, Standardabweichung, durchschnittlicher Radius und
% Z-Wert und p-Wert des Rayleigh Tests bestimmt. Anschließend werden diese
% Werte fuer jeden Fisch ueber alle Beobachtungen berechnet und in einer
% weiteren Tabelle (T_fisch_mean) abgespeichert.
%
% Syntax: 
%        [T_fisch_alle,T_fisch_mean] = read_tab_single_phi()
%
% Parameter:
%               art  Bestimmt, ob von der Fischart 'Dodos' oder von der
%                    Fischart 'Cyanos' die Daten eingelesen werden und die
%                    jeweiligen zwei Tabellen erstellt werden 
%
%      T_fisch_alle  Tabelle mit allen Beobachtungen. Fuer jede Beobachtung 
%                    wird der durchschnittliche Winkel phi, der 
%                    durchschnittliche Radius, Varianz der Winkel phi
%                    Standardabweichung, Z-wert und p-Wert des Rayleigh
%                    Tests der Winkel phi pro Beobachtung
%
%     T_fisch_mean   Tabelle mit den durchschnittlichen Winkel phi,
%                    durchschnittlicher Radius, Varianz,
%                    Standardabweichung, Z-wert und p-Wert des Rayleigh
%                    Tests der Winkel phi ueber alle Beobachtungen pro Fisch
%
% Nele, Schuff, 25-08-2023

%% Einladen der Tabelle für die Art der ...
cd 'Daten'\Dodos\
T_data = importdata('Dodos_angles_all_Gabi.xlsx');
cd ..\..

% In Zeile  11 sind die Fischnummern fuer alle Beobachtungen (1 1 1 1 2 2 2) 
% aufgelistet
fischnummern_daten = T_data.data(11,:);   

% Bestimmen der Fischnummern, so dass sie nur einmal aufgelistet werden
fischnummern = unique(T_data.data(11,:));


% In Spalte 12 bis 51 sind alle Positionswinkel phi für jeden Fisch und 
% alle Beobachtungen aufgelistet
phi_rad_vec = circ_ang2rad(T_data.data(12:51,:)); 

% Anzahl aller Beobachtungen fuer alle Fischen 
n_beobachtungen = size(T_data.data,2);

%% Fuer jeden Fisch und fuer jede Beobachtung mit dem Fisch werden mit den 
%  in den Zeilen angegebenen phi-Werte statistische Groessen berechnet und
%  dann fuer den Fisch und der jeweiligen Beobachtung in der Tabelle 
%  abgespeichert

T_fisch_alle = table('Size',[n_beobachtungen,9],'VariableTypes',{'double',...
    'double','double','double','double','double','double','double','double'},...
    'VariableNames',{'Fischnummer','Beobachtung',...
    'Bewoelkung','Winkel phi','r_strich','Varianz','Std','Rayleigh','p'});

% Variablen fuer die Fischnummer und Beobachtungen, die sich in der 
% for-Schleife erhoehen sollen (Counter)
c_fischnummer = 1;           
c_beob = 0;                 

for i_b = 1: n_beobachtungen

    % Wenn es noch die gleiche Fischnummer ist, ...
    if fischnummern_daten(i_b) == fischnummern(c_fischnummer)
        % ... wird die Beobachtungsnummer um eins erhoeht
        c_beob = c_beob+1;                         
    else
        % und wenn nicht, dann wird die Fischnummer um eins erhoeht
        c_beob = 1;
        c_fischnummer = c_fischnummer+1;
    end

    % Abspeichern der folgenden Werte in der Tabelle

    % Fischnummer
    T_fisch_alle{i_b,1} = fischnummern_daten(i_b);                   
    % Aktuelle Beobachtungsnummer des jeweiligen Fisches
    T_fisch_alle{i_b,2} = c_beob;       
    % Bewoelkung
    T_fisch_alle{i_b,3} = T_data.data(9,i_b);  
    % Durchschnitt ueber die Winkel phi 
    T_fisch_alle{i_b,4} = circ_mean(phi_rad_vec(:,i_b));              
    % Ordnungsparameter fuer die Winkel phi
    T_fisch_alle{i_b,5} = circ_r(phi_rad_vec(:,i_b)); 
    % Varianz der Winkel phi
    T_fisch_alle{i_b,6} = 1 - circ_r(phi_rad_vec(:,i_b));   
    % Standardabweichung der Winkel phi in Grad
    T_fisch_alle{i_b,7} = sqrt(2*(1-circ_r(phi_rad_vec(:,i_b)))); 
    % Rayleigh Test Z-Wert der Winkel phi 
    [~,T_fisch_alle{i_b,8}] = circ_rtest(phi_rad_vec(:,i_b));         
    % Rayleigh Test p-Wert der Winkel phi
    [T_fisch_alle{i_b,9},~] = circ_rtest(phi_rad_vec(:,i_b));         
end

%% Zusammenfassung der Beobachtungen pro Fisch
% Es wird jeweils der Durchschnitt von allen Beobachtungen fuer einen Fisch
% berechnet. Mit dem Ziel, dass in dieser Tabelle fuer jedenen einzelnen
% Fisch die durchschnittlichen Werte ueber alle Beobachtungen vorhanden
% sind
T_fisch_mean = table('Size',[length(fischnummern),8],'VariableTypes',...
    {'double','double','double','double','double','double','double','double'...
    },'VariableNames',{'Fischnummer',...
    'Anzahl der Beobachtungen','phi_mean','r_strich_mean','var_mean','std_mean',...
    'Rayleigh_mean','p_mean'});

for i_f = 1 : length(fischnummern)

    % Bestimmen der Indexe fuer die Fischnummer = fischnummern(k)
    ind_f = find((T_fisch_alle{:,1}) == fischnummern(i_f));
    
    % Berechnung der Durchschnitte ueber alle Beobachtung fuer jeden Fisch
    % und abspeichern dieser Werte in einer Tabelle
    % Fischnummer
    T_fisch_mean{i_f,1} = T_fisch_alle{ind_f(1),1};  
    % Anzahl der Beobachtungen pro Fisch
    T_fisch_mean{i_f,2} = T_fisch_alle{ind_f(end),2};                    
    % Durchschnittlicher Winkel phi der durchschnittlichen Winkel pro 
    % Beobachtung  -> Ordnungspräferenz 
    T_fisch_mean{i_f,3} = circ_mean(T_fisch_alle{ind_f,4});  
    % Durchschnittliche Ordnungsparameter der durchschnittlichen Winkel phi   
    T_fisch_mean{i_f,4} = circ_r(T_fisch_alle{ind_f,4});
    % Durchschnittliche Varianz der durchschnittlichen Winkel pro
    % Beobachtung
    T_fisch_mean{i_f,5} = 1-circ_r(T_fisch_alle{ind_f,4});   
    % Durchschnittliche Standardabweichung der durchschnittlichen Winkel pro
    % Beobachtung
    T_fisch_mean{i_f,6} = circ_std(T_fisch_alle{ind_f,4});
    %  Rayleigh Test Z-Wert von den durchschnittlichen Winkel phi pro 
    % Beobachtung
    [~,T_fisch_mean{i_f,7}] = circ_rtest(T_fisch_alle{ind_f,4});          
    %  Rayleigh Test p-Wert von den durchschnittlichen Winkel phi pro 
    % Beobachtung
    [T_fisch_mean{i_f,8},~] = circ_rtest(T_fisch_alle{ind_f,4});          
end
end
