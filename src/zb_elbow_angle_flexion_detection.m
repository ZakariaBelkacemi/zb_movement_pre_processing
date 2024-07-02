

function[movement_onsets,movement_durations] = zb_elbow_angle_flexion_detection(angle_vector)

temps = 1:numel(angle_vector);
angles = angle_vector;

% Variables de sortie
debut = [];
duree = [];

% Seuil pour détecter le mouvement
seuil = 100;

% Variables de contrôle
en_mouvement = false;
debut_mouvement = 0;
dernier_pic = 0;
lock_pic = 0;
fin_mouvement_en_attente = 0;

% Parcourir la série temporelle
for i = 2:length(angles)
    % Détecter un pic local comme point de départ potentiel
    if lock_pic == 0 && angles(i) < angles(i-1)        
        if fin_mouvement_en_attente == 1
            fin_mouvement = temps(i);
            duree_mouvement = fin_mouvement - debut(end); % Calculer la durée à partir du dernier début détecté
            duree = [duree; duree_mouvement]; % Ajouter la durée du mouvement
            fin_mouvement_en_attente = 0;
        end
        
        dernier_pic = i - 1; % Enregistrer le pic avant la baisse
        lock_pic = 1;

%         disp(['dernier pic ' num2str(dernier_pic)])
    elseif angles(i) > angles(i-1)
        lock_pic = 0;
    end
    
    % Détecter le début du mouvement en dessous du seuil
    if ~en_mouvement && angles(i) < seuil
        % Début d'un mouvement confirmé
        en_mouvement = true;
        debut_mouvement = temps(dernier_pic);
        debut = [debut; debut_mouvement]; % Ajouter le début exact du mouvement
%         disp(['debut_mouvement ' num2str(debut_mouvement)])
    elseif en_mouvement && angles(i) >= seuil
        % Fin d'un mouvement
        en_mouvement = false;
        fin_mouvement = temps(i);
        fin_mouvement_en_attente = 1;
%         duree_mouvement = fin_mouvement - debut(end); % Calculer la durée à partir du dernier début détecté
%         duree = [duree; duree_mouvement]; % Ajouter la durée du mouvement
    end
end

movement_onsets = debut;
movement_durations = duree;

%le dernier mouvement peut ne pas ?tre termin?, ce qui fera un mouvement
%onsets de plus que les movement durations donc on harmonise :
movement_onsets = movement_onsets(1:numel(movement_durations));

end


