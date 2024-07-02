

function[movement_onsets,movement_durations] = zb_elbow_angle_flexion_detection(angle_vector)

temps = 1:numel(angle_vector);
angles = angle_vector;

% Variables de sortie
debut = [];
duree = [];

% Seuil pour d�tecter le mouvement
seuil = 100;

% Variables de contr�le
en_mouvement = false;
debut_mouvement = 0;
dernier_pic = 0;
lock_pic = 0;
fin_mouvement_en_attente = 0;

% Parcourir la s�rie temporelle
for i = 2:length(angles)
    % D�tecter un pic local comme point de d�part potentiel
    if lock_pic == 0 && angles(i) < angles(i-1)        
        if fin_mouvement_en_attente == 1
            fin_mouvement = temps(i);
            duree_mouvement = fin_mouvement - debut(end); % Calculer la dur�e � partir du dernier d�but d�tect�
            duree = [duree; duree_mouvement]; % Ajouter la dur�e du mouvement
            fin_mouvement_en_attente = 0;
        end
        
        dernier_pic = i - 1; % Enregistrer le pic avant la baisse
        lock_pic = 1;

%         disp(['dernier pic ' num2str(dernier_pic)])
    elseif angles(i) > angles(i-1)
        lock_pic = 0;
    end
    
    % D�tecter le d�but du mouvement en dessous du seuil
    if ~en_mouvement && angles(i) < seuil
        % D�but d'un mouvement confirm�
        en_mouvement = true;
        debut_mouvement = temps(dernier_pic);
        debut = [debut; debut_mouvement]; % Ajouter le d�but exact du mouvement
%         disp(['debut_mouvement ' num2str(debut_mouvement)])
    elseif en_mouvement && angles(i) >= seuil
        % Fin d'un mouvement
        en_mouvement = false;
        fin_mouvement = temps(i);
        fin_mouvement_en_attente = 1;
%         duree_mouvement = fin_mouvement - debut(end); % Calculer la dur�e � partir du dernier d�but d�tect�
%         duree = [duree; duree_mouvement]; % Ajouter la dur�e du mouvement
    end
end

movement_onsets = debut;
movement_durations = duree;

%le dernier mouvement peut ne pas ?tre termin?, ce qui fera un mouvement
%onsets de plus que les movement durations donc on harmonise :
movement_onsets = movement_onsets(1:numel(movement_durations));

end


