

function[movement_onsets,movement_durations] = zb_elbow_angle_flexion_detection(angle_vector)

% angle_vector = movresult.croped_pos;
% velocity_vector = movresult.croped_vel;



[NegPValue,NegP] = findpeaks(-angle_vector,'MinPeakDistance',300,'MinPeakHeight',-mean(angle_vector));

% [VelPValue,VelP] = findpeaks(velocity_vector,'MinPeakDistance',50,'MinPeakHeight',1.025*mean(velocity_vector));

%plot(angle_vector)
%hold on
%plot(tvector(NegP),angle_vector(NegP),'*r')

%%

% % base_position = 0.95*mean(angle_vector(1:2500));%position de repos, environ 160°
%1:2500 = au milieu du bloc de repos de début (1:3000)
%facteur 0.9 pour rester un peu au dessous et détecter le seuil initial
% % max_position = 0.75*mean(angle_vector(3000:6000));%position de mouvement effectué
%3000:6000 = 1er bloc d'action

movement_onsets = [];
movement_durations = [];
ascending = 2;
for i=11:numel(angle_vector)%on commence à 11 pour évité d'être embêté par
    %le "angle_vector(i-2)" alors qu'on est en phase
    %de repos de toutes manières
%     disp(['movement_onsets = ' movement_onsets])
%     disp(['i = ' num2str(i)])
    if ascending == 1%si c'est une montée de signal
        if angle_vector(i-1)-angle_vector(i)>0.25 %si on détecte une descente
            movement_durations(end+1) = i-movement_onsets(end);
            movement_onsets(end+1) = i;
            ascending = 0;
        elseif angle_vector(i)-angle_vector(i-1)<0.25 %si on détecte un plat
            if angle_vector(i)-angle_vector(i-50)<0.25 %on vérifie que c'est vraiment plat
                ascending = 2;
                movement_durations(end+1) = i-movement_onsets(end);
            end
        end
    elseif ascending == 2%si c'est plat
        if angle_vector(i-1)-angle_vector(i)>0.25 %si on détecte une descente
            movement_onsets(end+1) = i;
            ascending = 0;
            disp('+++++++++++++++++++++++++++')
            disp(movement_onsets)
        end
    elseif ascending == 0%si c'est une descente de signal
        if angle_vector(i)-angle_vector(i-1)>0.25 %si on détecte une montée
            ascending = 1;
        end
    end
    
    %le dernier mouvement peut ne pas être terminé, ce qui fera un mouvement
    %onsets de plus que les movement durations donc on harmonise :
    movement_onsets = movement_onsets(1:numel(movement_durations));
    
    movement_ends = movement_onsets+movement_durations;
    
end

end


