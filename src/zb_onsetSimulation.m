
function[vector] = zb_onsetSimulation(block_number,block_duration_scan)

% echant = 0.8;
% block_number = 3;
% block_duration_scan = 38;

size = block_number*block_duration_scan*2;
vector = ones(size,1);

for iscan = 0:2:block_number*2-1
    index1 = block_duration_scan*iscan;
    if index1 == 0
        index1 = 1;
    end
    index2 = block_duration_scan*(iscan+1);
    vector(index1:index2) = 0;
end

end

% PREVIOUS VERSION :
% function[vector] = zb_onsetSimulation(echant,size)
% vector = 1:1:size;
% vector(:) = 1;
% 
% for i = 0:2:5
%     index1 = (round(20/echant))*i;
%     if index1 == 0
%         index1 = 1;
%     end
%     index2 = (round(20/echant))*(i+1);
%     vector(index1:index2) = 0;
% end
% 
% end