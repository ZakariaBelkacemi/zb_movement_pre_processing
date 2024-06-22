
function[movement_onsets,movement_durations] = zb_manual_mov_segmentation(position_vector)

movement_onsets = [];
movement_ends = [];
movement_durations = [];

fig_mov_segment = figure;
plot(position_vector);
hold on;
zoom(fig_mov_segment,'on');

while true
    title('movement begining');
    %asking the user to place a point at the begining of the movement
    point = ginput(1);%gives a tuple that only the first element interest us
    point = point(1);%taking the first element which is the x position
    if point<0%by clicking on the left of the GUI, the while loop ends
        break
    end
    movement_onsets(end+1) = point;
    plot([point point],[min(position_vector) max(position_vector)],':','color','blue');
    
    title('movement end');
    %asking the user to place a point at the begining of the movement
    point = ginput(1);%gives a tuple that only the first element interest us
    point = point(1);%taking the first element which is the x position
    movement_ends(end+1) = point;
    movement_durations(end+1) = movement_ends(end) - movement_onsets(end);
    plot([point point],[min(position_vector) max(position_vector)],':','color','red');
end
close(fig_mov_segment)

end
    