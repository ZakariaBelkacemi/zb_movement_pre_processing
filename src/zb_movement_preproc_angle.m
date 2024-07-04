

function[movresult] = zb_movement_preproc_angle(tsvfilename,taskName,varargin)

%crop synchro application

%output : asked marker positions croped

ignoredScan = 1;%defaultscanNumber = 229;%default

repetitionTime = 0.8;%default
blockDurationByScan = 38;%default
blockDurationBySec = blockDurationByScan * repetitionTime;
movSamplFreq = 100;%default
blockDurationByMovsampl = blockDurationBySec * movSamplFreq;
keepOldMovSegment = 0;
%markers of interest, depends on the task
marker1 = 'shoulderR';
marker2 = 'elbowR';
marker3 = 'handR';
markerOfInterest = 'handR';
scanNumber = 229;

for i = 1:2:numel(varargin)
    switch(varargin{i})
        case 'markerOfInterest'
            markerOfInterest = varargin{i+1};
        case 'otherMarker'
            otherMarker = varargin{i+1};
        case 'SamplFreq'
            movSamplFreq = varargin{i+1};
        case 'onsets_time'
            onsets_time = varargin{i+1};
        case 'repetition_time'
            repetitionTime = varargin{i+1};
        case 'scan_number'
            scanNumber = varargin{i+1};
        case 'block_duration'
            blockDurationByScan = varargin{i+1};
        case 'ignoredScan'
            ignoredScan = varargin{i+1};
        case 'keepOldMovSegment' % 0 or 1
            keepOldMovSegment = varargin{i+1};
        case 'oldMovOnsets' % give the upsampled one
            oldMovOnsets = varargin{i+1};
        case 'oldMovDurations' % give the upsampled one
            oldMovDurations = varargin{i+1};          
        otherwise
    end
end


%%
%===============================================================
%VARIABLES INITIALIZATION
%===============================================================

MRISamplFreq = 1/repetitionTime;

disp('extracting QTM.tsv file');
%The zb_qtm_tsvread function assumes that the header takes 12 lines.
%The new synchronisation script places an event in QTM at the first and
%the second s. Having 2 event make the header taking 13 lines.
try
    [qtmdata,qtmheader]=zb_qtm_tsvread(tsvfilename); %take ~20sec in FRAICHE
    tsv_header_length = 12;
catch
    [qtmdata,qtmheader]=zb_qtm_tsvread(tsvfilename,'header_length',13); %take ~20sec in FRAICHE
    tsv_header_length = 13;
    %si erreur, peut-être que le tsvfilename est vide et donc pb dans
    %l'appel de la fonction
end

%Extracting movement sampling frequency
movSamplFreq = qtmheader{4};
movSamplFreq = movSamplFreq{1};
movSamplFreq = str2num(movSamplFreq(2));
disp(['Movement sampling frequency = ' num2str(movSamplFreq) 'Hz']);

%Extracting frame of the 2nd acquired fMRI vol
if tsv_header_length == 12
    MRIbegin = qtmheader{10};
elseif tsv_header_length == 13 % if the header contains 2 event, the second
                               %corresponds to the second event we want,
                               %which is the 11th element of the header
    MRIbegin = qtmheader{11};
end
MRIbegin = MRIbegin{1};
MRIbegin = str2num(MRIbegin(3)); % = frame of the 2nd s !
%if we ignore the 9 first scans for the MRI stabilisation (238 scans),
%we have to crop the movement data not from the 2nd mri vol signal received
%but at the 9th
if ignoredScan == 9
    shift = MRIbegin + 7*(1/MRISamplFreq)*movSamplFreq;%+7 car MRI begin est déjà le 2e vol
elseif ignoredScan == 1
    shift = MRIbegin;
end
    
%by default, MRIbegin is the event placed at the 2nd volume
%and the fMRI Preprocessing delete the first calibration volume

%Extracting list of available markers
if tsv_header_length == 12
    markerlist = qtmheader(11);
    markerlist = markerlist{1};
elseif tsv_header_length == 13 % if the header contains 2 event, the second
                               %corresponds to the second event we want,
                               %which is the 11th element of the header
                               %decaying all the next elements
    markerlist = qtmheader{12};
end

markerlist = markerlist{:};

numberOfCameras = qtmheader(2);
numberOfCameras = numberOfCameras{1};
numberOfCameras = numberOfCameras{1};
movresult.movementjson.numberOfCameras = numberOfCameras{2};

numberOfMarker = qtmheader(3);
numberOfMarker = numberOfMarker{1};
numberOfMarker = numberOfMarker{1};
movresult.movementjson.numberOfMarkers = numberOfMarker{2};

motionCaptureTime = qtmheader(8);
motionCaptureTime = motionCaptureTime{1};
motionCaptureTime = motionCaptureTime{1};
movresult.movementjson.timestamp = motionCaptureTime{2};

movresult.movementjson.movSamplFreq = movSamplFreq;
movresult.movementjson.markerlist = markerlist;

%%
%angle between marker1, marker2 and marker3

%in qtmdata : one column per component x, y and z following the order of markers
marker_idx = find(markerlist == convertCharsToStrings(markerOfInterest))-1;%-1 because the first element is the title of the line
x = qtmdata(:,3*(marker_idx-1)+1);
y = qtmdata(:,3*(marker_idx-1)+2);
z = qtmdata(:,3*(marker_idx-1)+3);

marker_idx = find(markerlist == convertCharsToStrings(marker1))-1;%-1 because the first element is the title of the line
x1 = qtmdata(:,3*(marker_idx-1)+1);
y1 = qtmdata(:,3*(marker_idx-1)+2);
z1 = qtmdata(:,3*(marker_idx-1)+3);

marker_idx = find(markerlist == convertCharsToStrings(marker2))-1;%-1 because the first element is the title of the line
x2 = qtmdata(:,3*(marker_idx-1)+1);
y2 = qtmdata(:,3*(marker_idx-1)+2);
z2 = qtmdata(:,3*(marker_idx-1)+3);

marker_idx = find(markerlist == convertCharsToStrings(marker3))-1;%-1 because the first element is the title of the line
x3 = qtmdata(:,3*(marker_idx-1)+1);
y3 = qtmdata(:,3*(marker_idx-1)+2);
z3 = qtmdata(:,3*(marker_idx-1)+3);

A_point = [x1,y1,z1];
B_point = [x2,y2,z2];
C_point = [x3,y3,z3];

AB_segment = B_point - A_point;
CB_segment = B_point - C_point;

angles = zeros(size(A_point,1),1);
for iframe = 1:size(A_point,1)
    ABdotCB = dot(AB_segment(iframe,:),CB_segment(iframe,:));%dot product of AB and CB
    norm_AB = norm(AB_segment(iframe,:));
    norm_CB = norm(CB_segment(iframe,:));
    cos_theta = ABdotCB / (norm_AB * norm_CB);
    theta = acos(cos_theta);
    angles(iframe) = rad2deg(theta);
end

%%
%===============================================================
%MOVEMENT SYNCHRONIZATION
%importing data and synchronizing the movement vector with the MRI begin
%===============================================================

%si erreur de shift undefined, check si ignorescan == 1 ou 9
movresult.angular_amplitude = angles(shift:shift+(scanNumber*repetitionTime*movSamplFreq));%in degree

movresult.x_position = x(shift:shift+(scanNumber*repetitionTime*movSamplFreq));%in mm
movresult.x_position = movresult.x_position/1000;%in meter, international system unit
movresult.y_position = y(shift:shift+(scanNumber*repetitionTime*movSamplFreq));%in mm
movresult.y_position = movresult.y_position/1000;%in meter, international system unit
movresult.z_position = z(shift:shift+(scanNumber*repetitionTime*movSamplFreq));%in mm
movresult.z_position = movresult.z_position/1000;%in meter, international system unitmovresult.position = sqrt(x.^2+y.^2+z.^2);%in meter
movresult.position = sqrt(movresult.x_position.^2+movresult.y_position.^2+movresult.z_position.^2);%in meter

%%
%===============================================================
%MOVEMENT CLEANING
%gap filling eventual lack of data and filtering outliers
%===============================================================

%to keep in memory brut vectors :
movresult.brut_angular_amplitude = movresult.angular_amplitude;

%low-pass filter to avoid velocity peaks due to masked-unmasked markers
movresult.angular_amplitude = zb_butterworth(movresult.angular_amplitude,10,movSamplFreq);

%to display original vector before the filter
filterview = figure('units','normalized','outerposition',[0 0 1 1]);
plot(movresult.brut_angular_amplitude);
%displaying and saving butterworth filter result
hold on;
plot(movresult.angular_amplitude);
pause(3);
F = getframe(gcf);
fig_im = frame2im(F);
movresult.movfiltered = frame2im(F);
pause(1);
close(filterview);

%%
%===============================================================
%PARAMETERS CALCULATION
%===============================================================

%1) Velocity and Acceleration
%--------

%calculating tangential velocity
dangle = gradient(movresult.angular_amplitude);%in degree
dt = gradient((1:numel(movresult.angular_amplitude))');%in cs
movresult.angular_velocity = dangle./dt;%in degree/centisecond
movresult.angular_velocity = movresult.angular_velocity*100;%in degree/sec

%calculating acceleration
dangular_velocity = gradient(dangle);
movresult.angular_acceleration = dangular_velocity./(dt.^2);%in degree/cs^-2
movresult.angular_acceleration = movresult.angular_acceleration * (100^2);%in degree/s^-2

%NOTE : possibilité de ne garder que les vitesses angulaires positives pour avoir
%les descentes et les vitesses angulaires négatives pour avoir les montées

%NOTE : beaucoup plus simple de segmenter le mouvement car on a l'angle
%correspondant à la position neutre (160° environ)

%%
%2) Moving average velocity
%--------

%smooth by using a moving average : 
movresult.averaged_angular_amplitude = zb_butterworth(abs(movresult.angular_amplitude),repetitionTime/4,movSamplFreq);
movresult.averaged_angular_amplitude(movresult.averaged_angular_amplitude<0)=0;
movresult.averaged_angular_velocity = zb_butterworth(abs(movresult.angular_velocity),repetitionTime/4,movSamplFreq);
movresult.averaged_angular_velocity(movresult.averaged_angular_velocity<0)=0;
movresult.averaged_angular_acceleration = zb_butterworth(abs(movresult.angular_acceleration),repetitionTime/4,movSamplFreq);
movresult.averaged_angular_acceleration(movresult.averaged_angular_acceleration<0)=0;

%NOTE : je pense que l'enveloppe du signal n'a d'intérêt que si l'on prend
%montée et descente, donc valeur aboslue du signal, car si l'on prend que
%les montées ça se smoothera un peu et je ne sais pas si ça sert d'avoir la
%variation moyenne que des montées sachant qu'on sait qu'au milieu il y a
%des descentes de coude

%%
%3) Undersampling of continuous parameters
%--------
%Undersampling from motion capture sampling frequency to fMRI sampling frequency
%interpolation method : https://www.mathworks.com/help/matlab/ref/interp1.html#btwp6lt-1-method
movresult.interpolated_angular_amplitude = interp1([1:numel(movresult.angular_amplitude)],movresult.angular_amplitude,[1:numel(movresult.angular_amplitude)/scanNumber:numel(movresult.angular_amplitude)],'pchip');
movresult.interpolated_averaged_angular_amplitude = interp1([1:numel(movresult.averaged_angular_amplitude)],movresult.averaged_angular_amplitude,[1:numel(movresult.averaged_angular_amplitude)/scanNumber:numel(movresult.averaged_angular_amplitude)],'pchip');
movresult.interpolated_angular_velocity = interp1([1:numel(movresult.angular_velocity)],movresult.angular_velocity,[1:numel(movresult.angular_velocity)/scanNumber:numel(movresult.angular_velocity)],'pchip');
movresult.interpolated_averaged_angular_velocity = interp1([1:numel(movresult.averaged_angular_velocity)],movresult.averaged_angular_velocity,[1:numel(movresult.averaged_angular_velocity)/scanNumber:numel(movresult.averaged_angular_velocity)],'pchip');
movresult.interpolated_angular_acceleration = interp1([1:numel(movresult.angular_acceleration)],movresult.angular_acceleration,[1:numel(movresult.angular_acceleration)/scanNumber:numel(movresult.angular_acceleration)],'pchip');
movresult.interpolated_averaged_angular_acceleration = interp1([1:numel(movresult.averaged_angular_acceleration)],movresult.averaged_angular_acceleration,[1:numel(movresult.averaged_angular_acceleration)/scanNumber:numel(movresult.averaged_angular_acceleration)],'pchip');
movresult.interpolated_position = interp1([1:numel(movresult.position)],movresult.position,[1:numel(movresult.position)/scanNumber:numel(movresult.position)],'pchip');

%% 
% %4) Part with a code i took from outside and that i do not completely understand
% %Localizing pos and neg peak which are during action instruction
% %--------
% 
% movbloc = blockDurationByScan*repetitionTime*movSamplFreq;%duration of a bloc
% TmovOverfilt = [[movbloc:movbloc*2],[movbloc*3:movbloc*4],[movbloc*5:movbloc*6]];%Times index of action instruction
% %time in seconds or in 100Hz ?
% 
% 
% scaledPosVector = movresult.angular_amplitude - min(movresult.angular_amplitude);%scaling near to 0
% 
% % Ts = 1/movSamplFreq;
% 
% movDuration = numel(scaledPosVector)/100;
% %pourquoi /100 ?

%%
%5) Part with a code i took from outside and that i do not completely understand
% Velocity profile/ lifting movement
% max speed, time to peak,%longeuer trajectoire couvert
%--------
       
% [movresult.amplitudeMean, movresult.vMean, movresult.vMax, movresult.TTP, movresult.TLLift, movresult.DistToPeak, movresult.PercTrajPeak]=VelocityProfile(fivEP,T,movPosP,movNegP);



%%
%===============================================================
%VISUALIZATION
%plotting and saving the movement with some parameters
%===============================================================

onsetsNumber = 3*2; %3 blocks with 2 phases : action and rest
%in e-prime, we add an additional scan for an eventual music delay, if
%there is music in the protocol
% if contains(protocol,'music')
%     scanNumber = scanNumber-1;
% end

if exist('blockDurationByScan')%varargin not necessary defined
    onsetsDuration = blockDurationByScan;
else
    onsetsDuration = scanNumber/onsetsNumber; %unit : scan (repetition time)
end

onsets = zeros(scanNumber,1);
onsets = zb_onsetSimulation(onsetsNumber/2,onsetsDuration);
mov_visu = figure('units','normalized','outerposition',[0 0 1 1]);
plot(upfirdn(onsets,1,80)*max(abs(movresult.angular_velocity)),'y');
% plot(onsets*max(movresult.interpolate_croped_vel),'b');
hold on
plot(abs(movresult.angular_velocity),'r');
% plot(movresult.croped_envelope);
plot(movresult.averaged_angular_velocity,'black');
title('velocity in degree.s^-1 (absolute value)');
pause(3);
F = getframe(gcf);
movresult.block_visu = frame2im(F);
close(mov_visu);

%%
%===============================================================
%CALCULATION OF OTHER ONSETS
%===============================================================

%%
%1) real_onsets
%corresponds to onsets of the movement effectively performed by the subject
%--------

absolute_vel = abs(movresult.interpolated_angular_velocity);
absolute_acc = abs(movresult.interpolated_angular_acceleration);

movresult.real_onsets_vector = zeros(numel(absolute_vel),1);
baseline = mean(absolute_vel);
for iframe = 1:numel(absolute_vel)
    if absolute_vel(iframe)>baseline
         movresult.real_onsets_vector(iframe) = 1;
    end
end

nb_zeros = 0; % Counter pour consecutiv zeros

iframe = 0;
for iframe = 1:numel(movresult.real_onsets_vector)
    if movresult.real_onsets_vector(iframe) == 0
        nb_zeros = nb_zeros + 1;
    elseif nb_zeros >= 10
        nb_zeros = 0;
    elseif nb_zeros >= 1 && movresult.real_onsets_vector(iframe) == 1
        movresult.real_onsets_vector(iframe-nb_zeros:iframe-1) = 1; % Modification of zeros in ones
        nb_zeros = 0; % Reinitialization of the zeros counter
    else
        nb_zeros = 0;
    end    
end

% real onsets time and duration :
movresult.real_onsets = [];
movresult.real_onsets_durations = [];
lock_onsets = 0;
onset_duration = 0;

for iframe = 1:numel(movresult.real_onsets_vector)
    if movresult.real_onsets_vector(iframe) == 1 && lock_onsets == 0
        movresult.real_onsets(end+1) = iframe;
        lock_onsets = 1;
        onset_duration = onset_duration + 1;
    elseif movresult.real_onsets_vector(iframe) == 1 && lock_onsets == 1
        onset_duration = onset_duration + 1;
    elseif movresult.real_onsets_vector(iframe) == 0
        if lock_onsets == 1
            movresult.real_onsets_durations(end+1) = onset_duration;
            onset_duration = 0;
        end
        lock_onsets = 0;
    end
end

% Verify if an onset is occuring at the end of the vector
if lock_onsets == 1
    movresult.real_onsets_durations(end+1) = onset_duration;
end
        
%%
%2) initiation onset
% onsets for the time between the displaying of the instruction and the
% effective movement
%--------

%initiation : 
%onset time = onset time instruction
%onset duration = onset time real onset - onset time instruction
onsets_time = [blockDurationByScan blockDurationByScan*3 blockDurationByScan*5];
movresult.initiation_onsets_durations = movresult.real_onsets - onsets_time;
disp(['instruction onsets : ' num2str(onsets_time)]);
disp(['movement onsets : ' num2str(movresult.real_onsets)]);
disp(['difference : ' num2str(movresult.initiation_onsets_durations)]);

%transforming negative values to zeros
if movresult.initiation_onsets_durations(1) < 0
    movresult.initiation_onsets_durations(1) = 0;
end
if movresult.initiation_onsets_durations(2) < 0
    movresult.initiation_onsets_durations(2) = 0;
end
if movresult.initiation_onsets_durations(3) < 0
    movresult.initiation_onsets_durations(3) = 0;
end

%%
%3) Onsets vector creation
%converting onsets time and duration to vectors with zeros and ones
%--------

movresult.initiation_vector = zeros(1,scanNumber);
for i=1:numel(onsets_time)
    movresult.initiation_vector(onsets_time(i):onsets_time(i)+movresult.initiation_onsets_durations(i))=1;
end

movresult.onsets_vector = zeros(1,scanNumber);
for i=1:numel(onsets_time)
    movresult.onsets_vector(onsets_time(i):onsets_time(i)+onsetsDuration)=1;
end

%%
%4) Movement segementation

% scaledPosVector = movresult.croped_pos - min(movresult.croped_pos);%scaling near to 0
upsampled_onsets = upfirdn(movresult.real_onsets_vector,1,80);
upsampled_onsets = zb_butterworth(upsampled_onsets,0.2,100);
upsampled_onsets(upsampled_onsets<=0)=0;
upsampled_onsets(upsampled_onsets>mean(upsampled_onsets))=1;
upsampled_onsets(end:numel(movresult.angular_amplitude))=1;

%%
%finding peaks and movement detection : third method
% [PosPValue,PosP] = findpeaks(movresult.croped_pos,'MinPeakDistance',50,'MinPeakHeight',1.025*mean(movresult.croped_pos));
% % [PosPValue,PosP] = findpeaks(abs(movresult.angular_amplitude),'MinPeakDistance',150,'MinPeakHeight',1.07*mean(abs(movresult.angular_amplitude)));
% % %grâce au seuil, il s'avère non nécessaire de vérifier que les pics sont
% % %dans une période de mouvement
% % movPosP = PosP;
% % movresult.peakNumber = numel(movPosP);
% % movresult.movPosP = movPosP;


if keepOldMovSegment == 0
%     [movresult.each_mov_onsets_upsample,movresult.each_mov_duration_upsample] = zb_elbow_flexion_detection(movresult.croped_pos,movresult.croped_vel);
%     [movresult.each_mov_onsets_upsample,movresult.each_mov_duration_upsample] = zb_movement_detection(movresult.croped_pos,1.025*mean(movresult.croped_pos));
%     [movresult.each_mov_onsets_upsample,movresult.each_mov_duration_upsample] = zb_hand_opening_detection(movresult.croped_pos,movresult.croped_vel);
%     [movresult.each_mov_onsets_upsample,movresult.each_mov_duration_upsample] = zb_elbow_flexion_detection2(movresult.angular_amplitude);
%     [movresult.each_mov_onsets_upsample,movresult.each_mov_duration_upsample] = zb_manual_mov_segmentation(movresult.angular_amplitude);
    [movresult.each_mov_onsets_upsample,movresult.each_mov_duration_upsample] = zb_elbow_angle_flexion_detection(movresult.angular_amplitude);
else
    movresult.each_mov_onsets_upsample = oldMovOnsets;
    movresult.each_mov_duration_upsample = oldMovDurations;
end

movresult.each_mov_onsets_upsample_second = movresult.each_mov_onsets_upsample/100;%in second
movresult.each_mov_duration_upsample_second = movresult.each_mov_duration_upsample/100;%in second

%affichage
peak_visu = figure('units','normalized','outerposition',[0 0 1 1]);
plot(movresult.angular_amplitude)
hold on
tvector = 1:1:numel(movresult.angular_amplitude);
% plot(tvector(movresult.movPosP),movresult.angular_amplitude(movresult.movPosP),'*r')
% plot([1 numel(movresult.angular_amplitude)],[1.005*mean(movresult.angular_amplitude) 1.005*mean(movresult.angular_amplitude)])
for imov = 1:numel(movresult.each_mov_onsets_upsample)
    plot([movresult.each_mov_onsets_upsample(imov) movresult.each_mov_onsets_upsample(imov)],[min(movresult.angular_amplitude) max(movresult.angular_amplitude)],'-','color','blue');
%     plot([movresult.each_mov_onsets_upsample(imov)+movresult.each_mov_duration_upsample(imov) movresult.each_mov_onsets_upsample(imov)+movresult.each_mov_duration_upsample(imov)],[min(movresult.angular_amplitude) max(movresult.angular_amplitude)],'-','color','red');
    p1 = [movresult.each_mov_onsets_upsample(imov) mean(movresult.angular_amplitude)];% First Point
    p2 = [movresult.each_mov_onsets_upsample(imov)+movresult.each_mov_duration_upsample(imov) mean(movresult.angular_amplitude)];% Second Point
    dp = p2-p1;% Difference
    quiver(p1(1),p1(2),dp(1),dp(2),0,'r')
end
F = getframe(gcf);
pause(1);
movresult.mov_segmentation = frame2im(F);
pause(2);
close(peak_visu);
%Warning : not already downsampled, to permit a brut parameter calculaion

%downsampling produit en croix
%si deux mouvements correspondent au même point à la fréquence
%d'échantillonage de l'IRMf, prendre la moyenne des deux
movresult.each_mov_onsets_downsample = round((movresult.each_mov_onsets_upsample * scanNumber) / numel(movresult.angular_amplitude));
movresult.each_mov_duration_downsample = round((movresult.each_mov_duration_upsample * scanNumber) / numel(movresult.angular_amplitude));

%%

movresult.parametric_sparc = [];
movresult.parametric_velocity = [];
movresult.parametric_acceleration = [];
movresult.parametric_amplitude = [];
movresult.paramectric_nTL = [];%normalized Trajectory Length : accuracy
xpos = movresult.x_position; ypos = movresult.y_position; zpos = movresult.y_position;%for nTL calculation
for iMov = 1:numel(movresult.each_mov_onsets_upsample)
    movbeg = movresult.each_mov_onsets_upsample(iMov);
    movend = movresult.each_mov_onsets_upsample(iMov)+movresult.each_mov_duration_upsample(iMov);
    %anticiper le fait que la fin du dernier mouvement soit en dehors du
    %vecteur si on a cliqué trop loin
    if movend>numel(movresult.angular_velocity)
        movresult.parametric_sparc(iMov)=SpectralArcLength(movresult.angular_velocity(movbeg:end),1/movSamplFreq);%sparc de chaque mouvement, de la partie du vecteur de vitesse de début de mouvement à cette partie + durée du mouvement
        movresult.parametric_velocity(iMov)=mean(abs(movresult.angular_velocity(movbeg:end)))/movresult.each_mov_duration_upsample(iMov);
        movresult.parametric_acceleration(iMov)=mean(abs(movresult.angular_acceleration(movbeg:end)))/movresult.each_mov_duration_upsample(iMov);
        %amplitude angulaire : max-min pour chaque mouvement
        movresult.parametric_amplitude(iMov)=abs(max(movresult.angular_amplitude(movbeg:end))-min(movresult.angular_amplitude(movbeg:end)));
        %accuracy : normalized Trajectory Length, chemin vs baseline =
        %amplitude en position
        ideal_trajectory_length = max(movresult.position(movbeg:end))-min(movresult.position(movbeg:end));
        movresult.parametric_nTL(iMov) = sum(abs(diff(movresult.position(movbeg:end))))/ideal_trajectory_length;   
    else
        movresult.parametric_sparc(iMov)=SpectralArcLength(movresult.angular_velocity(movbeg:movend),1/movSamplFreq);%sparc de chaque mouvement, de la partie du vecteur de vitesse de début de mouvement à cette partie + durée du mouvement
        movresult.parametric_velocity(iMov)=mean(abs(movresult.angular_velocity(movbeg:movend)))/movresult.each_mov_duration_upsample(iMov);
        movresult.parametric_acceleration(iMov)=mean(abs(movresult.angular_acceleration(movbeg:movend)))/movresult.each_mov_duration_upsample(iMov);
        movresult.parametric_amplitude(iMov)=abs(max(movresult.angular_amplitude(movbeg:movend))-min(movresult.angular_amplitude(movbeg:movend)));
        %accuracy : normalized Trajectory Length, chemin vs baseline =
        %amplitude en position
        ideal_trajectory_length = max(movresult.position(movbeg:movend))-min(movresult.position(movbeg:movend));
        movresult.parametric_nTL(iMov) = sum(abs(diff(movresult.position(movbeg:movend))))/ideal_trajectory_length;   
    end
end

%efficiency
% movresult.mean_movement_time =
% mean(movresult.each_mov_duration_downsample); corrigé :
movresult.mean_movement_time = mean(movresult.each_mov_duration_upsample)/100;

%vecteur paramétrique par bloc
upsampled_onsets_time = onsets_time * 100/1.25;
upsampled_onsets_duration = 38 * 100/1.25;
velocity_block1 = mean(abs(movresult.angular_velocity(upsampled_onsets_time(1):upsampled_onsets_time(1)+upsampled_onsets_duration)));
velocity_block2 = mean(abs(movresult.angular_velocity(upsampled_onsets_time(2):upsampled_onsets_time(2)+upsampled_onsets_duration)));
velocity_block3 = mean(abs(movresult.angular_velocity(upsampled_onsets_time(3):upsampled_onsets_time(3)+upsampled_onsets_duration)));
movresult.velocity_block = zeros(scanNumber,1);
movresult.velocity_block(onsets_time(1):onsets_time(1)+blockDurationByScan) = velocity_block1;
movresult.velocity_block(onsets_time(2):onsets_time(2)+blockDurationByScan) = velocity_block2;
movresult.velocity_block(onsets_time(3):onsets_time(3)+blockDurationByScan) = velocity_block3;

%moyenne du sparc des mouvements entre un onset time et la fin de l'onset
%time
sparc_block1 = mean(movresult.parametric_sparc(movresult.each_mov_onsets_upsample>upsampled_onsets_time(1) & movresult.each_mov_onsets_upsample<upsampled_onsets_time(1)+upsampled_onsets_duration));
sparc_block2 = mean(movresult.parametric_sparc(movresult.each_mov_onsets_upsample>upsampled_onsets_time(2) & movresult.each_mov_onsets_upsample<upsampled_onsets_time(2)+upsampled_onsets_duration));
sparc_block3 = mean(movresult.parametric_sparc(movresult.each_mov_onsets_upsample>upsampled_onsets_time(3) & movresult.each_mov_onsets_upsample<upsampled_onsets_time(3)+upsampled_onsets_duration));
movresult.neg_sparc_block = zeros(scanNumber,1);
movresult.neg_sparc_block(onsets_time(1):onsets_time(1)+blockDurationByScan) = -sparc_block1;
movresult.neg_sparc_block(onsets_time(2):onsets_time(2)+blockDurationByScan) = -sparc_block2;
movresult.neg_sparc_block(onsets_time(3):onsets_time(3)+blockDurationByScan) = -sparc_block3;

%acceleration
acceleration_block1 = mean(movresult.parametric_acceleration(movresult.each_mov_onsets_upsample>upsampled_onsets_time(1) & movresult.each_mov_onsets_upsample<upsampled_onsets_time(1)+upsampled_onsets_duration));
acceleration_block2 = mean(movresult.parametric_acceleration(movresult.each_mov_onsets_upsample>upsampled_onsets_time(2) & movresult.each_mov_onsets_upsample<upsampled_onsets_time(2)+upsampled_onsets_duration));
acceleration_block3 = mean(movresult.parametric_acceleration(movresult.each_mov_onsets_upsample>upsampled_onsets_time(3) & movresult.each_mov_onsets_upsample<upsampled_onsets_time(3)+upsampled_onsets_duration));
movresult.acceleration_block = zeros(scanNumber,1);
movresult.acceleration_block(onsets_time(1):onsets_time(1)+blockDurationByScan) = acceleration_block1;
movresult.acceleration_block(onsets_time(2):onsets_time(2)+blockDurationByScan) = acceleration_block2;
movresult.acceleration_block(onsets_time(3):onsets_time(3)+blockDurationByScan) = acceleration_block3;

%nTL
nTL_block1 = mean(movresult.parametric_nTL(movresult.each_mov_onsets_upsample>upsampled_onsets_time(1) & movresult.each_mov_onsets_upsample<upsampled_onsets_time(1)+upsampled_onsets_duration));
nTL_block2 = mean(movresult.parametric_nTL(movresult.each_mov_onsets_upsample>upsampled_onsets_time(2) & movresult.each_mov_onsets_upsample<upsampled_onsets_time(2)+upsampled_onsets_duration));
nTL_block3 = mean(movresult.parametric_nTL(movresult.each_mov_onsets_upsample>upsampled_onsets_time(3) & movresult.each_mov_onsets_upsample<upsampled_onsets_time(3)+upsampled_onsets_duration));
movresult.nTL_block = zeros(scanNumber,1);
movresult.nTL_block(onsets_time(1):onsets_time(1)+blockDurationByScan) = nTL_block1;
movresult.nTL_block(onsets_time(2):onsets_time(2)+blockDurationByScan) = nTL_block2;
movresult.nTL_block(onsets_time(3):onsets_time(3)+blockDurationByScan) = nTL_block3;

%vecteur paramétrique crénelé
movresult.crenelated_velocity = zeros(scanNumber,1);
for i = 1:numel(movresult.each_mov_onsets_downsample)
    movresult.crenelated_velocity(movresult.each_mov_onsets_downsample(i):movresult.each_mov_onsets_downsample(i)+movresult.each_mov_duration_downsample(i))=movresult.parametric_velocity(i);
end
movresult.crenelated_velocity = movresult.crenelated_velocity(1:scanNumber);%si la segmentation du mouvement dépasse la taille du vecteur : ne garder que la partie dans la taille du nombre de volumes acquis
%raccorder les valeures nulles intra block
%si 0, prend précédente valeur différente de 0, mais va pas chercher plus
%loin en arrière que la durée d'un bloc/2 :
for i = 1:numel(movresult.crenelated_velocity)
    if movresult.crenelated_velocity(i) == 0
        for j = 1:round(blockDurationByScan/2)
            if movresult.crenelated_velocity(j) ~= 0
                movresult.crenelated_velocity(i) = movresult.crenelated_velocity(j);
            end
        end
    end
end

movresult.crenelated_neg_sparc = zeros(scanNumber,1);
for i = 1:numel(movresult.each_mov_onsets_downsample)
    movresult.crenelated_neg_sparc(movresult.each_mov_onsets_downsample(i):movresult.each_mov_onsets_downsample(i)+movresult.each_mov_duration_downsample(i))=movresult.parametric_sparc(i);
end
%raccorder les valeures nulles intra block
%si 0, prend précédente valeur différente de 0, mais va pas chercher plus
%loin en arrière que la durée d'un bloc/2 :
for i = 1:numel(movresult.crenelated_neg_sparc)
    if movresult.crenelated_neg_sparc(i) == 0
        for j = 1:round(blockDurationByScan/2)
            if movresult.crenelated_neg_sparc(j) ~= 0
                movresult.crenelated_neg_sparc(i) = movresult.crenelated_neg_sparc(j);
            end
        end
    end
end

movresult.crenelated_acceleration = zeros(scanNumber,1);
for i = 1:numel(movresult.each_mov_onsets_downsample)
    movresult.crenelated_acceleration(movresult.each_mov_onsets_downsample(i):movresult.each_mov_onsets_downsample(i)+movresult.each_mov_duration_downsample(i))=movresult.parametric_acceleration(i);
end
%raccorder les valeures nulles intra block
%si 0, prend précédente valeur différente de 0, mais va pas chercher plus
%loin en arrière que la durée d'un bloc/2 :
for i = 1:numel(movresult.crenelated_acceleration)
    if movresult.crenelated_acceleration(i) == 0
        for j = 1:round(blockDurationByScan/2)
            if movresult.crenelated_acceleration(j) ~= 0
                movresult.crenelated_acceleration(i) = movresult.crenelated_acceleration(j);
            end
        end
    end
end

movresult.crenelated_amplitude = zeros(scanNumber,1);
for i = 1:numel(movresult.each_mov_onsets_downsample)
    movresult.crenelated_amplitude(movresult.each_mov_onsets_downsample(i):movresult.each_mov_onsets_downsample(i)+movresult.each_mov_duration_downsample(i))=movresult.parametric_amplitude(i);
end
%raccorder les valeures nulles intra block
%si 0, prend précédente valeur différente de 0, mais va pas chercher plus
%loin en arrière que la durée d'un bloc/2 :
for i = 1:numel(movresult.crenelated_amplitude)
    if movresult.crenelated_amplitude(i) == 0
        for j = 1:round(blockDurationByScan/2)
            if movresult.crenelated_amplitude(j) ~= 0
                movresult.crenelated_amplitude(i) = movresult.crenelated_amplitude(j);
            end
        end
    end
end

movresult.crenelated_nTL = zeros(scanNumber,1);
for i = 1:numel(movresult.each_mov_onsets_downsample)
    movresult.crenelated_nTL(movresult.each_mov_onsets_downsample(i):movresult.each_mov_onsets_downsample(i)+movresult.each_mov_duration_downsample(i))=movresult.parametric_nTL(i);
end
%raccorder les valeures nulles intra block
%si 0, prend précédente valeur différente de 0, mais va pas chercher plus
%loin en arrière que la durée d'un bloc/2 :
for i = 1:numel(movresult.crenelated_nTL)
    if movresult.crenelated_nTL(i) == 0
        for j = 1:round(blockDurationByScan/2)
            if movresult.crenelated_nTL(j) ~= 0
                movresult.crenelated_nTL(i) = movresult.crenelated_nTL(j);
            end
        end
    end
end


%%
%sinusoid_block : une série temporelle simulant à l'aide d'une sinusoïde,
%un paradigme en bloc avec des variations ressemblant aux variations d'un
%paramètre cinématique dans les blocs d'action, afin de vérifier que c'est
%bien la variation du signal qui fait foncitonner l'analyse irmf et non
%juste l'aspect sinusoïdal du vecteur de vitesse.

movresult.sinusoid_block = 1:1:scanNumber;
movresult.sinusoid_block(:) = 0;
for i = 1:2:5
    movresult.sinusoid_block(blockDurationByScan*i:blockDurationByScan*(i+1))=1;
end

movresult.sinusoid_block([blockDurationByScan:blockDurationByScan*2-1])=1/20*sin(2*pi*0.3*[1:blockDurationByScan*2-blockDurationByScan] + 1)+1;
movresult.sinusoid_block([blockDurationByScan*3:blockDurationByScan*4-1])=1/20*sin(2*pi*0.3*[1:blockDurationByScan*2-blockDurationByScan] + 1)+1;
movresult.sinusoid_block([blockDurationByScan*5:blockDurationByScan*6-1])=1/20*sin(2*pi*0.3*[1:blockDurationByScan*2-blockDurationByScan] + 1)+1;

disp('done movement preprocessing')
%maintenant, downsample les each_mov onsets, duration et les parametrics
%penser peut-être à augmenter les valeurs, pour l'accélération et la
%vitesse, faire que la valeur max avoisine le 1.10^1

%et param sur moyenne par block, voir si pas déjà fait,
%le faire sur real onsets * 100/1.25, et parametric modulation du coup

%%
% movresult.response_delay_to_instruction = 'not already programmed';

% movname = categorical({'vMean','vMax','peakNumber/100','TimeToPeak','DistToPeak'});
% bar(movname,[movresult.vMean, movresult.vMax, movresult.peakNumber/100, movresult.TTP, movresult.DistToPeak]);
% pause(1);
% F = getframe(gcf);
% movresult.barplotim = frame2im(F);
% close all;


end

