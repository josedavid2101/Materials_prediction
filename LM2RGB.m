clear
clc

p = genpath('D:\Github\second\SGT_2D\get_SGT\Data\LM_64_120');
addpath(p) 

%% Load Data

% load('D:\Github\second\SGT_2D\GG_Simulations\Simulations_results\LM_RS_1000_256_800steps_2\LM_RS_1000_256_1.mat')
for i=1:50

    LM_name = ['LM_64_120_',num2str(i),'.mat'];
    load(LM_name)
    
    %% Obtain LM image
    
    LM = LM_continue.Label_matrices{200, 1};
    LM = LM(:,:,1);
    % Image must be single precision.
    LM = single(LM);
    R = rescale(LM,0,1);
    
    RGB(:,:,1) = R;
%     RGB(:,:,2) = R;
%     RGB(:,:,3) = R;
    
    % Create png file
    fileName = ['LM_64_120_',num2str(i),'.png'];
    imwrite(RGB,fileName);
    
    % Display it.
%     figure    
%     imshow(RGB, 'InitialMagnification', 100)
    
end
