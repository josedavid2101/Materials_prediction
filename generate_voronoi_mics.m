clear
clc

p = genpath('D:\Github\second\SGT_2D\Reconstruction');
addpath(p) 

%% Load Data

ng = 20;
dims = [64 64 2];

for i=1:500
   
    %% Obtain LM image    
    LM = get_mic(ng,dims);
    
    % Image must be single precision.
    LM = single(LM);
    R = rescale(LM,0,1);
    
    RGB(:,:,1) = R;
    RGB(:,:,2) = R;
    RGB(:,:,3) = R;
    
    % Create png file
    fileName = ['LM_64_120_',num2str(i),'.png'];
    imwrite(RGB,fileName);
    
    % Display it.
%     figure    
%     imshow(RGB, 'InitialMagnification', 100)
    
end

function [LM] = get_mic(ng,dims)

grains1 = voronoidata3dcolumnar(ng,dims);  %create initial grains
[LM] = label_matrixGB(grains1,dims); LM = LM(:,:,1);

end