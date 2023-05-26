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
    % Obtain GB pixels
    R = get_gbpixels(LM);

%     R = rescale(LM,0,1);
    
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

%% Functions

function R = get_gbpixels(LM)

LM_p = padding(LM);
R = zeros(size(LM));

for i = 2:length(LM_p)-1
    for j = 2:length(LM_p)-1

        N_ids = get_neighbors(i,j,LM_p);

        if N_ids > 0

            R(i-1,j-1) = 1;

        end

    end
end
end

function LM_p = padding(LM)

LM_p = [LM(end,:); LM; LM(1,:)];

column_l = [0; LM(:,end); 0];
column_r = [0; LM(:,1); 0];

LM_p = [column_l LM_p column_r];

end

function N_ids = get_neighbors(i,j,LM_p)

ids(1) = LM_p(i-1,j); %up
ids(2) = LM_p(i+1,j); %down
ids(3) = LM_p(i,j-1); %left
ids(4) = LM_p(i,j+1); %rigth

id = LM_p(i,j);

N_ids = sum(ids~=id);

end

