clear all
 clc

addpath(genpath('D:\Github\second\Large_files\426_2D_10000_800_simulations'));


%% Loop over simulations

for i=1:1
    i  
    %% Load data
    
    LMn = ['LM_10000_800_426_',num2str(i),'.mat'];
    load(LMn)

    %% Get descriptors

    no_grains = get_no_grains(LM_continue.Label_matrices); % Number of grains at each time step
%     tic
%     grain_size = get_grain_size(LM_continue.Label_matrices); % Average grain size at each time step
%     toc
    tic
    no_tjs = get_no_tjs(LM_continue.Label_matrices); % Number of triple junctions at each time step  
    toc
    
    %% Save data
    LM = ['LM_10000_800_426_',num2str(i),'.mat'];
    save(LM,'LM_continue','-v7.3')

end


%% FUNCTIONS


function [no_grains] = get_no_grains(LMs)
    
    no_grains = zeros(length(LMs),1);
    for i=1:length(LMs) %Loop over microstructures
        
        LM = LMs{i,1};
        ids = unique(LM);
        no_grains(i,1) = length(ids);

    end

end

function [grain_size] = get_grain_size(LMs)
    
    grain_size = zeros(length(LMs),1);
    for i=1:length(LMs) %Loop over microstructures
        
        LM = LMs{i,1};
        ids = unique(LM);
        sizes = zeros(length(ids),1);

        for j=1:length(ids)
            id = ids(j,1);
            find_ids = sum(sum(LM==id));
            sizes(j,1) = find_ids;
%             sizes(j,1) = length(find(LM==id));
        end

        grain_size(i,1) = mean(sizes);

    end

end

function no_tjs = get_no_tjs(LMs)
    no_tjs = zeros(length(LMs),1);
    for ii=1:length(LMs)

        LM = LMs{ii,1};
        no_tj = 0; %initialize the tj counter
        
        m = [LM LM(:,1)];
    
        lm_first_row = [LM(1,:) LM(1,1)];
    
        LM_ext = [m; lm_first_row];
    
        for i=1:length(LM_ext)-1 % loop over m
            for j=1:length(LM_ext)-1
                
                m = LM_ext(j:j+1,i:i+1); %Obtain 2x2 sample
                ids = unique(m); %get ids of the sample
                if length(ids)>2 % condition for a tj
                    no_tj = no_tj + 1;
                end
                
            end
        end
        no_tjs(ii,1) = no_tj;
    end
end