%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script to plot ODFs of a microstructure and calculate the Fourier
% coefficients to plot in the texture hull
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
clc

%% Add paths and Load data

q = genpath('D:\Github\second\SGT_2D\Functions');
addpath(q)

r = genpath('D:\Github\second\SGT_2D\Microstructure_hull\Fourier_basis');
addpath(r)

p = genpath('D:\Github\second\SGT_2D\Microstructure_hull\Fourier_basis\mtex-5.10.0');
addpath(p)

startup_mtex;

rng default

addpath(genpath('D:\Github\second\Large_files\426_2D_10000_800_simulations'));

for i=2:426 %Loop over simulations
    tic
    LMn = ['LM_10000_800_426_',num2str(i),'.mat'];
    load(LMn)
    ori_all = LM_continue.ori;
    times = LM_continue.timestosave;

    for jj=1:length(times)
        jj
        j = times(jj);
    
        %% Generate ODF
        ids = LM_continue.ids{j,1};
        ori = ori_all(ids,:);
        LM = LM_continue.Label_matrices{jj,1};
                
        %% Calculate
        
%         plot_odf(ori)

        f_coeff(:,jj) = get_odf_fourier(ori,LM,ids);

    end
    f_name = ['f_10000_800_426_',num2str(i),'.mat'];
    save(f_name,'f_coeff','-v7.3')
    clearvars f_coeff
    toc
end

%% Functions

function plot_odf(ori)

cs = crystalSymmetry('cubic');
rot = rotation(quaternion(ori(:,1),ori(:,2),ori(:,3),ori(:,4)));
ori = orientation(rot);

% the fixed crystal directions (100)
h = Miller({1,0,0},cs);

figure %Plot dots
plotPDF(ori,Miller(0,0,1,ori.CS),'antipodal','all');

figure %Plot ODF
plotPDF(ori,Miller(0,0,1,ori.CS),'contourf','all');

set(gca,'FontSize',20)

c = colorbar;
w = c.LineWidth;
c.LineWidth = 1.5;
c.FontSize = 30;
CLim(gcm,'equal');
% print(gcf,'BRK_Random_initial.png','-dpng','-r500');
% mtexColorbar
a = gca;
idx2kill = 1;
delete(a(idx2kill))



end

function [fhat_sym_Arbitrary2] = get_odf_fourier(ori,LM,ids)

rot = rotation(quaternion(ori(:,1),ori(:,2),ori(:,3),ori(:,4)));
N_grains = length(ori);
o_grains = orientation(rot);
V_grains = get_grains_area(LM,ids);

%% Generate Fundamental Orientations

% Define crystal and sample symmetry
CS = crystalSymmetry('Oh');
SS = specimenSymmetry('1');

%% Get Symmetrizing Coefficients

smax = 32;
[A_C,A_S] = symmetrizingCoefficientsWignerD(0:smax,CS,SS);

%% compute coefficients in symmetrized basis for each GRAIN orientation (takes about 1.19 seconds)

Mu = cellfun(@(X) size(X,1),A_C);
Nu = cellfun(@(X) size(X,2),A_S);
fhat_sym_Arbitrary2 = nan(sum(Mu.*Nu),N_grains);
for s = 0:smax
    workbar((s+1)/(smax+1),'Computing coefficients of symmetrized basis...')
    if Mu(s+1)*Nu(s+1) > 0 % only evaluate if there are symmetric basis functions for this order
        fhat_sym_tmp_Arbitrary2 = nan(Mu(s+1)*Nu(s+1),N_grains);
        for k = 1:N_grains
            % get D matrix
            D = sqrt((2*s)+1)*Wigner_D(s,o_grains(k)); % NOTE: we evaluate the basis function for each GRAIN orientation

            % multiply by symmetrizing coefficients
            D_sym = A_C{s+1}*D*A_S{s+1};

            fhat_sym_tmp_Arbitrary2(:,k) = conj(D_sym(:)); % note complex conjugate & normalization (the inner product has been defined as (1/(8*pi^2))*int(conj(f)*g) so we don't divide by 8*pi^2 here)
        end
        istart = sum(Mu(1:s).*Nu(1:s))+1;
        iend = sum(Mu(1:s+1).*Nu(1:s+1));
        fhat_sym_Arbitrary2(istart:iend,:) = fhat_sym_tmp_Arbitrary2;
    end
end

fhat_sym_Arbitrary2 = (fhat_sym_Arbitrary2*V_grains(:))./sum(V_grains); % NOTE: this is the weighted mean


end

function V_grains = get_grains_area(LM,ids)

    V_grains = zeros(length(ids),1);

    for j=1:length(ids)
        id = ids(j,1);
        find_ids = sum(sum(LM==id));
        V_grains(j,1) = find_ids;
    end

end