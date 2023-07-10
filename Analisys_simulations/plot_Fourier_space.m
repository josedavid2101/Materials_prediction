clear
clc


%% Add paths and Load data

% q = genpath('D:\Github\second\SGT_2D\Functions');
% addpath(q)

r = genpath('coefficients_426simulations');
addpath(r)

pp = genpath('D:\Github\second\SGT_2D\Microstructure_hull\Fourier_basis\mtex-5.10.0');
addpath(pp)

addpath(genpath('D:\Github\second\Large_files\426_2D_10000_800_simulations'));

startup_mtex;

load('coeff_426_10000.mat')

rng default


%% Plots
for i=1:426
    i
    f_name = ['f_10000_800_426_',num2str(i),'.mat'];
    load(f_name)
    LMn = ['LM_10000_800_426_',num2str(i),'.mat'];
    load(LMn)

    %% Plot ODF

%     figure
%     plot_odf_from_pvalues(idSample,fhat,p)
%     ori_all = LM_continue.ori;
ids1 = LM_continue.ids{1,/
%     ori1 = ori_all(ids1,:);
%     ori2 = ori_all(ids2,:);

%     plot_odf(ori1) % Plot initial ODF from orientations
% 
%     plot_odf(ori2) % Plot final ODF from orientations
    

    %% Plot point in Texture Hull
%     figure
%     hS = trisurf(coeff.tri,real(coeff.fhat_sym(coeff.id1,:)),imag(coeff.fhat_sym(coeff.id1,:)),real(coeff.fhat_sym(coeff.id2,:)));
%     axis equal tight vis3d
%     hS.FaceColor = 'g';
%     hS.EdgeColor = 'none';
%     hS.FaceAlpha = 0.5;
%     hold on
%     camlight
%     lighting gouraud
%     
%     
    times = [1 121];
    for j=1:length(times)
        jj = times(j);
        fhat_sym_Arbitrary2 = f_coeff(:,jj);

        p(j,1) = real(fhat_sym_Arbitrary2(coeff.id1));
        p(j,2) = imag(fhat_sym_Arbitrary2(coeff.id1));
        p(j,3) = real(fhat_sym_Arbitrary2(coeff.id2));
%         if j == 1  dot_color = 'b', else dot_color = 'r', end
% 
%         scatter3(p(j,1),p(j,2),p(j,3),30,dot_color,'filled','markeredgecolor','none')
% 
    end

    %% Get distance
    
    p1 = p(1,:); p2 = p(2,:);

    D(i,1) = norm(p1-p2);

end

%% Plot Histogram

D_sorted = sort(D);

figure
histogram(D)

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

% print(gcf,'BRK_Random_initial.png','-dpng','-r500');
% mtexColorbar
a = gca;
idx2kill = 1;
delete(a(idx2kill))

c = colorbar;
w = c.LineWidth;
c.LineWidth = 1.5;
c.FontSize = 30;
CLim(gcm,'equal');

end

function plot_odf_from_pvalues(idSample,fhat,p)

    CS = crystalSymmetry('Oh');
    smax = 32;
    smn = nan(0,3);
    for s = 0:smax
        [n,m] = ndgrid(-s:s,-s:s); % NOTE the order so that m is row index and n is column index
        smn = [smn; s*ones(numel(m),1) m(:) n(:)];
    end
    
    % make ODF
    odf = SO3FunHarmonic(((1-(smn(:,1)/(smax+1))).^2).*fhat*p(:,idSample),CS);
    
    plotPDF(odf,Miller(1,0,0,crystalSymmetry('1'))); 
    mtexColorbar

end