clc
clear
cd '/data/scripts/spinTest'

%%
coords = dlmread('sphere_HCP.txt');
coords_r = coords(1:180,:);
coords_l = coords(181:360,:);

% uncomment for generating new null distribution

% null_distr = rotate_parcellation(coords_l,coords_r,1000);
% save('null_distr.mat','null_distr')
load 'null_distr.mat'

%%

evo_sq = zeros(4,1);
%%%%%%%%%%%%%%%%%%%%%
x = csvread('/data/allometricScaling_regional.txt',0);
%y = csvread('/data/distortion_regional.txt',0);
y = csvread('/data/compressionEfficiency_send.txt');
evo_sq(1,:) = perm_sphere_p(x,y(:,2),null_distr,'spearman')

x = csvread('/data/myelin_regional.txt',0);
y = csvread('/data/compressionEfficiency_send.txt');
evo_sq(2,:) = perm_sphere_p(x,y(:,2),null_distr,'spearman')

x = csvread('/data/allometricScaling_regional.txt',0);
y = csvread('/data/compressionEfficiency_receive.txt');
evo_sq(3,:)= perm_sphere_p(x,y(:,2),null_distr,'spearman')

x = csvread('/data/myelin_regional.txt',0);
y = csvread('/data/compressionEfficiency_receive.txt');
evo_sq(4,:)= perm_sphere_p(x,y(:,2),null_distr,'spearman')

