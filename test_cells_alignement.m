%% This code is designed to identify bioarchitectural parameters that shape 
%% the internal and spatial organization of tumours.
%% 
%% This code contains 3 tutorials:
%% 
%% - "test_cells_alignement.m": In this tutorial, each cell is characterized 
%% using its main axis. To this end, a Principal Component Analysis (PCA) is 
%% applied on the voxel coordinates within the segmented mask of a structure 
%% of interest. The main axes of cells is used to calculate the best 2D 
%% “alignment plane” (in the least-squares sense: the sum of squared differences 
%% is minimized between the observed main axes of cells, and the fitted value 
%% provides by a 2D plane equation). Angles between the main axis of each 
%% cell and the alignment plane are then calculated.
%% 
%% - "test_cells_polarisation.m": In this tutorial, each cell is characterized 
%% using its main axis. To this end, a Principal Component Analysis (PCA) is 
%% applied on the voxel coordinates within the segmented mask of a structure 
%% of interest. 3D virtual rays are emitted by each cell along its main axis 
%% in both directions (the ray radius is a user-defined input parameter for 
%% the algorithm). The accumulated beam intensity is then calculated on a 
%% voxel-by-voxel basis.
%% 
%% - 'test_volume_versus_vessel.m': In this tutorial, the size of different 
%% cell components is compared to the distance to a blood capillary.
%% 
%% This code has been written by Baudouin Denis de Senneville.
%% Institut de Mathématiques de Bordeaux, 
%% UMR 5231 CNRS,
%% Université de Bordeaux 351, cours de la Libération - F 33 405 TALENCE
%% 
%% This code was developed under the commercial software Matlab ©1994-2021 
%% The MathWorks, Inc.
%% 
%% Please cite the following paper if you are using this code:
%% 
%% Baudouin Denis de Senneville, Fatma Zohra Khoubai, Marc Bevilacqua, 
%% Alexandre Labedade, Kathleen Flosseau, Christophe Chardot, Sophie 
%% Branchereau, Jean Ripoche, Stefano Cairo, Etienne Gontier, 
%% Christophe F. Grosset,  Deciphering Tumour Tissue Organization by 3D 
%% Electron Microscopy and machine learning, bioRxiv 2021.06.15.446847; 
%% doi: https://doi.org/10.1101/2021.06.15.446847. 
%% https://www.biorxiv.org/content/10.1101/2021.06.15.446847v1.abstract

clc
close all;
clear all;

%% Load cell labels
load('./segmentation_data/label_cells.mat');
[dimx,dimy,dimz] = size(mask_cells);

%% Display cell labels
delta_display_slice = floor(dimz/9);
subplotNb = 1;
for num_slice=round(delta_display_slice/2):delta_display_slice:dimz-round(delta_display_slice/2)
    figure(1);
    subplot(3,3,subplotNb)
    imagesc(squeeze(mask_cells(:,:, num_slice))'); axis off; axis image;
    title(strcat('Slice #',int2str(num_slice)));
    subplotNb = subplotNb+1;
end
pause(0.1);

%% Compute the label Id of each cell
label_Id = unique(mask_cells(:));

%% Remove the background label
label_Id(label_Id==0) = [];

%% Get the number of cell
nb_cells = numel(label_Id);

%% Compute the eigen vector of each cell
for i = 1:nb_cells
    
    fprintf('Computing eigenvector of cell %d/%d\n',i,nb_cells);
    
    %% Get the mask of the current cell
    idx = find(mask_cells==label_Id(i));    
    current_mask_cell = false(dimx,dimy,dimz);
    current_mask_cell(idx) = true;
    
    %% Compute eigen vectors of the current cell
    var_for_pca = [];
    for z=1:dimz
        for y=1:dimy
            for x=1:dimx
                if (current_mask_cell(x,y,z)>0)
                    var_for_pca(:,end+1) = [ x, y, z ];
                end
            end
        end
    end
    [coeff,score,latent] = pca(var_for_pca');    
    
    %% Store the eigen vector of the current cell
    cell_eigen_vec(:,:,i) = coeff;
    
    clear var_for_pca;
    clear current_mask_cell;
    
end

%% Least Square fit of a 2D plan
A=double(squeeze(cell_eigen_vec(1:2,1,:)))';
b=double(squeeze(cell_eigen_vec(3,1,:)));
x=inv(A'*A)*(A'*b);
plane_coeff = [x(1), x(2), -1];

%% Compute cell alignment (i.e., angle of each cell with the fitted 2D plane)
for i=1:nb_cells
    V1 = double(squeeze(cell_eigen_vec(:,1,i)))';
    V2 = plane_coeff;
    cell_alignement(i) = abs(V1*V2')/(norm(V1)*norm(V2));
    cell_alignement(i) = 90-acos(cell_alignement(i))*180/pi;
end        

%% Display histogram of cell alignments
figure(2);
set(gcf, 'defaultAxesColorOrder', [[0 0 0];[0 0 1]]);
yyaxis left;
h=histogram(squeeze(cell_alignement),18,'BinLimits',[0,90]);
set(h(1), 'FaceColor', [0.33 0.33 0.33]);
xlim([0,90]);
grid on;xlabel('Angle Cells/Best alignment plane [degree]');
ylabel('Occurence [#]');
set(h(1),'linewidth',1);
hold on;
yyaxis right;
h3=plot(h.BinEdges,[0,100*cumsum(h.Values)/sum(h.Values)], 'b-', 'LineWidth', 1);
set(h3(1),'linewidth',2);
ylabel('Cumulated occurence [%]', 'Color', 'b');
ylim([0,100]);
grid on;
