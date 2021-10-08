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
load('./segmentation_data/label_cells_complete.mat');
[dimx,dimy,dimz] = size(mask_cells);

%% Load vessel mask
load('./segmentation_data/mask_vessel.mat');

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

%% Display vessel mask
subplotNb = 1;
for num_slice=round(delta_display_slice/2):delta_display_slice:dimz-round(delta_display_slice/2)
    figure(2);
    subplot(3,3,subplotNb)
    imagesc(squeeze(mask_vessels(:,:, num_slice))'); axis off; axis image;
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

for i = 1:nb_cells
    
    fprintf('Computing size (in voxels) of cell %d/%d\n',i,nb_cells);
    
    %% Get the mask of the current cell
    idx = find(mask_cells==label_Id(i));    
    current_mask_cell = false(dimx,dimy,dimz);
    current_mask_cell(idx) = true;
    
    %% Compute gravity center coordinate (in voxels)
    gravity_center_x = 0; gravity_center_y = 0; gravity_center_z = 0; count = 0;
    for z=1:dimz
        if (max(max(current_mask_cell(:,:,z)))>0)
            for y=1:dimy
                if (max(max(current_mask_cell(:,y,z)))>0)
                    for x=1:dimx
                        if (current_mask_cell(x,y,z)>0)
                            gravity_center_x = gravity_center_x + x;
                            gravity_center_y = gravity_center_y + y;
                            gravity_center_z = gravity_center_z + z;
                            count = count + 1;
                        end
                    end
                end
            end
        end
    end
    gravity_center(:,i) = [gravity_center_x/count, gravity_center_y/count, gravity_center_z/count];
    
    %% Compute region volume (in voxels)
    idx = find(current_mask_cell>0);
    region_volume(i) = numel(idx);
    
end


%% Compute distance to vessel map
[distance_vessels, closest_voxel_idx] = bwdist(mask_vessels);

%% Compute distance between vessel and gravity center of each cell
for i=1:nb_cells
    distance_vessels_GC(i) = distance_vessels(round(gravity_center(1,i)), ...
        round(gravity_center(2,i)), ...
        round(gravity_center(3,i)));
end

%% Display distance from each cell to vessel
figure(3);
plot(distance_vessels_GC,region_volume(:),'bo');
[p,s] = polyfit(distance_vessels_GC,region_volume,1);
hold on;
plot(distance_vessels_GC,polyval(p,distance_vessels_GC),'k-')
xlabel('Distance to vessel [Voxels]');
ylabel('Cell volume [Voxels]');
xlim([min(distance_vessels_GC(:)), max(distance_vessels_GC(:))]);
grid on;
h=legend('Measures', 'Linear fit', 'Location', 'NorthEast');
