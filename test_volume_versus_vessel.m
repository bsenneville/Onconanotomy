

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
