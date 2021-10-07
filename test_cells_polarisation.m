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

%% Compute the eigen vector and gravity center of each cell
for i = 1:nb_cells
    
    fprintf('Computing eigenvector and gravity center of cell %d/%d\n',i,nb_cells);
    
    %% Get the mask of the current working ROI
    idx = find(mask_cells==label_Id(i));
    
    num_roi_cur = num_roi_cur+1;
    
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

%% Compute cumulated ray map
ray_cumul_map = zeros(dimx,dimy,dimz);
ray_radius = round(0.1*dimx); 
ray_length = 0.2;
ray_precision = 100; %% Number of observed location for each ray

%% Loop over each ray
for num_ray = 1:nb_cells
    
    fprintf('Cumulating ray %d/%d\n',num_ray,nb_cells);
    
    %% Initialize the current ray map
    ray_map = zeros(dimx,dimy,dimz);
    
    %% Compute ray tips
    radius_size = ray_length*[dimx,dimy,dimz]'.*squeeze(cell_eigen_vec(:,1,num_ray));
    P1 = gravity_center(:,num_ray)-radius_size;
    P2 = gravity_center(:,num_ray)+radius_size;
    
    %% Compute ray path (observed locations)
    for i=1:ray_precision
        ray_path_X(i) = P1(1)+(P2(1)-P1(1))*i/ray_precision;
        ray_path_Y(i) = P1(2)+(P2(2)-P1(2))*i/ray_precision;
        ray_path_Z(i) = P1(3)+(P2(3)-P1(3))*i/ray_precision;
    end
    ray_path_X = round(ray_path_X);
    ray_path_Y = round(ray_path_Y);
    ray_path_Z = round(ray_path_Z);
    
    %% Define a brush to draw the ray path
    [X,Y,Z] = meshgrid(1:2*ray_radius,1:2*ray_radius,1:2*round(ray_radius));
    my_brush = (ray_radius-sqrt((X-ray_radius).^2+(Y-ray_radius).^2+(Z-ray_radius).^2))/ray_radius;
    my_brush(find(my_brush<=0)) = 0;
    
    %% Draw the ray path using the brush
    for i=1:numel(ray_path_X)
        min_x = ray_path_X(i)-ray_radius;
        max_x = ray_path_X(i)+ray_radius-1;
        min_y = ray_path_Y(i)-ray_radius;
        max_y = ray_path_Y(i)+ray_radius-1;
        min_z = ray_path_Z(i)-round(ray_radius)+1;
        max_z = ray_path_Z(i)+round(ray_radius);
        delta_minX = 0;
        if (min_x<1)
            delta_minX = 1-min_x; min_x = 1;
        end
        delta_minY = 0;
        if (min_y<1)
            delta_minY = 1-min_y; min_y = 1;
        end
        delta_minZ = 0;
        if (min_z<1)
            delta_minZ = 1-min_z; min_z = 1;
        end
        delta_maxX = 0;
        if (max_x>dimx)
            delta_maxX = max_x-dimx; max_x = dimx;
        end
        delta_maxY = 0;
        if (max_y>dimy)
            delta_maxY = max_y-dimy; max_y = dimy;
        end
        delta_maxZ = 0;
        if (max_z>dimz)
            delta_maxZ = max_z-dimz; max_z = dimz;
        end
        
        aux = ray_map(min_x:max_x,min_y:max_y,min_z:max_z);
        my_brush_crop = my_brush(delta_minX+1:end-delta_maxX,delta_minY+1:end-delta_maxY,delta_minZ+1:end-delta_maxZ);
        index = find(aux<=my_brush_crop);
        aux(index) = my_brush_crop(index);
        ray_map(min_x:max_x,min_y:max_y,min_z:max_z) = aux;
    end
    
    %% Cumul current ray map
    ray_cumul_map = ray_cumul_map + ray_map;
    
end

%% Save RST map
subplotNb = 1;
for num_slice=round(delta_display_slice):delta_display_slice:dimz
    figure(2);
    subplot(3,3,subplotNb)
    imagesc(squeeze(ray_cumul_map(:,:,num_slice))',[0,max(ray_cumul_map(:))]); axis off; axis image;
    mask_slice = mask_cells(:,:,num_slice);
    region_id_cells = unique(mask_cells(:));
    for i = 1:nb_cells
        
        working_mask = false(size(mask_slice));
        idx = find(mask_slice==label_Id(i));
        if (numel(idx)>0)
            working_mask(idx) = true;
            hold on
            visboundaries(working_mask','LineStyle',':','LineWidth',0.2,'Color','k');
        end
    end
    working_mask = false(size(mask_slice));
    working_mask(2:end-1,2:end-1) = true;
    visboundaries(working_mask','LineStyle','-','LineWidth',1,'Color','k');
    subplotNb = subplotNb+1;
    pause(0.1);
end
