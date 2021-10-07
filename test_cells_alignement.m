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
