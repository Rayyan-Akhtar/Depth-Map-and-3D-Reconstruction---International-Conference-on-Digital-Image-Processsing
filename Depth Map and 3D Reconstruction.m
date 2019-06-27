% This code is copyright to the authors Mohd. Rayyan Akhtar from South China University of Technology.
% Department of Information and Communication Engineering

clc; clear; close all;

load logitech640x480.mat;   % Camera Calibration process file

% Path of the image file
path = '****************************************************';
fileName = '********.jpg';
Im = imread(strcat(path,fileName)); 
G = rgb2gray(Im);
figure, imshow(Im), title('\bf Input Image');

% Rotation matrix and translation vector
R = [0.0588503089179202, -0.998247005053399, 0.00628951843818058;
     -0.00500082358206347, -0.00659516430690762, -0.999965747199006;
     0.998254292704903, 0.0588168403578750, -0.00538018462042074];

t = [-0.00418039909319398;
     0.0693266037616307;
     -0.0669253851347902];

pointCl = [];

% Path of the LiDAR CSV files.
% Name file with offset like 0-4May, 2-4May, so on.
% It will read the file sequentially.
path = '********************************';
for i=0:2:150
    fileName = strcat(int2str(i),'-4May.csv');
    A = offsetZ(strcat(path,fileName),i);
    pointCl = [pointCl; A];
    fprintf('Read file no. %d\n',i/2);
end

x = pointCl(:,1);
y = pointCl(:,2);
z = pointCl(:,3);

ptCloud = pointCloud(pointCl);
figure, pcshow(ptCloud), title('\bf Point Cloud');
xlabel('\bf X(m)');
ylabel('\bf Y(m)');
zlabel('\bf Z(m)');

uv = worldToImage(logitech640x480, R', t, [x,y,z]);

i=1;
len = 0;
while i<=length(uv) 
    if uv(i,1)>=1 && uv(i,1)<=640 && uv(i,2)>=1 && uv(i,2)<=480
       len=len+1;
    end
    i=i+1;
end
uv_filtered = zeros(len,3);
r_filtered = zeros(len,1);
M = zeros(len,6);
n=1; i=1;
while i<=length(uv) 
    if uv(i,1)>=1 && uv(i,1)<=640 && uv(i,2)>=1 && uv(i,2)<=480
       uv_filtered(n,1:3) = [uv(i,1),uv(i,2),i];
       r_filtered(n,1) = x(i,1);
       M(n,1:6) = [x(i,1),y(i,1),z(i,1),double(Im(int16(uv(i,2)),int16(uv(i,1)),1)),...
               double(Im(int16(uv(i,2)),int16(uv(i,1)),2)),...
               double(Im(int16(uv(i,2)),int16(uv(i,1)),3))];
       n=n+1;
    end
    i=i+1;
end

% Saves the Point cloud in a text file in current folder
save('xyzrgb_4_May.txt','M','-ascii');

menu = menu('Select algorithm', 'GPR', 'Proposed', 'Run one after another');
list = {'10', '16', '20', '30', '32', '40', '48'};
patch_size_menu = listdlg('ListString', list);
switch(patch_size_menu)
    case 1
        patch_size = 10;
    case 2
        patch_size = 16;
    case 3
        patch_size = 20;
    case 4
        patch_size = 30;
    case 5
        patch_size = 32;
    case 6
        patch_size = 40;
    case 7
        patch_size = 48;
end

% ************** GPR Algorithm ************* %
if menu == 1
    empty_image = uint8(ones(480,640,3)*255);
    depth_intensity = zeros(480,640);
    uv_scatter = int16(uv_filtered);
    for i=1:length(uv_scatter)
       empty_image(uv_scatter(i,2),uv_scatter(i,1),1:3) = ...
           Im(uv_scatter(i,2),uv_scatter(i,1),1:3);
       depth_intensity(uv_scatter(i,2),uv_scatter(i,1)) = r_filtered(i,1);
    end
    
    figure, imshow(empty_image), title('\bf Projected Point Cloud');
    gpr_image = empty_image;

    counter=0;

    tic;
    c=0; d=1; a=0; b=1;
    for outer=1:480/patch_size  %24 %16 for 30
        for inner=1:640/patch_size  %32 %21 for 30
            fprintf('Outer: %d Inner: %d\n',outer,inner);
              empty=0; non_empty = 0;

            for row = 1+c*patch_size:d*patch_size
                for col=1+a*patch_size:b*patch_size
                   if depth_intensity(row,col) == 0
                     empty = empty+1;
                   else
                     non_empty = non_empty+1;
                   end
                end    
            end
            depth_empty = zeros(empty,3);
            depth_non_empty=zeros(non_empty,4);
            empty=1; non_empty=1;
            for row = 1+c*patch_size:d*patch_size
                for col=1+a*patch_size:b*patch_size
                   if depth_intensity(row,col) == 0
                     depth_empty(empty,1:2) = [row, col];
                     depth_empty(empty,3) = G(row,col);
                     empty = empty+1;
                  else
                     depth_non_empty(non_empty,1:2) = [row, col];
                     depth_non_empty(non_empty,3) = G(row, col);
                     depth_non_empty(non_empty,4) = depth_intensity(row, col);
                     non_empty = non_empty+1;
                   end
                end    
            end
              if isempty(depth_non_empty)

              else
                  K = calculateCovarienceMatrix(depth_non_empty,depth_non_empty,2,2);
                  K_star = calculateCovarienceMatrix(depth_non_empty,depth_empty,2,2);

                  mu = K_star'*(K\depth_non_empty(:,4));

                    [val, ind] = max(depth_non_empty(:,4));
                    [val_min, ind_min] = min(depth_non_empty(:,4));
                    if (val - val_min)>1.5
                        val = mean(depth_non_empty(:,4));
                    end
                  for i=1:length(mu)
                       [~, ind] = max(K_star(:,i));
                       if mu(i,1)<=val && mu(i,1)>val-1.2
                           depth_intensity(depth_empty(i,1),depth_empty(i,2)) = mu(i,1);
                           gpr_image(depth_empty(i,1),depth_empty(i,2),:) = ...
                           gpr_image(depth_non_empty(ind,1),depth_non_empty(ind,2),:);
                           counter = counter + 1;
                       end
                  end
              end
            a=a+1; b=b+1;
        end
       a=0; b=1;
       c=c+1; d=d+1;
    end
    
    fprintf('Pixels with depth: %f\n', length(uv_scatter));
    fprintf('Pixels with estimated depth: %f\n',counter);
    len_uv = length(uv_scatter);
    save('GPR Result.txt', 'patch_size', 'len_uv', 'counter', '-append', '-ascii');
    
    figure, imshow(gpr_image), title(sprintf('Reconstructed Image, Patch: %d',patch_size));
    time_gpr = toc;
    fprintf('time_gpr: %f\n',time_gpr);

    depth_image = uint8(zeros(480,640));

    gap_ = 0.02:0.02:4.4;

    color_depth = (uint8(linspace(255,32,221)))';
    gap_ = gap_';
    for row=1:480
       for col=1:640
           if(depth_intensity(row,col)==0)
              depth_image(row,col)=255; 
           else 
              for i=2:length(gap_)
                  if depth_intensity(row,col)>=gap_(i-1,1) && depth_intensity(row,col)<gap_(i,1)
                     depth_image(row,col) = color_depth(i,1);
                  end
              end
           end

       end
    end

    figure, imshow(depth_image), title(sprintf('Depth Image, Patch: %d',patch_size));
end

% ************** Proposed Algorithm *************** %
if menu == 2
    empty_image = uint8(ones(480,640,3)*255);
    depth_intensity = zeros(480,640);
    uv_scatter = int16(uv_filtered);
    for i=1:length(uv_scatter)
       empty_image(uv_scatter(i,2),uv_scatter(i,1),1:3) = ...
           Im(uv_scatter(i,2),uv_scatter(i,1),1:3);
       depth_intensity(uv_scatter(i,2),uv_scatter(i,1)) = r_filtered(i,1);
    end
    
    figure, imshow(empty_image), title('\bf Projected Point Cloud');
    gpr_image = empty_image;
    
    counter=0;
%     patch_size = 16;
    tic;
    c=0; d=1; a=0; b=1;
    for outer=1:480/patch_size  %24 %16 for 30
        for inner=1:640/patch_size  %32 %21 for 30
            fprintf('Outer: %d Inner: %d\n',outer,inner);
               empty=0; non_empty = 0;

            for row = 1+c*patch_size:d*patch_size
                for col=1+a*patch_size:b*patch_size
                   if depth_intensity(row,col) == 0
                     empty = empty+1;
                   else
                     non_empty = non_empty+1;
                   end
                end    
            end
            depth_empty = zeros(empty,3);
            depth_non_empty=zeros(non_empty,4);
            empty=1; non_empty=1;
            for row = 1+c*patch_size:d*patch_size
                for col=1+a*patch_size:b*patch_size
                   if depth_intensity(row,col) == 0
                     depth_empty(empty,1:2) = [row, col];
                     depth_empty(empty,3) = G(row,col);
                     empty = empty+1;
                  else
                     depth_non_empty(non_empty,1:2) = [row, col];
                     depth_non_empty(non_empty,3) = G(row, col);
                     depth_non_empty(non_empty,4) = depth_intensity(row, col);
                     non_empty = non_empty+1;
                   end
                end
            end
              if isempty(depth_non_empty)
                  %continue;
              else
   
                  K_star = calculateCovarienceMatrix(depth_non_empty,depth_empty,100,2);
   
                    K_star_size = size(K_star);
    
                  for i=1:K_star_size(1,2)
                       depth_intensity(depth_empty(i,1),depth_empty(i,2)) = mu(i,1);

                       [val, ind] = max(K_star(:,i));
                       depth_intensity(depth_empty(i,1),depth_empty(i,2))=depth_non_empty(ind,4);
                       gpr_image(depth_empty(i,1),depth_empty(i,2),:) = ...
                       gpr_image(depth_non_empty(ind,1),depth_non_empty(ind,2),:);
                       counter = counter + 1;
                  end
              end
            a=a+1; b=b+1;
        end
       a=0; b=1;
       c=c+1; d=d+1;
    end
    fprintf('Pixels with depth: %f\n', length(uv_scatter));
    fprintf('Pixels with estimated depth: %f\n',counter);
    len_uv = length(uv_scatter);
    save('My Algo Result.txt', 'patch_size', 'len_uv', 'counter', '-append', '-ascii');
    
    figure, imshow(gpr_image), title(sprintf('Reconstructed Image, Patch: %d',patch_size));
    time_gpr = toc;
    fprintf('time_gpr: %f\n',time_gpr);

    depth_image = uint8(zeros(480,640));

    gap_ = 0.02:0.02:4.4;

    color_depth = (uint8(linspace(255,50,221)))';
    gap_ = gap_';
    for row=1:480
       for col=1:640
           if(depth_intensity(row,col)==0)
              depth_image(row,col)=255; 
           else 
              for i=2:length(gap_)
                  if depth_intensity(row,col)>=gap_(i-1,1) && depth_intensity(row,col)<gap_(i,1)
                     depth_image(row,col) = color_depth(i,1);
                  end
              end
           end

       end
    end

    figure, imshow(depth_image), title(sprintf('Depth Image, Patch: %d',patch_size));
end

% *************** Do not change after this line *************** %
function [A] = offsetZ(fileName, offset_in_mm)
    B = readLidarDataAll(fileName,13);
    B(:,3) = B(:,3)-offset_in_mm/1000;
    B(:,6) = B(:,3);
    A(:,1) = B(:,1);
    A(:,2) = B(:,2);
    A(:,3) = B(:,3);
end

