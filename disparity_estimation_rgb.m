%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Trung-Hieu Tran @ IPVS - Uni Stuttgart
% Variational disparity estimation framework for plenoptic images
% ICME 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%define lightfield data name
dataname = 'buddha';

view_skip = 0; %skip some border view for faster processing.

lf_file = sprintf('inputs/%s.mat',dataname);
data = load(lf_file);
lf = double(data.LF);
lf = lf/max(lf(:));
lf = lf(1+view_skip:end-view_skip,1+view_skip:end-view_skip,:,:,:);

% prepare paramam
[sy,sx,ny,nx,nc] = size(lf);
cs =floor((sx+1)/2);
ct =floor((sy+1)/2);
Ic = squeeze(lf(ct,cs,:,:,:));
disp(sprintf('      working with LF size %d x %d x %d x %d',sy,sx,ny,nx));

% paramters
hx = 1;
hy = 1;
alpha = 1;
sor = 1.78;
iteration = 100;
sigma = 0.3;
warp_max_level = 11;
warp_eta = 0.8;
gamma = 1;
warp_sigma = sqrt(2)/4/warp_eta;


% prepare smooth kernel.
ker_size = ceil(3*sigma)*2-1;
gauss_kernel = fspecial('gaussi',ker_size,sigma);

ker_size = ceil(3*warp_sigma)*2-1;
warp_gauss_kernel = fspecial('gaussi',ker_size,warp_sigma);

% compute alpha for each level.
warp_idx = warp_max_level:-1:0;
alpha_warps = alpha*nx ./ floor(nx*warp_eta.^warp_idx);

if(exist('targetMaxLevel')==0 || targetMaxLevel~=warp_max_level )
   % compute resolution for each level.
   warp_idx = warp_max_level:-1:0;
   nx_warps = floor(nx*warp_eta.^warp_idx);
   ny_warps = floor(ny*warp_eta.^warp_idx);
   res_warps = [ny_warps; nx_warps]';

   % compute units for each level
   nx_warps = hx * nx ./ floor(nx*warp_eta.^warp_idx);
   ny_warps = hy * ny ./ floor(ny*warp_eta.^warp_idx);
   unit_warps = [ny_warps; nx_warps]';


   % downscale lf and store in cell first
   lf_cells = cell(1,warp_max_level+1);
   disp(sprintf('... PREPARE LF data ...'));
   tmpLF = lf;
   for j=1:warp_max_level
      i = warp_max_level-j+1;
      ny_fine = res_warps(i+1,1);
      nx_fine = res_warps(i+1,2);
      ny_coarse = res_warps(i,1);
      nx_coarse = res_warps(i,2);
      disp(sprintf('      downscale from res %d %d to res: ny nx = %d %d ...',ny_fine,nx_fine,ny_coarse,nx_coarse));
      restmpLF = zeros(sy,sx,ny_coarse,nx_coarse,3);

      % filter with kernel
      display('      Filter LF before downscale');
      tmpLFBlur = zeros(size(tmpLF));
      for s=1:sx
          for t=1:sy
             for c=1:3
                tmpLFBlur(t,s,:,:,c) = imfilter(squeeze(tmpLF(t,s,:,:,c)),warp_gauss_kernel,'replicate');
             end
          end
       end
       display('      Downscaling LF');
      for c=1:3
         tmpimg = squeeze(tmpLFBlur(:,:,:,:,c));
         restmpLF(:,:,:,:,c) = lf_downscale([sx sy nx_fine ny_fine],[nx_coarse ny_coarse],tmpimg);
      end
      tmpLF = restmpLF;
      lf_cells(1,i) = {tmpLF};
   end
   lf_cells(1,warp_max_level+1) = {lf};
   clear tmpLF;

   targetMaxLevel = warp_max_level;
end


fig_commons = figure;
figure(fig_commons);
subplot(1,2,1); imagesc(Ic);
for i=1:warp_max_level+1
   disp(sprintf('... WARP LEVEL %d ...',i));
   ny_fine = res_warps(i,1);
   nx_fine = res_warps(i,2);
   hy_fine = unit_warps(i,1);
   hx_fine = unit_warps(i,2);

   if(exist('tmpLF')==1)
      clear tmpLF;
   end
   tmpLF = cell2mat(lf_cells(1,i));

   disp(sprintf('      working with resolution: ny nx = %d %d, units %f %f...',ny_fine,nx_fine,hy_fine, hx_fine));
   if(i==1)
      disp(sprintf('      init z with zeros ...'));
      z_init = zeros(ny_fine,nx_fine);
      %z_init = zzxx
      lf_warp =tmpLF;
   else
      ny_coarse = res_warps(i-1,1);
      nx_coarse = res_warps(i-1,2);
      disp(sprintf('      interpolate z from (%d,%d) to (%d,%d) ...',ny_coarse,nx_coarse,ny_fine,nx_fine));
      fz = z;

      z_init = img_resample([nx_coarse ny_coarse],[nx_fine ny_fine],fz);
      %warp lf with new result
      disp(sprintf('      warp lf with interpolate result'));
      if(exist('lf_warp')==1)
         clear lf_warp;
      end
      lf_warp =zeros(sy,sx,ny_fine,nx_fine,3);
      for c=1:3
         lf_warp(:,:,:,:,c) = lf_warp_with_z([sx sy nx_fine ny_fine],[hx_fine hy_fine],squeeze(tmpLF(:,:,:,:,c)),z_init);
      end
   end


   % blur before compute motion tensor
   disp(sprintf('      fillter LF with sigma %f size %d',sigma,ker_size));
   %lfblur = lf_warp; %zeros(size(lf_warp));
   % filter with kernel
   lfblur = zeros(size(lf_warp));
   for s=1:sx
      for t=1:sy
         lfblur(t,s,:,:) = imfilter(squeeze(lf_warp(t,s,:,:)),gauss_kernel,'replicate');
      end
   end


   disp(sprintf('      compute gray motion tensor'));
   gJ11 = zeros(ny_fine,nx_fine,3);
   gJ12 = zeros(ny_fine,nx_fine,3);
   gJ22 = zeros(ny_fine,nx_fine,3);
   %lfblur = lfblur*255;
   clf = squeeze(lfblur(:,:,:,:,1));
   [J11 J22 J12] = comp_j_gray_hue(clf,[hx_fine hy_fine]);
   gJ11(:,:,1) = J11; gJ12(:,:,1) = J12; gJ22(:,:,1) = J22;
   disp('... ... channel 1 ok ');
   clf = squeeze(lfblur(:,:,:,:,2));
   [J11 J22 J12] = comp_j_gray(clf,[hx_fine hy_fine]);
   gJ11(:,:,2) = J11; gJ12(:,:,2) = J12; gJ22(:,:,2) = J22;
   disp('... ... channel 2 ok ');
   clf = squeeze(lfblur(:,:,:,:,3));
   [J11 J22 J12] = comp_j_gray(clf,[hx_fine hy_fine]);
   gJ11(:,:,3) = J11; gJ12(:,:,3) = J12; gJ22(:,:,3) = J22;
   disp('... ... channel 3 ok ');
   gJ11 = sum(gJ11(:,:,:),3)/3;
   gJ22 = sum(gJ22(:,:,:),3)/3;
   gJ12 = sum(gJ12(:,:,:),3)/3;


   disp(sprintf('      compute gradient motion tensor'));
   GJ11 = zeros(ny_fine,nx_fine,3);
   GJ12 = zeros(ny_fine,nx_fine,3);
   GJ22 = zeros(ny_fine,nx_fine,3);
   %lfblur = lfblur*255;
   clf = squeeze(lfblur(:,:,:,:,1));
   [J11 J22 J12] = comp_j_gradient_hue(clf,[hx_fine hy_fine]);
   GJ11(:,:,1) = J11; GJ12(:,:,1) = J12; GJ22(:,:,1) = J22;
   disp('... ... channel 1 ok ');
   clf = squeeze(lfblur(:,:,:,:,2));
   [J11 J22 J12] = comp_j_gradient(clf,[hx_fine hy_fine]);
   GJ11(:,:,2) = J11; GJ12(:,:,2) = J12; GJ22(:,:,2) = J22;
   disp('... ... channel 2 ok ');
   clf = squeeze(lfblur(:,:,:,:,3));
   [J11 J22 J12] = comp_j_gradient(clf,[hx_fine hy_fine]);
   GJ11(:,:,3) = J11; GJ12(:,:,3) = J12; GJ22(:,:,3) = J22;
   disp('... ... channel 3 ok ');
   GJ11 = sum(GJ11(:,:,:),3)/3;
   GJ22 = sum(GJ22(:,:,:),3)/3;
   GJ12 = sum(GJ12(:,:,:),3)/3;

   clear J11; clear J12; clear J22;

   cur_alpha =  alpha_warps(i);
   disp(sprintf('      iterative solver for alpha %f sor %f iteration %d',cur_alpha, sor, iteration));
   dz = warp_solver_gs_robust_data_flow_driven([nx_fine ny_fine],[hx_fine hy_fine],[cur_alpha gamma 0.0000001 sor iteration],gJ11, gJ22, gJ12, GJ11, GJ22, GJ12, z_init);

   z = dz+z_init;

   figure(fig_commons);
   subplot(1,2,2);
   imagesc(z); colorbar;title(sprintf('Estimated solution at %d level',i));

   tw=waitforbuttonpress;
end

figure; imagesc(z); colorbar;
