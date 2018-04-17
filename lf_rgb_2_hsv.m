function hsvlf = lf_rgb_2_hsv(lf)

   sy = size(lf,1);
   sx = size(lf,2);
   ny = size(lf,3);
   nx = size(lf,4);
   hsvlf=zeros(sy,sx,ny,nx,3);
   for s=1:sx
      for t=1:sy
         f0 = squeeze(lf(t,s,:,:,:));
         hsvlf(t,s,:,:,:) = rgb2hsv(f0);
      end
   end

end
