% compute motion tensor for single channel
function [J11 J22 J12]=comp_j_gray(lf,units)

   sy = size(lf,1);
   sx = size(lf,2);
   ny = size(lf,3);
   nx = size(lf,4);

   cs = floor((sx+1)/2);
   ct = floor((sy+1)/2);

   xi  = 0.0001;


   J11 = zeros(ny,nx);
   J12 = zeros(ny,nx);
   J22 = zeros(ny,nx);

   %center view extraction
   f0 = squeeze(lf(ct,cs,:,:));
   [f0x f0y] = img_gradient_xy(f0,units);
   %W = f0x.^2 + f0y.^2 +xi;

   weight = 1.0/(sx*sy-1);

   for s=1:sx
      for t=1:sy
         vs = s - cs;
         vt = t - ct;
         if( vs ==0 && vt ==0)
            continue;
         end
         %extract view;
         f1 = squeeze(lf(t,s,:,:));
         [f1x f1y] = img_gradient_xy(f1,units);
         fx = (f0x+f1x)/2.0;
         fy = (f0y+f1y)/2.0;

         ft = (f1-f0);% TODO
         % compute dataterm weight;
         W = fx.^2 + fy.^2 +xi;
         %W = ones(size(fx));

         tmp = fx*vs + fy*vt;
         J11 = J11 + tmp.*tmp ./ W;
         J12 = J12 + ft.*tmp ./ W ;
         J22 = J22 + ft .* ft ./W;
      end
   end

   J11 =  J11*weight;
   J12 =  J12*weight;
   J22 =  J22*weight;

end
