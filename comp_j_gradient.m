% compute motion tensor for single channel
function [J11 J22 J12]=comp_j_gradient(lf,units)

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
   %[fxx fxy] = img_gradient_xy(f0x,units);
   %[fxy fyy] = img_gradient_xy(f0y,units);
   %Wx = fxx.^2 + fxy.^2+xi;
   %Wy = fyy.^2 + fxy.^2+xi;

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

         %compute first derivative
         [f1x f1y] = img_gradient_xy(f1,units);
         fx = (f0x+f1x)/2.0;
         fy = (f0y+f1y)/2.0;
         ft = (f1-f0);% TODO

         %compute second derivative
         [fxx fxy] = img_gradient_xy(fx,units);
         [fxy fyy] = img_gradient_xy(fy,units);
         [fxt fyt] = img_gradient_xy(ft,units);


         % compute dataterm weight;
         %W = fx.^2 + fy.^2 +xi;
         %W = ones(size(f0));
         Wx = fxx.^2 + fxy.^2+xi;
         %Wx =ones(size(f0));
         Wy = fyy.^2 + fxy.^2+xi;
         %Wy =ones(size(f0));

         tmpx = fxx*vs + fxy*vt;
         tmpy = fxy*vs + fyy*vt;

         J11 = J11 + (tmpx.*tmpx ./Wx  + tmpy.*tmpy ./Wy);
         J12 = J12 + (tmpx.*fxt  ./Wx  + tmpy.*fyt  ./Wy);
         J22 = J22 + (fxt .*fxt  ./Wx  + fyt .*fyt  ./Wy);
      end
   end

   J11 =  J11*weight;
   J12 =  J12*weight;
   J22 =  J22*weight;

end
