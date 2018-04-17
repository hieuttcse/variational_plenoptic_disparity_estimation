% compute motion tensor for single channel
function [J11 J22 J12]=comp_j_gradient_hue(lf,units)

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


   coslf = cos(lf);
   sinlf = sin(lf);

   %center view extraction
   cf0 = squeeze(coslf(ct,cs,:,:));
   [cf0x cf0y] = img_gradient_xy(cf0,units);

   sf0 = squeeze(sinlf(ct,cs,:,:));
   [sf0x sf0y] = img_gradient_xy(sf0,units);

   weight = 1.0/(sx*sy-1);

   for s=1:sx
      for t=1:sy
         vs = s - cs;
         vt = t - ct;
         if( vs ==0 && vt ==0)
            continue;
         end


         %extract view;
         cf1 = squeeze(coslf(t,s,:,:));

         %compute first derivative
         [cf1x cf1y] = img_gradient_xy(cf1,units);
         cfx = (cf0x+cf1x)/2.0;
         cfy = (cf0y+cf1y)/2.0;
         cft = (cf1-cf0);% TODO

         %compute second derivative
         [cfxx cfxy] = img_gradient_xy(cfx,units);
         [cfxy cfyy] = img_gradient_xy(cfy,units);
         [cfxt cfyt] = img_gradient_xy(cft,units);
         % compute dataterm weight;
         %W = fx.^2 + fy.^2 +xi;
         %W = ones(size(f0));
         cWx = cfxx.^2 + cfxy.^2;
         %Wx =ones(size(f0));
         cWy = cfyy.^2 + cfxy.^2;
         %Wy =ones(size(f0));

         %extract view;
         sf1 = squeeze(sinlf(t,s,:,:));

         %compute first derivative
         [sf1x sf1y] = img_gradient_xy(sf1,units);
         sfx = (sf0x+sf1x)/2.0;
         sfy = (sf0y+sf1y)/2.0;
         sft = (sf1-sf0);% TODO

         %compute second derivative
         [sfxx sfxy] = img_gradient_xy(sfx,units);
         [sfxy sfyy] = img_gradient_xy(sfy,units);
         [sfxt sfyt] = img_gradient_xy(sft,units);
         % compute dataterm weight;
         %W = fx.^2 + fy.^2 +xi;
         %W = ones(size(f0));
         sWx = sfxx.^2 + sfxy.^2;
         %Wx =ones(size(f0));
         sWy = sfyy.^2 + sfxy.^2;
         %Wy =ones(size(f0));

         Wx = sWx + cWx + xi;
         Wy = sWy + cWy + xi;

         ctmpx = cfxx*vs + cfxy*vt;
         ctmpy = cfxy*vs + cfyy*vt;

         stmpx = sfxx*vs + sfxy*vt;
         stmpy = sfxy*vs + sfyy*vt;

         J11 = J11 + ( ( stmpx.*stmpx + ctmpx.*ctmpx) ./Wx  + (ctmpy.*ctmpy + stmpy.*stmpy) ./Wy);
         J12 = J12 + ( ( ctmpx.*cfxt  + stmpx.*sfxt ) ./Wx  + (ctmpy.*cfyt  + stmpy.*sfyt ) ./Wy);
         J22 = J22 + ( ( cfxt .*cfxt + sfxt .*sfxt ) ./Wx  +  (cfyt .*cfyt  + sfyt .*sfyt ) ./Wy);
      end
   end

   J11 =  J11*weight;
   J12 =  J12*weight;
   J22 =  J22*weight;

end
