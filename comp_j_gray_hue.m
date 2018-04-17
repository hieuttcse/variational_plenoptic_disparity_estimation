% compute motion tensor for single channel
function [J11 J22 J12]=comp_j_gray_hue(lf,units)

   sy = size(lf,1);
   sx = size(lf,2);
   ny = size(lf,3);
   nx = size(lf,4);

   cs = floor((sx+1)/2);
   ct = floor((sy+1)/2);

   xi  = 0.0001;


   coslf = cos(lf);
   sinlf = sin(lf);

   J11 = zeros(ny,nx);
   J12 = zeros(ny,nx);
   J22 = zeros(ny,nx);

   %center view extraction
   cf0 = squeeze(coslf(ct,cs,:,:));
   [cf0x cf0y] = img_gradient_xy(cf0, units);

   sf0 = squeeze(sinlf(ct,cs,:,:));
   [sf0x sf0y] = img_gradient_xy(sf0, units);
   %W = f0x.^2 + f0y.^2 +xi;

   weight = 1.0/(sx*sy-1);
%for cos
   for s=1:sx
      for t=1:sy
         vs = s - cs;
         vt = t - ct;
         if( vs ==0 && vt ==0)
            continue;
         end

         %extract view;
         cf1 = squeeze(coslf(t,s,:,:));
         [cf1x cf1y] = img_gradient_xy(cf1,units);
         cfx = (cf0x+cf1x)/2.0;
         cfy = (cf0y+cf1y)/2.0;

         cft = (cf1-cf0);% TODO
         % compute dataterm weight;
         cW = cfx.^2 + cfy.^2;
         %W = ones(size(fx));
         ctmp = cfx*vs + cfy*vt;

         %extract view;
         sf1 = squeeze(sinlf(t,s,:,:));
         [sf1x sf1y] = img_gradient_xy(sf1,units);
         sfx = (sf0x+sf1x)/2.0;
         sfy = (sf0y+sf1y)/2.0;

         sft = (sf1-sf0);% TODO
         % compute dataterm weight;
         sW = sfx.^2 + sfy.^2;
         %W = ones(size(fx));
         stmp = sfx*vs + sfy*vt;

         W = cW + sW + xi;

         J11 = J11 + (ctmp.*ctmp + stmp.*stmp) ./ W;
         J12 = J12 + (cft.*ctmp + sft.*stmp) ./ W ;
         J22 = J22 + (cft .* cft + sft .*sft) ./W;
      end
   end

   J11 =  J11*weight;
   J12 =  J12*weight;
   J22 =  J22*weight;

%clear to save mem.
   clear coslf sinlf;
   clear cf0 cf0x cf0y cf1 cf1x cf1y cfx cfy cft;
   clear sf0 sf0x sf0y sf1 sf1x sf1y sfx sfy sft;
   clear ctmp stmp cW sW weight;
end
