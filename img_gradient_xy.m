function [imgx imgy] = img_gradient_xy(img,units)
   hx = units(1);
   hy = units(2);

   nx = size(img,2);
   ny = size(img,1);

   f = zeros(ny+2,nx+2);
   fx = zeros(size(f));
   fy = zeros(size(f));
   f(2:end-1,2:end-1) = img;
   f = mirror_boundary(f,1,1);

   fy(2:end-1,:) = (f(3:end,:) - f(1:end-2,:))/2/hy;
   fx(:,2:end-1) = (f(:,3:end) - f(:,1:end-2))/2/hx;

   imgx = fx(2:end-1,2:end-1);
   imgy = fy(2:end-1,2:end-1);
end
