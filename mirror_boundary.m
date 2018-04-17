% dim = [10 10 1 1];
% scale = [1.0 1.0];
% params = [100 1.5 100];
% J11 = ones(12,12);
% J22 = ones(12,12)*2;
% J12 = ones(12,12)*3;

function f0 = mirror_boundary(f,bx,by)
   f0 = f;
   if(by>0)
      f0(1:by,:)=repmat(f0(by+1,:),by,1);
      f0(end-by+1:end,:)=repmat(f0(end-by,:),by,1);
   end
   if(bx>0)
      f0(:,1:bx)=repmat(f0(:,bx+1),1,bx);
      f0(:,end-bx+1:end)=repmat(f0(:,end-bx),1,bx);
   end
end
