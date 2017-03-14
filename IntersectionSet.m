function [index,Lvec,linearInd]=IntersectionSet(Source,Detector,xbox,ybox,theta)
global  omega m plotTravel BeforeEmit dz
global   Tol
x_coor=linspace(omega(1),omega(2),m(1)+1);
y_coor=linspace(omega(3),omega(4),m(2)+1);
slope=(Detector(2)-Source(2))/(Detector(1)-Source(1));
cut=(Source(2)*Detector(1)-Detector(2)*Source(1))/(Detector(1)-Source(1));

if(theta==pi/2)
    Q=[repmat(Source(1),size(y_coor')),y_coor'];
elseif(theta==0 | theta==2*pi)
    Q=[x_coor',repmat(Source(2),size(x_coor'))];
elseif(theta==pi)
    Q=[x_coor(end:-1:1)',repmat(Source(2),size(x_coor'))];
elseif(theta==3*pi/2)
    Q=[repmat(Source(1),size(y_coor')),y_coor(end:-1:1)'];
    
else
    Q=[[x_coor', slope*x_coor'+cut];...
        [(y_coor'-cut)./slope,y_coor']];
end

indx=find(Q(:,1)-xbox(1)<-1e-6*Tol |Q(:,1)-xbox(3)>1e-6*Tol); 
indy=find(Q(:,2)-ybox(1)<-1e-6*Tol |Q(:,2)-ybox(2)>1e-6*Tol);
dis=sqrt(sum(bsxfun(@minus,Q,Source).^2,2));
[~,InterOrder]=sort(dis);
Q=Q(InterOrder,:);
Lvec=sqrt(bsxfun(@minus,Q(2:end,1),Q(1:end-1,1)).^2+bsxfun(@minus,Q(2:end,2),Q(1:end-1,2)).^2);
%%%%%%%%%================================================================
QC=(Q(1:end-1,:)+Q(2:end,:))/2;
index=floor([(QC(:,1)-omega(1))/dz(1)+1, (QC(:,2)-omega(3))/dz(2)+1]);
indInside=find(index(:,1)>0 & index(:,1)<=m(2)& index(:,2)<=m(1) & index(:,2)>0);
index=index(indInside,:);
[~,subInd]=unique(index,'rows');
index=index(sort(subInd),:);
Lvec=Lvec(indInside);
Lvec=Lvec(sort(subInd));

linearInd=sub2ind(m,index(:,2),index(:,1));

