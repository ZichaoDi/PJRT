#include "utils.h"

[intersects] = intersectLinePolygon([Source(1) Source(2) Detector(1)-Source(1) Detector(2)-Source(2)], [xbox',ybox']);
Ax=intersects(:,1);Ay=intersects(:,2);
if(isempty(Ax) | length(unique(Ax))==1 & length(unique(Ay))==1)
    % fprintf('no intersection \n')
    index=[];
    Lvec=[];
    linearInd=[];
else
    A=unique([Ax,Ay],'rows');Ax=A(:,1);Ay=A(:,2);
    if(theta==pi/2)
        Q=[repmat(Ax(1),size(y_coor')),y_coor'];
    elseif(theta==0 | theta==2*pi)
        Q=[x_coor',repmat(Ay(1),size(x_coor'))];
    elseif(theta==pi)
        Q=[x_coor(end:-1:1)',repmat(Ay(1),size(x_coor'))];
    elseif(theta==3*pi/2)
        Q=[repmat(Ax(1),size(y_coor')),y_coor(end:-1:1)'];
        
    else
        Q=[[x_coor', (Ay(2)-Ay(1))/(Ax(2)-Ax(1)).*x_coor'+(Ay(1)*Ax(2)-Ay(2)*Ax(1))/(Ax(2)-Ax(1))];...
            [(y_coor'-(Ay(1)*Ax(2)-Ay(2)*Ax(1))/(Ax(2)-Ax(1)))./((Ay(2)-Ay(1))/(Ax(2)-Ax(1))),y_coor']];
    end
    
    indx=find(Q(:,1)-xbox(1)<-1e-6*Tol |Q(:,1)-xbox(3)>1e-6*Tol); 
    indy=find(Q(:,2)-ybox(1)<-1e-6*Tol |Q(:,2)-ybox(2)>1e-6*Tol);
    Q=setdiff(Q,Q([indx;indy],:),'rows');
    Q=unique(Q,'rows');
    if(BeforeEmit)
        dis=sqrt(sum(bsxfun(@minus,Q,Source).^2,2));
        [~,InterOrder]=sort(dis);
        Q=Q(InterOrder,:);
    else
        Q=unique([Q;A],'rows');
    end
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
    %%%%%%%%%================================================================
    if plotTravel
        if(BeforeEmit)
            
            figure(finalfig)
            subplot(1,2,1);
            set(fig2,'visible','off');
            drawnow;
            if(~isempty(index))
                fig2=plot((index(:,1)-1/2)*dz(1)-abs(omega(1)),(index(:,2)-1/2)*dz(2)-abs(omega(3)),'r*');%,Q(:,1),Q(:,2),'g-');
            end
        else
            figure(finalfig)
            subplot(1,2,1);
            set(fig5,'visible','off');
            drawnow;
            if(~isempty(index))
                fig5=plot((index(:,1)-1/2)*dz(1)-abs(omega(1)),(index(:,2)-1/2)*dz(2)-abs(omega(3)),'bo',Q(:,1),Q(:,2),'g-');
            end
            pause;
        end
    end
    %%%%%%%%%================================================================
    
end
