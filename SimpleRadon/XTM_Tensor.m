%%%Simulate XRT of a given object with predifined detector and beam
global x y omega m dz Tol
N=[30];%[33 17 9 5 3];% 17 9];%[129 65  9 5];%
numThetan=60; % number of scanning angles/projections
Tol=1e-2; 
omega=[-2     2    -2     2].*Tol;
m=[N N]; %Numerical Resolution
dz=[(omega(2)-omega(1))/m(2) (omega(4)-omega(3))/m(1)];
%%% +++set up experimental configuration+++++++++++++
alpha=atan((omega(4)-omega(3))/(omega(2)-omega(1)));
Tau= omega(2)-omega(1);
nTau=ceil(sqrt(2*prod(m)))-1;%m(1)-1;% % number of discrete beam%nTau;%
tol1=1/2*N(1);
%=============initiate transmission detector location
detS0=[Tau/2*tan(alpha)+tol1*dz(1),-Tau/2-tol1*dz(1)]; 
detE0=[Tau/2*tan(alpha)+tol1*dz(1),Tau/2+tol1*dz(1)];
% dTau=abs(-Tau-2*tol1*dz(1))/(nTau+1);%%% width of each discrete beam
knot=linspace(detS0(2),detE0(2),nTau+1)';
DetKnot0=[repmat(detS0(1),size(knot)),knot];%% transmission detectorlet knot points
SourceS0=[-Tau/2*tan(alpha)-tol1*dz(1),-Tau/2-tol1*dz(1)];%initiate beam source
SourceE0=[-Tau/2*tan(alpha)-tol1*dz(1),Tau/2+tol1*dz(1)];

knot=linspace(SourceS0(2),SourceE0(2),nTau+1)';
SourceKnot0=[repmat(SourceS0(1),size(knot)),knot];%% source knot points
% DetKnot0=DetKnot0(end:-1:1,:);
% SourceKnot0=SourceKnot0(end:-1:1,:);
%%%=========== Assign Projection Angles;
thetan=linspace(1,360,numThetan);% must be positive.
%%%+++++++++++++++++++++++++++++++++++++++++++++++++++++
%%%++++++++ Set up sample +++++++++++++++++++++++++++
load PeriodicTable
x=linspace(omega(1),omega(2),m(1)+1);
y=linspace(omega(3),omega(4),m(2)+1);

%%%=========== assign weight matrix for each element in each pixel
sample = 'Phantom';
Z = [14];
NumElement=length(Z);
W=sum(phantom3d(N),3);
E0=12.1;
load(['xRayLib',num2str(E0),'.mat'])
load AtomicWeight
%%%%% =================== Attenuation Matrix at beam energy
MUe=reshape(CS_TotalBeam(Z(1),1),1,1,NumElement);
MU_XTM=sum(W.*repmat(MUe,[m(1),m(2),1]),3);
%% +++++++++++++++++++++++++++++++++++++++++++++++++++++++
eX=ones(m(1),1);
eY=ones(m(2),1);
DisR_Simulated=zeros(numThetan,nTau+1);
L=sparse(numThetan*(nTau+1),prod(m));
for n=1:numThetan
    %% =============== No Probe Drift
    theta = thetan(n)/180*pi;
    TransMatrix=[cos(theta) sin(theta);-sin(theta) cos(theta)];
    DetKnot=DetKnot0*TransMatrix;
    SourceKnot=SourceKnot0*TransMatrix;
    %% =========================================
    Rdis=ones(1,nTau+1);
    xbox=[omega(1) omega(1) omega(2) omega(2) omega(1)];
    ybox=[omega(3) omega(4) omega(4) omega(3) omega(3)];
    for i=1:nTau+1 %%%%%%%%%========================================================
        %=================================================================
        [index,Lvec,linearInd]=IntersectionSet(SourceKnot(i,:),DetKnot(i,:),xbox,ybox,theta);
        %%%%%%%%================================================================
        if(~isempty(index)& norm(Lvec)>0)
            currentInd=sub2ind(m,index(:,2),index(:,1));
            L(sub2ind([numThetan,nTau+1],n,i),currentInd)=Lvec;
            Rdis(i)=eX'*(MU_XTM.*reshape(L(sub2ind([numThetan,nTau+1],n,i),:),m))*eY;
        end

    end
    DisR_Simulated(n,:)=Rdis;
end
L=L*sqrt(N(1));
%%==============================================================
