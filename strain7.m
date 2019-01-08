%smaller step size

%define variables
x0=40; y0=40; z0=40; %size of space in nm
cx=0.25; cy=0.25; cz=0.25; %step size/ unit cubic size
cx1=0.25e-9; cy1=0.25e-9; cz1=0.25e-9;
nx=x0/cx; ny=y0/cy; nz=z0/cz;
sub=15; height=6; base=12; wet=0.5;
mismatch=(6.0583-5.65325)*1e-10; %300K
tol=0.15; steps=100; factor=0.1;

%initialization of 4D displacement field u and strain S
u=zeros(nx,ny,nz,3);
c_g=1e10*[11.9;5.34;5.96]; %dyn/cm2=0.1N/m2
c_i=1e10*[8.34;4.54;3.95];
c11=zeros(nx,ny,nz); c12=zeros(nx,ny,nz); c44=zeros(nx,ny,nz);
cinv11=zeros(nx,ny,nz); cinv12=zeros(nx,ny,nz); cinv44=zeros(nx,ny,nz);
E=zeros(steps,1);


%hessian
m=zeros(nx,ny,nz);
H.g=[c_g(1) c_g(2) c_g(2) 0 0 0; c_g(2) c_g(1) c_g(2) 0 0 0; c_g(2) c_g(2) c_g(1) 0 0 0; 0 0 0 4*c_g(3) 0 0; 0 0 0 0 4*c_g(3) 0; 0 0 0 0 0 4*c_g(3)];
inv1=inv(H.g);

H.i=[c_i(1) c_i(2) c_i(2) 0 0 0; c_i(2) c_i(1) c_i(2) 0 0 0; c_i(2) c_i(2) c_i(1) 0 0 0; 0 0 0 4*c_i(3) 0 0; 0 0 0 0 4*c_i(3) 0; 0 0 0 0 0 4*c_i(3)];
inv2=inv(H.i);

%substrate 1:60
for z=1:sub/cz
    m(:,:,z)=1;
end
%last layer of substrate %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for y=1:ny
    for x=1:nx
        u(x,y,sub/cz+1,3)=-0.5*mismatch;
        if x<=nx/2
            u(x,y,sub/cz+1,1)=-0.5*mismatch;
        else
            u(x,y,sub/cz+1,1)=0.5*mismatch;
        end
        if y<=ny/2
            u(x,y,sub/cz+1,2)=-0.5*mismatch;
        else
            u(x,y,sub/cz+1,2)=0.5*mismatch;
        end
    end
end
u(:,:,sub/cz,3)=-0.5*mismatch; 

%wetting layer 0.5nm 61:62
for z=sub/cz+1:(sub+wet)/cz
    m(:,:,z)=2;
    u(:,:,sub/cz+1,3)=-0.5*mismatch;
    u(:,:,(sub+wet)/cz,3)=0.5*mismatch; 
    for y=1:ny
        for x=1:nx
            if x<=nx/2
                u(x,y,z,1)=mismatch;
            else
                u(x,y,z,1)=-1*mismatch;
            end
            if y<=ny/2
                u(x,y,z,2)=mismatch;
            else
                u(x,y,z,2)=-1*mismatch;
            end
        end
    end
end

%layer after wetting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for y=1:ny
    for x=1:nx
        u(x,y,(sub+wet)/cz+1,3)=0.5*mismatch;
        if x<=nx/2
            u(x,y,(sub+wet)/cz+1,1)=-0.5*mismatch;
        else
            u(x,y,(sub+wet)/cz+1,1)=0.5*mismatch;
        end
        if y<=ny/2
            u(x,y,(sub+wet)/cz+1,2)=-0.5*mismatch;
        else
            u(x,y,(sub+wet)/cz+1,2)=0.5*mismatch;
        end
    end
end
u(:,:,(sub+wet)/cz+1,3)=0.5*mismatch;

%pyramid 6nm %63:86
checkx=((sub+wet+height)/cz+1);
for z=(sub+wet)/cz+1:(sub+wet+height)/cz
    for y=1:ny
        for x=1:nx
            if abs(y-(ny/2+0.5))<(checkx-z) && abs(x-(nx/2+0.5))<(checkx-z) 
                u(x,y,z,1)=mismatch;
                u(x,y,z,2)=mismatch; %%%%%%%%%%%%%%%%%%%%%%%%%%%
                u(x,y,z,3)=mismatch;
                m(x,y,z)=2;
            else
                m(x,y,z)=1;
            end
            if x>nx/2
                u(x,y,z,1)=-u(x,y,z,1);
            end
            if y>ny/2
                u(x,y,z,2)=-u(x,y,z,2);
            end
            if z>74 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                u(x,y,z,1)=-u(x,y,z,1);
                u(x,y,z,2)=-u(x,y,z,2);
            end
        end
    end
end

%capping 87:160
for z=(sub+wet+height)/cz+1:nz
    m(:,:,z)=1;
end

for x=1:nx
    for y=1:ny
        for z=1:nz
            if m(x,y,z)==1
                c11(x,y,z)=c_g(1); %used in grad
                c12(x,y,z)=c_g(2);
                c44(x,y,z)=c_g(3);
                cinv11(x,y,z)=inv1(1,1); %components of inverse of hessian, used in minimization
                cinv12(x,y,z)=inv1(1,2); 
                cinv44(x,y,z)=inv1(4,4);
            else
                c11(x,y,z)=c_i(1);
                c12(x,y,z)=c_i(2);
                c44(x,y,z)=c_i(3);
                cinv11(x,y,z)=inv2(1,1); 
                cinv12(x,y,z)=inv2(1,2); 
                cinv44(x,y,z)=inv2(4,4);
            end
        end
    end
end

x_f=[2:nx,nx]; y_f=[2:ny,ny]; z_f=[2:nz,nz];
x_b=[1,1:nx-1]; y_b=[1,1:ny-1]; z_b=[1,1:nz-1];
S1=(u(x_f,:,:,1)-u(x_b,:,:,1))/cx1/2;
S2=(u(:,y_f,:,2)-u(:,y_b,:,2))/cy1/2;
S3=(u(:,:,z_f,3)-u(:,:,z_b,3))/cz1/2;
S6=0.5*((u(:,y_f,:,1)-u(:,y_b,:,1))/cy1/2+(u(x_f,:,:,2)-u(x_b,:,:,2))/cx1/2);
S5=0.5*((u(:,:,z_f,1)-u(:,:,z_b,1))/cz1/2+(u(x_f,:,:,3)-u(x_b,:,:,3))/cx1/2);
S4=0.5*((u(:,:,z_f,2)-u(:,:,z_b,2))/cz1/2+(u(:,y_f,:,3)-u(:,y_b,:,3))/cx1/2);


E0=sum(sum(sum(0.5*c11.*(S1.*S1+S2.*S2+S3.*S3)+c12.*(S1.*S2+S2.*S3+S1.*S3)+2*c44.*(S6.*S6+S5.*S5+S4.*S4))));
tic
for n=1:steps

gradE1=c11.*S1+c12.*(S2+S3);
gradE2=c11.*S2+c12.*(S1+S3);
gradE3=c11.*S3+c12.*(S1+S2);
gradE6=4*c44.*S6; %xy
gradE5=4*c44.*S5; %xz
gradE4=4*c44.*S4; %yz

dS1=-1*factor*(cinv11.*gradE1+cinv12.*gradE2+cinv12.*gradE3);
dS2=-1*factor*(cinv12.*gradE1+cinv11.*gradE2+cinv12.*gradE3);
dS3=-1*factor*(cinv12.*gradE1+cinv12.*gradE2+cinv11.*gradE3);
dS6=-1*factor*cinv44.*gradE6;
dS5=-1*factor*cinv44.*gradE5;
dS4=-1*factor*cinv44.*gradE4;

u(:,y_f,:,1)=u(:,y_f,:,1)+0.5*2*cy1*dS6;
u(:,y_b,:,1)=u(:,y_b,:,1)-0.5*2*cy1*dS6;
u(x_f,:,:,2)=u(x_f,:,:,2)+0.5*2*cx1*dS6;
u(x_b,:,:,2)=u(x_b,:,:,2)-0.5*2*cx1*dS6;

u(x_f,:,:,3)=u(x_f,:,:,3)+0.5*2*cx1*dS5;
u(x_b,:,:,3)=u(x_b,:,:,3)-0.5*2*cx1*dS5;
u(:,:,z_f,1)=u(:,:,z_f,1)+0.5*2*cz1*dS5;
u(:,:,z_b,1)=u(:,:,z_b,1)-0.5*2*cz1*dS5;

u(:,:,z_f,2)=u(:,:,z_f,2)+0.5*2*cz1*dS4;
u(:,:,z_b,2)=u(:,:,z_b,2)-0.5*2*cz1*dS4;
u(:,y_f,:,3)=u(:,y_f,:,3)+0.5*2*cy1*dS4;
u(:,y_b,:,3)=u(:,y_b,:,3)-0.5*2*cy1*dS4;

u(x_f,:,:,1)=u(x_f,:,:,1)+0.5*2*cx1*dS1;
u(x_b,:,:,1)=u(x_b,:,:,1)-0.5*2*cx1*dS1;
u(:,y_f,:,2)=u(:,y_f,:,2)+0.5*2*cy1*dS2;
u(:,y_b,:,2)=u(:,y_b,:,2)-0.5*2*cy1*dS2;
u(:,:,z_f,3)=u(:,:,z_f,3)+0.5*2*cz1*dS3;
u(:,:,z_b,3)=u(:,:,z_b,3)-0.5*2*cz1*dS3;

S1=(u(x_f,:,:,1)-u(x_b,:,:,1))/cx1/2;
S2=(u(:,y_f,:,2)-u(:,y_b,:,2))/cy1/2;
S3=(u(:,:,z_f,3)-u(:,:,z_b,3))/cz1/2;
S6=0.5*((u(:,y_f,:,1)-u(:,y_b,:,1))/cy1/2+(u(x_f,:,:,2)-u(x_b,:,:,2))/cx1/2);
S5=0.5*((u(:,:,z_f,1)-u(:,:,z_b,1))/cz1/2+(u(x_f,:,:,3)-u(x_b,:,:,3))/cx1/2);
S4=0.5*((u(:,:,z_f,2)-u(:,:,z_b,2))/cz1/2+(u(:,y_f,:,3)-u(:,y_b,:,3))/cx1/2);

E(n)=sum(sum(sum(0.5*c11.*(S1.*S1+S2.*S2+S3.*S3)+c12.*(S1.*S2+S2.*S3+S1.*S3)+2*c44.*(S6.*S6+S5.*S5+S4.*S4))));

fprintf('%i energy is %d\n',n,E(n));
end
toc

%iso and biaxial strain
I=S1+S2+S3;
B=sqrt((S1-S2).*(S1-S2)+(S2-S3).*(S2-S3)+(S3-S1).*(S3-S1));

A1=zeros(nz,1);
B1=zeros(nz,1);
C1=zeros(nz,1);
D1=zeros(nz,1);
E1=zeros(nz,1);
for a=1:nz
    A1(a)=I(80,80,a);
    B1(a)=B(80,80,a);
    C1(a)=S1(80,80,a);
    D1(a)=S2(80,80,a);
    E1(a)=S3(80,80,a);
end

figure;
ax1=axes;
hold on
plot (A1); plot(B1);
legend('I','B');
title(['n=',num2str(steps),' f=',num2str(factor), ' Xinout,Zup,layer b/a wetting'])
lh1=line([61 61], get(ax1,'Ylim'), 'LineStyle', '--', 'LineWidth', 0.3, 'Color','k');
% delete(lh1);
%lh1.Color='m'
lh2=line([63 63], get(ax1,'Ylim'), 'LineStyle', '--', 'LineWidth', 0.3, 'Color','k');
lh3=line([86 86], get(ax1,'Ylim'), 'LineStyle', '--', 'LineWidth', 0.3, 'Color','k');

figure;
ax2=axes;
hold on
plot (C1,'b'); plot(D1,'g');plot(E1,'r');
legend('S1','S2','S3');
title(['n=',num2str(steps),' f=',num2str(factor), ' Xinout,Zup,layer b/a wetting'])
lh4=line([61 61], get(ax2,'Ylim'), 'LineStyle', '--', 'LineWidth', 0.3, 'Color','k');
lh5=line([63 63], get(ax2,'Ylim'), 'LineStyle', '--', 'LineWidth', 0.3, 'Color','k');
lh6=line([86 86], get(ax2,'Ylim'), 'LineStyle', '--', 'LineWidth', 0.3, 'Color','k');

