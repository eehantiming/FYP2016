%piezo
eps0=8.854187e-12;
epsR_i=15.15; epsR_g=12.9; %static dielectric constant
e14_i=-4.5e-2; e14_g=-0.16; %Cm-2
steps=100; factor=1;

e14=zeros(nx,ny,nz);
epsR=zeros(nx,ny,nz);
GG=zeros(steps,1);

V=zeros(nx,ny,nz);

for x=1:nx
    for y=1:ny
        for z=1:nz
            if m(x,y,z)==1
                e14(x,y,z)=e14_g;
                epsR(x,y,z)=epsR_g;
                %V(x,y,z)=30*1.6e-22;
            else
                e14(x,y,z)=e14_i;
                epsR(x,y,z)=epsR_i;
            end
        end
    end
end

rho=-1*10000*e14.*((S4(x_f,:,:)-S4(x_b,:,:))/cx1/2+(S5(:,y_f,:)-S5(:,y_b,:))/cy1/2+(S6(:,:,z_f)-S6(:,:,z_b))/cz1/2);
FV=-2*(1/cx1^2+1/cy1^2+1/cz1^2);
tic
for n=1:steps
    lapV=(V(x_f,:,:)-2*V+V(x_b,:,:))/(cx1*2)+(V(:,y_f,:)-2*V+V(:,y_b,:))/(cy1*2)+(V(:,:,z_f)-2*V+V(:,:,z_b))/(cz1*2);
    F=eps0*epsR.*lapV-rho;
    G=sum(sum(sum(F)));
    GG(n)=G;
    cV=F/FV;
    V=V-factor*cV;
    fprintf('%i F is %d\n',n,G);
end
toc
figure
plot(GG);

F1=zeros(nx,ny);
F2=zeros(nx,ny);F3=zeros(nx,ny);F4=zeros(nx,ny);F5=zeros(nx,ny);
material=zeros(nx,ny);
for x=1:nx
    for y=1:ny
F1(x,y)=V(x,y,63);
F2(x,y)=V(x,y,123);
F3(x,y)=V(x,y,82);
F4(x,y)=V(x,y,142);
F5(x,y)=V(80,x,y);
material(x,y)=m(80,x,y);
    end
end
figure
surf(F1);title('z=base of lower pyramid')
figure
surf(F2);title('z=base of upper pyramid')
figure
surf(F3);title('z=top of lower pyramid')
figure
surf(F4);title('z=top of upper pyramid')
figure
surf(F5);title('y-z plane through pyramid')
figure
surf(material);