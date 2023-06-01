%Example-Interior Laplace Dirichlet problem (ILD)

clearvars
format long
format compact	

%Ifleft_Smooth parameters
maxIt = 40;
anglim = pi/8;

nPan = 100; %Number of panels
ngl = 16; %Number of G-L nodes per panel
if ngl == 16
    [GP,GW] = GaussTW_16(); %Very accurate
else
    [GP,GW]=GaussTW_gen(ngl); %Slightly less accurate
end

a = 0.3; %Boundary parameter
gf = @(z) real(z); %Dirichlet condition
%Boundary functions
zf = @(p) (1+a*cos(5*p)).*exp(1i*p); 
zpf = @(p) -5*a*sin(5*p).*exp(1i*p) + 1i*(1+a*cos(5*p)).*exp(1i*p);
zppf = @(p) -25*a*cos(5*p).*exp(1i*p) + -2*1i*5*a*sin(5*p).*exp(1i*p) - ...
        (1+a*cos(5*p)).*exp(1i*p);

%Initiate boundary
[z,zp,zpp,w,wzp,awzp,pts,nz,zPan] = zinit(nPan,zf,zpf,zppf,GW,GP,ngl);
g = gf(z);
N = nPan*ngl;

%System matrix
NP = NP_init(z,zp,zpp,w,wzp,N);

%Solve iteratively
[mu,nIt] = GMRESF(-NP,g,N,100,eps);

%Create test points to evaluate U
%Figure out which test points are inside boundary
ngrid = 300; %Square grid side length
gridX = [-1.5,1.5];
gridY = [-1.5,1.5];
xylim = [-1.5,1.5,-1.5,1.5];
u = testGrid(gridX,gridY,ngrid);
ifleft = ifleft_Smooth(zf,zpf,pts,u,maxIt,anglim);

%Field evaluation
U = FieldComp(mu,z,wzp,awzp,zPan,u,nPan,ngl,ifleft);

%Plots
fieldplot(U,u,xylim,z,ngrid)  
Uref=real(u);
Uref(~ifleft)=0;
errplot(Uref-U,u,xylim,ngrid)


function [z,zp,zpp,w,wzp,awzp,pts,nz,zPan] = zinit(nPan,zf,zpf,zppf,GW,GP,ngl)
    N = ngl*nPan;   %Total number of points
    dp = 2*pi/nPan*ones(nPan,1);   %Length of panel in parameter
    panStart = linspace(-pi,pi,nPan+1);     %Starting point of panels
    pts = zeros(N,1); %Scaled Gauss points
    w = zeros(N,1); %Scaled Gauss weights
    for i = 1:nPan
        pts((i-1)*ngl+1:i*ngl) = (panStart(i)+panStart(i+1))/2+dp(i)/2*GP;
        w((i-1)*ngl+1:i*ngl) = GW*dp(i)/2;
    end
    z = zf(pts);
    zp = zpf(pts);
    zpp = zppf(pts);
    nz=-1i*zp./abs(zp); %Outwards unit normal at z
    zPan = zf(panStart).'; %Panel start-points in complex plane
    wzp=w.*zp; 
    awzp=w.*abs(zp);
end

function U = FieldComp(mu,z,wzp,awzp,zPan,u, ...
				   nPan,ngl,ifleft)
    %Figure out which test points are close to boundary
    disp('Screening starts')
    tic
    dlim = 1.1; %Screening parameter
    closeCell = Screening_dlim(z,awzp,u,nPan,ngl,dlim);
    uLeft = u(ifleft);
    disp(['Screening and sorting time = ',num2str(toc),' seconds'])

    disp('Regular evaluation starts')
    tic
    U = zeros(length(u),1);
    ULeft = zeros(length(uLeft),1);
    for i = 1:length(uLeft)
        ULeft(i) = real(MC_Targ(uLeft(i),z,wzp))*mu; %Standard G-L
    end
    U(ifleft) = ULeft;
    disp(['Regular evaluation time    = ', num2str(toc), ' seconds'])

    disp('Special evaluation starts')
    tic
    c=(1-(-1).^(1:ngl))./(1:ngl);
    for i = 1:nPan
        a=zPan(i);
        b=zPan(i+1);
        tempind=(i-1)*ngl+1:i*ngl;
        zsc=z(tempind);
        wzpsc=wzp(tempind);
        closeInds = closeCell{i};
        closeLeftInds = closeInds(ifleft(closeInds));

        %Caculate and add compensation term
        wCcmp=wCcmp_Field(a,b,c,u(closeLeftInds),zsc,wzpsc,ifleft(closeLeftInds));
        U(closeLeftInds)=U(closeLeftInds)+real(MC_CmpTarg(wCcmp))*mu(tempind);
    end
    disp(['Special evaluation time    = ', num2str(toc), ' seconds'])
end

function u = testGrid(gridX,gridY,ngrid)
    linX = linspace(gridX(1),gridX(2),ngrid);
    linY = linspace(gridY(1),gridY(2),ngrid);
    [X,Y] = meshgrid(linX,linY);
    u = X + 1i*Y;
    u = u(:);
end

function fieldplot(U,u,xylim,z,Ng)
    F1 = zeros(Ng);
    F1(:)=U;
    maxafield=max(max(abs(F1)));
    set(0,'DefaultAxesFontSize',12)
    figure()
    imagesc(real(u),imag(u),F1,[-maxafield,maxafield]);      
    colormap(jet)
    axis xy
    colorbar
    hold on
    ztot=[-1.5;z;z(1);-1.5;-1.5-1.5i;1.5-1.5i;1.5+1.5i;-1.5+1.5i;-1.5];
    fill(real(ztot),imag(ztot),'w','EdgeColor','w');
    xlabel('$x_1$','Interpreter','LaTeX','FontSize',17)      
    ylabel('$x_2$','Interpreter','LaTeX','FontSize',17)      
    axis square
    axis(xylim)
    text(-2.05,-1.85,'(a)','FontSize',16)      
    drawnow
end

function errplot(Udiff,u,xylim,Ng)
    F1 = zeros(Ng);
    F1(:)=abs(Udiff);
    F1=log10(F1);
    maxa_err=max(max(F1));
    disp(['max field error = ',num2str(maxa_err)])
    set(0,'DefaultAxesFontSize',12)
    figure()
    imagesc(real(u),imag(u),F1,[log10(eps) maxa_err]);      
    colormap(flipud(pink(256)))
    axis xy
    hold on
    colorbar
    xlabel('$x_1$','Interpreter','LaTeX','FontSize',17)      
    ylabel('$x_2$','Interpreter','LaTeX','FontSize',17)
    text(-2.05,-1.85,'(b)','FontSize',16)      
    axis square
    axis(xylim)
    drawnow
    title("Absolute error for ILD, log scale")
end