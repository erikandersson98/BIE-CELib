%Example-Laplace Transmission problem (LT)

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
    [GP,GW] = GaussTW_gen(ngl); %Slightly less accurate
end
%Field constants
sigma1 = 1;
sigma2 = 100;
lambda = (sigma2-sigma1)/(sigma2+sigma1);

a = 0.3; %Boundary parameter
ze = 1; 
e = ze/abs(ze); %Unit field
gf = @(z) -2*lambda*real(conj(e)*z); %RHS
%Boundary functions
zf = @(p) (1+a*cos(5*p)).*exp(1i*p); 
zpf = @(p) -5*a*sin(5*p).*exp(1i*p) + 1i*(1+a*cos(5*p)).*exp(1i*p);
zppf = @(p) -25*a*cos(5*p).*exp(1i*p) + -2*1i*5*a*sin(5*p).*exp(1i*p) - ...
        (1+a*cos(5*p)).*exp(1i*p);

for i = 1:2
    if i == 1
        disp('First solve')
    elseif i == 2
        disp('Second solve')
        nPan = 1.5*nPan;
    end
    %Initiate boundary
    [z,zp,zpp,w,wzp,awzp,pts,nz,zPan] = zinit(nPan,zf,zpf,zppf,GW,GP,ngl);
    g = gf(nz); %RHS at z 
    N = nPan*ngl; %Total number of nodes
    
    %System matrix
    SysMat = lambda*NPA_init(z,zp,zpp,w,awzp,N,nz);
    
    %Solve iteratively
    [mu,nIt] = GMRESF(SysMat,g,N,100,eps);
    
    %Create test points to evaluate U
    %Figure out which test points are inside boundary
    ngrid = 300; %Square grid side length
    gridX = [-1.5,1.5];
    gridY = [-1.5,1.5];
    xylim = [-1.5,1.5,-1.5,1.5];
    u = testGrid(gridX,gridY,ngrid);
    ifleft = ifleft_Smooth(zf,zpf,pts,u,maxIt,anglim);

    %Field evaluation
    [U,ifleft] = FieldComp(mu,z,awzp,nz,zPan,u, ...
				       nPan,e,ngl,ifleft);
    %Plots
    if i == 1
        U1 = U;
        fieldplot(U,u,xylim,ngrid);
    elseif i == 2
        errplot(U-U1,u,xylim,ngrid)
    end
end

function [z,zp,zpp,w,wzp,awzp,pts,nz,zPan] = zinit(nPan,zf,zpf,zppf,GW,GP,ngl)
    N = ngl*nPan; %Total number of points
    dp = 2*pi/nPan*ones(nPan,1);   %Length of panel in parameter
    panStart = linspace(-pi,pi,nPan+1); %Starting point of panels
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
    zPan = zf(panStart).'; %Panel Startpoints in complex plane
    wzp=w.*zp; 
    awzp=w.*abs(zp);
end

function [U,ifleft] = FieldComp(mu,z,awzp,nz,zPan,u, ...
				   nPan,e,ngl,ifleft)
    %Figure out which test points are close to boundary
    disp('Screening starts')
    tic
    dlim = 1.1;
    closeCell = Screening_dlim(z,awzp,u,nPan,ngl,dlim);
    disp(['Screening and sorting time = ',num2str(toc),' seconds'])

    %Regular eval
    disp('Regular evaluation starts')
    tic
    U = real(conj(e)*u); %Term from underlying unit field
    for i = 1:length(u)
        U(i) = U(i) - 1/(2*pi)*ML_Targ(u(i),z,awzp)*mu; %Standard G-L
    end
    disp(['Regular evaluation time    = ', num2str(toc), ' seconds'])

    %Close eval
    disp('Special evaluation starts')
    tic
    c=(1-(-1).^(1:ngl))./(1:ngl);
    for i = 1:nPan
        a=zPan(i);
        b=zPan(i+1);
        tempind=(i-1)*ngl+1:i*ngl;
        zsc=z(tempind); %Source points
        awzpsc=awzp(tempind);
        nzsc = nz(tempind);
        musc = mu(tempind);
        closeInds = closeCell{i}; %Target points close to panel i have value 1
        uClose = u(closeInds); %Close target points
        ifleftClose = ifleft(closeInds);

        %Caculate and add compensation term
        wLcorr = wLcorr_Field(a,b,c,uClose,zsc,awzpsc,ifleftClose,nzsc);
        U(closeInds)=U(closeInds)-1/(2*pi)*ML_CmpTarg(awzpsc,wLcorr)*musc;
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

function fieldplot(U,u,xylim,Ng)
    F1 = zeros(Ng);
    F1(:)=U;
    maxfield=max(max(F1));
    minfield = min(min(F1));
    set(0,'DefaultAxesFontSize',12)
    figure()
    imagesc(real(u),imag(u),F1,[minfield,maxfield]);      
    colormap(jet)
    axis xy
    colorbar
    xlabel('$x_1$','Interpreter','LaTeX','FontSize',17)      
    ylabel('$x_2$','Interpreter','LaTeX','FontSize',17)      
    title("Field evaluation, LT problem")
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
  title("Absolute error for LT, log scale")
end
