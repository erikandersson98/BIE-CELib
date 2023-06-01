%Exterior Helmholz Neumann problem (EHN)

clearvars
format long
format compact	

%Ifleft_Smooth parameters
maxIt = 40;
anglim = pi/8;

nPan = 50; %Number of panels
ngl = 16; %Number of G-L nodes per panel
if ngl == 16
    [GP,GW] = GaussTW_16(); %Very accurate
else
    [GP,GW] = GaussTW_gen(ngl); %Slightly less accurate
end

k = 10; %Wavenumber
a = 0.3; %Boundary parameter
gsc = 0.3+0.5i; %Source point
g = @(z,nz) k*besselh(1,k*abs(gsc-z)).*real(conj(nz).*(gsc-z))./abs(gsc-z);
ref = @(z) besselh(0,k*abs(gsc-z));

%Boundary functions
zf = @(p) (1+a*cos(5*p)).*exp(1i*p); 
zpf = @(p) -5*a*sin(5*p).*exp(1i*p) + 1i*(1+a*cos(5*p)).*exp(1i*p);
zppf = @(p) -25*a*cos(5*p).*exp(1i*p) + -2*1i*5*a*sin(5*p).*exp(1i*p) - ...
        (1+a*cos(5*p)).*exp(1i*p);

%Initiate boundary
[z,zp,zpp,w,wzp,awzp,pts,nz,zPan,pPan,dp] = zinit(nPan,zf,zpf,zppf,GW,GP,ngl);
gz = g(z,nz);
kadz = k*abs(z-z.');
N = nPan*ngl;
ki = 1i*k;
kiadz = 1i*kadz;

%Construct system matrix
CauC = CauC_ProdIntP(zp,zpp,GP,GW,w,dp,ngl,N,nPan);
LogC = LogC_ProdIntP(zp,awzp,GP,dp,ngl,N,nPan);
HypC = HypC_LocRegP(z,zp,zpp,w,wzp,GP,dp,zPan,ngl,N,nPan);
SkiBnd = Sk_init(awzp,N,LogC,ki,kiadz);
KkABnd = KkA_init(z,zp,zpp,w,wzp,nz,N,LogC,kadz);
TkBnd = Tk_init(z,wzp,awzp,nz,N,LogC,HypC,k,kadz);

%Solve iteratively
[mu,~] = GMRESF(-KkABnd - 1i*TkBnd*SkiBnd,-2*gz,N,100,eps);

%Create test points to evaluate U
%Figure out which test points are inside boundary
ngrid = 300; %Square grid side length
gridX = [-1.5,1.5];
gridY = [-1.5,1.5];
xylim = [-1.5,1.5,-1.5,1.5];
u = testGrid(gridX,gridY,ngrid);
ifleft = ifleft_Smooth(zf,zpf,pts,u,maxIt,anglim);

%Field evaluation
U = FieldComp(mu,z,wzp,awzp,nz,zPan,u, ...
		           nPan,k,SkiBnd,ngl,ifleft);
%Plots
Uref = zeros(length(U),1);
Uref(~ifleft) = ref(u(~ifleft));
fieldplot(real(U),u,xylim,ngrid);
errplot(abs(U-Uref),u,xylim,ngrid)



function [z,zp,zpp,w,wzp,awzp,pts,nz,zPan,pPan,dp] = zinit(nPan,zf,zpf,zppf,GW,GP,ngl)
    N = ngl*nPan; %Total number of points
    dp = 2*pi/nPan*ones(nPan,1);   %Length of one panel
    pPan = linspace(-pi,pi,nPan+1); %At index k, starting point of panel k
    pts = zeros(N,1);
    w = zeros(N,1);
    for i = 1:nPan
        pts((i-1)*ngl+1:i*ngl) = (pPan(i)+pPan(i+1))/2+dp(i)/2*GP;
        w((i-1)*ngl+1:i*ngl) = GW*dp(i)/2;
    end
    z = zf(pts);
    zp = zpf(pts);
    zpp = zppf(pts);
    nz=-1i*zp./abs(zp);
    zPan = zf(pPan).';
    wzp=w.*zp;
    awzp=w.*abs(zp);
end

function [U,ifleft] = FieldComp(mu,z,wzp,awzp,nz,zPan,u, ...
				   nPan,k,SkiBnd,ngl,ifleft)
    %Figure out which test points are close to boundary
    disp('Screening starts')
    tic
    dlim = 1.1; %Screening parameter
    closeCell = Screening_dlim(z,awzp,u,nPan,ngl,dlim);
    disp(['Screening and sorting time = ',num2str(toc),' seconds'])

    disp('Regular evaluation starts')
    tic
    U = zeros(length(u),1);
    rightN = sum(~ifleft);
    URight = zeros(rightN,1);
    uRight = u(~ifleft);
    mu2 = SkiBnd*mu;
    for i = 1:rightN
        kadz = k*abs(uRight(i)-z.');
        URight(i) = 1/2*(Sk_Targ(awzp,kadz)*mu+1i*Kk_Targ(uRight(i),z,wzp,kadz)*mu2);
    end
    U(~ifleft) = URight;
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
        awzpsc=awzp(tempind);
        nzsc = nz(tempind);
        musc = mu(tempind);
        closeInds = closeCell{i};
        closeRightInds = closeInds(~ifleft(closeInds));
        ztg = u(closeRightInds);

        %Caculate and add compensation term
        [wCcmp,wLcorr] = wCL_Field(a,b,c,ztg,zsc,wzpsc,awzpsc,nzsc,ifleft(closeRightInds));
        kadz = k*abs(ztg-zsc.');
        U(closeRightInds)=U(closeRightInds) + 1/2*Sk_CmpTarg(awzpsc,kadz,wLcorr)*musc+...
            1i/2*Kk_CmpTarg(ztg,zsc,wzpsc,kadz,wLcorr,wCcmp)*mu2(tempind);
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
  maxfield=max(max(abs(F1)));
  set(0,'DefaultAxesFontSize',12)
  figure()
  imagesc(real(u),imag(u),F1,[-maxfield,maxfield]);      
  colormap(jet)
  axis xy
  colorbar
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
  F1(F1<eps) = eps;
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
  title("Absolute error for EHN, log scale")
end
