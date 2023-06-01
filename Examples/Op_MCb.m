%Solving system of equations MC*MC*f = g

clearvars
format long
format compact	

nPanV = 5:100; %Vector, number of panels
L = length(nPanV);

errLocRegP = zeros(L,1);
errLocRegZ = zeros(L,1);
errProdIntP = zeros(L,1);
errProdIntZ = zeros(L,1);
errGlobRegZ = zeros(L,1);
npvec=zeros(L,1);

ngl = 16;
if ngl == 16
    [GP,GW] = GaussTW_16(); %Very accurate
else
    [GP,GW] = GaussTW_gen(ngl); %Slightly less accurate
end

a = 0.3; %Boundary parameter
zf = @(p) (1+a*cos(5*p)).*exp(1i*p); %Boundary functions
zpf = @(p) -5*a*sin(5*p).*exp(1i*p) + 1i*(1+a*cos(5*p)).*exp(1i*p);
zppf = @(p) -25*a*cos(5*p).*exp(1i*p) + -2*1i*5*a*sin(5*p).*exp(1i*p) - ...
        (1+a*cos(5*p)).*exp(1i*p);

f = @(z) z.^6+1./z.^6; %Test function
for i = 1:L
    nPan = nPanV(i);

    %Initiate boundary
    [z,zp,zpp,w,wzp,awzp,pts,nz,zPan,dp] = zinit(nPan,zf,zpf,zppf,GW,GP,ngl);
    N = nPan*ngl;
    fz = f(z);
    
    %Local regularization - p
    CauC = CauC_LocRegP(z,zp,w,wzp,GP,dp,zPan,ngl,N,nPan);
    MC = MC_init(z,wzp,N,CauC);
    [mu,~] = GMRESF(MC*MC-eye(N),fz,N,100,eps);
    errLocRegP(i) = relerr(mu,fz,awzp);

    %Local regularization - z
    CauC = CauC_LocRegZ(z,zp,w,wzp,zPan,ngl,N,nPan);
    MC = MC_init(z,wzp,N,CauC);
    [mu,~] = GMRESF(MC*MC-eye(N),fz,N,100,eps);
    errLocRegZ(i) = relerr(mu,fz,awzp);
    
    %Product integration - z
    CauC = CauC_ProdIntZ(z,wzp,zPan,ngl,N,nPan);
    MC = MC_init(z,wzp,N,CauC);
    [mu,~] = GMRESF(MC*MC-eye(N),fz,N,100,eps);
    errProdIntZ(i) = relerr(mu,fz,awzp);
    
    %Product integration - p
    CauC = CauC_ProdIntP(zp,zpp,GP,GW,w,dp,ngl,N,nPan);
    MC = MC_init(z,wzp,N,CauC);
    [mu,~] = GMRESF(MC*MC-eye(N),fz,N,100,eps);
    errProdIntP(i) = relerr(mu,fz,awzp);

    %Global regularization
    CauC = CauC_GlobRegZ(z,zp,w,wzp,zPan,ngl,N,nPan);
    MC = MC_init(z,wzp,N,CauC);
    [mu,~] = GMRESF(MC*MC-eye(N),fz,N,100,eps);
    errGlobRegZ(i) = relerr(mu,fz,awzp);

    npvec(i)=N;
    myplot(errLocRegP,errLocRegZ,errProdIntP,errProdIntZ,errGlobRegZ,npvec)
end


function rout=relerr(f,fref,awzp)
%Returns: rout - relative error
  numer=abs(f-fref).^2'*awzp;
  denom=abs(fref).^2'*awzp;
  rout=sqrt(numer/denom);
  rout(rout<eps)=eps;
end

function [z,zp,zpp,w,wzp,awzp,pts,nz,zPan,dp] = zinit(nPan,zf,zpf,zppf,GW,GP,ngl)
    N = ngl*nPan; %Total number of points
    dp = 2*pi/nPan*ones(nPan,1);   %Length of one panel
    pPan = linspace(-pi,pi,nPan+1); %Starting point of panels
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
    zPan = zf(pPan(1:end)).';
    wzp=w.*zp;
    awzp=w.*abs(zp);
end

function myplot(ErrorLocRegP,ErrorLocRegZ,ErrorProdIntP,ErrorProdIntZ,ErrorGlobRegZ,npvec)
  np=nnz(npvec);
  figure(1)
  x=logspace(log10(80),log10(1600));
  loglog(npvec(1:np),ErrorGlobRegZ(1:np),'go', ...
	 npvec(1:np),ErrorLocRegP(1:np),'r*', ...
     npvec(1:np),ErrorLocRegZ(1:np),'b|', ...
     npvec(1:np),ErrorProdIntP(1:np),'m^',...
         npvec(1:np),ErrorProdIntZ(1:np),'ks')   
  hold on
  loglog(x,3e31*x.^(-15),'k--','Linewidth',1.2)  
  legend('GR',...
	 'LR-p','LR-z','PI-p','PI-z','15th order', 'Location','southwest')
  grid on
  axis([80 1.6e3 1e-16 1])
  title('M_C^2 : Solving system of equations')
  xlabel('number of discretization points')
  ylabel('relative error')
  hold off
  drawnow
end

