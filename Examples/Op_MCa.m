%MC implementation, matrix vector multiplication

clearvars
format long
format compact	

nPanV = 5:100; %Vector, number of panels
L = length(nPanV);

errLocRegP = zeros(L,1);
errLocRegZ= zeros(L,1);
errProdIntZ = zeros(L,1);
errProdIntP = zeros(L,1);
errGlobRegZ = zeros(L,1);
npvec=zeros(L,1);

ngl = 16; %Number of G-L nodes per panel
if ngl == 16
    [GP,GW] = GaussTW_16(); %Very accurate
else
    [GP,GW] = GaussTW_gen(ngl); %Slightly less accurate
end

a = 0.3;  %Boundary parameter
%Boundary functions
zf = @(p) (1+a*cos(5*p)).*exp(1i*p); 
zpf = @(p) -5*a*sin(5*p).*exp(1i*p) + 1i*(1+a*cos(5*p)).*exp(1i*p);
zppf = @(p) -25*a*cos(5*p).*exp(1i*p) + -2*1i*5*a*sin(5*p).*exp(1i*p) - ...
        (1+a*cos(5*p)).*exp(1i*p);

f = @(z) z.^6+1./z.^6; %Test function
ref = @(z) z.^6-1./z.^6; %Reference solution
for i = 1:L
    nPan = nPanV(i);

    %Initiate boundary
    [z,zp,zpp,w,wzp,awzp,pts,nz,zPan,dp] = zinit(nPan,zf,zpf,zppf,GW,GP,ngl);
    N = nPan*ngl;
    fz = f(z);
    refz = ref(z);

    %Local regularization - p
    CauC = CauC_LocRegP(z,zp,w,wzp,GP,dp,zPan,ngl,N,nPan);
    MC = MC_init(z,wzp,N,CauC);
    errLocRegP(i) = relerr(MC*fz,refz,awzp);

    %Local regularization - z
    CauC = CauC_LocRegZ(z,zp,w,wzp,zPan,ngl,N,nPan);
    MC = MC_init(z,wzp,N,CauC);
    errLocRegZ(i) = relerr(MC*fz,refz,awzp);

    %Product integration - z
    CauC = CauC_ProdIntZ(z,wzp,zPan,ngl,N,nPan);
    MC = MC_init(z,wzp,N,CauC);
    errProdIntZ(i) = relerr(MC*fz,refz,awzp);

    %Product integration - p
    CauC = CauC_ProdIntP(zp,zpp,GP,GW,w,dp,ngl,N,nPan);
    MC = MC_init(z,wzp,N,CauC);
    errProdIntP(i) = relerr(MC*fz,refz,awzp);

    %Global regularization - z
    CauC = CauC_GlobRegZ(z,zp,w,wzp,zPan,ngl,N,nPan);
    MC = MC_init(z,wzp,N,CauC);
    errGlobRegZ(i) = relerr(MC*fz,refz,awzp);

    %Plots
    npvec(i)=N;
    myplot(errLocRegP,errLocRegZ,errProdIntZ,errProdIntP,errGlobRegZ,npvec)
end


function rout=relerr(f,fref,awzp)
%Returns: relative error
  numer=abs(f-fref).^2'*awzp;
  denom=abs(fref).^2'*awzp;
  rout=sqrt(numer/denom);
  rout(rout<eps)=eps;
end

function [z,zp,zpp,w,wzp,awzp,pts,nz,zPan,dp] = zinit(nPan,zf,zpf,zppf,GW,GP,ngl)
    N = ngl*nPan; %Total number of points
    dp = 2*pi/nPan*ones(nPan,1);   %Length of one panel
    pPan = linspace(-pi,pi,nPan+1); %At index k, starting point of panel k in parameter
    pts = zeros(N,1);
    w = zeros(N,1);
    for i = 1:nPan
        pts((i-1)*ngl+1:i*ngl) = (pPan(i)+pPan(i+1))/2+dp(i)/2*GP;
        w((i-1)*ngl+1:i*ngl) = GW*dp(i)/2;
    end
    z = zf(pts);
    zp = zpf(pts);
    zpp = zppf(pts);
    nz=-1i*zp./abs(zp); %Outward unit normal
    zPan = zf(pPan).';
    wzp=w.*zp;
    awzp=w.*abs(zp);
end

function myplot(ErrorLocRegP,ErrorLocRegZ,ErrorProdIntZ,ErrorProdIntP,ErrorGlobRegZ,npvec)
  np=nnz(npvec);
  figure(2)
  x=logspace(log10(80),log10(1700));
  loglog(npvec(1:np),ErrorGlobRegZ(1:np),'go', ...
	 npvec(1:np),ErrorLocRegP(1:np),'r*', ...
     npvec(1:np),ErrorLocRegZ(1:np),'b|', ...
     npvec(1:np),ErrorProdIntP(1:np),'m^',...
         npvec(1:np),ErrorProdIntZ(1:np),'ks')  
  hold on
  loglog(x,5e29*x.^(-15),'k--','Linewidth',1.2)  
  loglog(x,1e31*x.^(-16),'k-.','Linewidth',1.2)  
  legend('GR',...
	 'LR-p','LR-z','PI-p','PI-z','15th ord.','16th ord.')
  grid on
  axis([80 1.6e3 1e-16 1e0])
  title('M_C : matrix-vector multiplication')
  xlabel('number of discretization points')
  ylabel('relative error')
  hold off
  drawnow
end
