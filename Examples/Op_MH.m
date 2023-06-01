%MH implementation, testing matrix-vector multiplication

clearvars
format long
format compact	

nPanV = 5:100; %Vector, number of panels
L = length(nPanV);

errLocRegP = zeros(L,1);
errLocRegZ = zeros(L,1);
errProdIntZ = zeros(L,1);
npvec=zeros(L,1);

ngl = 16; %Number of G-L nodes per panel
if ngl == 16
    [GP,GW] = GaussTW_16; %Very accurate
     
else
    [GP,GW]=GaussTW_gen(ngl); %Slightly less accurate
end

a = 0.3; %Boundary parameter
%Boundary functions
zf = @(p) (1+a*cos(5*p)).*exp(1i*p);
zpf = @(p) -5*a*sin(5*p).*exp(1i*p) + 1i*(1+a*cos(5*p)).*exp(1i*p);
zppf = @(p) -25*a*cos(5*p).*exp(1i*p) -10*1i*a*sin(5*p).*exp(1i*p) - ... 
        (1+a*cos(5*p)).*exp(1i*p);

f = @(z) z.^6+1./z.^6; %Test function
ref = @(z) 6*z.^5+6./z.^7; %Reference function
for i = 1:L
    nPan = nPanV(i);

    %Initiate boundary
    [z,zp,zpp,w,wzp,awzp,pts,nz,zPan,pPan,dp] = zinit(nPan,zf,zpf,zppf,GW,GP,ngl);
    N = nPan*ngl;
    fz = f(z);
    refz = ref(z);

    %Local regularization - p
    HypC = HypC_LocRegP(z,zp,zpp,w,wzp,GP,dp,zPan,ngl,N,nPan);
    MH = MH_init(z,wzp,N,HypC);
    errLocRegP(i) = relerr(MH*fz,refz,awzp);

    %Local regularization - z
    HypC = HypC_LocRegZ(z,zp,w,wzp,zPan,ngl,N,nPan);
    MH = MH_init(z,wzp,N,HypC);
    errLocRegZ(i) = relerr(MH*fz,refz,awzp);

    %Product integration - z
    HypC = HypC_ProdIntZ(z,wzp,zPan,ngl,N,nPan);
    MH = MH_init(z,wzp,N,HypC);
    errProdIntZ(i) = relerr(MH*fz,refz,awzp);

    npvec(i)=N;
    myplot(errProdIntZ,errLocRegP,errLocRegZ,npvec);
end


function rout=relerr(f,fref,awzp)
%Returns: rout - Relative error
    numer=abs(f-fref).^2'*awzp;
    denom=abs(fref).^2'*awzp;
    rout=sqrt(numer/denom);
    rout(rout<eps)=eps;
end

function [z,zp,zpp,w,wzp,awzp,pts,nz,zPan,pPan,dp] = zinit(nPan,zf,zpf,zppf,GW,GP,ngl)
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
    zPan = zf(pPan).';
    wzp=w.*zp;
    awzp=w.*abs(zp);
end

function myplot(ErrorProdIntZ,ErrorLocRegP,ErrorLocRegZ,npvec)
    np=nnz(npvec);
    figure(1)
    x=logspace(log10(80),log10(1600));
    loglog(npvec(1:np),ErrorLocRegP(1:np),'r*',...
        npvec(1:np),ErrorLocRegZ(1:np),'b|',...
        npvec(1:np),ErrorProdIntZ(1:np),'sk')  
    hold on
    loglog(x,4e30*x.^(-15),'k--','Linewidth',1.2)  
    loglog(x,7e31*x.^(-16),'k-.','Linewidth',1.2)  
    legend('LR-p','LR-z','PI-z','15th order','16th order')
    grid on
    axis([80 1.6e3 1e-16 1e0])
    title('M_H : matrix-vector multiplication')
    xlabel('number of discretization points')
    ylabel('relative error')
    hold off
    drawnow
end
