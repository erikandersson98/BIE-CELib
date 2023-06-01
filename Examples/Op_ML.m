%ML implementation, testing matrix-vector multiplication

clearvars
format long
format compact	

nPanV = 5:130; %Vector, number of panels
L = length(nPanV);

saveProdIntZ = zeros(L,1);
saveProdIntP = zeros(L,1);
npvec=zeros(L,1);

ngl = 16; %Number of G-L nodes per panel
if ngl == 16
    [GP,GW] = GaussTW_16(); %Very accurate
else
    [GP,GW]=GaussTW_gen(ngl); %Slightly less accurate
end

a = 0.3; %Boundary parameter
%Boundary functions
zf = @(p) (1+a*cos(5*p)).*exp(1i*p); 
zpf = @(p) -5*a*sin(5*p).*exp(1i*p) + 1i*(1+a*cos(5*p)).*exp(1i*p);
zppf = @(p) -25*a*cos(5*p).*exp(1i*p) + -2*1i*5*a*sin(5*p).*exp(1i*p) - ...
        (1+a*cos(5*p)).*exp(1i*p);

f = @(z) abs(z.^6+1./z.^6); %Test function
for i = 1:L
    nPan = nPanV(i);

    %Initiate boundary
    [z,zp,zpp,w,wzp,awzp,pts,nz,zPan,pPan,dp] = zinit(nPan,zf,zpf,zppf,GW,GP,ngl);
    N = nPan*ngl;
    fz = f(z);

    %Product integration - z
    LogC = LogC_ProdIntZ(z,awzp,nz,zPan,ngl,N,nPan);
    ML = ML_init(z,awzp,N,LogC);
    saveProdIntZ(i) = abserr(ML*fz,awzp);

    %Product integration - p
    LogC = LogC_ProdIntP(zp,awzp,GP,dp,ngl,N,nPan);
    ML = ML_init(z,awzp,N,LogC);
    saveProdIntP(i) = abserr(ML*fz,awzp);
    npvec(i)=N;
end
%Plots
plot1 = abs(saveProdIntZ-saveProdIntZ(80))/saveProdIntZ(80);
plot2 = abs(saveProdIntP-saveProdIntP(80))/saveProdIntP(80);
myplot(plot1,plot2,npvec);


function rout=abserr(f,awzp)
%Returns: rout - Absolute error
    rout=sqrt(abs(f).^2'*awzp);
    rout(rout<eps)=eps;
end

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

function myplot(ErrorProdIntZ,ErrorProdIntP,npvec)
    np=nnz(npvec);
    figure(1)
    x=logspace(2,log10(2080));
    loglog(npvec(1:np),ErrorProdIntP(1:np),'m^',...
        npvec(1:np),ErrorProdIntZ(1:np),'sk')
    hold on
    loglog(x,5e30*x.^(-14),'k--','Linewidth',1.2)
    legend('PI-p','PI-z','14th ord.')
    grid on
    axis([80 2.08e3 1e-16 1e0])
    title('M_L : matrix-vector multiplication')
    xlabel('number of discretization points')
    ylabel('relative error')
    hold off
    drawnow
end
