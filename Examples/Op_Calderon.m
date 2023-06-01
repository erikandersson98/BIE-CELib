%Implementing Kk, KkA, Sk, Tk, testing Caldéron identities

clearvars
format long
format compact	

nPanV = 5:100; %Vector, number of panels
L = length(nPanV);

err1C1 = zeros(L,1);
err2C1 = zeros(L,1);
err1C2 = zeros(L,1);
err2C2 = zeros(L,1);
err1C3 = zeros(L,1);
err2C3 = zeros(L,1);
err1C4 = zeros(L,1);
err2C4 = zeros(L,1);
npvec=zeros(L,1);

ngl = 16; %Number of G-L nodes per panel
if ngl == 16
    [GP,GW] = GaussTW_16(); %Very accurate
else
    [GP,GW] = GaussTW_gen(ngl); %Slightly less accurate
end

k = 10; %Wavenumber
a = 0.3; %Boundary parameter
%Boundary functions
zf = @(p) (1+a*cos(5*p)).*exp(1i*p); 
zpf = @(p) -5*a*sin(5*p).*exp(1i*p) + 1i*(1+a*cos(5*p)).*exp(1i*p);
zppf = @(p) -25*a*cos(5*p).*exp(1i*p) + -2*1i*5*a*sin(5*p).*exp(1i*p) - ...
        (1+a*cos(5*p)).*exp(1i*p);

f = @(z) z.^6+1./z.^6; %Test function
for i = 1:L
    nPan = nPanV(i);
    
    %Initiate boundary
    [z,zp,zpp,w,wzp,awzp,pts,nz,zPan,pPan,dp] = zinit(nPan,zf,zpf,zppf,GW,GP,ngl);
    kadz = k*abs(z-z.');
    N = nPan*ngl;
    fz = f(z);

    for iCase = 1:4
        if iCase == 1
            %CauC = CauC_ProdIntP(zp,zpp,GP,GW,w,dp,ngl,N,nPan);
            LogC = LogC_ProdIntP(zp,awzp,GP,dp,ngl,N,nPan);
            HypC = HypC_ProdIntZ(z,wzp,zPan,ngl,N,nPan);
        elseif iCase == 2
            %CauC = CauC_ProdIntZ(z,wzp,zPan,ngl,N,nPan);
            LogC = LogC_ProdIntZ(z,awzp,nz,zPan,ngl,N,nPan);
            HypC = HypC_ProdIntZ(z,wzp,zPan,ngl,N,nPan);
        elseif iCase == 3
            %CauC = CauC_ProdIntP(zp,zpp,GP,GW,w,dp,ngl,N,nPan);
            LogC = LogC_ProdIntP(zp,awzp,GP,dp,ngl,N,nPan);
            HypC = HypC_LocRegP(z,zp,zpp,w,wzp,GP,dp,zPan,ngl,N,nPan);
        elseif iCase == 4
            %CauC = CauC_ProdIntZ(z,wzp,zPan,ngl,N,nPan);
            LogC = LogC_ProdIntP(zp,awzp,GP,dp,ngl,N,nPan);
            HypC = HypC_LocRegZ(z,zp,w,wzp,zPan,ngl,N,nPan);
        end
        SkBnd = Sk_init(awzp,N,LogC,k,kadz);
        KkBnd = Kk_init(z,zp,zpp,w,wzp,N,LogC,kadz);
        KkABnd = KkA_init(z,zp,zpp,w,wzp,nz,N,LogC,kadz);
        TkBnd = Tk_init(z,wzp,awzp,nz,N,LogC,HypC,k,kadz);

        if iCase == 1
            err1C1(i) = relerr((KkBnd*KkBnd-SkBnd*TkBnd)*fz,fz,awzp);
            err2C1(i) = relerr((KkABnd*KkABnd-TkBnd*SkBnd)*fz,fz,awzp);
        elseif iCase == 2
            err1C2(i) = relerr((KkBnd*KkBnd-SkBnd*TkBnd)*fz,fz,awzp);
            err2C2(i) = relerr((KkABnd*KkABnd-TkBnd*SkBnd)*fz,fz,awzp);
        elseif iCase == 3
            err1C3(i) = relerr((KkBnd*KkBnd-SkBnd*TkBnd)*fz,fz,awzp);
            err2C3(i) = relerr((KkABnd*KkABnd-TkBnd*SkBnd)*fz,fz,awzp);
        elseif iCase == 4
            err1C4(i) = relerr((KkBnd*KkBnd-SkBnd*TkBnd)*fz,fz,awzp);
            err2C4(i) = relerr((KkABnd*KkABnd-TkBnd*SkBnd)*fz,fz,awzp);
        end
        %Plots
        npvec(i)=N;
        myplot(err1C1,err1C2,err1C3,err1C4,npvec,1)
        myplot(err2C2,err2C2,err2C3,err2C4,npvec,2)
    end
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

function myplot(ErrorC1,ErrorC2,ErrorC3,ErrorC4,npvec,fig)
    np=nnz(npvec);
    figure(fig)
    x=logspace(log10(80),log10(1700));
    loglog(npvec(1:np),ErrorC1(1:np),'o', ...
	    npvec(1:np),ErrorC2(1:np),'*', ...
        npvec(1:np),ErrorC3(1:np),'sk',...
        npvec(1:np),ErrorC4(1:np),'gx')  
    hold on
    axis([80 1.6e3 1e-14 1e0])
    if fig == 1
        title('K_k*K_k-S_k*T_k : matrix-vector multiplication')
        loglog(x,9e26*x.^(-14),'k--','Linewidth',1.2)  
        loglog(x,3e38*x.^(-15),'k-.','Linewidth',1.2) 
        legend('C1',...
	        'C2','C3','C4','14th order','15th order','Location','southwest')
    elseif fig == 2
        title('K_k^A*K_k^A-T_k*S_k : matrix-vector multiplication')
        loglog(x,7e30*x.^(-14),'k--','Linewidth',1.2)  
        loglog(x,1e36*x.^(-14),'k-.','Linewidth',1.2) 
        legend('C1',...
	        'C2','C3','C4','14th order','14th order','Location','southwest')
    end
    xlabel('number of discretization points')
    ylabel('relative error')
    hold off
    drawnow
end
