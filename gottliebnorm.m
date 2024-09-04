function [g, dgdr] = gottliebnorm(mu, re, xin, c, s, nax, max)
%%% Initalize Norm Sizes for Speed
size=nax+1;
norm1(1,1:size)=0;
norm11(1,1:size)=0;
norm2(1,1:size)=0;
normn10(1,1:size)=0;
normn20(1,1:size)=0;
norm1m(1:size,1:size)=0;
norm2m(1:size,1:size)=0;
normn1(1:size,1:size)=0;
normn2(1:size,1:size)=0;

%%% function to return both gradient and the jacobian ...
for n = 2:size %RAE
    norm1(n) = sqrt((2*n+1)/(2*n-1)); %RAE   - norm in (3.1)
    norm2(n) = sqrt((2*n+1)/(2*n-3)); %RAE   - norm in (3.2)
    norm11(n) = sqrt((2*n+1)/(2*n))/(2*n-1); %RAE - norm in (3.3)
    normn10(n) = sqrt((n+1)*n/2); %RAE            - norm in (3.14)/ Gott 7.19
    normn20(n) = sqrt(0.5*n*(n-1)*(n+1)*(n+2)); %MM - norm in Gott 7.22
    for m = 1:n %RAE
        norm1m(n,m) = sqrt((n-m)*(2*n+1)/((n+m)*(2*n-1))); %RAE       (3.4)
        norm2m(n,m) = sqrt((n-m)*(n-m-1)*(2*n+1)/((n+m)*(n+m-1)*(2*n-3))); %RAE   (3.5)
        normn1(n,m) = sqrt((n+m+1)*(n-m)); %RAE            Gott 7.19
        normn2(n,m) = sqrt((n-m)*(n-m-1)*(n+m+1)*(n+m+2)); %MM - Gott 7.22
    end %RAE
end %RAE
x=xin; %RAE
r = sqrt(x(1)^2+x(2)^2+x(3)^2);
ri=1/r;
xor=x(1)*ri;
yor=x(2)*ri;
zor=x(3)*ri;
ep=zor;
reor=re*ri;
reorn=reor;
muor2=mu*ri*ri;
p(1,1) = 1; %RAE
p(1,2) = 0; %RAE
p(1,3) = 0; %RAE
p(2,2) = sqrt(3); %RAE %norm
p(2,3) = 0; %RAE
p(2,4) = 0; %RAE
for n = 2:nax %RAE
    ni = n+1; %RAE
    p(ni,ni) = norm11(n)*p(n,n)*(2*n-1); %RAE %norm     - norm in (3.15)
    p(ni,ni+1) = 0; %RAE
    p(ni,ni+2) = 0; %RAE
end
ctil(1)=1; %RAE
stil(1)=0; %RAE
ctil(2)=xor; %RAE      % c1 = xor, s1 = yor
stil(2)=yor; %RAE      % Eq 3-18 recursions of normalized ctil,stil
sumh=0;    % summation for H
sumgm=1;   % summation for Gamma
sumj=0;    % summation for J
sumk=0;    % summation for K

%%% new terms needed for second order partials
sumL =2; sumM = 0; sumN = 0; sumO =0;
sumP = 0; sumQ = 0; sumR = 0; sumS = 0; sumT = 0;
p(2,1) = sqrt(3)*ep; %RAE %norm
for n=2:nax
    ni=n+1; %RAE
    reorn=reorn*reor;
    n2m1=n+n-1;
    nm1=n-1;
    np1=n+1; nm2 = n-2; np2 = n+2;
    p(ni,n) = normn1(n,n-1)*ep*p(ni,ni); %RAE %norm
    p(ni,1) = (n2m1*ep*norm1(n)*p(n,1)-nm1*norm2(n)*p(nm1,1))/n; %RAE %norm   Eq. 3.9/3.15
    p(ni,2) = (n2m1*ep*norm1m(n,1)*p(n,2)-n*norm2m(n,1)*p(nm1,2))/(nm1); %RAE %norm Eq 3.9/3.15

    if(n == 2)
        p(ni,3) = p(3,3);    % already defined in line 37.
    else
        p(ni,3) = (n2m1*ep*norm1m(n,2)*p(n,3) - np1*norm2m(n,2)*p(nm1,3))/(nm2); %MM norm Eq 3.15
    end
    sumhn=normn10(n)*p(ni,2)*c(ni,1); %norm %RAE
    sumgmn=p(ni,1)*c(ni,1)*np1; %RAE


    sumLn = p(ni,1)*c(ni,1)*np1*np2;
    sumMn = normn20(n)*p(ni,3)*c(ni,1);                %MM - norm Gott 7.22
    sumPn = normn10(n)*p(ni,2)*c(ni,1)*np1;
    if (max>0)
        for m = 2:n-2
            mi = m+1; %RAE
            p(ni,mi) = (n2m1*ep*norm1m(n,m)*p(n,mi)-...
                (nm1+m)*norm2m(n,m)*p(nm1,mi))/(n-m); %RAE %norm - Eq 3-14 in report)
        end
        sumjn=0;
        sumkn=0;
        % MM
        sumNn = 0;
        sumOn = 0; sumSn = 0;
        sumQn = 0; sumRn = 0; sumTn = 0;
        ctil(ni)=ctil(2)*ctil(ni-1)-stil(2)*stil(ni-1); %RAE  % recursion in 3-18a, ctil(1) = c_0
        stil(ni)=stil(2)*ctil(ni-1)+ctil(2)*stil(ni-1); %RAE  % recursion in 3-18b  stil(1) = s_0
        if(n<max)
            lim=n;
        else
            lim=max;
        end
        for m=1:lim
            mi=m+1; %RAE
            mm1=mi-1; %RAE
            %%% additions by MM
            mm2 = mi - 2;     % MM addition... starts from 0
            mp1=mi+1; %RAE
            mxpnm=m*p(ni,mi); %RAE
            bnmtil=c(ni,mi)*ctil(mi)+s(ni,mi)*stil(mi); %RAE
            sumhn=sumhn+normn1(n,m)*p(ni,mp1)*bnmtil; %RAE %norm - since it isnt done prior
            sumgmn=sumgmn+(n+m+1)*p(ni,mi)*bnmtil; %RAE
            bnmtm1=c(ni,mi)*ctil(mm1)+s(ni,mi)*stil(mm1); %RAE
            anmtm1=c(ni,mi)*stil(mm1)-s(ni,mi)*ctil(mm1); %RAE
            sumjn=sumjn+mxpnm*bnmtm1;
            sumkn=sumkn-mxpnm*anmtm1;

            %%% Additions by MM 04/04/2020 - second order recursions
            if (mm2 == 0)
                bnmtm2 = 0; anmtm2 =0;  % MM
            else
                bnmtm2=c(ni,mi)*ctil(mm2)+s(ni,mi)*stil(mm2); %MM
                anmtm2=c(ni,mi)*stil(mm2)-s(ni,mi)*ctil(mm2); %MM
            end

            %%% additions by MM
            mp2 = mi+2; mxmm1xpnm = m*(m-1)*p(ni,mi);
            mxpnmp1 = m*normn1(n,m)*p(ni,mp1);
            mxmpnp1xpnm = m*(n+m+1)*p(ni,mi);
            sumLn = sumLn + (n+m+1)*(n+m+2)*p(ni,mi)*bnmtil;
            sumMn = sumMn + normn2(n,m)*p(ni,mp2)*bnmtil;   %MM - following Gott 7-22
            sumNn = sumNn + mxmm1xpnm*bnmtm2;
            sumOn = sumOn + mxmm1xpnm*anmtm2;
            sumPn = sumPn + (m+n+1)*bnmtil;
            sumQn = sumQn + mxpnmp1*bnmtm1;
            sumRn = sumRn - mxpnmp1*anmtm1;
            sumSn = sumSn + mxmpnp1xpnm*bnmtm1;
            sumTn = sumTn - mxmpnp1xpnm*anmtm1;
        end  % loop for m - inner summation
        sumj=sumj+reorn*sumjn;
        sumk=sumk+reorn*sumkn;
        sumN = sumN + reorn*sumNn; sumO = sumO + reorn*sumOn;
        sumR = sumR + reorn*sumRn; sumQ = sumQ + reorn*sumQn;
        sumS = sumS + reorn*sumSn; sumT = sumT + reorn*sumTn;
    end

    sumh = sumh+reorn*sumhn;
    sumgm = sumgm+reorn*sumgmn;

    sumL = sumL + reorn*sumLn;
    sumM = sumM + reorn*sumMn;
    sumP = sumP + reorn*sumPn;
end    % loop for n - outer summation
lambda=sumgm+ep*sumh;
g(1,1)=-muor2*(lambda*xor-sumj);
g(2,1)=-muor2*(lambda*yor-sumk);
g(3,1)=-muor2*(lambda*zor-sumh);
g=g; %RAE

%%% Jacobian matrix
alp = [0;0;1]; rhat = [xor; yor; zor];
Y = [sumS; sumT; 0]; Z = [sumQ; sumR; 0];
Mat = zeros(3);
Mat(1,1) = sumN; Mat(1,2) = -sumO; Mat(2,1) = -sumO; Mat(2,2) = -sumN;
Uxx1 = (sumL + sumM*ep^2+ 2*sumP*ep + sumgm + 3*sumh*ep)*(rhat*rhat')  ...
    -(sumM*ep + sumP + sumh)*(rhat*alp' + alp*rhat') + sumM*(alp*alp') ...
    -lambda*eye(3) - (Y*rhat' + rhat*Y') + Z*alp' + alp*Z' ...
    - ep*(Z*rhat' + rhat*Z') + Mat;
dgdr = (muor2*ri)*Uxx1;

end