function [Y,numIter,tElapsed]=kurnnls(XtX,AtA,AtX,option)
% Kernel UR-NNLS: phi(X)=A_phi Y, s.t. Y>0
% XtX, kernel matrix, k(X,X).
% AtA, kernel matrix, k(A,A).
% AtX, kernel matrix, k(A,X).
% option.iter, max number of interations
% option.kernel, string, kernel function
% option.param, scalar or column vector, the parameter of kernel
% option.dis, boolen scalar, false: (not display information) or true
% (display).
% option.residual, if the ||X-XfitThis||<=option.residual, halt.
% option.tof, if ||XfitPrevious-XfitThis||<=option.tof, halt.
% Y, the coefficient matrix.
% numIter, the number of total iterations.
% tElapsed, the time consumed in second.
% [1]\bibitem{cibcb2012}
%     Y. Li and A. Ngom,
%     ``A New Kernel Non-Negative Matrix Factorization and Its Application in Microarray Data Analysis,''
%     {\it CIBCB},
%     submited, 2012.
% Yifeng Li, Otc 28, 2011
% yifeng.li.cn@gmail.com

tStart=tic;
optionDefault.iter=1000;
optionDefault.kernel='linear';
optionDefault.param=[];
optionDefault.dis=1;
optionDefault.residual=1e-4;
optionDefault.tof=1e-4;
if nargin<4
   option=[]; 
end
option=mergeOption(option,optionDefault);

Y=pinv(AtA)*AtX;

% Y=A\X;
Y(Y<0)=0;
Ap=(abs(AtX)+AtX)./2;
An=(abs(AtX)-AtX)./2;
Bp=(abs(AtA)+AtA)./2;
Bn=(abs(AtA)-AtA)./2;
prevRes=Inf;
for i=1:option.iter
    Y=Y.*sqrt((Ap + Bn*Y)./(An + Bp*Y));
    if mod(i,100)==0 || i==option.iter
        if option.dis
            disp(['Iterating >>>>>> ', num2str(i),'th']);
        end
        YtApn=Y'*AtX;
        curRes=trace(XtX-YtApn-YtApn'+Y'*AtA*Y); % use abs because the trace maybe negative due to numerical reasons
        fitRes=prevRes-curRes;%abs(prevRes-curRes);
        prevRes=curRes;
        if option.tof>=fitRes || option.residual>=curRes || i==option.iter
            disp(['KURNNLS successes!, # of iterations is ',num2str(i),'. The final residual is ',num2str(curRes)]);
            numIter=i;
            finalResidual=curRes;
            break;
        end
    end
end
tElapsed=toc(tStart);
end
