function a=constantQ_GSLS(Q,f0,N)
freq=f0;
f_ratio=100;
QPval=Q;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%单独求取非弹性系数%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
wmean=2*pi*freq;
wmin=wmean/sqrt(f_ratio);
wmax=wmean*sqrt(f_ratio); 
xfrq=[wmin:(wmax-wmin)/(10000-1):wmax];%横坐标
all_relax_num=[N];

for iteration_num=1:length(all_relax_num)
    Relax_num=all_relax_num(iteration_num);
    kmax=2*Relax_num-1;

    for j=1:kmax
        w_disc(j)=exp(log(wmin)+(j-1)/(kmax-1)*log(f_ratio));
        s_disc(j)=1/w_disc(j);
    end
    AP=zeros(kmax,Relax_num);
    QP=ones(kmax,1)/QPval;
    tauP=zeros(Relax_num,1);
    for m=1:kmax
        for j=1:Relax_num
            AP(m,j)=(s_disc(2*j-1).*w_disc(m)-(s_disc(2*j-1).^2.*w_disc(m).^2)/QPval)./(1+s_disc(2*j-1).^2.*w_disc(m).^2);
        end
    end
    tauP=AP\(QP*Relax_num);
    for j=1:Relax_num
        tausigma_j(j)=s_disc(2*j-1);
    end
    
    %% P-wave Q continuous values
    numP=0;
    denP=Relax_num;
    for j=1:Relax_num
        numP=numP+(tauP(j)*s_disc(2*j-1)*xfrq(:))./(1+s_disc(2*j-1)^2*xfrq(:).^2);
        %***正下方-号改为+号，w_disc(2*j-1)改为xfrq(:).^2)
        denP=denP+(tauP(j)*s_disc(2*j-1).^2*xfrq(:).^2)./(1+s_disc(2*j-1)^2*xfrq(:).^2);
    end
    Q_contP=denP./numP;
    Q_all(iteration_num,:)=Q_contP;
%     a=strcat('The mechanisim Number = ',num2str(j))
    taoe=tausigma_j'.*(tauP+1);
   a=[taoe'; tausigma_j];

% tausigma_j
%     %% Computing fitting quality (RMS and maximum difference)
%     maxPdif=0;
%     for j=1:length(Q_contP)
%         tempP=abs(Q_contP(j)-QPval);
%         if tempP >= maxPdif
%             maxPdif=tempP;
%         end
%     end
%     disp(strcat('Maximum QP fitting error=  ',num2str(maxPdif/QPval*100),' %'));
%     

    
%      hold on; semilogx(xfrq/(2*pi),Q_contP);
%      legend(Relax_num); %legend(a) its run not right

% dt=1e-4;t=dt:dt:1;
% for jin=1:Relax_num
% yi(jin,:)=1/tausigma_j(jin)*(1-taoe(jin)/tausigma_j(jin))*exp(-t/tausigma_j(jin));
% end
% y=-(sum(yi));
% y=log(y);
end

% plot(t,y)
% fun1=@(beta1,x1_1)(beta1(1).*x1_1+beta1(2));
% beta1=lsqcurvefit(fun1,[0 0],x1_1,y1_1);
% curve1_1=beta1(1).*x1_1+beta1(2);
% plot(x1_1,curve1_1,'r-')
