function []=scatter_plot_sta(x,y)
X =x;
Y =y;
numbins = 100;
    [values, centers] = hist3([X Y], [numbins numbins]);
    centers_X = centers{1,1};
    centers_Y = centers{1,2};
    binsize_X = abs(centers_X(2) - centers_X(1)) / 2;
    binsize_Y = abs(centers_Y(2) - centers_Y(1)) / 2;
    bins_X = zeros(numbins, 2);
    bins_Y = zeros(numbins, 2);
    for i = 1:numbins
        bins_X(i, 1) = centers_X(i) - binsize_X;
        bins_X(i, 2) = centers_X(i) + binsize_X;
        bins_Y(i, 1) = centers_Y(i) - binsize_Y;
        bins_Y(i, 2) = centers_Y(i) + binsize_Y;
    end
    scatter_COL = zeros(length(X), 1);
    onepercent = round(length(X) / 100);
    for i = 1:length(X)
        if (mod(i,onepercent) == 0)
            fprintf('.');
        end            
        last_lower_X = NaN;
        last_higher_X = NaN;
        id_X = NaN;
        c_X = X(i);
        last_lower_X = find(c_X >= bins_X(:,1));
        if (~isempty(last_lower_X))
            last_lower_X = last_lower_X(end);
        else
            last_higher_X = find(c_X <= bins_X(:,2));
            if (~isempty(last_higher_X))
                last_higher_X = last_higher_X(1);
            end
        end
        if (~isnan(last_lower_X))
            id_X = last_lower_X;
        else
            if (~isnan(last_higher_X))
                id_X = last_higher_X;
            end
        end
        last_lower_Y = NaN;
        last_higher_Y = NaN;
        id_Y = NaN;
        c_Y = Y(i);
        last_lower_Y = find(c_Y >= bins_Y(:,1));
        if (~isempty(last_lower_Y))
            last_lower_Y = last_lower_Y(end);
        else
            last_higher_Y = find(c_Y <= bins_Y(:,2));
            if (~isempty(last_higher_Y))
                last_higher_Y = last_higher_Y(1);
            end
        end
        if (~isnan(last_lower_Y))
            id_Y = last_lower_Y;
        else
            if (~isnan(last_higher_Y))
                id_Y = last_higher_Y;
            end
        end
        scatter_COL(i) = values(id_X, id_Y);
    end  
scatter(x, y, 50, scatter_COL, '.' );
colormap('jet');
h = colorbar;
caxis([0 100000]);
% plot(x,y,'dk','MarkerSize',5,'MarkerFaceColor','k');
%  title('Product Comparison');
xlabel('Sequence identity','FontSize',12,'FontWeight','normal','Color','k');
ylabel('TM-score','FontSize',12,'FontWeight','normal','Color','k');
 xlim([-1 4]);ylim([-1 4])
hold on
N=length(x);
a=polyfit(x,y,1);
C=corrcoef(x,y);
R=C(1,2)^2;
XX = x;
YY = y;
  nb_obs = length(XX);
    obs = XX;
    theo = YY;    
    sum_obs = sum(obs);  %%%%%%%%XX
  sum_theo = sum(theo);  %%%%%%%%%YY
  sum_sq_obs = sum(obs.^2);
  sum_sq_theo = sum(theo.^2);
  buf = theo - obs;  
  sum_diff = sum(buf);
  buf = buf.^2;  
  sum_sq_diff = sum(buf);
  buf = theo .* obs;
  cov = sum(buf);
  rse = sqrt(sum_sq_diff/(nb_obs - 2));
    bias = sum_diff./nb_obs;
  avg_obs = sum(obs)./nb_obs;
  avg_theo = sum(theo)./nb_obs;
  cov = (cov/nb_obs) - (avg_obs .* avg_theo);
  stdv_obs = sqrt((sum_sq_obs - (sum_obs.^2./nb_obs))./nb_obs);  
  stdv_theo = sqrt((sum_sq_theo - (sum_theo.^2./nb_obs))./nb_obs);
  slope_t1 = cov./stdv_obs.^2;
    intercept_t1 = avg_theo - (slope_t1.*avg_obs);
    rsq = (cov./(stdv_obs .* stdv_theo)).^2;   
    MRE = 100*((sum_obs - sum_theo)./sum_theo);   %% (Mod-Ins)/In    
    RE = MRE./nb_obs;     
    shi= abs(MRE)./nb_obs;    
    % Wenzhao edition on the biases calculation
A = XX./YY;
A2 = A(A~=0 & isfinite(A));
% Other way to exclude NAN and Inf: 
% B = A( ~any( isnan( A ) | isinf( A ), 2 ),: )
% ψ_i
ratio= A2 -1;
N=length(ratio);
bias = 100*sum(ratio)./N;
bias_abs = 100*sum(abs(ratio))./N;
% MaxSP=max(x);MaxV=max(y);Maxi=1.1*max(MaxSP,MaxV);
Maxi= 2000;
ax=linspace(0,Maxi,2000);
%plot(ax,10*ax,'k');hold on
axis([0 1 0 1])
set(gca,'FontSize',12)
%set(gca,'XTick',(0:10:150))
%title('BSk')
axis square
text(1*Maxi/500,9.0*Maxi/10,['N = ',num2str(N)], 'FontSize',12)
text(1*Maxi/500,8.25*Maxi/10,['R^2 = ',num2str(round(1000*R)/1000)], 'FontSize',12)
text(1*Maxi/500,7.50*Maxi/10,['RMSE = ',num2str(round(1000*rse)/1000)], 'FontSize',12)
text(1*Maxi/500,6.75*Maxi/10,['\Psi= ',num2str(round(1000*bias)/1000),' %'], 'FontSize',12)
text(1*Maxi/500,6.0*Maxi/10,['|\Psi|= ',num2str(round(1000*bias_abs)/1000),' %'], 'FontSize',12)
%ay=a(1)*ax+a(2);%a=round(100*a)/100;disp(a)
%plot(ax,ay,'r')
%RegrssStr=['Y = ',num2str(round(a(1),2)),'*X +',num2str(round(a(2),2))];
%%% xlim([0 0.08]);ylim([0 0.08])
%legend({'Data points','Y = 10*X',RegrssStr},'Location','northeast','FontSize',12);
box on
set(gca,'LineWidth',1.2)
set(gcf,'color','white')
