list=textread('pair.lst','%s');
%num=0;
%sum_num=0;

for i=1:1764
    [tmscore,identity]=textread(strcat('dom20pairs\',char(list(i)),'.dom20pair'),'%f %f');
    x=find(identity <= 0.25);
    y=find(tmscore <= 0.8 );
    j=union(x,y);
    identity(j)=[];
    tmscore(j)=[];
    tmscore_com=cat(1,tmscore_com,tmscore);
    identity_com=cat(1,identity_com,identity);
    %file=textread(strcat('C:\Users\dell\Desktop\pairs\d1fjra_.pair'));
%     tmscore=file(:,1)
%     identity=file(:,2);
%     num=num+nnz(tmscore>2 & identity<0.25);
%     sum_num=sum_num+nnz(tmscore);
%     Locate1=find(tmscore==1);
%     Locate2=find(identity==1);
%     if Locate1 == Locate2
%        tmscore(Locate1)=[];
%        identity(Locate2)=[];
%     end
%     figure(i)
    plot(identity,tmscore,'.')
%     set(gca,'XDir','reverse')
%     set(gca,'yDir','reverse')
    xlabel('Sequence identity')
    ylabel('TM-score')
    hold on;
    %break
%     x=0.25;
%     y=0:0.001:10;
%     plot(x,y,'-')
%     hold on;
%     y=2;
%     x=0:0.001:1;
%     plot(x,y,'-')
end

cftool(identity_com,tmscore_com)
%plot(identity_com,tmscore_com,'.')
%xlabel('Sequence identity')
%ylabel('TM-score')

% ratio=num/sum_num
% sum_num
%ratio=num/39938724
%ratio=num/39942252
%ratio=num/67970448
% x=0.25;
% y=0:0.001:1;
% plot(x,y,'-')
% hold on;
% y=1;
% x=0:0.001:1;
% plot(x,y,'-')