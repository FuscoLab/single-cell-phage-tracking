%% analyzing gene expression data


%%read all gene expression traces
path_directory='/Users/diana1/Desktop/Projects/PhageMM/gene_expression';
stats_file=[path_directory '/analysis_03_correlations_all.csv'];
stats_table=readtable(stats_file);
deltat=0.5; %%(time interval in mins)

%%add columns to stat table with regards to fitting parameters
newTabCol = zeros(height(stats_table), 1);

stats_table.('lambda1')=newTabCol;
stats_table.('lambda2')=newTabCol;
stats_table.('t0')=newTabCol
stats_table.('t0fit')=newTabCol;
stats_table.('K0')=newTabCol;
stats_table.('prod_time')=newTabCol;


max_length=max(stats_table.cell_length_at_plateau);
min_length=min(stats_table.cell_length_at_plateau);
color=colormap(redblue);
original_files=dir([path_directory '/expression_raw_data/*.csv']);
for k=1:size(original_files)
 %for k=1:1
     filename=[path_directory '/expression_raw_data/' original_files(k).name];
    a=split(original_files(k).name,'_');
    b=split(a{end},'.');
    cell_id=str2num(b{1});
    for j=1:size(stats_table,1)
        if(cell_id==stats_table.idx(j))
            cell_length=stats_table.cell_length_at_plateau(j);
            break;
        end
    end
    n=floor((log(cell_length)-log(min_length))*255/(log(max_length)-log(min_length)))+1;

    M=readtable(filename);
    
    %%clean last point(s) if it drops 
    while(M.total_YFP_bg_corrected(end)<M.total_YFP_bg_corrected(end-1))
        M(end,:)=[];
    end


    figure(3)
    subplot(1,2,1)
    plot(M.time_min-M.production_start_time,M.total_YFP_bg_corrected,'LineWidth',2,'Color',color(n,:));
    %set(gca,'FontSize',20);
    xlabel('Time')
    ylabel('YFP Fluorescence Signal')

    hold on
    subplot(1,2,2)
    plot(M.total_YFP_bg_corrected(2:end),diff(M.total_YFP_bg_corrected)/deltat,'LineWidth',2,'Color',color(n,:));
    %set(gca,'FontSize',20);
    xlabel('YFP')
    ylabel('dYFP/dt')

    hold on
    
   

    %%parameter fitting for each gene-expression time-series (take only
    %%10 time points before t00)
    gene_expression{cell_id}=M;
    for i=1:size(M,1)
        if(M.time_min(i)<M.production_start_time(i)-5)
            gene_expression{cell_id}(1,:)=[];
        end
    end

    xdata=gene_expression{cell_id}.total_YFP_bg_corrected(2:end);
    ydata=diff(gene_expression{cell_id}.total_YFP_bg_corrected)/deltat;

    fun=@(lambda,xdata)lambda(1)*(1-exp(-lambda(2)*xdata));
    lambda0=[max(ydata),0.00001];
    lambda{cell_id}=lsqcurvefit(fun,lambda0,xdata,ydata);
    lambda{cell_id}(2)
    figure(4)
    times = linspace(xdata(1),xdata(end));
    plot(xdata,ydata,'o','MarkerFaceColor',color(n,:),'MarkerEdgeColor',color(n,:));
    hold on
    plot(times,fun(lambda{cell_id},times),'--','Color',color(n,:),'LineWidth',2);

    for j=1:size(stats_table,1)
        if(cell_id==stats_table.idx(j))
            stats_table.lambda1(j)=lambda{cell_id}(2)*lambda{cell_id}(1);
            stats_table.lambda2(j)=lambda{cell_id}(1);
            break;
        end
    end



end

%%analysis of the parameters found by the fitting procedure

sub_table=stats_table(:,[3,4,11,13,14]);
sub_table.Properties.VariableNames = ["t2t5","t3t5","Size","l1","l2"]

figure(5)
[~,~,h]=corrplot(sub_table);


%%use lambda1 and lambda2 to fit YFP as a function of time and estimate
%%time0 and K_0

for i=1:size(stats_table,1)
%for i=1:3
%for i=1:1
    %%parameter fitting for each gene-expression time-series
    id=stats_table.idx(i);
    lambda1=stats_table.lambda1(i);
    lambda2=stats_table.lambda2(i);

    xdata=gene_expression{id}.time_min;
    ydata=gene_expression{id}.total_YFP_bg_corrected;
    t00=gene_expression{id}.production_start_time(1);
    n=floor((log(stats_table.cell_length_at_plateau(i))-log(min_length))*255/(log(max_length)-log(min_length)))+1;

    fun2=@(pars,xdata)lambda2/lambda1*log(1+(exp(lambda1*(xdata-pars(1)))-1)/pars(2));

    %set initial condition for fit search
    p=polyfit(xdata(end-3:end),ydata(end-3:end),1);
    b=exp((-p(2)/lambda2-t00)*lambda1);
    pars0=[t00,30];
    pars{id}=lsqcurvefit(fun2,pars0,xdata,ydata);
    stats_table.t0(i)=t00;
    stats_table.t0fit(i)=pars{id}(1);
    stats_table.K0(i)=pars{id}(2);
    stats_table.prod_time(i)=stats_table.t3t5(i)+t00-pars{id}(1);

    %%plot fit to check
    figure(6)
    times = linspace(xdata(1),xdata(end));
    plot(xdata,ydata,'o','MarkerFaceColor',color(n,:),'MarkerEdgeColor',color(n,:));
    hold on
    plot(times,fun2(pars{id},times),'--','Color',color(n,:),'LineWidth',2);

     %%rescaled plot at time 0
     figure(7)
     plot(xdata-pars{id}(1),ydata,'.','MarkerFaceColor',color(n,:),'MarkerEdgeColor',color(n,:));
     hold on
     plot(times-pars{id}(1),fun2([0,pars{id}(2)],times-pars{id}(1)),'--','Color',color(n,:),'LineWidth',1);



end




sub_table=stats_table(:,[3,4,7,11,13,14,17,18]);
sub_table.Properties.VariableNames = ["t2t5","t3t5","YFP","Size","l1","l2","K0","time"]

figure(10)
[~,~,h]=corrplot(sub_table);

