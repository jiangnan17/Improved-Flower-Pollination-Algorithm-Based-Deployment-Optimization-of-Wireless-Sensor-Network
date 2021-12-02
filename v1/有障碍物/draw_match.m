%%画圆
function draw_match(first_init_wolf,best_indivi,sensor_s)
    global N;
    [match,~] = get_match_value(first_init_wolf,best_indivi,sensor_s);
    %把数据加载出来
    x_pos_init = first_init_wolf(1,:,1);  %后期解决这个警告
    y_pos_init = first_init_wolf(2,:,1);
    x_pos_best = best_indivi(1,:);%后期解决这个警告
    y_pos_best = best_indivi(2,:);
    
    for k=1:N
        plot(x_pos_init(match(k,1)),y_pos_init(match(k,1)),'.','MarkerSize',10);%画初始化的结点圆心
        plot(x_pos_best(match(k,2)),y_pos_best(match(k,2)),'.','MarkerSize',10)%画优化后的结点圆心
        
        line([x_pos_init(match(k,1)),x_pos_best(match(k,2))],[y_pos_init(match(k,1)),y_pos_best(match(k,2))]...
               ,'LineWidth',2);
        %右边
        if rem(k,4)==1
            %初始化的进行画位置
            text(x_pos_init(match(k,1))+0.5,y_pos_init(match(k,1)),'N','FontSize',15);%进行标记
            text(x_pos_init(match(k,1))+2,y_pos_init(match(k,1))-0.5,['',num2str(match(k,1)),'i'],'FontSize',10);%进行标记
            %优化完后的画位置
            text(x_pos_best(match(k,2))+0.5,y_pos_best(match(k,2)),'N','FontSize',15);%进行标记
            text(x_pos_best(match(k,2))+2,y_pos_best(match(k,2))-0.5,['',num2str(match(k,2)),'e'],'FontSize',10);%进行标记
            
            
        elseif rem(k,4)==2%右边
            text(x_pos_init(match(k,1))-2.5,y_pos_init(match(k,1)),'N','FontSize',15);%进行标记
            text(x_pos_init(match(k,1))-1,y_pos_init(match(k,1))-0.5,['',num2str(match(k,1)),'i'],'FontSize',10);%进行标记
            
            
            text(x_pos_best(match(k,2))-2.5,y_pos_best(match(k,2)),'N','FontSize',15);%进行标记
            text(x_pos_best(match(k,2))-1,y_pos_best(match(k,2))-0.5,['',num2str(match(k,2)),'e'],'FontSize',10);%进行标记
        elseif rem(k,4)==3
            %上边
            text(x_pos_init(match(k,1))- 0.5,y_pos_init(match(k,1))+ 1.5,'N','FontSize',15);%进行标记
            text(x_pos_init(match(k,1))+1,y_pos_init(match(k,1))+0.5,['',num2str(match(k,1)),'i'],'FontSize',10);%进行标记
        
        
            text(x_pos_best(match(k,2))- 0.5,y_pos_best(match(k,2))+ 1.5,'N','FontSize',15);%进行标记
            text(x_pos_best(match(k,2))+1,y_pos_best(match(k,2))+0.5,['',num2str(match(k,2)),'e'],'FontSize',10);%进行标记
        else%下边
            text(x_pos_init(match(k,1))-0.5,y_pos_init(match(k,1))-1.5,'N','FontSize',15);%进行标记
            text(x_pos_init(match(k,1))+1,y_pos_init(match(k,1))-2,['',num2str(match(k,1)),'i'],'FontSize',10);%进行标记
        
            
            text(x_pos_best(match(k,2))-0.5,y_pos_best(match(k,2))-1.5,'N','FontSize',15);%进行标记
            text(x_pos_best(match(k,2))+1,y_pos_best(match(k,2))-2,['',num2str(match(k,2)),'e'],'FontSize',10);%进行标记
        
        end
        
        hold on;%加上这句  不然少一个圆心
        axis([0,50,0,50]);%横纵坐标最小和最大值  
        set(gca,'xtick',(0:2:50));%设置x坐标步长为1
        set(gca,'ytick',(0:2:50));%设置y坐标步长为1   还有个问题 横坐标数字太紧凑了
        axis square;%使横纵比例为1  这样圆才像圆  不然像椭圆
        set(gca,'yminorgrid','on');%最小网格
        set(gca,'xminorgrid','on');
        grid on;%网格化
        hold on;%把图一直保存起来
    end
    
    % 把网格中心点的图画出来
    % for i=1:L/1
    %     for j=1:W/1
    %         plot(Grid_cen_x(i),Grid_cen_y(j),'.','MarkerFaceColor','r');
    %         hold on;%把图一直保存起来
    %     end
    % end

        %画菱形
        x1 = [25,35,25,15,25];
        y1 = [35,25,15,25,35];
        line(x1,y1,'LineWidth',1,'LineStyle','-','color','k');
        fill(x1,y1,[0.5 0.5 0.5]);
        hold on;

        %画四个三角形
        x1 = [0,0,15,0];
        y1 = [50,35,50,50];
        x2 = [0,0,15,0];
        y2 = [15,0,0,15];
        x3 = [35,50,50,35];
        y3 = [50,50,35,50];
        x4 = [35,50,50,35];
        y4 = [0,0,15,0];
        line(x1,y1,'LineWidth',1,'LineStyle','-','color','k');
        line(x2,y2,'LineWidth',1,'LineStyle','-','color','k');
        line(x3,y3,'LineWidth',1,'LineStyle','-','color','k');
        line(x4,y4,'LineWidth',1,'LineStyle','-','color','k');
        fill(x1,y1,[0.5 0.5 0.5]);
        fill(x2,y2,[0.5 0.5 0.5]);
        fill(x3,y3,[0.5 0.5 0.5]);
        fill(x4,y4,[0.5 0.5 0.5]);
        hold on;

    

    
end