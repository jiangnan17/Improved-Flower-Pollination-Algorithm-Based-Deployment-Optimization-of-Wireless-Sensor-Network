%%��Բ
function draw_match(first_init_wolf,best_indivi,sensor_s)
    global N;
    [match,~] = get_match_value(first_init_wolf,best_indivi,sensor_s);
    %�����ݼ��س���
    x_pos_init = first_init_wolf(1,:,1);  %���ڽ���������
    y_pos_init = first_init_wolf(2,:,1);
    x_pos_best = best_indivi(1,:);%���ڽ���������
    y_pos_best = best_indivi(2,:);
    
    for k=1:N
        plot(x_pos_init(match(k,1)),y_pos_init(match(k,1)),'.','MarkerSize',10);%����ʼ���Ľ��Բ��
        plot(x_pos_best(match(k,2)),y_pos_best(match(k,2)),'.','MarkerSize',10)%���Ż���Ľ��Բ��
        
        line([x_pos_init(match(k,1)),x_pos_best(match(k,2))],[y_pos_init(match(k,1)),y_pos_best(match(k,2))]...
               ,'LineWidth',2);
        %�ұ�
        if rem(k,4)==1
            %��ʼ���Ľ��л�λ��
            text(x_pos_init(match(k,1))+0.5,y_pos_init(match(k,1)),'N','FontSize',15);%���б��
            text(x_pos_init(match(k,1))+2,y_pos_init(match(k,1))-0.5,['',num2str(match(k,1)),'i'],'FontSize',10);%���б��
            %�Ż����Ļ�λ��
            text(x_pos_best(match(k,2))+0.5,y_pos_best(match(k,2)),'N','FontSize',15);%���б��
            text(x_pos_best(match(k,2))+2,y_pos_best(match(k,2))-0.5,['',num2str(match(k,2)),'e'],'FontSize',10);%���б��
            
            
        elseif rem(k,4)==2%�ұ�
            text(x_pos_init(match(k,1))-2.5,y_pos_init(match(k,1)),'N','FontSize',15);%���б��
            text(x_pos_init(match(k,1))-1,y_pos_init(match(k,1))-0.5,['',num2str(match(k,1)),'i'],'FontSize',10);%���б��
            
            
            text(x_pos_best(match(k,2))-2.5,y_pos_best(match(k,2)),'N','FontSize',15);%���б��
            text(x_pos_best(match(k,2))-1,y_pos_best(match(k,2))-0.5,['',num2str(match(k,2)),'e'],'FontSize',10);%���б��
        elseif rem(k,4)==3
            %�ϱ�
            text(x_pos_init(match(k,1))- 0.5,y_pos_init(match(k,1))+ 1.5,'N','FontSize',15);%���б��
            text(x_pos_init(match(k,1))+1,y_pos_init(match(k,1))+0.5,['',num2str(match(k,1)),'i'],'FontSize',10);%���б��
        
        
            text(x_pos_best(match(k,2))- 0.5,y_pos_best(match(k,2))+ 1.5,'N','FontSize',15);%���б��
            text(x_pos_best(match(k,2))+1,y_pos_best(match(k,2))+0.5,['',num2str(match(k,2)),'e'],'FontSize',10);%���б��
        else%�±�
            text(x_pos_init(match(k,1))-0.5,y_pos_init(match(k,1))-1.5,'N','FontSize',15);%���б��
            text(x_pos_init(match(k,1))+1,y_pos_init(match(k,1))-2,['',num2str(match(k,1)),'i'],'FontSize',10);%���б��
        
            
            text(x_pos_best(match(k,2))-0.5,y_pos_best(match(k,2))-1.5,'N','FontSize',15);%���б��
            text(x_pos_best(match(k,2))+1,y_pos_best(match(k,2))-2,['',num2str(match(k,2)),'e'],'FontSize',10);%���б��
        
        end
        
        hold on;%�������  ��Ȼ��һ��Բ��
        axis([0,50,0,50]);%����������С�����ֵ  
        set(gca,'xtick',(0:2:50));%����x���경��Ϊ1
        set(gca,'ytick',(0:2:50));%����y���경��Ϊ1   ���и����� ����������̫������
        axis square;%ʹ���ݱ���Ϊ1  ����Բ����Բ  ��Ȼ����Բ
        set(gca,'yminorgrid','on');%��С����
        set(gca,'xminorgrid','on');
        grid on;%����
        hold on;%��ͼһֱ��������
    end
    
    % ���������ĵ��ͼ������
    % for i=1:L/1
    %     for j=1:W/1
    %         plot(Grid_cen_x(i),Grid_cen_y(j),'.','MarkerFaceColor','r');
    %         hold on;%��ͼһֱ��������
    %     end
    % end

        %������
        x1 = [25,35,25,15,25];
        y1 = [35,25,15,25,35];
        line(x1,y1,'LineWidth',1,'LineStyle','-','color','k');
        fill(x1,y1,[0.5 0.5 0.5]);
        hold on;

        %���ĸ�������
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