%%��Բ
function draw_circle(x_pos,y_pos,sersors_r)
    global N;
    global L;
    global W;
    global r;
    global r1;%��İ뾶
    global Grid_cen_x;
    global Grid_cen_y;
    global ger;
    angle=0:pi/100:2*pi;%�Ƕ�
    for k=1:N
        %�����²�ͬ�뾶��ɫ��ͬ
        if sersors_r(1,k) == 7
            plot(x_pos(1,k),y_pos(1,k),'.','MarkerSize',15,'color','r');%��Բ��  %�Ȼ�Բ�� �ٻ�Բ  ��Ȼ����һ��  ��֪զ���
            hold on;
            plot(sersors_r(1,k)*cos(angle)+x_pos(1,k),sersors_r(1,k)*sin(angle)+y_pos(1,k),'color','r');%��Բ
            hold on;
        elseif sersors_r(1,k) == 6
            plot(x_pos(1,k),y_pos(1,k),'.','MarkerSize',15,'color','g');%��Բ��  %�Ȼ�Բ�� �ٻ�Բ  ��Ȼ����һ��  ��֪զ���
            plot(sersors_r(1,k)*cos(angle)+x_pos(1,k),sersors_r(1,k)*sin(angle)+y_pos(1,k),'color','g');%��Բ
            hold on;
        else
            plot(x_pos(1,k),y_pos(1,k),'.','MarkerSize',15,'color','b');%��Բ��  %�Ȼ�Բ�� �ٻ�Բ  ��Ȼ����һ��  ��֪զ���
            plot(sersors_r(1,k)*cos(angle)+x_pos(1,k),sersors_r(1,k)*sin(angle)+y_pos(1,k),'color','b');%��Բ
            hold on;
        end
        
        if rem(k,4)==1
            text(x_pos(1,k)+0.5,y_pos(1,k),'N','FontSize',15);%���б��
            text(x_pos(1,k)+2,y_pos(1,k)-0.5,num2str(k),'FontSize',10);%���б��
        elseif rem(k,4)==2%�ұ�
            text(x_pos(1,k)-2.5,y_pos(1,k),'N','FontSize',15);%���б��
            text(x_pos(1,k)-1,y_pos(1,k)-0.5,num2str(k),'FontSize',10);%���б��
        elseif rem(k,4)==3
            %�ϱ�
            text(x_pos(1,k)- 0.5,y_pos(1,k)+ 1.5,'N','FontSize',15);%���б��
            text(x_pos(1,k)+1,y_pos(1,k)+0.5,num2str(k),'FontSize',10);%���б��
        else%�±�
            text(x_pos(1,k)-0.5,y_pos(1,k)-1.5,'N','FontSize',15);%���б��
            text(x_pos(1,k)+1,y_pos(1,k)-2,num2str(k),'FontSize',10);%���б��
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
    
%     ���������ĵ��ͼ������
%     for i=1:L/1
%         for j=1:W/1
%             plot(Grid_cen_x(i),Grid_cen_y(j),'.','MarkerFaceColor','r');
%             hold on;%��ͼһֱ��������
%         end
%     end
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