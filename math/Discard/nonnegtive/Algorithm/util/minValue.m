% ���������������������������������
% x ������  scale ������������ number ������ ||b-Ax|| ����
function [minNum,v]=minValue(x,scale,number,A,b,f,c)
[m,n]=size(x);
fx=f(A,b,x);
v=zeros(m,number);
minNum=0;
index=0;
while 1
    %�����������
		% ���������
    xi=x+(2*rand(m,n)-ones(m,n))*scale;
		% ������������
    g=c(xi);
    %����Լ������
		% ��
    if g>0
				% ������		
        fxi=f(A,b,xi);
        %�ȸ�����Сֵ��Ҫ�a�
				% ������������� ������
        if fxi<fx
            minNum=minNum+1;
            v(:,minNum)=xi;
        end
        %ֻҪ����Լ����ֵ��������Чֵ
        index=index+1;
        %ʵ��ֵ�ȸ�����Ҫ�� ��ֹͣ�ж�
        if index>number
            break;
        end
    end
end; 
