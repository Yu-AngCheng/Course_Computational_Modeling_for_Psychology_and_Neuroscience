clear
tic
Value=zeros(6);
rewardtable =...
    [-1, -1, -1, -1, 0, -1;
    -1, -1, -1, 0, -1, 100;
    -1, -1, -1, 0, -1, -1;
    -1, 0, 0, -1, 0, -1;
    0, -1, -1, 0, -1, 100;
    -1, 0, -1, -1, 0, 100];
gamma=0.1;
alpha=0.8;
eps=0.1;
iter_num=100000;
all=1:6;
for i=1:iter_num % ѵ��iter_num��
    state=randi(6);% ÿ���������һ��λ��
    while state~=6 % ���û���յ㣬�ߵ���һ��λ�ã��������λ�õ�value
        % eps greedy
        toss=binornd(1,eps);
        if toss==1 % exploration
            temp=all(rewardtable(state,:)>=0);
            action=temp(randperm(length(temp),1));
        else
            tempmax=max(Value(state,:));
            tempindex=all((Value(state,:)==tempmax)&(rewardtable(state,:)>=0));
            action=tempindex(randperm(length(tempindex),1));
        end
        immediate_reward=rewardtable(state,action);% ��ʱreward
        next_state=action;% ��һ��λ��
        % eps greedy
        toss=binornd(1,eps);
        if toss==1 % exploration
            temp=all(rewardtable(next_state,:)>=0);
            next_action=temp(randperm(length(temp),1));
        else
            tempmax=max(Value(next_state,:));
            tempindex=all((Value(next_state,:)==tempmax)&(rewardtable(next_state,:)>=0));
            next_action=tempindex(randperm(length(tempindex),1));
        end
        prediction_error=immediate_reward+gamma*Value(next_state,next_action)-Value(state,action);
        Value(state,action)=Value(state,action)+alpha*prediction_error;
        state = next_state;
        action = next_action;
    end
end
for j=1:10
    state=randi(3);
    count=0;% ���˶��ٲ�
    path(j).statechain=[];% ��¼��ô�ߵ�
    path(j).statechain=[path(j).statechain state];
    while state~=6
        if count>5 % 5��û�߳�������ʧ��
            path(j).statechain=[];
            break
        end
        tempmax=max(Value(state,:));
        tempindex=all(Value(state,:)==tempmax);
        next_state=tempindex(randperm(length(tempindex),1));
        path(j).statechain=[path(j).statechain next_state];
        state=next_state;
        count=count+1;
        path(j).count=count;
    end
    
end
toc