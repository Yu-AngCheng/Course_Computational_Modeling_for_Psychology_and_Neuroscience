%% Temporal Difference Learning (Sarsa)
clear
tic
% value_table for 6 states
Value=zeros(6);
% -1 means impermissible transition
rewardtable =...
    [-1, -1, -1, -1, 0, -1;
    -1, -1, -1, 0, -1, 100;
    -1, -1, -1, 0, -1, -1;
    -1, 0, 0, -1, 0, -1;
    0, -1, -1, 0, -1, 100;
    -1, 0, -1, -1, 0, 100];
gamma=0.1;% temporal discounting
alpha=0.8;% learning_rate
eps=0.1;% epsilon
iter_num=100000;
all=1:6;

% train
for i=1:iter_num % train iter_num times
    state=randi(6);% randomly put it at the first iteration
    while state~=6 % loop until destination
        % epsilon-greedy
        toss=binornd(1,eps);
        % exploration
        % randomly choose the action
        if toss==1 
            temp=all(rewardtable(state,:)>=0);
            action=temp(randperm(length(temp),1));
        % exploitation
        % choose the action with the greatest reward
        else 
            tempmax=max(Value(state,:));
            tempindex=all((Value(state,:)==tempmax)&(rewardtable(state,:)>=0));
            action=tempindex(randperm(length(tempindex),1));
        end
        immediate_reward=rewardtable(state,action);
        next_state=action;
        % epsilon-greedy
        toss=binornd(1,eps);
        % exploration
        % randomly choose the next action
        if toss==1
            temp=all(rewardtable(next_state,:)>=0);
            next_action=temp(randperm(length(temp),1));
        else
            tempmax=max(Value(next_state,:));
            tempindex=all((Value(next_state,:)==tempmax)&(rewardtable(next_state,:)>=0));
            next_action=tempindex(randperm(length(tempindex),1));
        end
        % calculate the prediction error
        prediction_error=immediate_reward+gamma*Value(next_state,next_action)-Value(state,action);
        % update the value table
        Value(state,action)=Value(state,action)+alpha*prediction_error;
        state = next_state;
        action = next_action;
    end
end
% test
for j= 1:10
    state=randi(3);% the starting point can only be at state 1-3
    count=0;% the number of steps taken
    path(j).statechain=[];% record the path
    path(j).statechain=[path(j).statechain state];
    while state~=6
        if count > 5 % if more than 5 steps, then fail
            path(j).statechain=[];
            break
        end
        % greedy
        tempmax=max(Value(state,:));
        tempindex=all(Value(state,:)==tempmax);
        next_state=tempindex(randperm(length(tempindex),1));
        path(j).statechain=[path(j).statechain next_state]; % record the path
        state=next_state; % update the state
        count=count+1;
        path(j).count=count;
    end
    
end
toc