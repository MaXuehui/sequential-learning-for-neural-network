clear;
clc;

x = -4:0.01:4;

space_num = 100;  %申请空间
% u = -5:1:5;
u = 0*ones(space_num,1);


theta = 0.1*ones(2*space_num+1,1);
P = 1000*eye(2*space_num+1);
K = ones(2*space_num+1,1);

Q0 = 0.01*eye(2*space_num+1);

index = 0;

delta_k = 0.1;

epsilon_max = 0.5;
epsilon_min = 0.2;
epsilon = epsilon_max;

for i = 1:800

    f(i) = 1.1*(1-x(i)+2*x(i)*x(i))*exp(-1/2*x(i)*x(i));

    if index == 0
        index = index+1;
        u(index) = x(i);
        delta_2(index) = 2;
    else
        for j = 1:index
            nearst(j) = abs(x(i)-u(j));
        end
        nearst_pot = min(nearst);
    
        if nearst_pot>epsilon
            if abs(f_hat(i)-f(i))>0.05
                index = index+1;
                u(index) = x(i);
                delta_2(index) = abs(x(i)-nearst_pot);
                theta(index+1) = f_hat(i)-f(i);
                
                if epsilon>epsilon_min
                    epsilon = 0.9*epsilon;
                else
                    epsilon = epsilon_min;
                end

            end
        end
    end

    for k = 1:index
        phi(k,1) = exp(-(norm(x(i)-u(k)))^2/delta_2(k));
    end


    K(1:index+1) = inv(1+[1; phi]'*P(1:index+1,1:index+1)*[1; phi])*P(1:index+1,1:index+1)*[1; phi];
    theta(1:index+1) = theta(1:index+1) + K(1:index+1)*(f(i)-(theta(1:index+1))'*[1; phi]);
    P(1:index+1,1:index+1) = (eye(index+1)-K(1:index+1)*[1; phi]')*P(1:index+1,1:index+1)+Q0(1:index+1,1:index+1);

    f_hat(i) =   [1; phi]' * theta(1:index+1);
    for k = 1:index
        phi(k,1) = exp(-(norm(x(i+1)-u(k)))^2/delta_2(k));
    end
    f_hat(i+1) =   [1; phi]' * theta(1:index+1);


end

% for i = 1:800
% 
%     for k = 1:index
%         phi(k,1) = exp(-(norm(x(i)-u(k)))^2/2);
%     end  
% 
%     f_hat(i) = [1; phi(1:index)]' * theta(1:index+1);
% end

plot(x(1:800),f,'o');hold on;
plot(x(1:801),f_hat,'.')








