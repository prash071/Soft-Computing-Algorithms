% program for particle swarm optimization
% Stagg and el-abiad system

%initialise population parameters
c1 = 1.8;
c2 = 1.6;

% weight factor

w=1;

% create population 
%No. of particles - m 
m=10

%Number of components 
n=9
% no. of generators
n1=2
% no. of transmission lines
n2=7




%range of forced outage rate 
%(0,0.05) for generators
%(0,0.2) for transmission lines


d1=2;d2=7;
ub1=ones(1,d1)*0.05; %/*upper bounds of the parameters. */
lb1=ones(1,d1)*(0.0001);%/*lower bound of the parameters.*/

ub2=ones(1,d2)*0.2;
lb2=ones(1,d2)*0.0001;
Range1 = repmat((ub1-lb1),[m 1]);
Lower1 = repmat(lb1, [m 1]);
x1 = rand(m,d1) .* Range1 + Lower1;

Range2 = repmat((ub2-lb2),[m 1]);
Lower2 = repmat(lb2, [m 1]);
x2 = rand(m,d2) .* Range2 + Lower2;




x = [x1,x2];
%disp(x);
ide=ones(m,n);

u=ide-x;

frac=u./x;


% velocities of particles
v= -0.01+rand(m,n)*0.02;

% calculate the fitness of each particle

 install_cost=[0.20 0.12 0.10 0.10 0.10 0.10 0.10 0.10 0.10];
 
 component_cost=[0.032 0.032 0.042 0.042 0.042 0.042 0.042 0.042 0.042];
 
 system_cost=zeros(m,1);
 
 sum_icost=sum(install_cost);
 
 for(i=1:m)
 system_cost(i,1) = sum_icost + sum(component_cost.*frac(i,:));
 
 
 system_cost(i,1) = 0.0000265 * edns(x(i,:))*8760 + system_cost(i,1);
 
 end
 
 %for total cost
 
 %min fitness
 
 [value index] = min(system_cost);
 
 % particle with index is the best of this swarm
 gbest_fitness=value;
 gbest=x(index,:);
 
 % personal best of each particle
 pbest_fitness = system_cost;
 pbest = x;
 % iteration to edit particle position and velocity
 % v= v + c1 * r1(pbest - x) + c2 * r2 (gbest - x)
 % x = x + v
 
 
 newv = v;
 newx = x;
 
 storegbestparam=zeros(150,9);
 storegbest=zeros(1,150);
 iter=1;
 for(it=1:150)
 
 for(i=1:m)
   for(j=1:n)
       
     newv(i,j) = w*newv(i,j) + c1 * rand() * (pbest(i,j) - x(i,j)) + c2 * rand() * (gbest(1,j) - x(i,j));
     
        if(newv(i,j) > 0.02)
            newv(i,j)=0.02;
        elseif(newv(i,j) < -0.02)
            newv(i,j) = -0.02;
        end
     newx(i,j)= newx(i,j) + newv(i,j);
     
     if(newx(i,j) > 0.05  && j<3)
         newx(i,j) = 0.05;
     elseif(newx(i,j) <= 0 && j <3 )
         newx(i,j) = 0.0001;
     elseif(newx(i,j) > 0.2 && j>2)
         newx(i,j) = 0.2;
     elseif(newx(i,j) <= 0 && j >2)
         newx(i,j) = 0.0001;
     end
     
   end
 end
 
 newu = ide - newx;
 
 newfr= newu./newx;
 
 for(i=1:m)
 system_cost(i,1) = sum_icost + sum(component_cost.*newfr(i,:));
 system_cost(i,1) = 0.0000265 * edns(newx(i,:))*8760 + system_cost(i,1);
 end
 
 [value index] = min(system_cost);
 
 if(gbest_fitness> value)
    gbest_fitness=value;
    gbest=newx(index,:);
 end  
  for(i=1:m)
    if(pbest_fitness(i,1) > system_cost(i,1))
        pbest_fitness(i,1) = system_cost(i,1);
        pbest(i,:) = newx(i,:);
    end
  end
  
  %gbest
  gbest_fitness
 storegbestparam(iter,:)=gbest;
 storegbest(1,iter)=gbest_fitness;
 iter=iter+1;
 end
 
 