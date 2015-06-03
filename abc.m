%/* ABC algorithm coded using MATLAB language */

%/* Artificial Bee Colony (ABC) is one of the most recently defined algorithms by Dervis Karaboga in 2005, 
%motivated by the intelligent behavior of honey bees. */





clear all
close all
clc

%/* Control Parameters of ABC algorithm*/
NP=20; %/* The number of colony size (employed bees+onlooker bees)*/
FoodNumber=NP/2; %/*The number of food sources equals the half of the colony size*/
limit=90; %/*A food source which could not be improved through "limit" trials is abandoned by its employed bee*/
maxCycle=150; %/*The number of cycles for foraging {a stopping criteria}*/


%/* Problem specific variables*/
%objfun='Sphere'; %cost function to be optimized
D=9; %/*The number of parameters of the problem to be optimized*/

% create population 
%No. of particles - m 
m=10;

%Number of components 
n=9;
% no. of generators
n1=2;
% no. of transmission lines
n2=7;





d1=2;d2=7;
ub1=ones(1,d1)*0.05; %/*upper bounds of the parameters. */
lb1=ones(1,d1)*(0.0001);%/*lower bound of the parameters.*/

ub2=ones(1,d2)*0.2;
lb2=ones(1,d2)*0.0001;
runtime=1;%/*Algorithm can be run many times in order to see its robustness*/





%Foods [FoodNumber][D]; /*Foods is the population of food sources. Each row of Foods matrix is a vector holding D parameters to be optimized. The number of rows of Foods matrix equals to the FoodNumber*/
%ObjVal[FoodNumber];  /*f is a vector holding objective function values associated with food sources */
%Fitness[FoodNumber]; /*fitness is a vector holding fitness (quality) values associated with food sources*/
%trial[FoodNumber]; /*trial is a vector holding trial numbers through which solutions can not be improved*/
%prob[FoodNumber]; /*prob is a vector holding probabilities of food sources (solutions) to be chosen*/
%solution [D]; /*New solution (neighbour) produced by v_{ij}=x_{ij}+\phi_{ij}*(x_{kj}-x_{ij}) j is a randomly chosen parameter and k is a randomly chosen solution different from i*/
%ObjValSol; /*Objective function value of new solution*/
%FitnessSol; /*Fitness value of new solution*/
%neighbour, param2change; /*param2change corresponds to j, neighbour corresponds to k in equation v_{ij}=x_{ij}+\phi_{ij}*(x_{kj}-x_{ij})*/
%GlobalMin; /*Optimum solution obtained by ABC algorithm*/
%GlobalParams[D]; /*Parameters of the optimum solution*/
%GlobalMins[runtime]; /*GlobalMins holds the GlobalMin of each run in multiple runs*/

GlobalMins=zeros(1,runtime);
Storemins=zeros(1,maxCycle);
Storeglobalparams=zeros(maxCycle,D);
for r=1:runtime
  
% /*All food sources are initialized */
%/*Variables are initialized in the range [lb,ub]. If each parameter has different range, use arrays lb[j], ub[j] instead of lb and ub */

Range1 = repmat((ub1-lb1),[FoodNumber 1]);
Lower1 = repmat(lb1, [FoodNumber 1]);
Foods1 = rand(FoodNumber,d1) .* Range1 + Lower1;

Range2 = repmat((ub2-lb2),[FoodNumber 1]);
Lower2 = repmat(lb2, [FoodNumber 1]);
Foods2 = rand(FoodNumber,d2) .* Range2 + Lower2;

Foods = [Foods1, Foods2];

ide=ones(FoodNumber,D);

av=ide-Foods;

frac=av./Foods;


%ObjVal=feval(objfun,Foods);
%Fitness=calculateFitness(ObjVal);

% calculate the fitness of each particle

 install_cost=[0.20 0.12 0.10 0.10 0.10 0.10 0.10 0.10 0.10];
 
 component_cost=[0.032 0.032 0.042 0.042 0.042 0.042 0.042 0.042 0.042];
 
 system_cost=zeros(FoodNumber,1);
 
 sum_icost=sum(install_cost);
 
 for(i=1:FoodNumber)
 system_cost(i,1) = sum_icost + sum(component_cost.*frac(i,:));
 
 system_cost(i,1) = 0.0000265 * edns(Foods(i,:))*8760 + system_cost(i,1);
 
 end
 rev=1+system_cost;
Fitness=1./rev;


%min fitness
BestInd=find(system_cost==min(system_cost));
BestInd=BestInd(end);
GlobalMin=system_cost(BestInd);
GlobalParams=Foods(BestInd,:); 



 %[value index] = min(system_cost);
 
 % particle with index is the best of this swarm
 %gbest_fitness=value;
 %gbest=x(index,:);
 
 % personal best of each particle
 %pbest_fitness = system_cost;
 %pbest = x;


%reset trial counters
trial=zeros(1,FoodNumber);

%/*The best food source is memorized*/
%BestInd=find(ObjVal==min(ObjVal));
%BestInd=BestInd(end);
%GlobalMin=ObjVal(BestInd);
%GlobalParams=Foods(BestInd,:);





iter=1;
while ((iter <= maxCycle)),

%%%%%%%%% EMPLOYED BEE PHASE %%%%%%%%%%%%%%%%%%%%%%%%
    for i=1:(FoodNumber)
        
        %/*The parameter to be changed is determined randomly*/
        Param2Change=fix(rand*D)+1;
        
        %/*A randomly chosen solution is used in producing a mutant solution of the solution i*/
        neighbour=fix(rand*(FoodNumber))+1;
       
        %/*Randomly selected solution must be different from the solution i*/        
            while(neighbour==i)
                neighbour=fix(rand*(FoodNumber))+1;
            end;
        
       sol=Foods(i,:);
       %  /*v_{ij}=x_{ij}+\phi_{ij}*(x_{kj}-x_{ij}) */
       sol(Param2Change)=Foods(i,Param2Change)+(Foods(i,Param2Change)-Foods(neighbour,Param2Change))*(-1+2*rand);
            
       %  /*if generated parameter value is out of boundaries, it is shifted onto the boundaries*/
        ind=find(sol(1,1:2)<lb1);
        sol(ind)=lb1(ind);
        ind=find(sol(1,1:2)>ub1);
        sol(ind)=ub1(ind);
        
        ind=find(sol(1,3:9)<lb2);
        sol(ind)=lb2(ind);
        ind=find(sol(1,3:9)>ub2);
        sol(ind)=ub2(ind);
        
        %evaluate new solution
        %ObjValSol=feval(objfun,sol);
        %FitnessSol=calculateFitness(ObjValSol);
        
        nide=ones(1,D);
        newav = nide - sol;
 
 newfr= newav./sol;
 
 
 newsystem_cost = sum_icost + sum(component_cost.*newfr(1,:));
 total_cost= 0.0000265 * edns(sol)*8760 + newsystem_cost;
 
 revsol=1+total_cost;
 FitnessSol=1./revsol;
 
   
        
        
        
       % /*a greedy selection is applied between the current solution i and its mutant*/
       if (FitnessSol>Fitness(i)) %/*If the mutant solution is better than the current solution i, replace the solution with the mutant and reset the trial counter of solution i*/
            Foods(i,:)=sol;
            Fitness(i)=FitnessSol;
            system_cost(i)=total_cost;
            trial(i)=0;
        else
            trial(i)=trial(i)+1; %/*if the solution i can not be improved, increase its trial counter*/
       end;
         
         
    end;

%%%%%%%%%%%%%%%%%%%%%%%% CalculateProbabilities %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%/* A food source is chosen with the probability which is proportioal to its quality*/
%/*Different schemes can be used to calculate the probability values*/
%/*For example prob(i)=fitness(i)/sum(fitness)*/
%/*or in a way used in the method below prob(i)=a*fitness(i)/max(fitness)+b*/
%/*probability values are calculated by using fitness values and normalized by dividing maximum fitness value*/

prob=(0.9.*Fitness./max(Fitness))+0.1;
  
%%%%%%%%%%%%%%%%%%%%%%%% ONLOOKER BEE PHASE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

i=1;
t=0;
while(t<FoodNumber)
    if(rand<prob(i))
        t=t+1;
        %/*The parameter to be changed is determined randomly*/
        Param2Change=fix(rand*D)+1;
        
        %/*A randomly chosen solution is used in producing a mutant solution of the solution i*/
        neighbour=fix(rand*(FoodNumber))+1;
       
        %/*Randomly selected solution must be different from the solution i*/        
            while(neighbour==i)
                neighbour=fix(rand*(FoodNumber))+1;
            end;
        
       sol=Foods(i,:);
       %  /*v_{ij}=x_{ij}+\phi_{ij}*(x_{kj}-x_{ij}) */
       sol(Param2Change)=Foods(i,Param2Change)+(Foods(i,Param2Change)-Foods(neighbour,Param2Change))*(-1+2*rand);
        
       %  /*if generated parameter value is out of boundaries, it is shifted onto the boundaries*/
       ind=find(sol(1,1:2)<lb1);
        sol(ind)=lb1(ind);
        ind=find(sol(1,1:2)>ub1);
        sol(ind)=ub1(ind);
        
        ind=find(sol(1,3:9)<lb2);
        sol(ind)=lb2(ind);
        ind=find(sol(1,3:9)>ub2);
        sol(ind)=ub2(ind);
        
        %evaluate new solution
       % ObjValSol=feval(objfun,sol);
        %FitnessSol=calculateFitness(ObjValSol);
        
        nide1=ones(1,D);
        newav1 = nide1 - sol;
 
        newfr1= newav1./sol;
 
 
        new1system_cost = sum_icost + sum(component_cost.*newfr1(1,:));
        total_cost1= 0.0000265 * edns(sol)*8760 + new1system_cost;
        
        revsol1=1+total_cost1;
         FitnessSol=1./revsol1;
        
        
        
        
        
        
        
       % /*a greedy selection is applied between the current solution i and its mutant*/
       if (FitnessSol>Fitness(i)) %/*If the mutant solution is better than the current solution i, replace the solution with the mutant and reset the trial counter of solution i*/
            Foods(i,:)=sol;
            Fitness(i)=FitnessSol;
            system_cost(i)=total_cost1;
            trial(i)=0;
        else
            trial(i)=trial(i)+1; %/*if the solution i can not be improved, increase its trial counter*/
       end;
    end;
    
    i=i+1;
    if (i==(FoodNumber)+1) 
        i=1;
    end;   
end; 


%/*The best food source is memorized*/
          ind=find(system_cost==min(system_cost));
         ind=ind(end);
         if (system_cost(ind)<GlobalMin)
         GlobalMin=system_cost(ind);
         GlobalParams=Foods(ind,:);

        
         end;
         
         
%%%%%%%%%%%% SCOUT BEE PHASE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%/*determine the food sources whose trial counter exceeds the "limit" value. 
%In Basic ABC, only one scout is allowed to occur in each cycle*/

ind=find(trial==max(trial));
ind=ind(end);
if (trial(ind)>limit)
    trial(ind)=0;
    sol(1,1:2)=(ub1-lb1).*rand(1,d1)+lb1;
    sol(1,3:9)=(ub2-lb2).*rand(1,d2)+lb2;
    %ObjValSol=feval(objfun,sol);
    %FitnessSol=calculateFitness(ObjValSol);
     nide1=ones(1,D);
        newav1 = nide1 - sol;
 
        newfr1= newav1./sol;
 
 
        new1system_cost = sum_icost + sum(component_cost.*newfr1(1,:));
        total_cost1= 0.0000265 * edns(sol)*8760 + new1system_cost;
        revsol1=1+total_cost1;
         FitnessSol=1./revsol1;
        
    
    
    Foods(ind,:)=sol;
    Fitness(ind)=FitnessSol;
    system_cost(ind)=total_cost1;
end;

Storeglobalparams(iter,:)=GlobalParams;
Storemins(iter)=GlobalMin;
%fprintf('Ýter=%d ObjVal=%g\n',iter,GlobalMin);
iter=iter+1;

end % End of ABC

GlobalMins(r)=GlobalMin;

end; %end of runs

save all

