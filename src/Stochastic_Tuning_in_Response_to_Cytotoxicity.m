num_trials = 5; 
doses = [1,10,50,100,250,500,1000];
all_final_N = zeros(length(doses), num_trials);
%this takes a snapshot of the number of cells remaining after the
%count the total number of cells at the end of
%the simulation for each dose.

for drug_idx = 1:length(doses)
    for t_idx = 1:num_trials
    
    N=50;
    
    D_current = doses(drug_idx);

    kd = 20; %induction midpoint dose


    %assume that all genes are off
    x=zeros(1,N);  % initial state; they are all off; 0/1 boolean array
    RNA=zeros(1,N);%initialize # of mRNA molecules
    prot=zeros(1,N);%initialize # of proteins
    
    t=0;  % initial time
    tmax = 50; %final time
    tracked = 1; %track behavior of one cell

    RNA_tracked = [];
    prot_tracked = [];
    t_tracked = [];
    N_history = N; %list of possiblities to know which cell experienced what
  

    tList = [t];
    mean_RNA_list = [mean(RNA)];
    mean_prot_list = [mean(prot)];
    RNA_list = [sum(RNA)];
    prot_list = [sum(prot)];

    
    
    %make experimental snapshots
    tsample = 5:5:tmax;
    sample_idx=1;
    RNA_sample = cell(1,length(tsample));
    prot_sample = cell(1,length(tsample));

    %parameters
    drug_delta=0.01; %max death rate drug can achieve
    k_on=(drug_delta + 0.2*D_current)/(kd+D_current);
    k_off=100;
    K=0.001;%protein level to beat the effects of the drug
    ktransc = 100; %activation rate of resistance gene
    ktransl = 70; %translation rate of resistance gene 
    division = 0.02; %constant rate of division 

    
    mrna_decay=0.1; %rate to degrade RNA
    prot_degrad=0.5; %prot half-life in denominator; protein degradation
  

    n=2; %slope representing drug response (drug-protein affinity)

    max_steps = 1e6;
    step = 0;

   
    %main simulation while loop
    while t<tmax
            
        if N <= 0
           break;
        elseif N>10000
            break;
        end             
        
    %hill function that represents function of protein and drug
        A = D_current^n;
        B = (K .* prot).^n;   % or K*prot if true matrix multiply intended
        death = drug_delta * (A ./ (A + B));    %A/(A+B) is inhibition 
        
        
        %gene state parameters
        active = k_on * (1-x); %if I multiply by just x, I will always get 0 except for the chosen cell id.
        off = k_off * x;
     
        transc=ktransc .* x;
        transl=ktransl * RNA;
        mrna_death = mrna_decay * RNA; %sum(kdegrad*RNA) in planning notes
        prot_death = prot_degrad * prot;
        div = division * ones(1, N);
            
            
        
        %all possible events
        events = [death,active,off,transc,transl,mrna_death,prot_death,div];
        Rtot=sum(events); %this was every event in the mRNAexpression that affected the state transition
        if Rtot <= 0
             break;  % no more reactions possible
        end

        %pick the event that happened: Gillespie algorithm
        r = rand * Rtot;
        %logic about the incorporation of event_id and cell_id by Google Gemini and mRNAexpression.m
        %once true value is found, it should be passed into the events list to find
        %the reaction and cell it should be linked to 
        event_id = find(cumsum(events)>= r, 1); %find where rate is greater than or equal to the random reaction/event weighted by the rate 
     
        %which cell; always rotating through each reaction for each cell based
        %on the random number that is generated (circular array)
        rxn_type = 1 + floor((event_id-1)/N);
        cell_id  = mod(event_id-1, N) + 1;
        
        %update cell fate
        %update RNA and protein counts based on the selected event

        %death
        if (rxn_type == 1)
            if cell_id == tracked 
                prot_tracked(end+1) = prot(cell_id);
                RNA_tracked(end+1) = RNA(cell_id);
                t_tracked(end+1) = t;
                tracked = NaN;         
            end
            
            %cell died: empty and reduce population by 1
            x(cell_id)=[];
            RNA(cell_id) = [];
            prot(cell_id) = [];
            N = N-1;
        
        
        %gene turned on        
        elseif (rxn_type == 2)
            x(cell_id)=1;
            
        %gene turned off
        elseif (rxn_type == 3)
            x(cell_id)=0;
        
            
        %transcription
        elseif (rxn_type == 4)
            if cell_id >= 1 && cell_id <= numel(RNA)
            RNA(cell_id) = RNA(cell_id) + 1;
            end

        %translation
        elseif (rxn_type == 5)
            if cell_id >= 1 && cell_id <= numel(prot)
            prot(cell_id) = prot(cell_id) + 1;
            end

        %mrna degrades
        elseif (rxn_type == 6)
            if RNA(cell_id)>0
                RNA(cell_id) = RNA(cell_id)-1;
            end
        
       
        %protein degrades
        elseif (rxn_type == 7)
           if prot(cell_id)>0
               prot(cell_id) = prot(cell_id)-1;
           end
           
        
        
        %cell divides
        elseif (rxn_type==8)

        % Parent cell keeps half, daughter cell gets half
            RNA(cell_id) = ceil(RNA(cell_id)/2);
            prot(cell_id) = ceil(prot(cell_id)/2);

        % Append to the end
            x(N+1) = x(cell_id); 
            RNA(N+1) = RNA(cell_id);
            prot(N+1) = prot(cell_id);
            N = N + 1;

        end
       
    %time update
       dt = -log(rand())/Rtot;
        t = t + dt;
        tList(end+1) = t;
        RNA_list(end+1) = sum(RNA);
        prot_list(end+1) = sum(prot);
        mean_RNA_list(end+1) = mean(RNA);
        mean_prot_list(end+1) = mean(prot);
        N_history(end + 1) = N;
        
        % Record tracking
        if ~isnan(tracked) && tracked <= N && cell_id==tracked
            RNA_tracked(end+1) = RNA(tracked);
            prot_tracked(end+1) = prot(tracked);
            t_tracked(end+1) = t;
        end
        
        %Experimental sample collection
        while sample_idx <= length(tsample) && t>=tsample(sample_idx)
            if ~isnan(tracked) && tracked <= N
                RNA_sample{sample_idx} = RNA(tracked);
                prot_sample{sample_idx} = prot(tracked);
            else 
                RNA_sample{sample_idx} = NaN;
                prot_sample{sample_idx} = NaN;
            end
                
                sample_idx = sample_idx + 1;
        end
        
      
        step = step + 1;
        if step > max_steps
            break;
        end
    end

        % dose response storage
        all_final_N(drug_idx, t_idx) = N;
    end
      
end 



RNA_all = RNA_list;   

mu_rna = mean(RNA_all);
var_rna = var(RNA_all);

Fano_rna = var_rna / mu_rna;

disp(['Mean RNA = ', num2str(mu_rna)])
disp(['Variance RNA = ', num2str(var_rna)])
disp(['Fano factor for RNA = ', num2str(Fano_rna)])
data = RNA_all;



% population RNA time series
figure(1);
% histogram of simulation
histogram(data, 'Normalization', 'pdf');
hold on;

% Poisson fit using same mean
lambda = mean(data);
x = 0:max(data);

poisson_pdf = poisspdf(x, lambda);


plot(x, poisson_pdf, 'g-', 'LineWidth', 2);

xlabel('RNA level')
ylabel('Probability')
title('RNA distribution: simulation vs Poisson')
legend('Simulation', 'Poisson Expectation')
hold off

% RNA time series
figure(2);
plot(tList, RNA_list);
xlabel('Time')
ylabel('RNA')
title('Sum RNA vs. time')

% protein time series
figure(3);
plot(tList, prot_list);
xlabel('Time')
ylabel('Protein level')
title('Sum Protein vs time')

figure(4);
plot(t_tracked,RNA_tracked)
xlabel('Time')
ylabel('RNA')
title('Single-Cell RNA')


figure(5);
plot(t_tracked, prot_tracked)
xlabel('Time')
ylabel('Protein')
title('Single-Cell Protein')

figure(6);
plot(tList,mean_RNA_list)
xlabel('Time')
ylabel('mRNA number')
title('Average RNA Time Series')

figure(7);
plot(tList, mean_prot_list);
xlabel('Time')
ylabel('Protein')
title('Average Protein Time Series')

figure(8);
avg_N = mean(all_final_N, 2);
semilogx(doses, avg_N, '-ob', 'LineWidth', 2);
xlabel('Drug dose');
ylabel('Final cell count (N)');
title('Dose-response');
grid on;


%Help in generating curve from Gemini
% Calculate death rate based on  Hill Equation
ratio = logspace(-1, 3, 100);
n_val = 1; %current Hill coefficient
delta = 0.1; %drug_delta
theoretical_death = delta * (ratio.^n_val ./ (ratio.^n_val + 1));
grid on;


figure(9);
K_values = [1, 5, 10, 50];  % different K's

for i = 1:length(K_values)
    K = K_values(i);
    
    theoretical_death = ratio ./ (K + ratio);  % adjust if needed
    
    semilogx(ratio, theoretical_death, 'LineWidth', 2);
end

xlabel('Drug Ratio (D / K*P)');
ylabel('Death Rate');
title('Effect of K on Sigmoid Death Curve');

legend('K=1','K=5','K=10','K=50');

D = logspace(-1,2,100);

k = 5;                 % protein scaling
K0 = 10;
beta = 20;

P = D ./ k;            % protein = D/k

KP = K0 + beta * P;    % resistance increases with protein

death = 1 ./ (1 + D ./ KP);   % DECREASING sigmoid



% histogram of  simulation
prot_all = prot_list;   % population protein time series

mu_prot = mean(prot_all);
var_prot = var(prot_all);
Fano_prot = var_prot / mu_prot;

disp(['Mean Protein = ', num2str(mu_prot)])
disp(['Variance Protein = ', num2str(var_prot)])
disp(['Fano factor for Protein = ', num2str(Fano_prot)])
figure(11);

data = prot_all;
histogram(data, 'Normalization', 'pdf');
hold on;

% Poisson fit using same mean
lambda = mean(data);
x = 0:max(data);

poisson_pdf = poisspdf(x, lambda);


plot(x, poisson_pdf, 'r-', 'LineWidth', 2);

xlabel('Protein level');
ylabel('Probability');
title('Protein distribution: simulation vs Poisson');
legend('Simulation', 'Poisson expectation');
hold off;

%Cell survival time series
figure(12);
plot(tList, N_history);
xlabel('Time');
ylabel('Cell Count');
title('Cell Survival')

%Theoretical drug-protein-resistance dynamic
figure(13);
final_death_rates = drug_delta * (D_current^n ./ (D_current^n + (K .* prot).^n));
figure;
scatter(prot, final_death_rates, 'filled', 'MarkerFaceAlpha', 0.5);
xlabel('Protein Level in Surviving Cells');
ylabel('Death Rate');
title('Resistance Status of Surviving Cells');