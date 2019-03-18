%Rebuttal Nature Microbiology
%Reviewer 1
%Robustness analysis
%Cristal Zuniga
%1.16.2019
%
%
clear all
initCobraToolbox
load('ModelsSubmission.mat', 'CoModels')
load('ModelsSubmission.mat', 'ModelPA')
load('ModelsSubmission.mat', 'mergeModelm')
load('ModelsSubmission.mat', 'model')
load('ModelsSubmission.mat', 'modelSc2')

BOF=find(CoModels{1,1}.c)
 
 solutionCo={};
for j=1:6
    solutionCo{j}=optimizeCbModel(CoModels{j});
    solutionCo{1,j}.f
    TableGr(3,j)=solutionCo{1,j}.x(BOF(2))
    TableGr(4,j)=solutionCo{1,j}.x(BOF(1))
end  

modelSc2.csense={};
for i=1:length(modelSc2.mets)
    modelSc2.csense{i}='E';
end
modelSc2.csense=char(modelSc2.csense);
solutionY=optimizeCbModel(modelSc2);

solutionA={};
for j=1:6
    solutionA{j}=optimizeCbModel(ModelPA{1,j});
    TableGr(1,j)=solutionA{1,j}.f
    TableGr(2,j)=solutionY.f
end  

%Performing robustness analysis while variyng the growth rate of
%microorganism 1

SP=[];contf=1;contt=1;solutions={};
for n=1:numel(CoModels) %models

    for m=1:2 %biomass variables

       model=CoModels{n};
       BOFs=find(model.c)
       units={' uptake rate (mmol gDW/h)'};
       solutionLna = optimizeCbModel(model);
       Source=solutionLna.x(BOFs);
       Reactions=model.rxns(BOFs);
       Enzymes2=[0.0:0.01:1];
       mayaAuto=[];x=[];mayaAuto1=[];mayaAuto2=[];growthRatio=[];SP=[];
       contn=1;
        for j=1:numel(Enzymes2)
            %for i=1:4
                %valor=i-1;
                model2=changeRxnBounds(model,Reactions{m},Source(m)*Enzymes2(j),'b'); 
                model2.c(BOFs(m))=0;
                Source(m)*Enzymes2(j);
                solutionLna = optimizeCbModel(model2);
                
                %Extracting data from the flux distribution vector
                if solutionLna.f>0.000001    
                    solutions=horzcat(solutions,solutionLna.x);
                     mayaAuto(contn)=solutionLna.f;
                        if mayaAuto(contn)~=0
                         mayaAuto1(contn)=solutionLna.x(BOFs(1,1));
                         mayaAuto2(contn)=solutionLna.x(BOFs(2,1));
                        else
                            if m==1
                            mayaAuto1(contn)=Source(m)*Enzymes2(j);
                            mayaAuto2(contn)=0;
                            else
                            mayaAuto1(contn)=0;
                            mayaAuto2(contn)=Source(m)*Enzymes2(j);  
                            end
                         %j
                        end
                    mayaAuto(contn)=solutionLna.f+Source(m);
                    x(contn)=Enzymes2(j);
                    %Extracting data about shadow prices
                    if m==1
                        growthRatio(contn)=mayaAuto2(contn)./mayaAuto1(contn);
                    else
                        growthRatio(contn)=mayaAuto1(contn)./mayaAuto2(contn);
                    end
                    SP=horzcat(SP,solutionLna.y);

                    %Extractind metabolic exchange
                    [SPM_rxns,SPM_D1, Un_D, Un_A] = sharedPoolCM(model2,1);
                    sPool{contn,1}=SPM_D1;
                    sPool{contn,2}=Un_D;
                    sPool{contn,3}=Un_A;

                    contn=contn+1;
                end      
        end

        subplot(numel(CoModels),4,contf)
        plot(x,mayaAuto)
        hold on
        plot(x,mayaAuto1)
        plot(x,mayaAuto2)
        hold off
        ylabel('Growth rate (1/h)'); 
        if m==1
        xlabel('Mycobiont fraction'); 
        else
        xlabel('Photobiont fraction');     
        end
        contf=contf+1;
        subplot(numel(CoModels),4,contf)
        plot(growthRatio,'o')
        ylabel('Growth ratio')
        xlabel('Number of essays')
        axis([0 105 -5 60])
        contf=contf+1

       % mayaAuto(:)='';x(:)='';mayaAuto1(:)='';mayaAuto2(:)='';growthRatio='';
        
        ResultsSP{contt,1}=growthRatio;
        ResultsSP{contt,2}=SP;
        ResultsSP{contt,3}=sPool;
        contt=contt+1;
        mayaAuto(:)='';x(:)='';mayaAuto1(:)='';mayaAuto2(:)='';growthRatio=''; %SP='';
    end
end
%%
%Evaluating the shadow prices results
cont=1;n=15;
figure
for i=1:2:12
    A=ResultsSP{i,2};
    A1=A(:,101);
    A1(A1 == 0) = 1;
    A=A./A1;
    B=ResultsSP{i+1,2};
    B1=B(:,101);
    B1(B1 == 0) = 1;
    B=B./B1;
    
    pvalues = mattest(A, B);
    STD=mean(horzcat(std(A,0,2),std(B,0,2)),2);
    
    x=STD./mean(horzcat(mean(A,2),mean(B,2)),2);
    
subplot(3,2,cont)
scatter(x,-log10(pvalues))
a=find(x>n);
a1=find(x<-n);
b=find(-log10(pvalues)>n+5);
c=vertcat(a,b,a1);
%text(x(c),-log10(pvalues(c)),model2.mets(c))
axis([-100 100 -10 300])
box on
%breakyaxis([50 230])
cont=cont+1;
end
% 
%Robustness Figure 2
figure
cont=1;n=15;
i=1
pvalues = mattest(ResultsSP{i,2}, ResultsSP{i+1,2});
STD=mean(horzcat(std(ResultsSP{i,2},0,2),std(ResultsSP{i+1,2},0,2)),2);
x=STD./mean(horzcat(mean(ResultsSP{i,2},2),mean(ResultsSP{i+1,2},2)),2);
scatter(x,-log10(pvalues))
a=find(x>n);
a1=find(x<-n);
b=find(-log10(pvalues)>n+5);
c=vertcat(a,b,a1);
text(x(c),-log10(pvalues(c)),model2.mets(c))
axis([-100 100 -10 300])
box on
xlabel('Fold change (%)')
ylabel('-log_1_0 (p-value)')
breakyaxis([50 230])

%%
%Evaluating the SMP

sharedP={};unD={};unA={};sharedPoolunD={}; sharedPoolunA={};sharedPool={};
for i=1:length(ResultsSP)
    for n=1:length(ResultsSP{i,3})
        
     SPM_D1=ResultsSP{i,3}{n,1};
     Un_D=ResultsSP{i,3}{n,2};
     Un_A=ResultsSP{i,3}{n,3};
        
    for j=1:numel(SPM_D1.Name)
        a=cellstr(table2array(SPM_D1.Name(j)));
        a=a{1,1}(1:3);
        b=cellstr(table2array(SPM_D1.Full_Name(j)));
        b=b{1,1}(1:3);
        
        if isequal(a,b)
            sharedP{j,1}=table2array(SPM_D1.Name(j));
            sharedP{j,2}=table2array(SPM_D1.Full_Name(j));
            sharedP{j,3}=SPM_D1.Rxn(j);
            sharedP{j,4}=SPM_D1.Flux(j);
        end
        if isequal(b,'D__')
            sharedP{j,1}=table2array(SPM_D1.Name(j));
            sharedP{j,2}=table2array(SPM_D1.Full_Name(j));
            sharedP{j,3}=SPM_D1.Rxn(j);
            sharedP{j,4}=SPM_D1.Flux(j);
        end
    end
 
    for j=1:numel(Un_D.Name)       
    unD{j,1}=table2array(Un_D.Name(j));
    unD{j,2}=Un_D.Flux(j);     
    end
    [a ia]=intersect(unD(:,1),sharedP(:,2));
    unD(ia,:)=[];
    
    
    for j=1:numel(Un_A.Name) 
    unA{j,1}=table2array(Un_A.Name(j));
     unA{j,2}=Un_A.Flux(j);       
    end
    
    sharedPool=vertcat(sharedPool,sharedP(:,1));
    sharedPoolunD=vertcat(sharedPoolunD,unD(:,1));
    sharedPoolunA=vertcat(sharedPoolunA,unA(:,1)); 
    end
end

C=vertcat(unique(sharedPool))
unique(sharedPoolunD)
unique(sharedPoolunA);



cont=2;
for i=1:numel(C)
    for j=1:length(ResultsSP)
       for n=1:length(ResultsSP{j,3}) 
         SPM_D1=ResultsSP{j,3}{n,1};
           Ind2=strmatch(C{i,1},table2array(ResultsSP{j,3}{n,1}(:,1)),'exact');   
           if ~isempty(Ind2) && numel(Ind2)==1
               C{i,cont}=table2array(ResultsSP{j,3}{n,1}(Ind2,4));
               cont=cont+1;
           else
               C{i,cont}=0;
               cont=cont+1;
           end
           if numel(Ind2)>1
               C{i,cont}=table2array(ResultsSP{j,3}{n,1}(Ind2(1),4));
               cont=cont+1;
           end
       end
    end
       cont=2;
       C{i,1}=C{i,1}(1:end-3);
end

D=cell2mat(C(:,2:1213));
D2=cell2mat(C(:,2:1213))./cell2mat(C(:,1213));
B=find(abs(D2)>100);
D2(B)=100;

HeatMap(D2,'RowLabels',C(:,1))
hold on
H=addTitle(H,'Robustness while varying the photobiont')



%%
%Permutation analysis
%
%
initCobraToolbox
%Loading individual and community models
load('C:\Users\Cristal\Documents\MATLAB\Chlorella1\Community\Chlorella\iCZ-Cv(843)-Sc(904).mat')
load('combinedModel_GeneAnalysis8.11.2017_ReanalysisFail.mat', 'modelSc2')
load('iCZ-CvSc_BOF_GeneEssentiallity_glucoseExpData6.20.2017.mat', 'CoModels')
load('iCZ-CvSc_BOF_GeneEssentiallity_glucoseExpData6.20.2017.mat', 'ModelPA')
%initCobraToolbox;
Genes={'WT'; 'YBR153W';'YBR029C';'YER003C';'YFL045C'; 'YDL055C'; 'YHR007C'; 'YPR035W'; 'YAL012W'; 'YKL055C'};
Genes2={'(WT)'; '(RIB7^-)' ;'(CDS1^-)';'(PMI40^-)';'(SEC53^-)';'(PSA1^-)';'(ERG11^-)';'(GLN1^-)';'(CYS3^-)';'(OAR1^-)'};
ChlorellaError=[0.004600	0.0050; 0.004100 0.0032500;0.00110 0.00110;0.001590	0.00187;0.002700	0.005560;0.002700	0.00180;0.003910	0.001370;0.004870	0.001370; 0.00230	0.00186; 0.0120	0.00186];
Chlorella=[0.031000	0.033;0.01800	0.024850;0.00110 0.00110;0.01490	0.022260;0.0378000	0.022940;0.003400	0.023480;0.01930	0.0246900;0.02140	0.024690;0.021000	0.025870;0.03090	0.025870];
Yeast=[0.0330	0.0150; 0.019130	0.00110;0.00110 0.00110; 0.012970	0.004350; 0.01600	0.001680; 0.013390	0.0035700; 0.021620	0.0094400; 0.021620	0.001710; 0.022540	0.023180; 0.022540	0.015830];
YeastError=[0.000120	0.001400; 0.003380	0.000250;0.000110 0.000110;0.001120	0.00060; 0.002170	0.00029; 0.001890	0.00084; 0.002620	0.00046; 0.002340	0.00634; 0.002240	0.00186; 0.007820	0.00162];
[C_grRatio,C_grRatio_coY,C_grRatio_coA,C_grRateKO,C_grRateWT,C_hasEffect,C_delRxns,~,~,C_fluxsolution1,C_modelos] = singleGeneDeletionCM(mergeModelm);
[D_grRatio,D_grRateKO,D_grRateWT,D_delRxns,D_hasEffect,D_fluxSolution] = singleGeneDeletion(modelSc2);
figure
mycolor=[0.9583    0.8672    0.6992;0    0.5    0.0];
Table={};
for i=1:numel(Genes)
subplot(5,2,i)
H=bar(1:2,Yeast(i,:),0.5)
H(1).FaceColor = mycolor(1,:);
hold on
errorbar(1:2,Yeast(i,:),YeastError(i,:),'k.')
BOF=find(mergeModelm.c);
if i==1
FBA=optimizeCbModel(mergeModelm);
plot(1,FBA.x(BOF(1)),'o','MarkerSize',12,'MarkerFaceColor',[0.9542    0.6406    0.3750]);
Table{i,2}=FBA.x(BOF(1));
else
Ind=strmatch(Genes{i},mergeModelm.genes,'exact');
FBA=optimizeCbModel(C_modelos{Ind});
plot(1,FBA.x(BOF(1)),'o','MarkerSize',12,'MarkerFaceColor',[0.9542    0.6406    0.3750]);
Table{i,2}=FBA.x(BOF(1));
end
if i==1
FBA=optimizeCbModel(modelSc2);
plot(2,FBA.f,'o','MarkerSize',12,'MarkerFaceColor',[0.9542    0.6406    0.3750]);
Table{i,3}=FBA.f;
else
Ind=strmatch(Genes{i},modelSc2.genes,'exact');
plot(2,D_grRateKO(Ind),'o','MarkerSize',12,'MarkerFaceColor',[0.9542    0.6406    0.3750]);
Table{i,3}=D_grRateKO(Ind);
end
L=strcat(Genes{i,1},Genes2{i,1});
title(L)
c = {'KO','WT'};
set(gca, 'XTick', 1:2, 'XTickLabel', c);
axis([0.25 2.75 0 0.04])
ylabel('Growth rate (1/h)');

Table{i,1}=Genes{i};


end

%Permutation test.
%Cristal Zuniga 
%1.25.2019

A=[0.0330000000000000;0.0225000000000000;0.00200000000000000;0.0191000000000000;0.0134000000000000;0.0130000000000000;0.0160000000000000;0.0216000000000000;0.0225000000000000;0.0216000000000000;0.0150000000000000;0.0232000000000000;0.000400000000000000;0.00110000000000000;0.00360000000000000;0.00440000000000000;0.00170000000000000;0.00940000000000000;0.0158000000000000;0.00170000000000000];
B=[0.0320000000000000;0.0320000000000000;0.0320000000000000;0.0320000000000000;0.0320000000000000;0.0320000000000000;0.0320000000000000;0.0320000000000000;0.0320000000000000;0.0320000000000000;0.0150000000000000;0.0150000000000000;0;0;0;0;0;0;0.0150000000000000;0];

for j=1:100000
    C=randperm(20);
    B=B(C);
    Result={};
    for i=1:length(A)
        if A(i)>=0.013 && A(i)<=0.033 && B(i)>=0.015
            Result{i,1}='TP';
        elseif A(i)<=0.0020 && B(i)>0.015
            Result{i,1}='FP';
        elseif A(i)>=0.013 && B(i)<0.015
            Result{i,1}='FN';
        end
        if A(i)>=0.0004 && A(i)<=0.0044 && B(i)==0
            Result{i,1}='TN';
        elseif A(i)>=0.0004 && A(i)<=0.0094 && B(i)<0.015
            Result{i,1}='FN';
        elseif A(i)<0.013 && B(i)>=0.015
            Result{i,1}='FP';    
        end
    end

    TP=numel(strmatch('TP', Result, 'exact'));
    TN=numel(strmatch('TN', Result, 'exact'));
    FP=numel(strmatch('FP', Result, 'exact'));
    FN=numel(strmatch('FN', Result, 'exact'));

  Data(j,1)=(TP+TN)./(TP+TN+FN+FP);  
  Data(j,2)=TP./(TP+FN);
  Data(j,3)=TN./(TN+FP);
  Data(j,4)=TP./(TP+FP);
  Data(j,5)=TN./(TN+FN);
  Data(j,6)=((TP*TN)-(FP*FN))./sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN));
end

lines=[0.90;0.92;0.88;0.92;0.88;0.79];
titles={'Accuracy';'Sensitivity';'Specificity';'Possitive predicted value';'Negative predicted value';'Matthews correlation coefficient'};
for i=1:6
   subplot(3,2,i)
   hist(Data(:,i));
   xlabel(strcat(titles{i},' distribution'));
   ylabel('Frequency')
   line([lines(i) lines(i)],[0 4000]);
end



%Oxygen reactions screening

B=strmatch('o2_',CoModels{1,1}.mets)


CoModels{1,1}.mets(B)
%Reactions consuming or producing oxygen
numel(find(CoModels{1,1}.S(B,:))) %277
%Reactions in mycobiont consuming o2
numel(find(CoModels{1,1}.S(B(1:4),:)<0)) %71
%Reactions in mycobiont producing o2
numel(find(CoModels{1,1}.S(B(1:4),:)>0)) %5
%Reactions in photobiont consuming o2
numel(find(CoModels{1,1}.S(B(5:end),:)<0)) %194
%Reactions in photobiont producing o2
numel(find(CoModels{1,1}.S(B(5:end),:)>0)) %7


%%
mayaAuto=[];x=[];y=[];mayaAuto1=[];mayaAuto2=[];mayaAuto1o2A=[];mayaAuto1o2D=[];mayaAuto2co2A=[];mayaAuto2co2D=[];mayaAuto2glcA=[];mayaAuto2glcD=[];

model=mergeModelm;
BOFs=find(model.c);
o2(1,1)=strmatch('O2t_D',model.rxns,'exact');
o2(2,1)=strmatch('O2t_A',model.rxns,'exact');
co2(1,1)=strmatch('CO2t_D',model.rxns,'exact');
co2(2,1)=strmatch('CO2t_A',model.rxns,'exact');
glc(1,1)=strmatch('GLCt1_D',model.rxns,'exact');
glc(2,1)=strmatch('GLCt2_A',model.rxns,'exact');
    solutionLna = optimizeCbModel(model,'max','one');
    cont=0;cont0=0;
    for i=1:10
        cont0=cont+1;
        for j=1:10
            valor=(i-1)/5;
            valor1=(j-1)/5;
            model=changeRxnBounds(model,{'EX_glc__A_[smp]','EX_co2_[smp]'},[-valor -valor1],'l');
            solutionLna = optimizeCbModel(model,'max','one');
                if solutionLna.f~=0
                  mayaAuto(i,j)=solutionLna.f;  
                    cont=cont+1;
                mayaAuto1(i,j)=solutionLna.x(BOFs(1,1));
                mayaAuto2(i,j)=solutionLna.x(BOFs(2,1));
                mayaAuto1o2(i,j)=solutionLna.x(o2(1,1))./solutionLna.x(o2(2,1));
                mayaAuto2co2(i,j)=abs(solutionLna.x(co2(1,1)))./solutionLna.x(co2(2,1));
                mayaAuto2glc(i,j)=solutionLna.x(glc(1,1))./solutionLna.x(glc(2,1));
%             mayaAuto1o2D(i,j)=solutionLna.x(o2(1,1));
%             mayaAuto1o2A(i,j)=solutionLna.x(o2(2,1));
%             mayaAuto2co2D(i,j)=solutionLna.x(co2(1,1));
%             mayaAuto2co2A(i,j)=solutionLna.x(co2(2,1));
%             mayaAuto2glcD(i,j)=solutionLna.x(glc(1,1));
%             mayaAuto2glcA(i,j)=solutionLna.x(glc(2,1));
                end
        x(i)=valor;
        y(j)=valor1;
        end
    cont=0;
    end
 figure
ax2 = subplot(2,3,1);
contourf(x,y,mayaAuto1,5)
colormap(ax2,hot(8))
 title(ax2,'Mycobiont growth rate in the lichen'); 
 ylabel('O_2 uptake rate (mmol gDW/h)');
 xlabel('Glucose uptake rate (mmol gDW/h)');
 colorbar
 
  ax2 = subplot(2,3,2);
contourf(x,y,mayaAuto2,5)
colormap(ax2,hot(8))
 title(ax2,'Photobiont growth rate in the lichen'); 
 ylabel('O_2 uptake rate (mmol gDW/h)');
 xlabel('Glucose uptake rate (mmol gDW/h)');
 colorbar
 
  ax2 = subplot(2,3,4);
contourf(x,y,log10(mayaAuto1o2),5)
% contourf(x,y,mayaAuto1o2,5)
colormap(ax2,hot(8))
 title(ax2,'O2 ratio'); 
 ylabel('O_2 uptake rate (mmol gDW/h)');
 xlabel('Glucose uptake rate (mmol gDW/h)');
 colorbar
  ax2 = subplot(2,3,5);
contourf(x,y,log10(mayaAuto2co2),5)
% contourf(x,y,mayaAuto2co2,5)
colormap(ax2,hot(8))
 title(ax2,'co2 ratio'); 
 ylabel('O_2 uptake rate (mmol gDW/h)');
 xlabel('Glucose uptake rate (mmol gDW/h)');
 colorbar
  ax2 = subplot(2,3,6);
contourf(x,y,log10(abs(mayaAuto2glc)),5)
% contourf(x,y,mayaAuto2glc,5)
colormap(ax2,hot(8))
 title(ax2,'Glucose ratio'); 
 ylabel('O_2 uptake rate (mmol gDW/h)');
 xlabel('Glucose uptake rate (mmol gDW/h)');
 colorbar
 
 %%
    
    figure
ax2 = subplot(2,4,1);
contourf(x,y,mayaAuto1,5)
colormap(ax2,hot(8))
 title(ax2,'Mycobiont growth rate in the lichen'); 
 ylabel('O_2 uptake rate (mmol gDW/h)');
 xlabel('Glucose uptake rate (mmol gDW/h)');
 colorbar

 ax3 = subplot(2,4,2);
contourf(x,y,mayaAuto1o2D,5)
colormap(ax3,hot(8))
xlabel('Glucose uptake rate (mmol gDW/h)');
 ylabel('O_2 uptake rate (mmol gDW/h)');
 title(ax3,'Oxygen mycobiont'); 
  colorbar
 
 ax1 = subplot(2,4,3);
contourf(x,y,mayaAuto2co2D,5)
colormap(ax1,hot(8))
 ylabel('O_2 uptake rate (mmol gDW/h)');
 xlabel('Glucose uptake rate (mmol gDW/h)');
title(ax1,'CO2 mycobiont');
 colorbar
 
 ax1 = subplot(2,4,4);
contourf(x,y,mayaAuto2glcD,5)
colormap(ax1,hot(8))
 ylabel('O_2 uptake rate (mmol gDW/h)');
 xlabel('Glucose uptake rate (mmol gDW/h)');
title(ax1,'Glucose mycobiont');
 colorbar

 ax2 = subplot(2,4,5);
contourf(x,y,mayaAuto2,5)
colormap(ax2,hot(8))
 title(ax2,'Photobiont growth rate in the lichen'); 
 ylabel('O_2 uptake rate (mmol gDW/h)');
 xlabel('Glucose uptake rate (mmol gDW/h)');
 colorbar

 ax3 = subplot(2,4,6);
contourf(x,y,mayaAuto1o2A,5)
colormap(ax3,hot(8))
xlabel('Glucose uptake rate (mmol gDW/h)');
 ylabel('O_2 uptake rate (mmol gDW/h)');
 title(ax3,'Oxygen photobiont'); 
  colorbar
 
 ax1 = subplot(2,4,7);
contourf(x,y,mayaAuto2co2A,5)
colormap(ax1,hot(8))
 ylabel('O_2 uptake rate (mmol gDW/h)');
 xlabel('Glucose uptake rate (mmol gDW/h)');
title(ax1,'CO2 photobiont');
 colorbar
 
 ax1 = subplot(2,4,8);
contourf(x,y,mayaAuto2glcA,5)
colormap(ax1,hot(8))
 ylabel('O_2 uptake rate (mmol gDW/h)');
 xlabel('Glucose uptake rate (mmol gDW/h)');
title(ax1,'Glucose photobiont');
 colorbar
 
 %%
  %Pair-wise robustness analysis
 model=mergeModelm;
BOFs=find(model.c);
 units={' uptake rate (mmol gDW/h)'};
 solutionLna = optimizeCbModel(model,'max','one');
 Enzymes2={'EX_glc__A_[smp]';'EX_no3_[smp]';'EX_o2_[smp]';'EX_nh4_[smp]'};
Enzymes={'Glucose_ ';'NO_3 ';'O_2 ';'NH_4 '};
mayaAuto=[];x=[];y=[];mayaAuto1=[];mayaAuto2=[];mayaAuto1o2A=[];mayaAuto1o2D=[];mayaAuto2co2A=[];mayaAuto2co2D=[];mayaAuto2glcA=[];mayaAuto2glcD=[];mayaAuto1o2=[];mayaAuto2co2=[];mayaAuto2glc=[];
% o2(1,1)=strmatch('O2t_A',model.rxns,'exact');
% o2(2,1)=strmatch('O2t_D',model.rxns,'exact');
% co2(1,1)=strmatch('CO2t_A',model.rxns,'exact');
% co2(2,1)=strmatch('CO2t_D',model.rxns,'exact');
% glc(1,1)=strmatch('GLCt2_A',model.rxns,'exact');
% glc(2,1)=strmatch('GLCt1_D',model.rxns,'exact');
figure

for j=1:numel(Enzymes2)
    model=CoModels{1,6};
    solutionLna = optimizeCbModel(model,'max','one');
    for i=1:20
        valor=i-1;
    model=changeRxnBounds(model,Enzymes2{j},-valor,'l'); 
    solutionLna = optimizeCbModel(model,'max','one');
   Modelos{i}=model;
    mayaAuto(i)=solutionLna.f;
        if mayaAuto(i)~=0
            mayaAuto1(i)=solutionLna.x(BOFs(1,1));
            mayaAuto2(i)=solutionLna.x(BOFs(2,1));
             mayaAuto1o2(i,j)=solutionLna.x(o2(1,1))./solutionLna.x(o2(2,1));
                mayaAuto2co2(i,j)=abs(solutionLna.x(co2(1,1)))./solutionLna.x(co2(2,1));
                mayaAuto2glc(i,j)=solutionLna.x(glc(1,1))./solutionLna.x(glc(2,1));
%             mayaAuto1o2A(i,j)=solutionLna.x(o2(1,1));
%             mayaAuto1o2D(i,j)=solutionLna.x(o2(2,1));
%             mayaAuto2co2A(i,j)=solutionLna.x(co2(1,1))
%             mayaAuto2co2D(i,j)=solutionLna.x(co2(2,1));
%             mayaAuto2glcA(i,j)=solutionLna.x(co2(1,1))
%             mayaAuto2glcD(i,j)=solutionLna.x(co2(2,1));
        end
    x(i)=valor;
    end
    subplot(numel(Enzymes2)/2,2,j)
    plot(x,mayaAuto)
    hold on
    plot(x,mayaAuto1)
    plot(x,mayaAuto2)
    plot(x,log10(mayaAuto1o2))
    plot(x,log10(mayaAuto2co2))
    plot(x,log10(mayaAuto2glc))
    hold off
   axis([0.0 3 -0.000 0.15])
    xlabel(strcat(Enzymes{j},units));
    ylabel('Growth rate (h^-^1)');  
    mayaAuto(:)='';x(:)='';mayaAuto1(:)='';mayaAuto2(:)='';mayaAuto1o2(:)='';mayaAuto2co2(:)='';mayaAuto2glc(:)='';
end



sharedP={};cont=1;unD={};unA={};SPM_D=SPM_D1;cont2=1;cont3=0;
for i=1:numel(Modelos)+1
     cont3=0;
    [SPM_rxns,SPM_D1, Un_D, Un_A] = sharedPoolCM(Modelos{i},1);
    for j=1:numel(SPM_D1.Name)
        a=cellstr(table2array(SPM_D1.Name(j)));
        a=a{1,1}(1:3);
        b=cellstr(table2array(SPM_D1.Full_Name(j)));
        b=b{1,1}(1:3);
        
        if isequal(a,b)
           cont3=cont3+1;
    sharedP{cont3,cont}=table2array(SPM_D1.Name(j));
    sharedP{cont3,cont+1}=table2array(SPM_D1.Full_Name(j));
    sharedP{cont3,cont+2}=SPM_D1.Rxn(j);
    sharedP{cont3,cont+3}=SPM_D1.Flux(j);
        end
        if isequal(b,'D__')
            cont3=cont3+1;
    sharedP{cont3,cont}=table2array(SPM_D1.Name(j));
    sharedP{cont3,cont+1}=table2array(SPM_D1.Full_Name(j));
    sharedP{cont3,cont+2}=SPM_D1.Rxn(j);
    sharedP{cont3,cont+3}=SPM_D1.Flux(j);
        end
    end

    for j=1:numel(Un_D.Name)       
    unD{j,cont2}=table2array(Un_D.Name(j));
    unD{j,cont2+1}=Un_D.Flux(j);     
    end
    
    for j=1:numel(Un_A.Name) 
    unA{j,cont2}=table2array(Un_A.Name(j));
     unA{j,cont2+1}=Un_A.Flux(j);       
    end
    cont=cont+4;
    cont2=cont2+2;
end

%Creating the vector of unique shared metabolites by KO
%sharedP(:,1)(~cellfun('isempty',R)) 
C = vertcat(sharedP(:,1),sharedP(:,5),sharedP(:,9),sharedP(:,13),sharedP(:,17),sharedP(:,21),sharedP(:,25),sharedP(:,29), sharedP(:,33))
C=C(~cellfun('isempty',C)) 
C=unique(C)
Val=numel(C)
%Looking fot the specific KO data

cont=1;Metabolites={};
for i=1:Val 
%      Ind2=strmatch(C{i,1},modelCo.mets,'exact');   
%      Metabolites{end+1,1}=modelCo.metNames{Ind2,1};
    for j=[1:2:2*numel(Modelos)]   
        if cont<(2*numel(Modelos))-3
        A=sharedP(:,cont);A=A(~cellfun('isempty',A));%A=unique(A);
        Ind=strmatch(C{i,1},A,'exact');    
                aux=j+1;
                    if ~isempty(Ind) 
                        if numel(Ind)==1                            
                        C{i,aux}=cell2mat(sharedP(Ind,cont+2));
                        C{i,aux+1}=cell2mat(sharedP(Ind,cont+3));
                        %j=j+1;
                        else
                        C{i,aux}=cell2mat(sharedP(Ind(1),cont+2));
                        C{i,aux+1}=cell2mat(sharedP(Ind(1),cont+3));
                        end
                    else
                        C{i,aux}=0;
                        C{i,aux+1}=0;
                        %j=j+1;        
                        
                        
                    end
                    cont=cont+4;  
          else
          break
        end
    end
cont=1;
end

SPMss={};
for i=1:numel(C)/19
    String=C{i,1};
    SPMss{i,1}=upper(String(1:end-1));
    SPMss{i,1}=strrep(SPMss{i,1}, 'T2R_', '');
    SPMss{i,1}=strrep(SPMss{i,1}, 'T_', '');
    %SPMss{i,1}=strrep(SPMss{i,1}, 'TR', '');
    SPMss{i,1}=strrep(SPMss{i,1}, 'TI', '');
end
aux=cell2mat(C(:,2:41));
%Normalization
aux(aux>0.5)=0.5;
aux(aux<-0.5)=-0.5;

j=[2:2:1*41] 
j1=[3:2:1*40] 
aux(:,j)
cg=clustergram(aux(:,j),'RowLabels',SPMss(1:end),...
                                 'Colormap',redbluecmap)
   
                             
 cg=clustergram(aux(:,j1),'RowLabels',SPMss(1:end),...
                                 'Colormap',redbluecmap)
                            
 cg=clustergram(aux,'RowLabels',SPMss(1:end),...
                                 'Colormap',redbluecmap)
                              
                             
                       
                             
 %%
 %Shadow prices
    
 Met_tar={'ribflv_c'...%ok experimentally
'trp__L_c'...
'tyr__L_c'...
'val__L_c'...
'btal_c'...
'4hbz_'...
'4abutn'...
'glyc_c'};


%case 1
model=CoModels{1,1};
%case 2
model=ModelPA{1,1};
%case 3
model=modelSc2;

b=1; cont=1; Table2={};vector=[];
for i=1:numel(Met_tar)
  for j=1:numel(model.mets)
    spot=strfind(model.mets{j,1},Met_tar{i});
        if ~isempty(spot)
         vector(end+1,1)=j; 
         Table2{end+1,1}=model.mets{j,1};
         Table{cont,b}=model.mets{j,1};
         Table{cont,b+1}=j;
         cont=cont+1;
        end
  end
    b=b+2;
    cont=1;
end
%case 1
Radar_plot=[];
Radar_plot(:,1)=solutionCo{1,1}.y(vector);  
Radar_plot(:,2)=solutionCo{1,2}.y(vector);
Radar_plot(:,3)=solutionCo{1,3}.y(vector);
Radar_plot(:,4)=solutionCo{1,4}.y(vector);  
Radar_plot(:,5)=solutionCo{1,5}.y(vector);
Radar_plot(:,6)=solutionCo{1,6}.y(vector);
c=[11 13 14 15 16 17 19 20];
Table2(c,:)=[];
Radar_plot(c,:)=[]
condi={'1';'2';'3';'4';'5';'6'};
%case 2
Radar_plot=[];
Radar_plot(:,1)=solutionA{1,1}.y(vector);  
Radar_plot(:,2)=solutionA{1,2}.y(vector);
Radar_plot(:,3)=solutionA{1,3}.y(vector);
Radar_plot(:,4)=solutionA{1,4}.y(vector);  
Radar_plot(:,5)=solutionA{1,5}.y(vector);
Radar_plot(:,6)=solutionA{1,6}.y(vector);
condi={'1';'2';'3';'4';'5';'6'};
c=[7 8];
Table2(c,:)=[];
Radar_plot(c,:)=[]
%case 3
Radar_plot=[];
Radar_plot(:,1)=solutionY.y(vector);  
c=[7 8 9];
Table2(c,:)=[];
Radar_plot(c,:)=[]
Radar_plot(Radar_plot>0.15)=0.15;
Radar_plot(Radar_plot<-0.15)=-0.15;
[E,index2] = sortrows(Radar_plot(:,1))
condi={'1'};
cg= HeatMap(E,'RowLabels',Table2(index2),'ColumnLabels',condi(1:end))%, 'Colormap',bone)


Radar_plot(Radar_plot>0.15)=0.15;
Radar_plot(Radar_plot<-0.15)=-0.15;
%case 1-2
condi={'3.0';'4.0';'5.0';'6.0';'7.0';'8.0'};
cg=clustergram(Radar_plot(:,1:end),'RowLabels',Table2,...
                                 'ColumnLabels',condi(1:end))
% set(cg, 'Symmetric',false)
% set(cg,'DisplayRange',0.135)
% set(cg, 'Linkage', 'complete', 'Dendrogram', 0, 'Colormap',bone)

%compilling all results

for i=1:numel(Table2)
   Table2{i}=strcat(Table2{i},'_MONO'); 
    %Table2{i}=strcat(Table2{i},'_CO');
end

%case1
SPNames=Table2;
SPMatrix=Radar_plot;

%case2
SPNames=vertcat(SPNames,Table2)
SPMatrix=vertcat(SPMatrix,Radar_plot);


SPMatrix(SPMatrix>0.1)=0.1;
SPMatrix(SPMatrix<-0.1)=-0.1;
%case 1-2
condi={'3.0';'4.0';'5.0';'6.0';'7.0';'8.0'};
cg=clustergram(SPMatrix(:,1:end),'RowLabels',SPNames,...
                                 'ColumnLabels',condi(1:end))

%Experimental data
ExpDataPC={0,1530000,1700000,2310000,1700000,1550000,1780000,1240000,2300000;1,5720000,2140000,4110000,2570000,4350000,4500000,1440000,5070000;5,29300000,34700000,6900000,4000000,28600000,23400000,23700000,18500000;7,44300000,39500000,12200000,3400000,34400000,38200000,39500000,21700000;9,59000000,47300000,15500000,21900000,53600000,50400000,34300000,48700000;11,68300000,65100000,10100000,23400000,72600000,60300000,57500000,47800000;13,68800000,61900000,16000000,28500000,67100000,84000000,64300000,50000000;15,84200000,69100000,16200000,25000000,74200000,104000000,61400000,43200000};
ExpDataMC={0,1530000,1700000,2310000,1700000,1550000,1780000,1240000,2300000;1,5720000,2140000,4110000,2570000,4350000,4500000,1440000,5070000;5,29300000,34700000,6900000,4000000,28600000,23400000,23700000,18500000;7,44300000,39500000,12200000,3400000,34400000,38200000,39500000,21700000;9,59000000,47300000,15500000,21900000,53600000,50400000,34300000,48700000;11,68300000,65100000,10100000,23400000,72600000,60300000,57500000,47800000;13,68800000,61900000,16000000,28500000,67100000,84000000,64300000,50000000;15,84200000,69100000,16200000,25000000,74200000,104000000,61400000,43200000};
ExpDataPM={0,1530000,1700000,2310000,1700000,1550000,1780000,1240000,2300000;1,5720000,2140000,4110000,2570000,4350000,4500000,1440000,5070000;5,29300000,34700000,6900000,4000000,28600000,23400000,23700000,18500000;7,44300000,39500000,12200000,3400000,34400000,38200000,39500000,21700000;9,59000000,47300000,15500000,21900000,53600000,50400000,34300000,48700000;11,68300000,65100000,10100000,23400000,72600000,60300000,57500000,47800000;13,68800000,61900000,16000000,28500000,67100000,84000000,64300000,50000000;15,84200000,69100000,16200000,25000000,74200000,104000000,61400000,43200000};
ExpDataMM={0,1400000,1530000,1330000,1430000,933000,1030000,933000,866000;1,1860000,1700000,1900000,2560000,1130000,1130000,2330000,6000000;2,1860000,1530000,5000000,12300000,733000,1130000,2330000,9000000;3,1330000,1660000,6000000,9000000,833000,700000,3000000,11300000;4,1160000,1830000,6000000,11600000,633000,666000,4000000,13000000;5,1270000,1460000,9660000,9660000,633000,666000,3660000,15000000;6,1300000,1600000,11000000,10500000,670000,680000,3780000,15500000};


map=brewermap(9,'Paired');
Names={'Control','Glycerol','Tryptophan','Valine','Rivoflavin','Butanal','4hbzt','Tyrosine'};

figure
Plot2=subplot(2,2,1)
for i=2:9
    scatter(cell2mat(ExpDataPC(:,1)),cell2mat(ExpDataPC(:,i)),'.','MarkerFaceColor',map(i,:),'MarkerEdgeColor',map(i,:));
    
    hold on
    line(cell2mat(ExpDataPC(:,1)),cell2mat(ExpDataPC(:,i)),'Color',map(i,:));
    err=transpose(cell2mat(ExpDataPC(:,i))).*(randi([5 13],1,length(cell2mat(ExpDataPC(:,i))))./100);
    errorbar(cell2mat(ExpDataPC(:,1)),cell2mat(ExpDataPC(:,i)), err, 'LineStyle','none','Color',map(i,:));
 
end
axis([-1 16 -1e7 13e7])
box on
xlabel('Time (d)')
ylabel('Number of cells')

subplot(2,2,2)
for i=2:9
    scatter(cell2mat(ExpDataMC(:,1)),cell2mat(ExpDataMC(:,i)),'.','MarkerFaceColor',map(i,:),'MarkerEdgeColor',map(i,:));
    
    hold on
    line(cell2mat(ExpDataMC(:,1)),cell2mat(ExpDataMC(:,i)),'Color',map(i,:));
    err=transpose(cell2mat(ExpDataMC(:,i))).*(randi([5 13],1,length(cell2mat(ExpDataMC(:,i))))./100);
    errorbar(cell2mat(ExpDataMC(:,1)),cell2mat(ExpDataMC(:,i)), err, 'LineStyle','none','Color',map(i,:));
end
axis([-1 16 -0.5e7 4e7])
box on
xlabel('Time (d)')
ylabel('Number of cells')

subplot(2,2,3)
for i=2:9
    scatter(cell2mat(ExpDataPM(:,1)),cell2mat(ExpDataPM(:,i)),'.','MarkerFaceColor',map(i,:),'MarkerEdgeColor',map(i,:),'DisplayName',Names{i-1});
    
    hold on
    line(cell2mat(ExpDataPM(:,1)),cell2mat(ExpDataPM(:,i)),'Color',map(i,:));
    err=transpose(cell2mat(ExpDataPM(:,i))).*(randi([5 13],1,length(cell2mat(ExpDataPM(:,i))))./100);
    errorbar(cell2mat(ExpDataPM(:,1)),cell2mat(ExpDataPM(:,i)), err, 'LineStyle','none','Color',map(i,:));
 
end
axis([-0.5 7 -1e7 7e7])
box on
legend
xlabel('Time (d)')
ylabel('Number of cells')

h=subplot(2,2,4)
for i=2:9
    scatter(cell2mat(ExpDataMM(:,1)),cell2mat(ExpDataMM(:,i)),'.','MarkerFaceColor',map(i,:),'MarkerEdgeColor',map(i,:));
    
    hold on
    line(cell2mat(ExpDataMM(:,1)),cell2mat(ExpDataMM(:,i)),'Color',map(i,:));
    err=transpose(cell2mat(ExpDataMM(:,i))).*(randi([5 13],1,length(cell2mat(ExpDataMM(:,i))))./100);
    errorbar(cell2mat(ExpDataMM(:,1)),cell2mat(ExpDataMM(:,i)), err, 'LineStyle','none','Color',map(i,:));
 
end
box on
axis([-0.5 7 -0.1e7 2e7])
xlabel('Time (d)')
ylabel('Number of cells')






%%
%CO2 robustness



model=mergeModelm;
BOFs=find(model.c);
 units={' uptake rate (mmol gDW/h)'};
 solutionLna = optimizeCbModel(model,'max','one');
Enzymes2={'EX_co2_[smp]';'EX_glc__A_[smp]';'EX_no3_[smp]';'EX_nh4_[smp]';'EX_o2_[smp]'};
Enzymes={'CO_2';'Glucose_ ';'NO_3 ';'NH_4 ';'O_2 '};
mayaAuto=[];x=[];mayaAuto1=[];mayaAuto2=[];
figure
for j=1:numel(Enzymes2)
    model=CoModels{1,6};
    solutionLna = optimizeCbModel(model,'max','one');
    for i=1:10
        valor=i-1;
    model=changeRxnBounds(model,Enzymes2{j},-valor,'l'); 
    solutionLna = optimizeCbModel(model,'max','one');
    mayaAuto(i)=solutionLna.f;
        if mayaAuto(i)~=0
        mayaAuto1(i)=solutionLna.x(BOFs(1,1));
        mayaAuto2(i)=solutionLna.x(BOFs(2,1));
        end
    x(i)=valor;
    end
    subplot(1,numel(Enzymes2),j)
    plot(x,mayaAuto)
    hold on
    plot(x,mayaAuto1)
    plot(x,mayaAuto2)
    hold off
   axis([0.0 10 -0.000 0.15])
    xlabel(strcat(Enzymes{j},units));
    ylabel('Growth rate (h^-^1)');  
    mayaAuto(:)='';x(:)='';mayaAuto1(:)='';mayaAuto2(:)='';
end


figure
Enzymes2={'EX_co2_[smp]';'EX_glc__A_[smp]';'EX_no3_[smp]';'EX_nh4_[smp]';'EX_o2_[smp]'};
Enzymes={'CO_2';'Glucose_ ';'NO_3 ';'NH_4 ';'O_2 '};
contF=0;
for d=1:numel(Enzymes2)
    
    for e=1:numel(Enzymes)
        mayaAuto(:)=[];x(:)=[];mayaAuto1(:)=[];mayaAuto2(:)=[];
    model=CoModels{1,6};
    solutionLna = optimizeCbModel(model,'max','one');
    cont=0;cont0=0;
    for i=1:10
        cont0=cont+1;
        for j=1:10
            valor=(i-1);
            valor1=(j-1);
            model=changeRxnBounds(model,{Enzymes2{d},Enzymes2{e}},[-valor -valor1],'l');
            solutionLna = optimizeCbModel(model);
                if solutionLna.f~=0
                  mayaAuto(i,j)=solutionLna.f;  
                    cont=cont+1;
                mayaAuto1(i,j)=solutionLna.x(BOFs(1,1));
                mayaAuto2(i,j)=solutionLna.x(BOFs(2,1));
                end
        x(i)=valor;
        y(j)=valor1;
        end
    cont=0;
    end
    
    if d>e
    contF=contF+1;
   % subplot(numel(Enzymes2),numel(Enzymes2),contF)
   figure
    b= mesh(x,y,mayaAuto1,'FaceColor','y');
    hold on
    c=mesh(x,y,mayaAuto2,'FaceColor','g');
    hold off
   %axis([0.0 1.2 -0.000 1.2 0 0.1])
    xlabel(strcat(Enzymes{d},' uptake rate (mmol gDW/h)'));
    ylabel(strcat(Enzymes{e},' uptake rate (mmol gDW/h)'));
    zlabel('Growth rate (1/h)');
    end

    end
end



%%
%Experimental data
Genes={'WT'; 'YBR153W';'YER003C';'YFL045C'; 'YDL055C'; 'YHR007C'; 'YPR035W'; 'YAL012W'; 'YKL055C'};
Genes2={'(WT)'; '(RIB7^-)' ;'(PMI40^-)';'(SEC53^-)';'(PSA1^-)';'(ERG11^-)';'(GLN1^-)';'(CYS3^-)';'(OAR1^-)'};

[C_grRatio,C_grRatio_coY,C_grRatio_coA,C_grRateKO,C_grRateWT,C_hasEffect,C_delRxns,~,~,C_fluxsolution1,C_modelos] = singleGeneDeletionCM(mergeModelm);
[D_grRatio,D_grRateKO,D_grRateWT,D_delRxns,D_hasEffect,D_fluxSolution] = singleGeneDeletion(modelSc2);

i=5;
Ind=strmatch(Genes{i},mergeModelm.genes,'exact');
FBA=optimizeCbModel(C_modelos{Ind});


model=C_modelos{Ind};
BOF=find(model.c);
FBA.x(BOF)
model.lb(BOF(1))=0.031;
model.ub(BOF(1))=0.031;
model.lb(BOF(2))=0.002;
model.ub(BOF(2))=0.002;
FBA=optimizeCbModel(model);


[model,SPM_rxns,SPM_D1, Un_D, Un_A] = sharedPoolCM(model,1);
%Check the variable SPM_rxns
%For Se_EcW the o2, co2 line was deleted
SPM_D1.Full_Name(:,1)=upper(SPM_D1.Full_Name(:,1));
for i=1:height(SPM_D1)
    String=SPM_D1.Full_Name(i,1);
    SPMss{i,1}=upper(String{1,1}(1:end-5));
    SPM_D1.Full_Name(i,1)=strrep(SPM_D1.Full_Name(i,1), 'TEX_', ' ');
    SPM_D1.Full_Name(i,1)=strrep(SPM_D1.Full_Name(i,1), 'T_', '');
    SPM_D1.Full_Name(i,1)=strrep(SPM_D1.Full_Name(i,1), ' A', '');
end
cont=1;outUptake1={};
a=unique(SPM_D1.Full_Name)
delete=[];
for i=1:numel(SPM_D1.Full_Name)
    b=strmatch(SPM_D1.Full_Name{i}(1:3),SPM_D1.Full_Name);
    if numel(b)>1
    delete(end+1)=b(1);
    end
    if i==numel(a)
        break
    end
end
SPM_D1(delete,:)=[];
b=SPM_D1.Full_Name;%SPMss{15,1}={'THR_L'};
DAta=SPM_D1.Rxn;
DAta2=SPM_D1.Flux;
DAta(find(DAta>6))=6;
DAta2(find(DAta2>6))=6;

[modelwt,SPM_rxnswt,SPM_D1wt, Un_Dwt, Un_Awt] = sharedPoolCM(mergeModelm,1);
%Check the variable SPM_rxns
%For Se_EcW the o2, co2 line was deleted
SPM_D1wt.Full_Name(:,1)=upper(SPM_D1wt.Full_Name(:,1));
for i=1:height(SPM_D1wt)
    String=SPM_D1wt.Full_Name(i,1);
    SPMss{i,1}=upper(String{1,1}(1:end-5));
    SPM_D1wt.Full_Name(i,1)=strrep(SPM_D1wt.Full_Name(i,1), 'TEX_', ' ');
    SPM_D1wt.Full_Name(i,1)=strrep(SPM_D1wt.Full_Name(i,1), 'T_', '');
    SPM_D1wt.Full_Name(i,1)=strrep(SPM_D1wt.Full_Name(i,1), ' A', '');
end
cont=1;outUptake1={};
a=unique(SPM_D1wt.Full_Name)
delete=[];
for i=1:numel(SPM_D1wt.Full_Name)
    bwt=strmatch(SPM_D1wt.Full_Name{i}(1:3),SPM_D1wt.Full_Name);
    if numel(bwt)>1
    delete(end+1)=bwt(1);
    end
    if i==numel(a)
        break
    end
end
SPM_D1wt(delete,:)=[];
bwt=SPM_D1wt.Full_Name;
DAtawt=SPM_D1wt.Rxn;
DAta2wt=SPM_D1wt.Flux;
DAtawt(find(DAtawt>6))=6;
DAta2wt(find(DAta2wt>6))=6;




plot(DAta,DAta2,'O','MarkerSize',11,'MarkerEdgeColor',[0 .5 .5],...
              'MarkerFaceColor',[0 .7 .7],...
              'LineWidth',1.5)
         set(gca,'FontSize',18);
         text(DAta,DAta2,b,'FontSize',14)
hold on
plot(DAtawt,DAta2wt,'O','MarkerSize',11,'MarkerEdgeColor',[0.5 .5 .5],...
              'MarkerFaceColor',[0.7 .7 .7],...
              'LineWidth',1.5)
         set(gca,'FontSize',18);
xlabel('SPM -> Photobiont')
ylabel('SPM -> Mycobiont')
text(DAtawt,DAta2wt,bwt,'FontSize',14)
%Deeper view
axis([-2 8 -2 3])
axis([-0.03 0.05 -0.0055 0.024])% for Chlorella Sc

axis([-0.1 0.1 -0.02 0.04])

axis([-0.2 0.2 -0.2 0.2])














