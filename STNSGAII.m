function STNSGAII
    global Global;
    cd(fileparts(mfilename('fullpath')));
    addpath(genpath(cd));
    %% Parameter setting
    Global.CurTime= 0; 
    Global.D = 10;
    Global.N = 100;
    Global.nt = 10;
    ParameterInitial('F1')    
    Global.CurGen = 0;  
    Global.Saves  = struct();
    Global.IGD    = zeros(1, Global.Times);    
    KModel        = cell(1,Global.M);
    THETA         = 5.*ones(Global.M,Global.D);
    LastIniTHETA  = 5.*ones(Global.M,Global.D);
    LastFlag      = zeros(1,Global.M);
    THEChange = 0;
    %% Optimization
    while Terminate()
        if Global.CurGen == Global.Gen || Global.CurGen ==0 
            if Global.CurGen == Global.Gen
                disp(['**************Time step: ', num2str(Global.CurTime), '**************'])
                %%% The combination between RPop and KPop
                [FinalPop,~] = ResultSelectionSTNSGAII(KPopulation, RPopulation, AllMSE, KModel, Rnets);
                SavesDataMulti(FinalPop, KModel, Rnets, KPopulation, RPopulation);
            end       
            Global.CurGen = 1;
            Global.CurTime = Global.CurTime + 1;
            IniData = LHS_sam(Global.D*11-1);   
            
            %%%% Record Offline data IGD
            [EOff,~,~] = EnvironmentalSelection(IniData,Global.N);
            Global.IniIGD(Global.CurTime) = CalIGD(EOff);
            
            %%%Select the parameter of RBFN and Construct RBFN
            [Rnets,LastFlag] = DLastRBFN_construct(IniData,LastFlag);
            OldTHETA = THETA;
            %%% Kriging Construction
            IniObj = objs(IniData); 
            for i = 1 : Global.M
                dmodel     = dacefit(decs(IniData),IniObj(:,i),'regpoly1','corrgauss',THETA(i,:),1e-5.*ones(1,Global.D),10^5.*ones(1,Global.D));
                KModel{i}  = dmodel;
                THETA(i,:) = dmodel.theta;  
            end
            THEDis = OldTHETA-THETA;
            if sum(sum(THEDis))==0
                THEChange = THEChange+1;
                if THEChange >= 5
                    if LastIniTHETA(1) == 5
                        THETA = ones(Global.M,Global.D);
                    else
                        THETA = 5.*ones(Global.M,Global.D);
                    end
                    LastIniTHETA = THETA;
                    THEChange = 0;
                end               
            end
            RPopulation = IniData;
            KPopulation = IniData;
            [~,RFrontNo,RCrowdDis] = EnvironmentalSelection(RPopulation,length(IniData)); 
            [~,KFrontNo,KCrowdDis] = EnvironmentalSelection(KPopulation,length(IniData));
            AllMSE = zeros(length(IniData),Global.M);
        else
            %%% Rnets Optimization
            RMatingPool = TournamentSelection(2,length(RPopulation),RFrontNo,-RCrowdDis);
            ROffDec     = GA(decs(RPopulation(RMatingPool)));
            ROffObj     = zeros(size(ROffDec,1),Global.M);
            for j = 1:Global.M
                PartObj = RBFN_cal(ROffDec, Rnets(j).net);
                ROffObj(:,j) = PartObj(:,j);
            end
            ROffspring  = PopStruct(ROffDec, ROffObj);
            [RPopulation,RFrontNo,RCrowdDis] = EnvironmentalSelection([RPopulation,ROffspring],Global.N);

            %%% Kriging Optimization
            KMatingPool = TournamentSelection(2,length(KPopulation),KFrontNo,-KCrowdDis);
            KOffDec     = GA(decs(KPopulation(KMatingPool))); 
            KOffObj     = zeros(size(KOffDec,1),Global.M);  
            MSE         = zeros(size(KOffDec,1),Global.M);
            for i = 1: size(KOffDec,1)
                for j = 1 : Global.M
                    [KOffObj(i,j),~,MSE(i,j)] = predictor(KOffDec(i,:),KModel{j});
                end
            end
            AllMSE = [AllMSE;MSE];
            KOffspring  = PopStruct(KOffDec, KOffObj);
            [KPopulation,KFrontNo,KCrowdDis] = EnvironmentalSelection([KPopulation,KOffspring],Global.N);
            Next = EnvironmentalSelectionAdapt([KPopulation,KOffspring],Global.N);
            AllMSE = AllMSE(find(Next==1),:);
            
            Global.CurGen = Global.CurGen + 1;
        end
    end
    MIGD     = mean(Global.IGD);
    IGD      = Global.IGD;
end