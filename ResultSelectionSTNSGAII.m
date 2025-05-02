function [FinalPop,KL] = ResultSelectionSTNSGAII(KPopulation, RPopulation, AllMSE, KModel,AllRnets)
    global Global
    %
    KDec = decs(KPopulation);     KObj = objs(KPopulation); AllMSE(find(AllMSE<10^(-3))) = 0;
    MKObj= KObj+3*sqrt(AllMSE);   MKPopulation = PopStruct(KDec,MKObj);

    %
    RDec = decs(RPopulation);
    MSE  = zeros(size(RDec,1),Global.M);
    RObj = zeros(size(RDec,1),Global.M);
    for i = 1: size(RDec,1)
        for j = 1 : Global.M
            [RObj(i,j),~,MSE(i,j)] = predictor(RDec(i,:),KModel{j});
        end
    end   
    MSE(find(MSE<10^(-3))) = 0;
    RObj = RObj-3*sqrt(MSE); MRPopulation = PopStruct(RDec,RObj);

    %
    KRNext = EnvironmentalSelectionAdapt([MKPopulation,MRPopulation],Global.N);
    RSelectSize = Global.N-sum(KRNext(1:Global.N));

    %
    LeftKID = find(KRNext(1:Global.N)~=1);
    LeftPop = [KPopulation(LeftKID),RPopulation];
    LeftDec = decs(LeftPop);
    LeftObj = zeros(size(LeftDec,1),Global.M);
    for j = 1:Global.M
        PartObj = RBFN_cal(LeftDec, AllRnets(j).net);
        LeftObj(:,j) = PartObj(:,j);
    end
    LeftPop = [PopStruct(LeftDec,LeftObj),RPopulation];
    RNext= EnvironmentalSelectionAdapt(LeftPop, RSelectSize);
    KSID = find(KRNext(1:Global.N)==1);
    
        
    FinalPop = [KPopulation(KSID),LeftPop(RNext)];
%    disp(CalIGD(KPopulation));
%     disp(CalIGD(RPopulation));
    KL = length(KSID);
%     disp(KL);
end