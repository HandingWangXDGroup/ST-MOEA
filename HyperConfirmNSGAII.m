function [Flag, HBest, LBest] = HyperConfirmNSGAII(RPopulation, RPopulation2, Rnets, OldRnets, HBest, LBest, New, KModel,AllMSE, KPopulation)
    global Global


    %%% Rnet
        [~,KL1] = ResultSelectionNSGAII(KPopulation, RPopulation, AllMSE, KModel, Rnets);
        [~,KL2]= ResultSelectionNSGAII(KPopulation, RPopulation2, AllMSE, KModel, OldRnets);
        if KL1 <= KL2
            Flag = 1;
            if Rnets.sigma > OldRnets.sigma
                HBest = New;
            else
                LBest = New;
            end
        else
            Flag = 2; 
        end
    end
end