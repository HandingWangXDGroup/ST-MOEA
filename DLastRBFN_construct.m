function [AllRnets,LastFlag] = DLastRBFN_construct(IniData,LastFlag)
    global Global
    AllRnets = struct();
    IniDec   = decs(IniData);
    N        = size(IniDec,1);     
    center_num = ceil(sqrt(N)); 
    centers= get_center(IniDec, center_num);
    
    IniObj = objs(IniData); 
    
    
    Sepvec = linspace(0,1,N+1);
    MinObj = min(IniObj); MaxObj = max(IniObj);
    IniObj = (IniObj-MinObj)./repmat((MaxObj-MinObj),N,1);
    HyperSigma  =[];
    for j = 1:Global.M
        OneObj = IniObj(:,j);
%         RandID = randi([1,center_num]);
%         OriOneObj = OneObj;
%         OneObj = (OneObj-min(OneObj))./(max(OneObj)-min(OneObj));
        ArrNum = histcounts(OneObj,Sepvec);
%         ArrID  = discretize(OneObj,Sepvec);
        DifNum = abs(ArrNum-1);
        MeanDif= sum(DifNum)/N;



        if MeanDif <= 1
            if LastFlag(j) >= 0 
                HyperSigma = mean(pdist(centers))*2;
                % Global.alphanum(j,1) = Global.alphanum(j,1)+1;
            else
                HyperSigma = mean(pdist(centers));
                % Global.alphanum(j,2) = Global.alphanum(j,2)+1;
            end
            LastFlag(j) = 1;
        elseif MeanDif > 1
            if LastFlag(j) <=0
                HyperSigma = mean(pdist(centers))/2;
                % Global.alphanum(j,3) = Global.alphanum(j,3)+1;
            else
                HyperSigma = mean(pdist(centers));
                % Global.alphanum(j,1) = Global.alphanum(j,1)+1;
            end
            LastFlag(j) = -1;
        end

        AllRnets(j).net = construct_RBFNHyper(decs(IniData),objs(IniData), centers, HyperSigma);
    end
end