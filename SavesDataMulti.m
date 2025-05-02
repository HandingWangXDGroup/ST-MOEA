function SavesDataMulti(Population, KModel, Rnets, KPopulation, RPopulation)
    global Global
    Global.Saves(Global.CurTime).RPop = RPopulation;
    Global.Saves(Global.CurTime).KPop = KPopulation;
    Global.Saves(Global.CurTime).OutPop= Population;
    Global.IGD(Global.CurTime) = CalIGD(Population);
    
end