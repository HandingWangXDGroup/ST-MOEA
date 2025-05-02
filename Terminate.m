function TFlag = Terminate()
    global Global;
    if Global.CurTime<=Global.Times
        TFlag    = 1;
    else
        SavesData= Global.Saves;
        % HyperNum = Global.HyperNum;
        TFlag    = 0; 
    end
end