
def kaplan(kaplandata,kaplancen):

    xx= len(kaplandata)
    tot  = len(kaplandata)
    division = len(kaplandata)
    kapper = [];   
    kaplandata, kaplancen = (list(t) for t in zip(*sorted(zip(kaplandata, kaplancen))))
    idsdays = [];
    censoreds = [];
    censoredids = [];
    place = 0;
    place2 = 0;
    place3 = 0;
    idsdays = [];
    former = 1;

    while(xx>0):

        cen = [];
        kap = [];
        matches = [];
        deaths = 0;

        matches = [place for i in kaplandata if i == kaplandata[place]]

        s = len(matches)
        ss = matches[0]
        sss = 0;
        
        
        while(s>0):
            
            kap.insert(sss,kaplandata[ss])
            
            cen.insert(sss,kaplancen[ss])

            if(cen[sss] == 0):
                deaths+=1;
            
            s-=1;
            ss+=1;
            sss+=1;

        censord = sum(cen)
        ov = (deaths/tot)
        st = 1 - ov;    
        percent = former*st

        if(censord >= 1):
            censoreds.insert(place3,percent)
            censoredids.insert(place3,kap[0])
            place3+=1;

        kapper.insert(place2,former)
        idsdays.insert(place2,kap[0])
        place2+=1;
        kapper.insert(place2,percent)
        idsdays.insert(place2,kap[0])
        place2+=1;

        former =percent;
        tot-=len(matches);
        
        place+=len(matches);
        xx-=len(matches);
    [1]+kapper
    [0]+idsdays
    x = len(kapper)

    kapper.insert(x,0)
    idsdays.insert(x,max(idsdays))

    
    return kapper,idsdays,censoreds,censoredids;
    

    



