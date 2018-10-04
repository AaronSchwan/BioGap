import scipy
from scipy import stats

def kaplan(kaplandata1,kaplancen1,kaplandata2,kaplancen2):
    
    kaplandata1, kaplancen1 = (list(t) for t in zip(*sorted(zip(kaplandata1, kaplancen1))))
    kaplandata2, kaplancen2 = (list(t) for t in zip(*sorted(zip(kaplandata2, kaplancen2))))

    tot = len(kaplandata1)+len(kaplandata2)
    tot1 = len(kaplandata1)
    tot2 = len(kaplandata2)
    risk1 = [];
    risk2 = [];
    breaker = 0;
    censord1 = 0;
    censord2 = 0;

    day = 0;
    if(max(kaplandata1) > max(kaplandata2)):
        maxday = max(kaplandata1)

    if(max(kaplandata1) < max(kaplandata2)):
        maxday = max(kaplandata2)

    if(max(kaplandata1) == max(kaplandata2)):
        maxday = max(kaplandata2)

        

    while(day<=maxday):
        
        deaths = 0;

        matches1 = [i for i in kaplandata1 if i == day]
        matches2 = [i for i in kaplandata2 if i == day]

        totmatches = matches1+matches2;

        if(len(matches1) >= 1):
            kap1 =[];
            cen1 =[];

            ss = kaplandata1.index(day)
            s = len(matches1)
            sss = 0;
            
            
            while(s>0):

                kap1.insert(sss,kaplandata1[ss])
                
                cen1.insert(sss,kaplancen1[ss])

                if(cen1[sss] == 0):
                    deaths+=1;
                
                s-=1;
                ss+=1;
                sss+=1;

            censord1+=sum(cen1)
            
        if(len(matches2) >= 1):
            kap2 = [];
            cen2 = [];
            
            ss = kaplandata2.index(day)
            s = len(matches2)
            sss = 0;
            
            
            while(s>0):

                kap2.insert(sss,kaplandata2[ss])
                
                cen2.insert(sss,kaplancen2[ss])


                if(cen2[sss] == 0):
                    deaths+=1;
                
                s-=1;
                ss+=1;
                sss+=1;

            censord2+=sum(cen2)

        if(len(totmatches) != 0):
            
            tot = tot1+tot2;
            
            risk = deaths/tot;
            
          
            
            risk1.insert(breaker,tot1*risk)
            risk2.insert(breaker,tot2*risk)
            
            

            
            tot1-=len(matches1);
            tot2-=len(matches2);

            breaker+=1;

        day+=1;
   

    tot1expected = sum(risk1)
    tot2expected = sum(risk2)

   

    

    observed1 = len(kaplandata1)-sum(kaplancen1)
    observed2 = len(kaplandata2)-sum(kaplancen2)

    #print(observed1,observed2)

    

    ch2, pval = scipy.stats.chisquare(f_obs=[observed1,observed2], f_exp=[tot1expected,tot2expected])

    #print(pval)

    return pval;
    
