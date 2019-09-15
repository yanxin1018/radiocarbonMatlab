function s = convertAge(CalAge)
if CalAge<0
    s=[num2str(abs(CalAge)),' CalBC'];
    if CalAge<-50000
        s='<50000 CalBC';
    end
else if CalAge==0
        s=' CalAD';
    else
        s=[num2str(CalAge),' CalAD'];
        if CalAge>1950
            s='>1950 CalAD';
        end
    end
end
end