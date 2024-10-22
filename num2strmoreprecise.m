%General supporting function. Makes sure that the precision provided in a string, converted from a number, is high enough.

function startstring=num2strmoreprecise(number,precision)
    if number == 0
        startstring=['0'];
    else
    if number < 0
        negative = 1;
        number = abs(number);
    else
        negative = 0;
    end
    a=floor(log(number)/log(10));
    addednr=floor(number./10.^floor(log(number)/log(10)));
    startstring=[num2str(addednr) '.'];
    number=number-addednr.*10.^floor(log(number)/log(10));
    i=1;
    while i < precision
        if floor(log(number)/log(10))< a-i
            startstring = [startstring '0'];
        elseif i==precision-1;
            addednr=round(number./10.^floor(log(number)/log(10)));
            startstring = [startstring num2str(addednr)];
        else
            addednr=floor(number./10.^floor(log(number)/log(10)));
            number=number-addednr.*10.^floor(log(number)/log(10));
            startstring = [startstring num2str(addednr)];
            
        end
        i=i+1;
    end
    startstring = [startstring 'e' num2str(a)];
    if negative
        startstring=['-' startstring];
    end
    end