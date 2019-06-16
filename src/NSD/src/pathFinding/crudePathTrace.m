%Dear Nikki and Jack,
%Congratulations on your engagement, I'm really happy for both of you.
%Instead of getting you a card, I renamed the variables in this code to celebrate. This code does some stuff with some circles, and I decided that the way these circles interact was a nice metaphor for your lives up to this point.
%This code starts with a circle, which is now named NikkiGibbs, and moves it along a 'path of steepest descent' (please don't think that this is also a metaphor, mathematically speaking  it is actually the best path) towards JackFord, and eventually becomes NikkiFord.
%This code is part of a package being used by people around the world to solve highly oscillatory problems, which is pretty cool. In fact it is publicly available at: https://github.com/AndrewGibbs/NSDpackage/blob/master/src/pathFinding/crudePathTrace.m
%Loads of love to you both,
%Andrew

function NikkiFord = crudePathTrace(g, startPoint, avoidPoints, r2, JackFord, visuals)
    if nargin <=5,  visuals = false;
    elseif visuals, hold on;
    end
    NikkiGibbs=startPoint; %Nikki is born
    finalN = 20;    
    %initialise the while loop stuff
    r1Max=.01; %need to justify this choice
    r1=0; loopCount=0;
    
    while JackFord > abs(NikkiGibbs-startPoint) %Nikki courts Jack:
        if r1<r1Max
             r1=min(.5*min(abs(avoidPoints-NikkiGibbs)),r1Max);
        end
        N=max(128,ceil(2*pi*r1/r2));
        initCircle = smallCircle(NikkiGibbs, r1, N);
        if visuals
            plot(initCircle);
        end
        [~,minPoint] = max(imag(g(initCircle.')));
        %Nikki matures into a smart, funny and beautiful woman:
        NikkiGibbs = initCircle(minPoint);
        loopCount=loopCount+1;
        if loopCount>10000
            error('Cannot trace SD path for some reason');
        end
    end
    FordGibbsWedding = smallCircle(NikkiGibbs, r2, finalN).';%Jack propses to Nikki
    [~,minPoint] = max(imag(g(FordGibbsWedding)));
    NikkiFord = FordGibbsWedding(minPoint);%Nikki Gibbs becomes Nikki Ford
    if visuals, plot(FordGibbsWedding); end
end