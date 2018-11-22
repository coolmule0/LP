function LPtimingScript()
%Script to read in timing files of LP methods and construct graphs

    [timing, timingstd, size, batch, Title]= ReadInValues();
    
    %PlotSurf(timing, size, batch, Title);
    
    
    %figs
    close all;
    
    set = [1,2,4,5, 9];
    PlotCutthrough(timing, timingstd, size, batch, Title, true, 4, set, 10^2, 2*10^5, 1, inf)
    PlotCutthrough(timing, timingstd, size, batch, Title, true, 8, set, 10^2, 2*10^5, 1, inf)
    PlotCutthrough(timing, timingstd, size, batch, Title, true, 11, set, 10^2, 2*10^5, 1, inf)
    PlotCutthrough(timing, timingstd, size, batch, Title, false, 3, set, 10^2, 2*10^5, 1, inf)
    PlotCutthrough(timing, timingstd, size, batch, Title, false, 10, [1,2], 10^2, 2*10^5, 1, inf)

%     PlotCutthrough(timing, timingstd, size, batch, Title, false, 10, set, 10^2, 2*10^5, 1, inf)
%     PlotSurf(timing{7}./(timing{7}+timing{8}),size{7},batch{7}, "Proportional amount of time spent transfering memory")

    %PlotCutthrough(timing, timingstd, size, batch, Title, false, 10, set, 10^2, 2*10^5, 5, 2*10^5)
    %PlotCutthrough(timing, timingstd, size, batch, Title, true, 14, [1,2], 10, 10^4, -inf, 10^5)
    %PlotCutthrough(timing, timingstd, size, batch, Title, false, 5, [7,8], 500, inf, 0.8, inf)
    %PlotCutthrough(timing, timingstd, size, batch, Title, false, 12, [1,2,3,4], 10, inf, 10, inf)
    %CompareNewOld(timing, timingstd, size, batch, Title, true, 9, 3, inf, -inf, inf)
    %CompareNewOld(timing, timingstd, size, batch, Title, true, 14, 3, 4500, -inf, inf)

    %save figures
    for y=1:4
        saveas(figure(y),strcat("figures/",num2str(y),".pdf")) %make sure file doesnt exist before saving
    end
    
    %Get max difference between 2 different sets
    comp1 = 3;
    comp2 = 4;
    [max, element] = getMax(timing, comp1, comp2);
    elemx=mod(element,length(size{comp1}));
    elemy=1+floor(element/length(size{comp1}));
    fprintf('max between %s and %s is %f, at element (size: %i,batch: %i) i.e. size %i batch %i\n', Title(comp1), Title(comp2),  max, element(1), elemx, elemy, size{comp1}(elemx), batch{comp1}(elemy) );


end

%Gets the maximum speedup of the first when compared to the second method. Only works if
%they contain the same dimensions (and please ensure the elements
%correspont to the same size/batch amount)
function [Max, element] = getMax(timing, comparison1, comparison2)
    relMat = (timing{comparison2}-timing{comparison1})./timing{comparison1};
    relMat = relMat(isfinite(relMat));
    [Max,element] = max(relMat); 
    
end

function [timing, timingstd, size, batch, Title] = ReadInValues()
%Reads in the data from files to produce a matrix of LP size-batch-timing
%returns:   cells of matrix timing, the time taken to solve the LPs.
%           cells of matrix timingstd, the standard deviation of timings
%           cell of array size, the values of sizes used
%           cell of array batch, the amount of LPs solved
%           array Title, the legend names of the different results

    %Editable variables
    %Gurung and Ray timing file
    GRTiming="..\timings\GRtimings.txt";
    %RGB timing
    RGBTiming="..\timings\RGBtimings.txt";
    %CPLEX timing with my code
    CPLEXTimingMe="..\timings\CPLEXtimings.txt";
    %CPLEX timing with mps file
    CPLEXTimingMPS="..\timings\CPLEXtimingsMPS.txt";
    %mGLPK timing file
    mGLPKTiming="..\timings\GLPK_OMPtimings.txt";
    %GLPK timing file
    GLPKTiming="..\timings\GLPKtimings.txt";
    %RGB memory transfer time
    RGBMem = "..\timings\RGBMemtimings.txt";
    %RGB solving time (without memory transfer)
    RGBComp = "..\timings\RGBCalctimings.txt";
    %CLP timing
    CLP = "..\timings\CLP.txt";

    Title=[     "RGB",      "Gurung and Ray",   "GLPK",     "mGLPK",        "CLPEX",      "CPLEXmps",       "RGBMem",   "RGBComp",  "CLP"      ];
    TimeLoc=[   RGBTiming,  GRTiming,         GLPKTiming, mGLPKTiming,    CPLEXTimingMe,CPLEXTimingMPS,   RGBMem,     RGBComp,    CLP  ];
   
    for g = 1:length(TimeLoc)
        if ~ isfile(TimeLoc(g))
            fprintf("Cannot find file %s\n", TimeLoc(g));
            continue;
        end
        fid=fopen(TimeLoc(g),'r');
        sizeA = [3 Inf];
        A=fscanf(fid, '%i %i %f', sizeA);
        fclose(fid);

        Size=A(1, :);
        Batch=A(2, :);
        Time=A(3, :);

        %Combing repeated size-batch pairs into an average
        size{g}=unique(Size);
        batch{g}=unique(Batch);
        repeat=zeros(length(size{g}),length(batch{g}));
        for i = 1:length(Time)
            sizeIndex = find(size{g}==Size(i),1);
            batchIndex = find(batch{g}>=Batch(i),1);
            repeat(sizeIndex,batchIndex) = repeat(sizeIndex,batchIndex)+1;
            timingInt(sizeIndex, batchIndex,repeat(sizeIndex,batchIndex)) = Time(i);
        end

        timing{g} = mean(timingInt,3);
        timingstd{g} = std(timingInt,0,3);

        %clear intermediate vars
        clearvars A timingInt
    end

end

function PlotSurf(timing, size, batch, Title)
%Plot the various surface plots of the data

    %For each different measured data
    %log color
    CO=(timing);


    figure
    surf(batch,size,timing,CO);
    xlabel('batch');
    ylabel('size');
    zlabel('timing');
    title(Title);
    set(gca,'XScale', 'log');
    set(gca,'YScale', 'log');
    %set(gca,'ZScale', 'log');
    xlim([batch(1) batch(end)]);
    ylim([size(1) size(end)]);


% %Where CPU outperforms RGB
% ints = floor(timing{1});
% differences = ints -timing{3};
% C = (differences<=0);
% C2 = double(C);
% figure
% surf(batch{1},size{1},differences,C2);
% xlabel('batch');
% ylabel('size');
% zlabel('timing');
% title('Differences');
% set(gca,'XScale', 'log');
% set(gca,'YScale', 'log');

end

function CompareNewOld(timing, timingstd, size, batch, Title, axisSize, value, xmin, xmax, ymin, ymax)
%Compares the timing of the old and new RGB method and displays results in
%graphical form

    %Get memory transfer timing
    
    relperf = timing{6}./timing{5}; %difference proportionalised to timing
    relerr = relperf .* sqrt( (timingstd{6}./timing{6} ).^2 + (timingstd{5}./timing{5}).^2 );

%     figure
%     surf(batch{1},size{1},relperf);
%     xlabel('batch');
%     ylabel('size');
%     zlabel('Time (ms)');
%     set(gca,'XScale', 'log');
%     set(gca,'YScale', 'log');
    
    h=figure;
    hold all
    for i=6:6 %only plot certain values
        if axisSize
            errorbar(size{i},relperf(:,value),relerr(:,value)); 
        else
            errorbar(batch{i},relperf(value,:),relerr(value,:));
        end
    end
    if axisSize
        title(strcat('Fixed batch: ', num2str(batch{1}(value))));
        xlabel('size');
    else
        title(strcat('Fixed size: ', num2str(size{1}(value))));
        xlabel('batch');
    end
    ylabel('Relative time taken');
    xlim([xmin xmax]);
    ylim([ymin ymax]);
    set(gca,'XScale', 'log');
    
    set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(h,'filename','-dpdf','-r0')
    
    hold off
end

function PlotCutthrough(timing, timingstd, size, batch, Title, axisSize, value, endi, xmin, xmax, ymin, ymax)
%plot layered cutthrough of data
%PlotCutthrough(axisSize,value):
%axisSize={true or false}, whether to vary LP size (true) or batch amount (false)
%value={int} what element to plot 
%endi=which elements of the cell array to use
%xmin, xmax=x-axis limits

    h = figure;
    hold all
    for i=endi %only plot certain values
        if axisSize
            errorbar(size{i},timing{i}(:,value),timingstd{i}(:,value));
        else
            errorbar(batch{i},timing{i}(value,:),timingstd{i}(value,:));
        end
    end
    legendTextArr = Title(endi);
    legend(legendTextArr, 'location', 'northwest');
    if axisSize
        title(strcat('Fixed batch: ', num2str(batch{1}(value))));
        xlabel('size');
    else
        title(strcat('Fixed size: ', num2str(size{1}(value))));
        xlabel('batch');
    end
    ylabel('Time (ms)');
    xlim([xmin xmax]);
    ylim([ymin ymax]);
    set(gca,'XScale', 'log');
    set(gca,'YScale', 'log');
    
    set(h,'Units','Inches');
    pos = get(h,'Position');
    set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    print(h,'filename','-dpdf','-r0')

    hold off
    
end