% **************************************************************************
%                       D R I V E R  f o r						
%        U P E M D  A N A L Y S I S  O F  B A R O R E F L E X 			
% **************************************************************************
%												
% License:										
% Copyright (C) 2019 J. R. Geddes, J. Mehlsen, and M. S. Olufsen
%
% Code, in part, adapted from code supplied by Dr. Ben Randall, University
% of Michigan 
%
% Contact information:									
% Mette Olufsen (msolufse@ncsu.edu)
% North Carolina State University
% Raleigh, NC
% 
% Permission is hereby granted, free of charge, to any person obtaining a
% copy of this software and associated documentation files (the "Software"),
% to deal in the Software without restriction, including without limitation
% the rights to use, copy, modify, and merge the Software subject to the 
% following conditions:
% 
% The above copyright notice and this permission notice shall be included
% in all copies or substantial portions of the Software.
% 
% The authors shall be cited in any work resulting from the use of the 
% Software. The associated published article is
% https://doi.org/10.1109/TBME.2020.2974095. 
% 
% THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
% WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
% MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
% ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES 
% WHATSOEVER RESULTING FROM LOSS OF USE, OR DATA, WHETHER IN AN ACTION OF 
% CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF OR IN 
% CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE
%
%%********************************************************************************

%%********************************************************************************
% Data_pre_proc(data,pkprom,figureson)
%
% Preprocesses given data to extract Heart Rate and Systolic Blood Pressure
% 
%
% Input: 
% .txt file containing time, ECG, and blood pressure data sampled at 1000 Hz
% pkprom 
% pkprom (scalar) defining the minimum peak prominence, needed for findpeaks (Matlab function) default value 25 
% figureson (0 = no figures, 1= figures) 

%
% Dependencies: MATLAB signal processing toolbox

% Output: .txt file containing time, Heart Rate, and systolic blood pressure
% data sampled at 250 Hz 
%*********************************************************************************

%%********************************************************************************
% Acknowledgement
% Code is inspired by similar code that was supplied courtesy of Dr. Eric
% Benjamin Randall 
%%********************************************************************************


function [Prepped_data] = Data_pre_proc(data,pkprom,figureson)
%This function turns a matrix with time, ECG, and BP to a matrix with
%Time, ECG, HR, and Systolic Blood Pressure
%Code adapted from code supplied by Dr. Ben Randall

sub_figs_on = 0;


t = data(:,1)-data(1,1); %Make initial time 0
ECG = data(:,2)*1000; %Needs to be in mV
BP = data(:,3);

HR = ECG_to_HR(t,ECG,sub_figs_on);
SPdata = SBPcalc(t,BP,pkprom,sub_figs_on);


dat = [t,HR,SPdata];

%Sample the data: CHANGE FOR DIFFERENT Hz CONVERSION
sample = 1:4: (floor(length(data)/4)*4); %Go from 1/250 to 1/1000

Prepped_data = dat(sample,:);

%% Plotting data 

if figureson == 1
    figure
    subplot(3,1,1)
    plot(t,ECG)
    ylabel('ECG mV')
    title('Data to be analyzed')
    subplot(3,1,2)
    plot(t,HR)
    ylabel('HR (bpm)')
    subplot(3,1,3)
    hold on 
    plot(t,BP)
    plot(t,SPdata,'r','linewidth',2)
    hold off
    legend('BP','SBP','location','southeast')
    ylabel('Pressure (mmHg)')
    xlabel('Time (s)')
end


%% ECG to HR function

function [E_HR] = ECG_to_HR(Tdata,ECG,figureson)
    %{
    Makes continuous respiration signal. Loads in 
        Tdata - Time points in seconds (sec)
        ECG   - ECG data must be in millivolts (mV)
        figureson - If 1, it plots a bunch of figures 
    %}
%    otherfigureson = 0; %Turn on and off all graphs except hr vs time
    dt = mean(diff(Tdata)); 

    %Correct baseline of ECG signal with medfilt1

    %Filter out P waves and QRS complex with a window of 200 ms
    smoothECG = medfilt1(ECG,round(.2/dt)); 

    %Filter out T waves with a window of 600 ms 
    smoothECG2 = medfilt1(smoothECG,round(.6/dt)); 

    %Baseline corrected ECG signal
    BaselineCorrectedECG = ECG - smoothECG2; 

    % Savitsky-Golay Filter 

    %Savitsky-Golay Filter filters out VERY low frequency signals. The window
    %must be odd and within .15 sec 
    SVGwindow = round(.15/dt); 
    if mod(SVGwindow,2) == 0
        SVGwindow = SVGwindow + 1;
    end 
    %Check to ensure order is less than window size 
    if SVGwindow > 5
        SVGorder = 5; 
    else
        SVGorder = 3; 
    end 
    smoothECG3 = sgolayfilt(BaselineCorrectedECG,SVGorder,SVGwindow); 

    % Accentuate peaks to easily find them 

    %Can now extract Q and R peaks 
    accentuatedpeaks = BaselineCorrectedECG - smoothECG3; 

    %Finds Q and R points with minimum peak distance of 200 ms 
    [~,z] = findpeaks(accentuatedpeaks,'MinPeakDistance',round(.2/dt)); 
    zz = mean(accentuatedpeaks(z)); 
    [~,iR] = findpeaks(accentuatedpeaks,'MinPeakHeight',zz,'MinPeakDistance',round(.2/dt)); 

    RRint = diff(Tdata(iR));%make this a step function with nodes at T_RRint  
    T_RRint = Tdata(iR(1:end-1)); 

    HRi = 60./RRint; 
    %interpolate over step function and evaluate at Tdata 
    E_HR_Func=griddedInterpolant(T_RRint,HRi,'pchip');
    E_HR = E_HR_Func(Tdata);

    if figureson == 1
        %plot(Tdata,E_HR,'b')
        figure
        scatter(Tdata,E_HR,2)
        title('Estimated HR')
        xlabel('Time')
        ylabel('Beats per minute')
    end

end

%% BP to SBP Function

function [SPdata] = SBPcalc(Tdata,Pdata,pkprom,graphsYoN)
    %Determine SBP signal NOTE FOR CREDIT: CODE ADAPTED FROM BEN RANDELL 
    %Pkprm := MinPeakProminence, usually 25 is good but may need to go lower
    %Run code with output graphs first to eyeball, then comment out. 
    graphs = graphsYoN; %0= no graphs, 1=with graphs

    dt = mean(diff(Tdata));

    [~,sbploc] = findpeaks(Pdata,'MinPeakDistance',round(.25/dt),'MinPeakProminence',pkprom);

    T = [Tdata(1); Tdata(sbploc); Tdata(end)]; %includes first and last time point 
    P = [Pdata(sbploc(1)); Pdata(sbploc); Pdata(sbploc(end))]; 
    SP = griddedInterpolant(T,P,'pchip');
    SPdata = SP(Tdata);


    if graphs == 1
        figure
        plot(Tdata,Pdata)
        hold on
        plot(Tdata,SPdata,'r')
    end


end

end