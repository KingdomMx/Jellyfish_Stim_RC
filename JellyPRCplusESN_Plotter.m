function [] = JellyPRCplusESN_Plotter()
%%% Fuse the Jellyfish Sensor data with an ESN 

close all
try
    %%% If varables are in ws load them here (Saves a lot of time for repeated executions)
    Files= evalin('caller','Files');
    Folders = evalin('caller','Folders');
    ImportedData_Organized= evalin('caller','ImportedData_Organized');
    UniqueInputs= evalin('caller','UniqueInputs');
    InputTestCombo = evalin('caller','InputTestCombo');
    ByJFTestCombo = evalin('caller','ByJFTestCombo');
    TestPerformance = evalin('caller','TestPerformance');
catch
    %%% If those varables arent in ws put them there
   disp('Loading Wave Varables Into Workspace')
   load('Compiled_Stimulated_TestData.mat',"Files","Folders","ImportedData_Organized","UniqueInputs","InputTestCombo","ByJFTestCombo");
   assignin('caller','Files',Files);
   assignin('caller','Folders',Folders);
   assignin('caller','ImportedData_Organized',ImportedData_Organized);
   assignin('caller','UniqueInputs',UniqueInputs);
   assignin('caller','InputTestCombo',InputTestCombo);
   assignin('caller','ByJFTestCombo',ByJFTestCombo);
   %%%%%% Load Sensor Analysis
   load('./PRCandESNData/JellyOpenLoopPRCResults.mat','TestPerformance')
   assignin("base","TestPerformance",TestPerformance);
end

%%%%% Run the Organization of Target data Only once
try
    %%%% If varables are in ws load them here (Saves a lot of time for repeated executions)
    BodyData= evalin('caller','BodyData');
    TargetData = evalin('caller','TargetData');
catch
   %%%% If those varables arent in ws put them there
   disp('Organizing Whole Body Target and Body Data ')
   [BodyData,TargetData] = CollectTargets(ByJFTestCombo,InputTestCombo,Folders,Files,ImportedData_Organized);
   assignin('caller','BodyData',BodyData);
   assignin('caller','TargetData',TargetData);
end


%%%%%%%  Compare unstimulated results
load('./PRCandESNData/JellUnstimESNresults_Scaled.mat',...
    'ESN100N','PRC100N','Hybrid100N','Sequential100N','Rand100N')
NoStimHyb= Hybrid100N;
NoStimSeq =Sequential100N;
NoStimRand = Rand100N;
NoStimESN =ESN100N;
NoStimPrc =PRC100N;


%%%% Load R2 values fit when the input is scaled in proportion to the number of inputs
load('./PRCandESNData/JellFusedESNresults4_Scaled.mat',...
    'Uscale','ESN100N','PRC100N','Hybrid100N','Sequential100N');
ESN100NScaled = ESN100N;
PRC100NScaled = PRC100N;
Sequential100NScaled = Sequential100N;
Hybrid100NScaled = Hybrid100N;
ScalledU = Uscale;

%%%% Load R2 values with no special scaling on the leaky-integrator inputs
load('./PRCandESNData/JellFusedESNresults3.mat','PRC50N','Sequential50N',"Hybrid50N",'ConfigNames',...
    'Sensor_Config','Spect','Uscale','ESN50N','ESN100N','PRC100N','Hybrid100N','Sequential100N');


% VideoName = 'JF4150NVelPredictions.mp4';
VideoName = [];

Plot_Predictions = 0;
Plot_R2 = 0;
Plot_R2_ScaledU = 1;

%%%%%% Adjust how error is computed
if Plot_Predictions  == 1
    N = 50;
    Plot_State_Data(BodyData,TargetData,Sensor_Config(1,:),N,Spect,Uscale,121,VideoName)

    N = 100;
    Plot_State_Data(BodyData,TargetData,Sensor_Config(1,:),N,Spect,Uscale,121,VideoName)
end

%%%%%% Plot R2
if Plot_R2  == 1
    Estimate_Best_Types(ESN50N,PRC50N,Sequential50N,Hybrid50N,ConfigNames,50)
    Plot_R2_Results(ESN50N,PRC50N,Sequential50N,Hybrid50N,ConfigNames,50)

    Estimate_Best_Types(ESN100N,PRC100N,Sequential100N,Hybrid100N,ConfigNames,100)
    Plot_R2_Results(ESN100N,PRC100N,Sequential100N,Hybrid100N,ConfigNames,100)
end


%%%%%% Plot R2 for the Scaled Input
if Plot_R2_ScaledU  == 1
    Estimate_Best_Types(ESN100NScaled,PRC100NScaled,Sequential100NScaled,Hybrid100NScaled,ConfigNames,100)

    Plot_State_Data2(BodyData,TargetData,Sensor_Config(1,:),100,Spect,ScalledU,121)

    Plot_R2_NoStim_Results(NoStimESN,NoStimPrc,NoStimSeq,NoStimHyb,NoStimRand,ConfigNames,100)

    Plot_R2_Results(ESN100NScaled,PRC100NScaled,Sequential100NScaled,Hybrid100NScaled,ConfigNames,100)
end

end

function [BodyData,TargetData] = CollectTargets(ByJFTestCombo,InputTestCombo,Folders,Files,ImportedData_Organized)
%%% Just Compute the Pulsatile and Aggregate Data from the ESN 

fillmethod = "center"; %%% Used to fix outliers 

%%% Use the full pulse
Pulse_Targets = {'COM_X','COM_Y','COM_Z','COM_R1','COM_R2','COM_R3','BODY_VX','BODY_VY','BODY_VZ'};
TargetData.Pulse_Targets = Pulse_Targets;
InputPeriods = [0.5,1,1.5,2];
InputPulseLens = [31,61,91,121];

MRKFLDS = fields(ImportedData_Organized{1, 1}.MkLen);
BodyData.SensorNames = {MRKFLDS{:},'Inner','Outer','Input'};

%%%%% Pulsatile by JF and Input Type
for TND = 1:4 %%% all FRQ
     %%% all JF
    for HND = 1:length(Folders)
        %%%%% Target Data
        for QND = 1:length(Pulse_Targets)
            TMP = ByJFTestCombo{TND, HND}.StimPulse.(Pulse_Targets{QND}).Pulses;
            TMP1D = reshape(TMP,[],1);
            TMP = ByJFTestCombo{TND, HND}.StimPulseSS.(Pulse_Targets{QND}).Pulses;
            TMP1D_S = reshape(TMP,[],1);
            TMP = ByJFTestCombo{TND, HND}.StimPulse.(Pulse_Targets{QND}).Pulses;
            TMP1D_M = reshape(TMP,[],1);
            %%%%% Store By JF Data
            if TND == 1
                TargetData.(erase(Folders{HND},'/')).P.(Pulse_Targets{QND}) = TMP1D';
                TargetData.(erase(Folders{HND},'/')).SS.(Pulse_Targets{QND}) = TMP1D_S';
                TargetData.(erase(Folders{HND},'/')).M.(Pulse_Targets{QND}) = TMP1D_M';
            else
                TargetData.(erase(Folders{HND},'/')).P.(Pulse_Targets{QND}) = [TargetData.(erase(Folders{HND},'/')).P.(Pulse_Targets{QND}),TMP1D'];
                TargetData.(erase(Folders{HND},'/')).SS.(Pulse_Targets{QND}) = [TargetData.(erase(Folders{HND},'/')).SS.(Pulse_Targets{QND}),TMP1D_S'];
                TargetData.(erase(Folders{HND},'/')).M.(Pulse_Targets{QND}) = [TargetData.(erase(Folders{HND},'/')).M.(Pulse_Targets{QND}),TMP1D_M'];
            end

        end

        %%%%%% Body Data 
        for QND2 = 1:length(BodyData.SensorNames)
            TMP = ByJFTestCombo{TND, HND}.StimPulse.(BodyData.SensorNames{QND2}).Pulses;
            TMP1D = reshape(TMP,[],1);
            TMP = ByJFTestCombo{TND, HND}.StimPulseSS.(BodyData.SensorNames{QND2}).Pulses;
            TMP1D_S = reshape(TMP,[],1);
            TMP = ByJFTestCombo{TND, HND}.StimPulse.(BodyData.SensorNames{QND2}).Pulses;
            TMP1D_M = reshape(TMP,[],1);
            %%%%% Store By JF Data
            if TND == 1
                BodyData.(erase(Folders{HND},'/')).P.(BodyData.SensorNames{QND2}) = TMP1D';
                BodyData.(erase(Folders{HND},'/')).SS.(BodyData.SensorNames{QND2}) = TMP1D_S';
                BodyData.(erase(Folders{HND},'/')).M.(BodyData.SensorNames{QND2}) = TMP1D_M';
            else
                BodyData.(erase(Folders{HND},'/')).P.(BodyData.SensorNames{QND2}) = [BodyData.(erase(Folders{HND},'/')).P.(BodyData.SensorNames{QND2}),TMP1D'];
                BodyData.(erase(Folders{HND},'/')).SS.(BodyData.SensorNames{QND2}) = [BodyData.(erase(Folders{HND},'/')).SS.(BodyData.SensorNames{QND2}),TMP1D_S'];
                BodyData.(erase(Folders{HND},'/')).M.(BodyData.SensorNames{QND2}) = [BodyData.(erase(Folders{HND},'/')).M.(BodyData.SensorNames{QND2}),TMP1D_M'];
            end
        end

    end

    %%%% Store Data by Input Type
    for QND = 1:length(Pulse_Targets)
        TMP1D = reshape(InputTestCombo{TND}.StimPulse.(Pulse_Targets{QND}).Pulses,[],1);
        TMP1D_S = reshape(InputTestCombo{TND}.StimPulseSS.(Pulse_Targets{QND}).Pulses,[],1);
        TMP1D_M = reshape(InputTestCombo{TND}.StimPulse_M.(Pulse_Targets{QND}).Pulses,[],1);

        TargetData.(strcat('INPUT_',num2str(TND))).P.(Pulse_Targets{QND}) = TMP1D';
        TargetData.(strcat('INPUT_',num2str(TND))).SS.(Pulse_Targets{QND}) = TMP1D_S';
        TargetData.(strcat('INPUT_',num2str(TND))).M.(Pulse_Targets{QND}) = TMP1D_M';
    end

    for QND2 = 1:length(BodyData.SensorNames)
        TMP1D = reshape(InputTestCombo{TND}.StimPulse.(BodyData.SensorNames{QND2}).Pulses,[],1);
        TMP1D_S = reshape(InputTestCombo{TND}.StimPulseSS.(BodyData.SensorNames{QND2}).Pulses,[],1);
        TMP1D_M = reshape(InputTestCombo{TND}.StimPulse_M.(BodyData.SensorNames{QND2}).Pulses,[],1);

        BodyData.(strcat('INPUT_',num2str(TND))).P.(BodyData.SensorNames{QND2}) = TMP1D';
        BodyData.(strcat('INPUT_',num2str(TND))).SS.(BodyData.SensorNames{QND2}) = TMP1D_S';
        BodyData.(strcat('INPUT_',num2str(TND))).M.(BodyData.SensorNames{QND2}) = TMP1D_M';
    end
end

%%%% Store Data for Unstimulated
for QND = 1:length(Pulse_Targets)
    TMP1D = reshape(InputTestCombo{5}.PulseWise.(Pulse_Targets{QND}).Pulses,[],1);
    TargetData.Unstimulated.(Pulse_Targets{QND}) = TMP1D(~isnan(TMP1D))';
end

for QND2 = 1:length(BodyData.SensorNames)
    TMP1D = reshape(InputTestCombo{5}.PulseWise.(BodyData.SensorNames{QND2}).Pulses,[],1);
    BodyData.Unstimulated.(BodyData.SensorNames{QND2}) = TMP1D(~isnan(TMP1D))';
end

%%%%% Find out corresponding Jelly
CorrJelly=cell(size(ImportedData_Organized));
CorrStim =cell(size(ImportedData_Organized));
for IND = 1:length(Files)
    for INC = 1:length(Folders)
        if contains(Files{IND},Folders{INC})
        CorrJelly{IND} = erase(Folders{INC},'/');
        end
    end

    for INC = 1:length(InputPeriods)
        if ImportedData_Organized{IND}.InputPeriod ==InputPeriods(INC)
        CorrStim{IND} = strcat('INPUT_',num2str(INC));
        end
    end
end

TMPStimNames = {'INPUT_1','INPUT_2','INPUT_3','INPUT_4'};

%%%%% Get full time info here, also extract full time marker lengths by JF 
for IND = 1:length(ImportedData_Organized)
    %%%%% Store mark length data
    CURDATA= ImportedData_Organized{IND};
    PWranges = CURDATA.PwRanges;
    CurJFld = CorrJelly{IND};
    CurStim = CorrStim{IND};
    StimLoc = find(strcmp(TMPStimNames,CurStim));

    %%%%% Run through all Marker Lengths
    for JND = 1:length(MRKFLDS)
        CurMarkerLen = CURDATA.MkLen.(MRKFLDS{JND}).Length';
        CurMarkerLen = filloutliers(CurMarkerLen,fillmethod);
        
        %%%% Store body Marker data
        if IND  == 1
            BodyData.All.All.(MRKFLDS{JND}) = CurMarkerLen;
        else
            BodyData.All.All.(MRKFLDS{JND}) = [BodyData.All.All.(MRKFLDS{JND}),CurMarkerLen];
        end

        %%%% Store body Marker Data Sorted by JF
        if ~isfield(BodyData.All,CorrJelly{IND})
            BodyData.All.(CurJFld).(MRKFLDS{JND}) = CurMarkerLen;
        else
            if ~isfield(BodyData.All.(CurJFld),(MRKFLDS{JND}))
                BodyData.All.(CurJFld).(MRKFLDS{JND}) = CurMarkerLen;
            else
                BodyData.All.(CurJFld).(MRKFLDS{JND}) = [BodyData.All.(CurJFld).(MRKFLDS{JND}),CurMarkerLen];
            end
        end

      %%%% Store body Marker Data Sorted by Unstimulated Tests
      if isempty(CurStim)
          if ~isfield(BodyData.All,'NoStim')
              BodyData.All.NoStim.(MRKFLDS{JND}) = CurMarkerLen;
          else
              if ~isfield(BodyData.All.NoStim,(MRKFLDS{JND}))
                  BodyData.All.NoStim.(MRKFLDS{JND}) = CurMarkerLen;
              else
                  BodyData.All.NoStim.(MRKFLDS{JND}) = [BodyData.All.NoStim.(MRKFLDS{JND}),CurMarkerLen];
              end
          end

          NoStimJell = strcat('NoStim_',CorrJelly{IND});
          if ~isfield(BodyData.All,NoStimJell)
              BodyData.All.(NoStimJell).(MRKFLDS{JND}) = CurMarkerLen;
          else
              if ~isfield(BodyData.All.(NoStimJell),(MRKFLDS{JND}))
                  BodyData.All.(NoStimJell).(MRKFLDS{JND}) = CurMarkerLen;
              else
                  BodyData.All.(NoStimJell).(MRKFLDS{JND}) = [BodyData.All.(NoStimJell).(MRKFLDS{JND}),CurMarkerLen];
              end
          end
      end

    end

    %%%%% Get inner, outer, and input data
    InputStim = CURDATA.input';
    OutLen = filloutliers(CURDATA.Outer.AvgRadius',fillmethod);
    InLen= filloutliers(CURDATA.Inner.AvgRadius',fillmethod);
    %%%% Store body Data Sorted all together
    if IND  == 1
        BodyData.All.All.Inner = InLen;
        BodyData.All.All.Outer = OutLen;
        BodyData.All.All.Input = InputStim;
    else
        BodyData.All.All.Inner = [ BodyData.All.All.Inner,InLen];
        BodyData.All.All.Outer = [BodyData.All.All.Outer,OutLen];
        BodyData.All.All.Input = [BodyData.All.All.Input,InputStim];
    end

    %%%% Store body  Data Sorted by JF
    if ~isfield(BodyData.All.(CurJFld),'Inner')
        BodyData.All.(CurJFld).Inner = InLen;
        BodyData.All.(CurJFld).Outer = OutLen;
        BodyData.All.(CurJFld).Input = InputStim;
    else
        BodyData.All.(CurJFld).Inner = [ BodyData.All.(CurJFld).Inner,InLen];
        BodyData.All.(CurJFld).Outer = [BodyData.All.(CurJFld).Outer,OutLen];
        BodyData.All.(CurJFld).Input = [BodyData.All.(CurJFld).Input,InputStim];
    end

      %%%% Store body Marker Data Sorted by Unstimulated Tests
      if isempty(CurStim)
          if ~isfield(BodyData.All.NoStim,'Inner')
              BodyData.All.NoStim.Inner = InLen;
              BodyData.All.NoStim.Outer = OutLen;
              BodyData.All.NoStim.Input = InputStim;
          else
              BodyData.All.NoStim.Inner = [ BodyData.All.NoStim.Inner,InLen];
              BodyData.All.NoStim.Outer = [BodyData.All.NoStim.Outer,OutLen];
              BodyData.All.NoStim.Input = [BodyData.All.NoStim.Input,InputStim];
          end

          if ~isfield(BodyData.All.(NoStimJell),'Inner')
              BodyData.All.(NoStimJell).Inner = InLen;
              BodyData.All.(NoStimJell).Outer = OutLen;
              BodyData.All.(NoStimJell).Input = InputStim;
          else
              BodyData.All.(NoStimJell).Inner = [BodyData.All.(NoStimJell).Inner,InLen];
              BodyData.All.(NoStimJell).Outer = [BodyData.All.(NoStimJell).Outer,OutLen];
              BodyData.All.(NoStimJell).Input = [BodyData.All.(NoStimJell).Input,InputStim];
          end
      end

    %%%%% Getting full time target Velocity
   CurTestBodyData = CURDATA.BodyFrame;
   Vdata =  CurTestBodyData.Vel;
   Adata =  CurTestBodyData.AccSmooth;
    VX = Vdata(:,1)';
    VY = Vdata(:,2)';
    VZ = Vdata(:,3)';
    AX = Adata(:,1)';
    AY = Adata(:,2)';
    AZ = Adata(:,3)';

    MarkerX = [CurTestBodyData.Top.X,CurTestBodyData.Bottom.X];
    MarkerY = [CurTestBodyData.Top.Y,CurTestBodyData.Bottom.Y];
    MarkerZ = [CurTestBodyData.Top.Z,CurTestBodyData.Bottom.Z];
    TopCenter = CurTestBodyData.Top.Center;
    BottomCenter = CurTestBodyData.Bottom.Center;

    %%%%% Storing fulltime Target Velocity
    if IND  == 1
        TargetData.All.All.VX = VX;
        TargetData.All.All.VY = VY;
        TargetData.All.All.VZ = VZ;
    else
        TargetData.All.All.VX = [ TargetData.All.All.VX,VX];
        TargetData.All.All.VY = [TargetData.All.All.VY,VY];
        TargetData.All.All.VZ = [TargetData.All.All.VZ,VZ];
    end


    %%%%%% Store the Marker data Based on Stimulus and JF
    if ~isempty(CurStim) %%% Check that Stimulus occurs
        %%% Save by JF PW
        if ~isfield(TargetData.(CurJFld).P,'MarkerX')
            pRng=PWranges.StimPulse(1):PWranges.StimPulse(2);
            TargetData.(CurJFld).P.MarkerX =  AddRepData(MarkerX(pRng,:),InputPulseLens(StimLoc));
            TargetData.(CurJFld).P.MarkerY =  AddRepData(MarkerY(pRng,:),InputPulseLens(StimLoc));
            TargetData.(CurJFld).P.MarkerZ =  AddRepData(MarkerZ(pRng,:),InputPulseLens(StimLoc));
            TargetData.(CurJFld).P.CenT =  AddRepData(TopCenter(pRng,:),InputPulseLens(StimLoc));
            TargetData.(CurJFld).P.CenB =  AddRepData(BottomCenter(pRng,:),InputPulseLens(StimLoc));

            pRng=PWranges.StimPulseSS(1):PWranges.StimPulseSS(2);
            TargetData.(CurJFld).SS.MarkerX =  AddRepData(MarkerX(pRng,:),InputPulseLens(StimLoc));
            TargetData.(CurJFld).SS.MarkerY =  AddRepData(MarkerY(pRng,:),InputPulseLens(StimLoc));
            TargetData.(CurJFld).SS.MarkerZ =  AddRepData(MarkerZ(pRng,:),InputPulseLens(StimLoc));
            TargetData.(CurJFld).SS.CenT =  AddRepData(TopCenter(pRng,:),InputPulseLens(StimLoc));
            TargetData.(CurJFld).SS.CenB =  AddRepData(BottomCenter(pRng,:),InputPulseLens(StimLoc));

            pRng=PWranges.StimPulse_M(1):PWranges.StimPulse_M(2);
            TargetData.(CurJFld).M.MarkerX =  AddRepData(MarkerX(pRng,:),InputPulseLens(StimLoc));
            TargetData.(CurJFld).M.MarkerY =  AddRepData(MarkerY(pRng,:),InputPulseLens(StimLoc));
            TargetData.(CurJFld).M.MarkerZ =  AddRepData(MarkerZ(pRng,:),InputPulseLens(StimLoc));
            TargetData.(CurJFld).M.CenT =  AddRepData(TopCenter(pRng,:),InputPulseLens(StimLoc));
            TargetData.(CurJFld).M.CenB = AddRepData(BottomCenter(pRng,:),InputPulseLens(StimLoc));

        else
            pRng=PWranges.StimPulse(1):PWranges.StimPulse(2);
            TargetData.(CurJFld).P.MarkerX =  [TargetData.(CurJFld).P.MarkerX;AddRepData(MarkerX(pRng,:),InputPulseLens(StimLoc))];
            TargetData.(CurJFld).P.MarkerY =  [TargetData.(CurJFld).P.MarkerY;AddRepData(MarkerY(pRng,:),InputPulseLens(StimLoc))];
            TargetData.(CurJFld).P.MarkerZ =  [TargetData.(CurJFld).P.MarkerZ;AddRepData(MarkerZ(pRng,:),InputPulseLens(StimLoc))];
            TargetData.(CurJFld).P.CenT =  [TargetData.(CurJFld).P.CenT;AddRepData(TopCenter(pRng,:),InputPulseLens(StimLoc))];
            TargetData.(CurJFld).P.CenB =  [TargetData.(CurJFld).P.CenB;AddRepData(BottomCenter(pRng,:),InputPulseLens(StimLoc))];

            pRng=PWranges.StimPulseSS(1):PWranges.StimPulseSS(2);
            TargetData.(CurJFld).SS.MarkerX =  [TargetData.(CurJFld).SS.MarkerX;AddRepData(MarkerX(pRng,:),InputPulseLens(StimLoc))];
            TargetData.(CurJFld).SS.MarkerY =  [TargetData.(CurJFld).SS.MarkerY;AddRepData(MarkerY(pRng,:),InputPulseLens(StimLoc))];
            TargetData.(CurJFld).SS.MarkerZ =  [TargetData.(CurJFld).SS.MarkerZ;AddRepData(MarkerZ(pRng,:),InputPulseLens(StimLoc))];
            TargetData.(CurJFld).SS.CenT =  [TargetData.(CurJFld).SS.CenT;AddRepData(TopCenter(pRng,:),InputPulseLens(StimLoc))];
            TargetData.(CurJFld).SS.CenB =  [TargetData.(CurJFld).SS.CenB;AddRepData(BottomCenter(pRng,:),InputPulseLens(StimLoc))];

            pRng=PWranges.StimPulse_M(1):PWranges.StimPulse_M(2);
            TargetData.(CurJFld).M.MarkerX =  [TargetData.(CurJFld).M.MarkerX;AddRepData(MarkerX(pRng,:),InputPulseLens(StimLoc))];
            TargetData.(CurJFld).M.MarkerY =  [TargetData.(CurJFld).M.MarkerY;AddRepData(MarkerY(pRng,:),InputPulseLens(StimLoc))];
            TargetData.(CurJFld).M.MarkerZ =  [TargetData.(CurJFld).M.MarkerZ;AddRepData(MarkerZ(pRng,:),InputPulseLens(StimLoc))];
            TargetData.(CurJFld).M.CenT =  [TargetData.(CurJFld).M.CenT;AddRepData(TopCenter(pRng,:),InputPulseLens(StimLoc))];
            TargetData.(CurJFld).M.CenB =  [TargetData.(CurJFld).M.CenB;AddRepData(BottomCenter(pRng,:),InputPulseLens(StimLoc))];
        end

        %%% Save by INPUT PW
        if ~isfield(TargetData.(CurStim).P,'MarkerX')
            pRng=PWranges.StimPulse(1):PWranges.StimPulse(2);
            TargetData.(CurStim).P.MarkerX =  AddRepData(MarkerX(pRng,:),InputPulseLens(StimLoc));
            TargetData.(CurStim).P.MarkerY =  AddRepData(MarkerY(pRng,:),InputPulseLens(StimLoc));
            TargetData.(CurStim).P.MarkerZ =  AddRepData(MarkerZ(pRng,:),InputPulseLens(StimLoc));
            TargetData.(CurStim).P.CenT =  AddRepData(TopCenter(pRng,:),InputPulseLens(StimLoc));
            TargetData.(CurStim).P.CenB =  AddRepData(BottomCenter(pRng,:),InputPulseLens(StimLoc));

            pRng=PWranges.StimPulseSS(1):PWranges.StimPulseSS(2);
            TargetData.(CurStim).SS.MarkerX =  AddRepData(MarkerX(pRng,:),InputPulseLens(StimLoc));
            TargetData.(CurStim).SS.MarkerY =  AddRepData(MarkerY(pRng,:),InputPulseLens(StimLoc));
            TargetData.(CurStim).SS.MarkerZ =  AddRepData(MarkerZ(pRng,:),InputPulseLens(StimLoc));
            TargetData.(CurStim).SS.CenT =  AddRepData(TopCenter(pRng,:),InputPulseLens(StimLoc));
            TargetData.(CurStim).SS.CenB =  AddRepData(BottomCenter(pRng,:),InputPulseLens(StimLoc));

            pRng=PWranges.StimPulse_M(1):PWranges.StimPulse_M(2);
            TargetData.(CurStim).M.MarkerX =  AddRepData(MarkerX(pRng,:),InputPulseLens(StimLoc));
            TargetData.(CurStim).M.MarkerY =  AddRepData(MarkerY(pRng,:),InputPulseLens(StimLoc));
            TargetData.(CurStim).M.MarkerZ =  AddRepData(MarkerZ(pRng,:),InputPulseLens(StimLoc));
            TargetData.(CurStim).M.CenT =  AddRepData(TopCenter(pRng,:),InputPulseLens(StimLoc));
            TargetData.(CurStim).M.CenB =  AddRepData(BottomCenter(pRng,:),InputPulseLens(StimLoc));

        else
            pRng=PWranges.StimPulse(1):PWranges.StimPulse(2);
            TargetData.(CurStim).P.MarkerX =  [TargetData.(CurStim).P.MarkerX;AddRepData(MarkerX(pRng,:),InputPulseLens(StimLoc))];
            TargetData.(CurStim).P.MarkerY =  [TargetData.(CurStim).P.MarkerY;AddRepData(MarkerY(pRng,:),InputPulseLens(StimLoc))];
            TargetData.(CurStim).P.MarkerZ =  [TargetData.(CurStim).P.MarkerZ;AddRepData(MarkerZ(pRng,:),InputPulseLens(StimLoc))];
            TargetData.(CurStim).P.CenT =  [TargetData.(CurStim).P.CenT;AddRepData(TopCenter(pRng,:),InputPulseLens(StimLoc))];
            TargetData.(CurStim).P.CenB =  [TargetData.(CurStim).P.CenB;AddRepData(BottomCenter(pRng,:),InputPulseLens(StimLoc))];

            pRng=PWranges.StimPulseSS(1):PWranges.StimPulseSS(2);
            TargetData.(CurStim).SS.MarkerX =  [TargetData.(CurStim).SS.MarkerX;AddRepData(MarkerX(pRng,:),InputPulseLens(StimLoc))];
            TargetData.(CurStim).SS.MarkerY =  [TargetData.(CurStim).SS.MarkerY;AddRepData(MarkerY(pRng,:),InputPulseLens(StimLoc))];
            TargetData.(CurStim).SS.MarkerZ =  [TargetData.(CurStim).SS.MarkerZ;AddRepData(MarkerZ(pRng,:),InputPulseLens(StimLoc))];
            TargetData.(CurStim).SS.CenT =  [TargetData.(CurStim).SS.CenT;AddRepData(TopCenter(pRng,:),InputPulseLens(StimLoc))];
            TargetData.(CurStim).SS.CenB =  [TargetData.(CurStim).SS.CenB;AddRepData(BottomCenter(pRng,:),InputPulseLens(StimLoc))];

            pRng=PWranges.StimPulse_M(1):PWranges.StimPulse_M(2);
            TargetData.(CurStim).M.MarkerX =  [TargetData.(CurStim).M.MarkerX;AddRepData(MarkerX(pRng,:),InputPulseLens(StimLoc))];
            TargetData.(CurStim).M.MarkerY =  [TargetData.(CurStim).M.MarkerY;AddRepData(MarkerY(pRng,:),InputPulseLens(StimLoc))];
            TargetData.(CurStim).M.MarkerZ =  [TargetData.(CurStim).M.MarkerZ;AddRepData(MarkerZ(pRng,:),InputPulseLens(StimLoc))];
            TargetData.(CurStim).M.CenT =  [TargetData.(CurStim).M.CenT;AddRepData(TopCenter(pRng,:),InputPulseLens(StimLoc))];
            TargetData.(CurStim).M.CenB =  [TargetData.(CurStim).M.CenB;AddRepData(BottomCenter(pRng,:),InputPulseLens(StimLoc))];
        end
    else
        %%%% Store body  Data Sorted by JF
        if  ~isfield(TargetData.All,'NoStim')
            TargetData.All.NoStim.VX = VX;
            TargetData.All.NoStim.VY = VY;
            TargetData.All.NoStim.VZ = VZ;
            TargetData.All.NoStim.AX = AX;
            TargetData.All.NoStim.AY = AY;
            TargetData.All.NoStim.AZ = AZ;

            TargetData.All.NoStim.MarkerX =  MarkerX;
            TargetData.All.NoStim.MarkerY =  MarkerY;
            TargetData.All.NoStim.MarkerZ =  MarkerZ;
            TargetData.All.NoStim.CenT =  TopCenter;
            TargetData.All.NoStim.CenB =  BottomCenter;

        else
            if ~isfield(TargetData.All.NoStim,'VX')
                TargetData.All.NoStim.VX = VX;
                TargetData.All.NoStim.VY = VY;
                TargetData.All.NoStim.VZ = VZ;
                TargetData.All.NoStim.AX = AX;
                TargetData.All.NoStim.AY = AY;
                TargetData.All.NoStim.AZ = AZ;

                TargetData.All.NoStim.MarkerX =  MarkerX;
                TargetData.All.NoStim.MarkerY =  MarkerY;
                TargetData.All.NoStim.MarkerZ =  MarkerZ;
                TargetData.All.NoStim.CenT =  TopCenter;
                TargetData.All.NoStim.CenB =  BottomCenter;
            else
                TargetData.All.NoStim.VX = [ TargetData.All.NoStim.VX,VX];
                TargetData.All.NoStim.VY = [TargetData.All.NoStim.VY,VY];
                TargetData.All.NoStim.VZ = [TargetData.All.NoStim.VZ,VZ];
                TargetData.All.NoStim.AX = [ TargetData.All.NoStim.AX,AX];
                TargetData.All.NoStim.AY = [TargetData.All.NoStim.AY,AY];
                TargetData.All.NoStim.AZ = [TargetData.All.NoStim.AZ,AZ];

                TargetData.All.NoStim.MarkerX =  [TargetData.All.NoStim.MarkerX;MarkerX];
                TargetData.All.NoStim.MarkerY =  [TargetData.All.NoStim.MarkerY;MarkerY];
                TargetData.All.NoStim.MarkerZ =  [TargetData.All.NoStim.MarkerZ;MarkerZ];
                TargetData.All.NoStim.CenT =  [TargetData.All.NoStim.CenT;TopCenter];
                TargetData.All.NoStim.CenB =  [TargetData.All.NoStim.CenB;BottomCenter];
            end
        end

        if  ~isfield(TargetData.All,NoStimJell)
            TargetData.All.(NoStimJell).VX = VX;
            TargetData.All.(NoStimJell).VY = VY;
            TargetData.All.(NoStimJell).VZ = VZ;
            TargetData.All.(NoStimJell).AX = AX;
            TargetData.All.(NoStimJell).AY = AY;
            TargetData.All.(NoStimJell).AZ = AZ;

            TargetData.All.(NoStimJell).MarkerX =  MarkerX;
            TargetData.All.(NoStimJell).MarkerY =  MarkerY;
            TargetData.All.(NoStimJell).MarkerZ =  MarkerZ;
            TargetData.All.(NoStimJell).CenT =  TopCenter;
            TargetData.All.(NoStimJell).CenB =  BottomCenter;

        else
            if ~isfield(TargetData.All.(NoStimJell),'VX')
                TargetData.All.(NoStimJell).VX = VX;
                TargetData.All.(NoStimJell).VY = VY;
                TargetData.All.(NoStimJell).VZ = VZ;
                TargetData.All.(NoStimJell).AX = AX;
                TargetData.All.(NoStimJell).AY = AY;
                TargetData.All.(NoStimJell).AZ = AZ;

                TargetData.All.(NoStimJell).MarkerX =  MarkerX;
                TargetData.All.(NoStimJell).MarkerY =  MarkerY;
                TargetData.All.(NoStimJell).MarkerZ =  MarkerZ;
                TargetData.All.(NoStimJell).CenT =  TopCenter;
                TargetData.All.(NoStimJell).CenB =  BottomCenter;
            else
                TargetData.All.(NoStimJell).VX = [ TargetData.All.(NoStimJell).VX,VX];
                TargetData.All.(NoStimJell).VY = [TargetData.All.(NoStimJell).VY,VY];
                TargetData.All.(NoStimJell).VZ = [TargetData.All.(NoStimJell).VZ,VZ];
                TargetData.All.(NoStimJell).AX = [ TargetData.All.(NoStimJell).AX,AX];
                TargetData.All.(NoStimJell).AY = [TargetData.All.(NoStimJell).AY,AY];
                TargetData.All.(NoStimJell).AZ = [TargetData.All.(NoStimJell).AZ,AZ];

                TargetData.All.(NoStimJell).MarkerX =  [TargetData.All.(NoStimJell).MarkerX;MarkerX];
                TargetData.All.(NoStimJell).MarkerY =  [TargetData.All.(NoStimJell).MarkerY;MarkerY];
                TargetData.All.(NoStimJell).MarkerZ =  [TargetData.All.(NoStimJell).MarkerZ;MarkerZ];
                TargetData.All.(NoStimJell).CenT =  [TargetData.All.(NoStimJell).CenT;TopCenter];
                TargetData.All.(NoStimJell).CenB =  [TargetData.All.(NoStimJell).CenB;BottomCenter];
            end
        end


    end
    
    %%%% Store body  Data Sorted by JF
    if  ~isfield(TargetData.All,CurJFld)
        TargetData.All.(CurJFld).VX = VX;
        TargetData.All.(CurJFld).VY = VY;
        TargetData.All.(CurJFld).VZ = VZ;
        TargetData.All.(CurJFld).AX = AX;
        TargetData.All.(CurJFld).AY = AY;
        TargetData.All.(CurJFld).AZ = AZ;

        TargetData.All.(CurJFld).MarkerX =  MarkerX;
        TargetData.All.(CurJFld).MarkerY =  MarkerY;
        TargetData.All.(CurJFld).MarkerZ =  MarkerZ;
        TargetData.All.(CurJFld).CenT =  TopCenter;
        TargetData.All.(CurJFld).CenB =  BottomCenter;

    else
        if ~isfield(TargetData.All.(CurJFld),'VX')
            TargetData.All.(CurJFld).VX = VX;
            TargetData.All.(CurJFld).VY = VY;
            TargetData.All.(CurJFld).VZ = VZ;
            TargetData.All.(CurJFld).AX = AX;
            TargetData.All.(CurJFld).AY = AY;
            TargetData.All.(CurJFld).AZ = AZ;

            TargetData.All.(CurJFld).MarkerX =  MarkerX;
            TargetData.All.(CurJFld).MarkerY =  MarkerY;
            TargetData.All.(CurJFld).MarkerZ =  MarkerZ;
            TargetData.All.(CurJFld).CenT =  TopCenter;
            TargetData.All.(CurJFld).CenB =  BottomCenter;
        else
            TargetData.All.(CurJFld).VX = [ TargetData.All.(CurJFld).VX,VX];
            TargetData.All.(CurJFld).VY = [TargetData.All.(CurJFld).VY,VY];
            TargetData.All.(CurJFld).VZ = [TargetData.All.(CurJFld).VZ,VZ];
            TargetData.All.(CurJFld).AX = [ TargetData.All.(CurJFld).AX,AX];
            TargetData.All.(CurJFld).AY = [TargetData.All.(CurJFld).AY,AY];
            TargetData.All.(CurJFld).AZ = [TargetData.All.(CurJFld).AZ,AZ];

            TargetData.All.(CurJFld).MarkerX =  [TargetData.All.(CurJFld).MarkerX;MarkerX];
            TargetData.All.(CurJFld).MarkerY =  [TargetData.All.(CurJFld).MarkerY;MarkerY];
            TargetData.All.(CurJFld).MarkerZ =  [TargetData.All.(CurJFld).MarkerZ;MarkerZ];
            TargetData.All.(CurJFld).CenT =  [TargetData.All.(CurJFld).CenT;TopCenter];
            TargetData.All.(CurJFld).CenB =  [TargetData.All.(CurJFld).CenB;BottomCenter];
        end
    end
end

%%%%%%%% Normalize the Bodies internal Data and center  [-1,1]
bdFldsA = fields(BodyData.All);
bdFldsB = fields(BodyData.All.All);
for INC = 1:length(bdFldsA)
    for JNC = 1:length(bdFldsB)
        CurBDat = BodyData.All.(bdFldsA{INC}).(bdFldsB{JNC});
        BodyData.All_Norm.(bdFldsA{INC}).(bdFldsB{JNC}) =2*(CurBDat-range(CurBDat)/2-min(CurBDat))/range(CurBDat);
    end
end

end

function [Outputdata] = AddRepData(InputData,Lengths)
%%%%% Needed becuase of the way pulses reset in pulsatile data

CurDataLength = length(InputData(:,1));
numPulses = floor((CurDataLength-1)/(Lengths-1));
NewDataLength  = CurDataLength+numPulses;
Outputdata = zeros([NewDataLength-1,length((InputData(1,:)))]);

StartPts = 1:Lengths:NewDataLength;
StartPts2 = 1:(Lengths-1):CurDataLength;

Outputdata(1:Lengths,:) = (InputData(1:Lengths,:));

for IND = 2:length(StartPts(1:end-1))
    Outputdata(StartPts(IND)+(0:Lengths-1),:) = InputData(StartPts2(IND):StartPts2(IND+1),:);
end


end

function [MultData,varargout] = MultiplexData(Data,CycleLength)
%%%%% Store data that is shifted by time and the length
MultData = zeros([CycleLength,length(Data)-CycleLength]);
    for ind = 1:CycleLength
        %%% Cut off end
        MultData(ind,:) =  Data(ind:end-(CycleLength-(ind-1)));
    end

    if nargout >= 2
        TMP = MultData;
        TMP(1,:) = 0; 
        varargout{1}= cumsum(TMP,1);
    end

end

function [] = Plot_R2_Results(ESN,PRC,Sequential,Hybrid,ConfigNames,N)
load('MyColors.Mat', 'FavoriteColors');

Real2SampleTime = 2/120;
XData = (0:120)*Real2SampleTime;
YData = (0:120)*Real2SampleTime;

LabAll = strcat('All Data N',num2str(N));
LabPulse = strcat('Pulsatile N',num2str(N));
LabX = 'Prediction Future Horizon';
LabY = 'Leaky-Int. History Horizon';
setclim = [-0.07,0.07];

JFIDs = {'JF41','JF42','JF43','JF44','JF46'};
InputIDs = {'INPUT_1','INPUT_2','INPUT_3','INPUT_4'};
PulsR2Fields = fields(Hybrid{1}.JF41);
PulsPosFields =PulsR2Fields(1:3);
PulsRotFields =PulsR2Fields(4:6);
PulsVelFields =PulsR2Fields(7:9);

%%%%% Get Beakdown By Jellyfsih
Mean_BYJelly = cell(length(Hybrid),length(JFIDs));
Agg_BYJelly = cell(size(Hybrid)); %%% using data trained by JF for each sensor config
Agg_BYJelly_PRC = cell(size(Hybrid));
Agg_BYJelly_SEQ = cell(size(Hybrid));
for IND = 1:length(Hybrid)
    Agg_BYJelly{IND} = zeros([121,121]);
    for JND = 1:length(JFIDs)
            Mean_BYJelly{IND,JND} = zeros([121,121]);
        for KND = 1:length(PulsR2Fields)
            Mean_BYJelly{IND,JND} = Mean_BYJelly{IND,JND} + Hybrid{IND}.(JFIDs{JND}).(PulsR2Fields{KND});
        end
        Mean_BYJelly{IND,JND} = Mean_BYJelly{IND,JND}/length(PulsR2Fields);
        Agg_BYJelly{IND}  = Agg_BYJelly{IND} + Mean_BYJelly{IND,JND};
    end
    Agg_BYJelly{IND} = Agg_BYJelly{IND}/length(Hybrid);

    Agg_BYJelly_PRC{IND} = zeros([121,121]);
    for JND = 1:length(JFIDs)
            TMP= zeros([121,121]);
        for KND = 1:length(PulsR2Fields)
            TMP = TMP + PRC{IND}.(JFIDs{JND}).(PulsR2Fields{KND});
        end
        TMP = TMP/length(PulsR2Fields);
        Agg_BYJelly_PRC{IND}  = Agg_BYJelly_PRC{IND} + TMP;
    end
    Agg_BYJelly_PRC{IND} = Agg_BYJelly_PRC{IND}/length(Hybrid);

    Agg_BYJelly_SEQ{IND} = zeros([121,121]);
    for JND = 1:length(JFIDs)
            TMP= zeros([121,121]);
        for KND = 1:length(PulsR2Fields)
            TMP = TMP + Sequential{IND}.(JFIDs{JND}).(PulsR2Fields{KND});
        end
        TMP = TMP/length(PulsR2Fields);
        Agg_BYJelly_SEQ{IND}  = Agg_BYJelly_SEQ{IND} + TMP;
    end
    Agg_BYJelly_SEQ{IND} = Agg_BYJelly_SEQ{IND}/length(Hybrid);

end

Mean_BYJellyType = cell(3,length(JFIDs));
for JND = 1:length(JFIDs)
    Mean_BYJellyType{1,JND} = zeros([121,121]);
    for KND = 1:length(PulsPosFields)
        Mean_BYJellyType{1,JND} = Mean_BYJellyType{1,JND} + Hybrid{1}.(JFIDs{JND}).(PulsPosFields{KND});
    end
    Mean_BYJellyType {1,JND} = Mean_BYJellyType{1,JND}/3;

    Mean_BYJellyType{2,JND} = zeros([121,121]);
    for KND = 1:length(PulsRotFields)
        Mean_BYJellyType{2,JND} = Mean_BYJellyType{2,JND} + Hybrid{1}.(JFIDs{JND}).(PulsRotFields{KND});
    end
    Mean_BYJellyType {2,JND} = Mean_BYJellyType{2,JND}/3;


    Mean_BYJellyType{3,JND} = zeros([121,121]);
    for KND = 1:length(PulsVelFields)
        Mean_BYJellyType{3,JND} = Mean_BYJellyType{3,JND} + Hybrid{1}.(JFIDs{JND}).(PulsVelFields{KND});
    end
    Mean_BYJellyType {3,JND} = Mean_BYJellyType{3,JND}/3;
end

Mean_BYInputType = cell(3,length(InputIDs));
for JND = 1:length(InputIDs)
    Mean_BYInputType{1,JND} = zeros([121,121]);
    for KND = 1:length(PulsPosFields)
        Mean_BYInputType{1,JND} = Mean_BYInputType{1,JND} + Hybrid{1}.(InputIDs{JND}).(PulsPosFields{KND});
    end
    Mean_BYInputType {1,JND} = Mean_BYInputType{1,JND}/3;

    Mean_BYInputType{2,JND} = zeros([121,121]);
    for KND = 1:length(PulsRotFields)
        Mean_BYInputType{2,JND} = Mean_BYInputType{2,JND} + Hybrid{1}.(InputIDs{JND}).(PulsRotFields{KND});
    end
    Mean_BYInputType {2,JND} = Mean_BYInputType{2,JND}/3;


    Mean_BYInputType{3,JND} = zeros([121,121]);
    for KND = 1:length(PulsVelFields)
        Mean_BYInputType{3,JND} = Mean_BYInputType{3,JND} + Hybrid{1}.(InputIDs{JND}).(PulsVelFields{KND});
    end
    Mean_BYInputType {3,JND} = Mean_BYInputType{3,JND}/3;
end

%%%%%%%%%%
figure
subplot(3,4,1)
contourf(XData,YData,ESN.All.VX_R2)
title('ESN')
ylabel({'VX',LabY})
clim(setclim)
subplot(3,4,2)
contourf(XData,YData,PRC{1}.All.VX_R2)
title('PRC')
clim(setclim)
subplot(3,4,3)
contourf(XData,YData,Sequential{1}.All.VX_R2)
title('Sequential')
clim(setclim)
subplot(3,4,4)
contourf(XData,YData,Hybrid{1}.All.VX_R2)
title('Hybrid')
clim(setclim)

subplot(3,4,1+4)
contourf(XData,YData,ESN.All.VY_R2)
ylabel({'VY',LabY})
clim(setclim)
subplot(3,4,2+4)
contourf(XData,YData,PRC{1}.All.VY_R2)
clim(setclim)
subplot(3,4,3+4)
contourf(XData,YData,Sequential{1}.All.VY_R2)
clim(setclim)
subplot(3,4,4+4)
contourf(XData,YData,Hybrid{1}.All.VY_R2)
clim(setclim)

subplot(3,4,1+8)
contourf(XData,YData,ESN.All.VZ_R2)
ylabel({'VZ',LabY})
xlabel(LabX)
clim(setclim)
subplot(3,4,2+8)
contourf(XData,YData,PRC{1}.All.VZ_R2)
xlabel(LabX)
clim(setclim)
subplot(3,4,3+8)
contourf(XData,YData,Sequential{1}.All.VZ_R2)
xlabel(LabX)
clim(setclim)
subplot(3,4,4+8)
contourf(XData,YData,Hybrid{1}.All.VZ_R2)
xlabel(LabX)
clim(setclim)
colormap("jet")
sgtitle({LabAll,'Velocity'})
set(gcf,"Position",[239         164        1592         733])


%%%%%%%%%%%%%
setclim = [0,1];
figure
subplot(3,4,1)
contourf(XData,YData,ESN.JF41.P_COM_X_R2)
ylabel({'ComX',LabY})
title('ESN')
clim(setclim)
subplot(3,4,2)
contourf(XData,YData,PRC{1}.JF41.P_COM_X_R2)
title('PRC')
clim(setclim)
subplot(3,4,3)
contourf(XData,YData,Sequential{1}.JF41.P_COM_X_R2)
title('Sequential')
clim(setclim)
subplot(3,4,4)
contourf(XData,YData,Hybrid{1}.JF41.P_COM_X_R2)
title('Hybrid')
clim(setclim)

subplot(3,4,1+4)
contourf(XData,YData,ESN.JF41.P_COM_Y_R2)
ylabel({'ComY',LabY})
clim(setclim)
subplot(3,4,2+4)
contourf(XData,YData,PRC{1}.JF41.P_COM_Y_R2)
clim(setclim)
subplot(3,4,3+4)
contourf(XData,YData,Sequential{1}.JF41.P_COM_Y_R2)
clim(setclim)
subplot(3,4,4+4)
contourf(XData,YData,Hybrid{1}.JF41.P_COM_Y_R2)
clim(setclim)

subplot(3,4,1+8)
contourf(XData,YData,ESN.JF41.P_COM_Z_R2)
ylabel({'ComZ',LabY})
xlabel(LabX)
clim(setclim)
subplot(3,4,2+8)
contourf(XData,YData,PRC{1}.JF41.P_COM_Z_R2)
xlabel(LabX)
clim(setclim)
subplot(3,4,3+8)
contourf(XData,YData,Sequential{1}.JF41.P_COM_Z_R2)
xlabel(LabX)
clim(setclim)
subplot(3,4,4+8)
contourf(XData,YData,Hybrid{1}.JF41.P_COM_Z_R2)
xlabel(LabX)
clim(setclim)
colormap("jet")
sgtitle({LabPulse,'JF41 COM'})
set(gcf,"Position",[239         164        1592         733])


%%%%%%%%%%%%%%%%
setclim = [0,1];
figure
subplot(3,4,1)
contourf(XData,YData,ESN.JF41.P_BODY_VX_R2)
title('ESN')
ylabel({'VX',LabY})
clim(setclim)
subplot(3,4,2)
contourf(XData,YData,PRC{1}.JF41.P_BODY_VX_R2)
title('PRC')
clim(setclim)
subplot(3,4,3)
contourf(XData,YData,Sequential{1}.JF41.P_BODY_VX_R2)
title('Sequential')
clim(setclim)
subplot(3,4,4)
contourf(XData,YData,Hybrid{1}.JF41.P_BODY_VX_R2)
title('Hybrid')
clim(setclim)

subplot(3,4,1+4)
contourf(XData,YData,ESN.JF41.P_BODY_VY_R2)
ylabel({'VY',LabY})
clim(setclim)
subplot(3,4,2+4)
contourf(XData,YData,PRC{1}.JF41.P_BODY_VY_R2)
clim(setclim)
subplot(3,4,3+4)
contourf(XData,YData,Sequential{1}.JF41.P_BODY_VY_R2)
clim(setclim)
subplot(3,4,4+4)
contourf(XData,YData,Hybrid{1}.JF41.P_BODY_VY_R2)
clim(setclim)

subplot(3,4,1+8)
contourf(XData,YData,ESN.JF41.P_BODY_VZ_R2)
ylabel({'VZ',LabY})
xlabel(LabX)
clim(setclim)
subplot(3,4,2+8)
contourf(XData,YData,PRC{1}.JF41.P_BODY_VZ_R2)
xlabel(LabX)
clim(setclim)
subplot(3,4,3+8)
contourf(XData,YData,Sequential{1}.JF41.P_BODY_VZ_R2)
xlabel(LabX)
clim(setclim)
subplot(3,4,4+8)
contourf(XData,YData,Hybrid{1}.JF41.P_BODY_VZ_R2)
xlabel(LabX)
clim(setclim)
colormap("jet")
sgtitle({LabPulse,'JF41 Velocity'})
set(gcf,"Position",[239         164        1592         733])



figure
subplot(3,1,1)
plot(XData,ESN.JF41.P_COM_X_R2(1,:),'Color',FavoriteColors.PptFrm.DarkYellow,'LineWidth',2)
hold on
plot(XData,PRC{1}.JF41.P_COM_X_R2(1,:),'Color',FavoriteColors.PptFrm.DarkPink,'LineWidth',2)
plot(XData,Sequential{1}.JF41.P_COM_X_R2(1,:),'Color',FavoriteColors.PptFrm.DarkBlue,'LineWidth',2)
plot(XData,Hybrid{1}.JF41.P_COM_X_R2(1,:),'Color',FavoriteColors.PptFrm.Green,'LineWidth',2)
plot(XData,ESN.JF41.P_COM_X_R2(end,:),'--','Color',FavoriteColors.PptFrm.DarkYellow,'LineWidth',2)
plot(XData,PRC{1}.JF41.P_COM_X_R2(end,:),'--','Color',FavoriteColors.PptFrm.DarkPink,'LineWidth',2)
plot(XData,Sequential{1}.JF41.P_COM_X_R2(end,:),'--','Color',FavoriteColors.PptFrm.DarkBlue,'LineWidth',2)
plot(XData,Hybrid{1}.JF41.P_COM_X_R2(end,:),'--','Color',FavoriteColors.PptFrm.Green,'LineWidth',2)
ylim([0,1])
ylabel({'VX','R^2'})
subplot(3,1,2)
plot(XData,ESN.JF41.P_COM_Y_R2(1,:),'Color',FavoriteColors.PptFrm.DarkYellow,'LineWidth',2)
hold on
plot(XData,PRC{1}.JF41.P_COM_Y_R2(1,:),'Color',FavoriteColors.PptFrm.DarkPink,'LineWidth',2)
plot(XData,Sequential{1}.JF41.P_COM_Y_R2(1,:),'Color',FavoriteColors.PptFrm.DarkBlue,'LineWidth',2)
plot(XData,Hybrid{1}.JF41.P_COM_Y_R2(1,:),'Color',FavoriteColors.PptFrm.Green,'LineWidth',2)
plot(XData,ESN.JF41.P_COM_Y_R2(end,:),'--','Color',FavoriteColors.PptFrm.DarkYellow,'LineWidth',2)
plot(XData,PRC{1}.JF41.P_COM_Y_R2(end,:),'--','Color',FavoriteColors.PptFrm.DarkPink,'LineWidth',2)
plot(XData,Sequential{1}.JF41.P_COM_Y_R2(end,:),'--','Color',FavoriteColors.PptFrm.DarkBlue,'LineWidth',2)
plot(XData,Hybrid{1}.JF41.P_COM_Y_R2(end,:),'--','Color',FavoriteColors.PptFrm.Green,'LineWidth',2)
ylim([0,1])
ylabel({'VY','R^2'})
subplot(3,1,3)
plot(XData,ESN.JF41.P_COM_Z_R2(1,:),'Color',FavoriteColors.PptFrm.DarkYellow,'LineWidth',2)
hold on
plot(XData,PRC{1}.JF41.P_COM_Z_R2(1,:),'Color',FavoriteColors.PptFrm.DarkPink,'LineWidth',2)
plot(XData,Sequential{1}.JF41.P_COM_Z_R2(1,:),'Color',FavoriteColors.PptFrm.DarkBlue,'LineWidth',2)
plot(XData,Hybrid{1}.JF41.P_COM_Z_R2(1,:),'Color',FavoriteColors.PptFrm.Green,'LineWidth',2)
plot(XData,ESN.JF41.P_COM_Z_R2(end,:),'--','Color',FavoriteColors.PptFrm.DarkYellow,'LineWidth',2)
plot(XData,PRC{1}.JF41.P_COM_Z_R2(end,:),'--','Color',FavoriteColors.PptFrm.DarkPink,'LineWidth',2)
plot(XData,Sequential{1}.JF41.P_COM_Z_R2(end,:),'--','Color',FavoriteColors.PptFrm.DarkBlue,'LineWidth',2)
plot(XData,Hybrid{1}.JF41.P_COM_Z_R2(end,:),'--','Color',FavoriteColors.PptFrm.Green,'LineWidth',2)
ylim([0,1])
xlabel(LabX)
ylabel({'VZ','R^2'})
legend({'Pure u(t_0)','PRC u(t_0)','Sequential u(t_0)','Hybrid u(t_0)','Pure u(t_0-t_2)','PRC u(t_0-t_2)','Sequential u(t_0-t_2)','Hybrid u(t_0-t_2)'},'NumColumns',2,'Location','northeast')
sgtitle({LabPulse,'JF41 COM R^2'})
set(gcf,"Position",[ 455         131        1070         785])



figure
subplot(3,1,1)
plot(XData,ESN.JF41.P_BODY_VX_R2(1,:),'Color',FavoriteColors.PptFrm.DarkYellow,'LineWidth',2)
hold on
plot(XData,PRC{1}.JF41.P_BODY_VX_R2(1,:),'Color',FavoriteColors.PptFrm.DarkPink,'LineWidth',2)
plot(XData,Sequential{1}.JF41.P_BODY_VX_R2(1,:),'Color',FavoriteColors.PptFrm.DarkBlue,'LineWidth',2)
plot(XData,Hybrid{1}.JF41.P_BODY_VX_R2(1,:),'Color',FavoriteColors.PptFrm.Green,'LineWidth',2)
plot(XData,ESN.JF41.P_BODY_VX_R2(end,:),'--','Color',FavoriteColors.PptFrm.DarkYellow,'LineWidth',2)
plot(XData,PRC{1}.JF41.P_BODY_VX_R2(end,:),'--','Color',FavoriteColors.PptFrm.DarkPink,'LineWidth',2)
plot(XData,Sequential{1}.JF41.P_BODY_VX_R2(end,:),'--','Color',FavoriteColors.PptFrm.DarkBlue,'LineWidth',2)
plot(XData,Hybrid{1}.JF41.P_BODY_VX_R2(end,:),'--','Color',FavoriteColors.PptFrm.Green,'LineWidth',2)
ylim([0,1])
ylabel({'VX','R^2'})
subplot(3,1,2)
plot(XData,ESN.JF41.P_BODY_VY_R2(1,:),'Color',FavoriteColors.PptFrm.DarkYellow,'LineWidth',2)
hold on
plot(XData,PRC{1}.JF41.P_BODY_VY_R2(1,:),'Color',FavoriteColors.PptFrm.DarkPink,'LineWidth',2)
plot(XData,Sequential{1}.JF41.P_BODY_VY_R2(1,:),'Color',FavoriteColors.PptFrm.DarkBlue,'LineWidth',2)
plot(XData,Hybrid{1}.JF41.P_BODY_VY_R2(1,:),'Color',FavoriteColors.PptFrm.Green,'LineWidth',2)
plot(XData,ESN.JF41.P_BODY_VY_R2(end,:),'--','Color',FavoriteColors.PptFrm.DarkYellow,'LineWidth',2)
plot(XData,PRC{1}.JF41.P_BODY_VY_R2(end,:),'--','Color',FavoriteColors.PptFrm.DarkPink,'LineWidth',2)
plot(XData,Sequential{1}.JF41.P_BODY_VY_R2(end,:),'--','Color',FavoriteColors.PptFrm.DarkBlue,'LineWidth',2)
plot(XData,Hybrid{1}.JF41.P_BODY_VY_R2(end,:),'--','Color',FavoriteColors.PptFrm.Green,'LineWidth',2)
ylim([0,1])
ylabel({'VY','R^2'})
subplot(3,1,3)
plot(XData,ESN.JF41.P_BODY_VZ_R2(1,:),'Color',FavoriteColors.PptFrm.DarkYellow,'LineWidth',2)
hold on
plot(XData,PRC{1}.JF41.P_BODY_VZ_R2(1,:),'Color',FavoriteColors.PptFrm.DarkPink,'LineWidth',2)
plot(XData,Sequential{1}.JF41.P_BODY_VZ_R2(1,:),'Color',FavoriteColors.PptFrm.DarkBlue,'LineWidth',2)
plot(XData,Hybrid{1}.JF41.P_BODY_VZ_R2(1,:),'Color',FavoriteColors.PptFrm.Green,'LineWidth',2)
plot(XData,ESN.JF41.P_BODY_VZ_R2(end,:),'--','Color',FavoriteColors.PptFrm.DarkYellow,'LineWidth',2)
plot(XData,PRC{1}.JF41.P_BODY_VZ_R2(end,:),'--','Color',FavoriteColors.PptFrm.DarkPink,'LineWidth',2)
plot(XData,Sequential{1}.JF41.P_BODY_VZ_R2(end,:),'--','Color',FavoriteColors.PptFrm.DarkBlue,'LineWidth',2)
plot(XData,Hybrid{1}.JF41.P_BODY_VZ_R2(end,:),'--','Color',FavoriteColors.PptFrm.Green,'LineWidth',2)
ylim([0,1])
xlabel(LabX)
ylabel({'VZ','R^2'})
legend({'Pure u(t_0)','PRC u(t_0)','Sequential u(t_0)','Hybrid u(t_0)','Pure u(t_0-t_2)','PRC u(t_0-t_2)','Sequential u(t_0-t_2)','Hybrid u(t_0-t_2)'},'NumColumns',2,'Location','northeast')
sgtitle({LabPulse,'JF41 Velocity R^2'})
set(gcf,"Position",[ 455         131        1070         785])



figure
subplot(3,1,1)
plot(XData,ESN.INPUT_4.P_BODY_VX_R2(1,:),'Color',FavoriteColors.PptFrm.DarkYellow,'LineWidth',2)
hold on
plot(XData,PRC{1}.INPUT_4.P_BODY_VX_R2(1,:),'Color',FavoriteColors.PptFrm.DarkPink,'LineWidth',2)
plot(XData,Sequential{1}.INPUT_4.P_BODY_VX_R2(1,:),'Color',FavoriteColors.PptFrm.DarkBlue,'LineWidth',2)
plot(XData,Hybrid{1}.INPUT_4.P_BODY_VX_R2(1,:),'Color',FavoriteColors.PptFrm.Green,'LineWidth',2)
plot(XData,ESN.INPUT_4.P_BODY_VX_R2(end,:),'--','Color',FavoriteColors.PptFrm.DarkYellow,'LineWidth',2)
plot(XData,PRC{1}.INPUT_4.P_BODY_VX_R2(end,:),'--','Color',FavoriteColors.PptFrm.DarkPink,'LineWidth',2)
plot(XData,Sequential{1}.INPUT_4.P_BODY_VX_R2(end,:),'--','Color',FavoriteColors.PptFrm.DarkBlue,'LineWidth',2)
plot(XData,Hybrid{1}.INPUT_4.P_BODY_VX_R2(end,:),'--','Color',FavoriteColors.PptFrm.Green,'LineWidth',2)
ylim([0,1])
ylabel({'VX','R^2'})
subplot(3,1,2)
plot(XData,ESN.INPUT_4.P_BODY_VY_R2(1,:),'Color',FavoriteColors.PptFrm.DarkYellow,'LineWidth',2)
hold on
plot(XData,PRC{1}.INPUT_4.P_BODY_VY_R2(1,:),'Color',FavoriteColors.PptFrm.DarkPink,'LineWidth',2)
plot(XData,Sequential{1}.INPUT_4.P_BODY_VY_R2(1,:),'Color',FavoriteColors.PptFrm.DarkBlue,'LineWidth',2)
plot(XData,Hybrid{1}.INPUT_4.P_BODY_VY_R2(1,:),'Color',FavoriteColors.PptFrm.Green,'LineWidth',2)
plot(XData,ESN.INPUT_4.P_BODY_VY_R2(end,:),'--','Color',FavoriteColors.PptFrm.DarkYellow,'LineWidth',2)
plot(XData,PRC{1}.INPUT_4.P_BODY_VY_R2(end,:),'--','Color',FavoriteColors.PptFrm.DarkPink,'LineWidth',2)
plot(XData,Sequential{1}.INPUT_4.P_BODY_VY_R2(end,:),'--','Color',FavoriteColors.PptFrm.DarkBlue,'LineWidth',2)
plot(XData,Hybrid{1}.INPUT_4.P_BODY_VY_R2(end,:),'--','Color',FavoriteColors.PptFrm.Green,'LineWidth',2)
ylim([0,1])
ylabel({'VY','R^2'})
subplot(3,1,3)
plot(XData,ESN.INPUT_4.P_BODY_VZ_R2(1,:),'Color',FavoriteColors.PptFrm.DarkYellow,'LineWidth',2)
hold on
plot(XData,PRC{1}.INPUT_4.P_BODY_VZ_R2(1,:),'Color',FavoriteColors.PptFrm.DarkPink,'LineWidth',2)
plot(XData,Sequential{1}.INPUT_4.P_BODY_VZ_R2(1,:),'Color',FavoriteColors.PptFrm.DarkBlue,'LineWidth',2)
plot(XData,Hybrid{1}.INPUT_4.P_BODY_VZ_R2(1,:),'Color',FavoriteColors.PptFrm.Green,'LineWidth',2)
plot(XData,ESN.INPUT_4.P_BODY_VZ_R2(end,:),'--','Color',FavoriteColors.PptFrm.DarkYellow,'LineWidth',2)
plot(XData,PRC{1}.INPUT_4.P_BODY_VZ_R2(end,:),'--','Color',FavoriteColors.PptFrm.DarkPink,'LineWidth',2)
plot(XData,Sequential{1}.INPUT_4.P_BODY_VZ_R2(end,:),'--','Color',FavoriteColors.PptFrm.DarkBlue,'LineWidth',2)
plot(XData,Hybrid{1}.INPUT_4.P_BODY_VZ_R2(end,:),'--','Color',FavoriteColors.PptFrm.Green,'LineWidth',2)
ylim([0,1])
xlabel(LabX)
ylabel({'VZ','R^2'})
legend({'Pure u(t_0)','PRC u(t_0)','Sequential u(t_0)','Hybrid u(t_0)','Pure u(t_0-t_2)','PRC u(t_0-t_2)','Sequential u(t_0-t_2)','Hybrid u(t_0-t_2)'},'NumColumns',2,'Location','northeast')
sgtitle({LabPulse,'\tau = 2 Velocity R^2'})
set(gcf,"Position",[ 455         131        1070         785])


%%%%%%%%%%%%%%
figure
contourf(XData,YData,(Hybrid{1}.All.VX_R2+Hybrid{1}.All.VY_R2+Hybrid{1}.All.VZ_R2)/3)
xlabel(LabX)
clim([0,0.5])
colormap("jet")
sgtitle({LabAll,'Average R2 by Sensor Trained on all Vel. data'})

%%%%%%%%%%%%%%
setclim = [0,1];
figure
for IND = 1:length(Hybrid)
    subplot(1,length(Hybrid),IND)
    contourf(XData,YData,Agg_BYJelly{IND})
    xlabel(LabX)
    clim(setclim)
    title(ConfigNames{IND})
    if IND == 1
        ylabel(LabY)
    end
end
colormap("jet")
sgtitle({LabPulse,'Average R2 by Sensor Trained by JF:All'})
set(gcf,"Position",[104         406        1708         318])

%%%%%%%%%%%%%%%%
setclim = [0,1];
figure
for IND = 1:length(Hybrid)
    subplot(2,length(Hybrid),IND)
    contourf(XData,YData,Agg_BYJelly_PRC{IND})
    xlabel(LabX)
    clim(setclim)
    title(ConfigNames{IND})
    if IND == 1
        ylabel(LabY)
    end

    subplot(2,length(Hybrid),IND+length(Hybrid))
    contourf(XData,YData,Agg_BYJelly_SEQ{IND})
    xlabel(LabX)
    clim(setclim)
    if IND == 1
        ylabel(LabY)
    end
end
colormap("jet")
sgtitle({LabPulse,'Average R2 by Sensor Trained by JF:PRC and Sequential'})
set(gcf,"Position",[104         406-318/2        1708         318*2])


%%%%%%%%%%%%%%%%
setclim = [0,1];
figure
for JND = 1:length(JFIDs)
    subplot(3,length(JFIDs),JND)
    contourf(XData,YData,Mean_BYJellyType{1,JND})
    xlabel(LabX)
    clim(setclim)
    title(JFIDs{JND})
    if JND == 1
        ylabel({'Local Pos',LabY})
    end

    subplot(3,length(JFIDs),JND+length(JFIDs))
    contourf(XData,YData,Mean_BYJellyType{2,JND})
    xlabel(LabX)
    clim(setclim)
    if JND == 1
        ylabel({'Rotation',LabY})
    end

    subplot(3,length(JFIDs),JND+2*length(JFIDs))
    contourf(XData,YData,Mean_BYJellyType{3,JND})
    xlabel(LabX)
    clim(setclim)
    if JND == 1
        ylabel({'Velocity',LabY})
    end
end
sgtitle({LabPulse,'Average R2 Trained by JF'})
colormap("jet")
set(gcf,"Position",[135         133        1708         804])



setclim = [0,1];
figure
for JND = 1:length(InputIDs)
    subplot(3,length(InputIDs),JND)
    contourf(XData,YData,Mean_BYInputType{1,JND})
    xlabel(LabX)
    clim(setclim)
    title(InputIDs{JND})
    if JND == 1
        ylabel({'Loca lPos',LabY})
    end

   subplot(3,length(InputIDs),JND+4)
    contourf(XData,YData,Mean_BYInputType{2,JND})
    xlabel(LabX)
    clim(setclim)
    if JND == 1
        ylabel({'Rotation',LabY})
    end

   subplot(3,length(InputIDs),JND+8)
    contourf(XData,YData,Mean_BYInputType{3,JND})
    xlabel(LabX)
    clim(setclim)
    if JND == 1
        ylabel({'Velocity',LabY})
    end
end
sgtitle({LabPulse,'Average R2 Trained by Input Period'})
colormap("jet")
set(gcf,"Position",[239         164        1592         733])



%%%%%%%%%%%%%%%%
setclim = [0,0.2];
figure
for JND = 1:length(JFIDs)
    subplot(3,length(JFIDs)+1,JND)
    contourf(XData,YData,Hybrid{1}.All.(JFIDs{JND}).VX_R2)
    xlabel(LabX)
    clim(setclim)
    title(JFIDs{JND})
    if JND == 1
        ylabel({'Velocity X',LabY})
    end

    subplot(3,length(JFIDs)+1,JND+(length(JFIDs)+1))
    contourf(XData,YData,Hybrid{1}.All.(JFIDs{JND}).VY_R2)
    xlabel(LabX)
    clim(setclim)
    if JND == 1
        ylabel({'Velocity Y',LabY})
    end

    subplot(3,length(JFIDs)+1,JND+2*(length(JFIDs)+1))
    contourf(XData,YData,Hybrid{1}.All.(JFIDs{JND}).VZ_R2)
    xlabel(LabX)
    clim(setclim)
    if JND == 1
        ylabel({'Velocity Z',LabY})
    end
end
subplot(3,length(JFIDs)+1,JND+1)
contourf(XData,YData,Hybrid{1}.All.VX_R2)
xlabel(LabX)
clim(setclim)
title('All')
subplot(3,length(JFIDs)+1,JND+1+(length(JFIDs)+1))
contourf(XData,YData,Hybrid{1}.All.VY_R2)
xlabel(LabX)
clim(setclim)
subplot(3,length(JFIDs)+1,JND+1+2*(length(JFIDs)+1))
contourf(XData,YData,Hybrid{1}.All.VZ_R2)
xlabel(LabX)
clim(setclim)
sgtitle({LabAll,'R2 Trained by JF'})
colormap("jet")
set(gcf,"Position",[135         133        1708         804])



%%%%%%%%%%%%%% Plot Full Time Position Predictions
setclim = [0,1];
figure
for JND = 1:length(Hybrid)
    subplot(3,length(Hybrid),JND)
    contourf(XData,YData,Hybrid{JND}.All.X_R2)
    xlabel(LabX)
    clim(setclim)
    title(ConfigNames{JND})
    if JND == 1
        ylabel({'X',LabY})
    end

    subplot(3,length(Hybrid),JND+length(Hybrid))
    contourf(XData,YData,Hybrid{JND}.All.Y_R2)
    xlabel(LabX)
    clim(setclim)
    if JND == 1
        ylabel({'Y',LabY})
    end

    subplot(3,length(Hybrid),JND+2*length(Hybrid))
    contourf(XData,YData,Hybrid{JND}.All.Z_R2)
    xlabel('Dead-Reckoned Future Time')
    clim(setclim)
    if JND == 1
        ylabel({'Z',LabY})
    end
end
colormap("jet")
sgtitle({LabAll,'Dead-Reckoned Future Position Estimation'})
set(gcf,"Position",[135         133        1708         804])


%%%%%%%%%%%%%%% Plot Full Time Position Predictions
setclim = [0,1];
figure
for JND = 1:length(JFIDs)
    subplot(3,length(JFIDs),JND)
    contourf(XData,YData,Hybrid{1}.All.(JFIDs{JND}).X_R2)
    xlabel(LabX)
    clim(setclim)
    title(JFIDs{JND})
    if JND == 1
        ylabel({'X',LabY})
    end

    subplot(3,length(JFIDs),JND+length(JFIDs))
    contourf(XData,YData,Hybrid{1}.All.(JFIDs{JND}).Y_R2)
    xlabel(LabX)
    clim(setclim)
    if JND == 1
        ylabel({'Y',LabY})
    end

    subplot(3,length(JFIDs),JND+2*length(JFIDs))
    contourf(XData,YData,Hybrid{1}.All.(JFIDs{JND}).Z_R2)
    xlabel('Dead-Reckoned Future Time')
    clim(setclim)
    if JND == 1
        ylabel({'Z',LabY})
    end
end
colormap("jet")
sgtitle({LabAll,'JF Dead-Reckoned Future Position Estimation'})
set(gcf,"Position",[135         133        1708         804])


%%%%%%%%%%%%%%% 
Ptarg = 'INPUT_3'; % 'INPUT_3' 'JF41'
setclim = [0,1];
figure
tiledlayout(1,9,"TileSpacing","compact")
nexttile
contourf(XData,YData,Hybrid{1}.(Ptarg).P_COM_X_R2)
xlabel(LabX)
ylabel({'Hybrid',LabY})
title('Position X')
clim(setclim)
nexttile
contourf(XData,YData,Hybrid{1}.(Ptarg).P_COM_Y_R2)
xlabel(LabX)
title('Position Y')
clim(setclim)
nexttile
contourf(XData,YData,Hybrid{1}.(Ptarg).P_COM_Z_R2)
xlabel(LabX)
title('Position Z')
clim(setclim)
nexttile
contourf(XData,YData,Hybrid{1}.(Ptarg).P_COM_R1_R2)
xlabel(LabX)
title('Rotation 1')
clim(setclim)
nexttile
contourf(XData,YData,Hybrid{1}.(Ptarg).P_COM_R2_R2)
xlabel(LabX)
title('Rotation 2')
clim(setclim)
nexttile
contourf(XData,YData,Hybrid{1}.(Ptarg).P_COM_R3_R2)
xlabel(LabX)
title('Rotation 3')
clim(setclim)
nexttile
contourf(XData,YData,Hybrid{1}.(Ptarg).P_BODY_VX_R2)
xlabel(LabX)
title('Velocity X')
clim(setclim)
nexttile
contourf(XData,YData,Hybrid{1}.(Ptarg).P_BODY_VY_R2)
xlabel(LabX)
title('Velocity Y')
clim(setclim)
nexttile
contourf(XData,YData,Hybrid{1}.(Ptarg).P_BODY_VZ_R2)
xlabel(LabX)
clim(setclim)
title('Velocity Z')
colormap("jet")
colorbar()
set(gcf,"Position",[76         442        1790         228])



setclim = [0,1];
figure
tiledlayout(2,9,"TileSpacing","compact")
nexttile
contourf(XData,YData,Sequential{1}.(Ptarg).P_COM_X_R2)
ylabel({'ESN Body',LabY})
title('Position X')
clim(setclim)
nexttile
contourf(XData,YData,Sequential{1}.(Ptarg).P_COM_Y_R2)
title('Position Y')
clim(setclim)
nexttile
contourf(XData,YData,Sequential{1}.(Ptarg).P_COM_Z_R2)
title('Position Z')
clim(setclim)
nexttile
contourf(XData,YData,Sequential{1}.(Ptarg).P_COM_R1_R2)
title('Rotation 1')
clim(setclim)
nexttile
contourf(XData,YData,Sequential{1}.(Ptarg).P_COM_R2_R2)
title('Rotation 2')
clim(setclim)
nexttile
contourf(XData,YData,Sequential{1}.(Ptarg).P_COM_R3_R2)
title('Rotation 3')
clim(setclim)
nexttile
contourf(XData,YData,Sequential{1}.(Ptarg).P_BODY_VX_R2)
title('Velocity X')
clim(setclim)
nexttile
contourf(XData,YData,Sequential{1}.(Ptarg).P_BODY_VY_R2)
title('Velocity Y')
clim(setclim)
nexttile
contourf(XData,YData,Sequential{1}.(Ptarg).P_BODY_VZ_R2)
clim(setclim)
title('Velocity Z')
colormap("jet")
colorbar()

nexttile
contourf(XData,YData,PRC{1}.(Ptarg).P_COM_X_R2)
xlabel(LabX)
ylabel({'PRC',LabY})
clim(setclim)
nexttile
contourf(XData,YData,PRC{1}.(Ptarg).P_COM_Y_R2)
xlabel(LabX)
clim(setclim)
nexttile
contourf(XData,YData,PRC{1}.(Ptarg).P_COM_Z_R2)
xlabel(LabX)
clim(setclim)
nexttile
contourf(XData,YData,PRC{1}.(Ptarg).P_COM_R1_R2)
xlabel(LabX)
clim(setclim)
nexttile
contourf(XData,YData,PRC{1}.(Ptarg).P_COM_R2_R2)
xlabel(LabX)
clim(setclim)
nexttile
contourf(XData,YData,PRC{1}.(Ptarg).P_COM_R3_R2)
xlabel(LabX)
clim(setclim)
nexttile
contourf(XData,YData,PRC{1}.(Ptarg).P_BODY_VX_R2)
xlabel(LabX)
clim(setclim)
nexttile
contourf(XData,YData,PRC{1}.(Ptarg).P_BODY_VY_R2)
xlabel(LabX)
clim(setclim)
nexttile
contourf(XData,YData,PRC{1}.(Ptarg).P_BODY_VZ_R2)
xlabel(LabX)
clim(setclim)
colormap("jet")
colorbar()
set(gcf,"Position",[22         296        1872         381])


end

function [] = Plot_R2_NoStim_Results(ESN,PRC,Sequential,Hybrid,RandESN,ConfigNames,N)
load('MyColors.Mat', 'FavoriteColors');

Real2SampleTime = 2/120;
XData = (0:120)*Real2SampleTime;
YData = (0:120)*Real2SampleTime;

LabAll = strcat('All Data N',num2str(N));
LabPulse = strcat('Pulsatile N',num2str(N));
LabX = 'Prediction Future Horizon';
LabY = 'Leaky-Int. History Horizon';
setclim = [-0.07,0.07];

 UnStimAllFLDS = fields(Hybrid{1}.All.NoStim);
 UnStimPulseFlds = fields(Hybrid{1}.Unstimulated); 

setclim = [-0.2,1];
TMPclin = linspace(setclim(1),setclim(2),256);
PositveLen =  length(TMPclin(TMPclin>=0));
NegativeLen =  length(TMPclin(TMPclin<0));
CmapPositive =  jet(PositveLen);
CmapNegative = [linspace(0,0.5,NegativeLen)',linspace(0,0.5,NegativeLen)',linspace(0,0.5,NegativeLen)'];
Cmap = [CmapNegative;CmapPositive];


figure
tiledlayout(5,length(UnStimAllFLDS),"TileSpacing","compact")
for JND = 1:length(UnStimAllFLDS)
    nexttile(JND)
    contourf(XData,YData,RandESN.All.NoStim.(UnStimAllFLDS{JND}))
    clim(setclim)
    title(UnStimAllFLDS{JND})
    if JND == 1
        ylabel({'ESN Random',LabY})
    elseif JND == length(UnStimAllFLDS)
        colorbar()
    end

    nexttile(JND+length(UnStimAllFLDS))
    contourf(XData,YData,ESN.All.NoStim.(UnStimAllFLDS{JND}))
    clim(setclim)
    if JND == 1
        ylabel({'ESN Input',LabY})
    elseif JND == length(UnStimAllFLDS)
        colorbar()
    end

    nexttile(JND+2*length(UnStimAllFLDS))
    contourf(XData,YData,Sequential{1}.All.NoStim.(UnStimAllFLDS{JND}))
    clim(setclim)
    if JND == 1
        ylabel({'ESN Body',LabY})
    elseif JND == length(UnStimAllFLDS)
        colorbar()
    end


    nexttile(JND+3*length(UnStimAllFLDS))
    contourf(XData,YData,PRC{1}.All.NoStim.(UnStimAllFLDS{JND}))
    clim(setclim)
    if JND == 1
        ylabel({'PRC',LabY})
    elseif JND == length(UnStimAllFLDS)
        colorbar()
    end


    nexttile(JND+4*length(UnStimAllFLDS))
    contourf(XData,YData,Hybrid{1}.All.NoStim.(UnStimAllFLDS{JND}))
    xlabel(LabX)
    clim(setclim)
    if JND == 1
        ylabel({'Hybrid',LabY})
    elseif JND == length(UnStimAllFLDS)
        colorbar()
    end
end
sgtitle('All: Unstimulated Data')
colormap(Cmap)
set(gcf,"Position",[331          99        1313         844])



figure
tiledlayout(5,length(UnStimPulseFlds),"TileSpacing","compact")
for JND = 1:length(UnStimPulseFlds)
    nexttile(JND)
    contourf(XData,YData,RandESN.Unstimulated.(UnStimPulseFlds{JND}))
    clim(setclim)
    title(UnStimPulseFlds{JND})
    if JND == 1
        ylabel({'ESN Random',LabY})
    elseif JND == length(UnStimPulseFlds)
        colorbar()
    end

    nexttile(JND+length(UnStimPulseFlds))
    contourf(XData,YData,ESN.Unstimulated.(UnStimPulseFlds{JND}))
    clim(setclim)
    if JND == 1
        ylabel({'ESN Input',LabY})
    elseif JND == length(UnStimPulseFlds)
        colorbar()
    end

    nexttile(JND+2*length(UnStimPulseFlds))
    contourf(XData,YData,Sequential{1}.Unstimulated.(UnStimPulseFlds{JND}))
    clim(setclim)
    if JND == 1
        ylabel({'ESN Body',LabY})
    elseif JND == length(UnStimPulseFlds)
        colorbar()
    end


    nexttile(JND+3*length(UnStimPulseFlds))
    contourf(XData,YData,PRC{1}.Unstimulated.(UnStimPulseFlds{JND}))
    clim(setclim)
    if JND == 1
        ylabel({'PRC',LabY})
    elseif JND == length(UnStimPulseFlds)
        colorbar()
    end


    nexttile(JND+4*length(UnStimPulseFlds))
    contourf(XData,YData,Hybrid{1}.Unstimulated.(UnStimPulseFlds{JND}))
    xlabel(LabX)
    clim(setclim)
    if JND == 1
        ylabel({'Hybrid',LabY})
    elseif JND == length(UnStimPulseFlds)
        colorbar()
    end
end
sgtitle('Only Spontaneous Pulses: Unstimulated Data')
colormap(Cmap)
set(gcf,"Position",[76          99        1809         844])



end

function [] = Plot_State_Data(BodyData,TargetData,SenseConfig,N,Spect,Uscale,Tback,VideoName)
%%% Just Compute the Pulsatile Data from the ESN 

MPLX_Len = 121;  %%% 61
Dt = 0.0167;

TWo_Pulse = 1000; 
TWo_ALL = 10000; 
T_Train_Frac = 1/1; %%% Fraction of data used for training


%%%%%% Full time Targets
% [Targ_VX,Targ_X] = MultiplexData(TargetData.All.All.VX,MPLX_Len);
% [Targ_VY,Targ_Y] = MultiplexData(TargetData.All.All.VY,MPLX_Len);
% [Targ_VZ,Targ_Z] = MultiplexData(TargetData.All.All.VZ,MPLX_Len);
% Targ_X = Targ_X*Dt;
% Targ_Y = Targ_Y*Dt;
% Targ_Z = Targ_Z*Dt;

%%%% Check both for Pulsatile targets and for Whole time targets
PulsatileFlds = {'JF41','JF42','JF43','JF44','JF46','INPUT_1','INPUT_2','INPUT_3','INPUT_4'};
JFFlds = {'JF41','JF42','JF43','JF44','JF46'};
PulseTargets = TargetData.Pulse_Targets;
NumSensors = length(SenseConfig);

 FlippedIND = fliplr(1:MPLX_Len);

B_Fld1 = 'JF41';
B_Fld2 = 'P';

%%%% Run through ALL time Data for Sensors
Ut = 2*flipud(MultiplexData(BodyData.(B_Fld1).(B_Fld2).Input,MPLX_Len))-1;
BodySenseData = zeros(MPLX_Len*NumSensors,size(Ut,2));
for LND = 1:NumSensors
    %%%% Normalize sensor data [-1,1]
    TMPSens = BodyData.(B_Fld1).(B_Fld2).(SenseConfig{LND});
    TMPSens = 2*(TMPSens-range(TMPSens)/2-min(TMPSens))/range(TMPSens);
    BodySenseData(Tback*(LND-1)+(1:Tback),:) = flipud(MultiplexData(TMPSens,Tback));
end

T_Flds = {'BODY_VX','BODY_VY','BODY_VZ'};
for KND = 1:length(T_Flds)
    CurTargets.(T_Flds{KND}) = MultiplexData(TargetData.(B_Fld1).(B_Fld2).(T_Flds{KND}),MPLX_Len);
    [CurTargets.(T_Flds{KND}),PureESN.(T_Flds{KND}),PRC.(T_Flds{KND}),SequentialESN.(T_Flds{KND}),HybridESN.(T_Flds{KND})] = ComputeStateData(Ut,BodySenseData, CurTargets.(T_Flds{KND}),Tback,N,TWo_Pulse,Spect,Uscale);
end

%%%%% Basic Figure
figure
for KND = 1:length(T_Flds)
    subplot(length(T_Flds),1,KND)
    plot(CurTargets.(T_Flds{KND})(:,1),'linewidth',2,'Color','k')
    hold on
    plot(PureESN.(T_Flds{KND}).Approx(:,1),'linewidth',2,'Color','b')
    plot(PRC.(T_Flds{KND}).Approx(:,1),'linewidth',2,'Color','r')
    plot(SequentialESN.(T_Flds{KND}).Approx(:,1),'linewidth',2,'Color','c')
    plot(HybridESN.(T_Flds{KND}).Approx(:,1),'linewidth',2,'Color','m')
    xlim([20000,22000])
    if KND == 1
        legend({'Target','Pure ESN','PRC','Sequential','Hybrid'},"Orientation","horizontal")
    end
    ylabel(T_Flds{KND})
end
sgtitle('JF41 Pulsatile f(t_0-t_2) = V(t_0)')
set(gcf,"Position",[680   177   717   801])


%%%%% Basic Figure (last frame)
figure
for KND = 1:length(T_Flds)
    subplot(length(T_Flds),1,KND)
    plot(CurTargets.(T_Flds{KND})(:,end),'linewidth',2,'Color','k')
    hold on
    plot(PureESN.(T_Flds{KND}).Approx(:,end),'linewidth',2,'Color','b')
    plot(PRC.(T_Flds{KND}).Approx(:,end),'linewidth',2,'Color','r')
    plot(SequentialESN.(T_Flds{KND}).Approx(:,end),'linewidth',2,'Color','c')
    plot(HybridESN.(T_Flds{KND}).Approx(:,end),'linewidth',2,'Color','m')
    xlim([20000,22000])
    if KND == 1
        legend({'Target','Pure ESN','PRC','Sequential','Hybrid'},"Orientation","horizontal")
    end
    ylabel(T_Flds{KND})
end
sgtitle('JF41 Pulsatile f(t_0-t_2) = V(t_2)')
set(gcf,"Position",[680   177   717   801])




%%%%% Video Here with projections and movement
VideoSampleLength = 20000:2:26000;
State_History_Length = 60*6;

if  ~isempty(VideoName)  == 1
    v = VideoWriter(VideoName,'MPEG-4');
    v.Quality = 95;
    open(v)

    figure
    for IND = VideoSampleLength
        for KND = 1:length(T_Flds)
            subplot(length(T_Flds),1,KND)
            plot(IND+(-State_History_Length:0),CurTargets.(T_Flds{KND})(IND+(-State_History_Length:0),1),'linewidth',1,'Color','k')
            hold on
            plot(IND+(-120:0),CurTargets.(T_Flds{KND})(IND-120,:),'linewidth',2,'Color',[212,175,55]/255)
            plot(IND+(0:120),CurTargets.(T_Flds{KND})(IND,:),'linewidth',1,'Color',[0.5,0.5,0.5,0.5])
            plot(IND+(0:120),PureESN.(T_Flds{KND}).Approx(IND-120,:),'linewidth',2,'Color','b')
            plot(IND+(0:120),PRC.(T_Flds{KND}).Approx(IND,:),'linewidth',2,'Color','r')
            plot(IND+(0:120),SequentialESN.(T_Flds{KND}).Approx(IND,:),'linewidth',2,'Color','c')
            plot(IND+(0:120),HybridESN.(T_Flds{KND}).Approx(IND,:),'linewidth',2,'Color','m')
            xlim(IND+[-State_History_Length,+1.5*120])
            xline(IND,':')
            if KND == 2
                legend({'Original','ESN Horizon','Actual Future','Pure ESN Prediction','PRC Prediction','Sequential Prediction','Hybrid Prediction'},...
                    "Orientation","horizontal",'NumColumns',2,'Location','southwest')
            end
            ylabel(T_Flds{KND})
            hold off
        end
        set(gcf,"Position",[680   177   717   801])
        drawnow
        A = getframe(gcf);
        writeVideo(v,A)
    end

    close(v)
end


end

function [] = Plot_State_Data2(BodyData,TargetData,SenseConfig,N,Spect,Uscale,Tback)
%%% Just Compute the Pulsatile Data from the ESN 
load('MyColors.Mat', 'FavoriteColors');

MPLX_Len = 121;  %%% 61
Dt = 0.0167;
TWo_Pulse = 1000; 
TWo_ALL = 10000; 
T_Train_Frac = 1/1; %%% Fraction of data used for training

%%%% Check both for Pulsatile targets and for Whole time targets
PulsatileFlds = {'JF41','JF42','JF43','JF44','JF46','INPUT_1','INPUT_2','INPUT_3','INPUT_4','Unstimulated'};
JFFlds = {'JF41','JF42','JF43','JF44','JF46'};
PulseTargets = TargetData.Pulse_Targets;
NumSensors = length(SenseConfig);

 FlippedIND = fliplr(1:MPLX_Len);

B_Fld1 = {'All','All','JF41'};
B_Fld2 = {'NoStim_JF41','JF41','P'};
T_FldsA = {'VX','VY','VZ'};
T_FldsP = {'BODY_VX','BODY_VY','BODY_VZ'};

%%%% Run through ALL time Data for Sensors
Ut  = cell(size(B_Fld2));
BodySenseData = cell(size(B_Fld2));
for INC = 1:length(B_Fld1)
    Ut{INC} = 2*flipud(MultiplexData(BodyData.(B_Fld1{INC}).(B_Fld2{INC}).Input,MPLX_Len))-1;
    BodySenseData{INC} = zeros(MPLX_Len*NumSensors,size(Ut{INC},2));
    for LND = 1:NumSensors
        %%%% Normalize sensor data [-1,1]
        TMPSens = BodyData.(B_Fld1{INC}).(B_Fld2{INC}).(SenseConfig{LND});
        TMPSens = 2*(TMPSens-range(TMPSens)/2-min(TMPSens))/range(TMPSens);
        BodySenseData{INC}(Tback*(LND-1)+(1:Tback),:) = flipud(MultiplexData(TMPSens,Tback));
    end
    Ut{INC} = Ut{INC}/Tback;
end

CurTargets  = cell(size(B_Fld2));
PureESN = cell(size(B_Fld2));
SequentialESN  = cell(size(B_Fld2));
PRC = cell(size(B_Fld2));
HybridESN  = cell(size(B_Fld2));

for INC = 1:length(B_Fld1)
    if strcmp(B_Fld1{INC},'All')
        T_Flds=  T_FldsA;
    else
        T_Flds =   T_FldsP;
    end

    for KND = 1:length(T_Flds)
        if  strcmp(B_Fld1{INC},'All')
             CurTargets{INC}.(T_Flds{KND}) = MultiplexData(smoothdata(TargetData.(B_Fld1{INC}).(B_Fld2{INC}).(T_Flds{KND}),"movmean",4),MPLX_Len);
        else
            CurTargets{INC}.(T_Flds{KND}) = MultiplexData(TargetData.(B_Fld1{INC}).(B_Fld2{INC}).(T_Flds{KND}),MPLX_Len);
        end
        [CurTargets{INC}.(T_Flds{KND}),PureESN{INC}.(T_Flds{KND}),PRC{INC}.(T_Flds{KND}),SequentialESN{INC}.(T_Flds{KND}),HybridESN{INC}.(T_Flds{KND})] =...
            ComputeStateData(Ut{INC},BodySenseData{INC}, CurTargets{INC}.(T_Flds{KND}),Tback,N,TWo_Pulse,Spect,Uscale);
    end
end

%%%%%%%% The confusion  test
HybridESN_Confuse  = cell(length(B_Fld2)*[1,1]);
for INC = 1:length(B_Fld1)
    if strcmp(B_Fld1{INC},'All')
        T_Flds_a=  T_FldsA;
    else
        T_Flds_a =   T_FldsP;
    end

    for JNC = 1:length(B_Fld1)
        if strcmp(B_Fld1{JNC},'All')
            T_Flds_b =  T_FldsA;
        else
            T_Flds_b =   T_FldsP;
        end

        %%%% Only look at x-x, y-y, and z-z
        for KND = 1:length(T_Flds_b)
           Ctarg =  CurTargets{INC}.(T_Flds_a{KND});
           CTrained = HybridESN{INC}.(T_Flds_a{KND});
           CDiffTrained = HybridESN{JNC}.(T_Flds_b{KND});

           TMPconf.W = CDiffTrained.W;
           TMPconf.States = CTrained.States;
           TMPconf.Approx = CTrained.States*CDiffTrained.W;
           TMPconf.R2 = (1-((sum((TMPconf.Approx-Ctarg).^2,1))./(sum((Ctarg - mean(Ctarg)).^2,1))));
           
           HybridESN_Confuse{INC,JNC}.(T_FldsA{KND})  = TMPconf;
        end
    end
end


B_Fld3 = {'All','Unstimulated','INPUT_1','INPUT_2','INPUT_3','INPUT_4'};
B_Fld4 = {'NoStim','','P','P','P','P'};
%%%% Run through ALL time Data for Sensors
Ut  = cell(size(B_Fld3));
BodySenseData = cell(size(B_Fld3));
for INC = 1:length(B_Fld3)
    if strcmp(B_Fld3{INC},'Unstimulated')
        Ut{INC} = 2*flipud(MultiplexData(BodyData.(B_Fld3{INC}).Input,MPLX_Len))-1;
        CBD = BodyData.(B_Fld3{INC});
    else
        Ut{INC} = 2*flipud(MultiplexData(BodyData.(B_Fld3{INC}).(B_Fld4{INC}).Input,MPLX_Len))-1;
        CBD = BodyData.(B_Fld3{INC}).(B_Fld4{INC});
    end
    BodySenseData{INC} = zeros(MPLX_Len*NumSensors,size(Ut{INC},2));
    for LND = 1:NumSensors
        %%%% Normalize sensor data [-1,1]
        TMPSens = CBD.(SenseConfig{LND});
        TMPSens = 2*(TMPSens-range(TMPSens)/2-min(TMPSens))/range(TMPSens);
        BodySenseData{INC}(Tback*(LND-1)+(1:Tback),:) = flipud(MultiplexData(TMPSens,Tback));
    end
    Ut{INC} = Ut{INC}/Tback;
end

CurTargets_In  = cell(size(B_Fld4));
PureESN_In = cell(size(B_Fld4));
SequentialESN_In  = cell(size(B_Fld4));
PRC_In = cell(size(B_Fld4));
HybridESN_In  = cell(size(B_Fld4));

for INC = 1:length(B_Fld3)
    if strcmp(B_Fld3{INC},'All')
        T_Flds=  T_FldsA;
    else
        T_Flds =   T_FldsP;
    end

    for KND = 1:length(T_Flds)
        if  strcmp(B_Fld3{INC},'All')
             CurTargets_In{INC}.(T_Flds{KND}) = MultiplexData(smoothdata(TargetData.(B_Fld3{INC}).(B_Fld4{INC}).(T_Flds{KND}),"movmean",4),MPLX_Len);
        elseif strcmp(B_Fld3{INC},'Unstimulated')
            CurTargets_In{INC}.(T_Flds{KND}) = MultiplexData(smoothdata(TargetData.(B_Fld3{INC}).(T_Flds{KND}),"movmean",4),MPLX_Len);
        else
            CurTargets_In{INC}.(T_Flds{KND}) = MultiplexData(TargetData.(B_Fld3{INC}).(B_Fld4{INC}).(T_Flds{KND}),MPLX_Len);
        end

        [CurTargets_In{INC}.(T_Flds{KND}),PureESN_In{INC}.(T_Flds{KND}),PRC_In{INC}.(T_Flds{KND}),SequentialESN_In{INC}.(T_Flds{KND}),HybridESN_In{INC}.(T_Flds{KND})] =...
            ComputeStateData(Ut{INC},BodySenseData{INC}, CurTargets_In{INC}.(T_Flds{KND}),Tback,N,TWo_Pulse,Spect,Uscale);
    end
end

%%%%%%%% The confusion  test
HybridESN_Confuse_In  = cell(length(B_Fld3)*[1,1]);
for INC = 1:length(B_Fld3)
    if strcmp(B_Fld3{INC},'All')
        T_Flds_a=  T_FldsA;
    else
        T_Flds_a =   T_FldsP;
    end

    for JNC = 1:length(B_Fld3)
        if strcmp(B_Fld3{JNC},'All')
            T_Flds_b =  T_FldsA;
        else
            T_Flds_b =   T_FldsP;
        end

        %%%% Only look at x-x, y-y, and z-z
        for KND = 1:length(T_Flds_b)
           Ctarg =  CurTargets_In{INC}.(T_Flds_a{KND});
           CTrained = HybridESN_In{INC}.(T_Flds_a{KND});
           CDiffTrained = HybridESN_In{JNC}.(T_Flds_b{KND});

           TMPconf.W = CDiffTrained.W;
           TMPconf.States = CTrained.States;
           TMPconf.Approx = CTrained.States*CDiffTrained.W;
           TMPconf.R2 = (1-((sum((TMPconf.Approx-Ctarg).^2,1))./(sum((Ctarg - mean(Ctarg)).^2,1))));
           
           HybridESN_Confuse_In{INC,JNC}.(T_FldsA{KND})  = TMPconf;
        end
    end
end


%%%%%%%%%%%%%%%%%%%%% Figures %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% Basic Figure: JF41
figure
for KND = 1:length(T_FldsP)
    subplot(length(T_FldsP),1,KND)
    plot(CurTargets{3}.(T_FldsP{KND})(:,1),'linewidth',2,'Color','k')
    hold on
    plot(PureESN{3}.(T_FldsP{KND}).Approx(:,1),'linewidth',2,'Color','b')
    plot(PRC{3}.(T_FldsP{KND}).Approx(:,1),'linewidth',2,'Color','r')
    plot(SequentialESN{3}.(T_FldsP{KND}).Approx(:,1),'linewidth',2,'Color','c')
    plot(HybridESN{3}.(T_FldsP{KND}).Approx(:,1),'linewidth',2,'Color','m')
    xlim([20000,22000])
    if KND == 1
        legend({'Target','Pure ESN','PRC','Sequential','Hybrid'},"Orientation","horizontal")
    end
    ylabel(T_FldsP{KND})
end
sgtitle('JF41 Pulsatile f(t_0-t_2) = V(t_0)')
set(gcf,"Position",[680   177   717   801])


%%%%% Basic Figure (last frame): JF41
figure
for KND = 1:length(T_FldsP)
    subplot(length(T_FldsP),1,KND)
    plot(CurTargets{3}.(T_FldsP{KND})(:,end),'linewidth',2,'Color','k')
    hold on
    plot(PureESN{3}.(T_FldsP{KND}).Approx(:,end),'linewidth',2,'Color','b')
    plot(PRC{3}.(T_FldsP{KND}).Approx(:,end),'linewidth',2,'Color','r')
    plot(SequentialESN{3}.(T_FldsP{KND}).Approx(:,end),'linewidth',2,'Color','c')
    plot(HybridESN{3}.(T_FldsP{KND}).Approx(:,end),'linewidth',2,'Color','m')
    xlim([20000,22000])
    if KND == 1
        legend({'Target','Pure ESN','PRC','Sequential','Hybrid'},"Orientation","horizontal")
    end
    ylabel(T_Flds{KND})
end
sgtitle('JF41 Pulsatile f(t_0-t_2) = V(t_2)')
set(gcf,"Position",[680   177   717   801])

%%%%% Basic Figure (0.2s): JF41
figure
for KND = 1:length(T_FldsP)
    subplot(length(T_FldsP),1,KND)
    plot(CurTargets{3}.(T_FldsP{KND})(:,12),'linewidth',2,'Color','k')
    hold on
    plot(PureESN{3}.(T_FldsP{KND}).Approx(:,12),'linewidth',2,'Color','b')
    plot(PRC{3}.(T_FldsP{KND}).Approx(:,12),'linewidth',2,'Color','r')
    plot(SequentialESN{3}.(T_FldsP{KND}).Approx(:,12),'linewidth',2,'Color','c')
    plot(HybridESN{3}.(T_FldsP{KND}).Approx(:,12),'linewidth',2,'Color','m')
    xlim([20000,22000])
    if KND == 1
        legend({'Target','Pure ESN','PRC','Sequential','Hybrid'},"Orientation","horizontal")
    end
    ylabel(T_Flds{KND})
end
sgtitle('JF41 Pulsatile f(t_0-t_{0.2}) = V(t_2)')
set(gcf,"Position",[680   177   717   801])


%%%%% Basic Figure: Tau = 1.5s
figure
for KND = 1:length(T_FldsP)
    subplot(length(T_FldsP),1,KND)
    plot(CurTargets_In{5}.(T_FldsP{KND})(:,1),'linewidth',2,'Color','k')
    hold on
    plot(PureESN_In{5}.(T_FldsP{KND}).Approx(:,1),'linewidth',2,'Color','b')
    plot(PRC_In{5}.(T_FldsP{KND}).Approx(:,1),'linewidth',2,'Color','r')
    plot(SequentialESN_In{5}.(T_FldsP{KND}).Approx(:,1),'linewidth',2,'Color','c')
    plot(HybridESN_In{5}.(T_FldsP{KND}).Approx(:,1),'linewidth',2,'Color','m')
    xlim([20000,22000])
    if KND == 1
        legend({'Target','Pure ESN','PRC','Sequential','Hybrid'},"Orientation","horizontal")
    end
    ylabel(T_FldsP{KND})
end
sgtitle('\tau=1.5s  Pulsatile f(t_0-t_2) = V(t_0)')
set(gcf,"Position",[680   177   717   801])


%%%%% Basic Figure (0.2s): Tau = 1.5s
figure
for KND = 1:length(T_FldsP)
    subplot(length(T_FldsP),1,KND)
    plot(CurTargets_In{5}.(T_FldsP{KND})(:,12),'linewidth',2,'Color','k')
    hold on
    plot(PureESN_In{5}.(T_FldsP{KND}).Approx(:,12),'linewidth',2,'Color','b')
    plot(PRC_In{5}.(T_FldsP{KND}).Approx(:,12),'linewidth',2,'Color','r')
    plot(SequentialESN_In{5}.(T_FldsP{KND}).Approx(:,12),'linewidth',2,'Color','c')
    plot(HybridESN_In{5}.(T_FldsP{KND}).Approx(:,12),'linewidth',2,'Color','m')
    xlim([20000,22000])
    if KND == 1
        legend({'Target','Pure ESN','PRC','Sequential','Hybrid'},"Orientation","horizontal")
    end
    ylabel(T_FldsP{KND})
end
sgtitle('\tau=1.5s  Pulsatile f(t_0-t_{0.2}) = V(t_2)')
set(gcf,"Position",[680   177   717   801])


%%%%% Basic Figure (last frame): Tau = 1.5s
figure
for KND = 1:length(T_FldsP)
    subplot(length(T_FldsP),1,KND)
    plot(CurTargets_In{5}.(T_FldsP{KND})(:,end),'linewidth',2,'Color','k')
    hold on
    plot(PureESN_In{5}.(T_FldsP{KND}).Approx(:,end),'linewidth',2,'Color','b')
    plot(PRC_In{5}.(T_FldsP{KND}).Approx(:,end),'linewidth',2,'Color','r')
    plot(SequentialESN_In{5}.(T_FldsP{KND}).Approx(:,end),'linewidth',2,'Color','c')
    plot(HybridESN_In{5}.(T_FldsP{KND}).Approx(:,end),'linewidth',2,'Color','m')
    xlim([20000,22000])
    if KND == 1
        legend({'Target','Pure ESN','PRC','Sequential','Hybrid'},"Orientation","horizontal")
    end
    ylabel(T_FldsP{KND})
end
sgtitle('\tau=1.5s  Pulsatile f(t_0-t_2) = V(t_2)')
set(gcf,"Position",[680   177   717   801])




StartPts = [23460;102655;113700;10500];  %%[49200, 54600,106700 ,2000]
% StartPts = [49290;103100;113700;10500];  %%[49200, 54600,106700 ,2000]
WindowSize = 2000-200;
PredictionSep =120; %120 
Timeconversion = 1/60;
PredictionsSt = StartPts + (0:PredictionSep:WindowSize);
% TMcolor = FavoriteColors.Map.Sea(6);
% PredictionColor = TMcolor(2,:);
PredictionColor = FavoriteColors.PptFrm.Green;
PredStColor = FavoriteColors.PptFrm.DarkPink;
Mkr = 's'; % '|' 
Mkrsz = 7;

%%%%%% Ploting data from different periods
figure
tiledlayout(4,1,TileSpacing="tight",Padding="tight")
%%%% One plot without stimulus
nexttile
plot(Timeconversion*(1:length(CurTargets{1}.VZ(:,1))),CurTargets{1}.VZ(:,1),'k',LineWidth=2);
hold on
for JNC = 1:length(PredictionsSt(1,:))
    plot(Timeconversion*PredictionsSt(1,JNC),HybridESN{1}.VZ.Approx(PredictionsSt(1,JNC),1),Mkr,'Color',PredStColor,'MarkerFaceColor',PredStColor,'MarkerSize',Mkrsz,LineWidth=2);
    plot(Timeconversion*(PredictionsSt(1,JNC)+(0:120)),HybridESN{1}.VZ.Approx(PredictionsSt(1,JNC),:),'Color',PredictionColor,LineWidth=2);
end
xlim(Timeconversion*(StartPts(1)+[0,WindowSize]))
ylabel('Velocity Z')

%%%% One with Jellyfish but no stimulus and sponaneous pulse
nexttile
plot(Timeconversion*(1:length(CurTargets{2}.VZ(:,1))),CurTargets{2}.VZ(:,1),'k',LineWidth=2);
hold on
for JNC = 1:length(PredictionsSt(1,:))
    plot(Timeconversion*PredictionsSt(2,JNC),HybridESN{2}.VZ.Approx(PredictionsSt(2,JNC),1),Mkr,'Color',PredStColor,'MarkerFaceColor',PredStColor,'MarkerSize',Mkrsz,LineWidth=2);
    plot(Timeconversion*(PredictionsSt(2,JNC)+(0:120)),HybridESN{2}.VZ.Approx(PredictionsSt(2,JNC),:),'Color',PredictionColor,LineWidth=2);
end
xlim(Timeconversion*(StartPts(2)+[0,WindowSize]))
ylabel('Velocity Z')

%%%% One with Jellyfish with stimulus
nexttile
plot(Timeconversion*(1:length(CurTargets{2}.VZ(:,1))),CurTargets{2}.VZ(:,1),'k',LineWidth=2);
hold on
for JNC = 1:length(PredictionsSt(1,:))
    plot(Timeconversion*PredictionsSt(3,JNC),HybridESN{2}.VZ.Approx(PredictionsSt(3,JNC),1),Mkr,'Color',PredStColor,'MarkerFaceColor',PredStColor,'MarkerSize',Mkrsz,LineWidth=2);
    plot(Timeconversion*(PredictionsSt(3,JNC)+(0:120)),HybridESN{2}.VZ.Approx(PredictionsSt(3,JNC),:),'Color',PredictionColor,LineWidth=2);
end
xlim(Timeconversion*(StartPts(3)+[0,WindowSize]))
ylabel('Velocity Z')
ylim([-20,20])

%%%% One with stimulated Pulse
nexttile
plot(Timeconversion*(1:length(CurTargets{3}.BODY_VZ(:,1))),CurTargets{3}.BODY_VZ(:,1),'k',LineWidth=2);
hold on
for JNC = 1:length(PredictionsSt(1,:))
    plot(Timeconversion*PredictionsSt(4,JNC),HybridESN{3}.BODY_VZ.Approx(PredictionsSt(4,JNC),1),Mkr,'Color',PredStColor,'MarkerFaceColor',PredStColor,'MarkerSize',Mkrsz,LineWidth=2);
    plot(Timeconversion*(PredictionsSt(4,JNC)+(0:120)),HybridESN{3}.BODY_VZ.Approx(PredictionsSt(4,JNC),:),'Color',PredictionColor,LineWidth=2);
end
xlim(Timeconversion*(StartPts(4)+[0,WindowSize]))
ylim([-20,20])
xlabel('Time (s)')
ylabel('Velocity Z')
set(gcf, 'Position',[680   120   845   758])


figure
hold on
for INC = 1:length(B_Fld1)
    for JNC = 1:length(B_Fld1)
        if INC == JNC
            plot((0:120)/60,HybridESN_Confuse{INC,JNC}.VZ.R2,LineStyle="-")
        else
            plot((0:120)/60,HybridESN_Confuse{INC,JNC}.VZ.R2,LineStyle=":")
        end
    end
end
ylim([-1,1])
xlim([0,2])




%%%%%%%%%%%%%%%%%%%%%%%%%%% Confusion matrix
TmpInst = zeros(size(HybridESN_Confuse_In));
Tmpplus5 = zeros(size(HybridESN_Confuse_In));
Tmpplus10 = zeros(size(HybridESN_Confuse_In));
Tmpplus1s = zeros(size(HybridESN_Confuse_In));
Tmpplus2s = zeros(size(HybridESN_Confuse_In));

% ColorRange = FavoriteColors.Map.Sea(length(B_Fld3));
ColorRange = hsv(length(B_Fld3)+1);
ColorRange(3,:) = [];
Mkrs = {'o','^','s','d','p','h'};
Mkrszs = [9,9,11,9,14,14];
figure
hold on
Lentry = cell([1,length(B_Fld3)*length(B_Fld3)]);
Tnames = {'Nostim','Spon','0.5s','0.1s','1.5s','2.0s'};
cntr = 1;
for INC = 1:length(B_Fld3)
    for JNC = 1:length(B_Fld3)
        TmpInst(INC,JNC) = HybridESN_Confuse_In{INC,JNC}.VZ.R2(1);
        Tmpplus5(INC,JNC) = HybridESN_Confuse_In{INC,JNC}.VZ.R2(5);
        Tmpplus10(INC,JNC) = HybridESN_Confuse_In{INC,JNC}.VZ.R2(10);
        Tmpplus1s(INC,JNC) = HybridESN_Confuse_In{INC,JNC}.VZ.R2(61);
        Tmpplus2s(INC,JNC) = HybridESN_Confuse_In{INC,JNC}.VZ.R2(121);

        Lentry{cntr} = strcat('T:',Tnames{INC},' W:',Tnames{JNC});

        if INC == JNC
            % plot((0:120)/60,HybridESN_Confuse_In{INC,JNC}.VZ.R2,LineStyle="-",LineWidth=2,Color=ColorRange(JNC,:),Marker=Mkrs{INC},MarkerFaceColor=ColorRange(JNC,:),MarkerSize=2.5)
            plot((0:120)/60,HybridESN_Confuse_In{INC,JNC}.VZ.R2,LineStyle="-",LineWidth=2,Color=ColorRange(JNC,:))
            plot((0:3:120)/60,HybridESN_Confuse_In{INC,JNC}.VZ.R2(1:3:121),LineStyle="none",...
                LineWidth=2,Color=ColorRange(JNC,:),Marker=Mkrs{INC},MarkerFaceColor=ColorRange(JNC,:),MarkerEdgeColor="none",MarkerSize=Mkrszs(INC))
        else
            % plot((0:120)/60,HybridESN_Confuse_In{INC,JNC}.VZ.R2,LineStyle=":",Color=ColorRange(JNC,:),Marker=Mkrs{INC},MarkerFaceColor=ColorRange(JNC,:),MarkerSize=2.5)
            plot((0:120)/60,HybridESN_Confuse_In{INC,JNC}.VZ.R2,LineStyle=":",LineWidth=2,Color=ColorRange(JNC,:))
            plot((0:3:120)/60,HybridESN_Confuse_In{INC,JNC}.VZ.R2(1:3:121),LineStyle="none",...
                Color=ColorRange(JNC,:),Marker=Mkrs{INC},MarkerFaceColor=ColorRange(JNC,:),MarkerEdgeColor="none",MarkerSize=Mkrszs(INC))
        end
        cntr = cntr+1;
    end
end
ylim([0,1])
xlim([0,2])
legend(Lentry,"NumColumns",length(B_Fld3))
% set(gcf, 'Position',[136         275        1687         644])
set(gcf, 'Position',[136         494        1687         418])

ColorRange = hsv(7);
ColorRange(3,:) = [];
StartPts = [23460+260;4560-20;9900-10];  
WindowSize = 1000-200;
PredictionSep =115; %120   161
PredictionLength = 60;  % 120,60
Timeconversion = 1/60;
PredictionsSt = StartPts + (0:PredictionSep:WindowSize);
PredStColor = FavoriteColors.PptFrm.DarkPink;
Mkr = 's'; % '|' 
Mkrsz = 7;


%%%%% Plot future time behaviors for a set of 2 second pulses
cntr = 1;
figure
tiledlayout(3,1,TileSpacing="tight",Padding="tight")
for INC = [1,2,5]
    nexttile
    if strcmp(B_Fld3{INC},'All')
         plot(Timeconversion*(1:length(CurTargets_In{INC}.VZ(:,1))),CurTargets_In{INC}.VZ(:,1),'k',LineWidth=2);
    else
        plot(Timeconversion*(1:length(CurTargets_In{INC}.BODY_VZ(:,1))),CurTargets_In{INC}.BODY_VZ(:,1),'k',LineWidth=2);
    end
    
    hold on
    for KNC = 1:length(PredictionsSt(1,:))

        for JNC = 1:length(B_Fld3)
            plot(Timeconversion*(PredictionsSt(cntr,KNC)),HybridESN_Confuse_In{INC,JNC}.VZ.Approx(PredictionsSt(cntr,KNC),1),...
                    LineStyle="-",LineWidth=2,Color=ColorRange(JNC,:),Marker=Mkr,MarkerSize=Mkrsz,MarkerFaceColor=ColorRange(JNC,:))
            % if INC == JNC
            %     plot(Timeconversion*(PredictionsSt(cntr,KNC)+(0:120)),HybridESN_Confuse_In{INC,JNC}.VZ.Approx(PredictionsSt(cntr,KNC),:),...
            %         LineStyle="-",LineWidth=2,Color=ColorRange(JNC,:))
            % else
                plot(Timeconversion*(PredictionsSt(cntr,KNC)+(0:PredictionLength)),HybridESN_Confuse_In{INC,JNC}.VZ.Approx(PredictionsSt(cntr,KNC),1:PredictionLength+1),...
                    LineStyle="-",LineWidth=2,Color=ColorRange(JNC,:))
            % end
        end

    end
    xlim(Timeconversion*(StartPts(cntr)+[0,WindowSize]))
    if INC == 1
        ylabel('Velocity Z')
        ylim([-50,90])
    else
        ylim([-12,20])
    end
    xlabel('Time (s)')
    cntr  = cntr +1;
end
% set(gcf, 'Position',[193         434        1615         297])
set(gcf, 'Position',[193         452        1489         297])

%%%%%% Plot only 5 timesteps in the future
cntr = 1;
figure
tiledlayout(3,1,TileSpacing="tight",Padding="tight")
for INC = [2,5,6]
    nexttile
    hold on
    if strcmp(B_Fld3{INC},'All')
        plot(Timeconversion*(1:length(CurTargets_In{INC}.VZ(:,1))),CurTargets_In{INC}.VZ(:,1),'k',LineWidth=2.5);
    else
        plot(Timeconversion*(1:length(CurTargets_In{INC}.BODY_VZ(:,1))),CurTargets_In{INC}.BODY_VZ(:,1),'k',LineWidth=2.5);
    end


    for JNC = 1:length(B_Fld3)
        plot(Timeconversion*(StartPts(cntr) + (0:WindowSize)),HybridESN_Confuse_In{INC,JNC}.VZ.Approx(StartPts(cntr)-6 + (0:WindowSize),5),...
            LineStyle="-",LineWidth=2,Color=ColorRange(JNC,:))
    end


    xlim(Timeconversion*(StartPts(cntr)+[0,WindowSize]))
    if INC == 1
        ylabel('Velocity Z')
        ylim([-12,45])
    else
        ylim([-12,20])
    end
    xlabel('Time (s)')
    cntr  = cntr +1;
end
% set(gcf, 'Position',[193         434        1615         297])
set(gcf, 'Position',[193         452        1489        297])




%%%%%%% Video
% 
% StartPts = [23460+260;4560-20;9900-10;7900-10;9900-10;9900-10];  
% WindowSize = 500;
% PredictionLength = 120;  % 120,60
% Timeconversion = 1/60;
% PredictionsSt = StartPts;
% 
% vframerate = 2*(30*Timeconversion);
% Nframes = 800;
% AlphaProfile = ones(1,Nframes);
% AlphaProfile(1,401:450) =linspace(1,0,50);
% AlphaProfile(1,451:end) =0;
% 
% 
% v = VideoWriter('ConfuseEstimates2.mp4','MPEG-4');
% v.Quality = 95;
% open(v)
% 
% figure
% tiledlayout(6,1,TileSpacing="tight",Padding="tight")
% set(gcf, 'Position',[694    80   777   907])
% for VFR = 1:Nframes
%     StartPts = StartPts+vframerate;
%     PredictionsSt = StartPts;
% 
%     cntr = 1;
%     for INC = [1,2,3,4,5,6]
%         nexttile(cntr)
%         hold off
%         if strcmp(B_Fld3{INC},'All')
%             plot(Timeconversion*(1:length(CurTargets_In{INC}.VZ(:,1))),CurTargets_In{INC}.VZ(:,1),'k',LineWidth=2);
%         else
%             plot(Timeconversion*(1:length(CurTargets_In{INC}.BODY_VZ(:,1))),CurTargets_In{INC}.BODY_VZ(:,1),'k',LineWidth=2);
%         end
% 
%         hold on
%         for KNC = 1:length(PredictionsSt(1,:))
%             for JNC = 1:length(B_Fld3)
%                 if JNC~=2 && JNC~=INC
%                     Curalpha = AlphaProfile(VFR);
%                 else
%                     Curalpha = 1;
%                 end
%                 Curcolor = [ColorRange(JNC,:),Curalpha];
%                 % plot(Timeconversion*(PredictionsSt(cntr,KNC)),HybridESN_Confuse_In{INC,JNC}.VZ.Approx(PredictionsSt(cntr,KNC),1),...
%                 %     LineStyle="-",LineWidth=2,Color=Curcolor,Marker=Mkr,MarkerSize=Mkrsz,MarkerFaceColor=Curcolor)
%               s1 =  scatter(Timeconversion*(PredictionsSt(cntr,KNC)),HybridESN_Confuse_In{INC,JNC}.VZ.Approx(PredictionsSt(cntr,KNC),1),Mkrsz,...
%                     LineWidth=3,Color=ColorRange(JNC,:),Marker=Mkr);
%               s1.MarkerFaceColor = ColorRange(JNC,:);
%               s1.MarkerEdgeColor = ColorRange(JNC,:);
%               alpha(s1,Curalpha);
%                 plot(Timeconversion*(PredictionsSt(cntr,KNC)+(0:PredictionLength)),HybridESN_Confuse_In{INC,JNC}.VZ.Approx(PredictionsSt(cntr,KNC),1:PredictionLength+1),...
%                     LineStyle="-",LineWidth=2,Color=Curcolor)
%             end
% 
%         end
%         xlim(Timeconversion*(StartPts(cntr)+WindowSize*[-2/5,3/5]))
%         if INC == 1
%             ylabel('Velocity Z')
%             ylim([-50,90])
%         else
%             ylabel('Velocity Z')
%             ylim([-12,20])
%         end
%         xlabel('Time (s)')
%         cntr  = cntr +1;
%     end
%     drawnow
%     A = getframe(gcf);
%     writeVideo(v,A)
% end
% close(v)
% 
% 



% 
% StartPts = [23460+260;4560-20;9900-10;7900-10;9900-10;9900-10];
% WindowSize = 300;
% PredictionLength = 120;  % 120,60
% Timeconversion = 1/60;
% PredictionsSt = StartPts;
% 
% vframerate = 2*(30*Timeconversion);
% Nframes = 800;
% AlphaProfile = ones(1,Nframes);
% AlphaProfile(1,401:450) =linspace(1,0,50);
% AlphaProfile(1,451:end) =0;
% 
% v = VideoWriter('ConfuseEstimates3.mp4','MPEG-4');
% v.Quality = 95;
% open(v)
% 
% figure
% tiledlayout(6,3,TileSpacing="tight",Padding="tight")
% set(gcf, 'Position',[206    80   777*2   907])
% 
% for VFR = 1:Nframes
% 
%     StartPts = StartPts+vframerate;
%     PredictionsSt = StartPts;
% 
%     cntr2 = 1;
%     for FutureINDS = [1,7,13]
% 
%         cntr = 1;
%         for INC = 1:6
%             nexttile(3*(INC-1)+cntr2)
%             hold off
% 
%             if strcmp(B_Fld3{INC},'All')
%                 plot(Timeconversion*(1:length(CurTargets_In{INC}.VZ(:,1))),CurTargets_In{INC}.VZ(:,1),'k',LineWidth=2.5);
%             else
%                 plot(Timeconversion*(1:length(CurTargets_In{INC}.BODY_VZ(:,1))),CurTargets_In{INC}.BODY_VZ(:,1),'k',LineWidth=2.5);
%             end
%             hold on
%             for JNC = 1:length(B_Fld3)
%                 if JNC~=2 && JNC~=INC
%                     Curalpha = AlphaProfile(VFR);
%                 else
%                     Curalpha = 1;
%                 end
%                 Curcolor = [ColorRange(JNC,:),Curalpha];
%                 s1 =  scatter(Timeconversion*(StartPts(cntr) +FutureINDS+WindowSize/2),HybridESN_Confuse_In{INC,JNC}.VZ.Approx(StartPts(cntr)+WindowSize/2,FutureINDS),Mkrsz,...
%                     LineWidth=3,Color=ColorRange(JNC,:),Marker=Mkr);
%                 s1.MarkerFaceColor = ColorRange(JNC,:);
%                 s1.MarkerEdgeColor = ColorRange(JNC,:);
%                 alpha(s1,Curalpha);
% 
%                 plot(Timeconversion*(StartPts(cntr)+FutureINDS + (0:WindowSize/2)),HybridESN_Confuse_In{INC,JNC}.VZ.Approx(StartPts(cntr) + (0:WindowSize/2),FutureINDS),...
%                     LineStyle="-",LineWidth=2,Color=Curcolor)
%             end
% 
%             xline(Timeconversion*(StartPts(cntr)+WindowSize/2),':')
%             xlim(Timeconversion*(StartPts(cntr)+[0,WindowSize]))
%             if INC == 1
%                 ylabel('Velocity Z')
%                 ylim([-12,45])
%             else
%                 ylim([-12,20])
%             end
%             xlabel('Time (s)')
%             cntr  = cntr +1;
%         end
%         cntr2 = cntr2 +1;
%     end
%     drawnow
%     A = getframe(gcf);
%     writeVideo(v,A)
% end
% close(v)

end

function [] = Estimate_Best_Types(ESN,PRC,Sequential,Hybrid,ConfigNames,N)
load('MyColors.Mat', 'FavoriteColors');


JFIDs = {'JF41','JF42','JF43','JF44','JF46'};
InputIDs = {'INPUT_1','INPUT_2','INPUT_3','INPUT_4'};
PulsatileIDS = {'JF41','JF42','JF43','JF44','JF46','INPUT_1','INPUT_2','INPUT_3','INPUT_4'};
PulsR2Fields = fields(Hybrid{1}.JF41);
PosFields = PulsR2Fields(1:3);
RotationFields =  PulsR2Fields(4:6);
VelocityFields =  PulsR2Fields(7:9);

AllVFields ={'VX_R2','VY_R2','VZ_R2'};

Horizons = {1:61,62:121};

%%%%% Combine Data By Type for stimulated region
for YNC = 1:length(AllVFields)
     [AVGR2.ESN.(AllVFields{YNC}),STDR2.ESN.(AllVFields{YNC})] = AvgbyQuadrant(Horizons,ESN.All.(AllVFields{YNC}));
    for XNC = 1:length(ConfigNames)
        [AVGR2.(ConfigNames{XNC}).PRC.(AllVFields{YNC}),STDR2.(ConfigNames{XNC}).PRC.(AllVFields{YNC})] = AvgbyQuadrant(Horizons,PRC{XNC}.All.(AllVFields{YNC}));
        [AVGR2.(ConfigNames{XNC}).Sequential.(AllVFields{YNC}),STDR2.(ConfigNames{XNC}).Sequential.(AllVFields{YNC})] = AvgbyQuadrant(Horizons,Sequential{XNC}.All.(AllVFields{YNC}));
        [AVGR2.(ConfigNames{XNC}).Hybrid.(AllVFields{YNC}),STDR2.(ConfigNames{XNC}).Hybrid.(AllVFields{YNC})] = AvgbyQuadrant(Horizons,Hybrid{XNC}.All.(AllVFields{YNC}));
    end
end


%%%%% Combine Data By Type for stimulated region
for YNC = 1:length(PulsatileIDS)
     CurBase = ESN.(PulsatileIDS{YNC});
     Curdata = (CurBase.(PosFields{1})+CurBase.(PosFields{2})+CurBase.(PosFields{3}))/3;
     [AVGR2.ESN.(PulsatileIDS{YNC}).Pos,STDR2.ESN.(PulsatileIDS{YNC}).Pos] = AvgbyQuadrant(Horizons,Curdata);
     Curdata = (CurBase.(RotationFields{1})+CurBase.(RotationFields{2})+CurBase.(RotationFields{3}))/3;
     [AVGR2.ESN.(PulsatileIDS{YNC}).Rot,STDR2.ESN.(PulsatileIDS{YNC}).Rot] = AvgbyQuadrant(Horizons,Curdata);
     Curdata = (CurBase.(VelocityFields{1})+CurBase.(VelocityFields{2})+CurBase.(VelocityFields{3}))/3;
     [AVGR2.ESN.(PulsatileIDS{YNC}).Vel,STDR2.ESN.(PulsatileIDS{YNC}).Vel] = AvgbyQuadrant(Horizons,Curdata);

    for XNC = 1:length(ConfigNames)
        CurBase = PRC{XNC}.(PulsatileIDS{YNC});
        Curdata = (CurBase.(PosFields{1})+CurBase.(PosFields{2})+CurBase.(PosFields{3}))/3;
        [AVGR2.(ConfigNames{XNC}).PRC.(PulsatileIDS{YNC}).Pos,STDR2.(ConfigNames{XNC}).PRC.(PulsatileIDS{YNC}).Pos] = AvgbyQuadrant(Horizons,Curdata);
        Curdata = (CurBase.(RotationFields{1})+CurBase.(RotationFields{2})+CurBase.(RotationFields{3}))/3;
        [AVGR2.(ConfigNames{XNC}).PRC.(PulsatileIDS{YNC}).Rot,STDR2.(ConfigNames{XNC}).PRC.(PulsatileIDS{YNC}).Rot] = AvgbyQuadrant(Horizons,Curdata);
        Curdata = (CurBase.(VelocityFields{1})+CurBase.(VelocityFields{2})+CurBase.(VelocityFields{3}))/3;
        [AVGR2.(ConfigNames{XNC}).PRC.(PulsatileIDS{YNC}).Vel,STDR2.(ConfigNames{XNC}).PRC.(PulsatileIDS{YNC}).Vel] = AvgbyQuadrant(Horizons,Curdata);

        CurBase = Sequential{XNC}.(PulsatileIDS{YNC});
        Curdata = (CurBase.(PosFields{1})+CurBase.(PosFields{2})+CurBase.(PosFields{3}))/3;
        [AVGR2.(ConfigNames{XNC}).Sequential.(PulsatileIDS{YNC}).Pos,STDR2.(ConfigNames{XNC}).Sequential.(PulsatileIDS{YNC}).Pos] = AvgbyQuadrant(Horizons,Curdata);
        Curdata = (CurBase.(RotationFields{1})+CurBase.(RotationFields{2})+CurBase.(RotationFields{3}))/3;
        [AVGR2.(ConfigNames{XNC}).Sequential.(PulsatileIDS{YNC}).Rot,STDR2.(ConfigNames{XNC}).Sequential.(PulsatileIDS{YNC}).Rot] = AvgbyQuadrant(Horizons,Curdata);
        Curdata = (CurBase.(VelocityFields{1})+CurBase.(VelocityFields{2})+CurBase.(VelocityFields{3}))/3;
        [AVGR2.(ConfigNames{XNC}).Sequential.(PulsatileIDS{YNC}).Vel,STDR2.(ConfigNames{XNC}).Sequential.(PulsatileIDS{YNC}).Vel] = AvgbyQuadrant(Horizons,Curdata);

        CurBase = Hybrid{XNC}.(PulsatileIDS{YNC});
        Curdata = (CurBase.(PosFields{1})+CurBase.(PosFields{2})+CurBase.(PosFields{3}))/3;
        [AVGR2.(ConfigNames{XNC}).Hybrid.(PulsatileIDS{YNC}).Pos,STDR2.(ConfigNames{XNC}).Hybrid.(PulsatileIDS{YNC}).Pos] = AvgbyQuadrant(Horizons,Curdata);
        Curdata = (CurBase.(RotationFields{1})+CurBase.(RotationFields{2})+CurBase.(RotationFields{3}))/3;
        [AVGR2.(ConfigNames{XNC}).Hybrid.(PulsatileIDS{YNC}).Rot,STDR2.(ConfigNames{XNC}).Hybrid.(PulsatileIDS{YNC}).Rot] = AvgbyQuadrant(Horizons,Curdata);
        Curdata = (CurBase.(VelocityFields{1})+CurBase.(VelocityFields{2})+CurBase.(VelocityFields{3}))/3;
        [AVGR2.(ConfigNames{XNC}).Hybrid.(PulsatileIDS{YNC}).Vel,STDR2.(ConfigNames{XNC}).Hybrid.(PulsatileIDS{YNC}).Vel] = AvgbyQuadrant(Horizons,Curdata);

        [VELR2.(ConfigNames{XNC}).(VelocityFields{1}).(PulsatileIDS{YNC}),VELR2std.(ConfigNames{XNC}).(VelocityFields{1}).(PulsatileIDS{YNC})] = AvgbyQuadrant(Horizons,CurBase.(VelocityFields{1}));
        [VELR2.(ConfigNames{XNC}).(VelocityFields{2}).(PulsatileIDS{YNC}),VELR2std.(ConfigNames{XNC}).(VelocityFields{2}).(PulsatileIDS{YNC})] = AvgbyQuadrant(Horizons,CurBase.(VelocityFields{2}));
        [VELR2.(ConfigNames{XNC}).(VelocityFields{3}).(PulsatileIDS{YNC}),VELR2std.(ConfigNames{XNC}).(VelocityFields{3}).(PulsatileIDS{YNC})] = AvgbyQuadrant(Horizons,CurBase.(VelocityFields{3}));
    end
end



%%%%%%% Data By Sensor 
BestALLAVG = AvgAllSubFields(AVGR2.Best);
TopALLAVG = AvgAllSubFields(AVGR2.Top);
BottomALLAVG = AvgAllSubFields(AVGR2.Bottom);
RadALLAVG = AvgAllSubFields(AVGR2.Rad);
CrossALLAVG = AvgAllSubFields(AVGR2.Cross);

%%%%%%% Data By Sensor  for the Sequential Form
BestSeqAVG = AvgAllSubFields(AVGR2.Best.Sequential);
TopSeqAVG = AvgAllSubFields(AVGR2.Top.Sequential);
BottomSeqAVG = AvgAllSubFields(AVGR2.Bottom.Sequential);
RadSeqAVG = AvgAllSubFields(AVGR2.Rad.Sequential);
CrossSeqAVG = AvgAllSubFields(AVGR2.Cross.Sequential);

BestSeqSTD = AvgAllSubFields(STDR2.Best.Sequential);
TopSeqSTD = AvgAllSubFields(STDR2.Top.Sequential);
BottomSeqSTD = AvgAllSubFields(STDR2.Bottom.Sequential);
RadSeqSTD = AvgAllSubFields(STDR2.Rad.Sequential);
CrossSeqSTD = AvgAllSubFields(STDR2.Cross.Sequential);

%%%%%%% Data By Sensor  for the ESN Form
ESNALLAVG = AvgAllSubFields(AVGR2.ESN);
ESNALLSTD = AvgAllSubFields(STDR2.ESN);
BestPRCAVG = AvgAllSubFields(AVGR2.Best.PRC);
TopPRCAVG = AvgAllSubFields(AVGR2.Top.PRC);
BottomPRCAVG = AvgAllSubFields(AVGR2.Bottom.PRC);
RadPRCAVG = AvgAllSubFields(AVGR2.Rad.PRC);
CrossPRCAVG = AvgAllSubFields(AVGR2.Cross.PRC);
BestHybridAVG = AvgAllSubFields(AVGR2.Best.Hybrid);

BestPRCSTD = AvgAllSubFields(STDR2.Best.PRC);
TopPRCSTD = AvgAllSubFields(STDR2.Top.PRC);
BottomPRCSTD = AvgAllSubFields(STDR2.Bottom.PRC);
RadPRCSTD = AvgAllSubFields(STDR2.Rad.PRC);
CrossPRCSTD = AvgAllSubFields(STDR2.Cross.PRC);

%%%%%%% Data By Sensor  for the ESN Form
BestHybridVXAVG = AvgAllSubFields(VELR2.Best.P_BODY_VX_R2);
BestHybridVYAVG = AvgAllSubFields(VELR2.Best.P_BODY_VY_R2);
BestHybridVZAVG = AvgAllSubFields(VELR2.Best.P_BODY_VZ_R2);


OrganizationalType_JUSTV = [];
OrganizationalType = [];
STDType_JUSTV = [];
STDType = [];
for YNC = 1:length(PulsatileIDS)
       OrganizationalType_JUSTV = [OrganizationalType_JUSTV,reshape(AVGR2.Best.Hybrid.(PulsatileIDS{YNC}).Vel,[],1)];
       OrganizationalType = [OrganizationalType,reshape(AvgAllSubFields(AVGR2.Best.Hybrid.(PulsatileIDS{YNC})),[],1)];

      STDType_JUSTV = [STDType_JUSTV,reshape(STDR2.Best.Hybrid.(PulsatileIDS{YNC}).Vel,[],1)];
      STDType = [STDType,reshape(AvgAllSubFields(STDR2.Best.Hybrid.(PulsatileIDS{YNC})),[],1)];
end
BestHybridGlobalVelMean = (AVGR2.Best.Hybrid.VX_R2+AVGR2.Best.Hybrid.VY_R2+AVGR2.Best.Hybrid.VZ_R2)/3;
BestHybridGlobalVelSTD = (STDR2.Best.Hybrid.VX_R2+STDR2.Best.Hybrid.VY_R2+STDR2.Best.Hybrid.VZ_R2)/3;
OrganizationalType_JUSTV = [OrganizationalType_JUSTV,reshape(BestHybridGlobalVelMean,[],1)];
STDType_JUSTV = [STDType_JUSTV,reshape(BestHybridGlobalVelSTD,[],1)];

figure
tiledlayout(2,1,"TileSpacing","tight")
nexttile
Bob = bar([reshape(BestSeqAVG,[],1),reshape(TopSeqAVG,[],1),reshape(BottomSeqAVG,[],1),reshape(RadSeqAVG,[],1),reshape(CrossSeqAVG,[],1),reshape(ESNALLAVG,[],1)]);
Barnum = 6;
BarXdata = repmat((1:4)',[1,Barnum]);
BarXdata = BarXdata + (linspace(-(Bob(1).GroupWidth/2-Bob(1).BarWidth/(Barnum*2)),Bob(1).GroupWidth/2-Bob(1).BarWidth/(Barnum*2),Barnum));
hold on
errorbar(BarXdata,[reshape(BestSeqAVG,[],1),reshape(TopSeqAVG,[],1),reshape(BottomSeqAVG,[],1),reshape(RadSeqAVG,[],1),reshape(CrossSeqAVG,[],1),reshape(ESNALLAVG,[],1)],...
    [reshape(BestSeqSTD,[],1),reshape(TopSeqSTD,[],1),reshape(BottomSeqSTD,[],1),reshape(RadSeqSTD,[],1),reshape(CrossSeqSTD,[],1),reshape(ESNALLSTD,[],1)],'.k')
xticklabels({})
ylabel({'Sequential Type','Average R2'})
nexttile
bar([reshape(BestALLAVG,[],1),reshape(TopALLAVG,[],1),reshape(BottomALLAVG,[],1),reshape(RadALLAVG,[],1),reshape(CrossALLAVG,[],1),reshape(ESNALLAVG,[],1)])
xticklabels({'ShortH-ShortMux','ShortH-LongMux','LongH-ShortMux','LongH-LongMux'})
ylabel({'All Types','Average R2'})
set(gcf,"Position",[418         349-100        1047         1.5*420])
legend({ConfigNames{:},'Input'})
colororder(FavoriteColors.Map.Sea(6))



Sensorbars = [reshape(BestSeqAVG,[],1),reshape(BestPRCAVG,[],1),...
    reshape(TopSeqAVG,[],1),reshape(TopPRCAVG,[],1),...
    reshape(BottomSeqAVG,[],1),reshape(BottomPRCAVG,[],1),...
    reshape(RadSeqAVG,[],1),reshape(RadPRCAVG,[],1),...
    reshape(CrossSeqAVG,[],1),reshape(CrossPRCAVG,[],1),...
    reshape(ESNALLAVG,[],1)];
SensorEbars = [reshape(BestSeqSTD,[],1),reshape(BestPRCSTD,[],1),...
    reshape(TopSeqSTD,[],1),reshape(TopPRCSTD,[],1),...
    reshape(BottomSeqSTD,[],1),reshape(BottomPRCSTD,[],1),...
    reshape(RadSeqSTD,[],1),reshape(RadPRCSTD,[],1),...
    reshape(CrossSeqSTD,[],1),reshape(CrossPRCSTD,[],1),...
    reshape(ESNALLSTD,[],1)];
figure
Bob = bar(Sensorbars);
Barnum = 11;
BarXdata = repmat((1:4)',[1,Barnum]);
BarXdata = BarXdata + (linspace(-(Bob(1).GroupWidth/2-Bob(1).BarWidth/(Barnum*2)),Bob(1).GroupWidth/2-Bob(1).BarWidth/(Barnum*2),Barnum));
hold on
errorbar(BarXdata,Sensorbars,...
    SensorEbars,'k','Marker','none','LineStyle','none',LineWidth=1.5)
xticklabels({'ShortH-ShortMux','ShortH-LongMux','LongH-ShortMux','LongH-LongMux'})
ylabel('Average R2')
set(gcf,"Position",[397         400        1047         324])
CorderTmp = zeros([11,3]);
CorderTmp(1:2:10,:) = FavoriteColors.Map.Sea(5);
CorderTmp(2:2:10,:) = FavoriteColors.Map.Sea(5);
CorderTmp(11,:) = [121,25,100]/255;
colororder(CorderTmp)
for Tnc = 1:Barnum
    Bob(Tnc).LineWidth = 1.25;
    if mod(Tnc,2) == 0
        Bob(Tnc).LineStyle = '-.';
    end
end



figure
bar([reshape(ESNALLAVG,[],1),reshape(BestPRCAVG,[],1),reshape(BestSeqAVG,[],1),reshape(BestHybridAVG,[],1)])
xticklabels({'ShortH-ShortMux','ShortH-LongMux','LongH-ShortMux','LongH-LongMux'})
set(gcf,"Position",[418         349        1047         420])
legend({'Pure ESN','PRC','Sequential','Hybrid'})
ylabel({'All Data Best Sensor','Average R2'})
colororder(FavoriteColors.Map.Make2(FavoriteColors.AquaMarine,FavoriteColors.Pink,4))


figure
tiledlayout(2,1,"TileSpacing","tight")
nexttile
Bob = bar(OrganizationalType_JUSTV);
Barnum = size(OrganizationalType_JUSTV,2);
BarXdata = repmat((1:size(OrganizationalType_JUSTV,1))',[1,Barnum]);
BarXdata = BarXdata + (linspace(-(Bob(1).GroupWidth/2-Bob(1).BarWidth/(Barnum*2)),Bob(1).GroupWidth/2-Bob(1).BarWidth/(Barnum*2),Barnum));
hold on
errorbar(BarXdata,OrganizationalType_JUSTV,STDType_JUSTV,'k','Marker','none','LineStyle','none',LineWidth=1.5)
xticklabels({})
ylabel({'Velocity Matching During Stim.','Average R2'})
nexttile
 bar(OrganizationalType);
xticklabels({'ShortH-ShortMux','ShortH-LongMux','LongH-ShortMux','LongH-LongMux'})
ylabel({'All Matching During Stim','Average R2'})
set(gcf,"Position",[418         349-100        1047         1.5*420])
legend(PulsatileIDS,'NumColumns',2)
colororder([FavoriteColors.Map.Make3(FavoriteColors.PptFrm.DeepRed,FavoriteColors.PptFrm.DarkBlue,FavoriteColors.PptFrm.Green,length(PulsatileIDS));[0,0,0]])


figure
bar([reshape(BestHybridVXAVG,[],1),reshape(BestHybridVYAVG,[],1),reshape(BestHybridVZAVG,[],1)])
xticklabels({'ShortH-ShortMux','ShortH-LongMux','LongH-ShortMux','LongH-LongMux'})
set(gcf,"Position",[418         349        1047         420])
legend({'VX','VY','VZ'})
ylabel({'Pulsatile Best Sensor','Average R2'})
colororder([([220,20,60]/255);([0,128,128]/255);([127,255,0]/255)])

end

function [Target_Clean,PureESN,PRC,SequentialESN,HybridESN] = ComputeStateData(Ut,BodySensorData,TargetStateData,Tback,N,TWo,Spect,Uscale)
%%% Just Compute the Pulsatile Data from the ESN 

MPLX_Len = 121;  %%% 61
T_Train_Frac = 1/1; %%% Fraction of data used for training


%%%%%% Set up ESN
 PURE_ESN = Gen_Hybrid_ESN(N,Spect,Uscale,'Activate','Tanh','InpNum',size(Ut,1));
 Sequential_ESN = Gen_Hybrid_ESN(N,Spect,Uscale,'Activate','Tanh','InpNum',size(BodySensorData,1));

%%%%% Set the time and empty vectors
TTotal = length(Ut(1,:));
T = TTotal-TWo;
T_Train = floor((T-MPLX_Len)*T_Train_Frac);

%%%%% Check ESN
Pure_States = PURE_ESN.ESN(TTotal,Ut);
Sequential_States = Sequential_ESN.ESN(TTotal,BodySensorData);


%%% Align First Term and washout
Target_Clean  = Sequential_ESN.Washout_Targets(TWo,TargetStateData(:,Tback:end)');
Pure_States_Clean = PURE_ESN.Washout_Targets(TWo,Pure_States(1:(end-Tback+1),:));
Sequential_States_Clean  =Sequential_ESN.Washout(TWo, Sequential_States(1:(end-Tback+1),:));
PRC_States_Align = BodySensorData(:,1:(end-Tback+1))';
PRC_States_Clean =  [PRC_States_Align((TWo+1):end,:),ones([length(PRC_States_Align((TWo+1):end,1)),1])];
Hybrid_States_Clean = [Sequential_States_Clean,PRC_States_Align((TWo+1):end,:)];


%%% Store all of the states
PureESN.States = Pure_States_Clean;
PRC.States = PRC_States_Clean;
SequentialESN.States = Sequential_States_Clean;
HybridESN.States = Hybrid_States_Clean;

%%%% Get weights
PureESN.W = PURE_ESN.Train(T_Train,Pure_States_Clean,Target_Clean);
PRC.W = Sequential_ESN.Train(T_Train,PRC_States_Clean,Target_Clean);
SequentialESN.W = Sequential_ESN.Train(T_Train,Sequential_States_Clean,Target_Clean);
HybridESN.W = Sequential_ESN.Train(T_Train,Hybrid_States_Clean,Target_Clean);


%%%% Get approximate
PureESN.Approx = Pure_States_Clean*PureESN.W;
PRC.Approx = PRC_States_Clean*PRC.W;
SequentialESN.Approx= Sequential_States_Clean*SequentialESN.W;
HybridESN.Approx = Hybrid_States_Clean*HybridESN.W;

%%%% Get R2
PureESN.R2 = Sequential_ESN.R2_1D(PureESN.Approx,Target_Clean);
PRC.R2 = Sequential_ESN.R2_1D(PRC.Approx,Target_Clean);
SequentialESN.R2= Sequential_ESN.R2_1D(SequentialESN.Approx,Target_Clean);
HybridESN.R2 = Sequential_ESN.R2_1D(HybridESN.Approx,Target_Clean);


end

function [Curavg,varargout] = AvgbyQuadrant(Horizons,Curdata)

Curavg = zeros([length(Horizons),length(Horizons)]);
CurSTD = zeros([length(Horizons),length(Horizons)]);
for INC = 1:length(Horizons)
    for JNC = 1:length(Horizons)
        Curavg(INC,JNC) = mean(Curdata(Horizons{INC},Horizons{JNC}),'all');
        CurSTD(INC,JNC) = std(Curdata(Horizons{INC},Horizons{JNC}),[],'all');
    end
end

if nargout>1
    varargout{1} = CurSTD;
end

end

function [AVG] = AvgAllSubFields(CurSTR)
Counter = 0;  %%% Countes the number of elements in the avg
FldCounter = zeros([100,1]); %%% Counts the number of fields at a level
NumFields =  zeros([100,1]); %%% Stores the total number of fields at a level
CurLevel = 1;  %%% Sets the Current struct level

Curfields = fields(CurSTR);
NumFields(CurLevel) = length(Curfields); %%% Keeps track of the number of felds 
ParentStruct = cell([1,10]);
ParentStruct{1} = CurSTR;

AVG = [0,0;0,0];

%%%%% Keep Running until there is no more data in the structs
while CurLevel > 0
    FldCounter(CurLevel) = FldCounter(CurLevel) +1;
    lastStruct =  ParentStruct{CurLevel};
    Curfields = fields(lastStruct);

    %%%%% Keep Running while there are remaining fields to check at the current level
    while FldCounter(CurLevel) <= NumFields(CurLevel)
        CurData = lastStruct.(Curfields{FldCounter(CurLevel)});

        %%%% Keep lowering the level until it isnt a struct anymore
        while isstruct(CurData)
            Curfields = fields(CurData);
            CurLevel = CurLevel+1;
            NumFields(CurLevel) = length(Curfields);
            FldCounter(CurLevel) = 1;

            ParentStruct{CurLevel-1} = lastStruct;
            lastStruct = CurData;
            CurData = CurData.(Curfields{FldCounter(CurLevel)});
        end

        %%% Sum up the data
        AVG = AVG+CurData;

        Counter = Counter+1;
        FldCounter(CurLevel) = FldCounter(CurLevel) +1;
    end
    CurLevel = CurLevel-1;
end

%%%% Compute the Avgerage
AVG = AVG/Counter;

end

%%%%%%%%%%%%%% GENERAL ESN FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%
function [ESN_OBJ] = Gen_Hybrid_ESN(N,Spect,Input_Scale,varargin)
%Generates an ESN object
% N is number of nodes, Spect ois spectral raduis, 
% Input_scale is the scale of the builtin random input (only used if you create the noise here)
% ActivateFn is the activation function selected

%%%%%% Variable Input Arguments
if  any(strcmp(varargin,'RSeed')) %%% Random seed internal weights
    index = find(strcmp(varargin,'RSeed'));
    RFun_Seed = varargin{index+1} ;
else
    RFun_Seed = 1 ;
end

if  any(strcmp(varargin,'InpRSeed')) %%% Random seed input weights
    index = find(strcmp(varargin,'InpRSeed'));
    RFun_Seed2 = varargin{index+1} ;
else
    RFun_Seed2 = RFun_Seed+1 ;
end

if  any(strcmp(varargin,'InpNum')) %%% Random seed input weights
    index = find(strcmp(varargin,'InpNum'));
    N_Inputs = varargin{index+1} ;
else
      N_Inputs = 1;
end

if any(strcmp(varargin,'Activate')) %%% activation function
    index = find(strcmp(varargin,'Activate'));
    ActivateFn = varargin{index+1} ;
else
    ActivateFn = 'Tanh';
end

if any(strcmp(varargin,'Inbias')) %%% Input bias
    index = find(strcmp(varargin,'Inbias'));
    Offset = varargin{index+1} ;
else
    Offset = 0;
end

if any(strcmp(varargin,'Afun_Param'))  %%% Parameters for the activation function
    index = find(strcmp(varargin,'Afun_Param'));
    alpha = varargin{index+1} ;
else
   alpha = 1;
end

if any(strcmp(varargin,'Init'))  %%% Parameters for the inital condition
    index = find(strcmp(varargin,'Init'));
    ESN_OBJ.Init = varargin{index+1} ;
else
   ESN_OBJ.Init = ones([1,N]);
end

%%%%%%% Parameters related to external nodes %%%%%%%%%
if any(strcmp(varargin,'External_Size'))  %%% Parameter with number of external nodes
    index = find(strcmp(varargin,'External_Size'));
    ESN_OBJ.NEx = varargin{index+1};
    ESN_OBJ.SpectEx = Spect; 
else
   ESN_OBJ.NEx = 0;
   ESN_OBJ.SpectEx = 0;
end

if any(strcmp(varargin,'External_Spect'))  %%%  Parameter with external Spectral density
    index = find(strcmp(varargin,'External_Spect'));
    ESN_OBJ.SpectEx = varargin{index+1};
end

if  any(strcmp(varargin,'ExtRSeed')) %%% Random seed input weights
    index = find(strcmp(varargin,'ExtRSeed'));
    RFun_Seed3 = varargin{index+1} ;
else
    RFun_Seed3 = RFun_Seed2 + 1;
end


%%%%%%%%%   ESN Object  %%%%%%%%%%%

% Save number of Nodes
ESN_OBJ.N = N;
ESN_OBJ.AF = ActivateFn;

% Random value generation
ESN_OBJ.Rnd1D = @(Low,Upp,N) (Low+(Upp-Low)*rand([N,1]));
ESN_OBJ.Rnd2D = @(Low,Upp,N) (Low+(Upp-Low)*rand([N(1),N(2)]));

% Random input weights
rng(RFun_Seed2); % creates a reproducible random set
% ESN_OBJ.InputW = ESN_OBJ.Rnd1D(-Input_Scale,Input_Scale,N); % input weights
ESN_OBJ.InputW = ESN_OBJ.Rnd2D(-Input_Scale,Input_Scale,[N,N_Inputs]); % input weights %%%% For upt to N inputs
ESN_OBJ.InputBias = Offset;

% Set the spectral raduis of the internal weights
rng(RFun_Seed); % creates a reproducible random set
TMP_IW = ESN_OBJ.Rnd2D(-1,1,[N,N]).*ESN_OBJ.Rnd2D(0,1,[N,N]);  % Internal weights
Eval = eig(TMP_IW);
ESN_OBJ.Internal = Spect*TMP_IW /max(abs(Eval));  % Internal weights

% Set the spectral raduis of the Extenal Nodes
if ESN_OBJ.NEx >0 &&  ESN_OBJ.NEx ==N
    rng(RFun_Seed3);
    TMP_IW = ESN_OBJ.Rnd2D(-1,1,[N,ESN_OBJ.NEx]).*ESN_OBJ.Rnd2D(0,1,[ESN_OBJ.NEx,N]);  % Internal weights
    Eval = eig(TMP_IW);
    ESN_OBJ.External = ESN_OBJ.SpectEx*TMP_IW/max(abs(Eval));  % Internal weights
elseif ESN_OBJ.NEx >0
    rng(RFun_Seed3);
    ESN_OBJ.External = ESN_OBJ.SpectEx*ESN_OBJ.Rnd2D(0,1,[N,ESN_OBJ.NEx]);  % Internal weights
else
    ESN_OBJ.External = 0;  % Internal weights
end


ESN_OBJ.alpha  = alpha; 

%%% Run the ESN
if ESN_OBJ.NEx >0
    ESN_OBJ.ESN = @(Time,u,X_External)  (ForFN_ESN(Time,u,X_External,ESN_OBJ));
else
    ESN_OBJ.ESN = @(Time,u)  (ForFN_ESN(Time,u,0,ESN_OBJ));
end

%%% Washout the ESN and set up states for estimate
ESN_OBJ.Washout = @(TWo,X) ( [X((TWo+1):end,:),ones([length(X((TWo+1):end,1)),1])]);
ESN_OBJ.Washout_Targets = @(TWo,Y) (Y((TWo+1):end,:));


%%%% Linear Regression
ESN_OBJ.Train = @(T_Train,X,Y) (X(1:T_Train,:)\Y(1:T_Train,:));
ESN_OBJ.Train_Ridge = @(T_Train,X,Y,ESN_OBJ) ComputRidge(T_Train,X,Y,ESN_OBJ.N+ESN_OBJ.NEx);

ESN_OBJ.NRMSE = @(Y_hat,Y_Eval) (sqrt(mean((Y_hat-Y_Eval).^2))./std(Y_Eval));
ESN_OBJ.R2 = @(Estimate,Real) (1-((sum((Real-Estimate).^2))/(sum((Real - mean(Real)).^2))));
ESN_OBJ.R2_1D = @(Estimate,Real) (1-((sum((Estimate-Real).^2,1))./(sum((Real - mean(Real)).^2,1))));

end

function [Weights] = ComputRidge(T_Train,X,Y,N)
% Runs the ridge regression algorithm based on the 

% X_Train = X(1:T_Train,1:end-1); %%% training without Bias
X_Train = X(1:T_Train,:); %%% training without Bias
Y_Train = Y(1:T_Train,:);
lambda = eig(X_Train'*X_Train); %%% Eigenvalue
Sigma = svd(X_Train); %%% singular value decomposition


%%% Find beta
beta = 1e-13; %% initial conditions
Betas = zeros([1,N+1]);
for df = (N+1):-1:1 %% select df to be whole numbers 
    dF_dbeta = sum(lambda./((lambda+beta).^2));
    F = df-sum(lambda./(lambda+beta));
    Betas(df) = beta -  F/dF_dbeta;  %% Newtons methods 
    beta = Betas(df);
end



%%%  Compute the ridge regression
Weights = (X_Train'*X_Train+Betas*eye(N+1))\(X_Train'*Y_Train); % may need to remove bias though


end

function [X] = ForFN_ESN(Time,u,X_External,ESN_OBJ)
%Generate ESN STate For loop

X = ones([Time,ESN_OBJ.N]);

X(1,:) = ESN_OBJ.Init;

if size(u,1)<Time
    u = u';
end

%%%%% Adjusting size incase ther are no external Nodes
if all(size(X_External) == [1,1])
    X_External = zeros([Time,1]);
else
    X_External = X_External./(max(abs(X_External),[],1));
end

for inc = 2:Time
    Activate_Input = ((ESN_OBJ.Internal*X(inc-1,:)') + (ESN_OBJ.External*X_External(inc-1,:)') + (ESN_OBJ.InputW*u(inc-1,:)' + ESN_OBJ.InputBias))';
    
    %%%% Switches the activation function defaults to Tanh
    if strcmp(ESN_OBJ.AF,'Sigmoid')
        X(inc,:)= 1./(1+exp(-1*Activate_Input));
    elseif strcmp(ESN_OBJ.AF,'Atan')
        X(inc,:)= atan(Activate_Input);
    elseif strcmp(ESN_OBJ.AF,'Ln')
        X(inc,:)= log(1+Activate_Input);
   elseif strcmp(ESN_OBJ.AF,'ELU')
       if  Activate_Input < 0
           X(inc,:)= ESN_OBJ.alpha*exp(Activate_Input-1);
       else
           X(inc,:)= Activate_Input;
       end
    elseif strcmp(ESN_OBJ.AF,'Linear')
        X(inc,:)= Activate_Input;
    elseif strcmp(ESN_OBJ.AF,'Logistic')
        X(inc,:)= ESN_OBJ.alpha*Activate_Input.*(1-Activate_Input);
    else 
        X(inc,:)= tanh(Activate_Input);
    end
end

end




