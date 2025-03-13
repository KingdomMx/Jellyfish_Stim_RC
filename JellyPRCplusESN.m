function [] = JellyPRCplusESN()
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


%%%%%% Settings 
MaxInputSensors = length(TestPerformance.Input_ID.SensorPermuations);
ListOfSensors = TestPerformance.Input_ID.BodyData;
FileName = 'JellFusedESNresults4_Scaled.mat';
RunR2SetConfigs = 1; %%%Run all selected sensor Variations of the mux RC 
RunR2OnlyBest = 0; %%%Only the best PRC configurations with the ESN


%%%%% Test 3 things Individual data, aggregate input reconstruction, and all velocity/input reconstrcution 
%%%%% Add Prediction Task Too
Spect= 0.35;
Uscale = 1/4;
BodySpect = 1;
LeakyScaling = 1;  %%% 1 or zeros to scale the input with the number of leaky integrators

%%%%% List af archetypical Sensor Configurations
ConfigNames = {'Best','Top','Bottom','Cross','Rad'};
Sensor_Config =  {'Outer','R2toO2','Inner','Y2toO1';...
    'R2toY2','R2toB2','O2toB2','Y2toO2';...
    'R1toY1','R1toB1','O1toB1','Y1toO1';...
    'R1toO1','R2toO2','Y1toB1','Y2toB2';...
    'R1toR2','Y1toY2','O1toO2','B1toB2'};

%%%%% Check Baseline using Raw PRC and Raw ESN then Combined
if RunR2SetConfigs == 1
    warning('off','all')
    %%%%% Set Up 35 N
    PRC100N = cell(size(ConfigNames));
    Sequential100N = cell(size(ConfigNames));
    Hybrid100N = cell(size(ConfigNames));
    N = 100;
    for QND = 1:length(ConfigNames)
        disp(ConfigNames{QND})
        [PRC100N{QND},Sequential100N{QND},Hybrid100N{QND}] = ComputeFullESNComparison(BodyData,TargetData,Sensor_Config(QND,:),N,Spect,Uscale,LeakyScaling);
    end
    save(FileName ,'PRC100N','Sequential100N',"Hybrid100N",'ConfigNames','Sensor_Config','Spect','Uscale','LeakyScaling','-v7.3')

    %%% The pure esn does not change with sensor inputs
    disp('Pure ESN')
    ESN100N = ComputeFullPureESN(BodyData,TargetData,N,Spect,Uscale,LeakyScaling);
    save(FileName ,'ESN100N','-append')
    warning('on','all')
end

%%%%% Check Using just the best sensory configurations from the same input type
if RunR2OnlyBest == 1
    for IND =1:4
        [InputMatch,VelocityMatch] = Compute_INPUT_HybridPRC(TestPerformance.(strcat('INPUT',num2str(IND))),InputTestCombo{IND},N,Spect,Uscale,BodySpect);
        Hybrid_RC_Results.(strcat('INPUT',num2str(IND))).InputMatch = InputMatch;
        Hybrid_RC_Results.(strcat('INPUT',num2str(IND))).VelocityMatch = VelocityMatch;
    end
    save(FileName ,'Hybrid_RC_Results','-append')
end


Jelly23N = cell(size(ConfigNames));
Jelly35N = cell(size(ConfigNames));
for QND = 1:length(ConfigNames)
    N = 35;
    fprintf('%s %i Node \n',ConfigNames{QND},N)
    Jelly35N{QND} = ComputeOnboardJellySim(BodyData,TargetData,Sensor_Config(QND,:),N,Spect,Uscale);

    N = 23;
    fprintf('\n %s %i Node \n',ConfigNames{QND},N)
    Jelly23N{QND} = ComputeOnboardJellySim(BodyData,TargetData,Sensor_Config(QND,:),N,Spect,Uscale);
end

save('JellyViableRC_Scaled.mat' ,'Jelly35N','Jelly23N','ConfigNames','Sensor_Config','Spect','Uscale','-v7.3')


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

function [PureESNResults] = ComputeFullPureESN(BodyData,TargetData,N,Spect,Uscale,LeakyScaling)
MPLX_Len = 121;  %%% 61
Dt = 0.0167;

TWo_Pulse = 1000; 
TWo_ALL = 10000; 
T_Train_Frac = 1/1; %%% Fraction of data used for training

EmptyMat = nan(MPLX_Len,MPLX_Len); %%% Used to record Accuracy

%%%%%% Full time Targets
[Targ_VX,Targ_X] = MultiplexData(TargetData.All.All.VX,MPLX_Len);
[Targ_VY,Targ_Y] = MultiplexData(TargetData.All.All.VY,MPLX_Len);
[Targ_VZ,Targ_Z] = MultiplexData(TargetData.All.All.VZ,MPLX_Len);
Targ_X = Targ_X*Dt;
Targ_Y = Targ_Y*Dt;
Targ_Z = Targ_Z*Dt;


%%% Define soem empy verctors
PRCResults.All.X_R2 = EmptyMat;
PRCResults.All.Y_R2 = EmptyMat;
PRCResults.All.Z_R2 = EmptyMat;
PureESNResults.All.VX_R2 = EmptyMat;
PureESNResults.All.VY_R2 = EmptyMat;
PureESNResults.All.VZ_R2 = EmptyMat;

%%%% Check both for Pulsatile targets and for Whole time targets
PulsatileFlds = {'JF41','JF42','JF43','JF44','JF46','INPUT_1','INPUT_2','INPUT_3','INPUT_4'};
% PulsatileTPs = fields(BodyData.(PulsatileFlds{1}));
JFFlds = {'JF41','JF42','JF43','JF44','JF46'};
PulsatileTPs = {'P'};
PulseTargets = TargetData.Pulse_Targets;


 PURE_ESN = Gen_Hybrid_ESN(N,Spect,Uscale,'Activate','Tanh','InpNum',MPLX_Len);

 FlippedIND = fliplr(1:MPLX_Len);

%%%% Run through ALL time Data for Sensors
Ut = 2*flipud(MultiplexData(BodyData.All.All.Input,MPLX_Len))-1;

%%%%% Set the time and empty vectors
TTotal = length(Ut(1,:));
T = TTotal-TWo_ALL;
T_Train = floor((T-MPLX_Len)*T_Train_Frac);
CurU = zeros(size(Ut));


%%%%%%%%%%% Run Full Time Data Here %%%%%%%%%%%
 for IND = 1:MPLX_Len
     %%%% Set current observer Data
     disp(num2str(IND))
     CurU(FlippedIND(IND):end,:) = Ut(FlippedIND(IND):end,:);

     %%%%% Check Pure ESN
     if LeakyScaling == 1
         InputScale =1/IND;
     else
         InputScale =1;
     end
     Pure_States = PURE_ESN.ESN(TTotal,CurU*InputScale);

     %%% Align First Term and washout
     TargetXP_Clean  = PURE_ESN.Washout_Targets(TWo_ALL,Targ_X(:,IND:end)');
     TargetYP_Clean  = PURE_ESN.Washout_Targets(TWo_ALL,Targ_Y(:,IND:end)');
     TargetZP_Clean  = PURE_ESN.Washout_Targets(TWo_ALL,Targ_Z(:,IND:end)');
     TargetX_Clean  = PURE_ESN.Washout_Targets(TWo_ALL,Targ_VX(:,IND:end)');
     TargetY_Clean   = PURE_ESN.Washout_Targets(TWo_ALL,Targ_VY(:,IND:end)');
     TargetZ_Clean   = PURE_ESN.Washout_Targets(TWo_ALL,Targ_VZ(:,IND:end)');
     Pure_States_Clean = PURE_ESN.Washout_Targets(TWo_ALL,Pure_States(1:(end-IND+1),:));

     %%%% Get weights, Regress, and calculate error: X
     Approx_Pure = Pure_States_Clean*PURE_ESN.Train(T_Train,Pure_States_Clean,TargetX_Clean);
     PureESNResults.All.VX_R2(IND,:) = PURE_ESN.R2_1D(Approx_Pure,TargetX_Clean);

     %%%% Get weights, Regress, and calculate error: Y
     Approx_Pure = Pure_States_Clean*PURE_ESN.Train(T_Train,Pure_States_Clean,TargetY_Clean);
     PureESNResults.All.VY_R2(IND,:) = PURE_ESN.R2_1D(Approx_Pure,TargetY_Clean);

     %%%% Get weights, Regress, and calculate error: Z
     Approx_Pure = Pure_States_Clean*PURE_ESN.Train(T_Train,Pure_States_Clean,TargetZ_Clean);
     PureESNResults.All.VZ_R2(IND,:) = PURE_ESN.R2_1D(Approx_Pure,TargetZ_Clean);

     %%%%%%%% Dead Reconing Position
     %%%% Get weights, Regress, and calculate error: X
     Approx_Pure = Pure_States_Clean*PURE_ESN.Train(T_Train,Pure_States_Clean,TargetXP_Clean);
     PureESNResults.All.X_R2(IND,:) = PURE_ESN.R2_1D(Approx_Pure,TargetXP_Clean);

     %%%% Get weights, Regress, and calculate error: Y
     Approx_Pure = Pure_States_Clean*PURE_ESN.Train(T_Train,Pure_States_Clean,TargetYP_Clean);
     PureESNResults.All.Y_R2(IND,:) = PURE_ESN.R2_1D(Approx_Pure,TargetYP_Clean);

     %%%% Get weights, Regress, and calculate error: Z
     Approx_Pure = Pure_States_Clean*PURE_ESN.Train(T_Train,Pure_States_Clean,TargetZP_Clean);
     PureESNResults.All.Z_R2(IND,:) = PURE_ESN.R2_1D(Approx_Pure,TargetZP_Clean);

 end

 %%%%%% For all Time and All Jellyfish
 for KND = 1:length(JFFlds)
     %%%%% Set Up Matricies
     PRCResults.All.(JFFlds{KND}).X_R2 = EmptyMat;
     PRCResults.All.(JFFlds{KND}).Y_R2 = EmptyMat;
     PRCResults.All.(JFFlds{KND}).Z_R2 = EmptyMat;
     PureESNResults.All.(JFFlds{KND}).VX_R2 = EmptyMat;
     PureESNResults.All.(JFFlds{KND}).VY_R2 = EmptyMat;
     PureESNResults.All.(JFFlds{KND}).VZ_R2 = EmptyMat;
     [Targ_VX,Targ_X] = MultiplexData(TargetData.All.(JFFlds{KND}).VX,MPLX_Len);
     [Targ_VY,Targ_Y] = MultiplexData(TargetData.All.(JFFlds{KND}).VY,MPLX_Len);
     [Targ_VZ,Targ_Z] = MultiplexData(TargetData.All.(JFFlds{KND}).VZ,MPLX_Len);
     Targ_X = Targ_X*Dt;
     Targ_Y = Targ_Y*Dt;
     Targ_Z = Targ_Z*Dt;

     %%%% Run through ALL time Data for Sensors
     CURBODYDATA = BodyData.All.(JFFlds{KND});
     Ut = 2*flipud(MultiplexData(CURBODYDATA.Input,MPLX_Len))-1;
     TTotal = length(Ut(1,:));
     T = TTotal-TWo_ALL;
     T_Train = floor((T-MPLX_Len)*T_Train_Frac);
     CurU = zeros(size(Ut));
     for IND = 1:MPLX_Len
         disp(num2str(IND))
         %%%% Set current observer Data
         CurU(FlippedIND(IND):end,:) = Ut(FlippedIND(IND):end,:);

         %%%%% Check Pure ESN
         if LeakyScaling == 1
             InputScale =1/IND;
         else
             InputScale =1;
         end
         Pure_States = PURE_ESN.ESN(TTotal,CurU*InputScale);

         %%% Align First Term and washout
         TargetXP_Clean  = PURE_ESN.Washout_Targets(TWo_ALL,Targ_X(:,IND:end)');
         TargetYP_Clean  = PURE_ESN.Washout_Targets(TWo_ALL,Targ_Y(:,IND:end)');
         TargetZP_Clean  = PURE_ESN.Washout_Targets(TWo_ALL,Targ_Z(:,IND:end)');
         TargetX_Clean  = PURE_ESN.Washout_Targets(TWo_ALL,Targ_VX(:,IND:end)');
         TargetY_Clean   = PURE_ESN.Washout_Targets(TWo_ALL,Targ_VY(:,IND:end)');
         TargetZ_Clean   = PURE_ESN.Washout_Targets(TWo_ALL,Targ_VZ(:,IND:end)');
         Pure_States_Clean = PURE_ESN.Washout_Targets(TWo_ALL,Pure_States(1:(end-IND+1),:));

         %%%% Get weights, Regress, and calculate error: X
         Approx_Pure = Pure_States_Clean*PURE_ESN.Train(T_Train,Pure_States_Clean,TargetX_Clean);
         PureESNResults.All.(JFFlds{KND}).VX_R2(IND,:) = PURE_ESN.R2_1D(Approx_Pure,TargetX_Clean);

         %%%% Get weights, Regress, and calculate error: Y
         Approx_Pure = Pure_States_Clean*PURE_ESN.Train(T_Train,Pure_States_Clean,TargetY_Clean);
         PureESNResults.All.(JFFlds{KND}).VY_R2(IND,:) = PURE_ESN.R2_1D(Approx_Pure,TargetY_Clean);

         %%%% Get weights, Regress, and calculate error: Z
         Approx_Pure = Pure_States_Clean*PURE_ESN.Train(T_Train,Pure_States_Clean,TargetZ_Clean);
         PureESNResults.All.(JFFlds{KND}).VZ_R2(IND,:) = PURE_ESN.R2_1D(Approx_Pure,TargetZ_Clean);

         %%%%%%%% Dead Reconing Position
         %%%% Get weights, Regress, and calculate error: X
         Approx_Pure = Pure_States_Clean*PURE_ESN.Train(T_Train,Pure_States_Clean,TargetXP_Clean);
         PureESNResults.All.(JFFlds{KND}).X_R2(IND,:) = PURE_ESN.R2_1D(Approx_Pure,TargetXP_Clean);

         %%%% Get weights, Regress, and calculate error: Y
         Approx_Pure = Pure_States_Clean*PURE_ESN.Train(T_Train,Pure_States_Clean,TargetYP_Clean);
         PureESNResults.All.(JFFlds{KND}).Y_R2(IND,:) = PURE_ESN.R2_1D(Approx_Pure,TargetYP_Clean);

         %%%% Get weights, Regress, and calculate error: Z
         Approx_Pure = Pure_States_Clean*PURE_ESN.Train(T_Train,Pure_States_Clean,TargetZP_Clean);
         PureESNResults.All.(JFFlds{KND}).Z_R2(IND,:) = PURE_ESN.R2_1D(Approx_Pure,TargetZP_Clean);
     end
 end


%%%% For all of the Pulasatile Data
 for JND = 1:length(PulsatileFlds)
     for KND = 1:length(PulsatileTPs)
         CURBODYDATA = BodyData.(PulsatileFlds{JND}).(PulsatileTPs{KND});
         Ut = 2*flipud(MultiplexData(CURBODYDATA.Input,MPLX_Len))-1;

         %%%%% Cut off UT and based on number of states
         TTotal = length(Ut(1,:));
         T = TTotal-TWo_Pulse;
         T_Train = floor((T-MPLX_Len)*T_Train_Frac);

         for IND = 1:MPLX_Len  
             disp(num2str(IND))
             CurU = zeros(size(Ut));
             CurU(FlippedIND(IND):end,:) = Ut(FlippedIND(IND):end,:);


             %%%%% Check Pure ESN
             Pure_States = PURE_ESN.ESN(TTotal,CurU);

             %%%% Check ALL Targets
             CURTARGETDATA = TargetData.(PulsatileFlds{JND}).(PulsatileTPs{KND});
             for LND = 1:length(PulseTargets)
                 %%%%% Set up Storage
                 CurTarget = MultiplexData(CURTARGETDATA.(PulseTargets{LND}),MPLX_Len);

                 %%% Align First Term
                 CurTarget_Align  = CurTarget(:,IND:end)';
                 Pure_States_Align = Pure_States(1:(end-IND+1),:);


                 %%% Washout
                 Target_Clean  = PURE_ESN.Washout_Targets(TWo_Pulse,CurTarget_Align);
                 Pure_States_Clean = PURE_ESN.Washout(TWo_Pulse,Pure_States_Align);


                 %%%%% Training
                 W_Pure = PURE_ESN.Train(T_Train,Pure_States_Clean,Target_Clean);


                 %%%%% Compute Approximation
                 Approx_Pure = Pure_States_Clean*W_Pure;

                 %%%%% Compute NRMSE
                 R2_Pure = PURE_ESN.R2_1D(Approx_Pure,Target_Clean);


                 %%%%%% Seting Up Results
                 if IND == 1
                     PureESNResults.(PulsatileFlds{JND}).(strcat(PulsatileTPs{KND},'_',PulseTargets{LND},'_R2')) = EmptyMat;
                 end

                 %%%%%% Saving Results data
                 PureESNResults.(PulsatileFlds{JND}).(strcat(PulsatileTPs{KND},'_',PulseTargets{LND},'_R2'))(IND,:) = R2_Pure;
             end
         end
     end
 end

end

function [JellyResults] = ComputeOnboardJellySim(BodyData,TargetData,SenseConfig,N,Spect,Uscale)
%%% Just Compute the Pulsatile Data from the ESN 

MPLX_Len = 121;  %%% 61
TWo_Pulse = 1000; 
TWo_ALL = 10000; 
T_Train_Frac = 1/1; %%% Fraction of data used for training
Dt = 0.0167;

%%%% Check both for Pulsatile targets and for Whole time targets
JFFlds = {'JF41','JF42','JF43','JF44','JF46'};
PulsatileFlds = {'JF41','JF42','JF43','JF44','JF46','INPUT_1','INPUT_2','INPUT_3','INPUT_4'};
PulsatileTPs = {'P'};
NumSensors = length(SenseConfig);
FlippedIND = fliplr(1:MPLX_Len);
EmptyMat = nan(MPLX_Len,MPLX_Len); %%% Used to record Accuracy

Jelly_ESN = Gen_Hybrid_ESN(N,Spect,Uscale,'Activate','Tanh','InpNum',NumSensors);

%%%%%% For all Time and All Jellyfish
for KND = 1:length(JFFlds)
    %%%%% Set Up Matricies
    JellyResults.All.(JFFlds{KND}).X_R2 = EmptyMat;
    JellyResults.All.(JFFlds{KND}).Y_R2 = EmptyMat;
    JellyResults.All.(JFFlds{KND}).Z_R2 = EmptyMat;
    JellyResults.All.(JFFlds{KND}).X_Err = EmptyMat;
    JellyResults.All.(JFFlds{KND}).Y_Err = EmptyMat;
    JellyResults.All.(JFFlds{KND}).Z_Err = EmptyMat;
    JellyResults.All.(JFFlds{KND}).VX_R2 = EmptyMat;
    JellyResults.All.(JFFlds{KND}).VY_R2 = EmptyMat;
    JellyResults.All.(JFFlds{KND}).VZ_R2 = EmptyMat;

    %%%%% Set Up Targets
    [Targ_VX,Targ_X] = MultiplexData(TargetData.All.(JFFlds{KND}).VX,MPLX_Len);
    [Targ_VY,Targ_Y] = MultiplexData(TargetData.All.(JFFlds{KND}).VY,MPLX_Len);
    [Targ_VZ,Targ_Z] = MultiplexData(TargetData.All.(JFFlds{KND}).VZ,MPLX_Len);
    Targ_X = Targ_X*Dt;
    Targ_Y = Targ_Y*Dt;
    Targ_Z = Targ_Z*Dt;

    %%%%% Load Acceleration Data
    Accel_Body = [TargetData.All.(JFFlds{KND}).AX(1:end-MPLX_Len);...
        TargetData.All.(JFFlds{KND}).AY(1:end-MPLX_Len);TargetData.All.(JFFlds{KND}).AZ(1:end-MPLX_Len)]';

    %%%% Run through ALL time Data for Sensors
    CURBODYDATA = BodyData.All.(JFFlds{KND});
    Ut = 2*flipud(MultiplexData(CURBODYDATA.Input,MPLX_Len))-1;
    BodySenseData = zeros(MPLX_Len*NumSensors,size(Ut,2));
    ESN_INPUT = zeros(NumSensors,size(Ut,2));
    for LND = 1:NumSensors
        %%%% Normalize sensor data [-1,1]
        TMPSens = CURBODYDATA.(SenseConfig{LND});
        TMPSens = 2*(TMPSens-range(TMPSens)/2-min(TMPSens))/range(TMPSens);
        BodySenseData(MPLX_Len*(LND-1)+(1:MPLX_Len),:) = flipud(MultiplexData(TMPSens,MPLX_Len));
        ESN_INPUT(LND,:) = TMPSens(1:end-MPLX_Len);
    end

    %%%%% Set the time and empty vectors
    TTotal = length(Ut(1,:));
    T = TTotal-TWo_ALL;
    T_Train = floor((T-MPLX_Len)*T_Train_Frac);

    %%%%% Check ESN
    ESN_Base_States = Jelly_ESN.ESN(TTotal,ESN_INPUT);

    %%%%%%%%%%% Run Full Time Data Here %%%%%%%%%%%
  fprintf('\n %s \n',JFFlds{KND})
    for IND = 1:MPLX_Len
        fprintf('%i,',IND)
        %%%% Set current observer Data
        PRC_States = zeros(IND*NumSensors,size(Ut,2));
        for LND = 1:NumSensors
            PRC_States(IND*(LND-1)+(1:IND),:) = BodySenseData(MPLX_Len*(LND-1)+(FlippedIND(IND):MPLX_Len),:);
        end

        %%% Align First Term and washout
        TargetXP_Clean  = Jelly_ESN.Washout_Targets(TWo_ALL,Targ_X(:,IND:end)');
        TargetYP_Clean  = Jelly_ESN.Washout_Targets(TWo_ALL,Targ_Y(:,IND:end)');
        TargetZP_Clean  = Jelly_ESN.Washout_Targets(TWo_ALL,Targ_Z(:,IND:end)');
        TargetX_Clean  = Jelly_ESN.Washout_Targets(TWo_ALL,Targ_VX(:,IND:end)');
        TargetY_Clean  = Jelly_ESN.Washout_Targets(TWo_ALL,Targ_VY(:,IND:end)');
        TargetZ_Clean  = Jelly_ESN.Washout_Targets(TWo_ALL,Targ_VZ(:,IND:end)');

        ESN_Base_States_Clean  =Jelly_ESN.Washout(TWo_ALL, ESN_Base_States(1:(end-IND+1),:));
        PRC_States_Align = PRC_States(:,1:(end-IND+1))';
        Accel_Body_States = Accel_Body(1:length(ESN_Base_States(1:(end-IND+1),1)),:);
        Jelly_States_Clean = [ESN_Base_States_Clean,PRC_States_Align((TWo_ALL+1):end,:),Accel_Body_States((TWo_ALL+1):end,:)];

        %%%% Get weights, Regress, and calculate error: X
        Approx_Jelly = Jelly_States_Clean*Jelly_ESN.Train(T_Train,Jelly_States_Clean,TargetX_Clean);
        JellyResults.All.(JFFlds{KND}).VX_R2(IND,:) = Jelly_ESN.R2_1D( Approx_Jelly,TargetX_Clean);

        %%%% Get weights, Regress, and calculate error: Y
        Approx_Jelly = Jelly_States_Clean*Jelly_ESN.Train(T_Train,Jelly_States_Clean,TargetY_Clean);
        JellyResults.All.(JFFlds{KND}).VY_R2(IND,:) = Jelly_ESN.R2_1D( Approx_Jelly,TargetY_Clean);

        %%%% Get weights, Regress, and calculate error: Z
        Approx_Jelly = Jelly_States_Clean*Jelly_ESN.Train(T_Train,Jelly_States_Clean,TargetZ_Clean);
        JellyResults.All.(JFFlds{KND}).VZ_R2(IND,:) = Jelly_ESN.R2_1D( Approx_Jelly,TargetZ_Clean);

        %%%%%%%% Dead Reconing Position
        %%%% Get weights, Regress, and calculate error: X
        Approx_Jelly = Jelly_States_Clean*Jelly_ESN.Train(T_Train,Jelly_States_Clean,TargetXP_Clean);
        JellyResults.All.(JFFlds{KND}).X_R2(IND,:) = Jelly_ESN.R2_1D(Approx_Jelly,TargetXP_Clean);
        JellyResults.All.(JFFlds{KND}).X_Err(IND,:) = vecnorm(Approx_Jelly-TargetXP_Clean,2)/length(Approx_Jelly(:,1));

        %%%% Get weights, Regress, and calculate error: Y
        Approx_Jelly = Jelly_States_Clean*Jelly_ESN.Train(T_Train,Jelly_States_Clean,TargetYP_Clean);
        JellyResults.All.(JFFlds{KND}).Y_R2(IND,:) = Jelly_ESN.R2_1D(Approx_Jelly,TargetYP_Clean);
        JellyResults.All.(JFFlds{KND}).Y_Err(IND,:) = vecnorm(Approx_Jelly-TargetYP_Clean,2)/length(Approx_Jelly(:,1));

        %%%% Get weights, Regress, and calculate error: Z
        Approx_Jelly = Jelly_States_Clean*Jelly_ESN.Train(T_Train,Jelly_States_Clean,TargetZP_Clean);
        JellyResults.All.(JFFlds{KND}).Z_R2(IND,:) = Jelly_ESN.R2_1D(Approx_Jelly,TargetZP_Clean);
        JellyResults.All.(JFFlds{KND}).Z_Err(IND,:) = vecnorm(Approx_Jelly-TargetZP_Clean,2)/length(Approx_Jelly(:,1));
        
    end
end


%%%% For all of the Pulasatile Data
 fprintf('\n')
disp('Running Pulsatile Time Data')
for JND = 1:length(PulsatileFlds)
    for KND = 1:length(PulsatileTPs)
        CURBODYDATA = BodyData.(PulsatileFlds{JND}).(PulsatileTPs{KND});
        Ut = 2*flipud(MultiplexData(CURBODYDATA.Input,MPLX_Len))-1;
        BodySenseData = zeros(MPLX_Len*NumSensors,size(Ut,2));
        ESN_INPUT = zeros(NumSensors,size(Ut,2));
        for LND = 1:NumSensors
            %%%% Normalize sensor data [-1,1]
            TMPSens = CURBODYDATA.(SenseConfig{LND});
            TMPSens = 2*(TMPSens-range(TMPSens)/2-min(TMPSens))/range(TMPSens);
            BodySenseData(MPLX_Len*(LND-1)+(1:MPLX_Len),:) = flipud(MultiplexData(TMPSens,MPLX_Len));
            ESN_INPUT(LND,:) = TMPSens(1:end-MPLX_Len);
        end

        %%%%% Cut off UT and based on number of states
        TTotal = length(Ut(1,:));
        T = TTotal-TWo_Pulse;
        T_Train = floor((T-MPLX_Len)*T_Train_Frac);

        %%%%% Check ESN
        ESN_Base_States = Jelly_ESN.ESN(TTotal,ESN_INPUT);
        
        %%%%% Set Up Targets
        CURTARGETDATA = TargetData.(PulsatileFlds{JND}).(PulsatileTPs{KND});
        AccelData = [gradient(CURTARGETDATA.BODY_VX'),gradient(CURTARGETDATA.BODY_VY'),gradient(CURTARGETDATA.BODY_VZ')]/Dt;
        Accel_Body = AccelData(2:end-(MPLX_Len-1),:);

        [Targ_VX,Targ_X] = MultiplexData(CURTARGETDATA.BODY_VX,MPLX_Len);
        [Targ_VY,Targ_Y] = MultiplexData(CURTARGETDATA.BODY_VY,MPLX_Len);
        [Targ_VZ,Targ_Z] = MultiplexData(CURTARGETDATA.BODY_VZ,MPLX_Len);
        Targ_X = Targ_X*Dt;
        Targ_Y = Targ_Y*Dt;
        Targ_Z = Targ_Z*Dt;

        fprintf('\n %s \n',PulsatileFlds{JND})
        for IND = 1:MPLX_Len
            fprintf('%i,',IND)
            PRC_States = zeros(IND*NumSensors,size(Ut,2));
            for LND = 1:NumSensors
                PRC_States(IND*(LND-1)+(1:IND),:) = BodySenseData(MPLX_Len*(LND-1)+(FlippedIND(IND):MPLX_Len),:);
            end

            %%% Align First Term and Washout Targets
            VXTarget_Clean = Jelly_ESN.Washout_Targets(TWo_Pulse,Targ_VX(:,IND:end)');
            VYTarget_Clean = Jelly_ESN.Washout_Targets(TWo_Pulse,Targ_VY(:,IND:end)');
            VZTarget_Clean = Jelly_ESN.Washout_Targets(TWo_Pulse,Targ_VZ(:,IND:end)');
            XTarget_Clean = Jelly_ESN.Washout_Targets(TWo_Pulse,Targ_X(:,IND:end)');
            YTarget_Clean = Jelly_ESN.Washout_Targets(TWo_Pulse,Targ_Y(:,IND:end)');
            ZTarget_Clean = Jelly_ESN.Washout_Targets(TWo_Pulse,Targ_Z(:,IND:end)');

            %%% Align First Term and Washout RC
            ESN_Base_States_Clean  =Jelly_ESN.Washout(TWo_Pulse, ESN_Base_States(1:(end-IND+1),:));
            PRC_States_Align = PRC_States(:,1:(end-IND+1))';
            Accel_Body_States = Accel_Body(1:length(ESN_Base_States(1:(end-IND+1),1)),:);
            Jelly_States_Clean = [ESN_Base_States_Clean,PRC_States_Align((TWo_Pulse+1):end,:),Accel_Body_States((TWo_Pulse+1):end,:)];

            %%%%% Training
            W_VX = Jelly_ESN.Train(T_Train,Jelly_States_Clean,VXTarget_Clean);
            W_VY = Jelly_ESN.Train(T_Train,Jelly_States_Clean,VYTarget_Clean);
            W_VZ = Jelly_ESN.Train(T_Train,Jelly_States_Clean,VZTarget_Clean);
            W_X = Jelly_ESN.Train(T_Train,Jelly_States_Clean,XTarget_Clean);
            W_Y = Jelly_ESN.Train(T_Train,Jelly_States_Clean,YTarget_Clean);
            W_Z = Jelly_ESN.Train(T_Train,Jelly_States_Clean,ZTarget_Clean);

            %%%%% Compute Approximation
            Approx_VX = Jelly_States_Clean*W_VX;
            Approx_VY = Jelly_States_Clean*W_VY;
            Approx_VZ = Jelly_States_Clean*W_VZ;
            Approx_X = Jelly_States_Clean*W_X;
            Approx_Y = Jelly_States_Clean*W_Y;
            Approx_Z = Jelly_States_Clean*W_Z;

            %%%%% Compute NRMSE
            R2_VX = Jelly_ESN.R2_1D(Approx_VX,VXTarget_Clean);
            R2_VY = Jelly_ESN.R2_1D(Approx_VY,VYTarget_Clean);
            R2_VZ = Jelly_ESN.R2_1D(Approx_VZ,VZTarget_Clean);
            R2_X = Jelly_ESN.R2_1D(Approx_X,XTarget_Clean);
            R2_Y = Jelly_ESN.R2_1D(Approx_Y,YTarget_Clean);
            R2_Z = Jelly_ESN.R2_1D(Approx_Z,ZTarget_Clean);


            %%%%%% Seting Up Results
            if IND == 1
                JellyResults.(PulsatileFlds{JND}).X_R2 = EmptyMat;
                JellyResults.(PulsatileFlds{JND}).Y_R2 = EmptyMat;
                JellyResults.(PulsatileFlds{JND}).Z_R2 = EmptyMat;
                JellyResults.(PulsatileFlds{JND}).X_Err = EmptyMat;
                JellyResults.(PulsatileFlds{JND}).Y_Err = EmptyMat;
                JellyResults.(PulsatileFlds{JND}).Z_Err = EmptyMat;
                JellyResults.(PulsatileFlds{JND}).VX_R2 = EmptyMat;
                JellyResults.(PulsatileFlds{JND}).VY_R2 = EmptyMat;
                JellyResults.(PulsatileFlds{JND}).VZ_R2 = EmptyMat;
            end

            %%%%%% Saving Results data
            JellyResults.(PulsatileFlds{JND}).VX_R2(IND,:) = R2_VX;
            JellyResults.(PulsatileFlds{JND}).VY_R2(IND,:) = R2_VY;
            JellyResults.(PulsatileFlds{JND}).VZ_R2(IND,:) = R2_VZ;
            JellyResults.(PulsatileFlds{JND}).X_R2(IND,:) = R2_X;
            JellyResults.(PulsatileFlds{JND}).Y_R2(IND,:) = R2_Y;
            JellyResults.(PulsatileFlds{JND}).Z_R2(IND,:) = R2_Z;
            JellyResults.(PulsatileFlds{JND}).X_Err(IND,:) = vecnorm(Approx_X-XTarget_Clean,2)/length(Approx_X(:,1));
            JellyResults.(PulsatileFlds{JND}).Y_Err(IND,:) = vecnorm(Approx_Y-YTarget_Clean,2)/length(Approx_Y(:,1));
            JellyResults.(PulsatileFlds{JND}).Z_Err(IND,:) = vecnorm(Approx_Z-ZTarget_Clean,2)/length(Approx_Z(:,1));


        end
    end
end


end

function [PRCResults,SequentialESNResults,HybridESNResults] = ComputeFullESNComparison(BodyData,TargetData,SenseConfig,N,Spect,Uscale,LeakyScaling)
%%% Just Compute the Pulsatile Data from the ESN 

MPLX_Len = 121;  %%% 61

TWo_Pulse = 1000; 
TWo_ALL = 10000; 
T_Train_Frac = 1/1; %%% Fraction of data used for training
Dt = 0.0167;


%%%%%%%%%%%%%% Targets (t or t+1): Input Determination, Velocity 
%%%%%%%%%%%%%% Targets (C+1): Next Input Signal, Next Body Pulses COM, Motion, Velocity (with and without stim)

%%%%% Check PRC and Check ESN First Both have 1 input
%%%%% Check hybrid

EmptyMat = nan(MPLX_Len,MPLX_Len); %%% Used to record Accuracy

%%%%%% Full time Targets
[Targ_VX,Targ_X] = MultiplexData(TargetData.All.All.VX,MPLX_Len);
[Targ_VY,Targ_Y] = MultiplexData(TargetData.All.All.VY,MPLX_Len);
[Targ_VZ,Targ_Z] = MultiplexData(TargetData.All.All.VZ,MPLX_Len);
Targ_X = Targ_X*Dt;
Targ_Y = Targ_Y*Dt;
Targ_Z = Targ_Z*Dt;

%%% Define soem empy verctors
PRCResults.All.X_R2 = EmptyMat;
PRCResults.All.Y_R2 = EmptyMat;
PRCResults.All.Z_R2 = EmptyMat;
PRCResults.All.VX_R2 = EmptyMat;
PRCResults.All.VY_R2 = EmptyMat;
PRCResults.All.VZ_R2 = EmptyMat;
SequentialESNResults = PRCResults;
HybridESNResults = PRCResults;

%%%% Check both for Pulsatile targets and for Whole time targets
PulsatileFlds = {'JF41','JF42','JF43','JF44','JF46','INPUT_1','INPUT_2','INPUT_3','INPUT_4'};
JFFlds = {'JF41','JF42','JF43','JF44','JF46'};
% PulsatileTPs = fields(BodyData.(PulsatileFlds{1}));
PulsatileTPs = {'P'};
PulseTargets = TargetData.Pulse_Targets;
NumSensors = length(SenseConfig);

 Sequential_ESN = Gen_Hybrid_ESN(N,Spect,Uscale,'Activate','Tanh','InpNum',MPLX_Len*NumSensors);

 FlippedIND = fliplr(1:MPLX_Len);

%%%% Run through ALL time Data for Sensors
Ut = 2*flipud(MultiplexData(BodyData.All_Norm.All.Input,MPLX_Len))-1;
BodySenseData = zeros(MPLX_Len*NumSensors,size(Ut,2));
for LND = 1:NumSensors
    %%%% Normalize sensor data [-1,1]
    TMPSens = BodyData.All_Norm.All.(SenseConfig{LND});
    BodySenseData(MPLX_Len*(LND-1)+(1:MPLX_Len),:) = flipud(MultiplexData(TMPSens,MPLX_Len));
end

%%%%% Set the time and empty vectors
TTotal = length(Ut(1,:));
T = TTotal-TWo_ALL;
T_Train = floor((T-MPLX_Len)*T_Train_Frac);
CurSense = zeros(size(BodySenseData));

%%%%%%%%%%% Run Full Time Data Here %%%%%%%%%%%
 fprintf('\n')
disp('Running Full Time Data')
 for IND = 1:MPLX_Len
   fprintf('%i,',IND)
     %%%% Set current observer Data
     PRC_States = zeros(IND*NumSensors,size(Ut,2));
     for LND = 1:NumSensors
         PRC_States(IND*(LND-1)+(1:IND),:) = BodySenseData(MPLX_Len*(LND-1)+(FlippedIND(IND):MPLX_Len),:);
         CurSense(MPLX_Len*(LND-1)+(FlippedIND(IND):MPLX_Len),:) = BodySenseData(MPLX_Len*(LND-1)+(FlippedIND(IND):MPLX_Len),:);
     end

     %%%%% Check ESN
     if LeakyScaling == 1
         InputScale =1/IND;
     else
         InputScale =1;
     end
     Sequential_States = Sequential_ESN.ESN(TTotal,CurSense*InputScale);

     %%% Align First Term and washout
     TargetXP_Clean  = Sequential_ESN.Washout_Targets(TWo_ALL,Targ_X(:,IND:end)');
     TargetYP_Clean  = Sequential_ESN.Washout_Targets(TWo_ALL,Targ_Y(:,IND:end)');
     TargetZP_Clean  = Sequential_ESN.Washout_Targets(TWo_ALL,Targ_Z(:,IND:end)');
     TargetX_Clean  = Sequential_ESN.Washout_Targets(TWo_ALL,Targ_VX(:,IND:end)');
     TargetY_Clean  = Sequential_ESN.Washout_Targets(TWo_ALL,Targ_VY(:,IND:end)');
     TargetZ_Clean  = Sequential_ESN.Washout_Targets(TWo_ALL,Targ_VZ(:,IND:end)');
     Sequential_States_Clean  =Sequential_ESN.Washout(TWo_ALL, Sequential_States(1:(end-IND+1),:));
     PRC_States_Align =PRC_States(:,1:(end-IND+1))';
     PRC_States_Clean =  [PRC_States_Align((TWo_ALL+1):end,:),ones([length(PRC_States_Align((TWo_ALL+1):end,1)),1])];
     Hybrid_States_Clean = [Sequential_States_Clean,PRC_States_Align((TWo_ALL+1):end,:)];

     %%%%%%%% Velocity
     %%%% Get weights, Regress, and calculate error: X
     Approx_Seq = Sequential_States_Clean*Sequential_ESN.Train(T_Train,Sequential_States_Clean,TargetX_Clean);
     Approx_PRC = PRC_States_Clean*Sequential_ESN.Train(T_Train,PRC_States_Clean,TargetX_Clean);
     Approx_Hyb = Hybrid_States_Clean*Sequential_ESN.Train(T_Train,Hybrid_States_Clean,TargetX_Clean);
     SequentialESNResults.All.VX_R2(IND,:) = Sequential_ESN.R2_1D(Approx_Seq,TargetX_Clean);
     PRCResults.All.VX_R2(IND,:) = Sequential_ESN.R2_1D( Approx_PRC,TargetX_Clean);
     HybridESNResults.All.VX_R2(IND,:) = Sequential_ESN.R2_1D( Approx_Hyb,TargetX_Clean);

     %%%% Get weights, Regress, and calculate error: Y
     Approx_Seq = Sequential_States_Clean*Sequential_ESN.Train(T_Train,Sequential_States_Clean,TargetY_Clean);
     Approx_Hyb = Hybrid_States_Clean*Sequential_ESN.Train(T_Train,Hybrid_States_Clean,TargetY_Clean);
     Approx_PRC = PRC_States_Clean*Sequential_ESN.Train(T_Train,PRC_States_Clean,TargetY_Clean);
     SequentialESNResults.All.VY_R2(IND,:) = Sequential_ESN.R2_1D(Approx_Seq,TargetY_Clean);
     PRCResults.All.VY_R2(IND,:) = Sequential_ESN.R2_1D( Approx_PRC,TargetY_Clean);
     HybridESNResults.All.VY_R2(IND,:) = Sequential_ESN.R2_1D( Approx_Hyb,TargetY_Clean);

     %%%% Get weights, Regress, and calculate error: Z
     Approx_Seq = Sequential_States_Clean*Sequential_ESN.Train(T_Train,Sequential_States_Clean,TargetZ_Clean);
     Approx_PRC = PRC_States_Clean*Sequential_ESN.Train(T_Train,PRC_States_Clean,TargetZ_Clean);
     Approx_Hyb = Hybrid_States_Clean*Sequential_ESN.Train(T_Train,Hybrid_States_Clean,TargetZ_Clean);
     SequentialESNResults.All.VZ_R2(IND,:) = Sequential_ESN.R2_1D(Approx_Seq,TargetZ_Clean);
     PRCResults.All.VZ_R2(IND,:) = Sequential_ESN.R2_1D( Approx_PRC,TargetZ_Clean);
     HybridESNResults.All.VZ_R2(IND,:) = Sequential_ESN.R2_1D( Approx_Hyb,TargetZ_Clean);

     %%%%%%%% Dead Reconing Position
     %%%% Get weights, Regress, and calculate error: X
     Approx_Seq = Sequential_States_Clean*Sequential_ESN.Train(T_Train,Sequential_States_Clean,TargetXP_Clean);
     Approx_PRC = PRC_States_Clean*Sequential_ESN.Train(T_Train,PRC_States_Clean,TargetXP_Clean);
     Approx_Hyb = Hybrid_States_Clean*Sequential_ESN.Train(T_Train,Hybrid_States_Clean,TargetXP_Clean);
     SequentialESNResults.All.X_R2(IND,:) = Sequential_ESN.R2_1D(Approx_Seq,TargetXP_Clean);
     PRCResults.All.X_R2(IND,:) = Sequential_ESN.R2_1D( Approx_PRC,TargetXP_Clean);
     HybridESNResults.All.X_R2(IND,:) = Sequential_ESN.R2_1D( Approx_Hyb,TargetXP_Clean);

     %%%% Get weights, Regress, and calculate error: Y
     Approx_Seq = Sequential_States_Clean*Sequential_ESN.Train(T_Train,Sequential_States_Clean,TargetYP_Clean);
     Approx_Hyb = Hybrid_States_Clean*Sequential_ESN.Train(T_Train,Hybrid_States_Clean,TargetYP_Clean);
     Approx_PRC = PRC_States_Clean*Sequential_ESN.Train(T_Train,PRC_States_Clean,TargetYP_Clean);
     SequentialESNResults.All.Y_R2(IND,:) = Sequential_ESN.R2_1D(Approx_Seq,TargetYP_Clean);
     PRCResults.All.Y_R2(IND,:) = Sequential_ESN.R2_1D( Approx_PRC,TargetYP_Clean);
     HybridESNResults.All.Y_R2(IND,:) = Sequential_ESN.R2_1D( Approx_Hyb,TargetYP_Clean);

     %%%% Get weights, Regress, and calculate error: Z
     Approx_Seq = Sequential_States_Clean*Sequential_ESN.Train(T_Train,Sequential_States_Clean,TargetZP_Clean);
     Approx_PRC = PRC_States_Clean*Sequential_ESN.Train(T_Train,PRC_States_Clean,TargetZP_Clean);
     Approx_Hyb = Hybrid_States_Clean*Sequential_ESN.Train(T_Train,Hybrid_States_Clean,TargetZP_Clean);
     SequentialESNResults.All.Z_R2(IND,:) = Sequential_ESN.R2_1D(Approx_Seq,TargetZP_Clean);
     PRCResults.All.Z_R2(IND,:) = Sequential_ESN.R2_1D( Approx_PRC,TargetZP_Clean);
     HybridESNResults.All.Z_R2(IND,:) = Sequential_ESN.R2_1D( Approx_Hyb,TargetZP_Clean);

 end


%%%%%% For all Time and All Jellyfish
for KND = 1:length(JFFlds)
    %%%%% Set Up Matricies
    PRCResults.All.(JFFlds{KND}).X_R2 = EmptyMat;
    PRCResults.All.(JFFlds{KND}).Y_R2 = EmptyMat;
    PRCResults.All.(JFFlds{KND}).Z_R2 = EmptyMat;
    PRCResults.All.(JFFlds{KND}).VX_R2 = EmptyMat;
    PRCResults.All.(JFFlds{KND}).VY_R2 = EmptyMat;
    PRCResults.All.(JFFlds{KND}).VZ_R2 = EmptyMat;
    SequentialESNResults.All.(JFFlds{KND}) = PRCResults.All.(JFFlds{KND});
    HybridESNResults.All.(JFFlds{KND}) = PRCResults.All.(JFFlds{KND});
    [Targ_VX,Targ_X] = MultiplexData(TargetData.All.(JFFlds{KND}).VX,MPLX_Len);
    [Targ_VY,Targ_Y] = MultiplexData(TargetData.All.(JFFlds{KND}).VY,MPLX_Len);
    [Targ_VZ,Targ_Z] = MultiplexData(TargetData.All.(JFFlds{KND}).VZ,MPLX_Len);
    Targ_X = Targ_X*Dt;
    Targ_Y = Targ_Y*Dt;
    Targ_Z = Targ_Z*Dt;

    %%%% Run through ALL time Data for Sensors
    CURBODYDATA = BodyData.All.(JFFlds{KND});
    Ut = 2*flipud(MultiplexData(CURBODYDATA.Input,MPLX_Len))-1;
    BodySenseData = zeros(MPLX_Len*NumSensors,size(Ut,2));
    for LND = 1:NumSensors
        %%%% Normalize sensor data [-1,1]
        TMPSens = CURBODYDATA.(SenseConfig{LND});
        TMPSens = 2*(TMPSens-range(TMPSens)/2-min(TMPSens))/range(TMPSens);
        BodySenseData(MPLX_Len*(LND-1)+(1:MPLX_Len),:) = flipud(MultiplexData(TMPSens,MPLX_Len));
    end

    %%%%% Set the time and empty vectors
    TTotal = length(Ut(1,:));
    T = TTotal-TWo_ALL;
    T_Train = floor((T-MPLX_Len)*T_Train_Frac);
    CurSense = zeros(size(BodySenseData));

    %%%%%%%%%%% Run Full Time Data Here %%%%%%%%%%%
    fprintf('\n %s \n',JFFlds{KND})
    for IND = 1:MPLX_Len
        fprintf('%i,',IND)
        %%%% Set current observer Data
        PRC_States = zeros(IND*NumSensors,size(Ut,2));
        for LND = 1:NumSensors
            PRC_States(IND*(LND-1)+(1:IND),:) = BodySenseData(MPLX_Len*(LND-1)+(FlippedIND(IND):MPLX_Len),:);
            CurSense(MPLX_Len*(LND-1)+(FlippedIND(IND):MPLX_Len),:) = BodySenseData(MPLX_Len*(LND-1)+(FlippedIND(IND):MPLX_Len),:);
        end

        %%%%% Check ESN
        Sequential_States = Sequential_ESN.ESN(TTotal,CurSense);

        %%% Align First Term and washout
        TargetXP_Clean  = Sequential_ESN.Washout_Targets(TWo_ALL,Targ_X(:,IND:end)');
        TargetYP_Clean  = Sequential_ESN.Washout_Targets(TWo_ALL,Targ_Y(:,IND:end)');
        TargetZP_Clean  = Sequential_ESN.Washout_Targets(TWo_ALL,Targ_Z(:,IND:end)');
        TargetX_Clean  = Sequential_ESN.Washout_Targets(TWo_ALL,Targ_VX(:,IND:end)');
        TargetY_Clean  = Sequential_ESN.Washout_Targets(TWo_ALL,Targ_VY(:,IND:end)');
        TargetZ_Clean  = Sequential_ESN.Washout_Targets(TWo_ALL,Targ_VZ(:,IND:end)');
        Sequential_States_Clean  =Sequential_ESN.Washout(TWo_ALL, Sequential_States(1:(end-IND+1),:));
        PRC_States_Align =PRC_States(:,1:(end-IND+1))';
        PRC_States_Clean =  [PRC_States_Align((TWo_ALL+1):end,:),ones([length(PRC_States_Align((TWo_ALL+1):end,1)),1])];
        Hybrid_States_Clean = [Sequential_States_Clean,PRC_States_Align((TWo_ALL+1):end,:)];

        %%%% Get weights, Regress, and calculate error: X
        Approx_Seq = Sequential_States_Clean*Sequential_ESN.Train(T_Train,Sequential_States_Clean,TargetX_Clean);
        Approx_PRC = PRC_States_Clean*Sequential_ESN.Train(T_Train,PRC_States_Clean,TargetX_Clean);
        Approx_Hyb = Hybrid_States_Clean*Sequential_ESN.Train(T_Train,Hybrid_States_Clean,TargetX_Clean);
        SequentialESNResults.All.(JFFlds{KND}).VX_R2(IND,:) = Sequential_ESN.R2_1D(Approx_Seq,TargetX_Clean);
        PRCResults.All.(JFFlds{KND}).VX_R2(IND,:) = Sequential_ESN.R2_1D( Approx_PRC,TargetX_Clean);
        HybridESNResults.All.(JFFlds{KND}).VX_R2(IND,:) = Sequential_ESN.R2_1D( Approx_Hyb,TargetX_Clean);

        %%%% Get weights, Regress, and calculate error: Y
        Approx_Seq = Sequential_States_Clean*Sequential_ESN.Train(T_Train,Sequential_States_Clean,TargetY_Clean);
        Approx_Hyb = Hybrid_States_Clean*Sequential_ESN.Train(T_Train,Hybrid_States_Clean,TargetY_Clean);
        Approx_PRC = PRC_States_Clean*Sequential_ESN.Train(T_Train,PRC_States_Clean,TargetY_Clean);
        SequentialESNResults.All.(JFFlds{KND}).VY_R2(IND,:) = Sequential_ESN.R2_1D(Approx_Seq,TargetY_Clean);
        PRCResults.All.(JFFlds{KND}).VY_R2(IND,:) = Sequential_ESN.R2_1D( Approx_PRC,TargetY_Clean);
        HybridESNResults.All.(JFFlds{KND}).VY_R2(IND,:) = Sequential_ESN.R2_1D( Approx_Hyb,TargetY_Clean);

        %%%% Get weights, Regress, and calculate error: Z
        Approx_Seq = Sequential_States_Clean*Sequential_ESN.Train(T_Train,Sequential_States_Clean,TargetZ_Clean);
        Approx_PRC = PRC_States_Clean*Sequential_ESN.Train(T_Train,PRC_States_Clean,TargetZ_Clean);
        Approx_Hyb = Hybrid_States_Clean*Sequential_ESN.Train(T_Train,Hybrid_States_Clean,TargetZ_Clean);
        SequentialESNResults.All.(JFFlds{KND}).VZ_R2(IND,:) = Sequential_ESN.R2_1D(Approx_Seq,TargetZ_Clean);
        PRCResults.All.(JFFlds{KND}).VZ_R2(IND,:) = Sequential_ESN.R2_1D( Approx_PRC,TargetZ_Clean);
        HybridESNResults.All.(JFFlds{KND}).VZ_R2(IND,:) = Sequential_ESN.R2_1D( Approx_Hyb,TargetZ_Clean);


        %%%%%%%% Dead Reconing Position
        %%%% Get weights, Regress, and calculate error: X
        Approx_Seq = Sequential_States_Clean*Sequential_ESN.Train(T_Train,Sequential_States_Clean,TargetXP_Clean);
        Approx_PRC = PRC_States_Clean*Sequential_ESN.Train(T_Train,PRC_States_Clean,TargetXP_Clean);
        Approx_Hyb = Hybrid_States_Clean*Sequential_ESN.Train(T_Train,Hybrid_States_Clean,TargetXP_Clean);
        SequentialESNResults.All.(JFFlds{KND}).X_R2(IND,:) = Sequential_ESN.R2_1D(Approx_Seq,TargetXP_Clean);
        PRCResults.All.(JFFlds{KND}).X_R2(IND,:) = Sequential_ESN.R2_1D( Approx_PRC,TargetXP_Clean);
        HybridESNResults.All.(JFFlds{KND}).X_R2(IND,:) = Sequential_ESN.R2_1D( Approx_Hyb,TargetXP_Clean);

        %%%% Get weights, Regress, and calculate error: Y
        Approx_Seq = Sequential_States_Clean*Sequential_ESN.Train(T_Train,Sequential_States_Clean,TargetYP_Clean);
        Approx_Hyb = Hybrid_States_Clean*Sequential_ESN.Train(T_Train,Hybrid_States_Clean,TargetYP_Clean);
        Approx_PRC = PRC_States_Clean*Sequential_ESN.Train(T_Train,PRC_States_Clean,TargetYP_Clean);
        SequentialESNResults.All.(JFFlds{KND}).Y_R2(IND,:) = Sequential_ESN.R2_1D(Approx_Seq,TargetYP_Clean);
        PRCResults.All.(JFFlds{KND}).Y_R2(IND,:) = Sequential_ESN.R2_1D( Approx_PRC,TargetYP_Clean);
        HybridESNResults.All.(JFFlds{KND}).Y_R2(IND,:) = Sequential_ESN.R2_1D( Approx_Hyb,TargetYP_Clean);

        %%%% Get weights, Regress, and calculate error: Z
        Approx_Seq = Sequential_States_Clean*Sequential_ESN.Train(T_Train,Sequential_States_Clean,TargetZP_Clean);
        Approx_PRC = PRC_States_Clean*Sequential_ESN.Train(T_Train,PRC_States_Clean,TargetZP_Clean);
        Approx_Hyb = Hybrid_States_Clean*Sequential_ESN.Train(T_Train,Hybrid_States_Clean,TargetZP_Clean);
        SequentialESNResults.All.(JFFlds{KND}).Z_R2(IND,:) = Sequential_ESN.R2_1D(Approx_Seq,TargetZP_Clean);
        PRCResults.All.(JFFlds{KND}).Z_R2(IND,:) = Sequential_ESN.R2_1D( Approx_PRC,TargetZP_Clean);
        HybridESNResults.All.(JFFlds{KND}).Z_R2(IND,:) = Sequential_ESN.R2_1D( Approx_Hyb,TargetZP_Clean);

    end
end


%%%% For all of the Pulasatile Data
 fprintf('\n')
disp('Running Pulsatile Time Data')
 for JND = 1:length(PulsatileFlds)
     for KND = 1:length(PulsatileTPs)
         CURBODYDATA = BodyData.(PulsatileFlds{JND}).(PulsatileTPs{KND});
         Ut = 2*flipud(MultiplexData(CURBODYDATA.Input,MPLX_Len))-1;
         BodySenseData = zeros(MPLX_Len*NumSensors,size(Ut,2));
         for LND = 1:NumSensors
             %%%% Normalize sensor data [-1,1]
             TMPSens = CURBODYDATA.(SenseConfig{LND});
             TMPSens = 2*(TMPSens-range(TMPSens)/2-min(TMPSens))/range(TMPSens);
             %%%% Multiplex Input and flip so it shift backwards in time T+n to T
             BodySenseData(MPLX_Len*(LND-1)+(1:MPLX_Len),:) = flipud(MultiplexData(TMPSens,MPLX_Len));
         end

         %%%%% Cut off UT and based on number of states
         TTotal = length(Ut(1,:));
         T = TTotal-TWo_Pulse;
         T_Train = floor((T-MPLX_Len)*T_Train_Frac);

         fprintf('\n %s \n',PulsatileFlds{JND})
         for IND = 1:MPLX_Len  
              fprintf('%i,',IND)
             CurSense = zeros(size(BodySenseData));
             for LND = 1:NumSensors
                 CurSense(MPLX_Len*(LND-1)+(FlippedIND(IND):MPLX_Len),:) = BodySenseData(MPLX_Len*(LND-1)+(FlippedIND(IND):MPLX_Len),:);
             end

             %%%%% Check Pure ESN
             if LeakyScaling == 1
                 InputScale =1/IND;
             else
                 InputScale =1;
             end
             Sequential_States = Sequential_ESN.ESN(TTotal,CurSense*InputScale);

             %%%% Check ALL Targets
             CURTARGETDATA = TargetData.(PulsatileFlds{JND}).(PulsatileTPs{KND});
             for LND = 1:length(PulseTargets)
                 %%%%% Set up Storage
                 CurTarget = MultiplexData(CURTARGETDATA.(PulseTargets{LND}),MPLX_Len);

                 %%% Align First Term
                 CurTarget_Align  = CurTarget(:,IND:end)';
                 Sequential_States_Align  = Sequential_States(1:(end-IND+1),:);
                 PRC_States_Align = CurSense(:,1:(end-IND+1))';

                 %%% Washout
                 Target_Clean  = Sequential_ESN.Washout_Targets(TWo_Pulse,CurTarget_Align);
                 Sequential_States_Clean  = Sequential_ESN.Washout(TWo_Pulse,Sequential_States_Align);
                 PRC_States_Clean =  [PRC_States_Align((TWo_Pulse+1):end,:),ones([length(PRC_States_Align((TWo_Pulse+1):end,1)),1])];
                 Hybrid_States_Clean = [Sequential_States_Clean,PRC_States_Align((TWo_Pulse+1):end,:)];

                 %%%%% Training
                 W_Seq = Sequential_ESN.Train(T_Train,Sequential_States_Clean,Target_Clean);
                 W_PRC = Sequential_ESN.Train(T_Train,PRC_States_Clean,Target_Clean);
                 W_Hybrid = Sequential_ESN.Train(T_Train,Hybrid_States_Clean,Target_Clean);

                 %%%%% Compute Approximation
                 Approx_Seq = Sequential_States_Clean*W_Seq;
                 Approx_PRC = PRC_States_Clean*W_PRC;
                 Approx_Hyb = Hybrid_States_Clean*W_Hybrid;

                 %%%%% Compute NRMSE
                 R2_Seq = Sequential_ESN.R2_1D(Approx_Seq,Target_Clean);
                 R2_PRC = Sequential_ESN.R2_1D( Approx_PRC,Target_Clean);
                 R2_Hyb = Sequential_ESN.R2_1D( Approx_Hyb,Target_Clean);


                 %%%%%% Seting Up Results
                 if IND == 1
                     PRCResults.(PulsatileFlds{JND}).(strcat(PulsatileTPs{KND},'_',PulseTargets{LND},'_R2')) = EmptyMat;
                     SequentialESNResults.(PulsatileFlds{JND}).(strcat(PulsatileTPs{KND},'_',PulseTargets{LND},'_R2')) = EmptyMat;
                     HybridESNResults.(PulsatileFlds{JND}).(strcat(PulsatileTPs{KND},'_',PulseTargets{LND},'_R2')) = EmptyMat;
                 end

                 %%%%%% Saving Results data
                 PRCResults.(PulsatileFlds{JND}).(strcat(PulsatileTPs{KND},'_',PulseTargets{LND},'_R2'))(IND,:) = R2_PRC;
                 SequentialESNResults.(PulsatileFlds{JND}).(strcat(PulsatileTPs{KND},'_',PulseTargets{LND},'_R2'))(IND,:) = R2_Seq;
                 HybridESNResults.(PulsatileFlds{JND}).(strcat(PulsatileTPs{KND},'_',PulseTargets{LND},'_R2'))(IND,:) = R2_Hyb;
             end
         end
     end
 end

end

function [InputMatch,VelocityMatch] = Compute_INPUT_HybridPRC(Teststruct,TestData,N,Spect,Uscale,BodySpect)
%%% Just Compute the Pulsatile Data from the ESN 

Stimtype = 'StimPulse'; %%% ''StimPulse', 'StimPulseSS'
TrainType = 'Linear'; %%% 'Linear' 'Ridge' 

%%%% Determining Maximum Number of Sensors on Jelly
MaxInputSensors = length(Teststruct.SensorPermuations);
ListOfSensors = Teststruct.BodyData;

%%%%% Setting RC input and Time Variables
Ut = Teststruct.StimPulse.InputTarget';
TTotal = length(Ut);
TWo = 10; %floor(TTotal/1000)
T = TTotal-TWo;
T_Train = floor(T/10);

%%%%% Setting RC Targets
InputTarget =Teststruct.(Stimtype).InputTarget;
VelocityTarget = Teststruct.(Stimtype).VelocityTarget;

%%%%% Creating Storage Elements
InputMatch.T1.Weights = cell([1,MaxInputSensors+1]);
InputMatch.T1.CleanStates = cell([1,MaxInputSensors+1]);
InputMatch.T1.Approx = cell([1,MaxInputSensors+1]);
InputMatch.T1.NMSE = ones([1,MaxInputSensors+1]);
InputMatch.T1.R2 = ones([1,MaxInputSensors+1]);
VelocityMatch = InputMatch;


%%%% Get the baseline with no Body Data
ESN_OBJ = Gen_Hybrid_ESN(N+MaxInputSensors,Spect,Uscale,'Activate','Tanh');
States_Clean = ESN_OBJ.Washout(TWo,ESN_OBJ.ESN(TWo+T,Ut));
          %%%%%%%  Washing Targets 
Targets_Clean_IN = ESN_OBJ.Washout_Targets(TWo,InputTarget);
Targets_Clean_V = ESN_OBJ.Washout_Targets(TWo,VelocityTarget);

%%%% Finisheing Baseline ESN
W_IN = ESN_OBJ.Train(T_Train,States_Clean,Targets_Clean_IN);
W_V = ESN_OBJ.Train(T_Train,States_Clean,Targets_Clean_V);
ApproxTarget_IN = States_Clean*W_IN;
ApproxTarget_V = States_Clean*W_V;
InputMatch.T1.Weights{1} = W_IN;
InputMatch.T1.CleanStates{1} = Targets_Clean_IN;
InputMatch.T1.Approx{1} = ApproxTarget_IN;
InputMatch.T1.NMSE(1) =mean((ApproxTarget_IN-Targets_Clean_IN).^2)./std(Targets_Clean_IN.^2);
VelocityMatch.T1.Weights{1} = W_V;
VelocityMatch.T1.CleanStates{1} = InputMatch.T1.CleanStates{1};
VelocityMatch.T1.Approx{1} = ApproxTarget_V;
VelocityMatch.T1.NMSE(1) = mean(mean((ApproxTarget_V-Targets_Clean_V).^2)./std(Targets_Clean_V.^2));
InputMatch.T1.R2(1) =  ESN_OBJ.R2(Targets_Clean_IN,ApproxTarget_IN);
VelocityMatch.T1.R2(1) =  ESN_OBJ.R2(Targets_Clean_V,ApproxTarget_V);
InputMatch.T2 = InputMatch.T1;
VelocityMatch.T2 = VelocityMatch.T1;



%%%%% Storing Cleaned Targets     
InputMatch.Target  = Targets_Clean_IN;
VelocityMatch.Target  = Targets_Clean_V;

for INC = 1:MaxInputSensors
    %%%%%% Finding the current optimal sensor names
    CurBestSensor_IN = {ListOfSensors{Teststruct.(Stimtype).BestPerformance.Input{INC}}};
    CurBestSensor_V ={ListOfSensors{Teststruct.(Stimtype).BestPerformance.Velocity{INC}}};

    %%%%%% Writing the Optimal Sensor names as an external Input
    Xternal_IN = [];
    Xternal_V = [];
    for JNC = 1:INC
        Xternal_IN = [Xternal_IN,reshape(TestData.(Stimtype).(CurBestSensor_IN{JNC}).Pulses,[],1)];
        Xternal_V= [Xternal_V,reshape(TestData.(Stimtype).(CurBestSensor_V{JNC}).Pulses,[],1)];
    end
    Xternal_IN_Normal =Xternal_IN./(max(abs(Xternal_IN),[],1));
    Xternal_V_Normal = Xternal_V./(max(abs(Xternal_IN),[],1));

    %%%%% Creating ESN Object
    ESN_OBJ_T1 = Gen_Hybrid_ESN(N,Spect,Uscale,'Activate','Tanh','External_Size',INC,'External_Spect',BodySpect);
    ESN_OBJ_T2 = Gen_Hybrid_ESN(N,Spect,Uscale,'Activate','Tanh','InpNum',INC);

    %%%%%% Running ESN for both external Node Types and Defining full ESN States
    Hybrid_States_IN = [ESN_OBJ_T1.ESN(TWo+T,Ut,Xternal_IN),Xternal_IN];
    Hybrid_States_Vel = [ESN_OBJ_T1.ESN(TWo+T,Ut,Xternal_V),Xternal_V];
    Sequential_States_IN = [ESN_OBJ_T2.ESN(TWo+T,Xternal_IN_Normal),Xternal_IN];
    Sequential_States_Vel = [ESN_OBJ_T2.ESN(TWo+T,Xternal_V_Normal),Xternal_V];

    %%%%%%% Defining Washed States with Bias
    States_Clean_and_Bias_IN = ESN_OBJ_T1.Washout(TWo,Hybrid_States_IN);
    States_Clean_and_Bias_V = ESN_OBJ_T1.Washout(TWo, Hybrid_States_Vel);
    States_Clean_and_Bias_IN2 = ESN_OBJ_T2.Washout(TWo, Sequential_States_IN);
    States_Clean_and_Bias_V2 = ESN_OBJ_T2.Washout(TWo, Sequential_States_Vel);
    
    %%%%%%%  Getting Weights From Regression
    if strcmpi(TrainType,'Ridge')
        W_IN = ESN_OBJ_T1.Train_Ridge(T_Train,States_Clean_and_Bias_IN,Targets_Clean_IN,ESN_OBJ_T1);
        W_V = ESN_OBJ_T1.Train_Ridge(T_Train,States_Clean_and_Bias_V,Targets_Clean_V,ESN_OBJ_T1);
        W_IN2 = ESN_OBJ_T2.Train_Ridge(T_Train,States_Clean_and_Bias_IN2,Targets_Clean_IN,ESN_OBJ_T2);
        W_V2 = ESN_OBJ_T2.Train_Ridge(T_Train,States_Clean_and_Bias_V2,Targets_Clean_V,ESN_OBJ_T2);
    else
        W_IN = ESN_OBJ_T1.Train(T_Train,States_Clean_and_Bias_IN,Targets_Clean_IN);
        W_V = ESN_OBJ_T1.Train(T_Train,States_Clean_and_Bias_V,Targets_Clean_V);
        W_IN2 = ESN_OBJ_T2.Train(T_Train,States_Clean_and_Bias_IN2,Targets_Clean_IN);
        W_V2 = ESN_OBJ_T2.Train(T_Train,States_Clean_and_Bias_V2,Targets_Clean_V);
    end

    %%%%%%%  Approximating States
    ApproxTarget_IN = States_Clean_and_Bias_IN*W_IN;
    ApproxTarget_V = States_Clean_and_Bias_V*W_V;
    ApproxTarget_IN2 = States_Clean_and_Bias_IN2*W_IN2;
    ApproxTarget_V2 = States_Clean_and_Bias_V2*W_V2;

    %%%%%%% Recording the Results
    InputMatch.T1.Weights{INC+1} = W_IN;
    InputMatch.T1.CleanStates{INC+1} = States_Clean_and_Bias_IN;
    InputMatch.T1.Approx{INC+1} = ApproxTarget_IN;
    InputMatch.T1.NMSE(INC+1) =mean((ApproxTarget_IN-Targets_Clean_IN).^2)./std(Targets_Clean_IN.^2);
    InputMatch.T1.R2(INC+1) =  ESN_OBJ_T1.R2(Targets_Clean_IN,ApproxTarget_IN);
    InputMatch.T2.Weights{INC+1} = W_IN2;
    InputMatch.T2.CleanStates{INC+1} = States_Clean_and_Bias_IN2;
    InputMatch.T2.Approx{INC+1} = ApproxTarget_IN2;
    InputMatch.T2.NMSE(INC+1) =mean((ApproxTarget_IN2-Targets_Clean_IN).^2)./std(Targets_Clean_IN.^2);
    InputMatch.T2.R2(INC+1) =  ESN_OBJ_T2.R2(Targets_Clean_IN,ApproxTarget_IN2);


    VelocityMatch.T1.Weights{INC+1} = W_V;
    VelocityMatch.T1.CleanStates{INC+1} = States_Clean_and_Bias_V;
    VelocityMatch.T1.Approx{INC+1} = ApproxTarget_V;
    VelocityMatch.T1.NMSE(INC+1) = mean(mean((ApproxTarget_V-Targets_Clean_V).^2)./std(Targets_Clean_V.^2));
    VelocityMatch.T1.R2(INC+1) =  ESN_OBJ_T1.R2(Targets_Clean_V,ApproxTarget_V);
    VelocityMatch.T2.Weights{INC+1} = W_V2;
    VelocityMatch.T2.CleanStates{INC+1} = States_Clean_and_Bias_V2;
    VelocityMatch.T2.Approx{INC+1} = ApproxTarget_V2;
    VelocityMatch.T2.NMSE(INC+1) = mean(mean((ApproxTarget_V2-Targets_Clean_V).^2)./std(Targets_Clean_V.^2));
    VelocityMatch.T2.R2(INC+1) =  ESN_OBJ_T2.R2(Targets_Clean_V,ApproxTarget_V2);

end


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




