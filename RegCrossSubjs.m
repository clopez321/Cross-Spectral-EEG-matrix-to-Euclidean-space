clc, clear all, close all

% This Script opens a window to select a folder containing 
% all the subjects to be regularized. The results of the regularization
% goes to the folder "Regularized", located inside the original foder.
% Inside of "Regularized" you can find a copy of all subjects with all
% frequencies regularized.

mainFolder = uigetdir('C:\');%Search for the .mat files 
subjList = dir([mainFolder, '\*.mat']);% List all .mat files
NTapers=20;
NFreq=49;
Nchannels=19;
Nsubjects=size(subjList,1);
Mreg=zeros(Nchannels-1, Nchannels-1, Nsubjects);

NewFolder=fullfile(mainFolder,'Regularized');
mkdir(NewFolder);

for i = 1:Nsubjects
    load(fullfile(mainFolder, subjList(i).name));%Load one Subject
    MCrossH=data_struct.CrossM(:,:,:);%Read all frequencies of one Subj
    for k = 1:NFreq
        H=eye(Nchannels)-ones(Nchannels)/Nchannels;%Apply average reference
        MCrossH(:,:,k)=H*MCrossH(:,:,k)*H;
        MH=MCrossH(1:Nchannels-1,1:Nchannels-1,k);%delete one row and one column

        CrossR = regularizeHS(MH, NTapers);%Apply HS regularization method
        Mreg(:,:,k)=CrossR;
    end
    data_struct.CrossM=Mreg;%Update Subject with regularized Cross-Spectral Data
    save(fullfile(NewFolder,strcat('R_',subjList(i).name)), 'data_struct');%Save Subject in a .mat file on a new folder
end
disp('End of Regularization')
%% 
% Organize the Manifolds by frequencies, for example,
% the first Manifold is formed by the frequency 0.39 Hz of all subjects. 
subjList = dir([NewFolder, '\*.mat']);% List all .mat files of Regularized subjects
Nsubjects=size(subjList,1);
AllSubjFreqs=zeros(Nchannels-1, Nchannels-1, NFreq, Nsubjects);%4-D array that contains all frequencies of all subjects
EuclVects=zeros((Nchannels-1)*(Nchannels)/2,NFreq,Nsubjects);%Dimension of the vectors p(p+1)/2

for i = 1:Nsubjects
    load(fullfile(NewFolder, subjList(i).name));%Load one Subject
    MCrossH=data_struct.CrossM(:,:,:);%Read all frequencies of one Subj
    for k = 1:NFreq
        AllSubjFreqs(:,:,k,i)=MCrossH(:,:,k);
    end 
end
for i = 1:Nsubjects
    [EuclVects(:,:,i)] = Tangent_space(AllSubjFreqs(:,:,:,i)); %Calculate the Riemannian Mean and do the mapping to the tangent space
end
disp('End of Mapping')

