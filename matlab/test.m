% test code for NCCA

clear all;
close all;

% create data
Y = rand(100,10);
Z = rand(100,15);

% create model with standard settings
model = modelCreate('ncca',[],[],Y,Z);

% add dynamic model
options = gpOptions;
model = nccaAddDynamics(model,'gp','YZ',options);

% train model 
% modelOptimise seems to be broken I've been trying to fix it but
% it does seem to dislike being called multiple times
model = nccaOptimise(model,true,30);

% try different outputs
fprintf('Multiple Mode Out:\n');
modelOut(model,Y(1:1:3,:),'YtoZ');
fprintf('Dynamic Disambiguation:\n');
modelOut(model,Y(1:1:10,:),'YtoZdyn');

fprintf('Status:\tOK\n');
