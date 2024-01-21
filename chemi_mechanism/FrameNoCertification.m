%This function is used to determined the frame no obtained by species file
%exists in the *.bonds and *.lammpstrj files
function  tartrajectoryact_new = FrameNoCertification (tartrajectoryact,trajper,choi,outputdatanew)
if mod(tartrajectoryact,trajper(2)) ~= 0 || mod(tartrajectoryact,trajper(3)) ~= 0
    fprintf('\nFrame No %d by species file does not exist in the *.bonds or *.lammpstrj files，maybe caused by different output frequency',tartrajectoryact)
    fprintf('\nRefinement is needed to match the frame No in the three files. Please wait...')
    
    while mod(tartrajectoryact,trajper(2)) ~= 0 || mod(tartrajectoryact,trajper(3)) ~= 0
        if choi == 2
            tartrajectoryact = tartrajectoryact - trajper(1);
        elseif choi == 4
            tartrajectoryact = tartrajectoryact + trajper(1);
        end
        if tartrajectoryact < outputdatanew{2,1} || tartrajectoryact > outputdatanew{size(outputdatanew,1),1}
            error('Frame No %d by refinement operation exceeds the recorded values (in outputdatanew). Please check it!',tartrajectoryact)
        end
    end
    tartrajectoryact_new = tartrajectoryact;
else
    fprintf('\nFrame No by species file exists in the *.bonds or *.lammpstrj files, and no refinement operation is needed')
    tartrajectoryact_new = tartrajectoryact;
end
fprintf('\nFrame No determined by refinement operation is %d\n',tartrajectoryact_new)
end