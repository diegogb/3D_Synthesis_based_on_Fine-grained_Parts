%--------------------------------------------------------------------------
%Developed by: DIEGO GONZALEZ BECERRA
%Dev. Date: Oct/2016
%Funcion: fGetSymmetryChain
%Input:
%   ixSegStart:     first segment to check for symmetric segments
%   allSegments:    structure array containing the segments that were previously computed for all 
%                   the input shapes
%   shapeIndices:   indices of the segments of the Shape for which we are looking for symmetric groups    
%Output:
%   err:            -1 if some error ocurrs; 0, otherwise
%
% Function that obtains a symmetry group or chain, starting from the index
% given by ixSegStart
%
% (See also fSynthesizeWithClosestSegments, fGetClosestSegments)
%%%
%--------------------------------------------------------------------------
function [err, symGroupLabels, symGroupIndices] = fGetSymmetryChain(ixSegStart, allSegments, shapeIndices)
    %indexCloseToUse= 1;
    err= 0;
    try
        %%%
        %TODO: CHECK  THIS CODE: it produces the expected result, but it is not efficient
        %%%
        symGroupLabels= [];
        symGroupIndices= [];
        symm= allSegments(ixSegStart).symLabels;
        if (~isempty(symm))
            if (~exist('shapeIndices', 'var'))
                shapeIndices= fGetSegIndicesPerShape(allSegments, allSegments(ixSegStart).ixShape);
            end
            %Add first the current segment
            symGroupIndices= [symGroupIndices ixSegStart];
            symGroupLabels= [symGroupLabels allSegments(ixSegStart).segLabel];
            
            for s= 1 : length(symm)
                symGroupLabels= [symGroupLabels allSegments(shapeIndices(symm(s))).symLabels];
                symGroupLabels= unique(symGroupLabels);
                symGroupIndices= [symGroupIndices shapeIndices(symm(s))];
                symGroupIndices= unique(symGroupIndices);
                
                symNext= allSegments(shapeIndices(symm(s))).symLabels;
                for j= 1 : length(symNext)
                    symGroupLabels= [symGroupLabels allSegments(shapeIndices(symNext(j))).symLabels];
                    symGroupLabels= unique(symGroupLabels);
                    symGroupIndices= [symGroupIndices shapeIndices(symNext(j))];
                    symGroupIndices= unique(symGroupIndices);
                end                
            end
        end
        
    catch ME
        err= -1;
        sErr= ['Error saving to disk the synthesized shape: (error in fSynthesizeWithClosestSegments))' ME.message];
        errordlg(sErr);
    end
end