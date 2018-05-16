%--------------------------------------------------------------------------
%Developed by: DIEGO GONZALEZ BECERRA
%Dev. Date: Nov/2017 - XXX/2018
%Function: fGetClosesOBBFaces
%--------------------------------------------------------------------------
function [face1, face2] = fGetClosesOBBFaces(segMainAx, facets1, facets2)
    dMin= 999999;
    if (segMainAx == 1)
        %In this case, by construction of the Bound. Box, we know that the closest faces are 1 and 2 or viceversa
        fI1= reshape(facets1(1, :, :), [4, 3]);
        fI2= reshape(facets1(2, :, :), [4, 3]);

        fNh1= reshape(facets2(1, :, :), [4, 3]);
        fNh2= reshape(facets2(2, :, :), [4, 3]);

        dBts= diag(pdist2(fI1, fNh2));
        if (mean(dBts) < dMin)
            face1= 1;
            face2= 2;
            dMin= mean(dBts);
        end

        dBts= diag(pdist2(fI2, fNh1));
        if (mean(dBts) < dMin)
            face1= 2;
            face2= 1;
        end
    end
    
    if (segMainAx == 2)
        %In this case, by construction of the Bound. Box, we know that the closest faces are 3 and 4 or viceversa
        fI1= reshape(facets1(3, :, :), [4, 3]);
        fI2= reshape(facets1(4, :, :), [4, 3]);

        fNh1= reshape(facets2(3, :, :), [4, 3]);
        fNh2= reshape(facets2(4, :, :), [4, 3]);

        dBts= diag(pdist2(fI1, fNh2));
        if (mean(dBts) < dMin)
            face1= 3;
            face2= 4;
            dMin= mean(dBts);
        end

        dBts= diag(pdist2(fI2, fNh1));
        if (mean(dBts) < dMin)
            face1= 4;
            face2= 3;
        end
    end

    if (segMainAx == 3)
        %In this case, by construction of the Bound. Box, we know that the closest faces are 5 and 6 or viceversa
        fI1= reshape(facets1(5, :, :), [4, 3]);
        fI2= reshape(facets1(6, :, :), [4, 3]);

        fNh1= reshape(facets2(5, :, :), [4, 3]);
        fNh2= reshape(facets2(6, :, :), [4, 3]);

        dBts= diag(pdist2(fI1, fNh2));
        if (mean(dBts) < dMin)
            face1= 5;
            face2= 6;
            dMin= mean(dBts);
        end

        dBts= diag(pdist2(fI2, fNh1));
        if (mean(dBts) < dMin)
            face1= 6;
            face2= 5;
        end
    end
end
