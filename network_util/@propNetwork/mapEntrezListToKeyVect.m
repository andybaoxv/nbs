function [keyGeneVect,unmappedGeneVect] = mapEntrezListToKeyVect(netobj,inEntrezList,inValueList)
%
% Map entrez list to Key vect
% INPUT: netobj - d x d gene interaction network, such as STRING
%        inEntrzList - 1 x p* entrez ID of a patient's gene vector
%        inValueList - 1 x p* values for a patient's gene vector
% OUTPUT: keyGeneVect - 1 x d vector
%         unmappedGeneVect - 1 x p* vector

    if (exist('inValueList','var') == 0)
        inValueList = ones(size(inEntrezList));
    end
    % Map entrzID to Ensembl ID for a patient's gene vector
    % keyList: p* x 1 Ensembl ID 
    keyList = nanvalues(netobj.entrezIDtoKeyMap,inEntrezList');
    keyGeneVect = zeros(1,length(netobj.key));
    unmappedGeneVect = false(1,length(inEntrezList));
    
    if (iscell(keyList))
        for i = 1:length(keyList)
            if (isempty(keyList{i}))                
                unmappedGeneVect(i) = 1;
                continue;
            end
            % The  Ensembl ID of one gene that is nonzero for a patient
            keyListP = keyList(i);
            % kPos - Position of Ensembl ID in gene interaction network
            kPos = nanvalues(netobj.keyPosMap,keyListP);
            if (isnan(kPos))
                unmappedGeneVect(i) = 1;
                fprintf(1,'Why is this?\n');
                continue;
            end
            % If the Ensembl ID is among the gene in STRING
            if (length(kPos) == 1)
                keyGeneVect(kPos(1)) =  keyGeneVect(kPos(1)) + inValueList(i);
            else
                warn('Ambigious ID');
                keyGeneVect(kPos{1}{1}) =  keyGeneVect(kPos{1}{1}) + inValueList(i);
%                 for zpos = kPos
%                     keyGeneVect(zpos{:}) = keyGeneVect(zpos{:}) + 1/length(kPos);
%                 end
            end            
        end
    else   
        kPos = nanvalues(netobj.keyPosMap,keyList);

        kPosNan = isnan(kPos);
        
        if (all(kPosNan))                
            unmappedGeneVect(:) = 1;
            fprintf(1,'Why is this?\n');            
        else
            keyGeneVect(kPos(~kPosNan)) = 1;
            unmappedGeneVect(kPosNan) = 1;
        end
        
%         for i = 1:length(keyList)
%             if (isempty(keyList(i)) || isnan(keyList(i)) )
%                 unmappedGeneVect(i) = 1;                
%                 continue;
%             end
%             
%             keyListP = keyList(i);
%             
%             kPos = nanvalues(netobj.keyPosMap,keyListP);
%             
%             if (isnan(kPos))                
%                 unmappedGeneVect(i) = 1;
%                 fprintf(1,'Why is this?\n');
%                 continue;
%             end
%             
%             if (length(kPos) == 1)
%                 keyGeneVect(kPos) = keyGeneVect(kPos) + inValueList(i);
%             else
%                 warn('Ambigious ID');
%                 keyGeneVect(kPos(1)) = keyGeneVect(kPos(1)) + inValueList(i);
% %                 for zi = 1:length(kPos)
% %                     zpos = kPos(zi);
% %                     keyGeneVect(zpos) = keyGeneVect(zpos) + 1/length(kPos);
% %                 end
%             end            
%         end
    end
    
    

end 