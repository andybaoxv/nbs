function [fullKeyMat,unMappedGenoMat] = entrezMatToKeyMat(netobj,igeno,genoEntrezKey)
% INPUT: netobj - d x d gene interaction network, such as STRING
%        igeno - n x p data matrix(somatic mutation or gene expression)
%        genoEntrezKey - p x 1 Entrez ID for p genes considered
% OUTPUT: fullKeyMat - n x d data matrix. Note that only the genes that 
%                      appear in d x d gene interaction matrix are kept and
%                      they are arranged according to Ensembl ID instead of
%                      Entrez ID in the original data matrix. Other elments
%                      among these d genes that don't appear in data matrix
%                      are set to zeros. So if we use a small set of genes
%                      (p<<d), then this matrix will be very sparse.
%          unMappedGenoMat - n x p matrix. Collect genes in the data matrix
%                            that don't appear in gene interactio network.
    
    fullKeyMat = zeros(size(igeno,1),length(netobj.key));
    unMappedGenoMat = zeros(size(igeno));
    
    % Process each patient one by one, size(igeno,1) = n
    for i = 1:size(igeno,1)
        % nnzPos - nonzeros position of patient i
        % genoEntreKey(nnzPos) - Entrez Key of nonzero elements of 
        % patient i's genes
        nnzPos = logical(igeno(i,:));
        
        % nnzPosIdx - indices of nonzeros gene index
        nnzPosIdx = find(nnzPos);
        
        % Map a 1xp*(p*<p for somatic mutation data, p*=p for gene 
        % expression data) entrez ID to 1 x d vector.
        % Note that there're 12233 genes in the STRING network, when
        % working with small number of genes, most elements of the obtained
        % 1 x d vector are zeros
        %
        % INPUT: netobj - d x d gene interaction network, such as STRING
        %        inEntrzList - 1 x p* entrez ID of a patient's gene vector
        %        inValueList - 1 x p* values for a patient's gene vector
        % OUTPUT: keyGeneVect - 1 x d vector
        %         unmappedGeneVect - 1 x p* vector
        [tt,tunmapped] = netobj.mapEntrezListToKeyVect(genoEntrezKey(nnzPos),igeno(i,nnzPos));
        
        fullKeyMat(i,:) = tt;    
        unMappedGenoMat(i,nnzPosIdx(tunmapped)) = 1;
    end


end