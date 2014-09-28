function [expgeno,cover_fraction,concat_id] = networkPropgateEntrezQuick(netobj,genoEntrez,genoEntrezKey,alpha,isInclude,doNorm,name_map)
% Network progate a signal on a network
%
% INPUT: netobj - d x d gene interaction network A in NBS paper
%        genoEntrez - n x p original data matrix,
%        genoEntrezKey - p x 1 corresponding EnsemblID
% OUTPUT: expgeno - n x d matrix 
    zoption.tof = 1e-3;
    zoption.max_iteration = 100;
    zoption.outlevel = 5;      
    
    % Number of genes considered(p)
    if (size(genoEntrez,2) ~= length(genoEntrezKey))
        error('Key missmatch error');
    end
    
%     if (~isfield(netobj,'propVal') || netobj.propVal ~= alpha)
%         error('Error: run generateQuickPropmat()');
%     end
    
% INPUT: netobj - d x d gene interaction network, such as STRING
%        genoEntrez - n x p data matrix(somatic mutation or gene expression)
%        genoEntrezKey - p x 1 Entrez ID for p genes considered
% OUTPUT: igeno - n x d data matrix. Note that only the genes that 
%                      appear in d x d gene interaction matrix are kept and
%                      they are arranged according to Ensembl ID instead of
%                      Entrez ID in the original data matrix. Other elments
%                      among these d genes that don't appear in data matrix
%                      are set to zeros. So if we use a small set of genes
%                      (p<<d), then this matrix will be very sparse.
%         unMappedGenoMat - n x p matrix. Collect genes in the data matrix
%                            that don't appear in gene interactio network.
    [igeno,unMappedGenoMat] = entrezMatToKeyMat(netobj,genoEntrez,genoEntrezKey);
    
    % Percentage of genes covered by the gene interaction network
    cover_fraction = 1 - sum(unMappedGenoMat,2)./sum(genoEntrez,2);
    %fprintf(1,'%f\n',cover_fraction);
    
    itcnt = 1;             
    adjM = netobj.adj_mat_norm;
    adjMInv = full(netobj.adj_mat_norm_val);
  
    
    if (exist('doNorm','var') == 1 && doNorm == 1)     
        fprintf(1,'Normalizing by mutation count\n');  
        igeno = bsxfun(@times,igeno,max(sum(igeno,2),eps).^-1);  
    end
    expgeno = igeno;
    
    fprintf(1,'Fast prop\n');    
    tic;
    %expgeno = mtimesx(sparse(expgeno),adjMInv,'SPEEDOMP') + (1-alpha)*igeno;    
    expgeno = sparse(expgeno)*adjMInv+(1-alpha)*igeno;
    fprintf(1,'%f\n',toc);
    
    epgeno_prev = expgeno;    
    while (itcnt < zoption.max_iteration)
        
        %expgeno = alpha*(mtimesx(expgeno,adjM,'SPEEDOMP')) + (1-alpha)*igeno;        
        expgeno = alpha*(expgeno*adjM)+(1-alpha)*igeno;
        l1residual = norm(epgeno_prev - expgeno,1);
        epgeno_prev = expgeno;
        
        if (l1residual < zoption.tof)
            fprintf(1,'Converged %d\n',itcnt);            
            break;
        end
        itcnt = itcnt + 1;
        % disp(itcnt);
    end
        
    if (exist('isInclude','var') == 1 && isInclude == 1)
        gunmapped = sum(unMappedGenoMat) > 1;
        fprintf(1,'Concat %d unmapped gene values\n',sum(gunmapped));        
        % expgeno = [expgeno mean(expgeno(:))*double(genoEntrez(:,gunmapped)>0) ];        
        % zfact = median(max(expgeno));
        zfact = mean(expgeno(:));
        expgeno = [expgeno zfact*double(genoEntrez(:,gunmapped)>0) ];        
        % expgeno = [logit(expgeno,1) double(genoEntrez(:,gunmapped)>0) ];
        concat_id = genoEntrezKey(gunmapped);
    end
    %expgeno = expgeno';
end