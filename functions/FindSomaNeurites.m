function [Neurites, newSkel, axon, newEndPoints] = FindSomaNeurites(inputSkl, cBody, prevEndPoints, MIN_LEN)     %% ps = [endPoints; branchPoints]
%FindSomaNeurites finds all longest paths originating from soma and returns the longest as the axon


%% find starting points on Soma
% remove cell body
noBody = inputSkl  & ~imerode(cBody,  strel('disk',1)); % imerode to make sure dendrites get a start point

% more stable via bodyBorderPoints, it adds more computation but fails less
NeuriteStartPoints = find(bwmorph(cBody, 'endpoints'));

% remove start points outside the skeleton
NeuriteStartPoints = intersect(NeuriteStartPoints, find(noBody));


Neurites = cell(size(NeuriteStartPoints));
newEndPoints = zeros(size(NeuriteStartPoints));


%% for each starting point find longest path in skeleton
disp('finding neurites');
    for k = flip(1:nnz(NeuriteStartPoints)) % going backwards (in case one is removed)

        [Neurites{k}, newEndPoints(k), noBody] = findLongestConnected(NeuriteStartPoints(k), noBody);

        len = nnz(Neurites{k});
        
        % if length of Neurites{k} < MIN_LEN) -> remove path (& endpoint)
        if  len < MIN_LEN
            noBody = noBody & ~Neurites{k};
            Neurites{k} = [];
            % remove endpoints
            prevEndPoints = setdiff(prevEndPoints, newEndPoints(k));
            newEndPoints(k) = [];          
        end
         
        
    
    end
    
    %% Test for overlap and fix if there is
    newSkel = noBody;
    disp('refining neurites');
    [newSkel, Neurites, newEndPoints ] = fixOverlap(newSkel, NeuriteStartPoints, Neurites, newEndPoints, MIN_LEN );
    
        
    %% classify as axon & dendrites
    
    axon_len = 0;
    axon =  cell(1);
    
    for k = flip(1:size(Neurites)) % going backwards (in case one is removed)
        len = nnz(Neurites{k});
        if len > axon_len    % find the longest one
          axon{1} = Neurites{k};
          axon_len = len;
          axonIndex = k;
        end
        % also remove already classified parts from skeleton now
        if ~isempty(Neurites{k})
            newSkel = newSkel & ~Neurites{k};
        end
    end
    
    
    
    % remove axon from Neurites{} and treat rest as Dendrites
    Neurites{axonIndex} = [];
    Neurites = Neurites(~cellfun('isempty',Neurites));
    
    newEndPoints = setdiff(prevEndPoints, newEndPoints);

    

end