classdef circVBuf < handle
    %circVBuf class defines a circular double buffered vector buffer
    %   -- newest vector has the highest index
    %   -- oldest vector has the lowest index
    %   -- append adds new vectors at the 'right side' of the buffer
    %   -- current buffer always accessable as a complete single subarray buf(fst:lst)
    %
    %   special advantages:
    %   -- any matrix operations at any time possible !
    %   -- time required for appending is independent on buffer status !
    %   -- index 'new' indicates the start of last appended vectors    
    % 
    % View on buffer while appending new vetors
    % -----------------------------------------
    % 1.  |aaaaaaaaa.....AAAAAAAAA.....| append 9 vectors => fill first time
    %                    f-------l
    %                    n
    % 2.  |aaaaaaaaaaaaaaAAAAAAAAAAAAAA| append 5 vectors => almost cycle, no old vectors (a/A) overwritten yet
    %                    f---------|--l                      l-f stays constant from now on until buffer cleared
    %                              n   
    % 3.  |BaaaaaaaaaaaaabAAAAAAAAAAAAA| append 1 vector  => switch to left copy
    %       f------------l                                   first time a vector gets overwritten (b/B overwrites a/A)
    %                    n
    % 4.  |BBBBBBaaaaaaaabbbbbbAAAAAAAA| append 5 vectors => valid buffer moves from left to right
    %            f--------|---l
    %                     n 
    % 5.  |BBBBBBBBBBBBBBbbbbbbbbbbbbbb| append 8 vectors => almost cycle
    %                    f-----|------l                      vectors from first fill (a/A) completely overwritten.
    %                          n
    % 6.  |cBBBBBBBBBBBBBCbbbbbbbbbbbbb| append 1 vector  => c/C overwrites b/B
    %       f------------l
    %                    n
    

    %
    % Spezials :
    % ---------------------------------------------
    % a. cleared/empty buffer (if buffer is empty: l<f, l<n, n=f)
    %     |............................| 
    %                   lf 
    %                    n
    % b. no new vectors added (nonew() called) => l<n   
    %     |cBBBBBBBBBBBBBCbbbbbbbbbbbbb| 
    %       f            l
    %                     n (<-- might exceed buffer dimention by 1 !!)
    %
    % f - first index (obj.fst)
    % l - last index (obj.lst)
    % n - new index (obj.new)
    %
    % example 1 - loop over new vectors of a circVBuf (slow because of copy):
    %     for ix=circVBufObj.VBuf.new:circVBuf.VBuf.lst
    %        vec(:) = circVBuf.VBuf.raw(:,bId);
    %     end
    %
    % example 2 - direct array operation on new vectors in the buffer (no copy => fast)
    %     new = circVBuf.VBuf.new;
    %     lst = circVBuf.VBuf.lst;
    %     mean = mean(circVBuf.VBuf.raw(3:7,new:lst)); 
    %     
    %
    % TODO: rethink "unique life time id function": |123456789.....123456789.....| => fst always 1
    %       get() function returns part of the buffer as copy
    %       pos(rIdx) translates from raw index to position
    %       idx(pos) translates from position to raw index
    
    % $Author: jhgoebbert $    $Date: 2014/06/29 20:00:00 $    $Revision: 0.9.3 $    
    % Copyright 2014 Jens Henrik Goebbert <jens.henrik.goebbert@rwth-aachen.de>
    %#codegen
    
    % changes by Ivo Trowitzsch, 2014/09
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% properties

    %% --------------------------------------------------------------------     
    properties (SetAccess = private, GetAccess = public)
        dat                    % buffer (initialized in constructor)
        matSz@int64 = int64(0) % size of matrix to store (only change in constructor)
        bufSz@int64 = int64(0) % max number of vectors to store (only change in constructor)
        
        fst@int64   = int64(nan) % first index == position of oldest/first value in circular buffer
        new@int64   = int64(nan) %   new index == position of first new value added in last append()
        lst@int64   = int64(nan) %  last index == position of newest/last value in circular buffer
        
        newCnt@int64= int64(0)   % number of new values added lately (last append call).
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% constructor/destructor
    methods     

        %% ----------------------------------------------------------------
        function obj = circVBuf(bufSize,matSize)
%            assert(isa(bufSize,'int64'))                
%            assert(isa(matSize,'int64'))
                
            obj.setup(bufSize, matSize);
        end
        
        %% ----------------------------------------------------------------        
        function delete(obj)
            obj.dat = []; % FIXME: probably not required ?
        end
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% public 
    methods (Access=public)

        %% ----------------------------------------------------------------
        function setup(obj,bufSize,matSize)
%            assert(isa(bufSize,'int64'))                
%            assert(isa(matSize,'int64'))
            
            % buffer initialized once here
            obj.bufSz = int64(bufSize); % fixed values         
            obj.matSz = int64(matSize); % fixed values
            
            obj.dat  = nan([bufSize*2, matSize], 'double');
              
            obj.clear();
        end
        
        %% ----------------------------------------------------------------        
        function clear(obj)
            obj.fst = obj.bufSz+1;
            obj.lst = obj.bufSz;
            
            obj.new = obj.fst;
            obj.newCnt = int64(0);               
        end        
        
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % ----------------------------------------------------------------
        function cpSz = append(obj,vec)
            %assert(isa(vec,'double')) % disabled because it consumes time   
            
            % preload values == increase performance !?
            f = obj.fst;
            l = obj.lst;          
            
            % preload values == increase performance !?
            vSz  = size(vec,1);
            bSz  = obj.bufSz;
            
            % calc number of vectors to add to buffer and start position in vec            
            cpSz  = min(vSz, bSz);         % do not copy more vectors than buffer size
            cpSz1 = min(cpSz, (bSz*2 -l)); % no. vectors added on the right side (beginning with pos lst) 
            cpSz2 = cpSz -cpSz1;           % no. vectors added on left side (beginning with pos 1)
            
            vSt = max(1, vSz-cpSz+1);      % start position in input vector array (we might have to skip values if vSz>bSz)
 
            % add data after lst
            obj.dat(l+1    :l+cpSz1    ,:) = vec(vSt:vSt+cpSz1-1,:);
            obj.dat(l+1-bSz:l+cpSz1-bSz,:) = vec(vSt:vSt+cpSz1-1,:);
            
            % cpSz2: number of vectors to add at buffer begin
            if(cpSz2 == 0)
                % add |bbbbbbb|: cpSz1==7, cpSz2==0
                % |AAAaaaaaaaaaaAAAAAAA|  -->  |AAABBBBBBBaaabbbbbbb|
                %     f--------l                          f--------l 
                %     4        13                         11       20
                obj.fst = min(bSz+1, f+cpSz1); % until buffer is completly filled the first time min() is required
                obj.lst = l +cpSz1;
                obj.new = l + cpSz1 - cpSz + 1;
            else % called only on buffer cycle (performance irrelevant)
                % add |bbbbbbbb|: cpSz1==7, cpSz2==2
                % |AAAaaaaaaaaaaAAAAAAA|  -->  |bbABBBBBBBBBabbbbbbb|
                %     f--------l                  f--------l 
                %     4        13                 3        12               
                obj.dat(    1:cpSz2,    :) = vec(vSt+cpSz1:vSt+cpSz-1,:); % copy bb
                obj.dat(bSz+1:cpSz2+bSz,:) = vec(vSt+cpSz1:vSt+cpSz-1,:); % copy BB
                
                obj.fst = cpSz2 +1;
                obj.lst = cpSz2 +bSz;            
                obj.new = cpSz2 + bSz - cpSz + 1;
            end
      
            obj.newCnt = cpSz;
        end

        %% ----------------------------------------------------------------     
        function nonew(obj)          
            obj.new    = obj.lst+1;
            obj.newCnt = int64(0);
        end
                
        %% ----------------------------------------------------------------        
        function id = lifetimeId(obj,idx)
            % Calc life-time id of a buffer value independent of fst and lst.
            % This id survives the switch of double buffering.
            id = mod(idx, obj.bufSz);
        end 
        
    end    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% static
    methods(Static)
        
        %% ---------------------------------------------------------------- 
        function success = test(appendType)
            % TEST Function to test class.
            success = false;
            
            if (nargin < 1) || isempty(appendType)
               appendType = 0;
            end
            
            % setup                        
            bufferSz   = 1000;
            vectorLen  = 7;
            stepSz     = 10;
            steps      = 100;
            
            % create/setup test object
            testObj = circVBuf(int64(bufferSz),int64(vectorLen),appendType);
                        
            % fill circular buffer with steps*stepSz vectors
            vecs = zeros(stepSz,vectorLen,'double');
            %tic
            for i=0:steps-1 % no. steps
                for j=1:stepSz
                    vecs(j,:) = (i*stepSz)+j;
                end
                testObj.append(vecs);
            end
            %toc
            
            % check last bufSz vectors 
            cnt = steps*stepSz;
            for i=testObj.lst:-1:testObj.fst
               vec = testObj.dat(i,:);
               assert( mean(vec(:)) == cnt, 'TEST FAILED: mean(..) ~= cnt');
               cnt = cnt -1;
            end
            
            success = true;
        end
        
    end
    
end

