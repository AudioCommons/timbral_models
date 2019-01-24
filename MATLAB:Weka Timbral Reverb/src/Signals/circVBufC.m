classdef circVBufC < handle
    %circVBufC class defines a circular complex single vector ring buffer
    %   for details see circVBuf.m
    
    %% properties
    
    properties (SetAccess = private, GetAccess = public)
        dat                    % buffer (initialized in constructor)
        matSz@int64 = int64(0) % size of matrix to store (only change in constructor)
        bufSz@int64 = int64(0) % max number of vectors to store (only change in constructor)
        
        fst@int64   = int64(nan) % first index == position of oldest/first value in circular buffer
        new@int64   = int64(nan) % new index == position of first new value added in last append()
        lst@int64   = int64(nan) % last index == position of newest/last value in circular buffer
        
        newCnt@int64= int64(0)   % number of new values added lately (last append call).
    end
    
    
    %% methods
    
    methods   
        
        function obj = circVBufC(bufSize,matSize)                
            obj.setup(bufSize, matSize);
        end
        
        function delete(obj)
            obj.dat = [];
        end
        
    end
    
    %% public 
    
    methods (Access=public)

        function setup(obj,bufSize,matSize)            
            
            % buffer initialized once here
            obj.bufSz = int64(bufSize); % fixed values         
            obj.matSz = int64(matSize); % fixed values
            % CHANGE TO SINGLE PRECISION
            obj.dat = nan([bufSize*2, matSize], 'single');
            % CHANGE TO COMPLEX
            obj.dat = obj.dat * 1i;
              
            obj.clear();
        end
        
        function clear(obj)
            obj.fst = obj.bufSz+1;
            obj.lst = obj.bufSz;
            
            obj.new = obj.fst;
            obj.newCnt = int64(0);               
        end        
        
        function cpSz = append(obj,vec)
            
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

        function nonew(obj)          
            obj.new    = obj.lst+1;
            obj.newCnt = int64(0);
        end
                
        function id = lifetimeId(obj,idx)
            % Calc life-time id of a buffer value independent of fst and lst.
            % This id survives the switch of double buffering.
            id = mod(idx, obj.bufSz);
        end 
        
    end    
    
    %% static
    
    methods(Static)
        
        function success = test()
            % TEST Function to test class.
            success = false;
            
            % setup                        
            bufferSz   = 1000;
            vectorLen  = 7;
            stepSz     = 10;
            steps      = 100;
            
            % create/setup test object
            testObj = circVBufC(int64(bufferSz),int64(vectorLen));
                        
            % fill circular buffer with steps*stepSz vectors
            vecs_r = zeros(stepSz,vectorLen,'single');
            vecs_i = zeros(stepSz,vectorLen,'single');
            vecs = complex(vecs_r, vecs_i);
            
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

