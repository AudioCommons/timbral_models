% CIRCVBUFARRAYINTERFACE   This class is a class adding onto the functionality
%  of the circVBuf class. It wraps the circular buffer and provides an array-like 
%  interface for it.
%
%  Data can be accessed using indexing as with arrays, e.g.:
%        buf(1) 
%        buf(:) 
%        buf(end)
%        buf(2:5)
%  The reference (1) of those indices is the oldest data point in the buffer. 'end'
%  refers to the newest data point. Additionally the latest added "chunk" of data
%  can be read by using 'new' as index:
%        buf('new')
%
%  Construct the interface by handing over the circular buffer instance, then use
%  the buffer through the interface object.
%
%
% Ivo Trowitzsch, 04 Novembre 2014
% Technical University of Berlin
% ivot@ni.tu-berlin.de

classdef circVBufArrayInterface < handle

    properties (Access = private)
        cbuf
        nd
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% constructor/destructor
    methods     

        %% ----------------------------------------------------------------
        function obj = circVBufArrayInterface(cbuf)
            obj.cbuf = cbuf;
            obj.nd = ndims( obj.cbuf.dat );
        end
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% public 
    methods (Access=public)

        % -------------- array-like interface ----------------------
        function varargout = subsref( obj, S )
            if (length(S) == 1) && strcmp(S(1).type,'()')
                idxs = S.subs{1,1};
                if isa( idxs, 'char' ) % then it is ':'
                    switch idxs
                        case ':'
                            S.subs{1,1} = obj.cbuf.fst:obj.cbuf.lst;
                        case 'new'
                            S.subs{1,1} = obj.cbuf.new:obj.cbuf.lst;
                        otherwise
                            error( 'indexing not known' );
                    end
                else % it is indexes (1-based)
                    S.subs{1,1} = int64(idxs) + obj.cbuf.fst - int64(1);
                end
                for k = length(S.subs)+1:obj.nd
                    S.subs{1,k} = ':';
                end
                varargout{1:nargout} = obj.cbuf.dat(S.subs{:});
            else
                error( 'Only array-indexing access is allowed through this object!' );
            end
        end
        
        function obj = subsasgn( obj, S, val )
            error( 'Indexed assignment into buffer is not supported at the moment!' );
            if (length(S) == 1) && strcmp(S(1).type,'()')
%                 idxs = S.subs{1,1};
%                 if isa( idxs, 'char' ) % then it is ':'
%                     S.subs{1,1} = obj.fst:obj.lst;
%                 else % it is indexes (1-based)
%                     S.subs{1,1} = int64(idxs) + obj.fst - int64(1);
%                 end
%                 for k = length(S.subs)+1:ndims(obj.dat)
%                     S.subs{1,k} = ':';
%                 end
%                 obj = builtin( 'subsasgn', obj.dat, S, val );
            else
            end
        end
        
        function l = length( obj )
            l = max( 0, obj.cbuf.lst - obj.cbuf.fst + 1 );
        end
            
        function s = size( obj, dim )
            if nargin < 2
                s = size(obj.cbuf.dat);
                s(1) = length( obj );
            else
                s = size(obj.cbuf.dat, dim);
                if dim == 1
                    s(1) = length( obj );
                end
            end
        end
        
        function n = numel( obj )
            n = prod( size( obj ) );
        end
        
        function ind = end( obj, k, n )
            if k == 1
                ind = length( obj );
            else
                ind = builtin( 'end', obj.cbuf.dat, k, ndims(obj.cbuf.dat) );
            end
        end
        
        function ie = isempty( obj )
            ie = (obj.cbuf.lst < obj.cbuf.fst);
        end
        
    end    
    
    
end

