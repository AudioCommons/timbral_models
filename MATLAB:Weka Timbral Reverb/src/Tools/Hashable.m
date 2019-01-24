classdef (HandleCompatible) Hashable    
    % Interface for classes that provide a method to encode their status as 
    % a hashcode. 
    % Overwrite getHashObjects to specify what elements to use for hash calculation.
    % Can be superclass of handle or value class
    
    %%---------------------------------------------------------------------
    properties 
    end
    
    %%---------------------------------------------------------------------
    methods
        
        function hashcode = getHash( obj )
            hashMembers = obj.getHashObjects();
            hashcode = calcDataHash( hashMembers );
        end
        
        function hashMembers = getHashObjects( obj )
            hashMembers = getNonTransientObjectProps( obj );
        end
        
    end
    
    %%---------------------------------------------------------------------
    methods (Abstract)
        
        
    end
    
    %%---------------------------------------------------------------------
    

    
    %%---------------------------------------------------------------------
    
end

