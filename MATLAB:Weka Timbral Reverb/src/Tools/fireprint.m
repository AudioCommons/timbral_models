function cmap=fireprint(n,varargin)
%FIREPRINT Colormap that increases linearly in lightness (with colors)
%
%   Colormap that increases linearly in lightness (such as a pure black to white
%   map) but incorporates additional colors that help to emphasize the
%   transitions and hence enhance the perception of the data.
%   This colormap is designed to be printer-friendly both for color printers as
%   as well as B&W printers.
%
%	Written by Matthias & Stefan Geissbuehler - matthias.geissbuehler@a3.epfl.ch
%	June 2011
%
%   Credit: The idea of the passages over blue&red stems from ImageJ's LUT 'Fire'
%   Our colormap corrects the color-printout-problems as well as the
%   non-linearity in the fire-colormap which would make it incompatible
%   with a B&W printing.
%
%
%
%   Usage:
%   cmap = fireprint(n)
%
%   All arguments are optional:
%
%   n           The number of elements (256)
%
%   Further on, the following options can be applied
%     'minColor' The absolute minimum value can have a different color
%                ('none'), 'white','black','lightgray', 'darkgray'
%                or any RGB value ex: [0 1 0]
%     'maxColor' The absolute maximum value can have a different color
%     'invert'   (0), 1=invert the whole colormap
%
%
%   Examples:
%     figure; imagesc(peaks(200));
%     colormap(fireprint)
%     colorbar
%
%     figure; imagesc(peaks(200));
%     colormap(fireprint(256,'minColor','black','maxColor',[0 1 0]))
%     colorbar
%
%     figure; imagesc(peaks(200));
%     colormap(fireprint(256,'invert',1,'minColor','darkgray'))
%     colorbar
%

%   Copyright 2011 Matthias Geissbuehler - matthias.geissbuehler@a3.epfl.ch 
%   $Revision: 1.1 $  $Date: 2011/06/14 12:00:00 $
p=inputParser;
p.addParamValue('minColor','none');
p.addParamValue('maxColor','none');
p.addParamValue('invert',0, @(x)x==0 || x==1);

if nargin==1
    p.addRequired('n', @(x)x>0 && mod(x,1)==0);
    p.parse(n);
elseif nargin>1
    p.addRequired('n', @(x)x>0 && mod(x,1)==0);
    p.parse(n, varargin{:});
else
    p.addParamValue('n',256, @(x)x>0 && mod(x,1)==0);
    p.parse();
end
config = p.Results;
n=config.n;


%the ControlPoints
cP(:,1) = [0 0 0];
cP(:,2) = [0 0.6268 0.8926];       %cyan
cP(:,3) = [0.8408 0.1233 0.4041];  %redish-magenta
cP(:,4) = [1 0.9309 0];            %yellow
cP(:,5) = [1 1 1];

cmap = abs(interp1((1:5),double(cP'),linspace(1,5,n)));  % for non-iso points, a normal interpolation gives better results

checkIfAnyAbove1 = 1;
while checkIfAnyAbove1
    % Normalize the shades in order to be compatible with a B&W printout
    cmap(:,1)=cmap(:,1)./sqrt(sum(cmap.^2,2)).*linspace(0,max(sqrt(sum(cmap.^2,2))),n)';
    cmap(:,2)=cmap(:,2)./sqrt(sum(cmap.^2,2)).*linspace(0,max(sqrt(sum(cmap.^2,2))),n)';
    cmap(:,3)=cmap(:,3)./sqrt(sum(cmap.^2,2)).*linspace(0,max(sqrt(sum(cmap.^2,2))),n)';
    cmap(isnan(cmap))=0;
    
    % check if during normalization any value is now bigger than 1
    above1 = cmap>1;
    if sum(above1(:))
        if sum(above1(:,1))  % any R>1 ?
            myIndexes = find(above1(:,1));
            cmap(myIndexes,1) = 1;
            cmap(myIndexes,2) = 0.05 .* (1-cmap(myIndexes,2)) + cmap(myIndexes,2);  % add a little bit to other values
            cmap(myIndexes,3) = 0.05 .* (1-cmap(myIndexes,3)) + cmap(myIndexes,3);  % add a little bit to other values
        end
        if sum(above1(:,2))  % any G>1 ?
            myIndexes = find(above1(:,2));
            cmap(myIndexes,2) = 1;
            cmap(myIndexes,1) = 0.05 .* (1-cmap(myIndexes,1)) + cmap(myIndexes,1);  % add a little bit to other values
            cmap(myIndexes,3) = 0.05 .* (1-cmap(myIndexes,3)) + cmap(myIndexes,3);  % add a little bit to other values
        end
        if sum(above1(:,3))  % any B>1 ?
            myIndexes = find(above1(:,3));
            cmap(myIndexes,3) = 1;
            cmap(myIndexes,1) = 0.05 .* (1-cmap(myIndexes,1)) + cmap(myIndexes,1);  % add a little bit to other values
            cmap(myIndexes,2) = 0.05 .* (1-cmap(myIndexes,2)) + cmap(myIndexes,2);  % add a little bit to other values
        end
        checkIfAnyAbove1 = 1;
    else
        checkIfAnyAbove1 = 0;
    end
end

if config.invert
    cmap = flipud(cmap);
end

if ischar(config.minColor)
    if ~strcmp(config.minColor,'none')
        switch config.minColor
            case 'white'
                cmap(1,:) = [1 1 1];
            case 'black'
                cmap(1,:) = [0 0 0];
            case 'lightgray'
                cmap(1,:) = [0.8 0.8 0.8];
            case 'darkgray'
                cmap(1,:) = [0.2 0.2 0.2];
        end
    end
else
    cmap(1,:) = config.minColor;
end
if ischar(config.maxColor)
    if ~strcmp(config.maxColor,'none')
        switch config.maxColor
            case 'white'
                cmap(end,:) = [1 1 1];
            case 'black'
                cmap(end,:) = [0 0 0];
            case 'lightgray'
                cmap(end,:) = [0.8 0.8 0.8];
            case 'darkgray'
                cmap(end,:) = [0.2 0.2 0.2];
        end
    end
else
    cmap(end,:) = config.maxColor;
end