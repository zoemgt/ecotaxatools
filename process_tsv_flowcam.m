%process_tsv_zooscan_new.m allows to assemble the data from the tsv files into a database called "base".  

function [base,m,zoo_groups,Idlistshort]=process_tsv_flowcam

answer2 = questdlg('Flowcam images could be separated in two groups by size with a different sub-sampling (needs to be in two different folder), is it the case?', ...
    'concentration step', ...
    'Yes','No','No');

switch answer2

    case 'Yes'

        f = msgbox('Please select the folder containing the *.tsv files of the first fraction');
        uiwait(f)
        folder1=uigetdir;
        
        prompt = {'Fractionation operated? (in %)'};
        title = 'fraction1';
        dims = [1 35];
        definput = {'100'};
        answerfrac1 = inputdlg(prompt,title,dims,definput);
        answerfrac1=str2num(cell2mat(answerfrac1))/100;
        
        f = msgbox('Please select the folder containing the *.tsv files of the second fraction');
        uiwait(f)
        folder2=uigetdir;
        prompt = {'Fractionation operated? (in %)'};
        title = 'fraction2';
        dims = [1 35];
        definput = {'30'};
        answerfrac2 = inputdlg(prompt,title,dims,definput);
        answerfrac2=str2num(cell2mat(answerfrac2))/100;
        
        cd(folder1)
        A=dir('*.tsv');
        filenames={A.name};
        [n,m]=size(filenames);
        
    case 'No'

        A=dir('*.tsv');
        filenames={A.name};
        [n,m]=size(filenames);       
end

%% create the bases

base=[]; %create the base which gathers the tsv data
Idlist=[]; %create Idlist; Idlist corresponds to the plankton groups list with taxonomic annotation hierarchy included
Idlistshort=[]; %create Idlistshort; Idlistshort corresponds to the plankton groups list without taxonomic annotation hierarchy 

h=waitbar(0,'Load tsv files and create base...');

for i=1:m

waitbar(i / m, h, sprintf('Load tsv files and create base: %d of %d', i, m));

    switch answer2
        
        %combine the 30% and 100% fraction in a single table S with: "d1" for the fraction 100% and "d2" for the fraction 30%
        case 'Yes'
            
            cd(folder1)
            S=readtable(char(filenames(i)),'Filetype','text','ReadVariableNames',1);
            S.acq_id(:)={'d1'}; %100% 
            cd(folder2)
            S2=readtable(char(filenames(i)),'Filetype','text','ReadVariableNames',1);
            S2.acq_id(:)={'d2'}; %30
       
            S.object_link = num2cell(S.object_link);
            S2.object_link = num2cell(S2.object_link); 

            S=[S;S2];   
            
        case 'No'

            S=readtable(char(filenames(i)),'Filetype','text','ReadVariableNames',1);
    end 


sample=unique(S.sample_id);

%% cleaning empty cells
    
[no_use,n]=size(S);
varnames=S.Properties.VariableNames;
    
for j=21:n

    eval(['var_j=S.' char(varnames(j)) ';']);
    if iscellstr(var_j)
    index = (cellfun(@isempty,  var_j) ==1);
    if sum(index)>0
    var_j(index)={'NaN'};
    var_j=(cellfun(@str2num,  var_j));    
    eval(['S.' char(varnames(j)) '=var_j;']);
    end    
    end
        
end
 
%% homogenize the writing of the taxonomic annotation hierarchy

[n,no_use]=size(S);
    
for j=1:n

    f=find(S.object_annotation_hierarchy{j,:}=='-');
    S.object_annotation_hierarchy{j,1}(f)='_';
        
    f=find(S.object_annotation_hierarchy{j,:}=='>');
    S.object_annotation_hierarchy{j,1}(f)='_';

end

%% construction of the base: base
    
    base(i).sampleid=sample; %name of the sample
    
    base(i).ship=unique(S.sample_ship); %name of the ship

    base(i).date=unique(S.object_date); %UTC sampling date
    base(i).time=unique(S.object_time); %UTC sampling time
    base(i).latitude=unique(S.object_lat); %sampling latitude or latitude of the start of the sampling for horizontal towing [decimal degree]
    base(i).longitude=unique(S.object_lon); %sampling longitude or longitude of the start of the sampling for horizontal towing [decimal degree]
    if ismember('object_lat_end', S.Properties.VariableNames) %sampling latitude or latitude of the end of the sampling for horizontal towing [decimal degree]
    base(i).latitude_end = unique(S.object_lat_end);
    end 
    if ismember('object_lon_end', S.Properties.VariableNames) %sampling longitude or longitude of the end of the sampling for horizontal towing [decimal degree]
    base(i).longitude_end=unique(S.object_lon_end); 
    end
    base(i).zmax=unique(S.object_depth_min); %maximum depth of the net when deployed
    base(i).zmin=unique(S.object_depth_min); %minimum depth of the net when deployed

%% getting pixel size in micrometer and converting in mm

if isnumeric(unique(S.process_pixel)) %dimension of the side of a pixel in the scanned image [um]
pixelsize=(unique(S.process_pixel))/1000;  %[um], convert in mm/pixel
else
pixelsize=str2num((cell2mat(unique(S.process_pixel))))/1000;  %[um], convert in mm/pixel
end

base(i).pixelsize=pixelsize;

%% construction of the base

    % raw_image_total
    if isnumeric(unique(S.acq_raw_image_total))  %number of raw images recorded by visualspreadsheet during the process [#]
        base(i).raw_image_total=unique(S.acq_raw_image_total); 
    else
        base(i).raw_image_total=str2num(cell2mat(unique(S.acq_raw_image_total)));
    end
    
    % process_nb_images
    if isnumeric(unique(S.process_nb_images))
        base(i).process_nb_images=unique(S.process_nb_images); %number of images processed by Zooprocess
    else
        base(i).process_nb_images=str2num(cell2mat(unique(S.process_nb_images)));
    end
    
    % volconc in m3 
    if ismember('sample_conc_vol_ml', S.Properties.VariableNames) %different version before and after reprocess FlowCam 2024

    if isnumeric(unique(S.sample_conc_vol_ml)) %concentrated or diluted water volume (from sample_comment_or_volume) [mL]
        base(i).volconc=unique(S.sample_conc_vol_ml)/1000000; %convert in m3
    else
        base(i).volconc=str2num(cell2mat(unique(S.sample_conc_vol_ml)))/1000000;
    end
    
    else
 
    if isnumeric(unique(S.sample_volconc)) %concentrated or diluted water volume (from sample_comment_or_volume) [mL]
        base(i).volconc=unique(S.sample_volconc)/1000000; %convert in m3
    else
        base(i).volconc=str2num(cell2mat(unique(S.sample_volconc)))/1000000;
    end

    end
    
    % sample_initial_col_vol in m3
    if ismember('sample_initial_col_vol_m3', S.Properties.VariableNames) %different version before and after reprocess FlowCam 2024

        test=unique(S.sample_initial_col_vol_m3); %initial collected volume, (if nets : sum of the nets) [m3]

        if isnumeric(test)
    
            if length(test)==1
            base(i).sample_initial_col_vol=unique(S.sample_initial_col_vol_m3);
            else
            base(i).sample_initial_col_vol=unique(S.sample_initial_col_vol_m3);
            end
        
            if isempty(base(i).sample_initial_col_vol)==1
            base(i).sample_initial_col_vol=NaN;
            end
        
            if size(base(i).sample_initial_col_vol,2)>1
            base(i).sample_initial_col_vol=base(i).sample_initial_col_vol(1);
            end
       
        else
    
            if length(test)==1
            base(i).sample_initial_col_vol=str2num(cell2mat(unique(S.sample_initial_col_vol)));
            else
            base(i).sample_initial_col_vol=str2num(cell2mat(unique(S.sample_initial_col_vol)));
            end
        
            if isempty(base(i).sample_initial_col_vol)==1
            base(i).sample_initial_col_vol=NaN;
            end
        
            if size(base(i).sample_initial_col_vol,2)>1
            base(i).sample_initial_col_vol=base(i).sample_initial_col_vol(1);
            end
    
        end
    
    else 

        test=unique(S.sample_comment_or_volume); %initial collected volume, (if nets : sum of the nets) [m3]

        if isnumeric(test)
    
            if length(test)==1
            base(i).sample_initial_col_vol=unique(S.sample_comment_or_volume)/1000000; %convert in m3
            else
            base(i).sample_initial_col_vol=unique(S.sample_comment)/1000000; %convert in m3
            end
        
            if isempty(base(i).sample_initial_col_vol)==1
            base(i).sample_initial_col_vol=NaN;
            end
        
            if size(base(i).sample_initial_col_vol,2)>1
            base(i).sample_initial_col_vol=base(i).sample_initial_col_vol(1);
            end
       
        else
    
            if length(test)==1
            base(i).sample_initial_col_vol=str2num(cell2mat(unique(S.sample_comment_or_volume)))/1000000; %convert in m3
            else
            base(i).sample_initial_col_vol=str2num(cell2mat(unique(S.sample_comment)))/1000000; %convert in m3
            end
        
            if isempty(base(i).sample_initial_col_vol)==1
            base(i).sample_initial_col_vol=1;
            end
        
            if size(base(i).sample_initial_col_vol,2)>1
            base(i).sample_initial_col_vol=base(i).sample_initial_col_vol(1);
            end
    
        end

    end
    
    % fluid_vol_imaged in m3
    if isnumeric(unique(S.acq_fluid_volume_imaged))

    base(i).fluid_vol_imaged=unique(S.acq_fluid_volume_imaged)/1000000;  %flowcam total images volume [mL] and convert in m3

    else 

        temp=unique(S.acq_fluid_volume_imaged);
        f=strfind(temp,'_'); %supressing the "_ml"
        temp=char(temp);
        temp2=temp(1:cell2mat(f)-1);   
        base(i).fluid_vol_imaged=str2num(temp2)/1000000;

    end

% compute the volume_imaged_processed in m3 
base(i).volume_imaged_processed = (base(i).fluid_vol_imaged/ base(i).raw_image_total)* min(base(i).raw_image_total,base(i).process_nb_images); %compute the final volume imaged processed by the FlowCam [m3] 

%% FracId (d1: fraction 100% and d2: fraction 30%)
%add d1 (100%) and d2 (30%) fraction 

base(i).fracids=unique(S.acq_id); %d1 and d2 fraction: 100% and 30% fraction
[nfrac, no_use]=size(base(i).fracids);

for fracnb=1:nfrac 

    I=strcmp(S.acq_id,base(i).fracids(fracnb));

    eval(['base(i).d' num2str(fracnb) '.fracmin=unique(S.acq_min_esd);']); %minimum fraction size (µm)
    eval(['base(i).d' num2str(fracnb) '.fracsup=unique(S.acq_max_esd);']); %maximum fraction size (µm)
    eval(['base(i).d' num2str(fracnb) '.object_annotation_hierarchy=unique(S.object_annotation_hierarchy(I));']); %taxonomical annotation hierarchy
    eval(['base(i).d' num2str(fracnb) '.anotStatus=S.object_annotation_status(I);']); 
    
    %operation to create Idlistshort:
    [Idlist, ia, ic]=unique([Idlist; S.object_annotation_hierarchy(I)]);
    Idlistshort=[Idlistshort; S.object_annotation_category(I)];
    Idlistshort=Idlistshort(ia);
    
    eval(['base(i).d' num2str(fracnb) '.object_annotation_hierarchy=S.object_annotation_hierarchy(I);']); 
    I2=strfind(S.object_annotation_status(I),'validated');
    index = (cellfun(@isempty,  I2) ==0);
    eval(['base(i).d' num2str(fracnb) '.percentValidated=100*sum(index)/length(index);']); %percentage of validated vignettes in the sample

    if iscell(S.object_major)==0

    eval(['base(i).d' num2str(fracnb) '.major=S.object_major(I)*pixelsize;']); %length in pixels of the major axis of the best ellipsoid approximation (length of the organism)
    eval(['base(i).d' num2str(fracnb) '.minor=S.object_minor(I)*pixelsize;']); %length in pixels of the minor axis of the bestellipsoid approximation (width of the organism)
    eval(['base(i).d' num2str(fracnb) '.area_exc=S.object_area_exc(I)*(pixelsize^2);']); %ratio between the distance between the foci of the ellipse and the length of its major axis (elongation of the organism)
    eval(['base(i).d' num2str(fracnb) '.area=S.object_area(I)*(pixelsize^2);']); %number of pixels in the region (size of the organism)
    eval(['base(i).d' num2str(fracnb) '.area_origin=S.object_area(I);']); % ?? 
    eval(['base(i).d' num2str(fracnb) '.ESD=2*(((S.object_area(I))*(pixelsize^2)/pi).^0.5);']); %ESD: Equivalent Spherical Diameter (µm)
    eval(['base(i).d' num2str(fracnb) '.perimferet=S.object_feret(I)*pixelsize;']); %distance around the border of the region (size of the organism)
    
    else
    
    eval(['base(i).d' num2str(fracnb) '.major=cellfun(@str2num,S.object_major(I))*pixelsize;']); %length in pixels of the major axis of the best ellipsoid approximation (length of the organism)
    eval(['base(i).d' num2str(fracnb) '.minor=cellfun(@str2num,S.object_minor(I))*pixelsize;']); %length in pixels of the minor axis of the bestellipsoid approximation (width of the organism)
    eval(['base(i).d' num2str(fracnb) '.area_exc=cellfun(@str2num,S.object_area_exc(I))*(pixelsize^2);']); %ratio between the distance between the foci of the ellipse and the length of its major axis (elongation of the organism)
    eval(['base(i).d' num2str(fracnb) '.ESD=2*((cellfun(@str2num,S.object_area(I))*(pixelsize^2)/pi).^0.5);']); %ESD: Equivalent Spherical Diameter (µm)
    eval(['base(i).d' num2str(fracnb) '.area=cellfun(@str2num,S.object_area(I))*(pixelsize^2);']); %number of pixels in the region (size of the organism)
    eval(['base(i).d' num2str(fracnb) '.area_origin=cellfun(@str2num,S.object_area(I));']); % ??
    eval(['base(i).d' num2str(fracnb) '.perimferet=cellfun(@str2num,S.object_feret(I))*pixelsize;']); %distance around the border of the region (size of the organism)
    
    end
        
switch answer2

     case 'Yes'

     if fracnb==1
     falsemotoda=answerfrac1;
     elseif fracnb==2
     falsemotoda=answerfrac2;
     end
                        
     eval(['base(i).d' num2str(fracnb) '.conver=  base(i).volconc./(base(i).volume_imaged_processed.*base(i).sample_initial_col_vol.*falsemotoda);']); %[m3]
     
     case 'No'

     eval(['base(i).d' num2str(fracnb) '.conver=  base(i).volconc./(base(i).volume_imaged_processed.*base(i).sample_initial_col_vol);']); %[m3]

    
end

end 

end 

close(h);

[n,m]=size(base);
zoo_groups=Idlist;
