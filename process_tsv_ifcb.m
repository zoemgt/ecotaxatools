%process_tsv_ifcb.m allows to assemble the data from the tsv files into a database called "base".  

function [base,m,zoo_groups,Idlistshort]=process_tsv_ifcb

%% create the bases

A=dir('*.tsv');
filenames={A.name};
[n,m]=size(filenames);

base=[]; %create the base which gathers the tsv data
Idlist=[]; %create Idlist; Idlist corresponds to the plankton groups list with taxonomic annotation hierarchy included 
Idlistshort=[]; %create Idlistshort; Idlistshort corresponds to the plankton groups list without taxonomic annotation hierarchy 

h=waitbar(0,'Load tsv files and create base...');

for i=1:m
   
S=readtable(char(filenames(i)),'Filetype','text','ReadVariableNames',1);
sample=unique(S.sample_id); 

waitbar(i / m, h, sprintf('Load tsv files and create base: %d of %d', i, m));

%% homogenize the writing of the taxonomic annotation hierarchy
    
[n,no_use]=size(S);
    
    for j=1:n;
        f=find(S.object_annotation_hierarchy{j,:}=='-');
        S.object_annotation_hierarchy{j,1}(f)='_';
        
        f=find(S.object_annotation_hierarchy{j,:}=='>');
        S.object_annotation_hierarchy{j,1}(f)='_';
    end   
      
%% construction of the base
    
    base(i).sampleid=sample;

    base(i).date=unique(S.object_date); %sampling date
    base(i).time=unique(S.object_time); %sampling time 
    base(i).latitude=double(unique(S.object_lat)); %sampling latitude [decimal degree]
    base(i).longitude=double(unique(S.object_lon)); %sampling longitude [decimal degree]
    if ismember('object_lat_end', S.Properties.VariableNames) %sampling latitude or latitude of the end of the sampling for horizontal towing [decimal degree]
    base(i).latitude_end = unique(S.object_lat_end);
    end 
    if ismember('object_lon_end', S.Properties.VariableNames) %sampling longitude or longitude of the end of the sampling for horizontal towing [decimal degree]
    base(i).longitude_end=unique(S.object_lon_end); 
    end

    base(i).vol=unique(S.acq_volume_sampled)/1000000; %[mL], convert in m3
       
    %% getting pixel size in micrometer and converting in mm

    pixelsize=(1/((unique(S.acq_resolution_pixel_per_micron))))/1000;  %[um], convert in mm/pixel

%% FracId (d1, d2, ... or tot) 

base(i).fracids=unique(S.acq_id); %see if there is size fraction
[nfrac, no_use]=size(base(i).fracids);

    for fracnb=1:nfrac 
        %nfrac=1, no size fraction (tot)
        %nfrac=2, size fraction was used 

        I=strcmp(S.acq_id,base(i).fracids(fracnb));

        eval(['base(i).d' num2str(fracnb) '.object_annotation_hierarchy=unique(S.object_annotation_hierarchy(I));']); %taxonomical annotation hierarchy
        
        %operation to create Idlistshort: 
        [Idlist, ia, ic]=unique([Idlist; S.object_annotation_hierarchy(I)]); 
        Idlistshort=[Idlistshort; S.object_annotation_category(I)];
        Idlistshort=Idlistshort(ia);
        
        eval(['base(i).d' num2str(fracnb) '.object_annotation_hierarchy=S.object_annotation_hierarchy(I);']);
        I2=strfind(S.object_annotation_status(I),'validated');
        index = (cellfun(@isempty,  I2) ==0);
        eval(['base(i).d' num2str(fracnb) '.percentValidated=100*sum(index)/length(index);']); %percentage of validated vignettes in the sample

        eval(['base(i).d' num2str(fracnb) '.major=S.object_major_axis_length(I)*pixelsize;']); %length in pixels of the major axis of the best ellipsoid approximation (length of the organism)
        eval(['base(i).d' num2str(fracnb) '.minor=S.object_minor_axis_length(I)*pixelsize;']); %length in pixels of the minor axis of the bestellipsoid approximation (width of the organism) 
        eval(['base(i).d' num2str(fracnb) '.area=S.object_surface_area(I)*(pixelsize^2);']); %number of pixels in the region (size of the organism)
        eval(['base(i).d' num2str(fracnb) '.ESD=2*(((S.object_surface_area(I)*(pixelsize^2))/pi).^0.5);']); %ESD: Equivalent Spherical Diameter (µm)
        
        eval(['base(i).d' num2str(fracnb) '.summedbiovolume=S.object_summed_biovolume(I)*(pixelsize^3);'])
        eval(['base(i).d' num2str(fracnb) '.summedarea=S.object_summed_surface_area	(I)*(pixelsize^2);'])

        eval(['base(i).d' num2str(fracnb) '.conver=  1./base(i).vol;']) %[m3]

    end

end
      
close(h);

[n,m]=size(base);
zoo_groups=Idlist;  

end