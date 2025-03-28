%% EcoTaxa Tool Box pipeline 2025 
% 
% We present a pipeline, EcoTaxa Tools Box on Matlab (available on the GitHub platform), designed to efficiently process and explore 
% the vast amounts of annotated images on EcoTaxa, coming from various quantitative imaging instrument: UVP, IFCB, Flowcam, PlanktoScope and Zooscan.
% EcoTaxa Tool Box pipeline helps you to process raw results originating from quantitative
% imaging devices and originating from Ecotaxa (https://ecotaxa.obs-vlfr.fr/)
% ecotaxa tool includes the initial process of raw *.tsv files available from
% the export function in ecotaxa (!!!! Warning !!!! , please use the option "one file/sample")

% The resulting structured base includes the following variables as final results:
%     Ab  abundance per groups (ind.m-3)
%     Bv biovolume per groups (mm3.m-3)
%     Ybv_Plain_Area_BV_spectra NBSSbiovolume per groups per size class (mm3.mm-3.m-3)
%     Ybv_Riddled_Area_BV_spectra NBSSbiovolume per groups per size class (mm3.mm-3.m-3)
%     Ybv_Ellipsoid_BV_spectra NBSSbiovolume per groups per size class (mm3.mm-3.m-3)
%     X  is the Middle of each biovolume size class (caution it is in log here) (log(mm3)
%     X1 is the amplitude of each biovolume size class - used for de-normalizing NBSS (mm3)
%     ESD vector is the conversion of X in ESD (µm this time)

% Results are stored in a multilevel layer:
% d1 => dx fifferent fractions containing raw measurments of data
% tot => assemblage of the different fractions
% regroupped => additional groups (tot, living/notliving , pft groups or trophic groups) are added from additionning the previous taxa

% This pipeline is made up of several codes, one of which is the main one and is able to call up the other codes, this one EcoTaxa_Tool_box.m 
% The other codes are available in the Box folder and include a description of each function. 
% The first part consists of creating a matlab analysis database containing key information from the tsvs exported from EcoTaxa; abundances, biovolumes, 
% size spectra by different levels of taxonomic grouping (ecotaxa_tools.m). 
% The 2nd part of the code provides an exploratory analysis of these different variables, such as maps and histograms describing the plankton communities. 
% Some initial avenues for multivariate analysis are also proposed at the end of the code (normalisation, PCA, clustering). 

%%
clear all
close all
warning('off');

% add lib to the path automatically 
currentScriptPath = mfilename('fullpath');
[currentFolder, ~, ~] = fileparts(currentScriptPath);
libPath = fullfile(currentFolder, 'lib');
addpath(genpath(libPath));

%% Data to load 

list = {'EcoTaxa export (tsv files) to process','Matlab base.m (processed)'};
    
[datatoload,tf] = listdlg('PromptString','What data should be uploaded?','ListString',list,'SelectionMode','single','ListSize', [400, 200]);

if datatoload==1; 

ecotaxa_tools; %Process before the tsv files into base.mat

list2 = {'Stop the pipeline after the process','Continue to ecological analysis'};
[indx,tf] = listdlg('PromptString','What data should be uploaded?','ListString',list2,'SelectionMode','single','ListSize', [400, 200]);

if indx==1 %Stop the pipeline after the process
    return; 
elseif indx==2; %Continue to ecological analysis
clear all 
f = msgbox('Select the matlab base *.mat (already processed)');
uiwait(f)
[file,path]=uigetfile('*.mat');
load([path file]);
end

elseif datatoload==2;
f = msgbox('Select the matlab base *.mat (already processed)');
uiwait(f)
[file,path]=uigetfile('*.mat');
load([path file]);
end

[n,m]=size(base);

%% latitude and longitude of the samples

%extract longitude and latitude from the dataset
lat=[base.latitude];
lon=[base.longitude];

%extract sample lables from the dataset
sampleid=[]; 
for i=1:numel(base)  
    sampleid=[sampleid;base(i).sampleid];
end

%% Export and Save graphic ? 

list2 = {'Yes', 'No'};  
[save,tf] = listdlg('PromptString','Would you like to export all the graphs (.pdf)?','ListString',list2,'SelectionMode','single','ListSize', [400, 200]);

%% Map Station

%display the figures in full screen
scrsz = get(0,'ScreenSize');

%color of the land in the map
landcolor=[1 0.949019610881805 0.866666674613953];

%use min/max latitude and longitude of the samples to define the map

latmax=max(lat);
latmin=min(lat); 

lonmax=max(lon);
lonmin=min(lon);

lat_extension = 0.30 * (latmax - latmin);
lon_extension = 0.30 * (lonmax - lonmin);

%case of one localisation 
if (latmax - latmin) == 0
    lat_extension = 0.5; 
end
if (lonmax - lonmin) == 0
    lon_extension = 0.5; 
end

latmax = latmax + lat_extension;
latmin = latmin - lat_extension;
lonmax = lonmax + lon_extension;
lonmin = lonmin - lon_extension;

%map1: station map
figure('numbertitle','off','Position',[10 50 scrsz(3) scrsz(4)-150],'Color','white');
m_proj('mercator','lon',[lonmin lonmax],'lat',[latmin latmax]);
m_coast('patch',landcolor,'edgecolor',[0.6 0.6 0.6]);
m_grid('fontsize',8,'linestyle','none','meridianlabel','off','parallellabel','off','xtick',[],'ytick',[]);
hold on
m_plot(lon,lat,'o','MarkerFaceColor',[0.7 0.7 0.7],'MarkerEdgeColor',[0.5 0.5 0.5]);

%Export and save? 
if save==false;
elseif save==true;
exportgraphics(gcf,'Map of stations.pdf','ContentType','vector','BackgroundColor','none');
end

pause(2)
close gcf

%% Abundance and Biovolume analysis

% create table with abundance (ind/m3) and biovolume (mm3/mm3/m3) per plankton group:
    %number of line corresponds to the number of sample (sampleid)
    %number of columns corresponds to the number of plankton group (Zoo_groups)

Ab=[]; 
Bv=[];
plankton_groups=[];

for i=1:m  
    Ab=[Ab;base(i).regroupped.Ab]; 
    Bv=[Bv;base(i).regroupped.Bv];
    plankton_groups=base(i).regroupped.Zoo_groups;
    taxo_annotation_shortname=[base(i).regroupped.Idlistshort;base(i).regroupped.Zoo_groups((size(base(1).regroupped.Idlistshort))+1:end)];
end

Ab_tot=Ab(:,base(1).regroupped.livingplace(1));
Bv_tot=Bv(:,base(1).regroupped.livingplace(1));

%% Taxon, function or trophic groups to exclude for analysis (on abundance and biovolume)

%%taxonomic analysis 
        
data_ab_taxo=Ab(:,base(1).regroupped.originalplace(1):base(1).regroupped.originalplace(2));
data_bv_taxo=Bv(:,base(1).regroupped.originalplace(1):base(1).regroupped.originalplace(2));
plankton_groups_data_taxo=taxo_annotation_shortname(base(1).regroupped.originalplace(1):base(1).regroupped.originalplace(2));

    list2 = {'Exclude manually', 'Load the list', 'No'};  
    [indx,tf] = listdlg('PromptString','Do you want to exclude somes of taxonomic annotation for the analysis?','ListString',list2,'SelectionMode','single','ListSize',[400,200]);
        
            if indx==1;
        
            %case 'Yes, exclude taxa manually'
            
            %select manually taxon to exclude
            list=plankton_groups_data_taxo;
            [indx2,tf] = listdlg('PromptString',{'Select taxon to exclude, several can be selected',...
            '(press control -PC- or command -Mac-)',''},'ListString',list,'ListSize',[400,200]);   
            
            %create and save 'taxon_to_exclude_ab_taxo' 
            taxon_to_exclude_taxo=indx2';
            taxon_to_exclude_taxo=num2cell(taxon_to_exclude_taxo);
            taxon_to_exclude_taxo(:,2)=plankton_groups_data_taxo(indx2);

            tosave=cell2table(taxon_to_exclude_taxo);
            [file,path] = uiputfile('taxon_to_exclude_taxo.xlsx');
            filename = fullfile(path,file);
            writetable(tosave,filename,'WriteVariableNames',0);
            
            %delete selected taxon
            data_ab_taxo(:,indx2)=[];
            data_bv_taxo(:,indx2)=[];
            plankton_groups_data_taxo(indx2)=[];

            elseif indx==2;

            f = msgbox('Select the taxon_to_exclude_taxo.xlsx list');
            uiwait(f)

            %case 'Yes, please load the list'
            
            %load 'taxon_to_exclude_ab_taxo'
            [file,path] = uigetfile('*.xlsx');
            taxon_to_exclude_taxo=readtable([path file],'ReadVariableNames',false);  %table
            indx2=table2array(taxon_to_exclude_taxo(:,1));

            %delete selected taxon
            data_ab_taxo(:,indx2)=[];
            data_bv_taxo(:,indx2)=[];
            plankton_groups_data_taxo(indx2)=[];

            elseif indx==3;
            
            %case 'No'
        
            end 

            clear indx
           
%%groups analysis 

data_ab_group=Ab(:,base(1).regroupped.planktongroupplace(1):base(1).regroupped.planktongroupplace(2));
data_bv_group=Bv(:,base(1).regroupped.planktongroupplace(1):base(1).regroupped.planktongroupplace(2));
plankton_groups_data_group=taxo_annotation_shortname(base(1).regroupped.planktongroupplace(1):base(1).regroupped.planktongroupplace(2));    

            % to have detritus by default
            default = 'detritus';
            index_default = find(strcmp(plankton_groups_data_group, default));
            plankton_groups_data_group = [plankton_groups_data_group(index_default); setdiff(plankton_groups_data_group, default, 'stable')];
            
            col_order = [index_default, setdiff(1:length(plankton_groups_data_group), index_default)];
            data_ab_group = data_ab_group(:, col_order);
            data_bv_group = data_bv_group(:, col_order);

    list2 = {'Exclude manually', 'Load the list', 'No'};  
    [indx,tf] = listdlg('PromptString','Do you want to exclude somes of taxonomic annotation for the analysis?','ListString',list2,'SelectionMode','single','ListSize',[400,200]);
        
            if indx==1;
            
            %case 'Yes, exclude function groups manually'
                          
            %select manually taxon to exclude
            list=plankton_groups_data_group;
            [indx2,tf] = listdlg('PromptString',{'Select taxon to exclude, several can be selected',...
            '(press control -PC- or command -Mac-)',''},'ListString',list,'ListSize',[400,200]);   
                
            %create and save 'taxon_to_exclude_group' 
            taxon_to_exclude_group=indx2';
            taxon_to_exclude_group=num2cell(taxon_to_exclude_group);
            taxon_to_exclude_group(:,2)=plankton_groups_data_group(indx2);

            tosave=cell2table(taxon_to_exclude_group);
            [file,path] = uiputfile('taxon_to_exclude_group.xlsx');
            filename = fullfile(path,file);
            writetable(tosave,filename,'WriteVariableNames',0);
            
            %delete selected taxon
            data_ab_group(:,indx2)=[];
            data_bv_group(:,indx2)=[];
            plankton_groups_data_group(indx2)=[];

            elseif indx==2;

            f = msgbox('Select the taxon_to_exclude_group.xlsx list');
            uiwait(f)

            %case 'Yes, please load the list'
            
            %load 'taxon_to_exclude_ab_taxo'
            [file,path] = uigetfile('*.xlsx');
            taxon_to_exclude_group=readtable([path file],'ReadVariableNames',false);  %table
            indx2=table2array(taxon_to_exclude_group(:,1));
    
            %delete selected taxon
            data_ab_group(:,indx2)=[];
            data_bv_group(:,indx2)=[];
            plankton_groups_data_group(indx2)=[];
    
            elseif indx==3; 
            
            %case 'No'
            
            end

%%trophic analysis

data_ab_troph=Ab(:,base(1).regroupped.trophicplace(1):base(1).regroupped.trophicplace(2));
data_bv_troph=Bv(:,base(1).regroupped.trophicplace(1):base(1).regroupped.trophicplace(2));
plankton_groups_data_troph=taxo_annotation_shortname(base(1).regroupped.trophicplace(1):base(1).regroupped.trophicplace(2));

    list2 = {'Exclude manually', 'Load the list', 'No'};  
    [indx,tf] = listdlg('PromptString','Do you want to exclude somes of taxonomic annotation for the analysis?','ListString',list2,'SelectionMode','single','ListSize',[400,200]);
        
            if indx==1;

            %plankton_groups_data_troph_name=({'-1 do not feed';'1 phototrophs';'1.5 mixotrophs';'2 grazers';'2.5 omnivorous';'3 predators';'3.5 unknown trophic group'})
            
            %case 'Yes, exclude trophic groups manually'
                
            %select manually taxon to exclude
            list=plankton_groups_data_troph;
            [indx2,tf] = listdlg('PromptString',{'Select taxon to exclude, several can be selected',...
            '(press control -PC- or command -Mac-)',''},'ListString',list,'ListSize',[400,200]);   
                
            %create and save 'taxon_to_exclude_ab_troph' 
            taxon_to_exclude_troph=indx2';
            taxon_to_exclude_troph=num2cell(taxon_to_exclude_troph);
            taxon_to_exclude_troph(:,2)=plankton_groups_data_troph(indx2);

            tosave=cell2table(taxon_to_exclude_troph);
            [file,path] = uiputfile('taxon_to_exclude_troph.xlsx');
            filename = fullfile(path,file);
            writetable(tosave,filename,'WriteVariableNames',0);
            
            %delete selected taxon
            data_ab_troph(:,indx2)=[];
            data_bv_troph(:,indx2)=[];
            plankton_groups_data_troph(indx2)=[];
            
            elseif indx==2;

            %case 'Yes, please load the list'

            f = msgbox('Select the taxon_to_exclude_troph.xlsx list');
            uiwait(f)
            
            %load 'taxon_to_exclude_ab_troph'
            [file,path] = uigetfile('*.xlsx');
            taxon_to_exclude_troph=readtable([path file],'ReadVariableNames',false);  %table
            indx2=table2array(taxon_to_exclude_troph(:,1));
    
            %delete selected taxon
            data_ab_troph(:,indx2)=[];
            data_bv_troph(:,indx2)=[];
            plankton_groups_data_troph(indx2)=[];

            elseif indx==3;
    
            %case 'No'
            
            end

            %create a base with the name of the trophic group
            for i = 1:numel(plankton_groups_data_troph)
            current_value = strtrim(plankton_groups_data_troph{i});
            if isequal(current_value, '-1')
            plankton_groups_data_troph_name{i} = 'not feeding';
            elseif isequal(current_value, '1')
            plankton_groups_data_troph_name{i} = 'phototrophs';
            elseif isequal(current_value, '1.5')
            plankton_groups_data_troph_name{i} = 'mixotrophs';
            elseif isequal(current_value, '2')
            plankton_groups_data_troph_name{i} = 'grazers';
            elseif isequal(current_value, '2.5')
            plankton_groups_data_troph_name{i} = 'omnivorous';
            elseif isequal(current_value, '3')
            plankton_groups_data_troph_name{i} = 'predators';
            elseif isequal(current_value, '3.5')
            plankton_groups_data_troph_name{i} = 'unknown';
            end
            end
            clear current_value;

%% Plankton community description for each station

for i=1:size(base,2); 

    figure('numbertitle', 'off', 'Position', [10 50 scrsz(3) scrsz(4) - 150], 'Color', 'white');

    base(i).sampleid = cellfun(@(x) strrep(x, '_', ' '), base(i).sampleid, 'UniformOutput', false);
    base(i).sampleid = cellfun(@(x) strrep(x, '_', ' '), base(i).sampleid, 'UniformOutput', false);
    plankton_groups_data_group = cellfun(@(x) strrep(x, '_', ' '), plankton_groups_data_group, 'UniformOutput', false);

    %colors for each plankton groups 
    p_group_colors;

    %compute Shannon index on all taxonomic groups 
    total=sum(data_ab_taxo(i,:),2);
    relativeabund=100*data_ab_taxo(i,:)./total;
    H=-nansum((relativeabund/100).*log(relativeabund/100),2);

    %compute absolute abundance on plankton group
    total=sum(data_ab_group(i,:),2);
    relativeabund=100*data_ab_group(i,:)./total;
   
    subplot(3,4,1)

    %%%%%Name and tot characteristics of the station

    t=text(0,0.8,char(base(i).sampleid),'FontSize',11,'FontWeight','bold'); axis off
    hold on
    text(0,0.5,['Total absolute abundance (ind. m^{-3})' num2cell(Ab_tot(i))],'FontSize',10);
    hold on
    text(0,0.3,['Total absolute biovolume (mm^{3}. mm^{-3})' num2cell(Bv_tot(i))],'FontSize',10);
    hold on
    text(0,0.1,['Shannon Index' num2cell(H)],'FontSize',10); 

    subplot(3,4,2)
    ax2 = subplot(3,4,2);
    pos2 = get(ax2, 'Position');
    set(ax2, 'Position', [pos2(1) + 0.03, pos2(2), pos2(3), pos2(4)]);

    %%%%%NBSS size spectra on living 
    
    ESD=base(i).tot.ESDvector;
    X=base(i).tot.X;

    plot(ESD, base(i).regroupped.Ybv_Ellipsoid_BV_spectra(:,base(i).regroupped.livingplace), 'o', 'MarkerFaceColor', [0.9, 0.9, 0.9], 'MarkerEdgeColor', [0.5, 0.5, 0.5]);

    set(gca, 'YScale', 'log', 'XScale', 'log');

    xticks([0.1 0.2 0.5 1 2 5 10 20 50 100 200 500 1000 2000 5000 10000 20000 50000 100000]);
    xticklabels({'0.1', '0.2', '0.5', '1', '2', '5', '10', '20', '50', '100', '200', '500', '1000', '2000', '5000', '10000', '20000', '50000', '100000'});

    xlabel('size on ESD (µm)','fontsize',8)
    ylabel({'NBSS (mm^3/mm^3/m^3)';'equivalent abundance'},'fontsize',8)
    title('NBSS','fontsize',10)

    subplot(3,4,3)

    %%%%%Pie relative abundance 

    colormap(plankton_group_colors);    

    labels = repmat({''},size(relativeabund));
    pie(relativeabund,labels)
    title('Relative abundance','fontsize',10)    

    %%%%%Legend plankton groups for abundance

    l=legend(plankton_groups_data_group,'fontsize',8,'Location','west','NumColumns',2);
    rect = [.72, .45, .2, .41];
    set(l, 'Position', rect)
    title(l,'Plankton groups','fontsize',10)

    f=subplot(3,4,5);

    %%%%%Trophic analysis (pyramid) 

    %create a base of color depending on the trophic group 
    color_troph = zeros(numel(plankton_groups_data_troph), 3);
    for k = 1:numel(plankton_groups_data_troph)
    current_value = strtrim(plankton_groups_data_troph{k});
    if isequal(current_value, '-1')
    color_troph(k,:) = [1.00 1.00 1.00];
    elseif isequal(current_value, '1')
    color_troph(k,:) = [0.40 0.69 0.39];
    elseif isequal(current_value, '1.5')
    color_troph(k,:) = [0.63 0.77 0.53];
    elseif isequal(current_value, '2')
    color_troph(k,:) = [0.3	0.54 0.69];
    elseif isequal(current_value, '2.5')
    color_troph(k,:) = [1.00 0.93 0.66];
    elseif isequal(current_value, '3')
    color_troph(k,:) = [0.78 0.33 0.1490];
    elseif isequal(current_value, '3.5')
    color_troph(k,:) = [0.90 0.90 0.90];
    clear current_value 
    end
    end

    height=str2double(plankton_groups_data_troph);
    width=data_bv_troph(i,:);
    width=log(width+1);
    h = [];

    [n,m]=size(plankton_groups_data_troph);

    for j=1:n;
    x=[-width(j) width(j) width(j) -width(j)];
    y=[height(j)-0.25 height(j)-0.25 height(j)+0.25 height(j)+0.25];
    hold on
    h = [h,patch(x,y,color_troph(j,:))];
    end

    xlabel({'log biovolume +1 (mm^3/m^3)';'equivalent abundance'},'fontsize',8)
    title('Trophic pyramids','fontsize',10)

    %Legend trophic

    l=legend(plankton_groups_data_troph_name,'fontsize',8,'Location','west','NumColumns',2);
    rect = [.130,.150,.15,.15];
    set(l, 'Position', rect)
    title(l,'Trophic groups','fontsize',10)

    subplot(3,4,6)
    ax6 = subplot(3,4,6);
    pos6 = get(ax6, 'Position');
    set(ax6, 'Position', [pos6(1) + 0.03, pos6(2), pos6(3), pos6(4)]);

    %%%%%BSS size spectra on living 

    %outpout biovolume per size denormalize
    livingplace=base(i).regroupped.livingplace;
    Ybv_denorm=(base(i).regroupped.Ybv_Ellipsoid_BV_spectra(:,base(i).regroupped.livingplace).*base(i).tot.X1);

    plot(ESD, Ybv_denorm, 'o', 'MarkerFaceColor', [0.9, 0.9, 0.9], 'MarkerEdgeColor', [0.5, 0.5, 0.5]);

    set(gca, 'YScale', 'log', 'XScale', 'log');

    xticks([0.1 0.2 0.5 1 2 5 10 20 50 100 200 500 1000 2000 5000 10000 20000 50000 100000]);
    xticklabels({'0.1', '0.2', '0.5', '1', '2', '5', '10', '20', '50', '100', '200', '500', '1000', '2000', '5000', '10000', '20000', '50000', '100000'});

    xlim_bss = xlim(gca);
     
    xlabel('size on ESD (µm)','fontsize',8)
    ylabel({'BSS (mm^3/m^3)'},'fontsize',8)
    title('BSS','fontsize',10)

    subplot(3,4,7)

    %%%%%Pie relative biovolume

    colormap(plankton_group_colors);

    total=sum(data_bv_group(i,:),2);
    relativebv=100*data_bv_group(i,:)./total;
    labels = repmat({''},size(relativebv));
    pie(relativebv,labels)
    title('Relative biovolume','fontsize',10)

    subplot(3,4,10)
    ax10 = subplot(3,4,10);
    pos10 = get(ax10, 'Position');
    set(ax10, 'Position', [pos10(1) + 0.03, pos10(2), pos10(3), pos10(4)]);

    %%%%%Biovolume composition (in %) per size class (ESD) 

    %outpout biovolume per size denormalize
    planktongroupplace=base(i).regroupped.planktongroupplace;
    Ybv_denorm_group=(base(i).regroupped.Ybv_Ellipsoid_BV_spectra(:,[planktongroupplace(1):planktongroupplace(2)]).*base(i).tot.X1);
    Ybv_denorm_group = Ybv_denorm_group(:, col_order);

    %taxon (plankton group) to exclude
    if iscell(taxon_to_exclude_group);
    indx=cell2mat(taxon_to_exclude_group(:,1));
    elseif istable(taxon_to_exclude_group);
    indx=table2array(taxon_to_exclude_group(:,1));
    end 

    Ybv_denorm_group(:,indx)=[];

    %biovolume denormalize per size in percent: relative (%) 
    Ybv_denorm_tot_persize=nansum(Ybv_denorm_group,2);
    relative_Ybv=100*(Ybv_denorm_group./Ybv_denorm_tot_persize);

    %remove NaN on the variable for the plot 
    columns_with_nan = not(Ybv_denorm_tot_persize==0);
    relative_Ybv_without_nan=relative_Ybv(columns_with_nan,:);
    ESD_without_nan=ESD(columns_with_nan,:);

    h=bar(log10(ESD_without_nan),relative_Ybv_without_nan,'stacked','BarWidth', 1);
    for j = 1:size(plankton_group_colors, 1);
    h(j).FaceColor = plankton_group_colors(j, :);
    end
    set(h, 'EdgeColor', 'none');
   
    ylim([0 100]);
    xlim(log10(xlim_bss));

    x_ticks = [0.1 0.2 0.5 1 2 5 10 20 50 100 200 500 1000 2000 5000 10000 20000 50000 100000];
    set(gca, 'XTick', log10(x_ticks), 'XTickLabel', arrayfun(@num2str, x_ticks, 'UniformOutput', false));

    xlabel('size on ESD (µm)','fontsize',8)
    ylabel({'Biovolume (mm^3/m^3)'},'fontsize',8)
    title('Biovolume composition per size class','fontsize',10)

    %Map

    subplot(3,4,[11,12])
    m_proj('mercator','lon',[lonmin lonmax],'lat',[latmin latmax]);
    m_coast('patch',landcolor,'edgecolor',[0.6 0.6 0.6]);
    m_grid('fontsize',8,'linestyle','none','meridianlabel','off','parallellabel','off','xtick',[],'ytick',[]);
    hold on
    m_scatter(lon(i),lat(i),'o','MarkerFaceColor',[0.7 0.7 0.7],'MarkerEdgeColor',[0.5 0.5 0.5]);

%Export and save? 
if save==false;
elseif save==true;
outputDir = 'raw analysis per station';
mkdir(outputDir);
filename = fullfile(outputDir, sprintf('%s.pdf', char(base(i).sampleid)));
exportgraphics(gcf, filename, 'ContentType', 'vector', 'BackgroundColor', 'none');
end

pause(1)
close gcf

end

%% Plankton community description for all stations

list3 = {'Abundance', 'Biovolume'};  
[indx3,tf] = listdlg('PromptString','Data to plot for H-Bar and Map per station','ListString',list3,'SelectionMode','single','ListSize',[400,200]);

% H-bar per station for abundance and biovolume all taxonomic

    if indx3==1
    data=data_ab_taxo;
    elseif indx3==2 
    data=data_bv_taxo;
    end

[n,m]=size(data);

figure('numbertitle','off','Position',[10 50 scrsz(3) scrsz(4)-150],'Color','white');

%colors for each plankton groups 
p_group_colors;
num_groups = length(plankton_groups_data_taxo);
num_base_colors = size(plankton_group_colors_palet, 1);

cmap = zeros(num_groups, 3);
for i = 1:num_groups
    base_idx = floor((i-1) * (num_base_colors-1) / (num_groups-1)) + 1;
    base_idx_next = min(base_idx + 1, num_base_colors);
    alpha = (i - 1) * (num_base_colors-1) / (num_groups-1) - (base_idx - 1);
    cmap(i, :) = (1 - alpha) .* plankton_group_colors_palet(base_idx, :) + alpha .* plankton_group_colors_palet(base_idx_next, :);
end

h=bar(data,'stacked','FaceColor','flat');
for k = 1:size(data,2);
    h(k).CData = k;
end

colormap(cmap);
set(gca,'TickLength',[0.002 0.002]);
set(h,'Edgecolor',[0 0 0]);
set(gca,'xtick',1:n);
set(gca,'xticklabel',sampleid,'fontsize',8);
set(gca,'XTickLabelRotation',75);
set(gca,'xlim',[0 n+1],'ylim',[0 max(sum(data,2))]);

if indx3==1;
for_legend = "Abundance (ind.m^{-3}) "; 
elseif indx3==2;
for_legend = "Biovolume (mm.m^{-3}) ";
end
ylabel(for_legend,'fontsize', 10);

l=legend(h,plankton_groups_data_taxo,'fontsize',8,'Location','eastoutside','NumColumns',2);
title(l,'Plankton Groups','fontsize',10);

%Export and save? 
if save==false;
elseif save==true;
            outputDir = 'global raw analysis';
                mkdir(outputDir);
    if indx3==1
    filename = fullfile(outputDir, 'Hbar of abundance.pdf');
    exportgraphics(gcf, filename, 'ContentType', 'vector', 'BackgroundColor', 'none');
    elseif indx3==2
    filename = fullfile(outputDir, 'Hbar of biovolume.pdf');
    exportgraphics(gcf, filename, 'ContentType', 'vector', 'BackgroundColor', 'none');
    end 
end

pause(2)
close gcf

% H-bar per station for abundance and biovolume per group

    if indx3==1
    data=data_ab_group;
    elseif indx3==2 
    data=data_bv_group;
    end

[n,m]=size(data);

figure('numbertitle','off','Position',[10 50 scrsz(3) scrsz(4)-150],'Color','white');

%colors for each plankton groups 
p_group_colors;

h=bar(data,'stacked','FaceColor','flat');
for k = 1:size(data,2);
    h(k).CData = k;
end

colormap(plankton_group_colors);
set(gca,'TickLength',[0.002 0.002]);
set(h,'Edgecolor',[0 0 0]);
set(gca,'xtick',1:n);
set(gca,'xticklabel',sampleid,'fontsize',8);
set(gca,'XTickLabelRotation',75);
set(gca,'xlim',[0 n+1],'ylim',[0 max(sum(data,2))]);

if indx3==1;
for_legend = "Abundance (ind.m^{-3}) "; 
elseif indx3==2;
for_legend = "Biovolume (mm.m^{-3}) ";
end
ylabel(for_legend,'fontsize', 10);

l=legend(h,plankton_groups_data_group,'fontsize',8,'Location','eastoutside','NumColumns',1);
title(l,'Plankton Groups','fontsize',10);

%Export and save? 
if save==false;
elseif save==true;
            outputDir = 'global raw analysis';
                mkdir(outputDir);
    if indx3==1
    filename = fullfile(outputDir, 'Hbar of abundance.pdf');
    exportgraphics(gcf, filename, 'ContentType', 'vector', 'BackgroundColor', 'none');
    elseif indx3==2
    filename = fullfile(outputDir, 'Hbar of biovolume.pdf');
    exportgraphics(gcf, filename, 'ContentType', 'vector', 'BackgroundColor', 'none');
    end 
end

pause(2)
close gcf

% Map absolute abundance or biovolume

data=sum(data,2);

figure('numbertitle','off','Position',[10 50 scrsz(3) scrsz(4)-150],'Color','white');
m_proj('mercator','lon',[lonmin lonmax],'lat',[latmin latmax]);
m_coast('patch',landcolor,'edgecolor',[0.6 0.6 0.6]);
hold on 
m_grid('fontsize',8,'linestyle','none','meridianlabel','off','parallellabel','off','xtick',[],'ytick',[]);

log_data=log10(data);

h = m_scatter(lon, lat, 50, log_data, 'filled', 'markeredgecolor', [0.5 0.5 0.5]);

load('slanCM_Data.mat');
color=slandarerCM(6).Colors(1,9); %spectral (change manually)
color=color{1, 1};
c2 = colorbar;
colormap(color);
colormap(flipud(color));

current_ticks = c2.Ticks;
new_tick_labels = 10.^current_ticks;
is_power_of_ten = @(x) abs(log10(x) - round(log10(x))) < 1e-10;
rounded_tick_indices = is_power_of_ten(new_tick_labels); 

set(c2, 'Ticks', current_ticks(rounded_tick_indices), ... 
         'TickLabels', arrayfun(@num2str, new_tick_labels(rounded_tick_indices), 'UniformOutput', false));
num2str('%.0e')

if indx3==1; 
for_legend = "Abundance (ind.m^{-3}) "; 
elseif indx3==2
for_legend = "Biovolume (mm.m^{-3}) ";
end 

title(c2,for_legend,'fontsize', 11,'fontweight', 'bold')

%Export and save? 
if save==false;
elseif save==true;
            outputDir = 'global raw analysis';
                mkdir(outputDir);
    if indx3==1
    filename = fullfile(outputDir, 'Map of abundance.pdf');
    exportgraphics(gcf, filename, 'ContentType', 'vector', 'BackgroundColor', 'none');
    elseif indx3==2
    filename = fullfile(outputDir, 'Map of biovolume.pdf');
    exportgraphics(gcf, filename, 'ContentType', 'vector', 'BackgroundColor', 'none');
    end 
end

pause(2)
close gcf

% Map Shannon diversity 

[n,m]=size(base);

for i=1:numel(base);

    %compute Shannon index on all taxonomic groups 
    total=sum(data_ab_taxo(i,:),2);
    relativeabund=100*data_ab_taxo(i,:)./total;
    H=-nansum((relativeabund/100).*log(relativeabund/100),2);
    diversityH(i,:)=H;
    clear H

end

figure('numbertitle','off','Position',[10 50 scrsz(3) scrsz(4)-150],'Color','white');
m_proj('mercator','lon',[lonmin lonmax],'lat',[latmin latmax]);
m_coast('patch',landcolor,'edgecolor',[0.6 0.6 0.6]);
hold on 
m_grid('fontsize',8,'linestyle','none','meridianlabel','off','parallellabel','off','xtick',[],'ytick',[]);

color=slandarerCM(2).Colors(1,15); %YlGnBu 
color=color{1, 1};
colormap(color);
colormap(flipud(color));

h = m_scatter(lon, lat, 50, diversityH, 'filled', 'markeredgecolor', [0.5 0.5 0.5]);

c = colorbar;

title(c,'Shannon Index','fontsize', 11,'fontweight', 'bold')

%Export and save? 
if save==false;
elseif save==true;
            outputDir = 'global raw analysis';
            mkdir(outputDir);
filename = fullfile(outputDir, 'Map of diversity.pdf');
exportgraphics(gcf,filename,'ContentType','vector','BackgroundColor','none');
end

pause(2)
close gcf

%% Classic analysis per group 

list5 = {'Yes', 'No'};  
[indx3, tf] = listdlg(...
    'PromptString', 'Do you want to select just one taxon or group for analysis?', ...
    'ListString', list5, ...
    'SelectionMode', 'single', ...
    'ListSize', [400, 200]); 

if indx3==1;

% Choose one group

list3 = {'Abundance', 'Biovolume'};  
[indx3,tf] = listdlg('PromptString','Data to plot for H-Bar and Map for one group','ListString',list3,'SelectionMode','single','ListSize',[400,200]);

list1 = [{'one taxon'}; {'one plankton group'}; {'one trophic group'}; {'all taxon individually'}; {'all plankton group'}; {'all trophic group'}];
[indx1,tf] = listdlg('PromptString','Choose the plankton grouping level for figure','ListString',list1,'SelectionMode','single','ListSize',[400,200]);

if indx1==1;

    list2 = plankton_groups_data_taxo;
    [indx2,tf] = listdlg('PromptString','Choose the plankton grouping level for figure','ListString',list2,'SelectionMode','single','ListSize',[400,200]);
    
    plankton_group=list2{indx2};

    if indx3==1
    data=data_ab_taxo(:,indx2);
    elseif indx3==2 
    data=data_bv_taxo(:,indx2);
    end   

elseif indx1==2

    list2 = plankton_groups_data_group;
    [indx2,tf] = listdlg('PromptString','Choose the plankton grouping level for figure','ListString',list2,'SelectionMode','single','ListSize',[400,200]);
    
    plankton_group=list2{indx2};

    if indx3==1
    data=data_ab_group(:,indx2);
    elseif indx3==2 
    data=data_bv_group(:,indx2);
    end  

elseif indx1==3

    list2 = plankton_groups_data_troph_name;
    [indx2,tf] = listdlg('PromptString','Choose the plankton grouping level for figure','ListString',list2,'SelectionMode','single','ListSize',[400,200]);
    
    plankton_group=list2{indx2};

    if indx3==1
    data=data_ab_troph(:,indx2);
    elseif indx3==2 
    data=data_bv_troph(:,indx2);
    end  

elseif indx1==4 
    
    plankton_group=plankton_groups_data_taxo;
    if indx3==1
    data=data_ab_taxo;
    elseif indx3==2 
    data=data_bv_taxo;
    end   

elseif indx1==5
    
    plankton_group=plankton_groups_data_group;
    if indx3==1
    data=data_ab_group;
    elseif indx3==2 
    data=data_bv_group;
    end  

elseif indx1==6
    
    plankton_group=plankton_groups_data_troph_name;
    if indx3==1
    data=data_ab_troph;
    elseif indx3==2 
    data=data_bv_troph;
    end  

end

plankton_group = cellstr(plankton_group);

% H-bar per station for abundance and biovolume for one group

[n,m]=size(data);

%colors for each plankton groups 
p_group_colors;

for i=1:size(data,2);

figure('numbertitle','off','Position',[10 50 scrsz(3) scrsz(4)-150],'Color','white');

h=bar(data(:,i),'stacked','FaceColor','flat');
for k = 1:size(data(:,i),2);
    h(k).CData = k;
end

set(gca,'TickLength',[0.002 0.002]);
set(h,'Facecolor',[0.8 0.8 0.8]);
set(h,'Edgecolor',[0 0 0]);
set(gca,'xtick',1:n);
set(gca,'xticklabel',sampleid,'fontsize',8);
set(gca,'XTickLabelRotation',75);
set(gca,'xlim',[0 n+1],'ylim',[0 max(sum(data(:,i),2))]);

if indx3==1;
for_legend = "Abundance (ind.m^{-3}) "; 
elseif indx3==2;
for_legend = "Biovolume (mm.m^{-3}) ";
end
ylabel(for_legend,'fontsize', 10);

l=legend(h,plankton_group(i),'fontsize',8,'Location','eastoutside','NumColumns',1);
title(l,'Plankton Groups','fontsize',10);

%Export and save? 
if save==false;
elseif save==true;
    
    if indx3==1 %abundance
        
        parentDir = 'global raw analysis';
        outputDir = fullfile(parentDir, ['Hbar abundance for ' list1{indx1}]);        
            mkdir(outputDir);

    filename = fullfile(outputDir, sprintf('Hbar of abundance %s.pdf', char(plankton_group{i})));
    exportgraphics(gcf, filename, 'ContentType', 'vector', 'BackgroundColor', 'none');

    elseif indx3==2 %biovolume 

        parentDir = 'global raw analysis';
        outputDir = fullfile(parentDir, ['Hbar biovolume for ' list1{indx1}]);
            mkdir(outputDir);

    filename = fullfile(outputDir, sprintf('Hbar of biovolume %s.pdf', char(plankton_group{i})));
    exportgraphics(gcf, filename, 'ContentType', 'vector', 'BackgroundColor', 'none');

    end 

end 

pause(1)
close gcf

end 

% Map absolute abundance or biovolume for one group

for i=1:size(data,2);

datamap=sum(data(:,i),2);

figure('numbertitle','off','Position',[10 50 scrsz(3) scrsz(4)-150],'Color','white');

m_proj('mercator','lon',[lonmin lonmax],'lat',[latmin latmax]);
m_coast('patch',landcolor,'edgecolor',[0.6 0.6 0.6]);
hold on 
m_grid('fontsize',8,'linestyle','none','meridianlabel','off','parallellabel','off','xtick',[],'ytick',[]);

load('slanCM_Data.mat');
color=slandarerCM(6).Colors(1,9); %spectral (change manually)
color=color{1, 1};
colormap(color);
colormap(flipud(color));

min_size = 10; 
max_size = 100;

scaled_size = min_size + (datamap - min(datamap)) / (max(datamap) - min(datamap)) * (max_size - min_size);

h = m_scatter(lon, lat, scaled_size, datamap, 'filled', 'markeredgecolor', [0.5 0.5 0.5]);

load('slanCM_Data.mat');
color=slandarerCM(6).Colors(1,9); %spectral (change manually)
color=color{1, 1};
c2 = colorbar;
colormap(color);
colormap(flipud(color));

current_ticks = c2.Ticks;
new_tick_labels = 10.^current_ticks;
is_power_of_ten = @(x) abs(log10(x) - round(log10(x))) < 1e-10;
rounded_tick_indices = is_power_of_ten(new_tick_labels); 

set(c2, 'Ticks', current_ticks(rounded_tick_indices), ... 
         'TickLabels', arrayfun(@num2str, new_tick_labels(rounded_tick_indices), 'UniformOutput', false));
num2str('%.0e')

%Export and save? 
if save==false;
elseif save==true;
    
    if indx3==1 %abundance
        
        parentDir = 'global raw analysis';
        outputDir = fullfile(parentDir, ['Map abundance for ' list1{indx1}]);
        mkdir(outputDir);

    filename = fullfile(outputDir, sprintf('Map of abundance %s.pdf', char(plankton_group{i})));
    exportgraphics(gcf, filename, 'ContentType', 'vector', 'BackgroundColor', 'none');

    elseif indx3==2 %biovolume 

        parentDir = 'global raw analysis';
        outputDir = fullfile(parentDir, ['Map biovolume for ' list1{indx1}]);
        mkdir(outputDir);

    filename = fullfile(outputDir, sprintf('Map of biovolume %s.pdf', char(plankton_group{i})));
    exportgraphics(gcf, filename, 'ContentType', 'vector', 'BackgroundColor', 'none');

    end 

end 

pause(1)
close gcf

end

end

%% Multivariates analysis 

%Choose data to analyse
list3_a = {'Abundance per taxo','Abundance per plankton group','Abundance per trophic group','Biovolume per taxo','Biovolume per plankton group','Biovolume per trophic group'}; 
[indx3,tf] = listdlg('PromptString','What data should be used for analysis?','ListString',list3_a,'SelectionMode','single','ListSize',[400,200]);
list3 = {'data_ab_taxo','data_ab_group','data_ab_troph','data_bv_taxo','data_bv_group','data_bv_troph'}; 

    dataACP_name = list3{indx3};
    data = eval(dataACP_name);

    if indx3 == 1 || indx3 == 4
    group=plankton_groups_data_taxo;
    elseif indx3 == 2 || indx3 == 5
    group=plankton_groups_data_group;
    elseif indx3 == 3 || indx3 == 6
    group=plankton_groups_data_troph_name;
    end

    for i = 1:numel(group)
    group{i} = strrep(group{i}, '_', ' ');
    end

%Choose normalisation
list3 = {'No normalisation',  'Log x+1, more weight to rare', 'Double cube root, more weight to rare', 'Relative abundance', 'Hellinger transform, relative + more weight to rare', 'Centred and reduced, forced normalisation loss of relative strength', 'Box cox, crude forced normalisation loss of relative strength'}; 
[normalisation,tf] = listdlg('PromptString','With what normalisation?','ListString',list3,'SelectionMode','single','ListSize',[600,200]);

if normalisation==1; %No
data_stat=data';

elseif normalisation==2; %logx+1
data_stat=log(data+1)';

elseif normalisation==3;
data_stat=((data').^0.25);

elseif normalisation==4; %relative abundance
total=sum(data,2);
[n,m]=size(data);
total=total*ones(1,m);
relativeabund=100*data./total;
data_stat=relativeabund';

elseif normalisation==5; %hellinger transform
total=sum(data,2);
[n,m]=size(data);
total=total*ones(1,m);
relativeabund=100*data./total;
data_stat=(relativeabund').^0.5;

elseif normalisation==6;
data_stat=normalize(data);
data_stat=data_stat';

elseif normalisation==7;
for i=1:size(data,2);
[transdat,lambda] = boxcox(data(:,i)+0.0000001);
transdattemp(:,i)=transdat;
end
data_stat=normalize(transdattemp);
data_stat=data_stat';

end

%to save
parentDir = fullfile(pwd, 'multivariates analysis');
mkdir(parentDir);

%% PCA

%to save
subDir = 'PCA';
fullSubDirPath = fullfile(parentDir, subDir);
    mkdir(fullSubDirPath);

[pc,zscores,latent,tsquared,pcvars] = pca(data_stat','Algorithm','svd'); 

% Scale scores to fit within -1 to 1 bounds (after eq. 1.10 in L&L, 1998):
zscores = (zscores./repmat(max(max(abs(zscores))), size(zscores)));

% PCA plot: 

figure('numbertitle','off','Position',[10 50 scrsz(3) scrsz(4)-150],'Color','white');

scatter(zscores(:,1),zscores(:,2),50,'filled','MarkerFaceColor', [0.9 0.9 0.9],'markeredgecolor',[0.5 0.5 0.5])

% Plot in function of the number of groups 

if numel(group) <= 15

[p,q]=size(pc);
hold on
for i=1:p
    if ((pc(i,1)^2+pc(i,2)^2)^0.5)
    plot([0 pc(i,1)],[0 pc(i,2)],'k-')
    text(pc(i,1),pc(i,2),group(i),'fontsize',8);
    end
end

elseif numel(group) > 15

[p,q]=size(pc);
hold on

distances = sqrt(pc(:,1).^2 + pc(:,2).^2);
[~, sortedIndices] = sort(distances, 'descend');
topIndices = sortedIndices(1:round(0.3 * numel(distances)));

for i=topIndices'
    if distances(i)
    plot([0 pc(i,1)],[0 pc(i,2)],'k-')
    text(pc(i,1),pc(i,2),group(i),'fontsize',8);
    end
end

end

xlabel(['PCA 1 ' num2str(fix(pcvars(1))) '% variance'],'fontsize',10);
ylabel(['PCA 2 ' num2str(fix(pcvars(2))) '% variance'],'fontsize',10);

%Export and save? 
if save==false;
elseif save==true;

        outputDir = fullfile(parentDir, subDir);
            mkdir(outputDir);

    filename = fullfile(outputDir, sprintf(['PCA on ' char(dataACP_name) ' ' list3{normalisation} '.pdf']));
    exportgraphics(gcf, filename, 'ContentType', 'vector', 'BackgroundColor', 'none');
end
pause(1)
close gcf

% Representation with RGB colors 
% Extract RGB color on the PCA

    RGB=zscores(:,1:3);
    minRGB=min(RGB);
    RGB=RGB-ones(size(data_stat,2),1)*minRGB;
    maxRGB=max(RGB);
    RGB=RGB./(ones(size(data_stat,2),1)*maxRGB);
    RGB=1-RGB;

    axis1=1;
    axis2=2;
    axis3=3;

% Plot PCA with RGB colors 

figure('numbertitle','off','Position',[10 50 scrsz(3) scrsz(4)-150],'Color','white');

subplot('Position', [0.05 0.1 0.3 0.8])

scatter(zscores(:,axis1),zscores(:,axis2),50,RGB,'filled','markeredgecolor',[0.5 0.5 0.5])

[p,q]=size(pc);
pc2=pc;
hold on

if numel(group) <= 15

[p,q]=size(pc);
hold on
for i=1:p
    if ((pc(i,1)^2+pc(i,2)^2)^0.5)
    plot([0 pc(i,1)],[0 pc(i,2)],'k-')
    text(pc(i,1),pc(i,2),group(i),'fontsize',8);
    end
end

elseif numel(group) > 15

[p,q]=size(pc);
hold on

distances = sqrt(pc(:,1).^2 + pc(:,2).^2);
[~, sortedIndices] = sort(distances, 'descend');
topIndices = sortedIndices(1:round(0.3 * numel(distances)));

for i=topIndices'
    if distances(i)
    plot([0 pc(i,1)],[0 pc(i,2)],'k-')
    text(pc(i,1),pc(i,2),group(i),'fontsize',8);
    end
end

end

xlabel(['PCA 1 ' num2str(fix(pcvars(1))) '% variance'],'fontsize',10);
ylabel(['PCA 2 ' num2str(fix(pcvars(2))) '% variance'],'fontsize',10);

% Map with RGB colors 

subplot('Position', [0.4 0.1 0.55 0.8])

m_proj('mercator','lon',[lonmin lonmax],'lat',[latmin latmax]);
m_coast('patch',landcolor,'edgecolor',[0.6 0.6 0.6]);
hold on 
m_grid('fontsize',8,'linestyle','none','meridianlabel','off','parallellabel','off','xtick',[],'ytick',[]);

[n,m] = size(data_stat);

m_scatter(lon,lat, 50, RGB, 'filled', 'markeredgecolor', [0.5 0.5 0.5]);

%Export and save? 
if save==false;
elseif save==true;

        outputDir = fullfile(parentDir, subDir);
            mkdir(outputDir);

    filename = fullfile(outputDir, sprintf(['PCA and RGB color on map ' char(dataACP_name) ' ' list3{normalisation} '.pdf']));
    exportgraphics(gcf, filename, 'ContentType', 'vector', 'BackgroundColor', 'none');
end
pause(1)
close gcf

%% PCoA 

%to save
subDir = 'PCoA';
fullSubDirPath = fullfile(parentDir, subDir);
    mkdir(fullSubDirPath);

%Methods for matrice (lots of other choices in f_dis Fathom)
list6 = {'Euclidean, recommanded for Hellinger normalisation','Bray curtis, recommanded for others normalisation'}; 
[indx6,tf] = listdlg('PromptString','Distance methods for PCoA?','ListString',list6,'SelectionMode','single','ListSize',[400,200]);

if indx6==1
dis = f_dis(data_stat','euc');  
elseif indx6==2
dis = f_dis(data_stat','bc');  
end

pcoa = f_pcoa_modified(dis,1); 

% Scale scores to fit within -1 to 1 bounds (after eq. 1.10 in L&L, 1998):
scores = pcoa.scores;
expl   = pcoa.expl;
scores = (scores./repmat(max([max(abs(scores(:,1)));max(abs(scores(:,2)))]), size(scores)));

%PCOA plot (adapted from f_pcoaPlot (Fathom): 

figure('numbertitle','off','Position',[10 50 scrsz(3) scrsz(4)-150],'Color','white');

scatter(scores(:,1), scores(:,2), 50, 'filled', 'MarkerFaceColor', [0.9 0.9 0.9], 'markeredgecolor', [0.5 0.5 0.5]);
%text(scores(:,1),scores(:,2),sampleid,'HorizontalAlignment','center', 'VerticalAlignment','middle','Color','k','Interpreter','tex'); 

Y=data_stat';
ncY = size(Y,2);
s = size(scores,2);
vec = zeros(ncY,s);
   for i = 1:s
      for j = 1:ncY
       vec(j,i) = f_corr(Y(:,j),scores(:,i))*std(Y(:,j))/std(scores(:,i));
      end
   end
vec = (vec./repmat(max([max(abs(vec(:,1)));max(abs(vec(:,2)))]), size(vec)));

eDis = f_dis([0 0;vec(:,1:2)],'euc');
eDis = eDis(2:end,1);

% Set delta inversely proportional to length of longest vector:
ratio = (abs(eDis./ repmat(max(eDis),size(eDis))));
ratio = ones(size(ratio))./ratio;
delta = 1 + (1* repmat(0,size(eDis)) .* ratio);
    
hold on 

if numel(group) <= 15

       n = size(vec,1); % get # of variables
   for j = 1:n
      % Plot vectors:
      f_arrow([0 0],[vec(j,1) vec(j,2)],'size',0.125*0.75,'angle',20,'Color','k')
      
      % Label vectors:
      h = text(vec(j,1)*delta(j),vec(j,2)*delta(j),group(j),'fontsize',8);
      set(h,'FontSize',8,'HorizontalAlignment','center','Color','k','Interpreter','tex');
   end

elseif numel(group) > 15

           n = size(vec,1); % get # of variables

distances = sqrt(vec(:,1).^2 + vec(:,2).^2);
[~, sortedIndices] = sort(distances, 'descend');
topIndices = sortedIndices(1:round(0.4 * numel(distances)));

for j=topIndices'
    if distances(j)
      % Plot vectors:
      f_arrow([0 0],[vec(j,1) vec(j,2)],'size',0.125*0.75,'angle',20,'Color','k')
      
      % Label vectors:
      h = text(vec(j,1)*delta(j),vec(j,2)*delta(j),group(j),'fontsize',8);
      set(h,'FontSize',8,'HorizontalAlignment','center','Color','k','Interpreter','tex');
    end
end

end

axis1 = sprintf('%2.2f',expl(1,1));
axis2 = sprintf('%2.2f',expl(2,1));
xlabel(['PCoA 1 (' num2str(axis1) ' % variance)'],'fontsize',10);
ylabel(['PCoA 2 (' num2str(axis2) ' % variance)'],'fontsize',10);

%Export and save? 
if save==false;
elseif save==true;

        outputDir = fullfile(parentDir, subDir);
            mkdir(outputDir);

    filename = fullfile(outputDir, sprintf(['PCoA on ' char(dataACP_name) ' ' list3{normalisation} '.pdf']));
    exportgraphics(gcf, filename, 'ContentType', 'vector', 'BackgroundColor', 'none');
end

pause(1)
close gcf

% PCoA plot representation with RGB

%Representation with RGB colors 
%Extract RGB color on the PCoA

    RGB=scores(:,1:3);
    minRGB=min(RGB);
    RGB=RGB-ones(size(data_stat,2),1)*minRGB;
    maxRGB=max(RGB);
    RGB=RGB./(ones(size(data_stat,2),1)*maxRGB);
    RGB=1-RGB;

    axis1=1;
    axis2=2;
    axis3=3;

%Plot PCoA with RGB colors

figure('numbertitle','off','Position',[10 50 scrsz(3) scrsz(4)-150],'Color','white');

subplot('Position', [0.05 0.1 0.3 0.8])

scatter(scores(:,1), scores(:,2), 50, RGB, 'filled', 'markeredgecolor', [0.5 0.5 0.5]);
%text(scores(:,1),scores(:,2),sampleid,'HorizontalAlignment','center', 'VerticalAlignment','middle','Color','k','Interpreter','tex'); 

Y=data_stat';
ncY = size(Y,2);
s = size(scores,2);
vec = zeros(ncY,s);
   for i = 1:s
      for j = 1:ncY
       vec(j,i) = f_corr(Y(:,j),scores(:,i))*std(Y(:,j))/std(scores(:,i));
      end
   end
vec = (vec./repmat(max([max(abs(vec(:,1)));max(abs(vec(:,2)))]), size(vec)));

eDis = f_dis([0 0;vec(:,1:2)],'euc');
eDis = eDis(2:end,1);

% Set delta inversely proportional to length of longest vector:
ratio = (abs(eDis./ repmat(max(eDis),size(eDis))));
ratio = ones(size(ratio))./ratio;
delta = 1 + (1* repmat(0,size(eDis)) .* ratio);
    
hold on 

if numel(group) <= 15

       n = size(vec,1); % get # of variables
   for j = 1:n
      % Plot vectors:
      f_arrow([0 0],[vec(j,1) vec(j,2)],'size',0.125*0.75,'angle',20,'Color','k')
      
      % Label vectors:
      h = text(vec(j,1)*delta(j),vec(j,2)*delta(j),group(j),'fontsize',8);
      set(h,'FontSize',8,'HorizontalAlignment','center','Color','k','Interpreter','tex');
   end

elseif numel(group) > 15

           n = size(vec,1); % get # of variables

distances = sqrt(vec(:,1).^2 + vec(:,2).^2);
[~, sortedIndices] = sort(distances, 'descend');
topIndices = sortedIndices(1:round(0.3 * numel(distances)));

for j=topIndices'
    if distances(j)
      % Plot vectors:
      f_arrow([0 0],[vec(j,1) vec(j,2)],'size',0.125*0.75,'angle',20,'Color','k')
      
      % Label vectors:
      h = text(vec(j,1)*delta(j),vec(j,2)*delta(j),group(j),'fontsize',8);
      set(h,'FontSize',8,'HorizontalAlignment','center','Color','k','Interpreter','tex');
    end
end

end

axis1 = sprintf('%2.2f',expl(1,1));
axis2 = sprintf('%2.2f',expl(2,1));
xlabel(['PCoA 1 (' num2str(axis1) ' % variance)'],'fontsize',10);
ylabel(['PCoA 2 (' num2str(axis2) ' % variance)'],'fontsize',10);

%Map with RGB colors

subplot('Position', [0.4 0.1 0.55 0.8])

m_proj('mercator','lon',[lonmin lonmax],'lat',[latmin latmax]);
m_coast('patch',landcolor,'edgecolor',[0.6 0.6 0.6]);
hold on 
m_grid('fontsize',8,'linestyle','none','meridianlabel','off','parallellabel','off','xtick',[],'ytick',[]);

[n,m] = size(data_stat);

m_scatter(lon,lat, 50, RGB, 'filled', 'markeredgecolor', [0.5 0.5 0.5]);

%Export and save? 
if save==false;
elseif save==true;
        outputDir = fullfile(parentDir, subDir);
            mkdir(outputDir);

    filename = fullfile(outputDir, sprintf(['PCoA and RGB color on map ' char(dataACP_name) ' ' list3{normalisation} '.pdf']));
    exportgraphics(gcf, filename, 'ContentType', 'vector', 'BackgroundColor', 'none');
end

pause(1)
close gcf

%% CA 

%to save
subDir = 'CA';
fullSubDirPath = fullfile(parentDir, subDir);
    mkdir(fullSubDirPath);

% CA (inspired from Rencher, A. C. (2000))

%  To cite this file, this would be an appropriate format:
%  Trujillo-Ortiz, A., R. Hernandez-Walls and K. Barba-Rojo. (2008). corran:
%    Correspondence Analysis. A MATLAB file. [WWW document]. URL http://
%    www.mathworks.com/matlabcentral/fileexchange/loadFile.do?objectId=21505
%  Reference:
%  Rencher, A. C. (2000), Methods of Multivariate Analysis. 2nd ed. 
%       John Wiley & Sons. p. 514-525

[a, b] = size(data_stat); 
N = sum(data_stat(:)); % total of the matrix
P = data_stat / N; % Correspondence matrix (relative frequencies)
r = sum(P, 2); % Row totals (marginal frequencies)
c = sum(P, 1); % Column totals (marginal frequencies)
Dr = diag(r); % Diagonal matrix of row totals
Dc = diag(c); % Diagonal matrix of column totals
R = inv(Dr) * P; % Row profiles (standardized)
C = P * inv(Dc); % Column profiles (standardized)
Z = sqrt(inv(Dr)) * (P - r * c) * sqrt(inv(Dc)); % Matrix for SVD
[U, S, V] = svd(Z); % Singular Value Decomposition
k = rank(inv(Dr) * (P - r * c) * inv(Dc) * (P - r * c)'); % Maximum number of dimensions
SV = diag(S); %singular values
SV = SV(1:k); %singular value equal to zero is deleted
EV = SV.^2; %eigenvalues (inertia)
PI = EV./sum(EV)*100; %percent of inertia
PA = cumsum(PI); %cumulative percent of inertia
x2 = N*EV; %chi-squared
A = sqrt(Dr)*U;
B = sqrt(Dc)*V;
X2 = N*trace(inv(Dr)*(P-r*c)*inv(Dc)*(P-r*c)'); %total chi-squared
v = (a-1)*(b-1); %degrees of freedom
P = 1-chi2cdf(X2,v); %P-value

% Scores for rows and columns
X=inv(Dr)*A*S; %row coordinates
X=X(:,1:k);
Y=inv(Dc)*B*S'; %column coordinates
Y=Y(:,1:k);

% CA plot:

figure('numbertitle','off','Position',[10 50 scrsz(3) scrsz(4)-150],'Color','white');

% distances
distances_X = sqrt(X(:,1).^2 + X(:,2).^2);
distances_Y = sqrt(Y(:,1).^2 + Y(:,2).^2);

if numel(group) <= 15

    %scatter(X(:,1), X(:,2), 50, 'filled', 'MarkerFaceColor', [0.9 0.9 0.9], 'MarkerEdgeColor', [0.5 0.5 0.5]);
    text(X(:,1), X(:,2), group', 'Color', 'k', 'fontsize', 8);

    hold on
    for i=1:size(X,1)
    plot([0 X(i,1)],[0 X(i,2)],'k-')
    end

elseif numel(group) > 15

    [~, sortedIndices_X] = sort(distances_X, 'descend');
    topIndices_X = sortedIndices_X(1:round(0.3 * numel(distances_X)));
    
    %scatter(X(topIndices_X, 1), X(topIndices_X, 2), 50, 'filled', 'MarkerFaceColor', [0.9 0.9 0.9], 'MarkerEdgeColor', [0.5 0.5 0.5]);
    text(X(topIndices_X, 1), X(topIndices_X, 2), group(topIndices_X), 'Color', 'k', 'fontsize', 8);

    hold on
    for i = topIndices_X'
        plot([0 X(i,1)], [0 X(i,2)], 'k-');
    end

end

hold on;

scatter(Y(:,1), Y(:,2), 30, 'filled', 'MarkerFaceColor', [0.9 0.9 0.9], 'MarkerEdgeColor', [0.5 0.5 0.5]);

xlim([-max(abs([X(:,1); Y(:,1)])) * 1.1, max(abs([X(:,1); Y(:,1)])) * 1.1]);
ylim([-max(abs([X(:,2); Y(:,2)])) * 1.1, max(abs([X(:,2); Y(:,2)])) * 1.1]);

%text(Y(:,1) + .01, Y(:,2) + .1, sampleid', 'Color', 'k', 'fontsize', 8);

xlabel(['CA 1 (' num2str(fix(PI(1))) '% variance)']);
ylabel(['CA 2 (' num2str(fix(PI(2))) '% variance)']);

%Export and save? 
if save==false;
elseif save==true;
    outputDir = fullfile(parentDir, subDir);
    mkdir(outputDir);
    
    filename = fullfile(outputDir, sprintf(['CA on ' char(dataACP_name) ' ' list3{normalisation} '.pdf']));
    exportgraphics(gcf, filename, 'ContentType', 'vector', 'BackgroundColor', 'none');
end

pause(1);
close(gcf);

% Plot CA with RGB colors

%Representation with RGB colors 
%Extract RGB color on the CA

    RGB=Y(:,1:3);
    minRGB=min(RGB);
    RGB=RGB-ones(size(data_stat,2),1)*minRGB;
    maxRGB=max(RGB);
    RGB=RGB./(ones(size(data_stat,2),1)*maxRGB);
    RGB=1-RGB;

figure('numbertitle','off','Position',[10 50 scrsz(3) scrsz(4)-150],'Color','white');

subplot('Position', [0.05 0.1 0.3 0.8])

% distances
distances_X = sqrt(X(:,1).^2 + X(:,2).^2);
distances_Y = sqrt(Y(:,1).^2 + Y(:,2).^2);

if numel(group) <= 15

    %scatter(X(:,1), X(:,2), 30, 'filled', 'MarkerFaceColor', [0.9 0.9 0.9], 'MarkerEdgeColor', [0.5 0.5 0.5]);
    text(X(:,1), X(:,2), group', 'Color', 'k', 'fontsize', 8);

    hold on
    for i=1:size(X,1)
    plot([0 X(i,1)],[0 X(i,2)],'k-')
    end

else

    [~, sortedIndices_X] = sort(distances_X, 'descend');
    topIndices_X = sortedIndices_X(1:round(0.3 * numel(distances_X)));
    
    %scatter(X(topIndices_X, 1), X(topIndices_X, 2), 30, 'filled', 'MarkerFaceColor', [0.9 0.9 0.9], 'MarkerEdgeColor', [0.5 0.5 0.5]);
    text(X(topIndices_X, 1), X(topIndices_X, 2), group(topIndices_X), 'Color', 'k', 'fontsize', 8);

    hold on
    for i = topIndices_X'
        plot([0 X(i,1)], [0 X(i,2)], 'k-');
    end

end

hold on;

scatter(Y(:,1), Y(:,2), 50, RGB, 'filled', 'MarkerEdgeColor', [0.5 0.5 0.5]);

xlim([-max(abs([X(:,1); Y(:,1)])) * 1.1, max(abs([X(:,1); Y(:,1)])) * 1.1]);
ylim([-max(abs([X(:,2); Y(:,2)])) * 1.1, max(abs([X(:,2); Y(:,2)])) * 1.1]);

%text(Y(:,1) + .01, Y(:,2) + .1, sampleid', 'Color', 'k', 'fontsize', 8);

xlabel(['CA 1 (' num2str(fix(PI(1))) '% variance)']);
ylabel(['CA 2 (' num2str(fix(PI(2))) '% variance)']);

%Map with RGB colors

subplot('Position', [0.4 0.1 0.55 0.8])

m_proj('mercator','lon',[lonmin lonmax],'lat',[latmin latmax]);
m_coast('patch',landcolor,'edgecolor',[0.6 0.6 0.6]);
hold on 
m_grid('fontsize',8,'linestyle','none','meridianlabel','off','parallellabel','off','xtick',[],'ytick',[]);

[n,m] = size(data_stat);

m_scatter(lon,lat, 50, RGB, 'filled', 'markeredgecolor', [0.5 0.5 0.5]);

%Export and save? 
if save==false;
elseif save==true;
    outputDir = fullfile(parentDir, subDir);
    mkdir(outputDir);
    
    filename = fullfile(outputDir, sprintf(['CA and RGB color on map on ' char(dataACP_name) ' ' list3{normalisation} '.pdf']));
    exportgraphics(gcf, filename, 'ContentType', 'vector', 'BackgroundColor', 'none');
end

pause(1);
close(gcf);

%% Clusters

%to save
subDir = 'Clusters analysis';
fullSubDirPath = fullfile(parentDir, subDir);
    mkdir(fullSubDirPath);

% Dendrogram to define cut

% Dendrocut value suggestion
Z = linkage(data_stat', 'complete', 'euclidean'); 
incon=inconsistent(Z);
test=cumsum(flipud(incon(:,4)));
I=find(test==max(test));
dendrocut_auto=incon(I(1),1); %Suggestion

% Plot dendrogram with dendrocut_auto

if ~iscell(sampleid); sampleid = num2cell(sampleid); end
if ischar(sampleid) || isstring(sampleid); sampleid = cellstr(sampleid); end

figure('numbertitle','off','Position',[10 50 scrsz(3) scrsz(4)-150],'Color','white');

subplot(1,2,1)
plot(incon(:,1),cumsum(flipud(incon(:,4))),'Color', [0.5 0.5 0.5])
hold on
plot(incon(I(1),1),cumsum(flipud(incon(I(1),4))),'o','MarkerFaceColor', [0.93,0.69,0.13])
text(incon(I(1),1), cumsum(flipud(incon(I(1),4)))+5, sprintf('%.2f', dendrocut_auto), 'Color', [0.93,0.69,0.13], 'HorizontalAlignment', 'left');

dis2=pdist(data_stat', 'euclidean');
leafOrder = optimalleaforder(Z,dis2);

subplot(1,2,2)
[lineH,T, Perm] = dendrogram(Z,0, 'ColorThreshold',dendrocut_auto,'reorder',leafOrder,'Labels',sampleid);
colors = cell2mat(get(lineH,'Color'));
set(lineH,'Color', [0.5 0.5 0.5]);
set(gca, 'XTick', []);

hold on;
hline = refline(0,dendrocut_auto);
hline.Color = [0.93,0.69,0.13]; 

color=slandarerCM(6).Colors(1,9); %spectral 
color=color{1, 1};

%Export and save? 
if save==false;
elseif save==true;

        outputDir = fullfile(parentDir, subDir);
        mkdir(outputDir);

    filename = fullfile(outputDir, sprintf(['Dendrogram for clustering on ' char(dataACP_name) ' ' list3{normalisation} '.pdf']));
    exportgraphics(gcf, filename, 'ContentType', 'vector', 'BackgroundColor', 'none');
end
pause(1)
close gcf

% Or set the dendrocut value manually?
list2 = {'Keep', 'Change'};  
[dendro,tf] = listdlg('PromptString','Keep automatic dendrocut or select manually?','ListString',list2,'SelectionMode','single','ListSize',[400,200]);

if dendro==1;
    dendrocut = dendrocut_auto; 
elseif dendro==2;
        name = 'Dendrocut value';
        prompt = {'Enter the Dendrocut value:'};
        dims = [1 35];
        definput = {'10'};
        answerfrac1 = inputdlg(prompt, name, dims, definput);
        dendrocut = str2double(answerfrac1{1});
end

close gcf

% Create and Plot Clustergram

CGobj=clustergram(data_stat,'RowLabels',group','ColumnLabels',sampleid,'Linkage','complete','ColumnPDist', 'euclidean','ColumnPDist', 'euclidean','Colormap',flipud(color),'Standardize',3,'Symmetric','false','DisplayRatio',[0.2 0.000000000000001],'OptimalLeafOrder',1);
CGobj.Dendrogram = dendrocut;
CGobj.DisplayRange = max(max(data_stat));
A=CGobj.ColumnLabels;
B=0;
close force 

plot(CGobj)

%Export and save? 
if save==false;
elseif save==true;

        outputDir = fullfile(parentDir, subDir);
        mkdir(outputDir);

    filename = fullfile(outputDir, sprintf(['Clustergram on ' char(dataACP_name) ' ' list3{normalisation} '.pdf']));
    exportgraphics(gcf, filename, 'ContentType', 'vector', 'BackgroundColor', 'none');
end
pause(1)
close gcf

% Create clusters

T2 = cluster(Z, 'cutoff', dendrocut, 'criterion', 'distance');

uniqueT2 = unique(T2); 

uniqueT2 = unique(T2); 
palette_couleurs = parula(numel(uniqueT2));

for i = 1:numel(T2)
    indice_cluster = T2(i);
    couleur_cluster = palette_couleurs(indice_cluster, :);
    T2(i, 2:4) = couleur_cluster; 
end

% PCA plot with clusters

figure('numbertitle','off','Position',[10 50 scrsz(3) scrsz(4)-150],'Color','white');

subplot('Position', [0.05 0.1 0.3 0.8])

scatter(zscores(:,1),zscores(:,2),50,T2(:,2:4),'filled','markeredgecolor',[0.5 0.5 0.5])
[p,q]=size(pc);

hold on

if numel(group) <= 15

[p,q]=size(pc);
hold on
for i=1:p
    if ((pc(i,1)^2+pc(i,2)^2)^0.5)
    plot([0 pc(i,1)],[0 pc(i,2)],'k-')
    text(pc(i,1),pc(i,2),group(i),'fontsize',8);
    end
end

elseif numel(group) > 15

[p,q]=size(pc);
hold on

distances = sqrt(pc(:,1).^2 + pc(:,2).^2);
[~, sortedIndices] = sort(distances, 'descend');
topIndices = sortedIndices(1:round(0.3 * numel(distances)));

for i=topIndices'
    if distances(i)
    plot([0 pc(i,1)],[0 pc(i,2)],'k-')
    text(pc(i,1),pc(i,2),group(i),'fontsize',8);
    end
end

end

xlabel(['PCA 1 ' num2str(fix(pcvars(1))) '% variance'],'fontsize',10);
ylabel(['PCA 2 ' num2str(fix(pcvars(2))) '% variance'],'fontsize',10);

% Clusters on the map 

subplot('Position', [0.4 0.1 0.55 0.8])

m_proj('mercator','lon',[lonmin lonmax],'lat',[latmin latmax]);
m_coast('patch',landcolor,'edgecolor',[0.6 0.6 0.6]);
hold on 
m_grid('fontsize',8,'linestyle','none','meridianlabel','off','parallellabel','off','xtick',[],'ytick',[]);

[n,m] = size(data_stat);
hold on;
for i = 1:m
    m_scatter(lon(i),lat(i), 50, T2(i,2:4), 'filled', 'markeredgecolor', [0.5 0.5 0.5]);
end

%Export and save? 
if save==false;
elseif save==true;

        outputDir = fullfile(parentDir, subDir);
            mkdir(outputDir);

    filename = fullfile(outputDir, sprintf(['PCA and Map Clusters on ' char(dataACP_name) ' ' list3{normalisation} '.pdf']));
    exportgraphics(gcf, filename, 'ContentType', 'vector', 'BackgroundColor', 'none');
end

pause(1)
close gcf

% PCoA plot with clusters

figure('numbertitle','off','Position',[10 50 scrsz(3) scrsz(4)-150],'Color','white');

subplot('Position', [0.05 0.1 0.3 0.8])

scatter(scores(:,1), scores(:,2), 50, T2(:,2:4), 'filled', 'markeredgecolor', [0.5 0.5 0.5]);
%text(scores(:,1),scores(:,2),sampleid,'HorizontalAlignment','center', 'VerticalAlignment','middle','Color','k','Interpreter','tex'); 
    
hold on 

if numel(group) <= 15

       n = size(vec,1); % get # of variables
   for j = 1:n
      % Plot vectors:
      f_arrow([0 0],[vec(j,1) vec(j,2)],'size',0.125*0.75,'angle',20,'Color','k')
      
      % Label vectors:
      h = text(vec(j,1)*delta(j),vec(j,2)*delta(j),group(j),'fontsize',8);
      set(h,'FontSize',8,'HorizontalAlignment','center','Color','k','Interpreter','tex');
   end

elseif numel(group) > 15

           n = size(vec,1); % get # of variables

distances = sqrt(vec(:,1).^2 + vec(:,2).^2);
[~, sortedIndices] = sort(distances, 'descend');
topIndices = sortedIndices(1:round(0.3 * numel(distances)));

for j=topIndices'
    if distances(j)
      % Plot vectors:
      f_arrow([0 0],[vec(j,1) vec(j,2)],'size',0.125*0.75,'angle',20,'Color','k')
      
      % Label vectors:
      h = text(vec(j,1)*delta(j),vec(j,2)*delta(j),group(j),'fontsize',8);
      set(h,'FontSize',8,'HorizontalAlignment','center','Color','k','Interpreter','tex');
    end
end

end

axis1 = sprintf('%2.2f',expl(1,1));
axis2 = sprintf('%2.2f',expl(2,1));
xlabel(['PCoA 1 (' num2str(axis1) ' % variance)'],'fontsize',10);
ylabel(['PCoA 2 (' num2str(axis2) ' % variance)'],'fontsize',10);

% Clusters on the map 

subplot('Position', [0.4 0.1 0.55 0.8])

m_proj('mercator','lon',[lonmin lonmax],'lat',[latmin latmax]);
m_coast('patch',landcolor,'edgecolor',[0.6 0.6 0.6]);
hold on 
m_grid('fontsize',8,'linestyle','none','meridianlabel','off','parallellabel','off','xtick',[],'ytick',[]);

[n,m] = size(data_stat);
hold on;
for i = 1:m
    m_scatter(lon(i),lat(i), 50, T2(i,2:4), 'filled', 'markeredgecolor', [0.5 0.5 0.5]);
end

%Export and save? 
if save==false;
elseif save==true;

        outputDir = fullfile(parentDir, subDir);
            mkdir(outputDir);

    filename = fullfile(outputDir, sprintf(['PCoA and Map Clusters on ' char(dataACP_name) ' ' list3{normalisation} '.pdf']));
    exportgraphics(gcf, filename, 'ContentType', 'vector', 'BackgroundColor', 'none');
end

pause(1)
close gcf

% CA plot with clusters

figure('numbertitle','off','Position',[10 50 scrsz(3) scrsz(4)-150],'Color','white');

subplot('Position', [0.05 0.1 0.3 0.8])

% distances
distances_X = sqrt(X(:,1).^2 + X(:,2).^2);
distances_Y = sqrt(Y(:,1).^2 + Y(:,2).^2);

if numel(group) <= 15

    %scatter(X(:,1), X(:,2), 30, 'filled', 'MarkerFaceColor', [0.9 0.9 0.9], 'MarkerEdgeColor', [0.5 0.5 0.5]);
    text(X(:,1), X(:,2), group', 'Color', 'k', 'fontsize', 8);

    hold on
    for i=1:size(X,1)
    plot([0 X(i,1)],[0 X(i,2)],'k-')
    end

else

    [~, sortedIndices_X] = sort(distances_X, 'descend');
    topIndices_X = sortedIndices_X(1:round(0.3 * numel(distances_X)));
    
    %scatter(X(topIndices_X, 1), X(topIndices_X, 2), 30, 'filled', 'MarkerFaceColor', [0.9 0.9 0.9], 'MarkerEdgeColor', [0.5 0.5 0.5]);
    text(X(topIndices_X, 1), X(topIndices_X, 2), group(topIndices_X), 'Color', 'k', 'fontsize', 8);

    for i=topIndices'
    if distances(i)
    plot([0 pc(i,1)],[0 pc(i,2)],'k-')
    text(pc(i,1),pc(i,2),group(i),'fontsize',8);
    end
    end

end

hold on;

scatter(Y(:,1), Y(:,2), 50, T2(:,2:4), 'filled', 'MarkerEdgeColor', [0.5 0.5 0.5]);

xlim([-max(abs([X(:,1); Y(:,1)])) * 1.1, max(abs([X(:,1); Y(:,1)])) * 1.1]);
ylim([-max(abs([X(:,2); Y(:,2)])) * 1.1, max(abs([X(:,2); Y(:,2)])) * 1.1]);

%text(Y(:,1) + .01, Y(:,2) + .1, sampleid', 'Color', 'k', 'fontsize', 8);

xlabel(['CA 1 (' num2str(fix(PI(1))) '% variance)']);
ylabel(['CA 2 (' num2str(fix(PI(2))) '% variance)']);

%Map with RGB colors

subplot('Position', [0.4 0.1 0.55 0.8])

m_proj('mercator','lon',[lonmin lonmax],'lat',[latmin latmax]);
m_coast('patch',landcolor,'edgecolor',[0.6 0.6 0.6]);
hold on 
m_grid('fontsize',8,'linestyle','none','meridianlabel','off','parallellabel','off','xtick',[],'ytick',[]);

[n,m] = size(data_stat);

[n,m] = size(data_stat);
hold on;
for i = 1:m
    m_scatter(lon(i),lat(i), 50, T2(i,2:4), 'filled', 'markeredgecolor', [0.5 0.5 0.5]);
end

%Export and save? 
if save
    outputDir = fullfile(parentDir, subDir);
    mkdir(outputDir);
    
    filename = fullfile(outputDir, sprintf(['CA and Map Clusters on ' char(dataACP_name) ' ' list3{normalisation} '.pdf']));
    exportgraphics(gcf, filename, 'ContentType', 'vector', 'BackgroundColor', 'none');
end

pause(1)
close gcf

%% Analysis on the NBSS 

%to save
parentDir = 'NBSS analysis';
fullSubDirPath = fullfile(pwd, parentDir);  
    mkdir(fullSubDirPath);

%% NBSS Size Spectra 

% NBSS size spectra 

figure('numbertitle','off','Position',[10 50 scrsz(3) scrsz(4)-150],'Color','white');

for i = 1:numel(base)
    hold on
    plot(ESD, base(i).regroupped.NBSSliving, 'o', 'MarkerFaceColor', [0.9, 0.9, 0.9], 'MarkerEdgeColor', [0.5, 0.5, 0.5]);
end 

set(gca, 'YScale', 'log', 'XScale', 'log');


xticks([0.1 0.2 0.5 1 2 5 10 20 50 100 200 500 1000 2000 5000 10000 20000 50000 100000]);
xticklabels({'0.1', '0.2', '0.5', '1', '2', '5', '10', '20', '50', '100', '200', '500', '1000', '2000', '5000', '10000', '20000', '50000', '100000'});

xlabel('size on ESD (µm)', 'fontsize', 10);
ylabel({'NBSS (mm^3/mm^3/m^3)'; 'equivalent abundance'}, 'fontsize', 10);
title('NBSS of all stations', 'fontsize', 10);

%Export and save? 
if save==false;
elseif save==true;

        outputDir = fullfile(parentDir);
            mkdir(outputDir);

    filename = fullfile(outputDir, sprintf('NBSS Spectra.pdf'));
    exportgraphics(gcf, filename, 'ContentType', 'vector', 'BackgroundColor', 'none');
end

pause(1)
close gcf

% NBSS size spectra trunked

figure('numbertitle','off','Position',[10 50 scrsz(3) scrsz(4)-150],'Color','white');

for i = 1:numel(base)
    hold on
    plot(ESD, base(i).regroupped.NBSSliving_trunked, 'o', 'MarkerFaceColor', [0.9, 0.9, 0.9], 'MarkerEdgeColor', [0.5, 0.5, 0.5]);
end 

set(gca, 'YScale', 'log', 'XScale', 'log');

xticks([0.1 0.2 0.5 1 2 5 10 20 50 100 200 500 1000 2000 5000 10000 20000 50000 100000]);
xticklabels({'0.1', '0.2', '0.5', '1', '2', '5', '10', '20', '50', '100', '200', '500', '1000', '2000', '5000', '10000', '20000', '50000', '100000'});

xlabel('size on ESD (µm)', 'fontsize', 10);
ylabel({'NBSS (mm^3/mm^3/m^3)'; 'equivalent abundance'}, 'fontsize', 10);
title('NBSS of all stations', 'fontsize', 10);

%Export and save? 
if save==false;
elseif save==true;

        outputDir = fullfile(parentDir);
        mkdir(outputDir);

    filename = fullfile(outputDir, sprintf('NBSS Spectra trunked.pdf'));
    exportgraphics(gcf, filename, 'ContentType', 'vector', 'BackgroundColor', 'none');
end

pause(1)
close gcf

% NBSS Slopes 

for i = 1:numel(base);
slope_values(i) = [base(i).regroupped.NBSSliving_slope];
end

% PCA 

figure('numbertitle','off','Position',[10 50 scrsz(3) scrsz(4)-150],'Color','white');

scatter(zscores(:,1), zscores(:,2), 50, slope_values, 'filled', 'MarkerEdgeColor', [0.5 0.5 0.5]); 

[p,q]=size(pc);

hold on

if numel(group) <= 15

[p,q]=size(pc);
hold on
for i=1:p
    if ((pc(i,1)^2+pc(i,2)^2)^0.5)
    plot([0 pc(i,1)],[0 pc(i,2)],'k-')
    text(pc(i,1),pc(i,2),group(i),'fontsize',8);
    end
end

elseif numel(group) > 15

[p,q]=size(pc);
hold on

distances = sqrt(pc(:,1).^2 + pc(:,2).^2);
[~, sortedIndices] = sort(distances, 'descend');
topIndices = sortedIndices(1:round(0.3 * numel(distances)));

for i=topIndices'
    if distances(i)
    plot([0 pc(i,1)],[0 pc(i,2)],'k-')
    text(pc(i,1),pc(i,2),group(i),'fontsize',8);
    end
end

end

xlabel(['PCA 1 ' num2str(fix(pcvars(1))) '% variance'],'fontsize',10);
ylabel(['PCA 2 ' num2str(fix(pcvars(2))) '% variance'],'fontsize',10);

color=slandarerCM(6).Colors(1,9); %spectral (change manually)
color=color{1, 1};
colormap(color);
colormap(flipud(color));
c2 = colorbar;

c=colorbar;
title(c,'NBSS Slope','fontweight','bold');

%annotation('ellipse', [0.65 0.025 0.01 0.015], 'EdgeColor', [0.5 0.5 0.5], 'LineWidth', 0.5, 'FaceColor', 'none');
%annotation('textbox', [0.66 0.02 0.6 0.04], 'String', '< 5 points on NBSS, slope not reliable', 'EdgeColor', 'none', 'FontSize', 8, 'Color', 'k');

%Export and save? 
if save==false;
elseif save==true;

        outputDir = fullfile(parentDir);
        mkdir(outputDir);

    filename = fullfile(outputDir, sprintf(['PCA with NBSS slopes on ' char(dataACP_name) '.pdf']));
    exportgraphics(gcf, filename, 'ContentType', 'vector', 'BackgroundColor', 'none');

end

pause(1)
close gcf

% % PCoA
% 
% figure('numbertitle','off','Position',[10 50 scrsz(3) scrsz(4)-150],'Color','white');
% 
% scatter(scores(:,1), scores(:,2), 50, slope_values, 'filled', 'markeredgecolor', [0.5 0.5 0.5]);
% 
% hold on 
% 
% if numel(group) <= 10
% 
%        n = size(vec,1); % get # of variables
%    for j = 1:n
%       % Plot vectors:
%       f_arrow([0 0],[vec(j,1) vec(j,2)],'size',0.125*0.75,'angle',20,'Color','k')
%       
%       % Label vectors:
%       h = text(vec(j,1)*delta(j),vec(j,2)*delta(j),group(j),'fontsize',8);
%       set(h,'FontSize',8,'HorizontalAlignment','center','Color','k','Interpreter','tex');
%    end
% 
% elseif numel(group) > 10
% 
%            n = size(vec,1); % get # of variables
% 
% distances = sqrt(vec(:,1).^2 + vec(:,2).^2);
% [~, sortedIndices] = sort(distances, 'descend');
% topIndices = sortedIndices(1:round(0.3 * numel(distances)));
% 
% for j=topIndices'
%     if distances(j)
%       % Plot vectors:
%       f_arrow([0 0],[vec(j,1) vec(j,2)],'size',0.125*0.75,'angle',20,'Color','k')
%       
%       % Label vectors:
%       h = text(vec(j,1)*delta(j),vec(j,2)*delta(j),group(j),'fontsize',8);
%       set(h,'FontSize',8,'HorizontalAlignment','center','Color','k','Interpreter','tex');
%     end
% end
% 
% end
% 
% color=slandarerCM(6).Colors(1,9); %spectral (change manually)
% color=color{1, 1};
% colormap(color);
% colormap(flipud(color));
% c2 = colorbar;
% 
% axis1 = sprintf('%2.2f',expl(1,1));
% axis2 = sprintf('%2.2f',expl(2,1));
% xlabel(['PCoA 1 (' num2str(axis1) ' % variance)'],'fontsize',10);
% ylabel(['PCoA 2 (' num2str(axis2) ' % variance)'],'fontsize',10);
% 
% c=colorbar;
% title(c,'NBSS Slope','fontweight','bold');
% 
% %annotation('ellipse', [0.65 0.025 0.01 0.015], 'EdgeColor', [0.5 0.5 0.5], 'LineWidth', 0.5, 'FaceColor', 'none');
% %annotation('textbox', [0.66 0.02 0.6 0.04], 'String', '< 5 points on NBSS, slope not reliable', 'EdgeColor', 'none', 'FontSize', 8, 'Color', 'k');
% 
% %Export and save? 
% if save==false;
% elseif save==true;
% 
%         outputDir = fullfile(parentDir);
%             mkdir(outputDir);
% 
%     filename = fullfile(outputDir, sprintf(['PCoA with NBSS slopes on ' char(dataACP_name) '.pdf']));
%     exportgraphics(gcf, filename, 'ContentType', 'vector', 'BackgroundColor', 'none');
% end
% pause(1)
% close gcf

% Map 

figure('numbertitle','off','Position',[10 50 scrsz(3) scrsz(4)-150],'Color','white');

m_proj('mercator','lon',[lonmin lonmax],'lat',[latmin latmax]);
m_coast('patch',landcolor,'edgecolor',[0.6 0.6 0.6]);
hold on 
m_grid('fontsize',8,'linestyle','none','meridianlabel','off','parallellabel','off','xtick',[],'ytick',[]);

hold on;
m_scatter(lon,lat, 50, slope_values, 'filled', 'markeredgecolor', [0.5 0.5 0.5]);

color=slandarerCM(6).Colors(1,9); %spectral (change manually)
color=color{1, 1};
colormap(color);
colormap(flipud(color));
c2 = colorbar;

c=colorbar;
title(c,'NBSS Slope','fontweight','bold');
%annotation('textbox', [0.7 0.05 0.2 0.1], 'String', 'empty if points nb on NBSS are < 5 : slope not reliable', 'EdgeColor', 'none', 'FontSize', 8, 'Color', 'k');

%Export and save? 
if save==false;
elseif save==true;

        outputDir = fullfile(parentDir);
        mkdir(outputDir);

    filename = fullfile(outputDir, sprintf(['Map of NBSS Slope.pdf']));
    exportgraphics(gcf, filename, 'ContentType', 'vector', 'BackgroundColor', 'none');
end

pause(1)
close gcf