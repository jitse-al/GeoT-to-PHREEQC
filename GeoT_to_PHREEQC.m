%This is the main file for conversion of GeoT/TOUGHREACT input database files. Please specify your pathname, filename and filetype here. The converted file will be put in the same path as the databasefile.


%INPUT 1. Adapt this with the path-/filename that you want to convert
pathname= 'C:\Users\Jitse\Documents\University\PhD - Zaragoza\Geot\V2.1\Thermo-databases\';
filename= 'tk-slt.h06_jun16';
%filename= 'tk1-ympR5';
filetype= '.dat';




%Input 2. 
ignoredoublenames= 'n'; %Either y or n. Defines whether double reactions (i.e., redox reactions that have been defined twice) are separated or not. Separation: 'n'/no separation: 'y'
%Preferred species defines the redox pair converted to the database. Uncomment the lines that you want, comment the rest with "%". Other redox species can be used as well, if present in databases
preferredspecies={ 'O'};
preferredspecies2={ '''O2'''};
preferredspecies={ 'S'};
preferredspecies2={ 'HS-'};








%DEFINITIONS
%Code rewrites Chemical formula, necessary for definition of reactions. This is the order of elements necessary to provide the molecular formula. IUPAC-based (roughly, may not always correspond to commonly used mineral formulas)
elementalorder={ 'Rn' 'Xe' 'Kr' 'Ar' 'Ne' 'He' 'Fr' 'Cs' 'Rb' 'K' 'Na' 'Li' 'Ra' 'Ba' 'Sr' 'Ca' 'Mg' 'Be' 'Lr' 'Ac' 'Lu' 'La' 'Y' 'Se' 'Rf' 'Hf' 'Zr' 'Ti' 'Db' 'Ta' 'Nb' 'V' 'Sg' 'W' 'Mo' 'Cr' 'Bh' 'Re' 'Te' 'Mn' 'Hs' 'Os' 'Ru' 'Fe' 'Mt' 'Ir' 'Rh' 'Co' 'Ds' 'Pt' 'Pd' 'Ni' 'Rg' 'Au' 'Ag' 'Cu' 'Hg' 'Cd' 'Zn' 'Tl' 'In' 'Ga' 'Al' 'B' 'Pb' 'Sn' 'Ge' 'Si' 'C' 'Bi' 'Sb' 'As' 'P' 'N' 'H' 'Po' 'Te' 'Se' 'S' 'O' '(OH)' 'At' 'I' 'Br' 'Cl' 'F'};

%Selection of elemental weights necessary for PRIMARY SPECIES block in PHREEQC. Based on databases present in PHREEQC, version 3.4.0
solution_master_species_matrix={ 'Acetate','Acetate-','1','59.045','59.045';'Acetate','HAcetate','0','Acetate','59';'Ag','Ag+','0','107.868','107.868';'Ag','Ag++','0','Ag','107.868';'Al','Al+++','0','26.9815','26.9815';'Am','Am+++','0','Am','243';'Am','Am++','0','Am','243';'Am','Am++++','0','Am','243';'Am','AmO2+','0','Am','243';'Am','AmO2++','0','Am','243';'Amm','AmmH+','0','AmmH','17';'Ar','Ar','0','Ar','39.948';'As','AsO4---','0','As','74.9216';'As','H2AsO4-','0','As','74.9216';'As','H3AsO4','-1','74.9216','74.9216';'As','H2AsO3-','0','As','74.9216';'As','H3AsO3','0','74.9216','74.9216';'As','AsH3','0','As','74.9216';'Au','Au+','0','Au','196.9665';'Au','Au+++','0','Au','196.9665';'B','B(OH)3','0','B','10.81';'B','B(OH)4-','0','B','10.811';'B','H3BO3','0','10.81','10.81';'B','BH4-','0','B','10.811';'Ba','Ba++','0','Ba','137.327';'Be','Be++','0','9.0122','9.0122';'Benzoate','Benzoate-','0','121.116','121.116';'Br','Br-','0','79.904','79.904';'Br','Br2','0','Br','79.904';'Br(-03)','Br3-','0','Br','79.904';'Br','BrO-','0','Br','79.904';'Br','BrO3-','0','Br','79.904';'Br','BrO4-','0','Br','79.904';'Butanoate','Butanoate-','0','87.043','87.043';'Butylamine','Butylamine','1','73.138','73.138';'Butyrate','Butyrate-','1','87.098','87.098';'C','CO2','0','HCO3','12.0111';'C','CO3--','2','CO3','12.0111';'C','HCO3-','1','HCO3','12.0111';'C','CO','0','C','12.0111';'C','C2H4','0','C2H4','12.0111';'C','C2H6','0','C2H6','12.0111';'C','CH4','0','CH4','12.011';'Ca','Ca++','0','40.078','40.078';'Cd','Cd++','0','112.399','112.399';'Ce','Ce+++','0','Ce','140.115';'Ce','Ce++','0','Ce','140.115';'Ce','Ce++++','0','Ce','140.115';'Cit','Cit---','0','Cit','189.1013';'Citrate','Citrate---','2','189.102','189.102';'Cl','Cl-','0','35.453','35.453';'Cl','ClO-','0','Cl','35.453';'Cl','ClO2-','0','Cl','35.453';'Cl','ClO3-','0','Cl','35.453';'Cl','ClO4-','0','Cl','35.453';'Cm','Cm+++','0','Cm','247';'Cn','Cn-','0','Cn','26.018';'Co','Co++','0','Co','58.9332';'Co','Co+++','-1','Co','58.9332';'Cr','CrO4--','0','CrO4','51.9961';'Cr','Cr++','0','Cr','51.9961';'Cr','Cr(OH)2+','1','51.996','51.9961';'Cr','Cr+++','0','Cr','51.9961';'Cr','CrO4---','0','Cr','51.9961';'Cs','Cs+','0','Cs','132.9054';'Cu','Cu++','0','63.546','63.546';'Cu','Cu+','0','Cu','63.546';'Cyanate','Cyanate-','0','42.017','42.017';'Cyanide','Cyanide-','1','Cyanide','26.0177';'Diethylamine','Diethylamine','1','73.138','73.138';'Dimethylamine','Dimethylamine','1','45.028','45.028';'Dom_a','Dom_a','0','C','12.0111';'Dom_b','Dom_b','0','C','12.0111';'Dom_c','Dom_c','0','C','12.0111';'Dy','Dy+++','0','Dy','162.5';'Dy','Dy++','0','Dy','162.5';'E','e-','0','0','';'Edta','Edta----','0','Edta','288.2134';'Er','Er+++','0','Er','167.26';'Er','Er++','0','Er','';'Ethylene','Ethylene','0','Ethylene','28.0536';'Ethylenediamine','Ethylenediamine','2','60.099','60.099';'Eu','Eu+++','0','Eu','151.965';'Eu','Eu++','0','Eu','151.965';'F','F-','0','18.9984','18.9984';'Fe','Fe++','0','55.847','55.847';'Fe','Fe+++','-2','Fe','55.847';'Formate','Formate-','0','45.018','45.018';'Four_methylpyridine','Four_methylpyridine','1','94','94';'Four_picoline','Four_picoline','1','93.128','93.128';'Fulvate','Fulvate--','0','650','650';'Ga','Ga+++','0','Ga','69.723';'Gd','Gd+++','0','Gd','157.25';'Gd','Gd++','0','Gd','157.25';'Glu','HGlu-','0','Glu','195.1459';'Glutamate','Glutamate--','1','145.115','145.115';'Glycine','Glycine-','1','74.059','74.059';'H','H+','-1','H','1.0079';'H','H3O+','-1','H','1.008';'H','H2','0','H','1.0079';'He','He','0','He','4.0026';'Hexylamine','Hexylamine','1','101.192','101.192';'Hf','Hf++++','0','Hf','178.49';'Hg','Hg(OH)2','0','200.59','200.59';'Hg','Hg','0','200.59','200.59';'Hg','Hg2++','0','Hg2','401.18';'Ho','Ho+++','0','Ho','164.9303';'Ho','Ho++','0','Ho','164.9303';'Humate','Humate--','0','2000','2000';'Hydrogen_gas','Hydrogen_gas','0','Hydrogen_gas','2.016';'I(-03)','I3-','0','I','126.9045';'I','IO-','0','I','126.9045';'I','I-','0','I','126.9045';'I','IO3-','0','I','126.9045';'I','IO4-','0','I','126.9045';'In','In+++','0','In','114.82';'Isa','HIsa-','0','HIsa','180.1548';'Isobutyrate','Isobutyrate-','1','87.043','87.043';'Isophthalate','Isophthalate--','0','164.117','164.117';'Isopropylamine','Isopropylamine','1','59.111','59.111';'Isovalerate','Isovalerate-','1','101.125','101.125';'K','K+','0','39.102','39.102';'Kr','Kr','0','Kr','83.8';'La','La+++','0','La','138.9055';'La','La++','0','La','138.9055';'Li','Li+','0','6.939','6.939';'Lu','Lu+++','0','Lu','174.967';'Methylamine','Methylamine','1','31.018','31.018';'Mg','Mg++','0','Mg','24.305';'Mn','Mn++','0','54.938','54.938';'Mn','MnO4--','0','54.938','54.938';'Mn','MnO4-','0','54.938','54.938';'Mo','MoO4--','0','Mo','95.94';'Mtg','Mtg','0','Mtg','16.032';'N','NH3','1','N','14.0067';'N','NO3-','0','14.0067','14.0067';'N','N2','0','14.0067','14.0067';'N(-03)','N3-','0','N','14.0067';'N','NO2-','0','14.0067','14.0067';'N','NH4+','0','14.0067','14.0067';'Na','Na+','0','22.9898','22.9898';'Nb','Nb(OH)6-','0','Nb','92.9064';'Nbutylamine','Nbutylamine','1','73','73';'Nd','Nd+++','0','Nd','144.24';'Nd','Nd++','0','Nd','144.24';'Ne','Ne','0','Ne','20.1797';'Ni','Ni++','0','Ni','58.69';'Nitrogen_gas','Nitrogen_gas','0','Nitrogen_gas','28.0134';'Np','Np++++','0','Np','237.048';'Np','Np+++','0','Np','237.048';'Np','NpO2+','0','Np','237.048';'Npropylamine','Npropylamine','1','59.04','59.04';'Nta','Nta---','0','Nta','188.1165';'O','H2O','0','O','15.9994';'O','O2','0','O','15.9994';'O_phthalate','O_phthalate--','0','1','1';'Ox','Ox--','0','Ox','88.0196';'Oxg','Oxg','0','Oxg','32';'P','H2(PO4)-','1','P','30.9738';'P','HPO4--','2','P','30.9738';'P','PO4---','2','30.9738','30.9738';'P','PH4+','0','P','30.9738';'Pa','Pa++++','0','Pa','231.0359';'Pa','PaO2+','0','Pa','231.0359';'Para_acetate','Para_acetate-','1','134.14','134.14';'Pb','Pb++','0','207.19','207.19';'Pb','Pb++++','0','Pb','207.19';'Pd','Pd++','0','Pd','106.42';'Phenylacetate','Phenylacetate-','0','135.142','135.142';'Phthalate','Phthalate--','1','164.117','164.117';'Pm','Pm+++','0','Pm','147';'Pm','Pm++','0','Pm','147';'Pr','Pr+++','0','Pr','140.9076';'Pr','Pr++','0','Pr','140.9076';'Propanoate','Propanoate-','1','73.032','73.032';'Propionate','Propionate-','1','73.072','73.072';'Propylamine','Propylamine','1','59.111','59.111';'Pu','Pu++++','0','Pu','244';'Pu','Pu+++','0','Pu','244';'Pu','PuO2+','0','Pu','244';'Pyrophos','Pyrophos----','0','Pyrophos','173.95';'Ra','Ra++','0','Ra','226.025';'Rb','Rb+','0','85.4699','85.4699';'Re','ReO4-','0','Re','186.207';'Rn','Rn','0','Rn','222';'Ru','RuO4--','0','Ru','101.07';'Ru','Ru++','0','Ru','101.07';'Ru','Ru+++','0','Ru','101.07';'Ru','Ru(OH)2++','0','Ru','101.07';'Ru','RuO4-','0','Ru','101.07';'Ru','RuO4','0','Ru','101.07';'S','SO4--','0','96.0616','32.066';'S','S2O3--','0','S','32.066';'S','H2S','0','32.064','32.066';'S','HS-','1','S','32.066';'S','S2O4--','0','S2O4-2','32.066';'S','SO3--','0','S','32.066';'S','S2O5--','0','S','32.066';'S','S2O8--','0','S','32.066';'S','HSO5-','0','S','32.066';'Salicylate','Salicylate--','1','136.107','136.107';'Sb','Sb(OH)6-','0','Sb','121.76';'Sb','Sb(OH)3','0','Sb','121.76';'Sb','Sb(OH)5','0','Sb','121.76';'Sc','Sc+++','0','Sc','44.9559';'Scn','Scn-','0','Scn','58.084';'Se','SeO3--','0','Se','78.96';'Se','HSe-','0','Se','78.96';'Se','HSeO3-','0','Se','78.96';'Sg','H2Sg','1','H2Sg','34.08';'Si','H4(SiO4)','0','Si','28.0855';'Si','H4SiO4','0','60.0843','28.0843';'Si','SiO2','0','SiO2','28.0855';'Sm','Sm+++','0','Sm','150.36';'Sm','Sm++','0','Sm','150.36';'Sn','Sn(OH)6--','0','Sn','118.71';'Sn','Sn++','0','Sn','118.71';'Sn','Sn(OH)2','0','Sn','118.71';'Sn','Sn++++','0','Sn','118.71';'Sr','Sr++','0','87.62','87.62';'Tartarate','Tartarate--','0','148.072','148.072';'Tartrate','Tartrate--','0','148.09','148.09';'Tb','Tb+++','0','Tb','158.9253';'Tb','Tb++','0','Tb','158.9253';'Tc','TcO(OH)2','0','Tc','98';'Tc','Tc+++','0','Tc','98';'Tc','TcO++','0','Tc','98';'Tc','TcO4---','0','Tc','98';'Tc','TcO4--','0','Tc','98';'Te','HTe-','0','Te','127.6';'Th','Th++++','0','Th','232.0381';'Thiocyanate','Thiocyanate-','0','SCN','58';'Three_methylpyridine','Three_methylpyridine','1','94','94';'Three_picoline','Three_picoline','1','93.128','93.128';'Ti','Ti(OH)4','0','Ti','47.88';'Tl','Tl(OH)3','0','Tl','204.3833';'Tl','Tl+','0','Tl','204.3833';'Tl','Tl+++','0','Tl','204.3833';'Tm','Tm+++','0','Tm','168.9342';'Tm','Tm++','0','Tm','168.9342';'Tributylphosphate','Tributylphosphate','0','265.97','265.97';'Trimethylamine','Trimethylamine','1','59.111','59.111';'Two_methylpyridine','Two_methylpyridine','1','94','94';'Two_picoline','Two_picoline','1','93.128','93.128';'U','UO2++','0','238.029','238.0289';'U','U+++','0','238.0289','238.0289';'U','U++++','0','238.0289','238.0289';'U','UO2+','0','238.029','238.0289';'V','VO++','0','V','50.9415';'V','V++','0','50.94','50.94';'V','V+++','-3','50.94','50.94';'Valerate','Valerate-','1','101.125','101.125';'W','WO4--','0','W','183.85';'Xe','Xe','0','Xe','131.29';'Y','Y+++','0','Y','88.9059';'Yb','Yb+++','0','Yb','173.04';'Yb','Yb++','0','Yb','173.04';'Zn','Zn++','0','65.37','65.37';'Zr','Zr(OH)2++','0','Zr','91.224';'Zr','Zr++++','0','Zr','91.22'};
alkalinity={ 'CO3--','1','Ca0.5(CO3)0.5','50.046';'CO3--','2','HCO3','61.0173';'HCO3-','1','Ca0.5(CO3)0.5','50.05'};
lowestpriorityelements = { 'H', 'O', 'C'}; %Alkalinity is rewritten. Other definition in PHREEQC.

%Removes the following part from specie names:
removedpartsfromspecies={ ' aq.','(aq)'};

%Text copied to file with errors
unbalancedreactions={['The following (if not empty) reactions are not balanced properly. Please check:' '\n']};

%Text copied directly into the database
LLNL_AQUEOUS_PARAMETERS={ 'LLNL_AQUEOUS_MODEL_PARAMETERS' ' -temperatures ',' 0.01 25 50 75 100 125 150 175 200 225 250 275 300','',' #debye huckel a (adh) ',' -dh_a ',' 0.4907 0.5092 0.5336 0.5639 0.5998 0.6416 0.6898 0.7454 0.8099 0.886 0.9785 1.096 1.2555','',' #debye huckel b (bdh) ',' -dh_b ',' 0.3245 0.3283 0.3325 0.3371 0.3422 0.3476 0.3533 0.3592 0.3655 0.3721 0.3792 0.3871 0.3965','','#bdot and CCO2 values: taken from LLNL, since they were impossible to convert (different equations) ',' -bdot ',' 0.0374 0.041 0.043 0.044625 0.046 0.0465 0.047 0.047 0.047 0.0405 0.034 0.017 0','',' #cco2   (coefficients for the Drummond (1981) polynomial)',' -co2_coefs','        -1.0312              0.0012806','          255.9                 0.4445','      -0.001606'};

%Species to calculate activities with CO2_llnl_gamma
CO2_llnl_gamma_names={ 'CO2_aq.','H2_aq.,','H2(aq)' ,'CO2(aq)'};

%O2-reaction to be added as link to pe
O2_reaction_string = ['2H2O =  O2 + 4H+ + 4e-  ' '\n' '   -CO2_llnl_gamma' '\n' '   log_k       -85.9951' '\n' '   -delta_H    559.543 kJ/mol # Calculated enthalpy of reaction O2' '\n' '#  Enthalpy of formation:	-2.9 kcal/mol' '\n' '        -analytic   38.0229    7.99407E-03   -2.7655e+004  -1.4506e+001  199838.45' '\n' '#  Range:  0-300' '\n'];
O2_strings_needed_to_be_replaced = {'O2(aq)'};


%Various starting definitions, necessary to run code. Don't touch.
positionold=0;
go=1;
linenumber=1;
s=1;
samenames=[0];
reactionmatrix={};
outputstringall={};
nameall={};
solution_master_species={};
solution_species={};
errormatrix_mineralsgases_doublenames={};
database_cell={};
llnlgamma= '';
speciestype={ 'beginning','elements','aqueous','mineral','gas','endone','endtwo'};

%CODE
%Actual program starts here. Following opens database file and writes outputfiles
outputfile=[pathname filename '_' preferredspecies{1} '_output' filetype];
errorfile=[pathname filename '_' preferredspecies{1} '_error' filetype];
outputhandle=fopen(outputfile,'w');
fprintf(outputhandle,'#PLEASE NOTE: THIS DATABASE HAS BEEN CONVERTED FROM THE FORMAT FOR TOUGHREACT TO A SUITABLE FORMAT FOR PHREEQC. THE FOLLOWING NOTES ARE COPIED FROM THE TOUGHREACT DATABASE');
database_raw=fopen([pathname filename filetype]);

%definitions needed for code, extracted from data.
lel=length(elementalorder);
lsmsm=length(solution_master_species_matrix);


while go %Reads databasefile, line by line. Databasefile is then transferred into a cell structuree, where each line forms a cell containing a string
    string=fgetline(database_raw);
    database_cell=[database_cell; {string}];
    positionnew=ftell(database_raw);
    if positionnew == positionold
        go=0;
    else
        positionold = positionnew;
        if length(string)>5
            if  strcmp(string(1:6),'''null'''); %THE FOLLOWING PART ONLY ACTIVATES IN BETWEEN GEOT DATA BLOCKS. THE "null" line is used in GeoT/TOUGHREACT to indicate the start of a new block. The code also prints the converted parts to the new, output database
                s=s+1;
                switch speciestype{s}
                    case 'mineral'
                        elecpresent=0;
                        [mssize1,mssize2]=size(solution_master_species);
                        for i=1:mssize1; %Add elemental weight (important)
                            for j=1:lsmsm
                                if strcmp(solution_master_species{i,3},solution_master_species_matrix{j,2})
                                    solution_master_species{i,4}=solution_master_species_matrix{j,3};
                                end
                                if strcmp(solution_master_species{i,1},solution_master_species_matrix{j,1})
                                    solution_master_species{i,6}=solution_master_species_matrix{j,5};
                                end
                            end
                            if strcmp(lower(solution_master_species{i,3}),'e-')
                                elecpresent=1;
                            end
                        end
                        if elecpresent == 0
                            solution_master_species=[solution_master_species; { 'E' '-1' 'e-' '0' '0' '0'}];
                        end
                        [mssize1,mssize2]=size(solution_master_species);
                        [solution_master_species(:,1),index]=sort(solution_master_species(:,1)); %Begin sorting
                        for i=2:mssize2
                            solution_master_species(:,i)=solution_master_species(:,i)(index);
                        end %End sorting
                        eq_ms=zeros(mssize1,1);
                        for i=2:mssize1
                            if strcmp(solution_master_species{i,1},solution_master_species{i-1,1})
                                eq_ms(i)=1;
                            end
                        end %End finding same elements > 2nd of same elements gets assigned 1 - first element in series gets assigned 2
                        for i=2:mssize1
                            if eq_ms(i)==1 & eq_ms(i-1)==0
                                eq_ms(i-1)=2;
                            end
                        end
                        indices=[];
                        nrline=0;
                        for i=1:mssize1
                            if eq_ms(i) == 2
                                indices = [indices; i i];
                                nrline=nrline+1;
                            end
                            if eq_ms(i) == 1
                                indices(nrline, 2) = indices(nrline, 2) + 1;
                            end
                        end %index: left: same elements begin - right: same elements end
                        switch ignoredoublenames
                            case 'n'
                                samenames=determine_most_important_of_double_reactions(samenames,reactionmatrix,preferredspecies,preferredspecies2);
                                errormatrix_aqueous_doublenames=outputstringall(samenames>0);
                                outputstringall=outputstringall(samenames<1);
                                nameall=nameall(samenames<1);
                                reactionmatrix=reactionmatrix(samenames<1,:);
                                samenames=samenames(samenames<1);
                            case 'y'
                                erromatrix_aqeous_doublenames={{ 'Double names of aqueous species were ignored'}};
                        end
                        defined_aqueous_species={};
                        sortingmatrix={};
                        [lreactionmatrix,wreactionmatrix]=size(reactionmatrix);
                        for i=1:lreactionmatrix;
                            defined_aqueous_species=[defined_aqueous_species;reactionmatrix{i,1}{length(reactionmatrix{i,1})}];
                            for j=1:length(reactionmatrix{i,1})
                                reactionmatrix{i,1}{j}=reactionmatrix{i,1}{j}((reactionmatrix{i,1}{j} == '''') < 1);
                            end
                            [reactionmatrixi1,reactionmatrixindextemp]=sort(reactionmatrix{i,1});
                            reactionmatrix{i,1}=reactionmatrixi1;
                            reactionmatrix{i,2}=reactionmatrix{i,2}(reactionmatrixindextemp);
                            sortingstring= '';
                            for j=1:length(reactionmatrix{i,1})
                                sortingstring =[sortingstring reactionmatrix{i,1}{j}];
                            end
                            sortingmatrix=[sortingmatrix;sortingstring];
                        end
                        [sortingmatrix,sortingmatrixindex]=sort(sortingmatrix);
                        revertsorting=1:length(sortingmatrix);
                        revertsorting=revertsorting(sortingmatrixindex);
                        
                        reactionmatrix=reactionmatrix(sortingmatrixindex,:);
                        defined_aqueous_species=defined_aqueous_species(sortingmatrixindex);
                        outputstringall=outputstringall(sortingmatrixindex);
                        
                        sortingmatrix=removedoublesfromreactionmatrix(reactionmatrix);
                        
                        sortingmatrix2=1:length(sortingmatrix);
                        sortingmatrix3=sortingmatrix2(sortingmatrix<1);
                        sortingmatrix2=sortingmatrix2(sortingmatrix>0);
                        errormatrix_aqueous=outputstringall(sortingmatrix2);
                        defined_aqueous_species=defined_aqueous_species(sortingmatrix3);
                        outputstringall=outputstringall(sortingmatrix3); %FROM CASE 'MINERAL' TO HERE: REMOVE DOUBLE REACTIONS FROM AQEOUS SPECIES
                        revertsorting=revertsorting(sortingmatrix3);
                        [temprevsort,revertsorting]=sort(revertsorting);
                        outputstringall=outputstringall(revertsorting);
                        defined_aqueous_species=unique(defined_aqueous_species);
                        [indices2, solution_master_species]=sortmasterspeciesbychargeandOH(indices,solution_master_species,defined_aqueous_species,nrline,mssize2,mssize1);
                        
                        %PRINT THE MASTER SPECIES                        
                        
                        outputstring_sms={};
                        for i=1:mssize1;
                            if indices2(i) == 1
                                outputstring_sms= [outputstring_sms; solution_master_species{i,1} ' ' solution_master_species{i,3} ' ' solution_master_species{i,4} ' ' solution_master_species{i,5} ' ' solution_master_species{i,6}];
                                outputstring_sms= [outputstring_sms; solution_master_species{i,1} '(' solution_master_species{i,2} ') ' solution_master_species{i,3} ' ' solution_master_species{i,4} ' ' solution_master_species{i,5} ' ' solution_master_species{i,6}];
                            elseif eq_ms(i) == 0
                                outputstring_sms= [outputstring_sms; solution_master_species{i,1} ' ' solution_master_species{i,3} ' ' solution_master_species{i,4} ' ' solution_master_species{i,5} ' ' solution_master_species{i,6}];
                                if strcmp(solution_master_species{i,1},'H')
                                    outputstring_sms= [outputstring_sms; solution_master_species{i,1} '(' solution_master_species{i,2} ') ' solution_master_species{i,3} ' ' solution_master_species{i,4} ' ' solution_master_species{i,5} ' ' solution_master_species{i,6}];
                                end
                            else
                                outputstring_sms= [outputstring_sms; solution_master_species{i,1} '(' solution_master_species{i,2} ') ' solution_master_species{i,3} ' ' solution_master_species{i,4} ' ' solution_master_species{i,5} ' ' solution_master_species{i,6}];
                            end
                            outputstring_sms=sort(outputstring_sms);
                        end
                        for i=1:length(outputstring_sms);
                            fprintf(outputhandle,[outputstring_sms{i} '\n']);
                        end
                        
                        %PRINT THE AQUEOUS SPECIES                           
                        fprintf(outputhandle,'SOLUTION_SPECIES \n');
                        if elecpresent == 0
                            fprintf(outputhandle, ['e- = e-' '\n']);
                            fprintf(outputhandle, ['log_k 0' '\n']);
                        end
                        for i = 1:length(solution_species)
                            fprintf(outputhandle, [solution_species{i} '\n']);
                        end
                        fprintf(outputhandle, O2_reaction_string);
                        for i=1:length(outputstringall);
                            for j=1:length(outputstringall{i}) 
                                fprintf(outputhandle,[outputstringall{i}{j} '\n']);
                            end
                        end
                        nameall={};
                        reactionmatrix={};
                        samenames=[0];
                        outputstringall={};
                        fprintf(outputhandle,'PHASES \n');
                    case 'endone'
                            switch ignoredoublenames
                                case 'n'
                                    samenames=determine_most_important_of_double_reactions(samenames,reactionmatrix,preferredspecies,preferredspecies2);
                                    errormatrix_mineralsgases_doublenames=outputstringall(samenames>0);
                                    outputstringall=outputstringall(samenames<1);
                                case 'y'
                                    erromatrix_mineralsgases_doublenames={{ 'Double names of minerals and gases were ignored'}};
                                end
                                for i=1:length(outputstringall);
                                    for j=1:length(outputstringall{i}) 
                                        fprintf(outputhandle,[outputstringall{i}{j} '\n']);
                                    end
                                end
                                
                end
            elseif   length(string)>44
                if strcmp(string(1:45),'!end-of-header     Do not remove this record!');
                    s=s+1;
                    for i=1:length(LLNL_AQUEOUS_PARAMETERS)
                        fprintf(outputhandle,[llnlgamma LLNL_AQUEOUS_PARAMETERS{i} '\n']);
                    end
                    fprintf(outputhandle,'\nSOLUTION_MASTER_SPECIES \n');                    
                end
            end
        end
        %MAIN PART OF THE CODE. REWRITES AQUEOUS SPECIES, MINERALS AND GASES       
        if length(string)>8
            if string(1:9) == '# source:' %First, it identifies the last line of the block. First three lines are then extracted and rewritten
                string1=database_cell{linenumber-3};
                if string1(1)== '#'
                else
                    string2=database_cell{linenumber-2};
                    string3=database_cell{linenumber-1};
                    apostrophes=0;
                    i=0;
                    while apostrophes<2
                        if string1(i+1)== ''''
                            apostrophes=apostrophes+1;
                        end
                        i=i+1;
                    end
                    name=string1(1:i);
                    for i=1:length(name)
                        if name(i)== ' '
                            name(i) = '_';
                        end
                    end
                    string1=string1(i+1:length(string1));
                    string1=removeparts(string1,removedpartsfromspecies);
                    string2=string2(i+1:length(string2));
                    string3=string3(i+1:length(string3));
                    string1=separatestringby(string1,' ');
                    string2=separatestringby(string2,' ');
                    string3=separatestringby(string3,' ');
                    logkvalue=str2num(string3{1})*log(10);
                    stringlogkvalue=num2strmoreprecise(logkvalue, length(string3{1}));
                    adjustmentcurve=[string3(2:4) {stringlogkvalue} string3(5)];
                    adjustmentcurve2= '';
                    for i=1:length(adjustmentcurve);
                        if strcmp(speciestype{s},'aqueous')
                            adjustmentcurve2=[adjustmentcurve2 ' ' num2strmoreprecise(-str2num(adjustmentcurve{i}),length(adjustmentcurve{i}))];
                        else
                            adjustmentcurve2=[adjustmentcurve2 ' ' adjustmentcurve{i}];
                        end
                    end
                    adjustmentcurve2=deblank(adjustmentcurve2);
                    weight=string1{1};
                    switch speciestype{s} 
                        case 'aqueous'
                            startnr=5;
                            a0=str2num(string1{2});
                            charge=str2num(string1{3});
                            if charge>0
                                a0=2*(a0+1.81*abs(charge))/(abs(charge)+1);
                            elseif charge<0
                                a0=2*(a0+1.91*abs(charge))/(abs(charge)+1);
                            end                
                            a0=num2str(a0);
                        case 'mineral'
                            startnr=4;
                            vm=string1{2};
                        case 'gas'
                            startnr=4;
                            dmdiam=string1{2};
                    end
                    lefthandside= '';
                    righthandside= '';
                    elementstotal={};
                    numberstotal=[];
                    name2=name(2:length(name)-1);
                    %the following loop determines the elements/type of species from the first string in the TOUGHREACT description of minerals, gases and aqeuous species, and adds any charges
                    allspeciesfor1reaction={};
                    allnumbersfor1reaction=[];
                    for i=startnr:2:length(string1);
                        speciesnr=string1{i};
                        speciesname=string1{i+1};
                        speciesname=replaceplusminusequalsignsforspecies(speciesname);
                        if str2num(speciesnr)<0
                            righthandside = [righthandside ' + ' speciesnr(2:length(speciesnr)) ' ' speciesname(2:length(speciesname)-1)];
                        elseif str2num(speciesnr)>0
                            lefthandside = [lefthandside ' + ' speciesnr ' ' speciesname(2:length(speciesname)-1)];
                        end
                        [elements,numbers]=extractspeciesfromname(speciesname);
                        numbers2=zeros(size(numbers));
                        lengths=zeros(size(numbers));
                        for i=1:length(numbers);
                            numbers2(i)=str2num(numbers{i}).*str2num(speciesnr);
                        end
                        numberstotal=[numberstotal;numbers2];
                        elementstotal=[elementstotal;elements];
                        allspeciesfor1reaction=[allspeciesfor1reaction {speciesname}];
                        allnumbersfor1reaction=[allnumbersfor1reaction str2num(speciesnr)];
                    end
                    
                    [elementstotalnew, numberstotalnew]=correctOH(elementstotal, numberstotal);
                    
                    specieswarning=0;
                    for i=1:length(numberstotalnew)
                        if numberstotalnew(i)<0;
                            specieswarning=1;
                        end
                    end                    
                    if specieswarning
                        unbalancedreactions=[unbalancedreactions; database_cell{linenumber-3};database_cell{linenumber-2};database_cell{linenumber-1};[database_cell{linenumber} '\n']];
                    else                    
                        elementstotalnew=elementstotalnew(numberstotalnew>0);
                        numberstotalnew=numberstotalnew(numberstotalnew>0);
                        %Sorting algorithm to sort according to chemistry
                        elementsorting=zeros(size(elementstotalnew));
                        for i=1:length(elementstotalnew);
                            for j=1:lel
                                if strcmp(elementalorder{j},elementstotalnew{i})
                                    elementsorting(i)=j;
                                end
                            end
                        end
                        [elementsorting,elementsortingindex]=sort(elementsorting);
                        elementstotalnew=elementstotalnew(elementsortingindex);
                        numberstotalnew=numberstotalnew(elementsortingindex);
                        %End of sorting algorithm
                        
                        mainspecies= '';
                        for i=1:length(elementstotalnew);
                            if numberstotalnew(i)==1
                                mainspecies=[mainspecies elementstotalnew{i}];
                            else
                                mainspecies=[mainspecies elementstotalnew{i} num2str(numberstotalnew(i))];
                            end
                        end
                        if charge>0
                            for i=1:charge
                                mainspecies = [mainspecies '+'];
                            end
                        elseif charge<0
                            charge = - charge;
                            for i=1:charge
                                mainspecies = [mainspecies '-'];
                            end
                        end
                        allspeciesfor1reaction=[allspeciesfor1reaction mainspecies];
                        allnumbersfor1reaction=[allnumbersfor1reaction -1];                    
                        reactionmatrix=[reactionmatrix;{allspeciesfor1reaction} allnumbersfor1reaction];
                        lefthandside=deblank(lefthandside);
                        righthandside=deblank(righthandside);
                        if lefthandside(length(lefthandside)) == ''''
                            lefthandside=lefthandside(1:length(lefthandside)-1);
                        end
                        if isempty(righthandside)
                        else
                            if righthandside(length(righthandside)) == ''''
                                righthandside=righthandside(1:length(righthandside)-1);
                            end
                        end
                        if strcmp(speciestype{s},'gas')
                            name2=deblank(name2);
                            lengthname2=length(name2);
                            if lengthname2>4
                                if strcmp('_gas',name2(lengthname2-3:lengthname2)),
                                    name2=name2(1:lengthname2-4);
                                    lengthname2=lengthname2-4;
                                end
                            end
                            if lengthname2>3
                                if strcmp('(g)',name2(lengthname2-2:lengthname2))
                                    name2=name2(1:lengthname2-3);
                                end                                
                            end
                            name2=[name2 '(g)'];
                        end
                        
                        lefthandside = lefthandside(4:length(lefthandside));    
                        switch speciestype{s}
                            case 'aqueous'
                                reaction = [lefthandside ' = ' mainspecies righthandside];
                                if  str2num(a0) == 0
                                    if sum(strcmp(name2,CO2_llnl_gamma_names))>0
                                        outputstring={['#' name2] reaction [llnlgamma '     -CO2_llnl_gamma'] ['    log_k ' num2str(-str2num(string2{2}))] ['    -analytic ' adjustmentcurve2] ['#    weight= ' weight] string};
                                    else
                                        outputstring={['#' name2] reaction [llnlgamma '     -llnl_gamma 3'] ['    log_k ' num2str(-str2num(string2{2}))] ['    -analytic ' adjustmentcurve2] ['#    weight= ' weight] string};
                                    end
                                    if sum(strcmp(name2,O2_strings_needed_to_be_replaced))>0
                                        for opst=1:length(outputstring)
                                            outputstring{opst}=['#' outputstring{opst}];
                                        end
                                    end
                                    
                                else
                                    outputstring={['#' name2] reaction [llnlgamma '     -llnl_gamma ' a0] ['    log_k ' num2str(-str2num(string2{2}))] ['    -analytic ' adjustmentcurve2] ['#    weight= ' weight] string};
                                end
                           case 'mineral'
                               reaction = [mainspecies righthandside ' = ' lefthandside];
                               outputstring={name2 reaction ['    log_k ' string2{2}] ['    -Vm ' vm] ['    -analytic ' adjustmentcurve2]  ['#    weight= ' weight] string};
                           case 'gas'
                               reaction = [mainspecies righthandside ' = ' lefthandside];
                               outputstring={name2 reaction ['    log_k ' string2{2}] ['    -analytic ' adjustmentcurve2] ['#    molecular diameter (m) = ' dmdiam] ['#    weight= ' weight] string};
                        end
                        outputstringall=[outputstringall;{outputstring}];
                        nameall=[nameall;name2];
                        if length(nameall)>1;
                            if strcmp(nameall{length(nameall)},nameall{length(nameall)-1})
                                samenames=[samenames;1];
                            else
                                samenames=[samenames;0];
                            end
                        end
                    end 
                    
                end
            elseif strcmp(speciestype{s},'elements') & string(1)== '''' %Following loop extracts primary species
                if strcmp(string(1:21),'''Temperature points:''')
                    string=string(22:length(string));
                    temperaturepoints=separatestringby(string,' ');
                    temperaturepoints=temperaturepoints(2:length(temperaturepoints));
                    temp2=zeros(size(temperaturepoints));
                    for i=1:length(temperaturepoints)
                        temp2(i)=str2num(temperaturepoints{i});
                    end
                    temperaturepoints=temp2;
                else
                    apostrophes=0;
                    i=0;
                    while apostrophes<2
                        if string(i+1)== ''''
                            apostrophes=apostrophes+1;
                        end
                        i=i+1;
                    end
                    name=string(1:i);
                    string=string(i+1:length(string));
                    name=removeparts(name,removedpartsfromspecies);
                    string=separatestringby(string,' ');
                    a0=str2num(string{1});
                    charge=str2num(string{2});
                    if charge>0
                        a0=2*(a0+1.81*abs(charge))/(abs(charge)+1);
                    elseif charge<0
                        a0=2*(a0+1.91*abs(charge))/(abs(charge)+1);
                    end                
                    a0=num2str(a0);
                    weight=deblank(string{3});
                    name=name(2:length(name)-1);
                    nameremover= (name == '+') + (name == '-') + (name == '=');
                    for i=2:length(name)
                        if isdigit(name(i))
                            if nameremover(i-1) == 1;
                                nameremover(i) = 1;
                            end
                        end
                    end
                    name=name(nameremover<1);
                    possibleelements=name((isdigit(name) | name == '(' | name == ')')<1);
                    capitalpositions=1:length(possibleelements);
                    capitalpositions=capitalpositions((possibleelements == lower(possibleelements))<1);
                    possibleelements2=cell(size(capitalpositions));
                    for i=1:length(capitalpositions);
                        if i<length(capitalpositions)
                            possibleelements2{i}=possibleelements(capitalpositions(i):capitalpositions(i+1)-1);
                        else
                            possibleelements2{i}=possibleelements(capitalpositions(i):length(possibleelements));
                        end
                    end
                    elementremover=1;
                    while length(possibleelements2)>1;
                        possibleelementsexcluder=ones(size(possibleelements2));
                        for i=1:length(possibleelements2);
                            if strcmp(possibleelements2{i}, lowestpriorityelements{elementremover});
                                possibleelementsexcluder(i)=0;
                            end
                        end
                        possibleelements2=possibleelements2(possibleelementsexcluder>0);
                        elementremover=elementremover+1;
                    end
                    element=possibleelements2{1};
                    if charge>0
                        for i=1:charge
                            name = [name '+'];
                        end
                    elseif charge<0
                        for i=1:-charge
                            name = [name '-'];
                        end
                    elseif charge == 0
                        a0 = '3';
                    end
                    solution_master_species=[solution_master_species;{element num2str(charge) name '0' weight weight}];
                    solution_species=[solution_species; {[name ' = ' name]};{[llnlgamma '    -llnl_gamma ' a0]};{ '    log_k 0'}];
                end    
            end
        end
        if s == 1 | s>5
            if length(deblank(string))<1
                fprintf(outputhandle,string);
            else
                string = strrep(string,'%','%%');
                fprintf(outputhandle, ['#' string]);
            end    
        end
    end
    
    
    linenumber=linenumber+1;
end
fclose(outputhandle);
%END OF MAIN CODE. Database is read, and transferred to a file format for PHREEQC


%The following code prints the errorfile, using the various parts that could not be converted, were double, etc.
errorhandle=fopen(errorfile,'w');
for i=1:length(unbalancedreactions)
    fprintf(errorhandle,unbalancedreactions{i});
end
fprintf(errorhandle,['\n' 'The following reactions for aqueous species are double and may need to be defined in another way. Please check:' '\n']);
for i=1:length(errormatrix_aqueous);
    for j=1:length(errormatrix_aqueous{i}) 
        fprintf(errorhandle,[errormatrix_aqueous{i}{j} '\n']);
    end
end
fprintf(errorhandle,['\n' 'The following reactions for aqueous species have been defined in multiple ways, meaning one (or more) have been deleted:' '\n']);
for i=1:length(errormatrix_aqueous_doublenames);
    for j=1:length(errormatrix_aqueous_doublenames{i}) 
        fprintf(errorhandle,[errormatrix_aqueous_doublenames{i}{j} '\n']);
    end
end
fprintf(errorhandle,['\n' 'The following reactions for minerals and gases have been defined in multiple ways, meaning one (or more) have been deleted:' '\n']);
for i=1:length(errormatrix_mineralsgases_doublenames);
    for j=1:length(errormatrix_mineralsgases_doublenames{i}) 
        fprintf(errorhandle,[errormatrix_mineralsgases_doublenames{i}{j} '\n']);
    end
end
fclose(errorhandle);
