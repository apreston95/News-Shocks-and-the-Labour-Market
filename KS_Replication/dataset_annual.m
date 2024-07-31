%dataset_annual.m
%
%reads in data and transforms to desired dataset
%
%Last modified: Andre Kurmann, March 2016

function [Y, X, vars, lev]=dataset_annual(esty1,esty2,nlags)

%Step 1: load data
%-----------------
dataseries = xlsread('KurmannSims_Part2_annual.xlsx',1,'A1:R66');  %see xls sheet for details

if ismac
    dataseries = dataseries(2:end,2:end);    %in a mac environment, the range option of xlsread.m does not work and the entire sheet is loaded instead
                                             %for the present case, the "entire" sheet begins at row 2, column 1         
end

varnames =   ['  Gross domestic product   '      %1) real gdp
              '  Consumption expenditures '      %2) real c
              '  Gross private investment '      %3) real i
              '      Total hours          '      %4) total hours
              '   Fernald TFP (2007)      '      %5) Fernald's utilization-adjusted TFP series (2007 vintage)
              '   Fernald TFP (2016)      '      %6) Fernald's utilization-adjusted TFP series (2016 vintage)
              '      Inflation            '      %7)
              '     Spread (5yr-ffr)      '      %8)
              '     real S&P500 index     '      %9)
              '  Alexopoulos TECH         '      %10) Bowker new book titles in the technology subject 
              '  Alexopoulos TECH2        '      %11) Alexopoulos? broader indicator of new technology book titles from Library of Congress records 
              '  Alexopoulos COMP         '      %12) Alexopoulos? indicator of new titles on computer software and hardware 
              '  Alexopoulos TEL          '      %14) Alexopoulos? indicator of new titles on telecommunications 
              '  Alexopoulos HIS          '      %15)
              '  Alexopoulos SCI          '      %16)
              '  NSF SIRD R&D             '];    %17)
             
%rescaling some of the variables   
dataseries = 100*dataseries;
dataseries(:,[7 8]) = dataseries(:,[7 8])/100;

%Step 2: setup VAR and compute lags
%----------------------------------
%select variables in VAR
var_select=[5 11 2 7];
lev = [1 1 1 0];

%select sample period
daty1=1947;       % First Year of Data Set
n1 = (esty1-daty1)+1;
n2 = (esty2-daty1)+1;
     
data = dataseries(n1:n2,var_select);
vars = varnames(var_select,:);

%checking for out-of-sample values
if sum(isnan(data))>0 | sum(sum(data<-998))>0;
    data
    disp('dataseries out of range')
    return
end

[T,nvars]=size(data);
%compute nlags lags
for p=1:nlags;
    X(:,1+(p-1)*nvars:p*nvars)=data((nlags+1-p):(T-p),:);
    %gives (T-nlags) x nvars*nlags matrix of lags for all variables
    %first lag of all variables first, then second lag of all variables
    %and so on...                                               
end;

%rescaling variables since we loose nlags observations through the lagging 
Y=data((nlags+1):T,:);



