function q=QualityOfScopes(doplot)
% Performs Optinal Segmentation on the scans in every csv file residing in the current directory. Although multiple scopes may be in each file the
% triplet scans (ULLR, LLUR, Hor) must be in adjacent titled columns 2 to 3*n. The first column must be titled "ID" and hold values 1:Longest scan.
% Input parameters: 
% doplot: if empty all scans in every file are segmented and and plotted, otherwise only the numbered scans listed are used.
% Output:
% For a single scope results are displayed on the screnen.  For multiple scopes a csv file is output
% q{1} holds cell array output of example2
% q{2} holds scope names
% Bugs:
% Presence of outliers may lead to failure of search for 3 segment solution, which is indicated by last column of output, MaxNumSegs, being
% a number other than 3.  Solution: Re-run with outlier deleted.

if(nargin==0) doplot=1:1000;end;% allow a selection of plots to be run, if doplot absent all are.
[y z]=system('ls -1 *.csv');z=regexp(z,'\n','split');pltnum=0;
for i=1:1:length(z)-1;%Outer loop reads whole csv files. Last element of z is empty
    vn=importfile(z{i});if(i==1) bw=[];end
    for j=1:1:size(vn)-1
        ID=evalin('base','ID');pltnum=pltnum+1;
                    if(sum(doplot==pltnum)==1) r=Segment(ID,vn{j+1},pltnum,nargin);bw=strvcat(bw,vn{j+1});
                            if(exist('result')==0) result={r};else result=horzcat(result,{r});end;
                    end;
    end;
end;
q=QualityScore({bw;result});[x,idx]=sort(q{1}(:,2));q{1}=q{1}(idx,:);q{2}=q{2}(idx,:);qsz=size(idx);
% Output results, write file if more than one scope scored
header=strvcat('ScopeName','ScopeNumber','QualityScore ','Score','RMSI','LinMid','LinRight','FirstDif','ScndDif','MaxNumSegs');
if(qsz(1)>1)
            [x,f]=system('date');filename=strcat('ReportSortedScopes_',regexprep(f(f~=' '),':',''),'.csv');
            fid=fopen(filename,'w');fprintf(fid,'%s',header');fprintf(fid,'%s\n','');
            for hi=1:1:qsz(1) 
                fprintf(fid,'%s %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f\n', q{2}(hi,:),q{1}(hi,:));end;
            fclose(fid);
            disp(strcat('Score output written to file: ',filename));
else
    disp(sprintf('\n%s\t %s',header(1,:),q{2}(1,:)));
    for h=2:1:9
        disp(sprintf('%s\t %10.4f', header(h+1,:),q{1}(:,h)));end;
end;

function y = Segment(x,varname,pltnum,upnarg)
% Perform optinal segmentation with PELTMeanNorm on input microscope scan
% held in varname. The series is transformed by:trimming, scaling
% Input:
% x: holds first column of csv file,ID=1:Longest scan
% varname: variable holding a scopes's scan name
% pltnum: if pltnum=1 start new plot window otherwise add to growing 6x6 plot matrix in current window (for large numbers of scopes multiple windows)
% upnarg: if upnarg~=0 (when Scopes argument is non-empty) then individual plots are displayed unwindowed 
% Output: cell array
% segnum: the segment number [1,2,3]
% segindexes: segment start:segment end for each segment
% segvalues: intensity values corresponding segindex for each segment
% startend:segment start, segment end, segment intercept, segment slope for each segment
myseries=evalin('base',varname);
% Transformations: 
x=x(6:end-5,1)';myseries=myseries(6:end-5,1)';% trim for flyback
myseries=-255*(1-myseries/max(myseries));% y-axis Intensity is scaled to  -(0:255)
myzeros=myseries<-249;myseries=myseries(myzeros~=1);x=x(myzeros~=1);%remove 0, these values were added to make csv files square
x=715*x/max(x);rnote=strcat('Current Variable= ',varname);% scaled abscissa

complw='y=PELTMeanNorm(myseries, fac*log(10));z=y.ChngPts;z(1)=1;z(length(z))=length(x);w=[z(2:length(z)-1),y.N+1]-z(1:length(z)-1);lw=1+length(w);';
% NB. z(1)=1;z(length(z))=length(x); assures for rare cases the end segments encompass all data
fac=5000;maxit=200;it=0;eval(complw);% search for fac yielding 3 segments. Such a value doesn't always exist....
if(lw<4) while(lw<4) fac=fac/2;cc=fac;eval(complw);end
        up=2*fac;
else while(lw>4) fac=fac+100;cc=fac;eval(complw);end
        up=fac;fac=fac-100;
end;up=111/100;
while(~((maxit<it) || (lw==4))) fac=up*fac;eval(complw);it=it+1;cc=0;
     while((lw>4)&&(it<maxit)) it=it+1;fac=fac*up;cc=cc+1;eval(complw);end
     fac=fac/up;up=1-(1-up)/(cc+it);
     disp(strcat('Iteration :',num2str(it),'. Value of Fac :',num2str(fac),' length(w) :',num2str(length(w))));
end
  c.segnum='Seg1';c.segindexes=1:w(1);c.segindep=x([c.segindexes]);c.segvalues=myseries(1:w(1));
  m=polyfit(c.segindep,c.segvalues,1);model=polyval(m,c.segindep);c.startend=[1,(w(1)),m(1,[2 1])];
  for i=2:1:length(w)
      c(i).segnum=strcat('Seg',num2str(i));
      c(i).segindexes=z(i)+(1:w(i))-1;
      c(i).segindep=x([c(i).segindexes]);
      c(i).segvalues=myseries(c(i).segindexes);
      m=polyfit(c(i).segindep,c(i).segvalues,1);
      c(i).startend=horzcat(z(i),(z(i)+w(i)-1),m(1,[2 1]));
      model=horzcat(model,polyval(m,c(i).segindep));% needed for plotting
  end
  o=[];for i=1:1:length(w) o=vertcat(o,horzcat(c(i).startend));end
rnote=strvcat(rnote,'Start, End, B0, B1 of each segment:');
rnote=strvcat(rnote,sprintf(' %5i %5i %10.2f %10.2f \n',o'));disp(rnote);
% cax=newplot;set(cax,'NextPlot','new','ZGrid','on');% not needed use
% figure() instead
if(upnarg==0) 
   if((pltnum==1)||(mod(pltnum,37)==0)) fh=figure('NumberTitle','off');disp(strcat('Plot figure handle=',num2str(fh)));end;
                                        subplot(6,6,1+mod((pltnum-1),36));
end;
if(upnarg~=0) fh=figure('NumberTitle','off');end;
plot(x,myseries,'or',x,model,'-b');title(varname);axis([0 725 -255 0]);
if(upnarg~=0) fn=strcat('mlplots/',varname,'.png');fn2=strcat('mlplots/',varname,'.txt');dlmwrite(fn2,rnote,' ');print('-dpng',fn);end;
disp(strcat('Returning 6 length ', num2str(length(w)),' nested vectors: segnum, segindexes, segindep, segvalues, startend'));
y=c;

function ans = PELTMeanNorm (d, p)
%   call with (time series data, penalty)
a=PELT_mean_norm (d, p);                                                                              
if(0~=a.Error(1)) 
     sprintf('Error! Type=  %i', error);
     a.Error=strcat('Error! Type = ', num2str(a.Error));
else 
     a.ChngPts=sort(a.ChngPts(a.ChngPts>0));
     sprintf('%s','Returning Structure; 1) +\y*2; 2) +\y; 3) n; 4) penalty; 5) cut points; 5) error');
end                                                                  
ans=a; 

function y=mll_mean(x,x2,n) 
y=x2-(x*x)/n;

function y=min_which(array)
y=[min(array),find(array==min(array))];

function y=max_which(array)
y=[max(array),find(array==max(array))];

function a=PELT_mean_norm(ts,pen)
y2=[0,cumsum(ts.*ts)];y=[0,cumsum(ts)];error=0;n=length(ts);cptsout=repmat(0,1,n);                                                                                             
n2=(n+1)*2;if(0<n2) 
    lastchangecpts=repmat(0,1,n2);% for storage of (cp(tmin),tmin)
else
    error=[1,1]; 
end
np1=n+1;if(0<np1) 
    lastchangelike=repmat(0.0,1,np1);% this is F 
else
    error=[1,2];
end
if(0<n) 
    checklist=repmat(0,1,n);% this is R
else error=[1,3];
end
if(0<np1) 
    tmplike=repmat(0.0,1,np1);% for storage of likelihoods over which the min is put into F
else
    error=[1,4];
end
if(2==sum(size(error)))
    lastchangelike(1)=-pen;                   % Set F(0) to -B
  for i=1:1:3 
      lastchangelike(i+1)=mll_mean(y(i+1),y2(i+1),i);% computes relevant part of likelihood (+/y*2)-((+/y)*2)?n for y[2],..y[4] 
  end
  lastchangecpts([(1:4),n+(2:5)])=[0 0 0 0,1:4];                 % initialize cp to 4/0, tmin to 1:4
  nchecklist=2;checklist([1 2])=[1 3];                           % initialize R to 0,2 (insist on at least 2 obs per segment)

  tstar=5;minout=-99.0;whichout=-99.0;nchecktmp=-99.0;
               while(tstar<n+2)                                                                                                                                          
                    i=1;while(i<nchecklist+1)
                        tmplike(i)=lastchangelike(checklist(i)) + pen + mll_mean((y(tstar)-y(checklist(i))),(y2(tstar)-y2(checklist(i))),tstar-checklist(i));
                           i=i+1; % for t* in R[t*] compute F(t) + C(y[t+1],..,y[t*]) + B, where t*=tstar which is step 1
                    end
                                     
                    z=min_which(tmplike(1:nchecklist));minout=z(1);whichout=z(2);% min and which element is Step2
                    lastchangelike(tstar)=minout; lastchangecpts([0,n]+tstar)=[checklist(whichout),tstar-1];
                    % disp(sprintf('tstar = %d whichout = %d, checklist(whichout) = %d', tstar, whichout, checklist(whichout)));
                    % disp(sprintf('lastchangecpts ='));disp(sprintf(' %i ',lastchangecpts));
                    % Update checklist for next iteration, first element is next tau
                    i=1;nchecktmp=1;
                    while(i<nchecklist) 
                         if(tmplike(i)<=(lastchangelike(tstar)+pen))
                           checklist(nchecktmp)=checklist(i);
                           nchecktmp=nchecktmp+1;
                         end
                         i=i+1;
                    end 
              checklist(nchecktmp+1)=tstar-1;% enforce at least 2 obs per seg %%checklist(nchecktmp)=tstar-1;
              nchecktmp=nchecktmp+1;nchecklist=nchecktmp;tstar=tstar+1;
              end

    ncpts=1;last=n;% backtrack to backtrack to get final set of changepoints
         while(last~=0)
              cptsout(ncpts)=lastchangecpts(np1 +last);last=lastchangecpts(last);ncpts=ncpts+1;
              % disp(sprintf('cptsout(%i)= %i, last= %i', (ncpts-1), cptsout(ncpts-1), last));
         end
else 
    y2=-99;y=-99;n=-99;pen=-99;cptsout=-99;
end
a=struct('CumSS',y2,'CumSum',y,'N',n,'Penalty',pen,'ChngPts',cptsout,'Error',error);

function s=QualityScore(Y)
% Score = 100xmin(1,1-min(1,(([-13.05841809 0.7083353545 4.638056856 5.174318301 -57.83098138 44.9805868 -0.001361034229]+.x(Mean;RMSI;LM;LR;CL;CQ;RMSI*2))-24)/71))
% where RMSI=tot error from max intensity
%        LM = total adjusted abs middle slope from 3 segments totaled over 3 traces
%        LR = total adjusted abs right slope  from 3 segments totaled over 3 traces
%        CL = total abs I(t)-I(t-1), 
%        CQ = total abs I(t)-2I(t-1) + I(t-2)
% and total means sum over three traces, sum over t for CL and CQ.
RMSI=[];curve=[];CL=[];CQ=[];N=[];big=[];scopenames=[];szsave=[1 3];
szd=size(Y{2})/3;qn=Y{1};
for scopes=1:1:szd(2)
   for d=1:1:3
      q=cell2mat(Y{2}(d+3*(scopes-1)));sz={q.segnum};sz=size(sz);
      if(sz(2)>szsave) 
          szsave=sz;end;
      for i=1:1:sz(2)
        x=q(i).segindep;dat=q(i).segvalues;b=q(i).startend;
        curve=horzcat(curve,abs(b(4))*(1+b(2)-b(1)));
        RMSI=horzcat(RMSI,(sum(dat.^2)));
        cl=0;N=horzcat(N,length(dat));
        for j=2:1:length(dat)
            cl=cl+abs(dat(j)-dat(j-1));
        end
        CL=horzcat(CL,cl/length(dat));
        cq=0;
        for j=3:1:length(dat)
            cq=cq+abs(dat(j)+dat(j-2)-2*dat(j-1));
        end
        CQ=horzcat(CQ,cq/length(dat));
      end
   end
      RMSI=(sum(RMSI)/sum(N))^0.5;FirstDif=sum(CL);ScndDif=sum(CQ);LinMid=sum(curve(1,[2 5 8]));LinRight=sum(curve(1,1+[2 5 8]));% when Segment produces > 3 segments only these are counted.
% Scoring function:
% 100xmin(1,1-min(1,(([-13.05841809 0.7083353545 4.638056856 5.174318301 -57.83098138 44.9805868 -0.001361034229]+.x(Mean;RMSI;LM;LR;CL;CQ;RMSI*2))-24)/71))
% where RMSI=tot error from max intensity, LM=total adjusted abs middle slope, LR =total adjusted abs right slope, CL = total abs I(t)-I(t-1), CQ =  total abs I(t)-2I(t-1) + I(t-2)
% and total means sum over three traces and also t for CL and CQ.
      score=sum([-13.05841809 0.7083353545 4.638056856 5.174318301 -57.83098138 44.9805868 -0.001361034229].*[1,RMSI,log([LinMid,LinRight,FirstDif,ScndDif]),RMSI^2]);
      qualscore=100*min(1,(1-min(1,(score-24)/71)));
      big=vertcat(big,[scopes,qualscore,score,RMSI,LinMid,LinRight,FirstDif,ScndDif,szsave(2)]);% szsave holds max number of segments found over all traces.
      scopenames=vertcat(scopenames,regexprep(qn(d+3*(scopes-1),1:end),'ULLR|LLUR|HOR|Hor',''));%Assumes scopenames ends with these tags
RMSI=[];curve=[];CL=[];CQ=[];N=[];szsave=[1 3];
end
disp('Returning: Scope number,QualityScore,score,RMSI,LinMid,LinRight,FirstDif,ScndDif,NumSegments');
s={big;scopenames};

function vn=importfile(fileToRead1)
%IMPORTFILE(FILETOREAD1)
%  Imports data from the specified file
%  FILETOREAD1:  file to read

%  Auto-generated by MATLAB on 13-Jun-2013 16:55:56

% Import the file
newData1 = importdata(fileToRead1);

% Break the data up into a new structure with one field per column.
colheaders = genvarname(newData1.colheaders);
for i = 1:length(colheaders)
    dataByColumn1.(colheaders{i}) = newData1.data(:, i);
end

% Create new variables in the base workspace from those fields.
vars = fieldnames(dataByColumn1);
for i = 1:length(vars)
    assignin('base', vars{i}, dataByColumn1.(vars{i}));
end
vn=vars;



