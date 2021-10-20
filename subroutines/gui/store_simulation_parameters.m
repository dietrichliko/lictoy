function varargout=store_simulation_parameters(src,evt,parameterfile,usestrings);
h=get(src,'parent');
data=guidata(h);
Text=data.Text;
Values=data.Values;
hedit=data.hedit;
for i=1:length(Text)
    text=Text{i};
    if usestrings(i)
        %values=Values{i};
       text=[text,':'];
       str=get(hedit(i,1),'string');
       Values{i}=str;
    else
        for j=1:length(str2num(Values{i}))
            if j==1,text=[text,':'];end
            str=get(hedit(i,1),'string');
            if ~isempty(str)
                values=str2num(str);
            else
                values=get(hedit(i,1),'value');
            end
        end
        if ~isempty(Values{i})
            Values{i}=num2str(values);
        end
    end
    Text{i}=text;
end
Text=str2mat(Text);
fid=fopen(parameterfile,'w');
for i=1:size(Text,1)
    fprintf(fid,'%s %s\n',Text(i,:),Values{i});
end
fl=fclose(fid);
close(h);
return