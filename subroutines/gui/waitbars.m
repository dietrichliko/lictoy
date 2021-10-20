function waitbars

global N
global hw_track hw_event

if N.Event>1
        hw_event=waitbar(0,' ');
        set(hw_event,'name','Event progress...')
        %ch=get(hw_event,'children');
        %set(ch(2),'String','Interrupt')
        drawnow
    end
    hw_track=waitbar(0,' ','CreateCancelBtn','global interrupt;interrupt=1;');
    if N.Event>1
        pos=get(hw_event,'position');
        pos=pos+[50 -100 0 25];
        actpos=zeros(size(pos));
        count=0;
        maxcount=10;
        while any(actpos~=pos) & count<maxcount % Strange behaviour of MATLAB
            set(hw_track,'name','Track progress...','position',pos)
            drawnow
            actpos=get(hw_track,'position');
            count=count+1;
        end
    else
        set(hw_track,'name','Track progress...')
        drawnow
    end
    ch=get(hw_track,'children');
    set(ch(2),'String','Interrupt')