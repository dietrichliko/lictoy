function LDT_Octave

clear
global unit octave
warning off;
unit=1;
octave=1;
path(path,'subroutines');
%path(path,'subroutines/gui');
path(path,'subroutines/outputs');
path(path,'subroutines/simulation');
path(path,'subroutines/tests');
edit('simulation_parameters_Octave.txt')
disp(' ');
disp('             __    ___   ______             ');
disp('            (  )  (   \ (__  __)            ');
disp('             )(__  ) > )   )(    v2.0       ');
disp('            (____)(___/   (__)              ');
disp('         The Vienna Fast Simulation Tool    ');
disp(' ');
disp('Ready!');
disp(' ');