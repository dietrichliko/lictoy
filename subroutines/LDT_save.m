
function LDT_save(filename)

load results
save('-v7',filename);
disp(['Results written to ',filename]);

endfunction
