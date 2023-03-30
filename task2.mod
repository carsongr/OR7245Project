param num_beams; 
param num_rows >= 1, integer;
param num_cols >= 1, integer;
set ROWS := 1 .. num_rows; 
set COLUMNS := 1 .. num_cols; 
set BEAMS := 1 .. num_beams;


param beam_values {BEAMS, ROWS, COLUMNS} >= 0;

param tumor_values {ROWS, COLUMNS} >= 0;

set tumor_area := {j in ROWS, k in COLUMNS: tumor_values[j,k] > 0};

param critical_values {ROWS, COLUMNS} >= 0;

set critical_area := {j in ROWS, k in COLUMNS: critical_values[j,k] > 0};

param critical_max;

param tumor_min;
# binary parameter for critical region and tumor region
param a {j in ROWS, k in COLUMNS} = if critical_values[j,k] > 0 then 1 else 0;
param e {j in ROWS, k in COLUMNS} = if tumor_values[j,k] > 0 then 1 else 0;
# dosage scalar of each beam
var X {i in BEAMS} >= 0;
# slack variables
var S {j in ROWS, k in COLUMNS} >= 0;
var T {j in ROWS, k in COLUMNS} >= 0;

# minimize total dosage in critical area
minimize total_critical_dosage: sum {j in ROWS, k in COLUMNS} (S[j,k] + T[j,k])
;
# total dosage at each tumor location [j,k] should be >= min tumor dosage
subject to tumor_limit {j in ROWS, k in COLUMNS} : sum {i in BEAMS} X[i] * beam_values[i,j,k] * e[j,k] >= (tumor_min - T[j,k]) * e[j,k] ;

# total dosage at each critical location [j,k] should be <= max critical dosage 
subject to critical_limit {j in ROWS, k in COLUMNS} : sum {i in BEAMS} a[j,k] * X[i] * beam_values[i,j,k] <= critical_max + S[j,k];