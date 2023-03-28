param num_beams; # number of beams
param num_rows >= 1, integer; # number of rows 
param num_cols >= 1, integer; # number of columns


read {m in BEAMS, i in ROWS, j in COLUMNS} beam_values[m,i,j] < beam_raw.txt;

read {i in ROWS, j in COLUMNS} tumor_values[i,j] < tumor_raw.txt; # read in tumor data

read {i in ROWS, j in COLUMNS} critical_values[i,j] < critical_raw.txt; # read in critical area data
set ROWS := 1 .. num_rows; # set of rows 
set COLUMNS := 1 .. num_cols; # set of columns 
set BEAMS := 1 .. num_beams; # set of beams




param beam_values {BEAMS, ROWS, COLUMNS} >= 0;


param tumor_values {ROWS, COLUMNS} >= 0; # values of tumor




# define the tumor area
set tumor_area := {j in ROWS, k in COLUMNS: tumor_values[j,k] > 0};
# values of critical area
param critical_values {ROWS, COLUMNS} >= 0;
# define critical area
set critical_area := {j in ROWS, k in COLUMNS: critical_values[j,k] > 0};
# critical maximum dosage requirement
param critical_max;
# tumor minimum dosage requirement
param tumor_min;





# binary parameter for critical region and tumor region
param a {j in ROWS, k in COLUMNS} = if critical_values[j,k] > 0 then 1 else 0;
param e {j in ROWS, k in COLUMNS} = if tumor_values[j,k] > 0 then 1 else 0;


# dosage scalar of each beam
var X {i in BEAMS} >= 0;



# minimize total dosage in critical area
minimize total_critical_dosage: sum {i in BEAMS} sum {j in ROWS, k in COLUMNS}
   a[j,k] * X[i] * beam_values[i,j,k];
   
   
# total dosage at each tumor location [j,k] should be >= min tumor dosage with
   slack variable
subject to tumor_limit {j in ROWS, k in COLUMNS} : sum {i in BEAMS} X[i] *
   beam_values[i,j,k] >= tumor_min * e[j,k];
# total dosage at each critical location [j,k] should be <= max critical dosage
    with slack variable
subject to critical_limit {j in ROWS, k in COLUMNS} : sum {i in BEAMS} a[j,k] *
    X[i] * beam_values[i,j,k] <= critical_max;