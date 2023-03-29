from docplex.mp.model import Model
import os
# input file name is set as mat_raw.txt. remember to change it if you use a different file
class modeling_data():
    def __init__(self, file_path):
        self.data_file_name = file_path

        self.num_matrices = None
        self.num_rows = None
        self.num_columns = None

        self.matrices = [[]]

    def get_data_from_file(self):
        file = open(self.data_file_name, 'r')
        lines = file.readlines()
        for line in lines:
            l = line[:-1].split()
            if l:
                self.matrices[-1].append([])
                for i in l:
                    self.matrices[-1][-1].append(float(i))
            else:
                self.matrices.append([])

        self.num_matrices = len(self.matrices)
        self.num_rows = len(self.matrices[0])
        self.num_columns = len(self.matrices[0][0])

    #aggregates by sum multiple matrices into a single matrix
    def combine_matrices(self):
        result = [[]]
        for k in range(self.num_matrices):
            mat = self.matrices[k]
            if result == [[]]:
                result = mat
            else:
                # iterate through rows
                for i in range(len(mat)):  
                # iterate through columns
                    for j in range(len(mat[0])):
                        result[i][j] = mat[i][j] + result[i][j]

        self.matrices = [result]

        self.num_matrices = len(self.matrices)
        self.num_rows = len(self.matrices[0])
        self.num_columns = len(self.matrices[0][0])

        return self                    
         
"""
critical_region: matrix of critical region 
tumor_region: matrix of tumor_region
beam_region: matrix of beam region - the beam matrices should be combined into a single aggregate matrix 
"""
def build_model(critical_region, tumor_region, beam_region, lower_limit, upper_limit, critical_region_balancing_param=1, tumor_region_balancing_param=-1):
    model = Model(log_output=True)

    # Variables: x_i_j is the maximum value in the matrices
    x = model.continuous_var_matrix(keys1=beam_region.num_rows, keys2=beam_region.num_columns, name="x")
    print(beam_region.matrices)
    for i in range(beam_region.num_rows):
        for j in range(beam_region.num_columns):
            for k in range(beam_region.num_matrices):
                #constraints for critical region upper bound
                model.add_constraint(critical_region.matrices[k][i][j] * beam_region.matrices[k][i][j] * x[i,j] <= upper_limit)

                #constraints for tumor region lower bound
                model.add_constraint(tumor_region.matrices[k][i][j] * beam_region.matrices[k][i][j] * x[i,j] >= lower_limit)

    print(model.export_as_lp_string())

    e = model.linear_expr()
    for i in range(tumor_region.num_rows):
        for j in range(tumor_region.num_columns):
            e += (tumor_region_balancing_param * tumor_region.matrices[0][i][j] * x[i,j]) + (critical_region_balancing_param * critical_region.matrices[0][i][j] * x[i,j])
    model.minimize(e)

    return model

def print_result(filename, vars, data):
    file = open(filename, 'w')
    for i in range(data.num_rows):
        for j in range(data.num_columns):
            file.write(str(vars[i * data.num_rows + j].solution_value) + '\t')
        file.write('\n')


critical_region = modeling_data('/Users/carson/Documents/Courses/Advanced Optimization and Network Analysis/Project/Workspace/critical_raw.txt')
critical_region.get_data_from_file()

tumor_region = modeling_data('/Users/carson/Documents/Courses/Advanced Optimization and Network Analysis/Project/Workspace/tumor_raw.txt')
tumor_region.get_data_from_file()

beam_region = modeling_data('/Users/carson/Documents/Courses/Advanced Optimization and Network Analysis/Project/Workspace/beam_raw.txt')
beam_region.get_data_from_file()

agg_beam_region = beam_region.combine_matrices()

print(beam_region.matrices)
model = build_model(critical_region, tumor_region, agg_beam_region, lower_limit=0, upper_limit=20)

s = model.solve()
#x_vars = model.find_matching_vars(pattern="x_")
#print_result("results.out", x_vars, data)
