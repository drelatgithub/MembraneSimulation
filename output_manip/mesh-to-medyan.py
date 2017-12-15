import os

cur_dir = os.path.split(os.path.abspath(__file__))[0]

file_in_dir = cur_dir + "\\..\\build\\VC\\MeshGeneration\\"
file_in = [
    file_in_dir + "position.txt",
    file_in_dir + "neighbors.txt"
]

file_out_dir = cur_dir + "\\"
file_out = file_out_dir + "sphere.txt"

len_ref = 1e9
coo_mov = [1500, 1500, 1500]

with open(file_in[0], 'r') as pos, open(file_in[1], 'r') as nei, open(file_out, 'w') as out:
    for pos_line, nei_line in zip(pos, nei):
        for each_num, each_mov in zip(pos_line.strip().split('\t'), coo_mov):
            out.write(str(float(each_num) * len_ref + each_mov) + '\t')
        out.write(nei_line)
