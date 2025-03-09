import sys
from itertools import combinations
from os import listdir
from gurobipy import *
from scipy.special import comb
from math import log2
from heuristic_async import *
from os.path import splitext


steady_states=[]
trajectories=[]

def heuristic_callback(model, where):
       ss_values=[[] for x in steady_states]
       trj_values=[[[] for x in y] for y in trajectories]
       trj_vars=[]
       ss_vars=[]
       cur_edges={}
       if where == GRB.Callback.MIPNODE:
        if model.cbGet(GRB.Callback.MIPNODE_STATUS)==GRB.OPTIMAL:

            colsums = [0 for x in range(len(node_names))]
            numvals = 0

            for row in range(len(steady_states)):
                numvals+=1
                for col in range(len(steady_states[0])):
                    hv=model.getVarByName('SS_'+str(row)+'['+str(col)+']')
                    ss_values[row].append((round(model.cbGetNodeRel(hv))+steady_states[row][col])%2)
                    ss_vars.append(hv)
                    colsums[col] += ((round(model.cbGetNodeRel(hv)) + steady_states[row][col]) % 2)

            for number in range(len(trajectories)):
                for row in range(len(trajectories[number])):
                    numvals += 1
                    for col in range(len(trajectories[number][0])):
                         hv = model.getVarByName('TRJ_' + str(number) + '[' + str(row) + ','+str(col)+']')
                         trj_values[number][row].append((round(model.cbGetNodeRel(hv))+trajectories[number][row][col])%2)
                         trj_vars.append(hv)
                         colsums[col] += ((round(model.cbGetNodeRel(hv)) + trajectories[number][row][col]) % 2)

            for gene in node_names:   #for genes without regulators we want an empty list stored
                reg_idx=0
                gene_num=[i for i in range(len(node_names)) if node_names[i]==gene][0]
                if not gene in cur_edges.keys():
                    cur_edges[gene_num]=[]
                for nodenum,reg in enumerate(node_names):
                    if gene in edges.keys() and reg in edges[gene]:
                        hv=model.getVarByName('R_'+gene+'['+str(reg_idx)+']')
                        reg_idx+=1
                        if model.cbGetNodeRel(hv)>0.5:
                            cur_edges[gene_num].append(nodenum)

            no_reg_indices = []

            for col_itr, gene in enumerate(node_names):
                if not gene in cur_edges.keys():
                    no_reg_indices.append(col_itr)
                    colsums[col_itr] /= numvals
                    colsums[col_itr] = round(colsums[col_itr])

            for col in no_reg_indices:
                for number in range(len(trajectories)):
                  for row in range(len(trajectories[number])):
                        trj_values[number][row][col] = colsums[col]

            if len(ss_vars)>0 and len(trj_vars)>0:
                cur_logic=resolve_ss_discrepancies(ss_values,cur_edges,True)
                model.cbSetSolution(ss_vars,unpack_ss(ss_values,steady_states))
                resolve_trj_discrepancies(trj_values,cur_edges,cur_logic)
                model.cbSetSolution(trj_vars, unpack_trj(trj_values,trajectories))
            elif len(ss_vars)>0:
                cur_logic = resolve_ss_discrepancies(ss_values, cur_edges)
                model.cbSetSolution(ss_vars, unpack_ss(ss_values,steady_states))
            else:
                if different_length_trajectories:
                    resolve_trj_discrepancies(trj_values,cur_edges,dict())
                else:
                    resolve_trj_discrepancies(trj_values, cur_edges)
                model.cbSetSolution(trj_vars, unpack_trj(trj_values,trajectories))


#sys.argv[1] is the directory where the data files are, sys.argv[2] is the prefix of a data file name, sys.argv[3] is the path and name of
#the file that contains the names of the genes and possible edges,

m = Model("medsi")
#read nodes and edges
f=open(sys.argv[3],'r')
node_names=[x.strip() for x in f.readline().split()]
edges=dict()  #the key is the target, the value is a list of regulators
for line in f:
    if line != '\n':
        pair=line.split()
        regulator = pair[0].strip()
        target=pair[1].strip()
        if target in edges.keys():
            edges[target].append(regulator)
        else:
            edges[target]=[regulator]

f.close()
#read the data
data_dir=listdir(sys.argv[1])
trj_file_names=[]
for file in data_dir:
    if sys.argv[2] in file:
        f=open(sys.argv[1]+'/'+file,'r')
        line=f.readline()
        if 'steady state' in line:
            for line in f:
                if line!='\n':
                    steady_states.append([int(x) for x in line.split()])
        elif 'trajectory' in line:
            new_traj=[]
            for line in f:
                if line!='\n':
                    new_traj.append([int(x) for x in line.split()])
            trajectories.append(new_traj)
            trj_file_names.append(file)
        f.close()

different_length_trajectories=len(set(len(x) for x in trajectories))>1

B_vars=[]  #one 1D variable per steady state
B_vars_trajectories=[]  #one 2D variable per trajectory

for i,state in enumerate(steady_states):
    next_vars=m.addVars(len(state),vtype=GRB.BINARY,name='SS_'+str(i))
    B_vars.append(next_vars)

for i,trajectory in enumerate(trajectories):
    next_vars = m.addVars(len(trajectory),len(trajectory[0]), vtype=GRB.BINARY,name='TRJ_'+str(i))
    B_vars_trajectories.append(next_vars)
    

#add all constraints for steady states
I_vars=[]
R_vars=[]
V_vars=[]

if len(steady_states)>0:
    for target in edges.keys():
        w=[0 for i in range(len(edges[target]))]
        Is=[]
        ws=[]
        reg_indices=[i for i in range(len(node_names)) if node_names[i] in edges[target]]
        w_dict={reg_indices[i]:i for i in range(len(reg_indices))}  #map regulator indices in the input data to indices in w
        target_index=[i for i,x in enumerate(node_names) if x==target][0]
        for i in range(2**len(w)):
            Iw=m.addVar(name='I_'+target+str(w),vtype=GRB.BINARY)
            Is.append(Iw)
            ws.append(list(w))
            for k in range(len(steady_states)):
                C=steady_states[k]
                B=B_vars[k]

                m.addConstr(quicksum(C[j]*(w[w_dict[j]]+(1-2*w[w_dict[j]])*B[j])+(1-C[j])*(1-w[w_dict[j]]+(2*w[w_dict[j]]-1)*B[j]) for j in reg_indices) \
                    +C[target_index]*B[target_index]+(1-C[target_index])*(1-B[target_index])<=(2-Iw)*(len(w)+1)-1)

                m.addConstr(quicksum(C[j]*(w[w_dict[j]]+(1-2*w[w_dict[j]])*B[j])+(1-C[j])*((1-w[w_dict[j]])+(2*w[w_dict[j]]-1)*B[j]) for j in reg_indices) \
                    +C[target_index]*(1-B[target_index])+(1-C[target_index])*B[target_index]<=(Iw+1)*(len(w)+1)-1)

            o=0
            while o<len(w)-1 and w[o]==1:
                w[o]=0
                o+=1
            w[o]=1

        Rs=m.addVars(len(w),vtype=GRB.BINARY,name='R_'+target)
        R_vars.append(Rs)  #save the lists of R variables of each target for use in the objective
        pairs=combinations(range(len(Is)),2)
        for i,j in pairs:
            changingRs=[Rs[k] for k in range(len(w)) if ws[i][k]!=ws[j][k]]
            m.addConstr(quicksum(changingRs)>=Is[i]-Is[j])
            m.addConstr(quicksum(changingRs)>=Is[j]-Is[i])

        Vs=m.addVars(len(w),vtype=GRB.BINARY)
        for i in range(len(w)):
            if i==0:
                nom=0.5
            else:
                nom=i
            m.addConstr(Vs[i]>=quicksum(x for x in Rs.values())/len(w)-nom/len(w))

        V_vars.append(Vs)  #save the lists of V variables of each target for use in the objective
        I_vars.append(Is)

D_vars=[]
#add all constraints for trajectories - have to re-use I_vars and R_vars if their list length is greater than 0

if len(trajectories)>0:
    for target_itr,target in enumerate(edges.keys()):
        w=[0 for i in range(len(edges[target]))]
        Is=[]
        ws=[]
        reg_indices = [i for i in range(len(node_names)) if node_names[i] in edges[target]]
        w_dict={reg_indices[i]:i for i in range(len(reg_indices))}  #map regulator indices in the input data to indices in w
        target_index = [i for i, x in enumerate(node_names) if x == target][0]
        first_Ds = []
        save_Ds = {}
        for i in range(2**len(w)):
            if len(I_vars)==0:
              Iw=m.addVar(name='I_'+target+str(w),vtype=GRB.BINARY)
              Is.append(Iw)
            else:
              Iw=I_vars[target_itr][i]
            ws.append(list(w))
            for k in range(len(trajectories)):
                C = trajectories[k]
                B = B_vars_trajectories[k]

                for l in range(len(C)-1):
                    if i==0:
                      D_1 = m.addVar(name='D1_' + target + str(k)+'-'+str(l), vtype=GRB.BINARY)
                      D_vars.append(D_1)
                      save_Ds[(k,l)]=D_1
                    else:
                      D_1=save_Ds[(k,l)]
                    if l==0:
                        first_Ds.append(D_1)

                    m.addConstr(1 - (C[l + 1][target_index] * (1 - B[l + 1, target_index]) + ( 1 - C[l + 1][target_index]) * B[l + 1, target_index] - ( C[l][target_index] * (1 - B[l, target_index]) + (1 - C[l][target_index]) * B[l, target_index])) >= D_1)
                    m.addConstr(1 - (C[l][target_index] * (1 - B[l, target_index]) + (1 - C[l][target_index]) * B[l, target_index] - (C[l + 1][target_index] * (1 - B[l + 1, target_index]) + (1 - C[l + 1][target_index]) * B[l + 1, target_index])) >= D_1)

                    m.addConstr(quicksum(C[l][j] * (w[w_dict[j]] + (1 - 2 * w[w_dict[j]]) * B[l,j]) + (1 - C[l][j]) * (1 - w[w_dict[j]] + (2 * w[w_dict[j]] - 1) * B[l,j]) for j in reg_indices) \
                                + C[l+1][target_index] * B[l+1,target_index] + (1 - C[l+1][target_index]) * (1 - B[l+1,target_index]) <= (2 - Iw) * (len(w) + 1) - 1+D_1)

                    m.addConstr(quicksum(
                        C[l][j] * (w[w_dict[j]] + (1 - 2 * w[w_dict[j]]) * B[l,j]) + (1 - C[l][j]) * ((1 - w[w_dict[j]]) + (2 * w[w_dict[j]] - 1) * B[l,j]) for j in reg_indices) \
                                + C[l+1][target_index] * (1 - B[l+1,target_index]) + (1 - C[l+1][target_index]) * B[l+1,target_index] <= (Iw + 1) * (len(w) + 1) - 1+D_1)

            o = 0
            while o < len(w) - 1 and w[o] == 1:
                w[o] = 0
                o += 1
            w[o] = 1

        if len(steady_states)==0:
            Rs = m.addVars(len(w), vtype=GRB.BINARY,name='R_'+target)
            R_vars.append(Rs)  # save the lists of R variables of each target for use in the objective
            pairs = combinations(range(len(Is)), 2)
            for i, j in pairs:
                changingRs = [Rs[k] for k in range(len(w)) if ws[i][k] != ws[j][k]]
                m.addConstr(quicksum(changingRs) >= Is[i] - Is[j])
                m.addConstr(quicksum(changingRs) >= Is[j] - Is[i])
            for var in first_Ds:
                m.addConstr(var<=sum(Rs))
            Vs = m.addVars(len(w), vtype=GRB.BINARY)
            for i in range(len(w)):
                if i == 0:
                    nom = 0.5
                else:
                    nom = i
                m.addConstr(Vs[i] >= quicksum(x for x in Rs.values()) / len(w) - nom / len(w))

            V_vars.append(Vs)  # save the lists of V variables of each target for use in the objective

#add objective
def nonred(n):   #returns the log of the number of non-redundant Boolean functions in n variables
    if n==0:
        return 0
    if n<=5:
        return log2(sum(pow(-1,n-k)*comb(n,k)*pow(2,pow(2,k)) for k in range(n,0,-1)))
    return 2**n

def regbits(N,i):  #returns the number of bits to add to the i-th V variable for regulator edges, i starting from 1 (first V variable)

    if i==0:
        return 0

    return log2((N-i+1)/i)

m.setObjective(sum(x for y in B_vars for x in y.values())+sum(x for y in B_vars_trajectories for x in y.values())+\
               sum(D_vars[i] for i in range(len(D_vars)))+sum((nonred(i)-nonred(i-1)+regbits(len(node_names),i))*Vs[i-1] for Vs in V_vars for i in range(1,len(Vs)+1)), GRB.MINIMIZE)

m.setParam('OBBT',3)
m.optimize(heuristic_callback)

#print inferred edges#
chosen_edges=dict()  #like edges but only ones that were chosen

f=open('InferenceResults.txt','w')
for v in m.getVars():
    if 'R_' in v.VarName and v.X==1:
        target=v.VarName[2:].split('[')[0]
        reg_list=[x for x in node_names if x in edges[target]]
        regulator=reg_list[int(v.VarName[2:].split('[')[1].split(']')[0])]
        if target in chosen_edges.keys():
            chosen_edges[target].append(regulator)
        else:
            chosen_edges[target]=[regulator]

for target in chosen_edges.keys():
    for regulator in chosen_edges[target]:
        f.write(regulator + '\t' + target + '\n')

#print inferred logic for genes with at least one regulator
for target in chosen_edges.keys():
    f.write('Regulatory Logic for '+target+'\n')
    covered_inputs=[]
    w = [0 for i in range(len(edges[target]))]
    reg_indices = [i for i in range(len(node_names)) if node_names[i] in chosen_edges[target]]
    f.write('Regulators: '+'\t'.join([node_names[i] for i in reg_indices])+'\n')
    original_reg_indices = [i for i in range(len(node_names)) if node_names[i] in edges[target]]
    w_indices=[i for i in range(len(w)) if original_reg_indices[i] in reg_indices]
    for i in range(2**len(edges[target])):
        logic_output='I_' + target + str(w)
        next_input=[w[i] for i in w_indices]
        if not next_input in covered_inputs:
            f.write(''.join([str(next_input[i]) for i in range(len(next_input))])+'\t'+str(int(m.getVarByName(logic_output).X))+'\n')
            covered_inputs.append(next_input)
        o = 0
        while o < len(w) - 1 and w[o] == 1:
            w[o] = 0
            o += 1
        w[o] = 1

    f.write('\n')

f.close()



diffs=0
total=0

for l in range(len(trajectories)):
    f = open(splitext(trj_file_names[l])[0]+'_modeled.txt', 'w')
    for i in range(len(trajectories[l])):
        f.write('\n')
        for j in range(len(trajectories[l][i])):
            f.write(str((round(m.getVarByName('TRJ_'+str(l)+'['+str(i)+','+str(j)+']').X+trajectories[l][i][j])%2))+'\t')
            diffs+=round(m.getVarByName('TRJ_'+str(l)+'['+str(i)+','+str(j)+']').X)==1
            total+=1
    f.close()


print('Percentage mismatches:',diffs/total)

for l in range(len(steady_states)):
    for i in range(len(steady_states[l])):
        diffs += round(m.getVarByName('SS_' + str(l) + '[' + str(i) + ']').X) == 1
        total += 1

