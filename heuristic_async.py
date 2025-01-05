from sklearn.cluster import BisectingKMeans
import numpy as np
from scipy.spatial.distance import hamming

def are_trjs_discrepant(trajectories, edges):  #reg_nums may be an empty list, meaning the target has no regulators
    finished=True
    tmp_logic={}
    for trj in trajectories:
        for state_itr in range(len(trj) - 1):
            for target_num in edges.keys():
                reg_nums = edges[target_num]
                k=(target_num, tuple(trj[state_itr][i] for i in reg_nums))
                if k in tmp_logic.keys():
                    if tmp_logic[k] != trj[state_itr + 1][target_num]  and trj[state_itr+1][target_num]!=trj[state_itr][target_num]:
                        finished = False
                        break
                else:
                    tmp_logic[k] = trj[state_itr + 1][target_num]
            if not finished:
                break
        if not finished:
            break

    return not finished

def are_ss_discrepant(steady_states,edges):  #reg_nums may be an empty list, meaning the target has no regulators
    finished=True
    tmp_logic=dict()
    for s in steady_states:
         for target_num in edges.keys():
             reg_nums=edges[target_num]
             k=(target_num,tuple(s[i] for i in reg_nums))
             if k in tmp_logic.keys():
                 if tmp_logic[k]!=s[target_num]:
                     finished=False
                     break
             else:
                 tmp_logic[k]=s[target_num]
         if not finished:
            break

    return not finished

def add_ss_to_dict(steady_states,edges,res):
    for s in steady_states:
        for target_num in edges.keys():
            reg_nums = edges[target_num]
            k = (target_num, tuple(s[i] for i in reg_nums))
            res[k] = s[target_num]

def add_trjs_to_dict(trajectories,edges,res):
   for trj in trajectories:
       for state_itr in range(len(trj) - 1):
           for target_num in edges.keys():
               reg_nums = edges[target_num]
               k = (target_num, tuple(trj[state_itr][i] for i in reg_nums))
               res[k]=trj[state_itr + 1][target_num]

def is_ss_discrepant(s,edges,logic):  #reg_nums may be an empty list, meaning the target has no regulators

        for target_num in edges.keys():
            reg_nums=edges[target_num]
            k=(target_num,tuple(s[i] for i in reg_nums))
            if k in logic.keys():
                 if logic[k]!=s[target_num]:
                        return True
        return False

#discrepancies between two different trajectories, if run with the same trajectory as both trj1 and trj2 looks for discrepancies
#within a trajectory
def is_trj_discrepant(trj, edges,logic):  #reg_nums may be an empty list, meaning the target has no regulators

        own_logic=dict()
        for state_itr in range(len(trj) - 1):
            for target_num in edges.keys():
                reg_nums = edges[target_num]
                k=(target_num, tuple(trj[state_itr][i] for i in reg_nums))
                if k in logic.keys():
                    if logic[k] != trj[state_itr + 1][target_num] and trj[state_itr+1][target_num]!=trj[state_itr][target_num]:
                          return True
                elif k in own_logic.keys():
                    if own_logic[k] != trj[state_itr + 1][target_num] and trj[state_itr+1][target_num]!=trj[state_itr][target_num]:
                          return True
                else:
                    own_logic[k]=trj[state_itr + 1][target_num]

        return False


def cluster_trj(values,num_trj,trj_len,num_genes):  #return either a 2D or a 3D array of binary cluster centers, depending on the input
    bisect_means = BisectingKMeans(n_clusters=max([len(values)//4,1]), random_state=0).fit(np.array(values).reshape(num_trj,trj_len*num_genes))
    return [np.round(x).reshape(trj_len,num_genes).tolist() for x in bisect_means.cluster_centers_]

def cluster_ss(values):  #return either a 2D or a 3D array of binary cluster centers, depending on the input
    bisect_means = BisectingKMeans(n_clusters=max([len(values)//4,1]), random_state=0).fit(values)
    return [[np.round(x).tolist() for x in y] for y in bisect_means.cluster_centers_]

def resolve_trj_disc_recursive(trj_values, edges):

    if not are_trjs_discrepant(trj_values,edges):
        return

    if len(trj_values)==1:
        cur_logic=dict()
        trj=trj_values[0]
        for trj_itr in range(len(trj) - 1):
            for target_num in edges.keys():
                regulator_values = tuple(trj[trj_itr][i] for i in edges[target_num])
                if (target_num, regulator_values) in cur_logic.keys() :
                    if trj[trj_itr+1][target_num]!=trj[trj_itr][target_num]:
                        trj[trj_itr + 1][target_num] = cur_logic[(target_num, regulator_values)]
                elif trj[trj_itr+1][target_num]!=trj[trj_itr][target_num]:
                    cur_logic[(target_num, regulator_values)] = trj[trj_itr + 1][target_num]
        return trj_values

    cluster = cluster_trj(trj_values,len(trj_values),len(trj_values[0]),len(trj_values[0][0]))

    resolve_trj_disc_recursive(cluster, edges)
    logic={}
    add_trjs_to_dict(cluster, edges,logic)

    for trj in trj_values:
        c_idx=0
        min_dist=sum(hamming(trj[k],cluster[0][k]) for k in range(len(trj)))
        for i, c in enumerate(cluster[1:len(cluster)]):
            new_dist=sum(hamming(trj[k],c[k]) for k in range(len(trj)))
            if new_dist< min_dist:
                min_dist = new_dist
                c_idx = i
        not_discrepant=False

        for i in range(len(trj)-1):
            for j in range(len(trj[0])):
                    disc=is_trj_discrepant(trj, edges,logic)
                    if trj[i][j]!=cluster[c_idx][i][j] and disc:
                            trj[i][j] = cluster[c_idx][i][j]
                    elif not disc:
                        not_discrepant=True
                        cluster.append(trj)
                        add_trjs_to_dict([trj], edges, logic)
                        break
            if not_discrepant:
                break

def resolve_ss_disc_recursive(ss_values,edges):

    if not are_ss_discrepant(ss_values,edges):
        return

    cluster=cluster_ss(ss_values)

    resolve_ss_disc_recursive(cluster,edges)

    logic={}

    add_ss_to_dict(cluster,edges, logic)

    for s in ss_values:
        c_idx=0
        min_dist=hamming(s,cluster[0])
        for i,c in enumerate(cluster[1:len(cluster)]):
            new_dist=hamming(s,c)
            if new_dist<min_dist:
                min_dist=new_dist
                c_idx=i
        for i in range(len(s)):
            disc=is_ss_discrepant(s,edges,logic)
            if s[i]!=cluster[c_idx][i] and disc:
                s[i]=cluster[c_idx][i]
            elif not disc:
                cluster.append(s)
                add_ss_to_dict([s],edges,logic)
                break


def extract_logic(ss_values,edges):
    res=dict()
    for target in edges.keys():
        for s in ss_values:
            res[(target,tuple(s[i] for i in edges[target]))]=s[target]

    return res

def resolve_ss_discrepancies(ss_values,edges,get_logic=False):
    resolve_ss_disc_recursive(ss_values, edges)
    if get_logic:
       return extract_logic(ss_values,edges)


def resolve_trj_discrepancies(trj_values,edges,cur_logic=None):

    if not cur_logic is None:
        for trj in trj_values:
            for trj_itr in range(len(trj)-1):
                for target_num in edges.keys():
                    regulator_values=tuple(trj[trj_itr][i] for i in edges[target_num])
                    if (target_num,regulator_values) in cur_logic.keys():
                        if trj[trj_itr+1][target_num]!=trj[trj_itr][target_num]:
                            trj[trj_itr+1][target_num]=cur_logic[(target_num,regulator_values)]
                    elif trj[trj_itr+1][target_num]!=trj[trj_itr][target_num]:
                        cur_logic[(target_num,regulator_values)]=trj[trj_itr+1][target_num]
        return

    resolve_trj_disc_recursive(trj_values, edges)


def unpack_ss(ss_values,steady_states):
    res=[]
    for row in range(len(ss_values)):
        for col in range(len(ss_values[0])):
            res.append((ss_values[row][col]+steady_states[row][col])%2)

    return res


def unpack_trj(trj_values,trajectories):
    res=[]
    for number in range(len(trj_values)):
        for row in range(len(trj_values[number])):
            for col in range(len(trj_values[number][0])):
                    res.append((trj_values[number][row][col]+trajectories[number][row][col])%2)
    return res
