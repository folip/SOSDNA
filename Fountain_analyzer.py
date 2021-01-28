from Encode.Helper_Functions import dna_to_int_array, rs_decode, preprocess, load_dna
from Analysis.Analysis import dna_chunk, save_simu_result
from Encode.DNAFountain import DNAFountain, Glass
import numpy as np
from scipy.stats import gumbel_r, poisson
import matplotlib.pyplot as plt

def error_profile(out_dnas, rs = 2):
    lost_num = 0
    fail_num = 0
    mis_judge_num = 0
    ignore_index = []
    for i,dna in enumerate(out_dnas):
        
        if dna['num'] is 0: 
            lost_num += 1
            continue
        
        dc = dna_chunk(dna)
        re_dna = dc.voting_result()
        
        re_data = dna_to_int_array(re_dna)
        flag, data_corrected = rs_decode(re_data, rs = rs)

        if flag == -1:
            fail_num += 1
        else:
            ori_data = dna_to_int_array(dna['ori'])
            if ori_data[:-2] != data_corrected:
                mis_judge_num += 1
                fail_num += 1
                ignore_index.append(i)
    return lost_num, fail_num, mis_judge_num, ignore_index

class FT_Analyzer:
    def __init__(self, file_name, model, alpha, rs_length, chunk_size = 20, encode = True):
        # path
        file_name, file_type = file_name.split('.')
        self.ori_path = 'files/' + file_name + '.' + file_type
        self.encode_path = 'files/' + file_name + '.dna'
        self.simu_path = 'files/' + 'simu_' + file_name + '.dna'
        self.simu_decode_path = 'files/' + 'simu_re_' + file_name + '.' + file_type
        # data
        self.chunk_size = chunk_size
        self.data, self.pad = preprocess(self.ori_path, chunk_size = chunk_size)
        self.chunk_num = len(self.data)
        print(f'Data split into {len(self.data)} data chunks.')
        # channel model
        self.model = model
        # encoder-decoder
        self.alpha = alpha
        self.rs_length = rs_length
        # distribution
        self.loss_nums = []
        self.fail_nums = []
        self.decode_lines = []
        if encode: self.encode()

    def encode(self):
        self.f = DNAFountain(self.data, self.alpha, rs = self.rs_length)
        good, tries = self.f.encode()
        self.good = good
        print(f'Encoding process done. Generated {good} strands after {tries} tries.')
        print(f'Alpha: {self.alpha}, chunk_size: {self.chunk_size}, rs_length: {self.rs_length}, information_density: {round(self.chunk_size/(self.rs_length + self.chunk_size)/(1+self.alpha),4)}')
        self.f.save(self.encode_path)

    def simu(self, save = True):
        in_dnas = load_dna(self.encode_path)
        out_dnas = self.model(in_dnas, print_state = False)
        loss_num, fail_num, mis_judge_num, ignore_index = error_profile(out_dnas)
        self.loss_nums.append(loss_num)
        self.fail_nums.append(fail_num)
        if save: save_simu_result(out_dnas, self.simu_path, ignore_index)
        self.out_dnas = out_dnas
        return out_dnas
    
    def decode(self, save = True):
        self.g = Glass(self.simu_path, len(self.data), rs = self.rs_length)
        flag, solve_num, line, chunksDone, errors =  self.g.decode()
        print('Decoding...', end = '\r')
        print("After reading %d lines, %d chunks are done. So far: %d rejections (%f)" % (line, chunksDone, errors, errors/(line+0.0)))
        if flag == -1: 
            print('Falied.')
            return -1
        else:
            print('Success!')
            if save: self.g.save(self.simu_decode_path)
            return line - errors
    
    def run(self, plot_res = False):
        out_dnas = self.simu()
        line = self.decode()
        if line == -1: 
            return
        self.decode_lines.append(line)
    
    def compute_dist(self):
        NP = np.mean(self.loss_nums) + np.mean(self.fail_nums)
        self.p = NP / (1+self.alpha) / self.chunk_num
        self.loc, self.scale = gumbel_r.fit(self.decode_lines)
        
    def fail_prob(self,alpha,plot = False, show_E = True):
        self.compute_dist()
        X = np.arange(0,int(alpha * self.chunk_num))
        
        X_decode = int((1+alpha) * self.chunk_num) - X
        P_correct = gumbel_r.pdf(X_decode,self.loc, self.scale)
        P_correct = P_correct / np.sum(P_correct)

        u_error = int((1+alpha) * self.chunk_num) * self.p
        P_error = poisson.pmf(X,u_error)
        
        p_fail_to_correct = 0
        for i in X:
            p_fail_to_correct += P_error[i] * np.sum(P_correct[:i+1])

        if plot: 
            no_zeros = np.where( (P_correct + P_error) > 0.00001)
            s,e = np.min(no_zeros), np.max(no_zeros)
            X = X[s:e]
            P_error = P_error[s:e]
            P_correct = P_correct[s:e]
            if show_E:
                plt.plot(X,P_error,color = 'red',label = f'generated losts')
                plt.fill_between(X,0,P_error,alpha = 0.3, color = 'red')

            plt.plot(X,P_correct,color = 'blue',label = f'alpha = {alpha}\nfail prob = {round(p_fail_to_correct*100,2)}%')
            plt.fill_between(X,0,P_correct,alpha = 0.3, color = 'blue')

            plt.legend()
            # plt.show()

        return p_fail_to_correct
    
    def alpha_scan(self, alpha_list = None, plot = True, points = 15, color = 'blue'):
        if not alpha_list:
            alpha_list = np.linspace(self.alpha-0.15,self.alpha+0.1,points)
        p_fail_list = []
        # info_ratio_list = []
        for alpha in alpha_list:
            p_fail = self.fail_prob(alpha)
            # info_ratio = round(self.chunk_size/(self.rs_length + self.chunk_size)/(1+alpha),4)
            p_fail_list.append(p_fail) 
            # info_ratio_list.append(info_ratio)
        
        if plot:
            plt.plot(alpha_list,p_fail_list,label = 'fail prob',color = color)
            plt.fill_between(alpha_list,0,p_fail_list,alpha = 0.3, color = color)
            plt.xlabel('Alpha')
            # plt.plot(alpha_list,info_ratio_list,label = 'information density')
            plt.legend()
        
        return p_fail_list

        


        
    
    
    
        
        

            