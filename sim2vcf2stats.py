#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 13 17:25:43 2021

@author: xinjunzhang
"""

import msprime, pyslim, os, random, itertools,argparse
import numpy as np
import pandas as pd
import allel

parser = argparse.ArgumentParser(description="A script for computing summary statistics in 50kb windows across given chromosome, between modern human and archaic human.")
parser.add_argument('-s', '--sim', action="store", dest="sim_id",
                        help="which simulation batch, default: 1",
                        default=1, type=int)
parser.add_argument('-g', '--gene', action="store", dest="gene_id",
                        help="which genomic region, default: 1; range 1-1000",
                        default=1, type=int)
parser.add_argument('-d', '--dominance', action="store", dest="dominance_id",
                        help="which dominance to use, default: 1; range 1-3",
                        default=1, type=int)
args = parser.parse_args()

batch = str(args.sim_id) #now sim id is going to refer to a 
whichgene = int(args.gene_id) 
whichh = int(args.dominance_id) 

######################################################################
def simplify_tree(filename, sampsize, archaic_sample_time):
    '''
    this method subsamples a tree sequence to a given set of
    sample sizes. archaic sample is set to be 1 individual in SLiM
    '''
    #load tree sequence
    ts = pyslim.load(filename)
    indivs = [x for x in ts.individuals()]
    #get nodes (chromosomes) alive today
    nodes_today_p1 = [x.id for x in ts.nodes() if ((x.time == 0.0) and
                      (x.population == 1))]
    nodes_today_p4 = [x.id for x in ts.nodes() if ((x.time == 0.0) and
                      (x.population == 4))]
    nodes_today_p3 = [x.id for x in ts.nodes() if ((x.time == 0.0) and
                      (x.population == 3))]
    #get a list of the individuals alive today
    indivs_today_p1 = [x.id for x in indivs if x.nodes[0] in
                       nodes_today_p1]
    indivs_today_p4 = [x.id for x in indivs if x.nodes[0] in
                       nodes_today_p4]
    indivs_today_p3 = [x.id for x in indivs if x.nodes[0] in
                       nodes_today_p3]
    #subsample individuals to sample sizes
    indivs_sample_p1 = random.sample(indivs_today_p1, sampsize[0])
    indivs_sample_p3 = random.sample(indivs_today_p3, sampsize[1])
    indivs_sample_p4 = random.sample(indivs_today_p4, sampsize[2])
    #get their nodes
    nodes_sample_p1 = [x.nodes for x in indivs if x.id in
                       indivs_sample_p1]
    nodes_sample_p4 = [x.nodes for x in indivs if x.id in
                       indivs_sample_p4]
    nodes_sample_p3 = [x.nodes for x in indivs if x.id in
                       indivs_sample_p3]
    nodes_sample_p1 = list(itertools.chain.from_iterable(
                           nodes_sample_p1))
    nodes_sample_p4 = list(itertools.chain.from_iterable(
                           nodes_sample_p4))
    nodes_sample_p3 = list(itertools.chain.from_iterable(
                           nodes_sample_p3))
    #get archaic samples
    archaic = [x.id for x in ts.nodes() if x.time ==
               archaic_sample_time]
    indivs_archaic = [x.id for x in indivs if x.nodes[0] in archaic]
    nodes_sample_archaic = [x.nodes for x in indivs if x.id in
                            indivs_archaic]
    nodes_sample_archaic = list(itertools.chain.from_iterable(
                                nodes_sample_archaic))
    #subsample while retaining admixture recorded nodes
    samp = nodes_sample_archaic + nodes_sample_p1 + nodes_sample_p3+nodes_sample_p4
    ts_sample = ts.simplify(samples=samp, filter_populations=False)
    return ts_sample

def ns_vcf(ts_sample, archaic_sample_time,vcfoutpath):
    with open(vcfoutpath, "w") as vcf_file:
        #these_indivs = np.append(ts_sample.individuals_alive_at(archaic_sample_time),ts_sample.individuals_alive_at(0))
        ts_sample.write_vcf(vcf_file,ploidy=2)

def get_selcoeff(ts):
    coeff = []
    muts = []
    for mut in ts.mutations():
        site = mut.position
        md = pyslim.decode_mutation(mut.metadata)
        sel = [x.selection_coeff for x in md] 
        coeff.append(sel[0])
        muts.append(site)
    return coeff, muts
#coeffs = []
#for mut in ts.mutations():
#    md = pyslim.decode_mutation(mut.metadata)
#    sel = [x.selection_coeff for x in md] 

def remove_mutations(ts, start, end, proportion):
    '''
    This function will return a new tree sequence the same as the input,
    but after removing each non-SLiM mutation within regions specified in lists
    start and end with probability `proportion`, independently. So then, if we
    want to add neutral mutations with rate 1.0e-8 within the regions and 0.7e-8
    outside the regions, we could do
      ts = pyslim.load("my.trees")
      first_mut_ts = msprime.mutate(ts, rate=1e-8)
      mut_ts = remove_mutations(first_mut_ts, start, end, 0.3)
    :param float proportion: The proportion of mutations to remove.
    '''
    tables = ts.dump_tables()
    tables.sites.clear()
    tables.mutations.clear()
    for tree in ts.trees():
        for site in tree.sites():
            assert len(site.mutations) == 1
            mut = site.mutations[0]
            keep_mutation = True
            for i in range(len(start)):
                left = start[i]
                right = end[i]
                assert(left < right)
                if i > 0:
                    assert(end[i - 1] <= left)
                if left <= site.position < right:
                    keep_mutation = (random.uniform(0, 1) > proportion)
            if keep_mutation:
                site_id = tables.sites.add_row(
                    position=site.position,
                    ancestral_state=site.ancestral_state)
                tables.mutations.add_row(
                    site=site_id, node=mut.node, derived_state=mut.derived_state)
    return tables.tree_sequence()

def noncoding_vcf(ts, vcfoutpath, mu, start, end, archaic_sample_time,proportion=1.):
    ts = pyslim.SlimTreeSequence(msprime.mutate(ts, rate=float(1.8e-8),
         keep=False))
    proportion = 1.
    ts = remove_mutations(ts, start, end, proportion)
    #these_indivs = np.append(ts_sample.individuals_alive_at(archaic_sample_time),ts_sample.individuals_alive_at(0))
    with open(vcfoutpath, "w") as vcf_file:
        ts.write_vcf(vcf_file,ploidy=2)

def syn_vcf(ts, vcfoutpath, mu, start, end,archaic_sample_time):
    mu = 1/3.31*mu
    ts = pyslim.SlimTreeSequence(msprime.mutate(ts, rate=float(mu),
         keep=False))
    proportion = 1.0
    first = ts.first().interval[0]
    last = ts.last().interval[1]
    new_start = [first] + [x for x in end[0:(len(end))]]
    new_end = [x for x in start[0:(len(start))]] + [last]
    ts = remove_mutations(ts, new_start, new_end, proportion)
    #these_indivs = np.append(ts_sample.individuals_alive_at(archaic_sample_time),ts_sample.individuals_alive_at(0))
    with open(vcfoutpath, "w") as vcf_file:
        ts.write_vcf(vcf_file,ploidy=2)


def update_par_file(temp_par, new_par,whichrep,region_name,dominance): #dominance:1 = additive, 2=partial, 3=recessive
    dominance_line=6
    m1_line=11
    segment_line = 22
    #Tadm_start = 103
    #Tadm_end = 109
    outputline=162
    
    oldfile = open(temp_par)
    newfile = open(new_par,'w')
    line_counter=0
    for line_counter, line in enumerate(oldfile):
        fields = line.split()
        #print(fields)
        #print(line_counter)

        if line_counter==dominance_line:
            fields[1] = str(dominance)+');'
        elif line_counter==m1_line:
            if dominance==3:
                fields[1] = str(0.0)+','
            else:
                fields[1] = str(0.5)+','
        elif line_counter==segment_line:
            #print(fields)
            fields[2] = 'readFile("'+DIR_master+'1kseg_byexon/sim_seq_info_'+str(region_name)+'.txt");'
        elif line_counter==outputline:
            #print(fields)
            fields[0] = 'sim.treeSeqOutput("'+DIR_master+''+str(region_name)+'_" ' 
            fields[6] = 'asString('+str(whichrep)+')'
            #print(fields)
        #elif line_counter == Tadm_start:
        #    fields[0] = str(adm_time) 
        new_line=str()
        for item in fields:
            new_line = new_line+item+" "
        newfile.write(new_line+'\n') 
        #print(fields)
    newfile.close()
    oldfile.close()

def ancestry_p_varies(ts,pop,recipientpop,nsize,duration): #pop=source pop
    #nsize = 4108
    n=nsize*2
    mixtime=duration+1
    #times = [x.time for x in ts.nodes() if (x.population == 2)]
    #ids = [x.id for x in ts.nodes() if (x.population == 2)]
    #times = list(set(times))
    #times.sort()
    p = [x.id for x in ts.nodes() if ((x.population == int(pop)) and (x.time == mixtime))] #source pop
    today = [x.id for x in ts.nodes() if ((x.population == int(recipientpop)) and (x.time == 0))] #assuming p3 is recipient
    tree_p = [sum([t.num_tracked_samples(u) for u in p])/n
               for t in ts.trees(tracked_samples=today, sample_counts=True)]      
    starts=[]
    ends=[]
    for x in ts.trees():
        starts.append(x.interval[0])
        ends.append(x.interval[1])
    #sum(tree_p)/len(tree_p)    
    return tree_p


def calc_p1ancestry (treepath, admpop, recipientpop,popsize,t_sinceadm):
    #treepath = "test.trees"
    ts = pyslim.load(treepath)
    any_ancestry = ancestry_p_varies(ts,admpop,recipientpop,popsize,t_sinceadm)
    meanp1 = sum(any_ancestry)/len(any_ancestry)   
    starts=[]
    ends=[]
    for x in ts.trees():
        starts.append(x.interval[0])
        ends.append(x.interval[1])
    return meanp1,any_ancestry,starts,ends


def calc_ancestry_window (ancestry,start_pos,end_pos,len_genome,window_size):
    #ancestry_file = DIR_anc+ region_name + "_"+str(n) + '.ancestry'
        #print(line)
    allpos_bin = np.linspace(0,len_genome,int((len_genome+1)/window_size)+1) #windows of every 50kb
    endpos_digitized = np.digitize(end_pos, allpos_bin)
    end_pos = np.array(end_pos)
    ancestry = np.array(ancestry)    
    anc_window = []
    anc_pos = []    
    for w in range(1,int((len_genome+1)/window_size)+1):
        these_pos = end_pos[endpos_digitized==w]
        these_anc = ancestry[endpos_digitized==w]           
        if(len(these_pos))>0:
            anc_window.append(np.mean(these_anc))
            anc_pos.append(these_pos)
        else:
            anc_window.append(float('nan'))
            anc_pos.append(these_pos)        
    return anc_window    

def ancestry_position_writeout (ancestry,starts,ends,outfilename):
    #record positions
    #write ancestry information to a file    
    outfile = open(outfilename, 'w')
    outfile.write('start,end,ancestry\n')
    for start, end, anc in zip(starts, ends, ancestry):
        outfile.write('{0},{1},{2}\n'.format(start, end, anc))
    outfile.close()



def vcf2genos (vcfpath,p1size,p2size,p3size):   #sizes are int
    #vcfpath = "test_combined.vcf"
    #p1size,p2size,p3size = 2,100,100    
    vcf = open(vcfpath)
    start = 0
    startg = 0
    pos = []
    geno = []
    muttype = []
    for count,line in enumerate(vcf):
        if line[0:6]=="#CHROM":
            start +=1
            startg = count
        if (start == 1) & (startg!=0) & (count>startg):
            fields = line.split()
            if ("1|2" not in fields) and ("0|2" not in fields) and ("2|0" not in fields) and ("2|1" not in fields) and ("2|2" not in fields):
                pos.append(int(fields[1]))
                muttype.append(fields[7][3:])
                geno.append(fields[9:])
    vcf.close()
    len_genome = int(pos[len(pos)-1])-int(pos[0])    
    geno1 = []
    geno2 = []
    geno3 = []
    for n in range(len(geno)):
        geno1.append([geno[n][x] for x in range(p1size)])
        geno2.append([geno[n][x] for x in range(p1size,p1size+p2size)])
        geno3.append([geno[n][x] for x in range(p1size+p2size,p1size+p2size+p3size)])
    return geno1, geno2, geno3, len_genome, muttype, pos


def geno2hap (geno):
    #geno=geno1    
    num_pos = len(geno)
    num_ind = len(geno[1])    
    geno_ref = np.empty([num_pos,num_ind])
    geno_alt = np.empty([num_pos,num_ind])
    for n in range(num_pos):
        geno_ref[n] = [int(i.split('|', 1)[0]) for i in geno[n]]
        geno_alt[n] = [int(i.split('|', 1)[1]) for i in geno[n]]
    geno_ref = np.transpose(geno_ref)
    geno_alt = np.transpose(geno_alt)
    hapmat = np.append(geno_ref,geno_alt,axis=0)   
    return hapmat



#######################################################
#scikit allel part
def haplotype_diversity(hap):
    return allel.haplotype_diversity(np.transpose(hap))

def moving_haplotype_diversity(pos, hap, len_genome,wind_size):
    windowed_hap_diversity = allel.moving_haplotype_diversity(hap, size=wind_size, start=0, stop=len_genome+1)
    return windowed_hap_diversity
  

def heterozygosity_observed(hap):
    g = allel.HaplotypeArray(np.transpose(hap)).to_genotypes(ploidy=2)
    return allel.heterozygosity_observed(g)


def calc_freq (pop_hap):
    popfreq = np.sum(pop_hap, axis=0)
    popfreq = popfreq/ float(pop_hap.shape[0])
    return popfreq

def vSumFunc(archaic_hap,other_hap, currentArchi):
    current_hap = np.array([archaic_hap[currentArchi,]])
    div = np.zeros(other_hap.shape)
    ones = np.ones((other_hap.shape[0],1))
    current_hap = current_hap
    current_hap_extended = np.dot(ones, current_hap)
    div = np.logical_xor(current_hap_extended == 1, other_hap == 1)
    return np.add.reduce(div, 1)   

def my_stats(p1_hap,p2_hap,p3_hap,p1_pos,window_ranges,wind_size): #D,fD, Het
    
    p1_freq = calc_freq (p1_hap)
    p2_freq = calc_freq (p2_hap)
    p3_freq = calc_freq (p3_hap)
    
    Dstat = []
    fD = []
    Het = []
    RD = []
    Q95 = []
    U0 = []
    U20 = []
    U50 = []
    U80 = []
    
    seg_p1 = []
    seg_p2=[]
    seg_p3=[]
    seg_priv_p1=[]
    seg_priv_p2=[]
    seg_priv_p3=[]
    for n in range(0,len(window_ranges)):
        start=window_ranges[n][0]
        end = window_ranges[n][1]
        these = (p1_pos >= start) & (p1_pos <= end)
        p1_freq_these = p1_freq[these]
        p2_freq_these = p2_freq[these]
        p3_freq_these = p3_freq[these]
        
        p1_hap_these=p1_hap[:,these]
        p2_hap_these=p2_hap[:,these]
        p3_hap_these=p3_hap[:,these]

        abbavec = (1.0 - p2_freq_these)*p3_freq_these*p1_freq_these
        babavec = p2_freq_these*(1.0 - p3_freq_these)*p1_freq_these
        abba = np.sum(abbavec)
        baba = np.sum(babavec)

        if (abba + baba > 0):
            Dstat_these = (abba - baba) / (abba + baba)
        else:
            Dstat_these = float('nan')
        Dstat.append(Dstat_these)

        checkfd1 = (p3_freq_these > p1_freq_these)
        abbafd1 = (1.0 - p2_freq_these)*p3_freq_these*p3_freq_these
        babafd1 = p2_freq_these*(1.0 - p3_freq_these)*p3_freq_these
        checkfd2 = (p3_freq_these < p1_freq_these)
        abbafd2 = (1.0 - p2_freq_these)*p1_freq_these*p1_freq_these
        babafd2 = p2_freq_these*(1.0 - p1_freq_these)*p1_freq_these
        abbafd = checkfd1 * abbafd1 + checkfd2 * abbafd2
        babafd = checkfd1 * babafd1 + checkfd2 * babafd2
        abbafd = np.sum(abbafd)
        babafd = np.sum(babafd)
        if (abbafd + babafd > 0):
            fD_these = (abba - baba) / (abbafd - babafd)
        else:
            fD_these = float('nan')
        fD.append(fD_these)
        
        hetvec = 2 * p3_freq_these * (1.0 - p3_freq_these)
        Het_these = np.sum(hetvec) / wind_size
        Het.append(Het_these)

        divratio = []
        for archi in range(0, p1_hap_these.shape[0]): #iterate over 0-99 haps; 100 total
            divarchintro = vSumFunc(p1_hap_these,p3_hap_these, archi)
            divarchintro = divarchintro.astype("float")
            divarchnonintro = vSumFunc(p1_hap_these,p2_hap_these, archi)        
            divarchnonintro = divarchnonintro.astype("float") #took the inversion here so that the multiplying below is really actually dividing; probably not necessary       
            for comb in itertools.product(divarchintro,divarchnonintro): #pairwise combos of divarchintro/divarchnonintro
                if comb[1] != 0:
                    divratio.append(comb[0]/comb[1])
        divratioavg = float(sum(divratio)) / float(len(divratio)) #len(divratio) = 100* (100*100)
        RD.append(divratioavg)

        Arc100 = (p1_freq_these == 1)
        NonAdm1 = (p2_freq_these < 0.01) 
        Arc100NonAdm1 = (Arc100 & NonAdm1)
        Freqs_Arc100NonAdm1 = p3_freq_these[np.where(Arc100NonAdm1 == True)]
        if Freqs_Arc100NonAdm1.size > 0:
            Q_1_100_q95 = np.percentile(Freqs_Arc100NonAdm1,95)
        else:
            Q_1_100_q95 = float('nan')
        U_1_0_100 = ( Arc100NonAdm1 & (p3_freq_these > 0) )
        U_1_20_100 = ( Arc100NonAdm1 & (p3_freq_these > 0.2) )
        U_1_50_100 = ( Arc100NonAdm1 & (p3_freq_these > 0.5) )
        U_1_80_100 = ( Arc100NonAdm1 & (p3_freq_these > 0.8) )        
        U_1_0_100 = np.sum(U_1_0_100)
        U_1_20_100 = np.sum(U_1_20_100)
        U_1_50_100 = np.sum(U_1_50_100)
        U_1_80_100 = np.sum(U_1_80_100)
        Q95.append(Q_1_100_q95)
        U0.append(U_1_0_100)
        U20.append(U_1_20_100)
        U50.append(U_1_50_100)
        U80.append(U_1_80_100)
        
        sites_p1_these =np.sum(p1_hap_these, axis=0)
        sites_p2_these =np.sum(p2_hap_these, axis=0)
        sites_p3_these =np.sum(p3_hap_these, axis=0)
        seg_p1_these = sum(sites_p1_these>0)
        seg_p2_these = sum(sites_p2_these>0)
        seg_p3_these = sum(sites_p3_these>0)
        seg_p1.append(seg_p1_these)
        seg_p2.append(seg_p2_these)
        seg_p3.append(seg_p3_these)
        private_p1 = sum((sites_p1_these>0) & (sites_p2_these==0) & (sites_p3_these==0))
        private_p2 = sum((sites_p2_these>0) & (sites_p1_these==0) & (sites_p3_these==0))
        private_p3 = sum((sites_p3_these>0) & (sites_p2_these==0) & (sites_p1_these==0))
        seg_priv_p1.append(private_p1)
        seg_priv_p2.append(private_p2)
        seg_priv_p3.append(private_p3)
        
    Dstat=np.array(Dstat)
    fD=np.array(fD)
    RD=np.array(RD)
    Het=np.array(Het)
    Q95=np.array(Q95)
    U0=np.array(U0)
    U20=np.array(U20)
    U50=np.array(U50)
    U80=np.array(U80)
    seg_p1=np.array(seg_p1)
    seg_p2=np.array(seg_p2)
    seg_p3=np.array(seg_p3)
    seg_priv_p1=np.array(seg_priv_p1)
    seg_priv_p2=np.array(seg_priv_p2)
    seg_priv_p3=np.array(seg_priv_p3)
    
    return Dstat,fD,RD,Het,Q95,U0,U20,U50,U80,seg_p1,seg_p2,seg_p3,seg_priv_p1,seg_priv_p2,seg_priv_p3


def moving_patterson_f3(hapA, hapB, hapC, len_genome,wind_size):
    aca = allel.HaplotypeArray(np.transpose(hapA)).count_alleles()
    acb = allel.HaplotypeArray(np.transpose(hapB)).count_alleles()
    acc = allel.HaplotypeArray(np.transpose(hapC)).count_alleles()                                                                   
    return allel.moving_patterson_f3(acc, aca, acb, wind_size, start=0, stop=len_genome+1)

def mean_pairwise_difference(hap):
    ac = allel.HaplotypeArray(np.transpose(hap)).count_alleles()                                                                                       
    return allel.mean_pairwise_difference(ac)

def windowed_diversity(pos, hap, len_genome,wind_size):
    g = allel.HaplotypeArray(np.transpose(hap)).to_genotypes(ploidy=2)
    ac = g.count_alleles()
    pi, windows, n_bases, counts = allel.windowed_diversity(pos, ac, size=wind_size, start=0, stop=len_genome+1)
    windowed_diversity_stats = [pi, windows, n_bases, counts]
    return windowed_diversity_stats

def windowed_divergence(pos, hap1, hap2, len_genome,wind_size):
    ac1 = allel.HaplotypeArray(np.transpose(hap1)).count_alleles() 
    ac2 = allel.HaplotypeArray(np.transpose(hap2)).count_alleles()
    dxy, windows, n_bases, counts = allel.windowed_divergence(
        pos, ac1, ac2, size=wind_size, start=0, stop=len_genome+1)
    return [dxy, windows, n_bases, counts]

def windowed_watterson_theta(pos, hap,len_genome,wind_size):
    ac = allel.HaplotypeArray(np.transpose(hap)).count_alleles()                                                                           
    return allel.windowed_watterson_theta(pos, ac, size=wind_size,start=0, stop=len_genome+1)

def windowed_tajima_d(pos, hap, len_genome,wind_size):
    ac = allel.HaplotypeArray(np.transpose(hap)).count_alleles()
    D, windows, counts = allel.windowed_tajima_d(pos, ac, wind_size, start=0, stop=len_genome+1)
    return [D, windows, counts]

def windowed_df(pos, hap1, hap2, len_genome,wind_size):
    ac1 = allel.HaplotypeArray(np.transpose(hap1)).count_alleles()
    ac2 = allel.HaplotypeArray(np.transpose(hap2)).count_alleles()
    df, windows, n_bases, counts = allel.windowed_df(pos, ac1, ac2, size=wind_size, start=0, stop=len_genome+1)
    return [df, windows, n_bases, counts]


def moving_garud_h(hap, len_genome,wind_size):
    h1, h12, h123, h2_h1 = allel.moving_garud_h((np.transpose(hap)), size=wind_size, start=0, stop=len_genome+1)
    return [h1, h12, h123, h2_h1]


def compute_scikit_stats(pos_list,hapmat1,hapmat2,hapmat3,len_genome, wind_size):
    p1_pos=p2_pos=p3_pos=pos_list    
    p1_hap=hapmat1.astype('int') #archiac
    p2_hap=hapmat2.astype('int') #african
    p3_hap=hapmat3.astype('int') #european   
    #wind_size=window_size
#####diversity related 
    #calculate allel.windowed_diversity for 3 pops separately 
    window_ranges = windowed_diversity(p3_pos, p3_hap, len_genome,wind_size)[1]    
    #p3windowed_diversities = windowed_diversity(p3_pos, p3_hap, len_genome,wind_size)[0]
    #print(len(p3windowed_diversities))
    p3windowed_variant_count = windowed_diversity(p3_pos, p3_hap, len_genome,wind_size)[3]
    #print(len(p3windowed_variant_count))
    #calculate windowed divergence
    #pop1 = admixed; pop2 = source
    divergence_p3_p1 = windowed_divergence(p1_pos, p3_hap, p1_hap, len_genome,wind_size)[0]
    #print(len(divergence_p3_p1))
    #pop1 = admixed; pop2 = outgroup
    divergence_p3_p2 = windowed_divergence(p1_pos, p3_hap, p2_hap, len_genome,wind_size)[0]
    #print(len(divergence_p3_p2))
    watterson_theta_p3 = windowed_watterson_theta(p3_pos, p3_hap,len_genome,wind_size)[0]
    #print(len(watterson_theta_p3))
    #calculate windowed df
    #pop1 = admixed; pop2 = source
    df_p3_p1 = windowed_df(p3_pos, p3_hap, p1_hap, len_genome,wind_size)[0]
    #print(len(df_p3_p1))
    #pop1 = admixed; pop2 = outgroup
    df_p3_p2 = windowed_df(p3_pos, p3_hap, p2_hap, len_genome,wind_size)[0]
    #print(len(df_p3_p2))
#####selection related    
    hap_diversity_p3 = moving_haplotype_diversity(p3_pos, p3_hap, len_genome,wind_size)
    #print(len(hap_diversity_p3))
    #calculate windowed tajima d for 3 pops separately
    windowed_tajima_d_p3 = windowed_tajima_d(p3_pos, p3_hap, len_genome,wind_size)[0]
    #print(len(windowed_tajima_d_p3))
    #calculate moving_garud_h in admixed pop
    garud_h1 = moving_garud_h(p3_hap, len_genome,wind_size)[0]
    garud_h12 = moving_garud_h(p3_hap, len_genome,wind_size)[1]
    garud_h2_h1 = moving_garud_h(p3_hap, len_genome,wind_size)[2]
    #print(len(garud_h1))    
    patterson_f3 = moving_patterson_f3(p1_hap,p2_hap,p3_hap,len_genome,wind_size)
    #print(len(patterson_f3))    
    Dstat,fD,RD,Het,Q95,U0,U20,U50,U80,seg_p1,seg_p2,seg_p3,seg_priv_p1,seg_priv_p2,seg_priv_p3 = my_stats(p1_hap,p2_hap,p3_hap,p1_pos,window_ranges,wind_size)    
    #patterson_d = moving_patterson_d(p2_hap,p3_hap, p1_hap, wind_size) #(A,B; C,D), D = all 0s
    #calculate number of private segregating sites for each population  
    #stats = [haplotype_diversities, het_observed, all_sites, mean_pairwise_differences, windowed_diversities, divergence_adm_source, divergence_adm_out, windowed_tajima_d_results, df_adm_source, df_adm_out, garud, private_sites]        
    return window_ranges,p3windowed_variant_count,divergence_p3_p1,divergence_p3_p2,watterson_theta_p3,df_p3_p1,df_p3_p2,hap_diversity_p3,windowed_tajima_d_p3,garud_h1,garud_h12,garud_h2_h1,patterson_f3,Dstat,fD,RD,Het,Q95,U0,U20,U50,U80,seg_p1,seg_p2,seg_p3,seg_priv_p1,seg_priv_p2,seg_priv_p3








#######################################################

if __name__ == "__main__":
    DIR_master = "/Users/xinjunzhang/Desktop/Dominance_project/simulation/"
    DIR_master = "/u/scratch/x/xinjunzh/dominance/"
    Path_slim = "/usr/local/bin/slim"
    Path_slim="/u/home/x/xinjunzh/slim_build/slim"

    os.chdir(DIR_master)
    

    segs_path = DIR_master+"1kseg_byexon/"
    allfiles = [f for f in os.listdir(segs_path) if f.endswith("txt")]
    genes = [f.partition("info_")[2].partition(".txt")[0] for f in allfiles] #len(genes) = 26
    thisgene = genes[whichgene-1]
    thisfile = allfiles[whichgene-1]

    
    num_rep = 100#0 #make it 10
    #num_rep = 100
    #num_rep = 1
    dominance_list = ["additive","partial","recessive"]
    dominance = dominance_list[whichh-1]
    
    slim_template = "Gravel_Eurasian_varyh.slim"
    
    
    rep=0
    while rep < num_rep:
        print(rep)

#1. run slim, get tree file
        #update slim file
        #adm_time = random.choice(range(8697,8747)) 
        adm_time=8734 #generation of admixture on slim script
        
        slim_rep = "rep"+str(rep)+"_model"+str(whichh)+thisgene+slim_template
        slim_rep_out = "rep"+str(rep)+"_model"+str(whichh)+thisgene+slim_template+".output"
        update_par_file(slim_template, slim_rep, rep, thisgene,whichh)
        #os.system('/u/home/x/xinjunzh/slim_build/slim %s' %(slim_rep))
        os.system('%s %s > %s' %(Path_slim,slim_rep,slim_rep_out))
        print("slim done")

#2. from tree to vcf
        filename = thisgene+"_"+str(whichh)+"_rep"+str(rep)
        #filename = "slim3.2"
        scalingfactor = 10
        mu = 1.5e-8 * scalingfactor
        archaic_sample_time = 152.0

#number of individuals from [afr, eur,asn]
        sampsize=[30,30,30] #by default Neanderthal is 2, so not specified here #order of p1, p3, p4 = yri, ceu, chb

#get exon definitions
        exonfilepath = segs_path+thisfile #'sim_seq_info_EPAS1.txt'
        lines = open(exonfilepath, 'r').readlines()
        rec_lines = [x for x in lines if x.startswith('rec')]
        rec_lines = [x.strip('\n').split(' ') for x in rec_lines]
        lines = [x for x in lines if x.startswith('exon')]
        lines = [x.strip('\n').split(' ') for x in lines]
        
        rec_annot,rec_end,rec_rate = zip(*rec_lines)
        rec_end = [int(x) for x in rec_end]
        rec_rate = [float(x) for x in rec_rate]
        intervals = np.diff(np.array(rec_end))
        intervals = np.insert(intervals,0,rec_end[0])
        
        annot, start, end = zip(*lines)
        start = [int(x) for x in start]
        end = [int(x) for x in end]

        ts_sample= simplify_tree('{0}.trees'.format(filename), sampsize, archaic_sample_time)
        
        selcoeff, mut_pos = get_selcoeff(ts_sample)
        selcoeff = [x / scalingfactor for x in selcoeff]
        ns_info = [mut_pos,selcoeff] 
        with open("{0}_ns_selcoeff.out".format(filename), "w") as file:
            for x in zip(*ns_info):
                file.write("{0}\t{1}\n".format(*x))

        ns_vcf(ts_sample,archaic_sample_time,'{0}_ns.vcf'.format(filename))


        noncoding_vcf(ts_sample, '{0}_noncoding.vcf'.format(filename), mu, start, end,archaic_sample_time)
        syn_vcf(ts_sample, '{0}_syn.vcf'.format(filename), mu, start, end,archaic_sample_time)

        cmd = '''
        cat {0}_ns.vcf | grep "##" > {0}_unsorted.vcf;
        echo '##INFO=<ID=TT,Number=A,Type=String,Description="Annotation">' >> {0}_unsorted.vcf;
        cat {0}_ns.vcf | grep "#CHROM" >> {0}_unsorted.vcf;
        cat {0}_ns.vcf | grep -v "#" | awk '$8="TT=NS"' | tr ' ' '\t' >> {0}_unsorted.vcf;
        cat {0}_syn.vcf | grep -v "#" | awk '$8="TT=SYN"' | tr ' ' '\t' >> {0}_unsorted.vcf;
        cat {0}_noncoding.vcf | grep -v "#" | awk '$8="TT=NC"' | tr ' ' '\t' >> {0}_unsorted.vcf;
        bcftools sort {0}_unsorted.vcf -o {0}_combined.vcf
        '''.format(filename)
        
        os.system(cmd)
        print("vcf done")

        vcfpath = filename + "_combined.vcf"
        #os.system("cp "+vcfpath+" vcf/") #store vcf copy
        
#get ancestry per window        
        treepath = filename+".trees"
        #os.system("cp "+treepath+" vcf/")
        
        
        len_genome = 4999999
        #len_genome = 19999999 #in the test run, 20MB
        window_size=1000000#50000
        p1size,p2size,p3size = 2,sampsize[0],sampsize[1]
        
        p3_endsize = 4081
        #p4_endsize = 4108
        meanp1inp3,anc_p3,starts,ends = calc_p1ancestry(treepath,2,3,p3_endsize,8900-(adm_time+1)) #8900-8716 = duration, since admixture happened at the end of 8715
        #meanp1inp4,anc_p4,starts,ends = calc_p1ancestry(treepath,2,4,p4_endsize,8900-(adm_time+1)) #8900-8716 = duration, since admixture happened at the end of 8715

        p3anc_window = calc_ancestry_window(anc_p3,starts,ends,len_genome,window_size) #ancestry by window
        #meanp1inp3=np.nanmean(p3anc_window)        
        #p4anc_window = calc_ancestry_window(anc_p4,starts,ends,len_genome,window_size) #ancestry by window
  
#get stats per window    
#for now ignoring p4; add back later if needed
        geno1, geno2,geno3,len_genome_vcf,mut_types, pos_list = vcf2genos(vcfpath,p1size,p2size,p3size)    
        hapmat1 = geno2hap(geno1)
        hapmat2 = geno2hap(geno2)
        hapmat3 = geno2hap(geno3)   
        
        window_ranges,p3windowed_variant_count,divergence_p3_p1,divergence_p3_p2,watterson_theta_p3,df_p3_p1,df_p3_p2,hap_diversity_p3,windowed_tajima_d_p3,garud_h1,garud_h12,garud_h2_h1,patterson_f3,Dstat,fD,RD,Het,Q95,U0,U20,U50,U80,seg_p1,seg_p2,seg_p3,seg_priv_p1,seg_priv_p2,seg_priv_p3 = compute_scikit_stats(pos_list,hapmat1,hapmat2,hapmat3,len_genome, window_size)
        #meanp1inp3
        #p3anc_window
        window_starts = window_ranges[:,0]
        window_ends = window_ranges[:,1]
        
        exon_windows = []
        rec_windows = []
        for i in range(0,len(window_ranges)):
            window_start = window_ranges[i][0]
            window_end = window_ranges[i][1]
            exon_count = sum((end<=window_end) & (end>=window_start))
            exon_windows.append(exon_count)
            
            rec_these = (rec_end<=window_end) & (rec_end>=window_start)
            intervals_these = intervals[rec_these] 
            rate_these = np.array(rec_rate)[rec_these]
            rate = sum(intervals_these*rate_these)/(window_end-window_start)
            rec_windows.append(rate)
        rec_windows = np.array(rec_windows)
        exon_windows = np.array(exon_windows)
        
        p3anc_window = np.array(p3anc_window)
            
        exon_total = len(end)
        mean_rec = sum(intervals*rec_rate)/(len_genome+1)
        
        #make constants arrays
        exon_total=np.repeat(exon_total,len(window_ranges))
        mean_rec = np.repeat(mean_rec,len(window_ranges))
        meanp1inp3 = np.repeat(meanp1inp3,len(window_ranges))
        rep_id = np.repeat(rep,len(window_ranges))
        seg_id = np.repeat(thisgene,len(window_ranges))
        if whichh==1:
            dom = 0
        elif whichh==3:
            dom=1
        elif whichh==2:
            dom=2
        dom_class = np.repeat(dom,len(window_ranges))
        
        stats_tuple = (seg_id,rep_id,window_starts,window_ends,dom_class,exon_total,mean_rec,meanp1inp3,exon_windows,rec_windows,
                             p3anc_window,p3windowed_variant_count,divergence_p3_p1,divergence_p3_p2,watterson_theta_p3,df_p3_p1,
                             df_p3_p2,hap_diversity_p3,windowed_tajima_d_p3,garud_h1,garud_h12,garud_h2_h1,patterson_f3,Dstat,fD,
                             RD,Het,Q95,U0,U20,U50,U80,seg_p1,seg_p2,seg_p3,seg_priv_p1,seg_priv_p2,seg_priv_p3)
        stats_all = np.vstack(stats_tuple)
        
        stats_table = pd.DataFrame(np.transpose(stats_all))
        
        stats_table.columns = ["segment","rep","start","end","dominance","exon_density","mean_recrate","mean_introg_anc","exon_window","recrate_window",
                               "introg_anc_window","num_variant_window","divergence_p3_p1","divergence_p3_p2","watterson_theta_p3","df_p3_p1",
                               "df_p3_p2","hap_diversity_p3","windowed_tajima_d_p3","garud_h1","garud_h12","garud_h2_h1","patterson_f3","D","fD",
                               "RD","Het","Q95","U0","U20","U50","U80","num_seg_p1","num_seg_p2","num_seg_p3","num_private_seg_p1","num_private_seg_p2","num_private_seg_p3"]

        stats_table.to_csv("stat_1kseg_byexon/"+filename+"_1MB_stat.csv", index=False)
        #make first columns: segment, rep, start, end, dominance, exon, rec, meanp1anc
        #windowed columns: exon_window, rec_window,anc_window, all other stats

#start another replicate        
        rep+=1
         #remove all temp files
        os.remove('{0}_ns.vcf'.format(filename))
        os.remove('{0}_noncoding.vcf'.format(filename))
        os.remove('{0}_syn.vcf'.format(filename))
        os.remove('{0}_unsorted.vcf'.format(filename))
        os.remove(slim_rep)
        os.remove(slim_rep_out)
        os.remove(vcfpath)
        os.remove('{0}.trees'.format(filename))
        os.remove('{0}_ns_selcoeff.out'.format(filename))
        

       

    print("END OF SIMULATIONS")














