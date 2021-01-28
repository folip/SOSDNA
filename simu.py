import streamlit as st
import plotly.io as pio
import matplotlib.pyplot as plt
from Model.Model import * 
import plotly.express as px
from Analysis.Analysis import dna_chunk, plot_oligo_number_distribution, plot_error_distribution
pio.templates.default = "plotly_white"
import Model.config as config

def inspect(dnas, num_th = 1, inspect_index = 20):
    fig = plt.figure(figsize = (12,6))
    plt.subplot(1,2,1)
    plot_oligo_number_distribution(dnas)
    plt.subplot(1,2,2)
    plot_error_distribution(dnas,th = num_th)
    st.write(fig)
    dc = dna_chunk(dnas[inspect_index],'html')
    table = dc.plot_re_dnas()
    st.markdown(table, unsafe_allow_html = True)
    st.write(dc.plot_voting_result())

# assign parameters
arg = config.DEFAULT_PASSER

arg.syn_number = st.sidebar.slider('Syn number', min_value = 10, max_value = 50, value = 30)
arg.syn_sub_prob = st.sidebar.number_input('Syn Error rate', min_value = 0.0, max_value = 0.1, value = 0.001) / 3 # 3 kinds of substitutions
arg.syn_yield = st.sidebar.slider('Syn Yield', min_value = 0.98, max_value = 0.995, value = 0.99)

arg.pcrc = st.sidebar.slider('PCR cycle',min_value = 0, max_value =20,value = 12)
arg.pcrp = st.sidebar.number_input('PCR prob',min_value = 0.5, max_value = 1.0,value = 0.8)

arg.sam_ratio = st.sidebar.number_input('Sampling ratio',min_value = 0.0, max_value =1.0,value = 0.005)
arg.seq_depth = st.sidebar.slider('Seq Depth', min_value = 1, max_value = 100, value = 10)
seq_platform = st.sidebar.selectbox('Sequencing Platform',['NGS','Nanopore'])
if seq_platform == 'NGS':
    arg.seq_TM = config.TM_NGS
else:
    arg.seq_TM = config.TM_NNP
st.header('Simulation model for DNA data storage channel')
# st.image('model.png')
# load data
st.subheader('Load Data')
with open('files/out.dna') as f:
    dnas = f.readlines()
in_dnas = [dna.split('\n')[0] for dna in dnas]

'dnas loaded. ', len(in_dnas), ' strands of length ', len(in_dnas[0])
index = st.sidebar.slider('inspect index', max_value = len(in_dnas)-1, value = 0)

st.subheader('Synthesis')
SYN = Synthesizer(arg)
dnas_syn = SYN(in_dnas)
inspect(dnas_syn,inspect_index = index)

st.subheader('Decay')
DEC = Decayer(arg)
dnas_dec = DEC(dnas_syn)
inspect(dnas_dec,inspect_index = index)

st.subheader('PCR')
PCR = PCRer(arg = arg)
dnas_pcr = PCR(dnas_dec)
inspect(dnas_pcr,inspect_index = index)

st.subheader('Sampling')
SAM = Sampler(arg = arg)
dnas_sam = SAM(dnas_pcr)
inspect(dnas_sam,inspect_index = index)

st.subheader('Sequencing')
SEQ = Sequencer(arg)
dnas_seq = SEQ(dnas_sam)
inspect(dnas_seq,inspect_index = index)
