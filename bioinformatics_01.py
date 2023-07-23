# Import libraries

from Bio.Seq import Seq
import streamlit as st
from Bio import SeqIO
from io import StringIO
from Bio import pairwise2
import matplotlib.pyplot as plt

# Header Main Page

st.title('Bioinformatics Server v:1.00')
st.write("""
#### Welcome to the Bioinformatics Server
---        
""")

# Sidebar Page

st.sidebar.header("Available Operations🧬💻")
sidebar_options = st.sidebar.radio("Select Operations",['Sequence Operation','FASTA File Analyser','GENBANK File Analyser','Pairwise Alignment','Nucleotide Plots','About Us'])

# 1 - Basic sequence operations 

if sidebar_options == 'Sequence Operation':
    st.warning("Default value set below.Change as per your work! 😉")
    dna_query = st.text_area('Enter Query DNA Sequence',value='ATGC')
    dna_query = Seq(dna_query)
    transcribe_dna = dna_query.transcribe()
    translate_dna = dna_query.translate()
    nucleotide_base_counter = len(dna_query)
    content_GC = (dna_query.count('G') + dna_query.count('C'))/(dna_query.count('A')+dna_query.count('T')+dna_query.count('C')+dna_query.count('G'))

    if st.button('Get Results'):
        st.write(f"Transcribed DNA  : {transcribe_dna}")
        st.write(f"Translated DNA : {translate_dna}")
        st.write(f"Length of DNA : {nucleotide_base_counter} bp")
        st.write(f"GC Content of sequence : {content_GC*100} %")

# 2 - FASTA File Analyser

if sidebar_options == 'FASTA File Analyser':

    # We have to specifically decode these files to allow streamlit work

    def read_fasta_file(uploaded_file):
        fasta_content = uploaded_file.read().decode()
        return list(SeqIO.parse(StringIO(fasta_content), "fasta"))

    # Create function to display the biological information
    
    def display_sequences(sequence_records):
        if len(sequence_records) == 0:
            st.write("Sorry ! No Sequences Found.")
        else:
            for i, record in enumerate(sequence_records, start=1):
                st.write(f"Sequence {i} - ID: {record.id}")
                st.write(record.seq)
                st.write("Length:", len(record.seq))
                st.write("Description:", record.description)
                st.write("---")
    st.header("FASTA File Viewer")
    st.write("Upload a FASTA file and view its sequences:")

    uploaded_file = st.file_uploader("Upload your FASTA File", type=["fasta"])

    if uploaded_file is not None:
        sequence_records = read_fasta_file(uploaded_file)
        display_sequences(sequence_records)
        st.success("You can copy the sequence and perform sequence operations in the previous page.",icon="😍")


# 3 - GENBANK File Analyser

if sidebar_options == 'GENBANK File Analyser':
    # We have to specifically decode these files to allow streamlit work

    def read_genbank_file(uploaded_file_gb):
        genbank_content = uploaded_file_gb.read().decode()
        return list(SeqIO.parse(StringIO(genbank_content), "genbank"))
    
    def display_sequences_gb(sequence_records_gb):
        if len(sequence_records_gb) == 0:
            st.write("Sorry ! No Sequences Found.")
        else:
            for i, record in enumerate(sequence_records_gb, start=1):
                st.write(f"Sequence {i} - ID: {record.id}")
                st.write(record.seq)
                st.write("Length:", len(record.seq))
                st.write("Description:", record.description)
                st.write("---")
    st.header("GENBANK File Viewer")
    st.write("Upload a GENBANK file and view its sequences:")

    uploaded_file_gb = st.file_uploader("Upload your GENBANK File", type=["genbank"])

    if uploaded_file_gb is not None:
        sequence_records_gb = read_genbank_file(uploaded_file_gb)
        display_sequences_gb(sequence_records_gb)
        st.success("You can copy the sequence and perform sequence operations in the 1st page.",icon="😍")

# 4 - Pairwise Alignment

if sidebar_options == 'Pairwise Alignment':
    st.header("Pairwise Alignment Operations")
    st.write("Upload FASTA files to visualise pairwise alignment.")
    # We have to specifically decode these files to allow streamlit work

    # Dedicated function to access query fasta file

    def read_fasta_file_pairwise_query(uploaded_file_pairwise_query):
        fasta_content_pairwise_query = uploaded_file_pairwise_query.read().decode()
        return list(SeqIO.parse(StringIO(fasta_content_pairwise_query), "fasta"))
    
    # Dedicated function to access target fasta file

    def read_fasta_file_pairwise_target(uploaded_file_pairwise_target):
        fasta_content_pairwise_target = uploaded_file_pairwise_target.read().decode()
        return list(SeqIO.parse(StringIO(fasta_content_pairwise_target), "fasta"))
    
    # Dedicated function to operate query fasta file

    def display_sequences_pairwise_query(sequence_records_pair_query):
        if len(sequence_records_pair_query) == 0:
            st.write("Sorry ! No Sequences Found.")
        else:
            for i, record_query in enumerate(sequence_records_pair_query, start=1):
                st.write(f"Sequence Query : - ID: {record_query.id}")
                st.write(record_query.seq)
    
                
    # Dedicated function to operate target fasta file

    def display_sequences_pairwise_target(sequence_records_pair_target):
        if len(sequence_records_pair_target) == 0:
            st.write("Sorry ! No Sequences Found.")
        else:
            for i, record_target in enumerate(sequence_records_pair_target, start=1):
                st.write(f"Sequence Target - ID: {record_target.id}")
                st.write(record_target.seq)


    uploaded_file_pairwise_query = st.file_uploader("Upload your Query FASTA File", type=["fasta"])
    uploaded_file_pairwise_target = st.file_uploader("Upload your Target FASTA File", type=["fasta"])

    if uploaded_file_pairwise_query is not None:
        sequence_records_pair_query = read_fasta_file_pairwise_query(uploaded_file_pairwise_query)
        display_sequences_pairwise_query(sequence_records_pair_query)
        

    if uploaded_file_pairwise_target is not None:
        sequence_records_pair_target = read_fasta_file_pairwise_query(uploaded_file_pairwise_target)
        display_sequences_pairwise_target(sequence_records_pair_target)
    
    query_input = st.text_input("Enter Query",value='ATGC')
    target_input = st.text_input("Enter Target",value='AGCA')
    alignments = pairwise2.align.globalxx(query_input,target_input)
    
    align = pairwise2.format_alignment(*alignments[0]).splitlines()
    for line in align:
        st.write(line)

# 5 - Creating Dot Plots 

if sidebar_options == "Nucleotide Plots":

    def calculate_similarity(seq1, seq2):
        # Check the length of the sequences
        if len(seq1) != len(seq2):
            st.warning("Oh no!! please enter sequences of same length!")

        # Calculate nucleotide similarity
        similarity = [(i, j) for i, a in enumerate(seq1) for j, b in enumerate(seq2) if a == b]
        return similarity

    def plot_alignment(seq1, seq2):

        similarity = calculate_similarity(seq1, seq2)

        # Create a dot plot for alignment
        x, y = zip(*similarity)
        plt.scatter(x, y, color='blue')
        plt.xticks(range(len(seq1)), seq1)
        plt.yticks(range(len(seq2)), seq2)
        plt.xlabel('Sequence 1')
        plt.ylabel('Sequence 2')
        plt.title('Dot Plot for Sequence Alignment')
        plt.grid(True, linestyle='--', alpha=0.5)
        plt.tight_layout()
        st.pyplot(plt)


    st.title("Dot Plot for Sequence Alignment")

    # Text input for the first nucleotide sequence
    seq1 = st.text_input("Enter Sequence 1", value='ATGC')

    # Text input for the second nucleotide sequence
    seq2 = st.text_input("Enter Sequence 2", value='AGCA')

    if st.button("Plot Alignment"):
        try:
            plot_alignment(seq1, seq2)
        except ValueError as e:
            st.error(str(e))

# 6 - About Me
st.sidebar.header('Created with ❤️ by Rajarshi Ray')
st.sidebar.image('https://miro.medium.com/max/1200/1*7Upq9bHVVZ8f6IwGnHvxsA.jpeg')

if sidebar_options == 'About Us':
    st.header("""
        Hello everyone my name is Rajarshi Ray. I am a student of UEM Kolkata from Biotechnology.This is my first ever bioinformatics application build using python.
""")