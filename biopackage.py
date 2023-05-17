from Bio.Seq import Seq
import random

def reverse_complement(input):
    output = ''
    for letter in input:
        if letter == 'A':
            output += 'T'
        elif letter == 'T':
            output += 'A'
        elif letter == 'G':
            output += 'C'
        else:
            output += 'G'

    return(output[::-1])

def get_base_counts(dna):
    counts = {'A': 0, 'T': 0, 'G': 0, 'C': 0}
    for base in dna:
        counts[base] += 1
    return counts
'''def get_base_frequencies_v1(dna):
    counts = get_base_counts(dna)
    return {base: count*1.0/len(dna) for base, count in counts.items()}
'''
def get_base_frequencies_v2(dna):
    return {base: dna.count(base)/float(len(dna)) for base in 'ATGC'}
def format_frequencies(frequencies):
    return ', '.join(['%s: %.2f' % (base, frequencies[base]) for base in frequencies])
    print ("Base frequencies of sequence",dna ,"is" ,format_frequencies(frequencies))
def count_v11(dna, base):
    return len([i for i in range(len(dna)) if dna[i] == base])
def freq_lists(dna_list):
    n = len(dna_list[0])
    A = [0]*n
    T = [0]*n
    G = [0]*n
    C = [0]*n
    for dna in dna_list:
        for index, base in enumerate(dna):
            if base == 'A':
                A[index] += 1
            elif base == 'C':
                C[index] += 1
            elif base == 'G':
                G[index] += 1
            elif base == 'T':
                T[index] += 1
    return A, C, G, T          
def mutate_v1(dna):
        dna_list = list(dna)
        mutation_site = random.randint(0, len(dna_list) - 1)
        dna_list[mutation_site] = random.choice(list('ATCG'))
        return ''.join(dna_list)
   



dna = '''ATGTACTCATTCGTTTCGGAAGAGACAGGTACGTTAATAGTTAATAGCGTACTTCTTTTTCTTGCTTTCG
TGGTATTCTTGCTAGTTACACTAGCCATCCTTACTGCGCTTCGATTGTGTGCGTACTGCTGCAATATTGT
TAACGTGAGTCTTGTAAAACCTTCTTTTTACGTTTACTCTCGTGTTAAAAATCTGAATTCTTCTAGAGTT
CCTGATCTTCTGGTCTAA'''

dna_new = dna.replace("\n", "")

rna_map = dna_new.maketrans("ATGC", "UACG")
#protein = ""

rna = dna_new.translate(rna_map)

RNA_Codons = {
    # 'M' - START, '*' - STOP
    "GCU": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "UGU": "C", "UGC": "C",
    "GAU": "D", "GAC": "D",
    "GAA": "E", "GAG": "E",
    "UUU": "F", "UUC": "F",
    "GGU": "G", "GGC": "G", "GGA": "G", "GGG": "G",
    "CAU": "H", "CAC": "H",
    "AUA": "I", "AUU": "I", "AUC": "I",
    "AAA": "K", "AAG": "K",
    "UUA": "L", "UUG": "L", "CUU": "L", "CUC": "L", "CUA": "L", "CUG": "L",
    "AUG": "M",
    "AAU": "N", "AAC": "N",
    "CCU": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAA": "Q", "CAG": "Q",
    "CGU": "R", "CGC": "R", "CGA": "R", "CGG": "R", "AGA": "R", "AGG": "R",
    "UCU": "S", "UCC": "S", "UCA": "S", "UCG": "S", "AGU": "S", "AGC": "S",
    "ACU": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GUU": "V", "GUC": "V", "GUA": "V", "GUG": "V",
    "UGG": "W",
    "UAU": "Y", "UAC": "Y",
    "UAA": "*", "UAG": "*", "UGA": "*"
}

def butto1():
    print("---Representing the place where the Bases occured in Gene---\n")
    count=0
    dna = input("Enter DNA : ")
    base = input("Enter base : ")
    for ele in dna:
        if('A' in ele or 'T' in ele or 'G' in ele or 'C'in ele):
            count+=1
    if(count==len(dna)):
        indices = [i for i in range(len(dna)) if dna[i] == base]
        print("The position where the base occurs is:-- \n",base,indices)
        print("The total number of times it occurs is:--",count_v11(dna,base),"\n\n")
    else:
        print("Enter the correct input!!")

def butto2(): 
    flag=0
    print("---Frequency of the Bases occured in Gene---\n")
    dna_list=[]
    num=int(input("Enter number of dna samples :"))
    print("dna strings you are going to enter is to be of equal length\n")
    for i in range(num):
        b=input("dna string==> ")
        count=0
        for a in b:
            if('A' in a or 'T' in a or 'G' in a or 'C'in a):
                count+=1
        if(count==len(b)):
            dna_list.append(b)
        else:
            flag=1
            break
        #dna_list = ['GGTAG', 'GGTAC', 'GGTGC']
        
    if(flag==1):
        print("wrong input given\n")
    else:
        A, C, G, T = freq_lists(dna_list)
        print("The frequency matrix is:--\n")
        print (A,"\n")
        print (C,"\n")
        print (G,"\n")
        print (T,"\n")

def butto3():
    print("----Random generation of genes-----\n")
    dna = input("Enter the DNA sequence:--")
    frequencies = get_base_frequencies_v2(dna)
    print ("Base frequencies of sequence '%s':\n%s" % (dna, format_frequencies(frequencies)))
    print ('Starting DNA:', dna)
    print (format_frequencies(get_base_frequencies_v2(dna)))

    mutations = int(input("Enter the number of mutations:-- "))
    for i in range(mutations):
        dna = mutate_v1(dna)

    print ('DNA after %d mutations:' % mutations, dna)


def butto4():
    print("---DNA sequence to RNA sequence of homo sapiens---\n")
    DnaSeq ='''ACACTCGCTTCTGGAACGTCTGAGGTTATCAATAAGCTCCTAGTCCAGACGCCATGGGTCATTTCACAGA
            GGAGGACAAGGCTACTATCACAAGCCTGTGGGGCAAGGTGAATGTGGAAGATGCTGGAGGAGAAACCCTG
            GGAAGGTAGGCTCTGGTGACCAGGACAAGGGAGGGAAGGAAGGACCCTGTGCCTGGCAAAAGTCCAGGTC
            GCTTCTCAGGATTTGTGGCACCTTCTGACTGTCAAACTGTTCTTGTCAATCTCACAGGCTCCTGGTTGTC
            TACCCATGGACCCAGAGGTTCTTTGACAGCTTTGGCAACCTGTCCTCTGCCTCTGCCATCATGGGCAACC
            CCAAAGTCAAGGCACATGGCAAGAAGGTGCTGACTTCCTTGGGAGATGCCATAAAGCACCTGGATGATCT
            CAAGGGCACCTTTGCCCAGCTGAGTGAACTGCACTGTGACAAGCTGCATGTGGATCCTGAGAACTTCAAG
            GTGAGTCCAGGAGATGTTTCAGCACTGTTGCCTTTAGTCTCGAGGCAACTTAGACAACTGAGTATTGATC
            TGAGCACAGCAGGGTGTGAGCTGTTTGAAGATACTGGGGTTGGGAGTGAAGAAACTGCAGAGGACTAACT
            GGGCTGAGACCCAGTGGCAATGTTTTAGGGCCTAAGGAGTGCCTCTGAAAATCTAGATGGACAACTTTGA
            CTTTGAGAAAAGAGAGGTGGAAATGAGGAAAATGACTTTTCTTTATTAGATTTCGGTAGAAAGAACTTTC
            ACCTTTCCCCTATTTTTGTTATTCGTTTTAAAACATCTATCTGGAGGCAGGACAAGTATGGTCATTAAAA
            AGATGCAGGCAGAAGGCATATATTGGCTCAGTCAAAGTGGGGAACTTTGGTGGCCAAACATACATTGCTA
            AGGCTATTCCTATATCAGCTGGACACATATAAAATGCTGCTAATGCTTCATTACAAACTTATATCCTTTA
            ATTCCAGATGGGGGCAAAGTATGTCCAGGGGTGAGGAACAATTGAAACATTTGGGCTGGAGTAGATTTTG
            AAAGTCAGCTCTGTGTGTGTGTGTGTGTGTGTGCGCGCGTGTGTTTGTGTGTGTGTGAGAGCGTGTGTTT
            CTTTTAACGTTTTCAGCCTACAGCATACAGGGTTCATGGTGGCAAGAAGATAACAAGATTTAAATTATGG
            CCAGTGACTAGTGCTGCAAGAAGAACAACTACCTGCATTTAATGGGAAAGCAAAATCTCAGGCTTTGAGG
            GAAGTTAACATAGGCTTGATTCTGGGTGGAAGCTTGGTGTGTAGTTATCTGGAGGCCAGGCTGGAGCTCT
            CAGCTCACTATGGGTTCATCTTTATTGTCTCCTTTCATCTCAACAGCTCCTGGGAAATGTGCTGGTGACC
            GTTTTGGCAATCCATTTCGGCAAAGAATTCACCCCTGAGGTGCAGGCTTCCTGGCAGAAGATGGTGACTG
            GAGTGGCCAGTGCCCTGTCCTCCAGATACCACTGAGCTCACTGCCCATGATGCAGAGCTTTCAAGGATAG
            GCTTTATTCTGCAAGCAATCAAATAATAAATCTATTCTGCTAAGAGATCACACA'''
    script = DnaSeq.maketrans("ATGC","UACG")
    a=DnaSeq.translate(script)
    print("RNA sequence is ==>\n")
    print(a,"\n\n")    

def butto5():
    print("----RNA sequence to Protein sequence of homo sapiens-----\n")
    protein=""
    print(f"\nThe DNA sequence is: \n\n{dna_new}")
    print(f"\nThe RNA sequence is: \n\n{rna}")
    for i in range(0, len(rna), 3):
        codon = rna[i:i+3]
        protein += RNA_Codons[codon]
    print(f"\nThe protein sequence is: \n\n{protein}")

    print(f"\nThe RNA sequence is: \n\n{rna}")
    
def butto6():
    print("----To find the DNA complement-----\n")
    my_dna = input("Enter the DNA sequence:--")
    print("Original data: ",my_dna)
        #print(my_dna.complement())
    print("Reverse and Complemented data: ",reverse_complement(my_dna))

def butto0():
    print("\nApplication terminated...\n\n")
    
    exit(0)


from tkinter import *

root = Tk()
root.geometry("500x500")
root.title("Bioinformatics")

#bg=PhotoImage(file="image.png")
#label=Label(root,image=bg)
#label.place(x=0,y=0,relwidth=1,relheight=1)
label1=Label(root,text="Representing the place where the Bases occured in Gene",fg="blue")
label1.grid(row=1,column=6)

label2=Label(root,text="Frequency of the Bases occured in Gene",fg="blue")
label2.grid(row=2,column=6)

label3=Label(root,text="Random generation of genes",fg="blue")
label3.grid(row=3,column=6)

label4=Label(root,text="DNA sequence to RNA sequence of homo sapies",fg="blue")
label4.grid(row=4,column=6)

label5=Label(root,text="RNA sequence to Protein sequence of homo sapies",fg="blue")
label5.grid(row=5,column=6)

label6=Label(root,text="To find the DNA complement",fg="blue")
label6.grid(row=6,column=6)




button_1=Button(root,text="1",padx=40,pady=20,command=butto1)
button_2=Button(root,text="2",padx=40,pady=20,command=butto2)
button_3=Button(root,text="3",padx=40,pady=20,command=butto3)
button_4=Button(root,text="4",padx=40,pady=20,command=butto4)
button_5=Button(root,text="5",padx=40,pady=20,command=butto5)
button_6=Button(root,text="6",padx=40,pady=20,command=butto6)
button_0=Button(root,text="EXIT",fg="red",padx=40,pady=20,command=root.destroy)

button_1.grid(row=1,column=5)
button_2.grid(row=2,column=5)
button_3.grid(row=3,column=5)
button_4.grid(row=4,column=5)
button_5.grid(row=5,column=5)
button_6.grid(row=6,column=5)
button_0.grid(row=7,column=7)
root.mainloop()


