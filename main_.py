#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Nov 13 09:37:03 2022

@author: taniapeniche
"""

######################CONSTANTS#################
#INPUT MESSAGE
START_MSG = "Load input from file (Y/N)?"
INPUT_MSG = "> "
SOURCE_MSG = "Input source (F/T):"
D_STRUCT_MSG="How many structures (enter a positive number):"
D_IN_MSG = "Give two SMILES to compare:"
F = "F"
T = "T"
Y ="Y"
N = "N"
C = "C"
M = "M"
D = "D"
S = "S"
I = "I"
H = "H"
Q = "Q"

#CONSTANTS
simb = ["B", "C", "N", "O", "P", "S", "F", "I","H","c"]
BR = "Br"
CL = "Cl"
BR_L = ["br","Br","BR"]
CL_L = ["cl", "Cl","CL"]
char = ["-", "=", "#", ":", "/"]
char_ = "\\"
par = ["(", ")"]
qpar = ["[","]"]
TC = "@"
num = ["1","2","3","4","5","6"]
IS = " is "
CONTAINS = " contains "

#OUTPUT MESSAGE
MENU = """C: count the number of times each sub-string from an external list (given file) occurs in the SMILES strings of the list.
M: Count the number of times each atomic element occurs in the strings in the list and obtain the molecular formula (number of atoms of each element, e.g., C8NO2). The output of the command should appear in the terminal and be in lexicographic order.
D: compare a given pair of molecules from their SMILES representation (calculate their dissimilarity, i.e., sum of squared differences between the number of occurrences of the sub-strings in two SMILES).
I: input a new SMILES string to be added to the current list, if valid (if not, the application reports it found a problem and waits for the user’s to input a new command).
H: help – list all commands.
Q: quit – quit the application."""
D_OUT = "Dissimilarity degree between SMILES "
I_OUT = "SMILES list updated: ", " inserted"
Q_OUT = "Save SMILES list to file (Y/N)?"
GB = "Goodbye."
AND = " and "
TIME = " time"
TIMES = " times"
VI = ","
DP = ": "
WRT = " w.r.t "


#ERROR MESSAGE
ERROR = "ERROR"
INV_ANSWER = "Answer invalid"
INV_INPUT = "Input Invalid"
READ_ERROR = "Failed reading file "
EMPTY_LIST = "SMILES list empty"
INV_COMMAND = "Command invalid"
D_ERROR = "SMILES list empty or singular"
D_SMILES_ERROR="SMILES unknown"
D_VALID_ONE="input a valid one:"
I_NOT_VAL = "String ", " is not a valid SMILES"
I_EX = "SMILES "," already loaded"
WRITE_ERROR = "Failed writing file "
SM_UKN = "SMILES "," unknown; input a valid one:"


#####################FUNCTIONS#####################
def read_from_file(file_name):
    smiles_list = []
    try:
        file_handle = open(file_name, "r")
    except:
        return ERROR
    lines = file_handle.read().split()
    n_SMILES = int(lines.pop(0))
    for i in range(n_SMILES):
        s = str(lines.pop(0))
        if val_smiles(s) == True:
            smiles_list.append(s)
    file_handle.close()
    return smiles_list


#Write to a file
def write_to_file(SMILES_list, file_name):
    file_handle = open(file_name, "w")
    file_handle.write(str(len(SMILES_list))+"\n")
    for s in SMILES_list:
        file_handle.write(s + "\n")
    file_handle.close()


#Upper_Case
def upper_smiles(s):
    i = 0
    l = len(s)
    s_up = ""
    while i < l:
        if i < (l-1):
            if s[i]+ s[i+1]  in BR_L or s[i]+ s[i+1] in CL_L:
                s_up += s[i].upper()
                s_up += s[i+1].lower()
                i +=2
            else:
                s_up += s[i].upper()
                i +=1
        elif i == (l-1):
            s_up+=s[i].upper()
            i +=1
    return s_up


#Count the number of a specif character
def count(a,s):
    count = 0
    l = len(s)
    for i in range(l):
        if s[i] == a:
            count +=1
    return count


#Validate the SMILES
def val_smiles(s):
    l = len(s)
    c = 0
    start = is_first_atm(s)
    c += start
    i = start
    if start == 0:
        return I_NOT_VAL[0] + s + I_NOT_VAL[1]
    elif start == l:
        return c == l
    elif is_last_atm(l, s) == False:
        c = 0
    else: 
        while i < l:
            if i<(l-1) and is_doub_simb(i, s) == True:
                i +=2
                c +=2
            elif is_simb(i, s) == True:
                i +=1
                c +=1
            elif is_par(i, s) == True:
                if val_par(i, s)==True:
                    i +=1
                    c +=1
                else:
                    break
            elif is_qpar(i, s) == True:
                if val_qpar(i, s)==True:
                    i +=1
                    c +=1
                else:
                    break
            elif is_char(i , s) == True:
                if val_char(i, s) == True:
                    i +=1
                    c +=1
                else:
                    break
            elif is_char_(i, s) == True:
                if val_char_(i, s)== True:
                    i +=1
                    c +=1
                else:
                    break           
            elif is_TC(i, s) == True:
                if val_TC(i, s) == True:
                    i +=1
                    c +=1
                else:
                    break
            elif is_num(i, s) == True:
                if val_num(i, s,l) == True:
                    i +=1
                    c +=1
                else:
                    break
            else:
                i = l 
    return l == c
        

#check if the first character is valid
def is_first_atm(s):
    if s[0] not in simb and s[0]+s[1] != BR and s[0]+s[1] != CL:
        i=0
    else:
        if s[0]+s[1] == BR or s[0]+s[1] == CL:
            i=2
        else:
            i=1
    return i


#check if the last character is valid
def is_last_atm(l,s):
    if s[l-1] not in simb and s[l-1] not in num and  s[l-2]+s[l-1] != BR and s[l-2]+s[l-1] != CL and s[l-1] != par[1]:
        return False

#check if the character is an atom
def is_simb(i,s):
    return s[i] in simb

#check if the character is an atom with two characters
def is_doub_simb(i,s):
    if s[i]+s[i+1] == BR or s[i]+s[i+1] == CL:
        return True

#count the # of [] or ()
def find_tp(s, tp):
    c = 0
    c_i = 0
    for i in range(len(s)):
        if s[i] == tp[0]:
            c +=1
        elif s[i] == tp[1]:
            c_i +=1
    if c != c_i:
        return False
    else:
        return True
    
#check if the character is a ()
def is_par(i,s):
    return s[i] in par

#check if the () is valid
def val_par(i,s):
    if find_tp(s,par)==True:
        if (s[i] == par[0] and s[i+1] == par[1]):
            return False
        return True

#check if the character is a []
def is_qpar(i,s):
    return s[i] in qpar

#check if the [] is valid
def val_qpar(i,s):
    if find_tp(s,qpar)==True:
        if (s[i] == qpar[0] and s[i+1] == qpar[1]):
            return False
        return True

#check if the character is a simbol
def is_char(i,s):
    return s[i] in char

#check if the simbol in the position is valid
def val_char(i,s):
    if s[i+1] in char or s[i+1] in TC:
        return False
    elif s[i-1] in char or s[i-1] in TC:
        return False
    return True

#check if the character is \
def is_char_(i,s):
    return s[i] == char_

#check if the \ in the position is valid
def val_char_(i,s):
    if s[i+1] == char_:
        return False
    else:
        return True

#Check if is an @
def is_TC(i,s):
    return s[i]==TC

#check if the @ in the position is valid
def val_TC(i,s):
    c = count(TC, s)
    if c > 2 and s[i+1]==TC and s[i-1]==TC:
        return False
    else:
        return True

#Check if is an number
def is_num(i,s):
    return s[i] in num

#check if the quantities of number and if number is valid
def val_num(i,s,l):
    c = count(s[i],s)
    if c%2 != 0 and l<8:
        return False    
    return True

#Count if subtring in SMILES sequence
def count_s(s,sb):
    """ structures: str,list-> str
    Description: 
    Example: 
    """
    ls = len(s)
    lsb = len (sb)
    lsb_ = lsb-1
    i = 0
    count = 0
    if ls < lsb:
        return str(count) 
    else:
        while i+lsb_ <= ls:
            if s[i:lsb] == sb:
                count+=1
            i+=1
            lsb +=1
        return str(count)

#Returns the molecular formule
def molecular_formule(s):
    form = ""
    form_c = create_elements(s)
    elem_sorted = sorted(form_c)
    for i in elem_sorted:
        if form_c[i] == 1 :
            form += i 
        else:
            form += i + str(form_c[i])
    return form

#creates a dicionary with element and is respective count
def create_elements(s):
    form = {}
    i = 0
    l =len(s)
    while i<(l-1):
        if s[i] not in form and s[i]+s[i+1] == BR:
            n = count_s(s, BR)
            form[BR] = int(n)
            i+=2
        elif s[i] not in form and s[i]+s[i+1] == CL:
            n = count_s(s, CL)
            form[CL] = int(n)
            i+=2
        elif s[i] not in form and s[i] in simb:
            n = count(s[i], s)
            form[s[i]] = n
            i+=1
        else:
            i+=1
    if i == (l-1) and s[i] in simb:
        n = count(s[i], s)
        form[s[i]] = n
        i+=1
    return form

#the number of times that sb occurs in s
def count_smiles_split(s,sb):
    para = slice_par(s)
    ss=slice_other(s, para)
    para.append(ss)
    l_smiles = out_bonds(para)
    c = 0
    for i in l_smiles:
        c += int(count_s(i, sb))
    return str(c)

#positions of (),[]
def find_par(s,para):
    l = len(s)
    pos_par0 = []
    pos_par1 = []
    for i in range(l):
        if s[i] == para[0]:
            pos_par0.append(i)
        elif s[i] == para[1]:
            pos_par1.append(i)
    return pos_par0, pos_par1

#Give us what is btw (),[]
def slice_par(s):
    pp0, pp1 = find_par(s, par)
    pq0, pq1 = find_par(s, qpar)
    l = len(pp0)
    lq = len(pq0)
    para = []
    for i in range(l):
         ss = s[pp0[i]:(pp1[i]+1)]
         para.append(ss)
    for j in range(lq):
        ss = s[pq0[i]:(pq1[i]+1)]
        para.append(ss)
    return para

#Give us the strings withou the things btw ()[]
def slice_other(s,para):
    s_spliced = s
    for i in para:
        s_spliced=s_spliced.replace(i,"")
    return s_spliced

#take the bonds to allow us to count
def out_bonds(para):
    smile = []
    ss =""
    for s in para:
        for i in s:
            if i not in char[0:3]:
                ss+=i
        smile.append(ss)
    return smile

#Counts the number of each subtring occurs in each smile, returns a string with all the information
def count_all_s(s,sbs):
    """ structures: str,list-> str
    Description: 
    Example: 
    """
    a = ""
    l = len(sbs)
    for j in s:
      sub = j + CONTAINS
      for i in range(0,l):
          count = count_smiles_split(j,sbs[i])
          if i == (l-1):
            if int(count) == 1 :
                sub += AND + sbs[i] +" "+ count + TIME
            else:
                sub += AND + sbs[i] +" "+ count + TIMES
          else:
            if int(count) != 1:
                sub += sbs[i] +" "+ count + TIMES
            else:
                sub += sbs[i] +" "+ count + TIME
      if j == s[len(s)-1]:
        a += sub
      else:
        a += sub + "/n"
    return sorted(a.split("/n"))

#Define the count in case of single,double or triple bounds:
def check_bonds_count(s,ss):
    """ check_bonds: str,list-> str
    Description: 
    Example: 
    """
    l = len(s)
    pos = []
    count = count_smiles_split(s,ss)
    for i in range(l):
        if s[i] in char[0:3]:
            if s[i-1] in ss and s[i+1] in ss:
                pos.append(i)
    lp = len(pos)
    if lp > 0:
        count = count_bonds(s, pos, count)
    return count

#give us the weight of the bonds
def count_bonds (s,pos, count):
    count = int(count)
    for i in pos:
      if s[i] == char[1]:
          count +=2
      elif s[i] == char[2]:
          count +=3
    return count

#Count how many times x sub-structures are in the 2 smiles structures (S1,S2):
def count_differences(s,sst):
    """ count_differences: str,list-> list
    Description: 
    Example: 
    """
    count=[]
    for ss in sst:
      c=check_bonds_count(s, ss)
      count.append(c)
    return count

#Calculate the dissimilarity degree between SMILES s1 and s2 w.r.t. ss1 and ss2:
#n é o resultado do cálculo
#count the occurence of each sub-string in the SMILES
def differences(s1,s2,sst):
    """ differences: str,str,list-> int
    Description: 
    Example: 
    """
    sst_s1=count_differences(s1,sst)
    sst_s2=count_differences(s2,sst)
    n = 0
    for i in range(len(sst)):
        n += (int(sst_s1[i])-int(sst_s2[i]))**2 
    sst_s=""
    l = len(sst)-1
    sst = sorted(sst)
    for i in sst:
      if i != sst[l]:
        sst_s += i + AND
      else:
        sst_s +=i
    s = sorted([s1,s2])
    ss =""
    for i in s:
      if i != s[len(s)-1]:
        ss += i + AND
      else:
        ss +=i 
    return D_OUT + ss + WRT + sst_s + DP + str(n)




############Interections with the user###########

#Help Command/Start:
def help_M():
    print(MENU)

#Next Command
def next_command():
    return str(input(INPUT_MSG)).upper()

#Input/ Create Smiles List
def input_from():
    smiles_list = []
    sel = str(input(START_MSG)).upper()
    while sel != Y and sel != N:
        print(INV_ANSWER)
        sel = str(input(START_MSG)).upper()
    if sel == Y:
        file_name = str(input(INPUT_MSG))
        out = read_from_file(file_name)
        if out == ERROR:
            print(READ_ERROR + file_name)
            return smiles_list
        else:
            smiles_list = out
            if smiles_list == []:
                print(EMPTY_LIST)
            else:
                smiles_list=sorted(smiles_list)
                print(smiles_list)
            return smiles_list
    else:    
        return smiles_list
    

#Command Q
def quit_app(SMILES_list):
    op = str(input(Q_OUT)).upper()
    while op != Y and op != N:
        print(INV_ANSWER)
        op = str(input(Q_OUT)).upper()
    if op == Y:
        file_name = str(input(INPUT_MSG))
        try:
            write_to_file(SMILES_list, file_name)
        except:
            print(WRITE_ERROR + file_name)
        print(GB)
    else:
        print(GB)


#Comnmand I
def insert_new_smile(SMILES_list):
    s = upper_smiles(str(input()))
    if s not in SMILES_list:
        if val_smiles(s) == True:
            SMILES_list.append(s)
            print(I_OUT[0] + s + I_OUT[1])
        else:
            print(I_NOT_VAL[0] + s + I_NOT_VAL[1])
    else:
        print(I_EX[0]+ s + I_EX[1])


#Command M
def molecular_formule_io (smiles_list):
    form = []
    if len(smiles_list)==0:
        print(EMPTY_LIST)
    else:
        for s in smiles_list:
            form.append(s+IS+molecular_formule(s))
        for i in form:
            print(i)
            
#Command C
#Number of substring to compare in smiles:
def number_ssb_smiles_io(smiles_list):
    command_app=str(input(SOURCE_MSG)).upper()
    while command_app!=F and command_app!=T:
        print(INV_INPUT)
        command_app=str(input(SOURCE_MSG)).upper()
    if command_app==F:
        file_name = str(input(INPUT_MSG))
        subs=read_from_file(file_name)
    else:
        subs=read_sub_from_terminal()
    c = count_all_s(smiles_list, subs)
    for i in c:
        print(i)
    
    
def read_sub_from_terminal():
    num_subs=int(input())
    subs=[]
    for i in range(0,num_subs):
        sub=str(input())
        subs.append(sub)
    return subs
    

#Command D
def dissimilarity_io(smiles_list):
    if len(smiles_list)<2:
        print(D_ERROR)
    else:
        num_struc=int(input(D_STRUCT_MSG))
        while num_struc<=0:
            num_struc=int(input(D_STRUCT_MSG))
        if num_struc>0:
            sst=structures(num_struc)
            print(D_IN_MSG)
            s1=s1_(smiles_list)
            s2=s2_(smiles_list)
            print(differences(s1,s2,sst))

#Input the respective sub-structures to search in the list of smiles
def structures(num_struc):
    """ structures: int -> list
    Description: Consoant the number of sub-structures to compare into smiles, 
    inputs valid sub-structures into a list sst until there aren't no many 
    substructures to input, otherwise prints "input a valid one:"
    to input a valid substructure.
    Example: structure(4)->CC,CBr,C,Cl->[CC,CBr,C,Cl]
    """
    sst=[]
    while num_struc>0:
        c_smiles=str(input()).upper()
        if val_smiles(c_smiles)==True:
            sst.append(c_smiles)
            num_struc-=1
        else:
            print(D_VALID_ONE)
            c_smiles=str(input()).upper()
    return sst

def s1_(smiles_list):
    """ s1_s2: list -> str,str
    Description: Input 2 smiles variables to compare with eachother, but if s1 is not present in the list of SMILES, print 
    "SMILES s1 unknown;input a valid one:" and submit a valid s1 SMILE, or s2 is not present also in the list of smiles, 
    print "SMILES s2 unknown;input a valid one:" and submit a valid s2 SMILE, otherwise just return the 2 valid SMILES, s1 and s2.
    Example:
    """
    s1=upper_smiles(str(input()))
    while s1 not in smiles_list:
       print(SM_UKN[0]+s1+SM_UKN[1])
       s1=upper_smiles(str(input()))
    return s1

def s2_(smiles_list):
    s2=upper_smiles(str(input()))
    while s2 not in smiles_list:
        print(SM_UKN[0]+s2+SM_UKN[1])
        s2=upper_smiles(str(input()))
    return s2



######################Main#####################
def main():
    smiles_list=input_from()
    command = help_M()
    while command != Q:
        if command == H:
            help_M()
            command = next_command()
        elif command==C:
            number_ssb_smiles_io(smiles_list)
            command = next_command()
        elif command == M:
            molecular_formule_io(smiles_list)
            command = next_command()
        elif command == D:
            dissimilarity_io(smiles_list)
            command = next_command()
        elif command == I:
            insert_new_smile(smiles_list)
            command = next_command()
        else:
            print(INV_COMMAND)
            command = next_command()
    quit_app(smiles_list)
    

#################Execute#############
main()








