#!/usr/bin/python
import numpy as np
import random as rd
import os
import copy as cp

'''
vasp POSCAR manipulate
#example usage:

#alloy_layers('POSCAR_out','log',4.1625)
#alloy('POSCAR_out','log',4.1625,'1.vasp')
#get_layers('Alloy.vasp')
#stacking_faults('111.vasp','test.vasp',4.05)
#sf('1152.vasp')
#silica_surfaces('POSCAR',3)
#add_vac('POSCAR')
#remove_vac('CONTCAR')
#get_inter('POSCAR_out','Ag.vasp')
#get_surface_energy()
#fix_atom('POSCAR_out',lst,'POSCAR_out')  eg. where lst=[[1,2],[5,7]]

#get_inter('POSCAR_out','POSCAR_0',3.6)
#get_inter_vac('ZnO.vasp','Ag.vasp',2.3)   #2.3 distance between top ZnO atom and bot Ag atom
#zcenter('POSCAR_out','mass')  # move mass center of slab to 0.5 in z direction,mass file needed
#shift_z_align('POSCAR_out',7,9)  # align two atoms in xy plane(share same x,y coordinates)
#z_reverse('ZnO.vasp')
'''

def calc_dist(lst,n1,n2,H):
    if  lst[n2][0]-lst[n1][0]>=(0.5*H[0][0]):
        a=(lst[n1][0]-lst[n2][0]+H[0][0])
    elif  lst[n2][0]-lst[n1][0]<(-0.5*H[0][0]):
        a=(lst[n1][0]-lst[n2][0]-H[0][0])
    else:
        a=(lst[n1][0]-lst[n2][0])

    if  lst[n2][1]-lst[n1][1]>=(0.5*H[1][1]):
        b=(lst[n1][1]-lst[n2][1]+H[1][1])
    elif  lst[n2][1]-lst[n1][1]<(-0.5*H[1][1]):
        b=(lst[n1][1]-lst[n2][1]-H[1][1])
    else:
        b=(lst[n1][1]-lst[n2][1])

    if  lst[n2][2]-lst[n1][2]>=(0.5*H[2][2]):
        c=(lst[n1][2]-lst[n2][2]+H[2][2])
    elif  lst[n2][2]-lst[n1][2]<(-0.5*H[2][2]):
        c=(lst[n1][2]-lst[n2][2]-H[2][2])
    else:
        c=(lst[n1][2]-lst[n2][2])
    return np.sqrt(a**2+b**2+c**2)

def get_lmp(cell,str1,str2,data,filename='out.lmp'):
    lst_a=[int(s) for s in str2.split()]
    tot_num=int(sum(lst_a))
    with open(filename,'w') as fout:
        fout.write('comments\n\n    '+str(tot_num)+'    atoms\n\n   2   atom types\n\n')
        fout.write('0.000   '+str(cell[0][0])+'  xlo xhi\n')
        fout.write('0.000   '+str(cell[1][1])+'  ylo yhi\n')
        fout.write('0.000   '+str(cell[2][2])+'  zlo zhi\n'+str("%12.6f"% cell[1][0])+'  '+str("%12.6f" % cell[2][0])+'  '+str("%12.6f" % cell[1][2])+'   xy xz yz\n\n    Masses\n\n   1   26.9815\n   2   47.8670\n\nAtoms\n\n')
        nm=1
        for i in range(lst_a[0]):
            fout.write('    '+str(nm)+'     1       '+str(data[i][0])+'    '+str(data[i][1])+'    '+str(data[i][2])+'\n')
            nm+=1
        for i in range(lst_a[0],tot_num):
            fout.write('    '+str(nm)+'     2       '+str(data[i][0])+'    '+str(data[i][1])+'    '+str(data[i][2])+'\n')
            nm+=1
    return

def get_xsf(cell,str1,str2,data,filename='out.xsf'):
    lst_a=[int(s) for s in str2.split()]
    tot_num=int(sum(lst_a))
    with open(filename,'w') as fout:
        fout.write('CRYSTAL\nPRIMVEC\n')
        for i in range(3):
            fout.write(str(str( cell[i][0])+' '+str( cell[i][1])+' '+str( cell[i][2])+'\n'))
        fout.write(str('PRIMCOORD\n'+str(tot_num)+' 1\n'))
        for i in range(lst_a[0]):
            fout.write('Al      '+str(data[i][0])+'    '+str(data[i][1])+'    '+str(data[i][2])+'\n')
        for i in range(lst_a[0],tot_num):
            fout.write('Ti      '+str(data[i][0])+'    '+str(data[i][1])+'    '+str(data[i][2])+'\n')
    return

def get_POSCAR(cell,str1,str2,data,flag1=True,filename='POSCAR_out'):
    #flag1 determines whether Selective dynamics
    lst_a=[int(s) for s in str2.split()]
    tot_num=int(sum(lst_a))
    with open(filename,'w') as fout:
        fout.write('generated POSCAR\n')
        fout.write('1.0\n')
        for i in range(3):
            fout.write(str(str(cell[i][0])+' '+str(cell[i][1])+' '+str(cell[i][2])+'\n'))
        fout.write(str1)
        fout.write(str2)
        if flag1=='True':
            fout.write('Selective\nCart\n')
            for i  in range(tot_num):
                fout.write('    '+str(data[i][0])+'    '+str(data[i][1])+'    '+str(data[i][2]))
                fout.write('    T   T   T\n')
        elif flag1=='False':
            fout.write('Cart\n')
            for i  in range(tot_num):
                fout.write('    '+str(data[i][0])+'    '+str(data[i][1])+'    '+str(data[i][2]))
                fout.write('\n')
    return

def combine_str(str1,str2):
    lst_1=[str(s) for s in str1.split()]
    lst_2=[str(s) for s in str2.split()]
    tot_lst= lst_1 + lst_2
    mystring =' '.join(tot_lst)
    mystring += '\n'
    return mystring

def get_energy(filename='OUTCAR'):
    nrg=[]
    with open(filename,'r') as fin:
        for i in fin:
            if 'entropy=' in i:
                nrg.append(float(i.split()[6]))
    print nrg[-1]
    return float(nrg[-1])

def get_surface_energy(filename):
    cwd=os.getcwd()
    os.chdir(cwd+'/bulk')
    e1=get_energy()
    os.chdir(cwd+'/'+filename)
    e2=get_energy()
    cell = get_data('POSCAR')[0]
    area=cell[0][0]*cell[1][1]
    se=(e2-e1)*16/float(area) /2.0
    with open('DATA','a') as fin:
        fin.write(str(filename)+'energy is : ' + str(se) +'\n')
    return se

def get_low_high(cell,data):
    high= data.max(axis=0)[2]
    low= data.min(axis=0)[2]
    return low,high

def get_data(filename,f1=True):   #return cell,data,str1,str2  f1=='TRUE'?retrun cart;return direct
    with open(filename,'r') as fin:
        lines=fin.readlines()
    cell=np.zeros((3,3))
    for i in range(3):
        for j in range(3):
            cell[i][j]=lines[2+i].split()[j]
    lst_a=[int(s) for s in lines[6].split()]
    tot_num=int(sum(lst_a))
    data=np.zeros((tot_num,3))
    if lines[7][0] == 'C'  :
        for i in range(tot_num):
            for j in range(3):
                data[i][j]=float(lines[8+i].split()[j])
    elif lines[7][0] == str('D')  :
        data1=np.zeros((tot_num,3))
        for i in range(tot_num):
            for j in range(3):
                data1[i][j]=float(lines[8+i].split()[j])
        if f1=='True':

            #data1=np.dot(cell,np.transpose(data1))
            #data=np.transpose(data1)
            data = np.dot(data1,cell)
        else:
            data = data1

    elif lines[7][0]== 'S'  :
        if lines[8][0]== 'C' :
            for i in range(tot_num):
                for j in range(3):
                    data[i][j]=float(lines[9+i].split()[j])
        elif lines[8][0]== 'D'  :
            data1=np.zeros((tot_num,3))
            for i in range(tot_num):
                for j in range(3):
                    data1[i][j]=float(lines[9+i].split()[j])
            if f1=='True':
                #data1=np.dot(cell,np.transpose(data1))
                #data=np.transpose(data1)
                data = np.dot(data1,cell)
            else:
                data = data1
    return cell,data,lines[5],lines[6]

def add_vac(filename):
    cell=get_data(filename)[0]
    data=get_data(filename)[1]
    str1=get_data(filename)[2]
    str2=get_data(filename)[3]
    cell[2][2]+=5.0
    lst_a=[int(s) for s in str2.split()]
    tot_num=int(sum(lst_a))
    for i in range(tot_num):
        data[i][2]+=2.5
    get_POSCAR(cell,str1,str2,data,'False','POSCAR_out')
    return

def remove_vac(filename):
    cell=get_data(filename)[0]
    data=get_data(filename)[1]
    str1=get_data(filename)[2]
    str2=get_data(filename)[3]
#    cell1=get_low_high(cell,data)[0]
#    data1=get_low_high(cell,data)[1]
    low = get_low_high(cell,data)[0]
    high = get_low_high(cell,data)[1]
    for i in range(len(data)):
        data[i][2]-=float(low)
        data[i][2]+=0.1
    tmp1=high - low + 0.2
    cell[2][2]=tmp1
    get_POSCAR(cell,str1,str2,data,'False','POSCAR_out')
    return
 
def get_inter_vac(filename1,filename2,dist):
    cell_1=get_data(filename1)[0]
    data_1=get_data(filename1)[1]
    str1_1=get_data(filename1)[2]
    str2_1=get_data(filename1)[3]
    cell_2=get_data(filename2)[0]
    data_2=get_data(filename2)[1]
    str1_2=get_data(filename2)[2]
    str2_2=get_data(filename2)[3]
    str1=combine_str(str1_1,str1_2)
    str2=combine_str(str2_1,str2_2)
    high1=get_low_high(cell_1,data_1)[1]
    low1=get_low_high(cell_1,data_1)[0]
    high2=get_low_high(cell_2,data_2)[1]
    low2=get_low_high(cell_2,data_2)[0]
    lst_2=[int(s) for s in str2_2.split()]
    tot_num_2=int(sum(lst_2))
    for i in range(tot_num_2):
        data_2[i][2]+=high1
        data_2[i][2]+=dist
    data=np.concatenate((data_1,data_2),axis=0)
    cell_1[2][2]+=12.0
    cell_1[2][2]+=cell_2[2][2]
    get_POSCAR(cell_1,str1,str2,data,'False','POSCAR_out')
    return

def get_inter(filename1,filename2,a):
    cell_1=get_data(filename1)[0]
    data_1=get_data(filename1)[1]
    str1_1=get_data(filename1)[2]
    str2_1=get_data(filename1)[3]
    cell_2=get_data(filename2)[0]
    data_2=get_data(filename2)[1]
    str1_2=get_data(filename2)[2]
    str2_2=get_data(filename2)[3]
    str1=combine_str(str1_1,str1_2)
    str2=combine_str(str2_1,str2_2)
    high1=get_low_high(cell_1,data_1)[1]
    low1=get_low_high(cell_1,data_1)[0]
    high2=get_low_high(cell_2,data_2)[1]
    low2=get_low_high(cell_2,data_2)[0]
    lst_2=[int(s) for s in str2_2.split()]
    tot_num_2=int(sum(lst_2))
    for i in range(tot_num_2):
        data_2[i][2]+=high1
        data_2[i][2]+=2.0
    data=np.concatenate((data_1,data_2),axis=0)
    cell_1[2][2]=high1-low1+high2-low2 +float(a)
    get_POSCAR(cell_1,str1,str2,data,'False','POSCAR_out')
    return

def get_layers(filename,layers=3):
    cell=get_data(filename)[0]
    data=get_data(filename)[1]
    str1=get_data(filename)[2]
    str2=get_data(filename)[3]
    lst_a=[int(s) for s in str2.split()]
    tot_num=int(sum(lst_a))
#    layers=int(tot_num/16)
    print layers
    num=0
    n1=cell[2][2]/float(layers)
    data_2=np.zeros((tot_num))
    for i in range(layers):
        for j in range(tot_num):
            n1=cell[2][2]/float(layers)
            data[j][2]+=n1
            if data[j][2]>cell[2][2]: data[j][2]-=cell[2][2]
        out_name='layer_'+str(num)+'.vasp'
        get_POSCAR(cell,str1,str2,data,'False',out_name)
        num+=1

#        for j in range(tot_num):
#            data_2[j]=-data[j][2]
#            if data_2[i]<0: data_2[i]+=cell[2][2]
#        out_name='alloy_'+str(num)+'.vasp'
#        get_POSCAR(cell,str1,str2,data,'False',out_name)
#        num+=1
    return

def silica_surfaces(filename,n):
    cell=get_data(filename)[0]
    data=get_data(filename)[1]
    str1=get_data(filename)[2]
    str2=get_data(filename)[3]
    lst_a=[int(s) for s in str2.split()]
    tot_num=int(sum(lst_a))
    layers=int(n)
    print layers
    num=1
    n1=cell[2][2]/float(layers)
    data_2=np.zeros((tot_num))
    for i in range(layers):

        for j in range(tot_num):
            n1=cell[2][2]/float(layers)
            data[j][2]+=n1
            if data[j][2]>cell[2][2]: data[j][2]-=cell[2][2]
        keepgoing = True
        while keepgoing:
            get_POSCAR(cell,str1,str2,data,'False','POSCAR_tmp')
            add_vac('POSCAR_tmp')
            cell_t=get_data('POSCAR_out')[0]
            data_t=get_data('POSCAR_out')[1]
            str1_t=get_data('POSCAR_out')[2]
            str2_t=get_data('POSCAR_out')[3]
            lst_a_t=[int(st) for st in str2_t.split()]
            tot_num_t=int(sum(lst_a_t))
            at1=int(lst_a_t[0])
            at2=int(lst_a_t[1])
            high=get_low_high(cell_t,data_t)[1]
            low=get_low_high(cell_t,data_t)[0]
            
            flag1,flag2=True,True
            for ii in range(at1):
                if (data_t[ii][2]>=(high-2.0)) or (data_t[ii][2] <=(low+2.0)):
                    f1=0
                    for jj in range(at1,tot_num):
                        if calc_dist(data_t,ii,jj,cell_t) <= 2.0 : f1+=1
                    if f1==0 : 
                        print ii
                        flag1=False
                        pass
            for kk in range(at1,tot_num):
                if (data_t[kk][2]>=(high-2.0)) or (data_t[kk][2] <=(low+2.0)):
                    f2=0
                    for ll in range(at1):
                        if calc_dist(data_t,kk,ll,cell_t) <= 2.0 : f2+=1
                    if f2 <=1 :
                        print kk
                        flag2=False
                        pass
            print flag1 
            print flag2
            if flag1 and flag2 :
                filename_out='silica_'+str(num)+'.vasp'
                os.system('mv POSCAR_out '+str(filename_out))
                keepgoing = False
            else:
                for ii in range(tot_num):
                    data[ii][2]+=0.1
                    if data[ii][2]>cell[2][2]: data[ii][2]-=cell[2][2]
        num+=1
    return

def fix_atom(in1,lst_in,filename='POSCAR_sd'):
    cell=get_data(in1)[0]
    data=get_data(in1)[1]
    str1=get_data(in1)[2]
    str2=get_data(in1)[3]
    lst_a=[int(s) for s in str2.split()]
    tot_num=int(sum(lst_a))

    with open(filename,'w') as fout:
        fout.write('generated POSCAR\n')
        fout.write('1.0\n')
        for i in range(3):
            fout.write(str(str(cell[i][0])+' '+str(cell[i][1])+' '+str(cell[i][2])+'\n'))
        fout.write(str1)
        fout.write(str2)
        fout.write('Selective\nCart\n')

        for i  in range(tot_num):
            fout.write('    '+str(data[i][0])+'    '+str(data[i][1])+'    '+str(data[i][2]))
            a='    T   T   T\n' 
            for j in range(len(lst_in)):
                if data[i][2]<=lst_in[j][1] and data[i][2]>=  lst_in[j][0]:
                    a='    F   F   F\n'
            fout.write(a)
    return

   

def stacking_faults(in1,out,latt):
# latt for lattice constant 4.05 for Al in lammps 
    cell=get_data(in1)[0]
    data=get_data(in1)[1]
    str1=get_data(in1)[2]
    str2=get_data(in1)[3]

    for i in range(len(data)):
        if data[i][2]>=cell[2][2]/2.0:
            data[i][1]-=2*latt/np.sqrt(6)

    out1=out+'_sf.xsf'
    get_xsf(cell,str1,str2,data,out1)
    out2=out+'_sf.lmp'
    get_lmp(cell,str1,str2,data,out2)
    return 
def solid_solution(filename,conc):
    cell=get_data(filename)[0]
    data=get_data(filename)[1]
    str1=get_data(filename)[2]
    str2=get_data(filename)[3]
    lst_a=[int(s) for s in str2.split()]
    tot_num=int(sum(lst_a))
    atomtype=[0]*tot_num
    dop=int(conc*tot_num)
    for i in range(dop):
        atomtype[i]=1
    rd.shuffle(atomtype)
    data1=np.zeros((tot_num,3))
    ii=0
    for i in range(tot_num):
        if atomtype[i]==0:
            data1[ii][0]=data[i][0]
            data1[ii][1]=data[i][1]
            data1[ii][2]=data[i][2]
            ii+=1
    for i in range(tot_num):
        if atomtype[i]==1:
            data1[ii][0]=data[i][0]
            data1[ii][1]=data[i][1]
            data1[ii][2]=data[i][2]
            ii+=1
    str1='  Al Ti\n'
    str2='  '+str(tot_num-dop)+'    '+str(dop)+'\n'
    get_POSCAR(cell,str1,str2,data1,'False','ss.vasp')
    return

 
def sf(filename):
    solid_solution(filename,0.4) # which outputs a ss.vasp file
    get_layers('ss.vasp',12) # which outputs layer_*.vasp file
    os.system('rm ss.vasp')
    for i in range(12):
        in1='layer_'+str(i)+'.vasp'
        out1=str(i)+'.xsf'
        out2=str(i)+'.lmp'
        add_vac(in1)
        cell=get_data('POSCAR_out')[0]
        data=get_data('POSCAR_out')[1]
        str1=get_data('POSCAR_out')[2]
        str2=get_data('POSCAR_out')[3]
        get_xsf(cell,str1,str2,data,out1)
        get_lmp(cell,str1,str2,data,out2)
        stacking_faults('POSCAR_out',str(i),4.05)

        os.system('rm POSCAR_out')
        os.system('rm '+str(in1))
    os.system('mkdir xsf lmp_input')
    os.system('mv *.xsf xsf')
    os.system('mv *.lmp lmp_input')

    return

def alloy(filename,log,LC,out):
    cell=get_data(filename)[0]
    data=get_data(filename)[1]
    str1=get_data(filename)[2]
    str2=get_data(filename)[3]
    lst_a=[int(s) for s in str2.split()]
    tot_num=int(sum(lst_a))
    data[:,0] /= cell[0][0]
    data[:,1] /= cell[1][1]
    data[:,2] /= cell[2][2]
    #atomtype=[0]*lst_a[0]+[1]*lst_a[1]+[3]*lst_a[2]
    #print atomtype
    with open(log,'r') as fin2:
        mylines2=fin2.readlines()
    data_tmp=np.zeros((96,4))
    c1=int(lst_a[0]+lst_a[1])
    print c1
    for i in range(96):
        data_tmp[i][0]=data[i+c1][0]
        data_tmp[i][1]=data[i+c1][1]
        data_tmp[i][2]=data[i+c1][2]
        data_tmp[i][3]=int(mylines2[i].split()[0])
    
    data_tmp = sorted(data_tmp, key=lambda x:x[3], reverse=True)

    a,b=LC*4/np.sqrt(2),LC*2*np.sqrt(3)/np.sqrt(2)
    c=LC/np.sqrt(3)*6
    x=cell[2][2]-4.1625*6/np.sqrt(3)
    c+=x
    
    f=open(out,'w')
    f.write('#most\n1.0\n%10.6f    0.0 0.0\n0.0    %10.6f     0.0\n0.0    0.0     %10.6f\n' %(a,b,c))
    f.write('O Si Al Ag\n 64 32 10 86\n')
    f.write('Direct\n')
    for i in range(c1):
        f.write(str(data[i][0])+ '   '+str(data[i][1])+'   '+str(data[i][2])+'\n')
    for i in range(96):
        f.write(str(data_tmp[i][0])+ '   '+str(data_tmp[i][1])+'   '+str(data_tmp[i][2])+'\n')
    f.close()
    return

def alloy_layers(filename,log,LC):

    with open(log,'r') as fin2:
        mylines1=fin2.readlines()
    for i in range(6):
        outname='log'+str(i)
        num=i*16
        with open(outname,'w') as fin3:
            for j in range(96-num):
                fin3.write(mylines1[num+j])
            for j in range(num):
                fin3.write(mylines1[j])
        outname2=str(i)+'.vasp'
        alloy(filename,outname,4.1625,outname2)
    return

def zcenter(filename,mass):
    cell=get_data(filename)[0]
    data=get_data(filename)[1]
    str1=get_data(filename)[2]
    str2=get_data(filename)[3]
    lst_a=[int(s) for s in str2.split()]
    tot_num=int(sum(lst_a))


    with open(mass,'r') as fin2:
        mylines1=fin2.readlines()
    mass = []
    tot_m=0
    for i in range(len(lst_a)):
        mass.append(float(mylines1[i].split()[0]))
        tot_m += mass[i]*lst_a[i]
    print tot_m
    tmp1 = 0

    c=0
    for i in range(len(lst_a)):
        for j in range(lst_a[i]):
            tmp1+=data[c][2]*mass[i]
            c+=1
    print c, tmp1
    tmp2=tmp1/float(tot_m)
    print tmp2
    dist=cell[2][2]/2.0-tmp2
    print dist
    data[:,2]+=dist


    get_POSCAR(cell,str1,str2,data,'False','POSCAR_out')
    return


def shift_z_align(filename,at1,at2):
    # move at2 to at1's (x,y)
    cell=get_data(filename)[0]
    data=get_data(filename)[1]
    str1=get_data(filename)[2]
    str2=get_data(filename)[3]
    lst_a=[int(s) for s in str2.split()]
    tot_num=int(sum(lst_a))
    dx = data[at1][0] - data[at2][0]
    dy = data[at1][1] - data[at2][1]
    for i in range(lst_a[0]+lst_a[1],tot_num):
        data[i][0] +=dx
        data[i][1] +=dy
    get_POSCAR(cell,str1,str2,data,'False','POSCAR_out')
    return
def tear_interface(filename,num1):
    cell=get_data(filename)[0]
    data=get_data(filename)[1]
    str1=get_data(filename)[2]
    str2=get_data(filename)[3]
    lst_a=[int(s) for s in str2.split()]
    tot_num=int(sum(lst_a))
    tot_num1 = 0
    for i in range(num1):
        tot_num1+=lst_a[i]
    tot_num2=tot_num-tot_num1
    data_1=np.zeros((tot_num1,3))
    data_2=np.zeros((tot_num2,3))
    for i in range(tot_num1):
        data_1[i][0] = data[i][0]
        data_1[i][1] = data[i][1]
        data_1[i][2] = data[i][2]
    for i in range(tot_num1,tot_num):
        data_2[i-tot_num1][0] = data[i][0]
        data_2[i-tot_num1][1] = data[i][1]
        data_2[i-tot_num1][2] = data[i][2]
    tmp = []
    for i in range(num1):
        tmp.append(str1.split()[i])
    str1_1 = '\t'.join(str(x) for x in tmp)+'\n'
    tmp = []
    for i in range(num1,len(str1.split())):
        tmp.append(str1.split()[i])
    str1_2 = '\t'.join(str(x) for x in tmp)+'\n'
    tmp = []
    for i in range(num1):
        tmp.append(str2.split()[i])
    str2_1 = '\t'.join(str(x) for x in tmp)+'\n'
    tmp = []
    for i in range(num1,len(str2.split())):
        tmp.append(str2.split()[i])
    str2_2 = '\t'.join(str(x) for x in tmp)+'\n'

    get_POSCAR(cell,str1_1,str2_1,data_1,'False','1.vasp')
    get_POSCAR(cell,str1_2,str2_2,data_2,'False','2.vasp')
    return

def z_reverse(filename):
    cell=get_data(filename)[0]
    data=get_data(filename)[1]
    str1=get_data(filename)[2]
    str2=get_data(filename)[3]
    lst_a=[int(s) for s in str2.split()]
    tot_num=int(sum(lst_a))
    for i in range(tot_num):
        data[i][2]=cell[2][2]-data[i][2]
    get_POSCAR(cell,str1,str2,data,'False','reverse.vasp')
    return

def z_rotate(filename,rg,tht,box):
    cell=get_data(filename,'True')[0]
    data=get_data(filename,'True')[1]
    str1=get_data(filename,'True')[2]
    str2=get_data(filename,'True')[3]
    lst_a=[int(s) for s in str2.split()]
    tot_num=int(sum(lst_a))
    radi_tht=tht*np.pi/180.0
    R_matrix=[[np.cos(radi_tht),-np.sin(radi_tht),0],
              [np.sin(radi_tht),np.cos(radi_tht),0],
              [               0,               0,1]]
    print R_matrix
    print data
    #data1=cp.copy(data)
    #data1=np.dot(R_matrix,np.transpose(data1))
    #data1=np.transpose(data1)
    for i in range(tot_num):
        if data[i][2]<=rg[1] and data[i][2]>=rg[0] :
            data[i]=np.dot(R_matrix,np.transpose(data[i]))
            data[i]=np.transpose(data[i])
    if box == 'True':
        cell1=np.dot(R_matrix,np.transpose(cell))
        cell1=np.transpose(cell1)
        get_POSCAR(cell1,str1,str2,data,'False','rotated.vasp')
    else:
        get_POSCAR(cell,str1,str2,data,'False','rotated.vasp')
    return

def get_H_int(filename,lst,r,out):
    cell=get_data(filename,'True')[0]
    data=get_data(filename,'True')[1]
    str1=get_data(filename,'True')[2]
    str2=get_data(filename,'True')[3]
    lst_a=[int(s) for s in str2.split()]
    tot_num=int(sum(lst_a))
    f=open(out,'w')
    for i in range(len(lst)):
        for j in range(i+1,len(lst)):
            tmp_d=calc_dist(data,lst[i],lst[j],cell)
            if tmp_d< r:
                f.write('%12.6f  %12.6f  %12.6f\n'%((data[lst[i]][0]+data[lst[j]][0])/2.0,(data[lst[i]][1]+data[lst[j]][1])/2.0,data[lst[i]][2]-1.2))
    f.close()
    return
#get_H_int('full.vasp',[61,85,49,73,60,84],3,'out')

def add_dist_type(in1,lst_in,filename='POSCAR_sd'):
    cell=get_data(in1)[0]
    data=get_data(in1)[1]
    str1=get_data(in1)[2]
    str2=get_data(in1)[3]
    lst_a=[int(s) for s in str2.split()]
    tot_num=int(sum(lst_a))

    with open(filename,'w') as fout:
        fout.write('generated POSCAR\n')
        fout.write('1.0\n')
        for i in range(3):
            fout.write(str(str(cell[i][0])+' '+str(cell[i][1])+' '+str(cell[i][2])+'\n'))
        fout.write(str1)
        fout.write(str2)
        fout.write('Selective\nCart\n')

        for i  in range(tot_num):
            for j in range(len(lst_in)):
                if data[i][2]<=lst_in[j][1] and data[i][2]>=  lst_in[j][0]:
                    fout.write('    '+str(data[i][0])+'    '+str(data[i][1])+'    '+str(data[i][2]+0.4))
                else:
                    fout.write('    '+str(data[i][0])+'    '+str(data[i][1])+'    '+str(data[i][2]))
                fout.write('\n')

    return




add_dist_type('POSCAR',[[63.8,75]],'POSCAR_out')
#z_reverse('ZnO.vasp')
#get_inter_vac('reverse.vasp','Ag.vasp',2.3)
#zcenter('POSCAR_out','mass')
#fix_atom('POSCAR_out',[[10,12],[16,17]],'POSCAR_out')
#rotate('a.vasp',(0,0,1),45,(1,1,0),'True')
#z_rotate('slab.vasp',[16,30],30,'False')

#shift_z_align('POSCAR_out',4,8)
#tear_interface('last.vasp',2)
