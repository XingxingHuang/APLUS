from bpz_tools import *
import sys,os,shutil

"""
   Usage:
   python gen_trans.py filter_name filter.dat ccd.dat atm.dat
   Generates a file called filter_name.res which contains the 
   normalized product of all the transmissions filter.dat ccd.dat atm.dat
"""

try:
    filter_name=sys.argv[1]
    ntrans=len(sys.argv[2:])
except:
    print """
Usage:
    python gen_trans.py filter_name filter.dat ccd.dat atm.dat

    Generates a file called filter_name.res which contains the 
    normalized product of all the transmissions filter.dat ccd.dat atm.dat
    """ 
    sys.exit()
    
l_t=[[]]*ntrans
f_t=[[]]*ntrans
trans=[[]]*ntrans

for i in range(ntrans):
    trans[i]=sys.argv[2+i]
    l_t[i],f_t[i]=get_data(trans[i],range(2))
    #Make a security copy of the filter if there isn't one already
    if not os.path.exists(trans[i]+'.bak'):
	shutil.copy(trans[i],trans[i]+'.bak')
    if not ascend(l_t[i]):
	print 'File ',trans[i]
	print 'wavelengths not ordered'
	print 'Ordering them properly'
	ind=argsort(l_t[i])
	l_t[i]=take(l_t[i],ind)
	f_t[i]=take(f_t[i],ind)

#Check that the wavelenght ranges and units are right
for i in range(ntrans):
    print 'File: ',trans[i]
    print l_t[i][0],' < lambda < ',l_t[i][-1]
    if i==0: lmin,lmax=l_t[i][0],l_t[i][-1]
    if not ask('Are the units right(y/n)?\n'):
	factor=raw_input('Factor to multiply lambda by?\n')
	l_t[i]=l_t[i]*float(factor)

#Normalize and Multiply the transmissions
f=ones((len(l_t[0]),))
for i in range(ntrans):
    f_t[i]=f_t[i]/max(f_t[i])
    #Save the final transmissions
    put_data(trans[i],(l_t[i],f_t[i]))
    header=get_header(trans[i])
    if len(header)>0: put_header(trans[i],header)
    #Find the total throughput
    f=f*match_resol(l_t[i],f_t[i],l_t[0])

#Get rid of negative/zero values in the filter transmission
#(Saves time when estimating the object fluxes)
f=clip(f,0.,1e40)
high=greater(f,0.)
i1=searchsorted(high,1)
i2=searchsorted(high[::-1],1)
if i2: put_data(filter_name+'.res',(l_t[0][i1:-i2],f[i1:-i2]))
else: put_data(filter_name+'.res',(l_t[0][i1:],f[i1:]))
shutil.copy(filter_name+'.res',fil_dir)

#Plot them
buffer="""del 0-1000
device xterm
expand 1.1
lw 1
erase limits """+str(lmin)+" "+str(lmax)+" 0 1.05 \n"
for i in range(ntrans):
    buffer=buffer+("data "+trans[i]+'\n'+
		   "read {l"+str(i)+" 1 f"+str(i)+" 2}\n"+
		   "ltype "+str(i+1)+" con l"+str(i)+" f"+str(i)+"\n"+
		   "ltype 0 box\n")
buffer=buffer+"""data """+filter_name+""".res
read {l 1 f 2}
lw 2 ltype 0
con l f box"""
open('plot_trans.sm','w').write(buffer)

print "Super mongo commands to plot the transmissions are written in plot_trans.sm"

#os.system('sm execute plot_trans.sm&')
