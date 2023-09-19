_author__    = "Dr. Pham Quang Phuong","Cao Anh Quoc Hung"
__copyright__ = "Copyright 2023"
__license__   = "All rights reserved"
__email__     = "phuong.phamquang@hust.edu.vn","caoanhquochung@gmail.com"
__status__    = "Released"
__version__   = "1.1"

# IMPORT -----------------------------------------------------------------------
import sys,os,time
import winreg
import tkinter as tk
import tkinter.filedialog as tkf
import tkinter.messagebox as tkm
from tkinter import ttk
import traceback
import subprocess,threading
import new_kernel as KERNEL
import warnings
warnings.simplefilter("ignore", UserWarning)
import new_optimum
import argparse
PARSER_INPUTS = argparse.ArgumentParser(epilog= "")
PARSER_INPUTS.usage = 'BK Power System Analysis'
ARGVS = PARSER_INPUTS.parse_known_args()[0]
#
class WIN_REGISTRY:
    def __init__(self,path,keyUser,nmax):
        #
        self.nmax = nmax
        #
        if keyUser == "LOCAL_MACHINE":
            self.Registry = winreg.ConnectRegistry(None, winreg.HKEY_LOCAL_MACHINE)
        else:#"CURRENT_USER":
            self.Registry = winreg.ConnectRegistry(None, winreg.HKEY_CURRENT_USER)
        #
        try:
            self.RawKey  = winreg.OpenKey(self.Registry, path)
        except:
            self.createKey(path)
            self.RawKey  = winreg.OpenKey(self.Registry, path)
        #
        self.reg_key = winreg.OpenKey(self.Registry,path,0, winreg.KEY_SET_VALUE)
    #
    def createKey(self,path):
        patha = str(path).split("\\")
        path1 = patha[0]
        for i in range(1,len(patha)):
            p2 = path1 +"\\"+ patha[i]
            try:
                access_key = winreg.OpenKey(self.Registry,p2)
            except:
                access_key = winreg.OpenKey(self.Registry,path1)
                winreg.CreateKey(access_key,patha[i])
            #
            path1 = p2
    #
    def getAllNameValue(self):
        i = 0
        name,vala = [],[]
        while True:
            try:
                a = winreg.EnumValue(self.RawKey, i)
                name.append(a[0])
                vala.append(a[1])
            except:
                break
            i+=1
        return name,vala
    #
    def getAllValue(self):
        return self.getAllNameValue()[1]
    #
    def getValue0(self):
        try:
            return self.getAllNameValue()[1][0]
        except:
            return ""
    #
    def appendValue(self,val):
        name,vala = self.getAllNameValue()
        if len(vala)>0 and val==vala[0]:
            return False
        #
        for n1 in name:
            winreg.DeleteValue(self.reg_key, n1)
        #
        r1 = [val]
        for ri in vala:
            if len(r1)>=self.nmax:
                break
            if ri!=val :
                r1.append(ri)
        #
        for i in range(len(r1)):
            winreg.SetValueEx(self.reg_key, "File"+str(i+1), 0, winreg.REG_SZ, r1[i])
        #
        return True
    #
    def deleteValue(self,val):
        name,vala = self.getAllNameValue()
        #
        for i in range(len(name)):
            if vala[i]==val:
                winreg.DeleteValue(self.reg_key, name[i])
#
class TraceThread(threading.Thread):
    def __init__(self, *args, **keywords):
        threading.Thread.__init__(self, *args, **keywords)
        self.killed = False
    def start(self):
        self._run = self.run
        self.run = self.settrace_and_run
        threading.Thread.start(self)
    def settrace_and_run(self):
        sys.settrace(self.globaltrace)
        self._run()
    def globaltrace(self, frame, event, arg):
        return self.localtrace if event == 'call' else None
    def localtrace(self, frame, event, arg):
        if self.killed and event == 'line':
            raise SystemExit()
        return self.localtrace
#
class MainGUI(tk.Frame):
    def __init__(self,master):
        tk.Frame.__init__(self, master=master)
        self.sw = self.master.winfo_screenwidth()
        self.sh = self.master.winfo_screenheight()
        w,h = 900,500 #
        self.master.geometry("{0}x{1}+{2}+{3}".format(w,h,int(self.sw/2-w/2),int(self.sh/2-h/2)))
        self.master.resizable(0,0)# fixed size
        self.master.wm_title("BKGA")
        #
        try:
            if os.path.isfile('bk.ico'):
                ico = 'bk.ico'
            
            self.master.wm_iconbitmap(ico)
        except:
            pass
        # registry
        self.reg = WIN_REGISTRY(path = "SOFTWARE\BKGA\RECLOSEROPTIM",keyUser="",nmax =1)
        self.initGUI()

   #
    def write(self, txt):
        self.text1.insert(tk.INSERT,txt)
        #
    def flush(self):
        pass
    #
    def close_buttons(self):
        self.master.destroy()
    #
    def clearConsol(self):
        self.text1.delete(1.0,tk.END)
    #
    def editInput(self):
        sfile = self.ipf_v.get()
        os.system('start "excel.exe" "{0}"'.format(sfile))
    #
    def initGUI(self):
        sys.stdout = self
        #
        fileFrame = tk.LabelFrame(self.master, text = "Files")
        fileFrame.grid(row=0, column=0, padx=10, pady=10)
        ipf = tk.Label(fileFrame, text="Input File : ")
        ipf.grid(row=0, column=0, sticky='E', padx=5, pady=5)

        self.ipf_v = tk.StringVar()
        self.opf_v = tk.StringVar()
        try:
            ft1 = self.reg.getAllValue()[0]
            if os.path.isfile(ft1):
                self.ipf_v.set(ft1)
        except:
            pass
        ipfTxt = tk.Entry(fileFrame,width= 103,textvariable=self.ipf_v)
        ipfTxt.grid(row=0, column=1, sticky="W",padx=5, pady=5)
        #
        ipf_b1 = tk.Button(fileFrame, text="...",width= 6,relief= tk.GROOVE,command=self.selectInputFile)
        ipf_b1.grid(row=0, column=2, sticky='W', padx=5, pady=5)
        ipf_b2 = tk.Button(fileFrame, text="Edit in Excel",width= 10,relief= tk.GROOVE,command=self.editInput)
        ipf_b2.grid(row=0, column=3, sticky='W', padx=5, pady=5)
        #
        opf = tk.Label(fileFrame, text="Output File : ")
        opf.grid(row=2, column=0, sticky='E', padx=5, pady=5)
        opfTxt = tk.Entry(fileFrame,width= 103,textvariable=self.opf_v)
        opfTxt.grid(row=2, column=1, sticky="W",padx=5, pady=5)
        opf_b1 = tk.Button(fileFrame, text="...",width= 6,relief= tk.GROOVE,command=self.selectOutputFile)
        opf_b1.grid(row=2, column=2, sticky='W', padx=5, pady=5)
        #
        csFrame = tk.LabelFrame(self.master, relief= tk.GROOVE, bd=0)
        #
        self.text1 = tk.Text(csFrame,wrap = tk.NONE,width=108,height=20)#,
        # yScroll
        yscroll = tk.Scrollbar(csFrame, orient=tk.VERTICAL, command=self.text1.yview)
        self.text1['yscroll'] = yscroll.set
        yscroll.pack(side=tk.RIGHT, fill=tk.Y)
        #
        # xScroll
        xscroll = tk.Scrollbar(csFrame, orient=tk.HORIZONTAL, command=self.text1.xview)
        self.text1['xscroll'] = xscroll.set
        xscroll.pack(side=tk.BOTTOM, fill=tk.X)
        #
        self.text1.pack(fill=tk.BOTH, expand=tk.Y)

        #
        btFrame = tk.Frame(self.master, relief= tk.GROOVE, bd=0)
        #
        self.run_b = tk.Button(btFrame, text="Launch",relief= tk.GROOVE,width= 10, command=self.run1)
        self.run_b.grid(row=0,column=0, padx=0, pady=5)
        #
        close_b = tk.Button(btFrame, text="Exit",width= 10,relief= tk.GROOVE, command=self.close_buttons)
        close_b.grid(row=0,column=1, padx=35, pady=5)
        #
        clearCs_b = tk.Button(btFrame, text="Clear console",width= 10,relief= tk.GROOVE, command=self.clearConsol)
        clearCs_b.grid(row=0,column=2, padx=0, pady=5)
        #
        self.pgFrame = tk.Frame(self.master, relief= tk.GROOVE, bd=0)
        self.var1 = tk.StringVar()
        self.percent =  ttk.Label(self.pgFrame , textvariable=self.var1,width= 18) #
        self.percent.grid(row=0,column=0, padx=10, pady=5)
        #
        self.progress = ttk.Progressbar(self.pgFrame,orient='horizontal',length=250,maximum=100, mode='determinate')
        self.progress.grid(row=0,column=1, padx=10, pady=5)
        #
        self.stop_b = tk.Button(self.pgFrame, text="Stop",relief= tk.GROOVE,width= 10, command=self.stop_progressbar)
        self.stop_b.grid(row=0,column=2, padx=10, pady=5)
        self.stop_b['state']='disabled'
        #
        #tab1_file.grid(row=, sticky='W', padx=10, pady=5, ipadx=5, ipady=5)

        csFrame.grid(row=3, column=0, padx=10, pady=5)
        btFrame.grid(row=4, column=0, padx=10, pady=2)
        #
        sys.stdout = self
        #color
        self.master.configure(bg='#BF2135')
        fileFrame.configure(bg="#F6E9E9", fg="#C72C41")
        ipf.configure(bg="#F6E9E9", fg="#C72C41")
        ipf_b1.configure(bg="#F6E9E9", fg="#C72C41")
        ipf_b2.configure(bg="#F6E9E9", fg="#C72C41")
        opf.configure(bg="#F6E9E9", fg="#C72C41")
        opf_b1.configure(bg="#F6E9E9", fg="#C72C41")
        self.text1.configure(bg="#F6E9E9", fg="#C72C41")
        self.run_b.configure(bg="#F6E9E9", fg="#C72C41")
        close_b.configure(bg="#F6E9E9", fg="#C72C41")
        clearCs_b.configure(bg="#F6E9E9", fg="#C72C41")
        self.stop_b.configure(bg="#F6E9E9", fg="#C72C41")
        btFrame.configure(bg='#C72C41')
        
    #   
    def selectInputFile(self):
        # Excell/CSV file
        v1 = tkf.askopenfilename(filetypes=[('Excel/csv file','*.xlsx *.xls .*csv'), ("All Files", "*.*")],title='Select Input file')
        if v1!='':
            self.ipf_v.set(v1)
            self.reg.appendValue(v1)
     #
    def selectOutputFile(self):
        v1 = tkf.asksaveasfile(defaultextension=".csv", filetypes=(("CSV file", "*.csv"),("All Files", "*.*") ))
        try:
            self.opf_v.set(v1.name)
        except:
            pass
    #
    def simulate(self):
        self.pgFrame.grid(row=3,column=0, padx=0, pady=0)
        self.stop_b['state']='disabled'
        self.var1.set("Reading data")
        try:
            self.stop_b['state']='active'
            self.var1.set("Running")
            self.progress['value'] = 0
            fi = self.ipf_v.get()
            for char_idx in range(len(fi)):
                if fi[char_idx] == '/':
                    found = char_idx
            fi = fi[found+1:]
            fo = self.opf_v.get()
            if fo:
                for char_idx in range(len(fo)):
                    if fo[char_idx] == '/':
                        found = char_idx
                fo = fo[found+1:]
                fo = 'res\\' + fo
            else:
                fo ='res\\resOptim12.csv'
            ga = new_optimum.GA(fi,fo)
            ga.run([],[])
            #os.system('start "excel.exe" "{0}"'.format(fo))
        except Exception as err:
            print(err)
        #
        self.finish()
    #
    def finish(self):
        self.run_b['state']='active'
        self.stop_b['state']='disabled'
        self.progress.stop()
        self.pgFrame.grid_forget()
    #
    def stop_progressbar(self):
        self.finish()
        self.t.killed = True
    #
    def run1(self):
        ft = self.ipf_v.get()
        if os.path.isfile(ft):
            self.reg.appendValue(ft)
        else:
            print('Input file not found:\n%s'%ft)
            return
        #
        try:
            self.progress.start()
            self.t = TraceThread(target=self.simulate)
            self.t.start()
            #
        except Exception as e:
            print(str(e))

        self.run_b['state']='disabled'
        self.stop_b['state']='active'
        #
def main():
    
    root = tk.Tk()
    feedback = MainGUI(root)
    root.mainloop()
#
if __name__ == '__main__':
    main()

