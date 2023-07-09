
import tkinter as tk
from tkinter import ttk
from tkinter import scrolledtext
from tkinter import Menu
from tkinter import messagebox as mBox
from tkinter import Spinbox
from tkinter import filedialog as fd
import Analysisprotein as ap
import bp
from collapsible import CollapsiblePane as cp

def on_configure(event):
    # update scrollregion after starting 'mainloop'
    # when all widgets are in canvas
    mCanvas.configure(scrollregion=mCanvas.bbox('all'))

# Create instance
win = tk.Tk()
menubar = Menu(win)
win.config(menu=menubar)
win.geometry("1000x850")

# Add a title
win.title("SBA GUI")

def a_env():
    #add new window for credential input
    cred_win = tk.Toplevel(win)
    cred_win.title("Credentials")
    cred_win.geometry("400x150")
    cred_win.resizable(False, False)
    #add labels and entries
    #email
    email_label = ttk.Label(cred_win, text="Email:")
    email_label.grid(column=0, row=0,sticky='NW', padx=10, pady=10)
    email_entry = ttk.Entry(cred_win, width=25)
    email_entry.grid(column=1, row=0,sticky='NW', padx=10, pady=10)
    #api key
    api_label = ttk.Label(cred_win, text="API key (optional):")
    api_label.grid(column=0, row=1,sticky='NW', padx=10, pady=10)
    api_entry = ttk.Entry(cred_win, width=25)
    api_entry.grid(column=1, row=1,sticky='NW', padx=10, pady=10)
    #button
    def save_cred():
        email = email_entry.get()
        api_key = api_entry.get()
        if email == '':
            mBox.showinfo('Error', 'Please enter a valid email')
        else:
            bp.save_credentials(email, api_key)
            cred_win.destroy()
    save_button = ttk.Button(cred_win, text="Save", command=save_cred )
    save_button.grid(column=1, row=2,sticky='NW', padx=10, pady=10)

    #exit button
    def exit_cred():
        cred_win.destroy()
    exit_button = ttk.Button(cred_win, text="Exit", command=exit_cred )
    exit_button.grid(column=0, row=2,sticky='NW', padx=10, pady=10)

#open credentials on startup if they don't exist
try:
    bp.get_credentials()
except:
    a_env()


cred_menu=Menu(menubar)
menubar.add_cascade(label ='Setup', menu = cred_menu)
cred_menu.add_command(
    label='Credentials',
    command=a_env
)


#mCanvas= tk.Canvas(win)
#mCanvas.pack(side=tk.LEFT)

#scroll_bar = tk.Scrollbar(win,command=mCanvas.yview)
#scroll_bar.pack(side=tk.RIGHT, fill='y')

#mCanvas.configure(yscrollcommand = scroll_bar.set)
#mCanvas.bind('<Configure>', on_configure)

mframe = tk.Frame(win)
#mCanvas.create_window((0,0), window=mframe)
mframe.pack()

### Query Entrez ###
fasta_frame = tk.Frame(mframe)
#collapsible pane
cEntrez = cp(fasta_frame, 'Close Entrez', 'Open Entrez')
cEntrez.toggle()
cEntrez.grid(row = 0, column = 0)

# Query Entrez Label
query_label = ttk.Label(cEntrez.frame, text="Search Entrez:")
query_label.grid(column=0, row=0,sticky='NW', padx=10, pady=10)

# Query Entrez Entry
st = tk.StringVar()
query_entry = ttk.Entry(cEntrez.frame, width=50, text=st)
st.set("Dehalococcoides mccartyi[Organism] AND methionine synthase[Protein]")
query_entry.grid(column=1, row=0,sticky='NW', padx=10, pady=10)

# entrez query function
global last_accession
last_accession = None
global last_record
last_record = None

def entrez_query():
    analyze_box.configure(state='normal')
    analyze_box.delete(1.0, tk.END)
    analyze_box.configure(state='disabled')

    fasta_box.configure(state='normal')
    fasta_box.delete(1.0, tk.END)
    fasta_box.configure(state='disabled')

    save_button.configure(state='disabled')
    analyze_button.configure(state='disabled')

    res_list.delete(*res_list.get_children())

    align_text.configure(state='normal')
    align_text.delete('1.0', tk.END)
    align_text.configure(state='disabled')

    user_input = query_entry.get()

    entrtez_accession = bp.get_accession_number(user_input)
    entrez_fasta = bp.download_protein_sequence(entrtez_accession)

    fasta_box.configure(state='normal')
    fasta_box.insert(tk.INSERT, entrez_fasta)
    fasta_box.configure(state='disabled')
    save_button.configure(state='normal')
    analyze_button.configure(state='normal')

    global last_accession
    last_accession = entrtez_accession
    global last_record
    last_record = entrez_fasta

#query button
query_button = ttk.Button(cEntrez.frame, text="Search", command=entrez_query)
query_button.grid(column=2, row=0,sticky='SW', padx=10, pady=10)

#text output label
fasta_label = ttk.Label(cEntrez.frame, text="Fasta sequence:")
fasta_label.grid(column=0, row=1,sticky='SW', padx=10)

#Text output box
fasta_scroll = ttk.Scrollbar(cEntrez.frame)
fasta_box = scrolledtext.ScrolledText(cEntrez.frame, yscrollcommand=fasta_scroll.set, height=10, wrap=tk.WORD, )
fasta_box.grid(column=0, row=2,sticky='NW', padx=10, pady=10, columnspan=3)
fasta_box.configure(state='disabled')

#analyze output box
analyze_box = scrolledtext.ScrolledText(cEntrez.frame, width=35, height=10, wrap=tk.WORD, )
analyze_box.grid(column=3, row=2,sticky='NW', padx=10, pady=10, columnspan=3)
analyze_box.configure(state='disabled')

#save function
def save_fasta():
    filetypes = (
        ('fasta file', '*.fasta'),
        ('All files', '*.*')
    )

    filename = fd.asksaveasfilename(
        title='Save file as',
        initialdir='./',
        filetypes=filetypes,
        defaultextension='.fasta',
        
    )

    if filename:
        bp.write_fasta(filename, last_record)

def analyze_prot():
    analyze_box.configure(state='normal')
    analyze_box.delete(1.0, tk.END)
    analyze_box.configure(state='disabled')
    read = bp.read_fasta(last_record)
    analysis = ap.analyze(read)
    analyze_box.configure(state='normal')
    analyze_box.insert(tk.INSERT, analysis)
    analyze_box.configure(state='disabled')


#analyze button
analyze_button = ttk.Button(cEntrez.frame, text="Analyze", command=analyze_prot)
analyze_button.grid(column=3, row=1,sticky='NW', padx=10, pady=10)
analyze_button.configure(state='disabled')
#save button
save_button = ttk.Button(cEntrez.frame, text="Save", command=save_fasta )
save_button.grid(column=2, row=1,sticky='NW', padx=10, pady=10)
save_button.configure(state='disabled')

fasta_frame.pack() #pack the frame

### Search BLAST ###
# user options
user_options = tk.Frame(mframe)
#collapsible pane
cBlast = cp(user_options, 'Close BLAST options', 'Open BLAST options')
cBlast.grid(row = 0, column = 0)

# Search Label
search_label = ttk.Label(cBlast.frame, text="Filter BLAST:")
search_label.grid(column=0, row=3,sticky='NW', padx=10, pady=10)

#search Entry
search_entry = ttk.Entry(cBlast.frame, width=30)
search_entry.grid(column=1, row=3,sticky='NW', padx=10, pady=10)

# dropdown
search_var = tk.StringVar()
search_combobox = ttk.Combobox(cBlast.frame, width=12, textvariable=search_var, state='readonly')
search_combobox['values'] = ('All Fields','Accession', 'Organism', 'Protein Name')
search_combobox.current(0)
search_combobox.grid(column=2, row=3,sticky='NW', padx=10, pady=10)


# Add a label
res_size = ttk.Label(cBlast.frame, text="Result count:")
res_size.grid(column=0, row=4,sticky='W', padx=10, pady=10)

# Add a spinbox widget
spin_size = Spinbox(cBlast.frame, from_=1, to=100, values=(5,10,20,30,40,50,75,100), width=5, bd=8)
spin_size.grid(column=1, row=4,sticky='NW', padx=10, pady=10)

#add a label
threshold_label = ttk.Label(cBlast.frame, text="Threshold:")
threshold_label.grid(column=0, row=5,sticky='W', padx=10, pady=10)

#add a float entry
v = tk.StringVar()
threshold_entry = ttk.Entry(cBlast.frame, width=5,text=v)
v.set("12")
threshold_entry.grid(column=1, row=5,sticky='NW', padx=10, pady=10)

# Button Click Event Function
global b_r
b_r = []
def click_me():

    align_text.configure(state='normal')
    align_text.delete('1.0', tk.END)
    align_text.configure(state='disabled')
    save_selected_button.configure(state='disabled')
    search_query = search_entry.get()
    search_type = '' if search_query == '' else '[' + search_var.get() + ']'
    search = (search_query + search_type)#.replace(' ', '%20')           # replace spaces with %20 for url
    search_threshold = threshold_entry.get()
    print(search_query+search_type)
    res_list.delete(*res_list.get_children())
    protein_sequence = bp.read_fasta(last_record)
    blast_record = bp.blast_protein_sequence(protein_sequence, hl=int(spin_size.get()),entrez=search, threshold=float(search_threshold))
    global b_r
    b_r = blast_record
    save_selected_button.configure(state='normal')

    for a in blast_record:
        
        res_list.insert('', 'end', values=(a[1], a[2], a[0], a[4], a[6], a[3], a[5]))




# search button
search_button = ttk.Button(user_options, text="BLAST",command=click_me)
search_button.grid(column=0, row=5,sticky='NW', padx=10, pady=10)

user_options.pack()



# results
results = tk.Frame(mframe)
results.pack(fill='x')

result_label = ttk.Label(results, text="BLAST result:")
result_label.pack(side='top', anchor=tk.W, padx=10, pady=10)

res_scroll = tk.Scrollbar(results)
res_scroll.pack(side='right', fill='y')
res_list = ttk.Treeview(results,yscrollcommand=res_scroll.set, height=10, selectmode='extended')
res_list.pack(fill='x')
res_scroll.config(command=res_list.yview)


# sort columns in treeview on click
def treeview_sort_column(tv, col, reverse):
    l = [(tv.set(k, col), k) for k in tv.get_children('')]
    l.sort(reverse=reverse)

    # rearrange items in sorted positions
    for index, (val, k) in enumerate(l):
        tv.move(k, '', index)

    # reverse sort next time
    tv.heading(col, text=col, command=lambda _col=col: \
                 treeview_sort_column(tv, _col, not reverse))
    
columns = ('protein', 'organism', 'acc', 'e', 'score', 'length', 'align_length')
res_list ['columns'] = columns
for col in columns:
    res_list.heading(col, text=col, command=lambda _col=col: \
                     treeview_sort_column(res_list, _col, False))

# format our column
res_list.column("#0", width=0,  stretch='no')
res_list.column("protein",anchor='center', width=80)
res_list.column("organism",anchor='center',width=80)
res_list.column("acc",anchor='center',width=80)
res_list.column("e",anchor='center',width=50)
res_list.column("score",anchor='center',width=50)
res_list.column("length",anchor='center',width=50)
res_list.column("align_length",anchor='center',width=110)

#Create Headings 
res_list.heading("#0",text="",anchor='center')
res_list.heading("protein",text="Protein",anchor='center')
res_list.heading("organism",text="Organism",anchor='center')
res_list.heading("acc",text="Accession no.",anchor='center')
res_list.heading("e",text="e",anchor='center')
res_list.heading("score",text="Score",anchor='center')
res_list.heading("length",text="Length",anchor='center')
res_list.heading("align_length",text="Alignment Length",anchor='center')



    
# save seleced rows to file function
def save_selected():
    # get selected rows
    selected = res_list.selection()
    records = []
    for entry in selected:
        curIndex = res_list.index(entry)
        add_record = b_r[curIndex]
        print(add_record)
        records.append(add_record)

    # save to file as fasta dialog
    filetypes = (
        ('fasta file', '*.fasta'),
        ('All files', '*.*')
    )

    filename = fd.asksaveasfilename(
        title='Save file as',
        initialdir='./',
        filetypes=filetypes,
        defaultextension='.fasta',
        
    )

    if filename:
        bp.write_fasta_multi(filename, records)


# save selected button
save_selected_button = ttk.Button(results, text="Save selected", command=save_selected)
save_selected_button.pack(side='top', anchor=tk.E, padx=10, pady=10)
save_selected_button.config(state='disabled')

# alignment
align_frame = tk.Frame(mframe)
align_frame.pack()

# radio buttons global and local
radio_frame = tk.Frame(align_frame)
radio_frame.pack(side='top', padx=10, pady=10)
global_local_var = tk.IntVar()


def plot_alignment():
    global plot_stats
    reference = last_accession
    ap.dotplot(plot_stats[0],plot_stats[1],mode=plot_stats[3],xlabel=last_accession,ylabel=plot_stats[2])
    
# plot alignment button
global plot_stats
plot_button = ttk.Button(align_frame, text="Plot alignment", command=plot_alignment)
plot_button.pack(side='left', padx=10, pady=10)
plot_button.config(state='disabled')

align_label = ttk.Label(radio_frame, text="Alignment:")
align_label.pack(side='top', padx=10, pady=10)

align_scroll = tk.Scrollbar(align_frame, orient='horizontal')
align_scroll.pack(side='bottom', fill='x')
align_text = tk.Text(align_frame, height=10, width=100, wrap='none', xscrollcommand=align_scroll.set, state='disabled',font='TkFixedFont')
align_text.pack(side='left', fill='x')

align_scroll.config(command=align_text.xview)

#event on selecting a row in the treeview
def selectItem(a=None):
    plot_button.config(state='disabled')
    curItem = res_list.focus()
    curIndex = res_list.index(curItem)
    curQuery = b_r[curIndex][10]
    currAcc = b_r[curIndex][0]
    curAlign = bp.read_fasta(bp.download_protein_sequence(currAcc)) 

    #annotation counter
    sep_no=25
    char_count = [str(i) for i in range(25, len(curAlign), sep_no)]
    newline_str =  20*" " + "1" + 24*" " #fml
    for i in char_count:
        newline_str += i + (sep_no-len(i))*" " 
    # call the alignment function
    g_l = 'global' if global_local_var.get() else 'local'
    alignment, score = bp.align_sequences(curQuery, curAlign, g_l)
    global plot_stats
    plot_stats = [curQuery, curAlign, currAcc, g_l]

    align_text.configure(state='normal')
    align_text.delete('1.0', tk.END)
    align_text.insert(tk.END, alignment  )
    align_text.insert(tk.END, newline_str)
    align_text.configure(state='disabled')
    plot_button.config(state='normal')
    

    

tk.Radiobutton(radio_frame, text="global", variable=global_local_var, value=1, command=selectItem).pack(side='left')
tk.Radiobutton(radio_frame, text="local", variable=global_local_var, value=0, command=selectItem).pack(side='left')

res_list.bind('<<TreeviewSelect>>', selectItem)
radio_frame.bind('<Button-1>', selectItem)

if __name__ == '__main__':
    # Start GUI
    win.mainloop()