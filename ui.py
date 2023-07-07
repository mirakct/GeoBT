
import tkinter as tk
from tkinter import ttk
from tkinter import scrolledtext
from tkinter import Menu
from tkinter import messagebox as mBox
from tkinter import Spinbox
from tkinter import filedialog as fd
import Analysisprotein as ap
import bp

# Create instance
win = tk.Tk()
win.geometry("1000x850")

# Add a title
win.title("SBA GUI")



### Query Entrez ###
fasta_frame = tk.Frame(win)

# Query Entrez Label
query_label = ttk.Label(fasta_frame, text="Search Entrez:")
query_label.grid(column=0, row=0,sticky='NW', padx=10, pady=10)

# Query Entrez Entry
st = tk.StringVar()
query_entry = ttk.Entry(fasta_frame, width=50, text=st)
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
query_button = ttk.Button(fasta_frame, text="Search", command=entrez_query)
query_button.grid(column=2, row=0,sticky='SW', padx=10, pady=10)

#text output label
fasta_label = ttk.Label(fasta_frame, text="Fasta sequence:")
fasta_label.grid(column=0, row=1,sticky='SW', padx=10)

#Text output box
fasta_scroll = ttk.Scrollbar(fasta_frame)
fasta_box = scrolledtext.ScrolledText(fasta_frame, yscrollcommand=fasta_scroll.set, height=10, wrap=tk.WORD, )
fasta_box.grid(column=0, row=2,sticky='NW', padx=10, pady=10, columnspan=3)
fasta_box.configure(state='disabled')

#analyze output box
analyze_box = scrolledtext.ScrolledText(fasta_frame, width=35, height=10, wrap=tk.WORD, )
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
analyze_button = ttk.Button(fasta_frame, text="Analyze", command=analyze_prot)
analyze_button.grid(column=3, row=1,sticky='NW', padx=10, pady=10)
analyze_button.configure(state='disabled')
#save button
save_button = ttk.Button(fasta_frame, text="Save", command=save_fasta )
save_button.grid(column=2, row=1,sticky='NW', padx=10, pady=10)
save_button.configure(state='disabled')

fasta_frame.pack() #pack the frame

### Search BLAST ###
# user options
user_options = tk.Frame(win)
# Search Label
search_label = ttk.Label(user_options, text="Filter BLAST:")
search_label.grid(column=0, row=3,sticky='NW', padx=10, pady=10)

#search Entry
search_entry = ttk.Entry(user_options, width=30)
search_entry.grid(column=1, row=3,sticky='NW', padx=10, pady=10)

# dropdown
search_var = tk.StringVar()
search_combobox = ttk.Combobox(user_options, width=12, textvariable=search_var, state='readonly')
search_combobox['values'] = ('All Fields','Accession', 'Organism', 'Protein Name')
search_combobox.current(0)
search_combobox.grid(column=2, row=3,sticky='NW', padx=10, pady=10)


# Add a label
res_size = ttk.Label(user_options, text="Result count:")
res_size.grid(column=0, row=4,sticky='W', padx=10, pady=10)

# Add a spinbox widget
spin_size = Spinbox(user_options, from_=1, to=100, values=(5,10,20,30,40,50,75,100), width=5, bd=8)
spin_size.grid(column=1, row=4,sticky='NW', padx=10, pady=10)

#add a label
threshold_label = ttk.Label(user_options, text="Threshold:")
threshold_label.grid(column=0, row=5,sticky='W', padx=10, pady=10)

#add a float entry
v = tk.StringVar()
threshold_entry = ttk.Entry(user_options, width=5,text=v)
v.set("12")
threshold_entry.grid(column=1, row=5,sticky='NW', padx=10, pady=10)

# Button Click Event Function
global b_r
b_r = []
def click_me():

    align_text.configure(state='normal')
    align_text.delete('1.0', tk.END)
    align_text.configure(state='disabled')

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
    

    for a in blast_record:
        
        res_list.insert('', 'end', values=(a[1], a[2], a[0], a[4], a[6], a[3], a[5]))




# search button
search_button = ttk.Button(user_options, text="BLAST",command=click_me)
search_button.grid(column=2, row=5,sticky='NW', padx=10, pady=10)

user_options.pack()



# results
results = tk.Frame(win)
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



# alignment
align_frame = tk.Frame(win)
align_frame.pack()

# radio buttons global and local
radio_frame = tk.Frame(align_frame)
radio_frame.pack(side='top', padx=10, pady=10)
global_local_var = tk.IntVar()
tk.Radiobutton(radio_frame, text="global", variable=global_local_var, value=1).pack(side='left')
tk.Radiobutton(radio_frame, text="local", variable=global_local_var, value=0).pack(side='left')




align_label = ttk.Label(radio_frame, text="Alignment:")
align_label.pack(side='top', padx=10, pady=10)

align_scroll = tk.Scrollbar(align_frame, orient='horizontal')
align_scroll.pack(side='bottom', fill='x')
align_text = tk.Text(align_frame, height=10, width=100, wrap='none', xscrollcommand=align_scroll.set, state='disabled',font='TkFixedFont')
align_text.pack(side='left', fill='x')

align_scroll.config(command=align_text.xview)

#event on selecting a row in the treeview
def selectItem(a):
    curItem = res_list.focus()
    curIndex = res_list.index(curItem)
    curQuery = b_r[curIndex][10]
    curAlign = b_r[curIndex][12]
    #annotation counter
    sep_no=25
    char_count = [str(i) for i in range(25, len(curAlign), sep_no)]
    newline_str =  20*" " + "1" + 24*" " #fml
    for i in char_count:
        newline_str += i + (sep_no-len(i))*" " 
    # call the alignment function
    g_l = 'global' if global_local_var.get() else 'local'
    alignment, score = bp.align_sequences(curQuery, curAlign, g_l)

    align_text.configure(state='normal')
    align_text.delete('1.0', tk.END)
    align_text.insert(tk.END, alignment  )
    align_text.insert(tk.END, newline_str)
    align_text.configure(state='disabled')
    
res_list.bind('<<TreeviewSelect>>', selectItem)
radio_frame.bind('<Button-1>', selectItem)

if __name__ == '__main__':
    # Start GUI
    win.mainloop()