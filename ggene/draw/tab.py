


from tabulate import tabulate

from ggene.draw.colors import Colors





def plot_table(plot_tab, title = "", col_labels = [], row_labels = []):
    
    tab = []
    
    for i, plot_row in enumerate(plot_tab):
        
        row = []
        if row_labels:
            row.append(row_labels[i])
        
        for plot in plot_row:
            
            lines = plot.get_rows()
            row.append("\n".join(lines))
        tab.append(row)
    
    if title:
        print(f"{Colors.BOLD}{title}{Colors.RESET}")
    print(tabulate(tab, headers = col_labels))
    print()



