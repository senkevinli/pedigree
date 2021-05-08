from graphviz import Graph
h = Graph('html_table', format='png')

h.node('tab', label='''<<TABLE>
 <TR>
   <TD>left</TD>
   <TD>right</TD>
 </TR>
</TABLE>>''')

h.render('woah')