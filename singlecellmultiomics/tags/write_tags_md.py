import tabulate
import singlecellmultiomics.tags

table = []
for tag in singlecellmultiomics.tags.tags:
    table.append( [tag.tag, tag.humanName, ('No','Yes')[tag.isPhred], ('No','Yes')[tag.doNotWrite]] )

headers = ['Tag','Description','Is phred','Written by default']
with open('TAGS.MD','wt') as f:
    f.write( tabulate.tabulate(table, headers, tablefmt='github') )
