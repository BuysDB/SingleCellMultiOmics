import tabulate
import singlecellmultiomics.tags

if __name__ == '__main__':
    table = []
    for tag in singlecellmultiomics.tags.tags:
        table.append(
            (tag.tag, tag.humanName, ('No', 'Yes')[
                tag.isPhred], ('No', 'Yes')[
                tag.doNotWrite]))

    # Order alphabetically
    table = sorted(table)
    headers = ['Tag', 'Description', 'Is phred', 'Written by default']
    with open('TAGS.MD', 'wt') as f:
        f.write(tabulate.tabulate(table, headers, tablefmt='github'))
