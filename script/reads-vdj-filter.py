import sys,re

print_first_line_tag = 0

for line in sys.stdin:
    if not line: break
    llst = line.split('\t')
    if print_first_line_tag == 0 : 
        sys.stdout.write(line)
        print_first_line_tag = 1
        continue
    complete_vdj = llst[20].replace('\n', '')
    if complete_vdj == 'T' :
        sys.stdout.write(line)
