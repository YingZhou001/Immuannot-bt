import sys, re

#load id list

buf = ''
print_tag = False
for line in sys.stdin:
    if not line : break
    buf += line
    if line[0:2] == 'OS' :
        if 'human' in line : print_tag = True
    if '//' in line :
        if print_tag : sys.stdout.write(buf)
        buf = ''
        print_tag = False

