import re
with open('output_sequence_.txt', 'r') as f:
    for line in f:
        matchObj = re.match(r'(.*)\[(.*?)\].*', line, re.M|re.I)
        print(matchObj.group(2))
        next(f)
