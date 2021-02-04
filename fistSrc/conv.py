
files = "ck.fis", "dc.fis", "hl.p3fis", "sp.p3fis", "udm_hl.p3fis", "ct.fis", "ft.fis", "sp.p2fis", "tt.fis"
import re

operators = ":+.=-/*^&|<>#(),@\n;\"'{}[]"
tokens = set()
for f in files:
    print(f)
    for line in open(f, "r").readlines():

        m = re.match(r'^([^;]*);(.*)$', line)
        if m:  # The line contains a hash / comment
            line = m.group(1)
        for op in operators:
            line = line.replace(op, " ")
        token_list = line.split(" ")
        for t in token_list:
            tokens.add(t)

l = list(tokens)
l.sort()

with open("tokens.txt", "w") as f:
    for item in l:
        print("global ", item, file=f)
