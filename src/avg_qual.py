import sys

num = 0
dom = 0
for line in sys.stdin:
    for char in line:
        num += 10 ** (- ord(char)/ 10)
        dom += 1
print(num / dom)
