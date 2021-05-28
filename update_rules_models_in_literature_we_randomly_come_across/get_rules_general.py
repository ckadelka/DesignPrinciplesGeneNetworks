rules = '''
SHR = SHR & ! CYA | (SHR | ! CYA) & PHB 
SCR = SYS | JKD
JKD = ! PHB & (MGP | SYS) | JKD & SYS & MGP
MGP = JKD & ! WOX5 (! SYS | ! MGP) 
miRNA165 = SYS | miRNA165
PHB = ! miRNA165 & ! CYA
Auxin = Auxin & ! PHB
IAA5 = ! Auxin
WOX5 = PYI5 & ! CYA & MGP | CYA & MGP
CLE = (CLE | ! SHR) & ! IAA5
ACR = CLE
SYS = SHR & SCR
CYA = CLE & ACR
PYI5 = ! PHB & ! IAA5
'''

rules = rules.replace('(', ' ( ').replace('&', 'AND').replace('|', 'OR').replace('!', 'NOT')

print(rules)