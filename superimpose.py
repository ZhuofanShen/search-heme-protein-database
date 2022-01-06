import os
from pymol import *
import time

pdb = 'TIPE3AF_clean'

# heme = '1YNR_truncated'
heme = '2YL7_truncated'

aa = 43

cmd.load(pdb + '.pdb')
cmd.load(heme + '.pdb')

cmd.pair_fit('/' + heme + '//A/1/CA', '/' + pdb + '//A/' + str(aa) + '/CA', \
    '/' + heme + '//A/4/CA', '/' + pdb + '//A/' + str(aa + 3) + '/CA', \
    '/' + heme + '//A/5/CA', '/' + pdb + '//A/' + str(aa + 4) + '/CA')
cmd.save(str(aa) + '_heme.pdb', heme)
time.sleep(1)
hec_lines = list()
with open(str(aa) + '_heme.pdb', 'r') as pf:
    for line in pf:
        if line.startswith('HETATM') and line[17:26] == 'HEX A 999':
            hec_lines.append(line)
os.remove(str(aa) + '_heme.pdb')

cmd.pair_fit('/' + heme + '//A/1/CA', '/' + pdb + '//A/' + str(aa) + '/CA', \
    '/' + heme + '//A/1/N', '/' + pdb + '//A/' + str(aa) + '/N', \
    '/' + heme + '//A/1/C', '/' + pdb + '//A/' + str(aa) + '/C')
cmd.pair_fit('/' + heme + '//A/1/CA', '/' + pdb + '//A/' + str(aa) + '/CA', \
    '/' + heme + '//A/1/CB', '/' + pdb + '//A/' + str(aa) + '/CB')
cmd.pair_fit('/' + heme + '//A/1/CA', '/' + pdb + '//A/' + str(aa) + '/CA')
cmd.save(str(aa) + '_cyx1.pdb', heme)
time.sleep(1)
cyx1_lines = list()
with open(str(aa) + '_cyx1.pdb', 'r') as pf:
    for line in pf:
        if line.startswith('ATOM') and line[17:26] == 'CYX A   1' and \
            line[13:15] != 'N ' and line[13:16] != 'CA ' and line[13:15] != 'C ' and line[13:15] != 'O ':
            cyx1_lines.append(line)
os.remove(str(aa) + '_cyx1.pdb')

cmd.pair_fit('/' + heme + '//A/4/CA', '/' + pdb + '//A/' + str(aa + 3) + '/CA', \
    '/' + heme + '//A/4/N', '/' + pdb + '//A/' + str(aa + 3) + '/N', \
    '/' + heme + '//A/4/C', '/' + pdb + '//A/' + str(aa + 3) + '/C')
cmd.pair_fit('/' + heme + '//A/4/CA', '/' + pdb + '//A/' + str(aa + 3) + '/CA', \
    '/' + heme + '//A/4/CB', '/' + pdb + '//A/' + str(aa + 3) + '/CB')
cmd.pair_fit('/' + heme + '//A/4/CA', '/' + pdb + '//A/' + str(aa + 3) + '/CA')
cmd.save(str(aa) + '_cyx4.pdb', heme)
time.sleep(1)
cyx4_lines = list()
with open(str(aa) + '_cyx4.pdb', 'r') as pf:
    for line in pf:
        if line.startswith('ATOM') and line[17:26] == 'CYX A   4' and \
            line[13:15] != 'N ' and line[13:16] != 'CA ' and line[13:15] != 'C ' and line[13:15] != 'O ':
            cyx4_lines.append(line)
os.remove(str(aa) + '_cyx4.pdb')

cmd.pair_fit('/' + heme + '//A/5/CA', '/' + pdb + '//A/' + str(aa + 4) + '/CA', \
    '/' + heme + '//A/5/N', '/' + pdb + '//A/' + str(aa + 4) + '/N', \
    '/' + heme + '//A/5/C', '/' + pdb + '//A/' + str(aa + 4) + '/C')
cmd.pair_fit('/' + heme + '//A/5/CA', '/' + pdb + '//A/' + str(aa + 4) + '/CA', \
    '/' + heme + '//A/5/CB', '/' + pdb + '//A/' + str(aa + 4) + '/CB')
cmd.pair_fit('/' + heme + '//A/5/CA', '/' + pdb + '//A/' + str(aa + 4) + '/CA')
cmd.save(str(aa) + '_his.pdb', heme)
time.sleep(1)
his_lines = list()
with open(str(aa) + '_his.pdb', 'r') as pf:
    for line in pf:
        if line.startswith('ATOM') and line[17:26] == 'HIS A   5' and \
            line[13:15] != 'N ' and line[13:16] != 'CA ' and line[13:15] != 'C ' and line[13:15] != 'O ':
            his_lines.append(line)
os.remove(str(aa) + '_his.pdb')

remark_lines = list()
pdb_lines = list()
cyx1_cb = False
cyx4_cb = False
his_cb = False
with open(pdb + '.pdb', 'r') as pf:
    for line in pf:
        if line.startswith('ATOM') or line.startswith('HETATM'):
            idx = int(line[22:26].strip(' '))
            if idx == aa:
                if line[13:15] == 'N ' or line[13:16] == 'CA ' or line[13:15] == 'C ' or line[13:15] == 'O ':
                    pdb_lines.append('ATOM  ' + line[6:17] + 'CYX' + line[20:])
            elif idx == aa + 1:
                if not cyx1_cb:
                    if aa < 10:
                        aa_str = '  ' + str(aa)
                    elif aa < 100:
                        aa_str = ' ' + str(aa)
                    elif aa < 1000:
                        aa_str = str(aa)
                    for cyx1_line in cyx1_lines:
                        pdb_lines.append(cyx1_line[:23] + aa_str + cyx1_line[26:])
                    remark_lines.append('REMARK 666 MATCH TEMPLATE A HEX  999 MATCH MOTIF A CYX  ' + aa_str + '  1  1               \n')
                    cyx1_cb = True
                pdb_lines.append(line)
            elif idx == aa + 3:
                if line[13:15] == 'N ' or line[13:16] == 'CA ' or line[13:15] == 'C ' or line[13:15] == 'O ':
                    pdb_lines.append('ATOM  ' + line[6:17] + 'CYX' + line[20:])
            elif idx == aa + 4:
                if not cyx4_cb:
                    if aa + 3 < 10:
                        aa_str = '  ' + str(aa + 3)
                    elif aa + 3 < 100:
                        aa_str = ' ' + str(aa + 3)
                    elif aa + 3 < 1000:
                        aa_str = str(aa + 3)
                    for cyx4_line in cyx4_lines:
                        pdb_lines.append(cyx4_line[:23] + aa_str + cyx4_line[26:])
                    remark_lines.append('REMARK 666 MATCH TEMPLATE A HEX  999 MATCH MOTIF A CYX  ' + aa_str + '  2  1               \n')
                    cyx4_cb = True
                if line[13:15] == 'N ' or line[13:16] == 'CA ' or line[13:15] == 'C ' or line[13:15] == 'O ':
                    pdb_lines.append('ATOM  ' + line[6:17] + 'HIS' + line[20:])
            elif idx == aa + 5:
                if not his_cb:
                    if aa + 4 < 10:
                        aa_str = '  ' + str(aa + 4)
                    elif aa + 4 < 100:
                        aa_str = ' ' + str(aa + 4)
                    elif aa + 4 < 1000:
                        aa_str = str(aa + 4)
                    for his_line in his_lines:
                        pdb_lines.append(his_line[:23] + aa_str + his_line[26:])
                    remark_lines.append('REMARK 666 MATCH TEMPLATE A HEX  999 MATCH MOTIF A HIS  ' + aa_str + '  3  1               \n')
                    his_cb = True
                pdb_lines.append(line)
            else:
                pdb_lines.append(line)
pdb_lines.append('TER\n')
for line in hec_lines:
    pdb_lines.append(line)
pdb_lines.append('END\n')
with open(pdb[:5] + 'X' + str(aa) + 'X' + str(aa + 3) + '.pdb', 'w') as pf:
    pf.writelines(remark_lines)
    pf.writelines(pdb_lines)
