import argparse
import os
from pymol import *
import shutil
import time
import traceback

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input_files', type=str, default='pdb')
parser.add_argument('-o', '--processed_files', type=str, default='processed_pdb')
args = parser.parse_args()

for output_dir in ['HEM', 'HEM_symm', 'HEM_asymm', 'HEM_noREMARK350', \
        'HEA', 'HEA_symm', 'HEA_asymm', 'HEA_noREMARK350', \
        'HEB', 'HEB_symm', 'HEB_asymm', 'HEB_noREMARK350', \
        'HEC', 'HEC_symm', 'HEC_asymm', 'HEC_noREMARK350', \
        'ABAB', 'error', 'no_REMARK620', 'invalid_select', 'inequivalent_hemes', args.processed_files]:
    if not os.path.isdir(output_dir):
        os.mkdir(output_dir)
for pdb in filter(lambda x: x.endswith('.pdb'), os.listdir(args.input_files)):
    try:
        #print('\n' + pdb[:-4])
        log = list()
        # pymol load pdb file
        cmd.load(args.input_files + '/' + pdb)
        # Read the PDB file
        with open(args.input_files + '/' + pdb, 'r') as ppdb:
            e_coli = False
            total_biomolecules = 0
            assembly = None
            consistent_assembly = True
            main_chains = None
            main_heme_chain = None
            no_remark350 = False # if not main_chains: no_remark350 = True
            hemes = dict()
            standard_heme = None
            inequivalent_hemes = False
            in_heme_remark_620 = False
            proximal_residues = dict()
            proximal_ligands = dict()
            proximal_atom_heme_dihedral = None
            link_lines = list()
            ssbond_lines = list()
            for line in ppdb:
                if line.startswith('SOURCE') and line[11:46] == 'EXPRESSION_SYSTEM: ESCHERICHIA COLI':
                    e_coli = True
                elif line.startswith('REMARK 300 BIOMOLECULE: '):
                    total_biomolecules = len(line[24:-1].strip(' ').split(', '))
                elif line.startswith('REMARK 350 AUTHOR DETERMINED BIOLOGICAL UNIT: ') or \
                        line.startswith('REMARK 350 SOFTWARE DETERMINED QUATERNARY STRUCTURE: '):
                    current_assembly = line.strip('\n').strip(' ').split(' ')[-1]
                    if not assembly:
                        assembly = current_assembly
                    elif assembly != current_assembly:
                        consistent_assembly = False
                    total_biomolecules -= 1
                elif line.startswith('REMARK 350 APPLY THE FOLLOWING TO CHAINS: '):
                    if not assembly:
                        assembly = 'VOID'
                        current_assembly = 'VOID'
                    if not main_chains:
                        main_chains = line[42:-1].strip(' ').split(', ')
                        #print('main chains: ' + ' '.join(main_chains))
                        log.append('main chains: ' + ' '.join(main_chains) + '\n')
                    else:
                        current_main_chains = line[42:-1].strip(' ').split(', ')
                        if len(current_main_chains) > len(main_chains) and \
                                not set(main_chains).issubset(set(current_main_chains)):
                            assembly = current_assembly
                            main_chains = current_main_chains
                            #print('updated main chains: ' + ' '.join(main_chains))
                            log.append('updated main chains: ' + ' '.join(main_chains) + '\n')
                elif line.startswith('REMARK 620'):
                    if not main_chains:
                        no_remark350 = True
                        main_chains = list() # = ['A']
                        #print('Automatically generate main chains.')
                        log.append('Automatically generate main chains.\n')
                    if (line[43] in main_chains or no_remark350) and \
                            line.startswith('REMARK 620                           ') and \
                            line[39:42] in ['HEM', 'HEA', 'HEB', 'HEC']:
                        current_heme = (line[39:42], line[43], line[44:48].strip(' '))
                        if not main_heme_chain:
                            main_heme_chain = current_heme[1]
                        if no_remark350: # and current_heme[1] not in main_chains:
                            main_chains.append(current_heme[1])
                        #print('heme: ' + str(current_heme))
                        log.append(str(current_heme) + '\n')
                        if current_heme[1] not in hemes.keys():
                            in_heme_remark_620 = True
                            hemes[current_heme[1]] = current_heme
                            if not standard_heme:
                                standard_heme = current_heme
                            elif current_heme[0] != standard_heme[0] or current_heme[2] != standard_heme[2]:
                                inequivalent_hemes = True
                                break
                        elif line[44:48].strip(' ') != hemes[current_heme[1]][2]:
                            inequivalent_hemes = True
                            break
                    elif in_heme_remark_620 and not line.startswith('REMARK 620 N') and \
                            (line[17] == main_heme_chain and line[13:16] != hemes[line[17]][0] or line[17] != main_heme_chain) and \
                            line[13:16] + line[17] + line[18:22].strip(' ') not in proximal_residues and \
                            line[13:16] + line[17] + line[18:22].strip(' ') not in proximal_ligands:
                        # only record each proximal residue or ligand once for one heme
                        if line[13:16] in ['ALA', 'CYS', 'ASP', 'GLU', 'PHE', 'GLY', 'HIS', 'ILE', 'LYS', 'LEU', \
                                'MET', 'ASN', 'PRO', 'GLN', 'ARG', 'SER', 'THR', 'VAL', 'TRP', 'TYR']:
                            #print(line[13:16] + ': cmd.get_dihedral("/' + pdb[:-4] + '//' + current_heme[1] + '/' + current_heme[2] + '/NA", "/' + \
                                    #pdb[:-4] + '//' + line[17] + '/' + line[18:22].strip(' ') + '/' + line[24:28].strip(' ') + '", "/' + \
                                    #pdb[:-4] + '//' + current_heme[1] + '/' + current_heme[2] + '/FE", "/' + pdb[:-4] + '//' + current_heme[1] + '/' + current_heme[2] + '/NB")')
                            log.append(line[13:16] + ': cmd.get_dihedral("/' + pdb[:-4] + '//' + current_heme[1] + '/' + current_heme[2] + '/NA", "/' + \
                                    pdb[:-4] + '//' + line[17] + '/' + line[18:22].strip(' ') + '/' + line[24:28].strip(' ') + '", "/' + \
                                    pdb[:-4] + '//' + current_heme[1] + '/' + current_heme[2] + '/FE", "/' + pdb[:-4] + '//' + current_heme[1] + '/' + current_heme[2] + '/NB")' + '\n')
                            try:
                                proximal_atom_heme_dihedral = cmd.get_dihedral(\
                                        '/' + pdb[:-4] + '//' + current_heme[1] + '/' + current_heme[2] + '/NA', \
                                        '/' + pdb[:-4] + '//' + line[17] + '/' + line[18:22].strip(' ') + '/' + line[24:28].strip(' '), \
                                        '/' + pdb[:-4] + '//' + current_heme[1] + '/' + current_heme[2] + '/FE', \
                                        '/' + pdb[:-4] + '//' + current_heme[1] + '/' + current_heme[2] + '/NB')
                            except Exception as ex:
                                # traceback.print_exc()
                                log.append(traceback.format_exc())
                                proximal_atom_heme_dihedral = 'invalid'
                                break
                            if proximal_atom_heme_dihedral > 0:
                                position = 'BELOW'
                            else:
                                position = 'ABOVE'
                            proximal_residues[line[13:16] + line[17] + line[18:22].strip(' ')] = \
                                    (line[13:16], line[17], line[18:22].strip(' '), line[24:28].strip(' '), position)
                        else:
                            proximal_ligands[line[13:16] + line[17] + line[18:22].strip(' ')] = \
                                    (line[13:16], line[17], line[18:22].strip(' '), line[24:28].strip(' '), 'VOID')
                            # cmd.remove('resi ' + line[18:22].strip(' ') + ' & resn ' + line[13:16])
                    elif line.startswith('REMARK 620 N                    '):
                        in_heme_remark_620 = False
                elif line.startswith('LINK'): # and line[17:20] == heme[0] and not line.startswith('LINK        FE')
                    link_lines.append(line)
                elif line.startswith('SSBOND'):
                    ssbond_lines.append(line)
                elif line.startswith('ATOM'):
                    break
        # Exception part 1
        if proximal_atom_heme_dihedral == 'invalid':
            #print('Invalid dihedral atom select.')
            log.append('Invalid dihedral atom select.\n')
            with open('invalid_select/' + pdb[:-4] + '.log', 'w') as plog:
                plog.writelines(log)
            os.rename(args.input_files + '/' + pdb, 'invalid_select/' + pdb)
            continue
        elif inequivalent_hemes:
            #print('Inequivalent hemes exist in a single chain or multiple asymmetric chains.')
            log.append('Inequivalent hemes exist in a single chain or multiple asymmetric chains.\n')
            with open('inequivalent_hemes/' + pdb[:-4] + '.log', 'w') as plog:
                plog.writelines(log)
            os.rename(args.input_files + '/' + pdb, 'inequivalent_hemes/' + pdb)
            continue
        elif len(hemes) == 0:
            #print('No REMARK 620 lines for the main chain.')
            log.append('No REMARK 620 lines for the main chain.\n')
            with open('no_REMARK620/' + pdb[:-4] + '.log', 'w') as plog:
                plog.writelines(log)
            os.rename(args.input_files + '/' + pdb, 'no_REMARK620/' + pdb)
            continue
        # Otherwise we will have the "proximal_residues" and "proximal_ligands" information
        # Decide output directory.
        main_heme = hemes[main_heme_chain]
        if assembly == 'MONOMERIC': # if remark300, there must be remark350
            symmetry = 'SYMMETRIC'
            outdir = main_heme[0] + '/' + pdb[:-4]
            #print('SYMMETRIC')
        else:
            if no_remark350:
                symmetry = 'VOID'
                outdir = main_heme[0] + '_noREMARK350/' + pdb[:-4]
            elif len(main_chains) == len(hemes):
                outdir = main_heme[0] + '_symm/' + pdb[:-4]
                symmetry = 'SYMMETRIC'
                #print('SYMMETRIC')
            elif len(hemes) == 1:
                outdir = main_heme[0] + '_asymm/' + pdb[:-4]
                symmetry = 'ASYMMETRIC'
                #print('ASYMMETRIC')
            else: # Exception part 2
                #print('ABAB asymmetric hemoprotein.')
                log.append('ABAB asymmetric hemoprotein.\n')
                with open('ABAB/' + pdb[:-4] + '.log', 'w') as plog:
                    plog.writelines(log)
                os.rename(args.input_files + '/' + pdb, 'ABAB/' + pdb)
                continue
        os.mkdir(outdir)
        # check if proximal_residues
        if len(proximal_residues) == 0:
            with open(outdir + '/no_proximal_res', 'w') as _:
                pass
        # else:
            #print('proximal_residues: \n' + str(proximal_residues.values()))
        # E. coli
        if e_coli:
            with open(outdir + '/e_coli', 'w') as _:
                pass
        # Exception part 3
        if no_remark350: # if no REMARK 620, no_remark350 is always False. But the no REMARK 620 case is handled above.
            assembly = 'VOID'
            #print('No REMARK 350 main chains information.')
            log.append('No REMARK 350 main chains information.\n')
            with open(outdir + '/' + pdb[:-4] + '.log', 'w') as plog:
                plog.writelines(log)
        elif not consistent_assembly:
            #print(pdb[:-4] + ': Assembly information is not consistant among different molecules.')
            log.append(pdb[:-4] + ': Assembly information is not consistant among different molecules.\n')
            with open(outdir + '/' + pdb[:-4] + '.log', 'w') as plog:
                plog.writelines(log)
        # Write .out file and _clean.pdb file
        shutil.copy(args.input_files + '/' + pdb, outdir + '/' + pdb)
        with open(outdir + '/' + pdb[:-4] + '.out', 'w') as pout:
            pout.write(symmetry + ' ' + assembly + '\n')
            pout.write(' '.join(main_chains) + '\n')
            pout.write(' '.join(main_heme) + '\n')
            for proximal_residue in proximal_residues.values():
                pout.write(' '.join(proximal_residue) + '\n')
            for proximal_ligand in proximal_ligands.values():
                pout.write(' '.join(proximal_ligand) + '\n')
        main_chains_selection = ' | '.join(['chain ' + chain for chain in main_chains])
        cmd.create('clean', '(' + main_chains_selection + ') & (polymer.protein | resn ' + main_heme[0] + ')')
        cmd.save(outdir + '/' + pdb[:-4] + '_clean.pdb', 'clean')
        if assembly != 'MONOMERIC':
            cmd.delete('clean')
            cmd.create('mono', '(chain ' + main_heme_chain + ') & (polymer.protein | resn ' + main_heme[0] + ')')
            cmd.save(outdir + '/' + pdb[:-4] + '_mono.pdb', 'mono')
        cmd.delete('*')
        # write link and ssbond lines
        if len(link_lines) > 0 or len(ssbond_lines) > 0:
            time.sleep(5)
            with open(outdir + '/' + pdb[:-4] + '_clean.pdb', 'r+') as p_clean:
                clean_lines = p_clean.readlines()
                p_clean.seek(0)
                if len(ssbond_lines) > 0:
                    p_clean.writelines(ssbond_lines)
                if len(link_lines) > 0:
                    p_clean.writelines(list(filter(lambda x: x[21] in main_chains, link_lines)))
                p_clean.writelines(clean_lines)
                p_clean.truncate()
            if assembly != 'MONOMERIC':
                with open(outdir + '/' + pdb[:-4] + '_mono.pdb', 'r+') as p_mono:
                    mono_lines = p_mono.readlines()
                    p_mono.seek(0)
                    if len(ssbond_lines) > 0:
                        p_mono.writelines(ssbond_lines)
                    if len(link_lines) > 0:
                        p_mono.writelines(list(filter(lambda x: x[21] == main_heme_chain, link_lines)))
                    p_mono.writelines(mono_lines)
        os.rename(args.input_files + '/' + pdb, args.processed_files + '/' + pdb)
    except Exception as ex:
        # traceback.print_exc()
        shutil.rmtree(outdir)
        with open('error/' + pdb[:-4] + '.log', 'w') as plog:
            plog.writelines(log)
        with open('error/' + pdb[:-4] + '.err', 'w') as perr:
            perr.write(traceback.format_exc())
        os.rename(args.input_files + '/' + pdb, 'error/' + pdb)
