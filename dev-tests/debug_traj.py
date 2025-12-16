import mdtraj as md

traj_file = '/Users/fengyu/interaction_map_test/E1A_pat-summary/BB/S_0/__traj_0.xtc'
pdb_file = '/Users/fengyu/interaction_map_test/E1A_pat-summary/BB/S_0/__START_0.pdb'

print("Loading trajectory...")
t = md.load(traj_file, top=pdb_file)
print(f'Residues: {t.n_residues}')
print(f'Frames: {t.n_frames}')

print("\nSelecting pairs...")
indices = t.top.select_pairs('all', 'all')
print(f'Pairs: {len(indices)}')
print(f'Max atom index: {indices.max()}')

print("\nComputing contacts...")
try:
    distances, pairs = md.compute_contacts(t, indices)
    print('✓ Success!')
    print(f'Distances shape: {distances.shape}')
    print(f'Pairs shape: {pairs.shape}')
    print(f'Pairs max: {pairs.max()}')
except Exception as e:
    print(f'✗ Error: {e}')
