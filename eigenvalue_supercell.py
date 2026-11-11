import numpy as np
from scipy.sparse import lil_matrix


def get_eigenvalues_from_hoppings(k_points, hoppings_file, num_orbitals=16):

    data = np.loadtxt(hoppings_file)
    R_vecs = data[:, 0:3]
    orb_i = data[:, 3].astype(int)
    orb_j = data[:, 4].astype(int)
    t_vals = data[:, 5] + 1j * data[:, 6]
    
    all_eigvals = []

    for k in k_points:
        H_k = np.zeros((num_orbitals, num_orbitals), dtype=complex)
        
        # Convention 2 Phase Factor (from tbmodels documentation)
        phases = np.exp(2j * np.pi * np.dot(R_vecs, k))
        
        # Sum the reduced half-set of hoppings
        for idx in range(len(t_vals)):
            i, j = orb_i[idx], orb_j[idx]
            H_k[i, j] += t_vals[idx] * phases[idx]
            
       
        H_final = H_k + H_k.conj().T
        
        # Diagonalize
        evals = np.linalg.eigvalsh(H_final)
        all_eigvals.append(evals)
        
    return np.array(all_eigvals)

def build_supercell_matrix_fixed(hoppings_file, nx, ny, num_orbitals=16):
   
    data = np.loadtxt(hoppings_file)
    
    # Map R-vectors to hopping matrices
    hop_map = {}
    for row in data:
        R = tuple(row[0:3].astype(int))
        if R not in hop_map: 
            hop_map[R] = np.zeros((num_orbitals, num_orbitals), dtype=complex)
        hop_map[R][int(row[3]), int(row[4])] = row[5] + 1j * row[6]
    total_cells = nx * ny
    dim = total_cells * num_orbitals

    # Use lil_matrix for efficient construction
    H_super = lil_matrix((dim, dim), dtype=complex)

    # Create coordinates for all cells in the supercell
    cells = [(ix, iy, 0) for ix in range(nx) for iy in range(ny)]

    # Tile with Modulo and Addition
    for i, C1 in enumerate(cells):

        for R, matrix in hop_map.items():
            
            # This 'folds' the hoppings into the supercell dimensions
            target_ix = (C1[0] + R[0]) % nx
        
            target_iy = (C1[1] + R[1]) % ny
            target_iz = 0 
            # Find the index of that target cell
            j = cells.index((target_ix, target_iy, target_iz))
            # Define the block position in the big matrix
            row_s, col_s = i * num_orbitals, j * num_orbitals
            
            H_super[row_s:row_s+num_orbitals, col_s:col_s+num_orbitals] += matrix
            
    H_final = H_super + H_super.conj().T

    return H_final.tocsr()

