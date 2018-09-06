import pandas as pd

kcal_j = 4.184
ang_to_nm = 0.1
inv_ang_square_to_inv_nm_square = 100

angle_names = ["n2-ce-c3", "n2-ce-ce", "n-n2-ce", "hn-n-n2", "c-n-n2", "ce-c3-c3"]
lmbd_amber = [122.64, 123.0, 118.09, 118.33, 120.59, 111.22]
k_lmbd_amber = [66.94, 68.52, 70.33, 49.62, 68.06, 63.65]
angle_df = pd.DataFrame(
    {
        "Angle": angle_names,
        "Lambda (Amber)": lmbd_amber,
        "k_lambda (Amber)": k_lmbd_amber,
    }
)
angle_df["Lambda (GROMACS)"] = angle_df["Lambda (Amber)"]
angle_df["k_lambda (GROMACS)"] = angle_df["k_lambda (Amber)"] * kcal_j

bond_names = ["n2-ce", "n2-n"]
b0_amber = [1.2790, 1.3710]
k_b_amber = [599.8, 499.7]
bond_df = pd.DataFrame(
    {"Bond": bond_names, "b0 (Amber)": b0_amber, "k_b (Amber)": k_b_amber}
)
bond_df["b0 (GROMACS)"] = bond_df["b0 (Amber)"] * ang_to_nm
bond_df["k_b (GROMACS)"] = (
    bond_df["k_b (Amber)"] * kcal_j * inv_ang_square_to_inv_nm_square
)

dihedral_names = ["n2-ce", "n2-n"]
theta_0 = [180.0, 0.0]
k_theta = [8.3, 0.8]
dihedral_df = pd.DataFrame(
    {"Dihedral": dihedral_names, "Theta (Amber)": theta_0, "k_theta (Amber)": k_theta}
)
dihedral_df["Theta (GROMACS)"] = dihedral_df["Theta (Amber)"]
dihedral_df["k_Theta (GROMACS)"] = dihedral_df["k_theta (Amber)"] * kcal_j

print(angle_df)
print(bond_df)
print(dihedral_df)
