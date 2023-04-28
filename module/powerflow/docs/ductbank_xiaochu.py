## ductbank calculation, Xiaochu Wang, 4/19/2023
import numpy as np

# thermal resistance: degree-cm/W
U = {
    "INSULATION": 550,
    "JACKET": 550,
    "DUCT": 480,
    "EARTH": 120,
    "CONCRETE": 85,
}


class DuctBank:
    def __init__(self, **kwargs):
        self.Ta = self.getarg(kwargs, "ambient_temperature")
        self.I_given = self.getarg(kwargs, "current")
        self.R_dc = self.getarg(kwargs, "conductor_resistance")
        self.Dc = self.getarg(kwargs, "segmental_conductor_diameter")
        self.Di = self.getarg(kwargs, "insulation_outer_diameter")
        self.Ds = self.getarg(kwargs, "sheath_outer_diameter")
        self.nn = self.getarg(kwargs, "conductors_per_cable")
        self.N = self.getarg(kwargs, "number_of_cables")
        self.id_duct = self.getarg(kwargs, "duct_inner_diameter")
        self.t_duct = self.getarg(kwargs, "duct_thickness")
        self.od_duct = self.id_duct + 2 * self.t_duct
        self.row = self.getarg(kwargs, "duct_bank_rows")
        self.col = self.getarg(kwargs, "duct_bank_columns")
        self.L = self.getarg(kwargs, "depth_of_reference_cable")
        self.d_vert = self.getarg(kwargs, "vertical_distance_between_duct_banks")
        self.d_hori = self.getarg(kwargs, "horizontal_distance_between_duct_banks")
        self.x = self.getarg(kwargs, "short_side_of_duct")
        self.y = self.getarg(kwargs, "long_side_of_duct")
        self.Lb = self.getarg(kwargs, "depth_to_center_of_duct_bank")
        self.LF = self.getarg(kwargs, "loss_factor")
        self.Dx = self.getarg(kwargs, "fictitious_diameter", 8.3)

        # Calculate thermal resistances
        self.calculate_thermal_resistances()

        # Calculate mutual heating factor
        self.calculate_mutual_heating_factor()

        # Calculate geometric factor
        self.calculate_geometric_factor()

        # Calculate total thermal resistance
        self.calculate_total_thermal_resistance()

        # Calculate cable core temperature
        self.calculate_cable_core_temperature()

    def getarg(self, args, name, default=None, required=True):
        value = args[name] if name in args else default
        if value is None and required:
            raise ValueError(f"{name} is required")
        setattr(self, name, value)
        return value

    def calculate_thermal_resistances(self):
        self.Ri = 0.012 * U["INSULATION"] * np.log10(self.Di / self.Dc)
        self.Rsd = self.nn * 4.6 / (self.Ds + 0.27)
        self.Rd = 0.0104 * U["DUCT"] * self.nn * self.t_duct / (self.od_duct - self.t_duct)

    def calculate_mutual_heating_factor(self):
        N = self.row * self.col
        d_diag = (self.d_vert ** 2 + self.d_hori ** 2) ** 0.5
        d_vec = np.array([1, self.d_vert, self.d_vert, self.d_hori, d_diag, d_diag])
        d_diag_short = (self.d_hori ** 2 + (2 * self.L - self.d_vert) ** 2) ** 0.5
        d_diag_mid = (self.d_hori ** 2 + (2 * self.L) ** 2) ** 0.5
        d_diag_long = (self.d_hori ** 2 + (2 * self.L + self.d_vert) ** 2) ** 0.5
        D_vec = np.array(
            [1, 2 * self.L - self.d_vert, 2 * self.L + self.d_vert, d_diag_short, d_diag_mid, d_diag_long]
        )
        product1 = np.prod(d_vec)
        product2 = np.prod(D_vec)
        self.F = round(product2 / product1)

    def calculate_geometric_factor(self):
        rb_c1 = self.x / self.y / 2 * (4 / np.pi - self.x / self.y)
        rb_c2 = np.log10(1 + self.y ** 2 / self.x ** 2)
        rb_c3 = np.log10(self.x / 2)
        rb = 10 ** (rb_c1 * rb_c2 + rb_c3)
        u = self.Lb / rb
        self.Gb = np.log10(u + (u ** 2 - 1) ** 0.5)
        self.Gb = round(self.Gb, 2)

    def calculate_total_thermal_resistance(self):
        Re_1 = np.log10(self.Dx / self.od_duct)
        Re_2 = self.LF * np.log10(4 * self.L / self.Dx * self.F)
        Re_corr = 0.012 * (U["EARTH"] - U["CONCRETE"]) * self.nn * self.N * self.LF * self.Gb
        self.Re = 0.012 * U["CONCRETE"] * self.nn * (Re_1 + Re_2) + Re_corr
        self.Rca = self.Ri + self.Rsd + self.Rd + self.Re

    def calculate_cable_core_temperature(self):
        self.Tc = self.Ta + (self.I_given / 1000) ** 2 * self.R_dc * self.Rca
        self.Tc = round(self.Tc, 1)


if __name__ == "__main__":
    specs = dict(
    	# cable details
    	segmental_conductor_diameter=1.543,
        insulation_outer_diameter=2.113,
        sheath_outer_diameter=2.373,
        # key parameters input
        ambient_temperature=25,
        current=700,
        conductor_resistance=8.6,
        conductors_per_cable=1,
        number_of_cables=6,
        duct_inner_diameter=5,
        duct_thickness=0.25,
        duct_bank_rows=3,
        duct_bank_columns=2,
        depth_of_reference_cable=43.5,
        vertical_distance_between_duct_banks=9,
        horizontal_distance_between_duct_banks=9,
        short_side_of_duct=18,
        long_side_of_duct=27,
        depth_to_center_of_duct_bank=43.5,
        loss_factor=0.8,
    )
    duct = DuctBank(**specs)
    print(f"Cable core temperature: {duct.Tc} (in deg C)")
