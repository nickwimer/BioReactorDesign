import math
import os
from functools import reduce
from pathlib import Path

import numpy as np

from bird.utilities.parser import parse_yaml

CUBIC_IN_TO_L = 0.0163871  # conversion factor from cubic inches to liters


class StirredTankReactor:
    def __init__(
        self,
        units,
        tank_diameter,
        impeller_tip_diameter,
        reactor_height,
        nimpellers,
        impeller_centers,
        blade_width,
        blade_length,
        inner_blade_length,
        baffle_width,
        hub_height_width,
        polyrad,
        reactor_bottom,
        nr,
        nz,
        n_poly,
        n_azimuth,
        nbaffles,
        n_fins_per_impeller,
        blade_pitch,
        impeller_scale,
        aspect_ratio,
        target_volume_L,
        round_bottom,
        bottom_inlet,
        mrf_vertical_buffer_cells=0,
        mrf_radial_buffer_cells=0,
    ):
        # Loop through params and setattr v to self.k
        for k, v in locals().items():
            if k != "self":
                setattr(self, k, v)

        # Convert degrees to radians for blade pitch
        self.blade_pitch = [angle * np.pi / 180.0 for angle in blade_pitch]

        # Solve for the tank parameters using the provided geometry or aspect ratio
        (
            self.tank_diameter,
            self.reactor_height,
            self.aspect_ratio,
            self.final_volume_L,
        ) = StirredTankReactor.solve_cylinder(
            target_volume_L=target_volume_L,
            aspect_ratio=aspect_ratio,
            tank_diameter=tank_diameter,
            reactor_height=reactor_height,
        )
        print(
            f"Tank diameter: {self.tank_diameter:.2f}, "
            f"Reactor height: {self.reactor_height:.2f}, "
            f"Aspect ratio: {self.aspect_ratio:.2f}, "
            f"Final volume (L): {self.final_volume_L:.2f}"
        )

        self.impeller_tip_diameter = 0.5 * self.tank_diameter
        self.baffle_width = 0.075 * self.tank_diameter

        # compute least common multiple of n_fins_per_impeller using greatest common divisor
        def _lcm(a, b):
            return a * b // math.gcd(a, b)

        self.hub_diameter = (
            self.impeller_tip_diameter - 2 * self.blade_length
        )  # (old Dh)
        self.mrf_region_diameter = (
            self.impeller_tip_diameter + self.tank_diameter - 2 * self.baffle_width
        ) / 2  # (old Dmrf)

        # Radial MRF buffer measured outward from the local blade tip radius
        stator_inner_radius = self.tank_diameter / 2.0 - self.baffle_width
        tank_radius = self.tank_diameter / 2.0
        radial_buffer_cells = max(0, int(self.mrf_radial_buffer_cells))
        radial_buffer_thickness = (
            radial_buffer_cells / float(nr) if radial_buffer_cells > 0 else 0.0
        )
        has_radial_buffer_interface = radial_buffer_cells > 0

        def _section_radii(scale_val):
            tip_radius = scale_val * self.impeller_tip_diameter / 2
            radii = [
                scale_val * (self.hub_diameter / 2 - self.inner_blade_length),
                scale_val * self.hub_diameter / 2,
                tip_radius,
            ]
            if has_radial_buffer_interface:
                # Keep a dedicated MRF outer interface beyond the tip
                mrf_buffer_outer_radius = min(
                    tip_radius + radial_buffer_thickness,
                    stator_inner_radius - 1e-12,
                )
                if mrf_buffer_outer_radius <= tip_radius:
                    mrf_buffer_outer_radius = tip_radius + 1e-12
                radii.append(mrf_buffer_outer_radius)
            radii.extend([stator_inner_radius, tank_radius])
            return np.array(radii)

        base_counts = [self.nbaffles] + self.n_fins_per_impeller
        n_base = reduce(_lcm, base_counts)
        print(f"Least common multiple of baffles and fins: {n_base}")

        self.nsplits = 2 * n_base  # we need twice the number of splits
        self.dangle = 2.0 * np.pi / float(self.nsplits)

        if self.round_bottom:
            self.curved_bottom_center = [
                0.0,
                0.0,
                self.reactor_bottom + self.reactor_height / 3,
            ]
            curved_bottom_edge = [
                self.tank_diameter / 2,
                0.0,
                self.reactor_bottom,
            ]
            self.curved_bottom_radius = np.sqrt(
                (curved_bottom_edge[0] - self.curved_bottom_center[0]) ** 2
                + (curved_bottom_edge[1] - self.curved_bottom_center[1]) ** 2
                + (curved_bottom_edge[2] - self.curved_bottom_center[2]) ** 2
            )

        self.circradii = _section_radii(self.impeller_scale[0])
        self.ncirc = len(self.circradii)
        self.hub_circ = 1
        self.inhub_circ = self.hub_circ - 1  # circle inside hub
        self.rot_circ = self.hub_circ + 1
        self.tank_circ = self.ncirc - 1
        self.mrf_circ = self.rot_circ + (1 if has_radial_buffer_interface else 0)

        self.reacthts = [reactor_bottom]
        self.baff_sections = []
        self.baff_volumes = []
        self.hub_volumes = []
        self.mrf_volumes = []
        count = 1
        self.angle_offsets = [0.0]

        # Buffer cells for MRF region above and below impeller (vertical direction)
        nbuf = int(self.mrf_vertical_buffer_cells)
        buffer_ht = nbuf / float(nz) if nbuf > 0 else 0.0

        for n_imp in range(self.nimpellers):
            pitch = self.blade_pitch[n_imp]
            tmp_len = self.blade_width
            tip_rad = self.impeller_scale[n_imp] * (self.impeller_tip_diameter / 2)
            dz = tmp_len * np.cos(pitch)
            # prevent blade from going below rotor hub
            dz_min = self.hub_height_width * 1.05
            if dz < dz_min:
                dz = dz_min
            dtheta = (tmp_len * np.sin(pitch)) / max(tip_rad, 1e-12)
            zc = self.reactor_bottom + self.impeller_centers[n_imp]

            def _theta_offset(z):
                return (z - zc) * dtheta / dz

            z0 = zc - dz / 2.0

            # Add buffer section below; make sure to not go below bot or touch next imp
            if buffer_ht > 0.0 and (z0 - buffer_ht) > self.reacthts[-1]:
                zb = z0 - buffer_ht
                self.reacthts.append(zb)
                self.circradii = np.append(
                    self.circradii,
                    _section_radii(self.impeller_scale[n_imp]),
                )
                self.mrf_volumes.append(count)
                self.angle_offsets.append(_theta_offset(z0))  # zero twist in buffer
                count += 1

            self.reacthts.append(z0)
            self.circradii = np.append(
                self.circradii,
                _section_radii(self.impeller_scale[n_imp]),
            )

            self.baff_sections.append(count)
            self.baff_volumes.append(count)
            self.angle_offsets.append(_theta_offset(z0))
            count = count + 1

            z1 = zc - self.hub_height_width / 2.0
            self.reacthts.append(z1)
            self.circradii = np.append(
                self.circradii,
                _section_radii(self.impeller_scale[n_imp]),
            )

            self.baff_sections.append(count)
            self.baff_volumes.append(count)
            self.hub_volumes.append(count)
            self.angle_offsets.append(_theta_offset(z1))
            count = count + 1

            z2 = zc + self.hub_height_width / 2.0
            self.reacthts.append(z2)
            self.circradii = np.append(
                self.circradii,
                _section_radii(self.impeller_scale[n_imp]),
            )

            self.baff_sections.append(count)
            self.baff_volumes.append(count)
            self.angle_offsets.append(_theta_offset(z2))
            count = count + 1

            z3 = zc + dz / 2.0
            self.reacthts.append(z3)
            self.circradii = np.append(
                self.circradii,
                _section_radii(self.impeller_scale[n_imp]),
            )
            self.baff_sections.append(count)
            self.angle_offsets.append(_theta_offset(z3))
            count = count + 1

            if buffer_ht > 0.0:
                zt = z3 + buffer_ht
                # Only add if it doesn't exceed next impeller or tank top
                if n_imp < self.nimpellers - 1:
                    next_zc = self.reactor_bottom + self.impeller_centers[n_imp + 1]
                    next_dz_val = self.blade_width * np.cos(self.blade_pitch[n_imp + 1])
                    if next_dz_val < self.hub_height_width * 1.05:
                        next_dz_val = self.hub_height_width * 1.05
                    next_z = next_zc - next_dz_val / 2.0
                else:
                    next_z = self.reactor_bottom + self.reactor_height

                if zt < next_z:
                    self.reacthts.append(zt)
                    self.circradii = np.append(
                        self.circradii,
                        _section_radii(self.impeller_scale[n_imp]),
                    )

                    self.angle_offsets.append(_theta_offset(z3))  # zero twist in buffer

                    self.mrf_volumes.append(count - 1)
                    count = count + 1

        self.reacthts.append(self.reactor_bottom + self.reactor_height)
        self.circradii = np.append(
            self.circradii,
            _section_radii(self.impeller_scale[-1]),
        )
        self.angle_offsets.append(0.0)
        self.nsections = len(self.reacthts)
        self.circradii = self.circradii.reshape(self.nsections, self.ncirc)
        self.nvolumes = self.nsections - 1
        self.meshz = nz * np.diff(self.reacthts)
        self.meshz = self.meshz.astype(int) + 1  # avoid zero mesh elements

        # mapping from section to impeller number
        self.section2imp = -1 * np.ones(self.nsections, dtype=int)
        # index 0-3 -> impeller0, 4-7 -> impeller1, etc.
        for j, sec in enumerate(self.baff_sections):
            self.section2imp[sec] = j // 4

        self.all_volumes = range(self.nvolumes)
        self.nonbaff_volumes = [
            sec for sec in self.all_volumes if sec not in self.baff_volumes
        ]
        self.nonstem_volumes = (
            list(range(self.hub_volumes[0])) if len(self.hub_volumes) > 0 else [0, 1]
        )

        # note: stem_volumes include hub volumes also
        # these are volumes where we miss out polygon block
        self.stem_volumes = [
            sec for sec in self.all_volumes if sec not in self.nonstem_volumes
        ]

        # removes hub_volumes here for declaring patches
        self.only_stem_volumes = [
            sec for sec in self.stem_volumes if sec not in self.hub_volumes
        ]

        # Define MRF volumes by merging baffled volumes with buffer volumes
        mrf_vols = set(self.baff_volumes) | set(self.mrf_volumes)
        self.mrf_volumes = sorted(list(mrf_vols))

        # # increase grid points in the impeller section
        # for i in self.baff_volumes:
        #     self.meshz[i] *= 2

        self.meshr = nr * np.diff(self.circradii)

        # adding polygon to hub mesh resolution
        self.meshr = np.append(nr * polyrad, self.meshr)
        self.meshr = self.meshr.astype(int)
        self.meshr += 1  # to avoid being zero

        self.centeroffset = 1  # one point on the axis
        self.polyoffset = self.nsplits  # number of points on polygon
        self.npts_per_section = (
            self.centeroffset + self.polyoffset + self.ncirc * self.nsplits
        )  # center+polygon+circles

    @classmethod
    def from_file(cls, yamlfile):
        if ".yaml" not in yamlfile:
            yamlfile += ".yaml"
        in_dict = parse_yaml(yamlfile)
        react_dict = {**in_dict["geometry"], **in_dict["mesh"]}
        return cls(**react_dict)

    def solve_cylinder(
        target_volume_L: float,
        aspect_ratio: float | None = None,  # AR = H/D
        tank_diameter: float | None = None,  # inches
        reactor_height: float | None = None,  # inches
    ):
        """
        Solves for the cylinder tank geometry using the aspect ratio, tank diameter, or
        height.
        """

        def _cylinder_volume_L(diameter_in: float, height_in: float) -> float:
            return (math.pi * (diameter_in / 2.0) ** 2 * height_in) * CUBIC_IN_TO_L

        if target_volume_L is None:
            assert (
                tank_diameter is not None and reactor_height is not None
            ), "Provide target_volume_L or both tank_diameter and reactor_height."
            aspect_ratio_final = reactor_height / tank_diameter

            water_level_height_final = 0.9 * reactor_height
            water_level_volume_L = _cylinder_volume_L(
                tank_diameter, water_level_height_final
            )
            print(
                f"Water level height for 90% volume: {water_level_height_final:.2f} inches, "
                f"which corresponds to {water_level_volume_L:.2f} L"
            )
            return (
                tank_diameter,
                reactor_height,
                aspect_ratio_final,
                _cylinder_volume_L(tank_diameter, reactor_height),
            )

        volume_in3 = target_volume_L / CUBIC_IN_TO_L

        if tank_diameter is not None and reactor_height is not None:
            raise ValueError(
                "Provide only one of tank_diameter or reactor_height (or neither)."
            )

        if tank_diameter is not None:
            tank_diameter_final = float(tank_diameter)
            reactor_height_final = 4.0 * volume_in3 / (math.pi * tank_diameter_final**2)
        elif reactor_height is not None:
            reactor_height_final = float(reactor_height)
            tank_diameter_final = math.sqrt(
                4.0 * volume_in3 / (math.pi * reactor_height_final)
            )
        else:
            if aspect_ratio is None:
                raise ValueError(
                    "If neither diameter nor height is provided, "
                    "aspect_ratio is required."
                )
            aspect_ratio_final = float(aspect_ratio)
            tank_diameter_final = (
                4.0 * volume_in3 / (math.pi * aspect_ratio_final)
            ) ** (1.0 / 3.0)
            reactor_height_final = aspect_ratio_final * tank_diameter_final

        aspect_ratio_final = reactor_height_final / tank_diameter_final
        final_volume_L = (
            math.pi * (tank_diameter_final / 2.0) ** 2 * reactor_height_final
        ) * CUBIC_IN_TO_L

        # Compute the height that yeils 90% of the target volume for water level line
        water_level_height_final = 0.9 * reactor_height_final
        water_level_volume_L = _cylinder_volume_L(
            tank_diameter_final, water_level_height_final
        )
        print(
            f"Water level height for 90% volume: {water_level_height_final:.2f} inches, "
            f"which corresponds to {water_level_volume_L:.2f} L"
        )

        return (
            tank_diameter_final,
            reactor_height_final,
            aspect_ratio_final,
            final_volume_L,
        )
