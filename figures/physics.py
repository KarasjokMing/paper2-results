"""
Physics computation: 2-atom 3-waveguide scattering amplitudes.
Angular inputs are pi-normalized; conversion to radians is internal.
GPU (torch) and CPU (numpy) backends.

Corrected formula: Gamma scaling factors from unified_scan/physics_core.py.
The bare coupling constants (2, 4, 14, 6) are scaled by the appropriate
power of Gamma to ensure probability conservation (unitarity).
"""
import math
import cmath
import numpy as np
from typing import Tuple, Optional

_TORCH = False
try:
    import torch
    _TORCH = True
except ImportError:
    pass

_PI = math.pi


def _to_rad(*angles_pi: float) -> tuple:
    return tuple(a * _PI for a in angles_pi)


# ═══════════════════════════════════════════════════════════════
#  Numpy CPU (scalar/vector)
# ═══════════════════════════════════════════════════════════════

def compute_numpy(Gamma: float, Delta: float,
                  th_a1_pi, th_a2_pi, th_b2_pi, th_b3_pi,
                  phi1_pi, phi2_pi, phi3_pi) -> np.ndarray:
    """Return [T1,R1,T2,R2,T3,R3] as float64 array."""
    th = _to_rad(th_a1_pi, th_a2_pi, th_b2_pi, th_b3_pi)
    ph = _to_rad(phi1_pi, phi2_pi, phi3_pi)
    th_a1, th_a2, th_b2, th_b3 = th
    phi1, phi2, phi3 = ph

    tau_a1 = 2 * Gamma * cmath.exp(1j * phi1) * math.cos(th_a1)
    tau_a2 = 2 * Gamma * cmath.exp(1j * phi2) * math.cos(th_a2)
    tau_b2 = 2 * Gamma * cmath.exp(1j * phi2) * math.cos(th_b2)
    tau_b3 = 2 * Gamma * cmath.exp(1j * phi3) * math.cos(th_b3)
    d = th_a2 - th_b2

    M = (-2 * Gamma ** 2 * math.cos(d)
         - 4 * Gamma ** 2 * cmath.exp(2j * phi2) * (math.sin(d / 2)) ** 2
         + tau_a1 * tau_b2 + (tau_a1 + tau_a2) * tau_b3
         + (2 * Gamma - 1j * Delta) * (tau_a2 + tau_b2)
         + (4 * Gamma - 1j * Delta) * (tau_a1 + tau_b3)
         + 14 * Gamma ** 2 - 8j * Gamma * Delta - Delta ** 2)

    e_ia1 = cmath.exp(1j * th_a1)
    e_ia2 = cmath.exp(1j * th_a2)
    e_ib2 = cmath.exp(1j * th_b2)
    e_ib3 = cmath.exp(1j * th_b3)
    e_ip2 = cmath.exp(1j * phi2)

    B = 1 + e_ia2 / e_ib2 + e_ip2 * e_ia2 + e_ip2 / e_ib2
    P1p = 1 + cmath.exp(1j * (phi1 + th_a1))
    P1m = 1 + cmath.exp(1j * (phi1 - th_a1))
    P3p = 1 + cmath.exp(1j * (phi3 + th_b3))
    P3m = 1 + cmath.exp(1j * (phi3 - th_b3))

    t1 = -(2 * Gamma ** 2 * math.cos(d)
           - 2j * Gamma * (tau_b2 + tau_b3 + 4 * Gamma - 1j * Delta) / e_ia1 * math.sin(phi1)
           + 4 * Gamma ** 2 * cmath.exp(2j * phi2) * (math.sin(d / 2)) ** 2
           - tau_a2 * tau_b3 + 1j * Delta * tau_b2
           - (2 * Gamma - 1j * Delta) * (tau_a2 + tau_b3)
           - (6 * Gamma ** 2 - 6j * Gamma * Delta - Delta ** 2)) / M
    r1 = -(P1m * P1p * Gamma * (tau_b2 + tau_b3 + 4 * Gamma - 1j * Delta)) / M
    t2 = -(P1p * Gamma * ((1 + e_ip2 * e_ia2) * (tau_b3 + 2 * Gamma - 1j * Delta)
                          + Gamma * (1 - cmath.exp(2j * phi2)) * (1 - e_ia2 / e_ib2))) / M
    r2 = -(P1p * Gamma * ((1 + e_ip2 / e_ia2) * (tau_b3 + 2 * Gamma - 1j * Delta)
                          + Gamma * (1 - cmath.exp(2j * phi2)) * (1 - e_ib2 / e_ia2))) / M
    t3 = (P1p * B * Gamma ** 2 * P3p) / M
    r3 = (P1p * B * Gamma ** 2 * P3m) / M

    return np.array([abs(t1)**2, abs(r1)**2, abs(t2)**2,
                     abs(r2)**2, abs(t3)**2, abs(r3)**2], dtype=np.float64)


# ═══════════════════════════════════════════════════════════════
#  Torch GPU (batched)
# ═══════════════════════════════════════════════════════════════

if _TORCH:

    @torch.inference_mode()
    def compute_torch(Gamma: torch.Tensor, Delta: torch.Tensor,
                      th_a1_pi, th_a2_pi, th_b2_pi, th_b3_pi,
                      phi1_pi, phi2_pi, phi3_pi,
                      device='cuda') -> Tuple[torch.Tensor, ...]:
        """Batch GPU: return (T1,R1,T2,R2,T3,R4) each [N] float32."""
        pi = torch.tensor(_PI, device=device, dtype=torch.float32)

        def _ei(c, s):
            return torch.complex(c, s)

        th_a1, th_a2 = th_a1_pi * pi, th_a2_pi * pi
        th_b2, th_b3 = th_b2_pi * pi, th_b3_pi * pi
        phi1, phi2, phi3 = phi1_pi * pi, phi2_pi * pi, phi3_pi * pi

        ca1, sa1 = torch.cos(th_a1), torch.sin(th_a1)
        ca2, sa2 = torch.cos(th_a2), torch.sin(th_a2)
        cb2, sb2 = torch.cos(th_b2), torch.sin(th_b2)
        cb3, sb3 = torch.cos(th_b3), torch.sin(th_b3)

        e_p1 = _ei(torch.cos(phi1), torch.sin(phi1))
        e_p2 = _ei(torch.cos(phi2), torch.sin(phi2))
        e_p3 = _ei(torch.cos(phi3), torch.sin(phi3))

        tau_a1 = 2 * Gamma * e_p1 * ca1
        tau_a2 = 2 * Gamma * e_p2 * ca2
        tau_b2 = 2 * Gamma * e_p2 * cb2
        tau_b3 = 2 * Gamma * e_p3 * cb3

        d = th_a2 - th_b2
        cd = torch.cos(d)
        sd2 = torch.sin(d / 2)
        e2p2 = _ei(torch.cos(2 * phi2), torch.sin(2 * phi2))

        M = (-2.0 * Gamma ** 2 * cd - 4.0 * Gamma ** 2 * e2p2 * sd2 ** 2
             + tau_a1 * tau_b2 + (tau_a1 + tau_a2) * tau_b3
             + (2.0 * Gamma - 1j * Delta) * (tau_a2 + tau_b2)
             + (4.0 * Gamma - 1j * Delta) * (tau_a1 + tau_b3)
             + 14.0 * Gamma ** 2 - 8j * Gamma * Delta - Delta ** 2)

        e_ia1 = _ei(ca1, sa1)
        e_ia2 = _ei(ca2, sa2)
        e_ib2 = _ei(cb2, sb2)
        e_ib3 = _ei(cb3, sb3)

        B = 1 + e_ia2 / e_ib2 + e_p2 * e_ia2 + e_p2 / e_ib2
        P1p = 1 + e_p1 * e_ia1
        P1m = 1 + e_p1 / e_ia1
        P3p = 1 + e_p3 * e_ib3
        P3m = 1 + e_p3 / e_ib3

        t1 = -(2.0 * Gamma ** 2 * cd
               - 2j * Gamma * (tau_b2 + tau_b3 + 4.0 * Gamma - 1j * Delta) / e_ia1 * torch.sin(phi1)
               + 4.0 * Gamma ** 2 * e2p2 * sd2 ** 2
               - tau_a2 * tau_b3 + 1j * Delta * tau_b2
               - (2.0 * Gamma - 1j * Delta) * (tau_a2 + tau_b3)
               - (6.0 * Gamma ** 2 - 6j * Gamma * Delta - Delta ** 2)) / M
        r1 = -(P1m * P1p * Gamma * (tau_b2 + tau_b3 + 4.0 * Gamma - 1j * Delta)) / M
        t2 = -(P1p * Gamma * ((1 + e_p2 * e_ia2) * (tau_b3 + 2.0 * Gamma - 1j * Delta)
                              + Gamma * (1 - e2p2) * (1 - e_ia2 / e_ib2))) / M
        r2 = -(P1p * Gamma * ((1 + e_p2 / e_ia2) * (tau_b3 + 2.0 * Gamma - 1j * Delta)
                              + Gamma * (1 - e2p2) * (1 - e_ib2 / e_ia2))) / M
        t3 = (P1p * B * Gamma ** 2 * P3p) / M
        r3 = (P1p * B * Gamma ** 2 * P3m) / M

        return (torch.abs(t1)**2, torch.abs(r1)**2, torch.abs(t2)**2,
                torch.abs(r2)**2, torch.abs(t3)**2, torch.abs(r3)**2)
