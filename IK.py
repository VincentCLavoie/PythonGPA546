import numpy as np
import sys

def Cinematique_Inverse(x, y, z, ex, ey, ez):
    ex = np.deg2rad(ex)
    ey = np.deg2rad(ey)
    ez = np.deg2rad(ez)
    H_outil_R = np.array([
        [np.cos(ez)*np.cos(ey), np.cos(ez)*np.sin(ey)*np.sin(ex)-np.sin(ez)*np.cos(ex), np.cos(ez)*np.sin(ey)*np.cos(ex)+np.sin(ez)*np.sin(ex), x],
        [np.sin(ez)*np.cos(ey), np.sin(ez)*np.sin(ey)*np.sin(ex)+np.cos(ez)*np.cos(ex), np.sin(ez)*np.sin(ey)*np.cos(ex)-np.cos(ez)*np.sin(ex), y],
        [-np.sin(ey), np.cos(ey)*np.sin(ex), np.cos(ey)*np.cos(ex), z],
        [0, 0, 0, 1]
    ])
    H0R = np.array([
        [0, -1, 0, 0],
        [1,  0, 0, -500],
        [0,  0, 1, 0],
        [0,  0, 0, 1]
    ])
    Ho4 = np.array([
        [1, 0, 0,   0],
        [0, 1, 0,   0],
        [0, 0, 1, 100],
        [0, 0, 0,   1]
    ])
    P = np.linalg.inv(H0R) @ H_outil_R @ np.linalg.inv(Ho4)
    nx, ny, nz = P[0, 0], P[1, 0], P[2, 0]
    ox, oy, oz = P[0, 1], P[1, 1], P[2, 1]
    ax, ay, az = P[0, 2], P[1, 2], P[2, 2]
    px, py, pz = P[0, 3], P[1, 3], P[2, 3]
    d1 = -(400 * az - pz + 215)
    s2 = -(400 * ay - py) / 800
    c2 = -(400 * ax - px) / 800
    theta2 = np.arctan2(s2, c2)
    s3 = ay / np.sin(theta2)
    c3 = -az
    theta3 = np.arctan2(s3, c3) - np.pi / 2
    denom = nx * oy - ny * ox
    s4 = (np.cos(theta2) * ox + oy * np.sin(theta2)) / denom
    c4 = -(np.cos(theta2) * nx + ny * np.sin(theta2)) / denom
    theta4 = np.arctan2(s4, c4)
    theta2 = np.rad2deg(theta2)
    theta3 = np.rad2deg(theta3)
    theta4 = np.rad2deg(theta4)
    return {
        'd1': d1,
        'theta2': theta2,
        'theta3': theta3,
        'theta4': theta4
    }

if __name__ == "__main__":
    if len(sys.argv) != 7:
        print(f"Usage: python {sys.argv[0]} x y z ex ey ez")
        sys.exit(1)
    x = float(sys.argv[1])
    y = float(sys.argv[2])
    z = float(sys.argv[3])
    ex = float(sys.argv[4])
    ey = float(sys.argv[5])
    ez = float(sys.argv[6])
    result = Cinematique_Inverse(x, y, z, ex, ey, ez)
    print("Inverse Kinematics Result:")
    for k, v in result.items():
        print(f"{k}: {v}")
