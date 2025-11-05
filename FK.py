import numpy as np
import sympy as sp

def dh(theta, d, a, alpha):
    ct = np.cos(theta)
    st = np.sin(theta)
    ca = np.cos(alpha)
    sa = np.sin(alpha)
    return sp.Matrix([
        [ct, -st * ca,  st * sa, a * ct],
        [st,  ct * ca, -ct * sa, a * st],
        [ 0,       sa,      ca,      d],
        [ 0,        0,       0,      1]
    ])

def Cinematique_Directe(d1, theta2, theta3, theta4):
    theta2 = np.deg2rad(theta2)
    theta3 = np.deg2rad(theta3)
    theta4 = np.deg2rad(theta4)

    H10 = dh(0, 215 + d1, 0, 0)
    H21 = dh(theta2, 0, 800, np.pi/2)
    H32 = dh(theta3 + np.pi/2, 0, 0, np.pi/2)
    H43 = dh(theta4, 400, 0, 0)

    H0R = sp.Matrix([
        [0, -1, 0, 0],
        [1,  0, 0, -500],
        [0,  0, 1, 0],
        [0,  0, 0, 1]
    ])

    Ho4 = sp.Matrix([
        [1, 0, 0,   0],
        [0, 1, 0,   0],
        [0, 0, 1, 100],
        [0, 0, 0,   1]
    ])

    # Full symbolic transformation chain
    H = H0R * H10 * H21 * H32 * H43 * Ho4

    if abs(H[2, 0]) == 1:
        ey = -H[2, 0] * 90
        ex = 0
        ez = np.arctan2(-H[2, 0]*H[1, 2], H[1, 1]) * 180/np.pi
    else:
        ey = np.arctan2(-float(H[2, 0]), np.sqrt(float(H[0, 0])**2 + float(H[1, 0])**2)) * 180/np.pi
        cp = np.cos(np.deg2rad(ey))
        ez = np.arctan2(float(H[1, 0])/cp, float(H[0, 0])/cp) * 180/np.pi
        ex = np.arctan2(float(H[2, 1])/cp, float(H[2, 2])/cp) * 180/np.pi

    pose = {
        'x': float(H[0, 3]),
        'y': float(H[1, 3]),
        'z': float(H[2, 3]),
        'ex': float(ex),
        'ey': float(ey),
        'ez': float(ez)
    }
    return pose

if __name__ == "__main__":
    import sys
    if len(sys.argv) != 5:
        print(f"Usage: python {sys.argv[0]} d1 theta2 theta3 theta4")
        sys.exit(1)
    # Parse inputs as floats
    d1 = float(sys.argv[1])
    theta2 = float(sys.argv[2])
    theta3 = float(sys.argv[3])
    theta4 = float(sys.argv[4])

    pose = Cinematique_Directe(d1, theta2, theta3, theta4)
    print("Pose Output:")
    for k, v in pose.items():
        print(f"{k}: {v}")
