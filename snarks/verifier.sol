// This file is LGPL3 Licensed

/**
 * @title Elliptic curve operations on twist points for alt_bn128
 * @author Mustafa Al-Bassam (mus@musalbas.com)
 */
library BN256G2 {
    uint256 internal constant FIELD_MODULUS = 0x30644e72e131a029b85045b68181585d97816a916871ca8d3c208c16d87cfd47;
    uint256 internal constant TWISTBX = 0x2b149d40ceb8aaae81be18991be06ac3b5b4c5e559dbefa33267e6dc24a138e5;
    uint256 internal constant TWISTBY = 0x9713b03af0fed4cd2cafadeed8fdf4a74fa084e52d1852e4a2bd0685c315d2;
    uint internal constant PTXX = 0;
    uint internal constant PTXY = 1;
    uint internal constant PTYX = 2;
    uint internal constant PTYY = 3;
    uint internal constant PTZX = 4;
    uint internal constant PTZY = 5;

    /**
     * @notice Add two twist points
     * @param pt1xx Coefficient 1 of x on point 1
     * @param pt1xy Coefficient 2 of x on point 1
     * @param pt1yx Coefficient 1 of y on point 1
     * @param pt1yy Coefficient 2 of y on point 1
     * @param pt2xx Coefficient 1 of x on point 2
     * @param pt2xy Coefficient 2 of x on point 2
     * @param pt2yx Coefficient 1 of y on point 2
     * @param pt2yy Coefficient 2 of y on point 2
     * @return (pt3xx, pt3xy, pt3yx, pt3yy)
     */
    function ECTwistAdd(
        uint256 pt1xx, uint256 pt1xy,
        uint256 pt1yx, uint256 pt1yy,
        uint256 pt2xx, uint256 pt2xy,
        uint256 pt2yx, uint256 pt2yy
    ) public pure returns (
        uint256, uint256,
        uint256, uint256
    ) {
        if (
            pt1xx == 0 && pt1xy == 0 &&
            pt1yx == 0 && pt1yy == 0
        ) {
            if (!(
                pt2xx == 0 && pt2xy == 0 &&
                pt2yx == 0 && pt2yy == 0
            )) {
                assert(_isOnCurve(
                    pt2xx, pt2xy,
                    pt2yx, pt2yy
                ));
            }
            return (
                pt2xx, pt2xy,
                pt2yx, pt2yy
            );
        } else if (
            pt2xx == 0 && pt2xy == 0 &&
            pt2yx == 0 && pt2yy == 0
        ) {
            assert(_isOnCurve(
                pt1xx, pt1xy,
                pt1yx, pt1yy
            ));
            return (
                pt1xx, pt1xy,
                pt1yx, pt1yy
            );
        }

        assert(_isOnCurve(
            pt1xx, pt1xy,
            pt1yx, pt1yy
        ));
        assert(_isOnCurve(
            pt2xx, pt2xy,
            pt2yx, pt2yy
        ));

        uint256[6] memory pt3 = _ECTwistAddJacobian(
            pt1xx, pt1xy,
            pt1yx, pt1yy,
            1,     0,
            pt2xx, pt2xy,
            pt2yx, pt2yy,
            1,     0
        );

        return _fromJacobian(
            pt3[PTXX], pt3[PTXY],
            pt3[PTYX], pt3[PTYY],
            pt3[PTZX], pt3[PTZY]
        );
    }

    /**
     * @notice Multiply a twist point by a scalar
     * @param s     Scalar to multiply by
     * @param pt1xx Coefficient 1 of x
     * @param pt1xy Coefficient 2 of x
     * @param pt1yx Coefficient 1 of y
     * @param pt1yy Coefficient 2 of y
     * @return (pt2xx, pt2xy, pt2yx, pt2yy)
     */
    function ECTwistMul(
        uint256 s,
        uint256 pt1xx, uint256 pt1xy,
        uint256 pt1yx, uint256 pt1yy
    ) public pure returns (
        uint256, uint256,
        uint256, uint256
    ) {
        uint256 pt1zx = 1;
        if (
            pt1xx == 0 && pt1xy == 0 &&
            pt1yx == 0 && pt1yy == 0
        ) {
            pt1xx = 1;
            pt1yx = 1;
            pt1zx = 0;
        } else {
            assert(_isOnCurve(
                pt1xx, pt1xy,
                pt1yx, pt1yy
            ));
        }

        uint256[6] memory pt2 = _ECTwistMulJacobian(
            s,
            pt1xx, pt1xy,
            pt1yx, pt1yy,
            pt1zx, 0
        );

        return _fromJacobian(
            pt2[PTXX], pt2[PTXY],
            pt2[PTYX], pt2[PTYY],
            pt2[PTZX], pt2[PTZY]
        );
    }

    /**
     * @notice Get the field modulus
     * @return The field modulus
     */
    function GetFieldModulus() public pure returns (uint256) {
        return FIELD_MODULUS;
    }

    function submod(uint256 a, uint256 b, uint256 n) internal pure returns (uint256) {
        return addmod(a, n - b, n);
    }

    function _FQ2Mul(
        uint256 xx, uint256 xy,
        uint256 yx, uint256 yy
    ) internal pure returns(uint256, uint256) {
        return (
            submod(mulmod(xx, yx, FIELD_MODULUS), mulmod(xy, yy, FIELD_MODULUS), FIELD_MODULUS),
            addmod(mulmod(xx, yy, FIELD_MODULUS), mulmod(xy, yx, FIELD_MODULUS), FIELD_MODULUS)
        );
    }

    function _FQ2Muc(
        uint256 xx, uint256 xy,
        uint256 c
    ) internal pure returns(uint256, uint256) {
        return (
            mulmod(xx, c, FIELD_MODULUS),
            mulmod(xy, c, FIELD_MODULUS)
        );
    }

    function _FQ2Add(
        uint256 xx, uint256 xy,
        uint256 yx, uint256 yy
    ) internal pure returns(uint256, uint256) {
        return (
            addmod(xx, yx, FIELD_MODULUS),
            addmod(xy, yy, FIELD_MODULUS)
        );
    }

    function _FQ2Sub(
        uint256 xx, uint256 xy,
        uint256 yx, uint256 yy
    ) internal pure returns(uint256 rx, uint256 ry) {
        return (
            submod(xx, yx, FIELD_MODULUS),
            submod(xy, yy, FIELD_MODULUS)
        );
    }

    function _FQ2Div(
        uint256 xx, uint256 xy,
        uint256 yx, uint256 yy
    ) internal pure returns(uint256, uint256) {
        (yx, yy) = _FQ2Inv(yx, yy);
        return _FQ2Mul(xx, xy, yx, yy);
    }

    function _FQ2Inv(uint256 x, uint256 y) internal pure returns(uint256, uint256) {
        uint256 inv = _modInv(addmod(mulmod(y, y, FIELD_MODULUS), mulmod(x, x, FIELD_MODULUS), FIELD_MODULUS), FIELD_MODULUS);
        return (
            mulmod(x, inv, FIELD_MODULUS),
            FIELD_MODULUS - mulmod(y, inv, FIELD_MODULUS)
        );
    }

    function _isOnCurve(
        uint256 xx, uint256 xy,
        uint256 yx, uint256 yy
    ) internal pure returns (bool) {
        uint256 yyx;
        uint256 yyy;
        uint256 xxxx;
        uint256 xxxy;
        (yyx, yyy) = _FQ2Mul(yx, yy, yx, yy);
        (xxxx, xxxy) = _FQ2Mul(xx, xy, xx, xy);
        (xxxx, xxxy) = _FQ2Mul(xxxx, xxxy, xx, xy);
        (yyx, yyy) = _FQ2Sub(yyx, yyy, xxxx, xxxy);
        (yyx, yyy) = _FQ2Sub(yyx, yyy, TWISTBX, TWISTBY);
        return yyx == 0 && yyy == 0;
    }

    function _modInv(uint256 a, uint256 n) internal pure returns(uint256 t) {
        t = 0;
        uint256 newT = 1;
        uint256 r = n;
        uint256 newR = a;
        uint256 q;
        while (newR != 0) {
            q = r / newR;
            (t, newT) = (newT, submod(t, mulmod(q, newT, n), n));
            (r, newR) = (newR, r - q * newR);
        }
    }

    function _fromJacobian(
        uint256 pt1xx, uint256 pt1xy,
        uint256 pt1yx, uint256 pt1yy,
        uint256 pt1zx, uint256 pt1zy
    ) internal pure returns (
        uint256 pt2xx, uint256 pt2xy,
        uint256 pt2yx, uint256 pt2yy
    ) {
        uint256 invzx;
        uint256 invzy;
        (invzx, invzy) = _FQ2Inv(pt1zx, pt1zy);
        (pt2xx, pt2xy) = _FQ2Mul(pt1xx, pt1xy, invzx, invzy);
        (pt2yx, pt2yy) = _FQ2Mul(pt1yx, pt1yy, invzx, invzy);
    }

    function _ECTwistAddJacobian(
        uint256 pt1xx, uint256 pt1xy,
        uint256 pt1yx, uint256 pt1yy,
        uint256 pt1zx, uint256 pt1zy,
        uint256 pt2xx, uint256 pt2xy,
        uint256 pt2yx, uint256 pt2yy,
        uint256 pt2zx, uint256 pt2zy) internal pure returns (uint256[6] memory pt3) {
            if (pt1zx == 0 && pt1zy == 0) {
                (
                    pt3[PTXX], pt3[PTXY],
                    pt3[PTYX], pt3[PTYY],
                    pt3[PTZX], pt3[PTZY]
                ) = (
                    pt2xx, pt2xy,
                    pt2yx, pt2yy,
                    pt2zx, pt2zy
                );
                return pt3;
            } else if (pt2zx == 0 && pt2zy == 0) {
                (
                    pt3[PTXX], pt3[PTXY],
                    pt3[PTYX], pt3[PTYY],
                    pt3[PTZX], pt3[PTZY]
                ) = (
                    pt1xx, pt1xy,
                    pt1yx, pt1yy,
                    pt1zx, pt1zy
                );
                return pt3;
            }

            (pt2yx,     pt2yy)     = _FQ2Mul(pt2yx, pt2yy, pt1zx, pt1zy); // U1 = y2 * z1
            (pt3[PTYX], pt3[PTYY]) = _FQ2Mul(pt1yx, pt1yy, pt2zx, pt2zy); // U2 = y1 * z2
            (pt2xx,     pt2xy)     = _FQ2Mul(pt2xx, pt2xy, pt1zx, pt1zy); // V1 = x2 * z1
            (pt3[PTZX], pt3[PTZY]) = _FQ2Mul(pt1xx, pt1xy, pt2zx, pt2zy); // V2 = x1 * z2

            if (pt2xx == pt3[PTZX] && pt2xy == pt3[PTZY]) {
                if (pt2yx == pt3[PTYX] && pt2yy == pt3[PTYY]) {
                    (
                        pt3[PTXX], pt3[PTXY],
                        pt3[PTYX], pt3[PTYY],
                        pt3[PTZX], pt3[PTZY]
                    ) = _ECTwistDoubleJacobian(pt1xx, pt1xy, pt1yx, pt1yy, pt1zx, pt1zy);
                    return pt3;
                }
                (
                    pt3[PTXX], pt3[PTXY],
                    pt3[PTYX], pt3[PTYY],
                    pt3[PTZX], pt3[PTZY]
                ) = (
                    1, 0,
                    1, 0,
                    0, 0
                );
                return pt3;
            }

            (pt2zx,     pt2zy)     = _FQ2Mul(pt1zx, pt1zy, pt2zx,     pt2zy);     // W = z1 * z2
            (pt1xx,     pt1xy)     = _FQ2Sub(pt2yx, pt2yy, pt3[PTYX], pt3[PTYY]); // U = U1 - U2
            (pt1yx,     pt1yy)     = _FQ2Sub(pt2xx, pt2xy, pt3[PTZX], pt3[PTZY]); // V = V1 - V2
            (pt1zx,     pt1zy)     = _FQ2Mul(pt1yx, pt1yy, pt1yx,     pt1yy);     // V_squared = V * V
            (pt2yx,     pt2yy)     = _FQ2Mul(pt1zx, pt1zy, pt3[PTZX], pt3[PTZY]); // V_squared_times_V2 = V_squared * V2
            (pt1zx,     pt1zy)     = _FQ2Mul(pt1zx, pt1zy, pt1yx,     pt1yy);     // V_cubed = V * V_squared
            (pt3[PTZX], pt3[PTZY]) = _FQ2Mul(pt1zx, pt1zy, pt2zx,     pt2zy);     // newz = V_cubed * W
            (pt2xx,     pt2xy)     = _FQ2Mul(pt1xx, pt1xy, pt1xx,     pt1xy);     // U * U
            (pt2xx,     pt2xy)     = _FQ2Mul(pt2xx, pt2xy, pt2zx,     pt2zy);     // U * U * W
            (pt2xx,     pt2xy)     = _FQ2Sub(pt2xx, pt2xy, pt1zx,     pt1zy);     // U * U * W - V_cubed
            (pt2zx,     pt2zy)     = _FQ2Muc(pt2yx, pt2yy, 2);                    // 2 * V_squared_times_V2
            (pt2xx,     pt2xy)     = _FQ2Sub(pt2xx, pt2xy, pt2zx,     pt2zy);     // A = U * U * W - V_cubed - 2 * V_squared_times_V2
            (pt3[PTXX], pt3[PTXY]) = _FQ2Mul(pt1yx, pt1yy, pt2xx,     pt2xy);     // newx = V * A
            (pt1yx,     pt1yy)     = _FQ2Sub(pt2yx, pt2yy, pt2xx,     pt2xy);     // V_squared_times_V2 - A
            (pt1yx,     pt1yy)     = _FQ2Mul(pt1xx, pt1xy, pt1yx,     pt1yy);     // U * (V_squared_times_V2 - A)
            (pt1xx,     pt1xy)     = _FQ2Mul(pt1zx, pt1zy, pt3[PTYX], pt3[PTYY]); // V_cubed * U2
            (pt3[PTYX], pt3[PTYY]) = _FQ2Sub(pt1yx, pt1yy, pt1xx,     pt1xy);     // newy = U * (V_squared_times_V2 - A) - V_cubed * U2
    }

    function _ECTwistDoubleJacobian(
        uint256 pt1xx, uint256 pt1xy,
        uint256 pt1yx, uint256 pt1yy,
        uint256 pt1zx, uint256 pt1zy
    ) internal pure returns(
        uint256 pt2xx, uint256 pt2xy,
        uint256 pt2yx, uint256 pt2yy,
        uint256 pt2zx, uint256 pt2zy
    ) {
        (pt2xx, pt2xy) = _FQ2Muc(pt1xx, pt1xy, 3);            // 3 * x
        (pt2xx, pt2xy) = _FQ2Mul(pt2xx, pt2xy, pt1xx, pt1xy); // W = 3 * x * x
        (pt1zx, pt1zy) = _FQ2Mul(pt1yx, pt1yy, pt1zx, pt1zy); // S = y * z
        (pt2yx, pt2yy) = _FQ2Mul(pt1xx, pt1xy, pt1yx, pt1yy); // x * y
        (pt2yx, pt2yy) = _FQ2Mul(pt2yx, pt2yy, pt1zx, pt1zy); // B = x * y * S
        (pt1xx, pt1xy) = _FQ2Mul(pt2xx, pt2xy, pt2xx, pt2xy); // W * W
        (pt2zx, pt2zy) = _FQ2Muc(pt2yx, pt2yy, 8);            // 8 * B
        (pt1xx, pt1xy) = _FQ2Sub(pt1xx, pt1xy, pt2zx, pt2zy); // H = W * W - 8 * B
        (pt2zx, pt2zy) = _FQ2Mul(pt1zx, pt1zy, pt1zx, pt1zy); // S_squared = S * S
        (pt2yx, pt2yy) = _FQ2Muc(pt2yx, pt2yy, 4);            // 4 * B
        (pt2yx, pt2yy) = _FQ2Sub(pt2yx, pt2yy, pt1xx, pt1xy); // 4 * B - H
        (pt2yx, pt2yy) = _FQ2Mul(pt2yx, pt2yy, pt2xx, pt2xy); // W * (4 * B - H)
        (pt2xx, pt2xy) = _FQ2Muc(pt1yx, pt1yy, 8);            // 8 * y
        (pt2xx, pt2xy) = _FQ2Mul(pt2xx, pt2xy, pt1yx, pt1yy); // 8 * y * y
        (pt2xx, pt2xy) = _FQ2Mul(pt2xx, pt2xy, pt2zx, pt2zy); // 8 * y * y * S_squared
        (pt2yx, pt2yy) = _FQ2Sub(pt2yx, pt2yy, pt2xx, pt2xy); // newy = W * (4 * B - H) - 8 * y * y * S_squared
        (pt2xx, pt2xy) = _FQ2Muc(pt1xx, pt1xy, 2);            // 2 * H
        (pt2xx, pt2xy) = _FQ2Mul(pt2xx, pt2xy, pt1zx, pt1zy); // newx = 2 * H * S
        (pt2zx, pt2zy) = _FQ2Mul(pt1zx, pt1zy, pt2zx, pt2zy); // S * S_squared
        (pt2zx, pt2zy) = _FQ2Muc(pt2zx, pt2zy, 8);            // newz = 8 * S * S_squared
    }

    function _ECTwistMulJacobian(
        uint256 d,
        uint256 pt1xx, uint256 pt1xy,
        uint256 pt1yx, uint256 pt1yy,
        uint256 pt1zx, uint256 pt1zy
    ) internal pure returns(uint256[6] memory pt2) {
        while (d != 0) {
            if ((d & 1) != 0) {
                pt2 = _ECTwistAddJacobian(
                    pt2[PTXX], pt2[PTXY],
                    pt2[PTYX], pt2[PTYY],
                    pt2[PTZX], pt2[PTZY],
                    pt1xx, pt1xy,
                    pt1yx, pt1yy,
                    pt1zx, pt1zy);
            }
            (
                pt1xx, pt1xy,
                pt1yx, pt1yy,
                pt1zx, pt1zy
            ) = _ECTwistDoubleJacobian(
                pt1xx, pt1xy,
                pt1yx, pt1yy,
                pt1zx, pt1zy
            );

            d = d / 2;
        }
    }
}


// This file is MIT Licensed.
//
// Copyright 2017 Christian Reitwiessner
// Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
// The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

pragma solidity ^0.5.0;
library Pairing {
    struct G1Point {
        uint X;
        uint Y;
    }
    // Encoding of field elements is: X[0] * z + X[1]
    struct G2Point {
        uint[2] X;
        uint[2] Y;
    }
    /// @return the generator of G1
    function P1() pure internal returns (G1Point memory) {
        return G1Point(1, 2);
    }
    /// @return the generator of G2
    function P2() pure internal returns (G2Point memory) {
        return G2Point(
            [11559732032986387107991004021392285783925812861821192530917403151452391805634,
             10857046999023057135944570762232829481370756359578518086990519993285655852781],
            [4082367875863433681332203403145435568316851327593401208105741076214120093531,
             8495653923123431417604973247489272438418190587263600148770280649306958101930]
        );
    }
    /// @return the negation of p, i.e. p.addition(p.negate()) should be zero.
    function negate(G1Point memory p) pure internal returns (G1Point memory) {
        // The prime q in the base field F_q for G1
        uint q = 21888242871839275222246405745257275088696311157297823662689037894645226208583;
        if (p.X == 0 && p.Y == 0)
            return G1Point(0, 0);
        return G1Point(p.X, q - (p.Y % q));
    }
    /// @return the sum of two points of G1
    function addition(G1Point memory p1, G1Point memory p2) internal returns (G1Point memory r) {
        uint[4] memory input;
        input[0] = p1.X;
        input[1] = p1.Y;
        input[2] = p2.X;
        input[3] = p2.Y;
        bool success;
        assembly {
            success := call(sub(gas, 2000), 6, 0, input, 0xc0, r, 0x60)
            // Use "invalid" to make gas estimation work
            switch success case 0 { invalid() }
        }
        require(success);
    }
    /// @return the sum of two points of G2
    function addition(G2Point memory p1, G2Point memory p2) internal pure returns (G2Point memory r) {
        (r.X[1], r.X[0], r.Y[1], r.Y[0]) = BN256G2.ECTwistAdd(p1.X[1],p1.X[0],p1.Y[1],p1.Y[0],p2.X[1],p2.X[0],p2.Y[1],p2.Y[0]);
    }
    /// @return the product of a point on G1 and a scalar, i.e.
    /// p == p.scalar_mul(1) and p.addition(p) == p.scalar_mul(2) for all points p.
    function scalar_mul(G1Point memory p, uint s) internal returns (G1Point memory r) {
        uint[3] memory input;
        input[0] = p.X;
        input[1] = p.Y;
        input[2] = s;
        bool success;
        assembly {
            success := call(sub(gas, 2000), 7, 0, input, 0x80, r, 0x60)
            // Use "invalid" to make gas estimation work
            switch success case 0 { invalid() }
        }
        require (success);
    }
    /// @return the result of computing the pairing check
    /// e(p1[0], p2[0]) *  .... * e(p1[n], p2[n]) == 1
    /// For example pairing([P1(), P1().negate()], [P2(), P2()]) should
    /// return true.
    function pairing(G1Point[] memory p1, G2Point[] memory p2) internal returns (bool) {
        require(p1.length == p2.length);
        uint elements = p1.length;
        uint inputSize = elements * 6;
        uint[] memory input = new uint[](inputSize);
        for (uint i = 0; i < elements; i++)
        {
            input[i * 6 + 0] = p1[i].X;
            input[i * 6 + 1] = p1[i].Y;
            input[i * 6 + 2] = p2[i].X[0];
            input[i * 6 + 3] = p2[i].X[1];
            input[i * 6 + 4] = p2[i].Y[0];
            input[i * 6 + 5] = p2[i].Y[1];
        }
        uint[1] memory out;
        bool success;
        assembly {
            success := call(sub(gas, 2000), 8, 0, add(input, 0x20), mul(inputSize, 0x20), out, 0x20)
            // Use "invalid" to make gas estimation work
            switch success case 0 { invalid() }
        }
        require(success);
        return out[0] != 0;
    }
    /// Convenience method for a pairing check for two pairs.
    function pairingProd2(G1Point memory a1, G2Point memory a2, G1Point memory b1, G2Point memory b2) internal returns (bool) {
        G1Point[] memory p1 = new G1Point[](2);
        G2Point[] memory p2 = new G2Point[](2);
        p1[0] = a1;
        p1[1] = b1;
        p2[0] = a2;
        p2[1] = b2;
        return pairing(p1, p2);
    }
    /// Convenience method for a pairing check for three pairs.
    function pairingProd3(
            G1Point memory a1, G2Point memory a2,
            G1Point memory b1, G2Point memory b2,
            G1Point memory c1, G2Point memory c2
    ) internal returns (bool) {
        G1Point[] memory p1 = new G1Point[](3);
        G2Point[] memory p2 = new G2Point[](3);
        p1[0] = a1;
        p1[1] = b1;
        p1[2] = c1;
        p2[0] = a2;
        p2[1] = b2;
        p2[2] = c2;
        return pairing(p1, p2);
    }
    /// Convenience method for a pairing check for four pairs.
    function pairingProd4(
            G1Point memory a1, G2Point memory a2,
            G1Point memory b1, G2Point memory b2,
            G1Point memory c1, G2Point memory c2,
            G1Point memory d1, G2Point memory d2
    ) internal returns (bool) {
        G1Point[] memory p1 = new G1Point[](4);
        G2Point[] memory p2 = new G2Point[](4);
        p1[0] = a1;
        p1[1] = b1;
        p1[2] = c1;
        p1[3] = d1;
        p2[0] = a2;
        p2[1] = b2;
        p2[2] = c2;
        p2[3] = d2;
        return pairing(p1, p2);
    }
}

contract Verifier {
    using Pairing for *;
    struct VerifyingKey {
        Pairing.G1Point a;
        Pairing.G2Point b;
        Pairing.G2Point gamma;
        Pairing.G2Point delta;
        Pairing.G1Point[] gammaABC;
    }
    struct Proof {
        Pairing.G1Point A;
        Pairing.G2Point B;
        Pairing.G1Point C;
    }
    function verifyingKey() pure internal returns (VerifyingKey memory vk) {
        vk.a = Pairing.G1Point(uint256(0x11b2f3ea59a41c913142c87dbed1f006b8f93661232fed3cac2a694902718ef0), uint256(0x2a27186304474c020275e5f116ce532da4a24ec38dd212e33f69601cabf9bdc3));
        vk.b = Pairing.G2Point([uint256(0x11155fcda2fff2a9a14feb3d1ad84b091b5894ead9cdc53ff7652db287552507), uint256(0x0bfacbcbd97c7333cff340ab7a15384a3303b1e782e0e1260d13eb75e001e94c)], [uint256(0x1e3e0368c3b9148a3bf3b42d0fe49327cffccd78c843a530a3f47fd4fbb7ceff), uint256(0x0d9d3f2628585b7f8b034024fd436cc0711f7ef51f92e592f5d075715ee30028)]);
        vk.gamma = Pairing.G2Point([uint256(0x09ff0f6b954e5e472a02c474be2ce3039d5e2ebe79b283944f95ff9685df4e66), uint256(0x2d2031d2ca59dbeebce5ce3721533339493103dbb4aa05b68ecb854a2196dc31)], [uint256(0x0aeee4b01f377163e0db5df283d639d634939bf524b40c1d16deffd9bbc5d508), uint256(0x24906da081510616199767bfc53305640dd501f9ede181aab2d4803a7e29722c)]);
        vk.delta = Pairing.G2Point([uint256(0x3056fbe13a558e3a9f0ec06c6ce761bc5a420fcde0450efa6ca354f5e8d8662e), uint256(0x05c254ce7de7a76eb0289f1df00cfcbeff92b522ce752b77006583e046d955c5)], [uint256(0x170ad2b8e48f7eccb1fdc73a1907a6d13ddd2d884e3b7ab688dc788aab69cac4), uint256(0x12c02e27cdd0b7578d70cd54472729fa96adf36697f28d041b08489295d99b7d)]);
        vk.gammaABC = new Pairing.G1Point[](522);
        vk.gammaABC[0] = Pairing.G1Point(uint256(0x23083f71a055e0297aef587245938df03096a9c510efcaa10879408a078ba062), uint256(0x2498f7d4ca7d92514b717af8ed6105fde7253cf7f4f9ae9ea02f4ef4928cc111));
        vk.gammaABC[1] = Pairing.G1Point(uint256(0x05281fe7d33f52a59bfe1f226fcfff60e9fadc34dcfe9f235f13b2beb6a909b0), uint256(0x195c7c4730dbc07dc4ac685ca6e210efc4284194aef716a97099a1a86086cd45));
        vk.gammaABC[2] = Pairing.G1Point(uint256(0x133f9dd3b4f56562f0c27e5e97e36bda8d6cbcf4f7a95a9df9f41108acf0c604), uint256(0x2fca04095df7832d3792c6897cb193f1e32976d8f408da98ca3e855fe7f8238d));
        vk.gammaABC[3] = Pairing.G1Point(uint256(0x183b733f279458293eb527dd6ea93af3fa5cc398d562364cd92699dd432629a7), uint256(0x03bba4a35a1b04fe26eac52433226c32e7729d0213307866acdb1604c183798b));
        vk.gammaABC[4] = Pairing.G1Point(uint256(0x21c6ed332f948f917f7ac0193adf37b0a9fc9cfef384c79fdbf7c79e469dfcd7), uint256(0x0376f50f566750dbe1073686d691694a6b76c02a7986cc453d48329e89fb8b86));
        vk.gammaABC[5] = Pairing.G1Point(uint256(0x0b64df0051a33d26ad1ae658eaa6e07b6207568c674215868b811918e48b8e76), uint256(0x054e598d8c8837acfd73eef25f0d6e8954caa3472dfa31c5da4f581b494b2532));
        vk.gammaABC[6] = Pairing.G1Point(uint256(0x2b93d2ab115fc3ea9fe7bcebe8f65d9429d7e6b3e2397d5ebe00bb5823e2f783), uint256(0x08ca827ce54541a78cd3a8f594746f9bf342931b1eb8a67c6de549f7e0fef8ae));
        vk.gammaABC[7] = Pairing.G1Point(uint256(0x1aa1f6db0ba88badeb85d6759a68146ca8b85a38abb8da9e2c373611f8943480), uint256(0x2ea32e643b860d19217a4ec3280e537b691feaa74648094e27100569ce57d4ad));
        vk.gammaABC[8] = Pairing.G1Point(uint256(0x14cf4a4efb357a4f1696f9165fb8d38a2e5c54a286a2fe3b6fd54d18f8b99036), uint256(0x11d8f2f49113343b5a1f8a384c579b518eab258704c32add30652e86bc89fc2a));
        vk.gammaABC[9] = Pairing.G1Point(uint256(0x257bce52522fb8ab52ba8b3eecf9311d935afea05450a15f8fc086b49f26f336), uint256(0x039dfcf5f35f03aee245b54a721a3de7dbbb6001bf2a8e9c48ca17d4afaea4d9));
        vk.gammaABC[10] = Pairing.G1Point(uint256(0x0dcd98909ae60a80d1959418439b0e8962a72882ca324ac7f38cd455654f43f9), uint256(0x0948eedcd0f9488a3b6923390fbc8ed149b2c29bb573d09a588637a5ad8d697d));
        vk.gammaABC[11] = Pairing.G1Point(uint256(0x2d30f621116582b2888d5d9727cb090e41d8385ab39d2264cb1e3bf0a12b49b2), uint256(0x1609505bedfbc536b98d236db926a1a3b8aee54d4c38f7e3c527e899f3a97588));
        vk.gammaABC[12] = Pairing.G1Point(uint256(0x13075cba0526a3a61e0487dc73adf929956a0aff91d8339260d27b37684743e9), uint256(0x0c7d6a49e044e124566ffc0313c17d43bddc6da80d1a8114eb99d4834b1c66d6));
        vk.gammaABC[13] = Pairing.G1Point(uint256(0x10e9b4162892309799111c15438c2c6601b6b84169a30af09442a1d276b5b2ed), uint256(0x01cd068ee94566caa7feb5364088fd250869905939315641208fb4c82f1d8522));
        vk.gammaABC[14] = Pairing.G1Point(uint256(0x1ee815cdedf8a4488d066d2c8302ace19f9764172968ac8a036fc947ed542f6a), uint256(0x14885b7425596ed662d8053339e5a06f096a73f9fde19143fe5bd3a8e219035b));
        vk.gammaABC[15] = Pairing.G1Point(uint256(0x198c51e2bc1045959325147707b4befca963d476d1847730c94fff4b62a31df0), uint256(0x0f5bf111bc58262d1ef165e856df88c6d321e2f15ad8642782efb73e4f4d9783));
        vk.gammaABC[16] = Pairing.G1Point(uint256(0x19a5c663b18bd73286fde091633edfe8ec0674e067a442b557c5c9a32e2946ca), uint256(0x05b9a23fcfb83b0abbba51af887693120f1bdb96bebf7726d9bd88f5536707b8));
        vk.gammaABC[17] = Pairing.G1Point(uint256(0x2e8d388e2d8efa004f76e3d48450b01841d57ef6eb7152d2a9b45fd331396a1c), uint256(0x2c70f73ec98c019994d07916ff6eab6ead848815c88607fe4c09b361493e9f96));
        vk.gammaABC[18] = Pairing.G1Point(uint256(0x20fcca68728140a804b340491810e04f58bbb04fe88c373766bf09fdc79d6a7b), uint256(0x14ef51ace20ff8f09e28cacf88aecd7b21970838f03a63b7e23c74c9d3883e57));
        vk.gammaABC[19] = Pairing.G1Point(uint256(0x0d8f4aefc3fb3a1bf1bea4202170cae69e825880f0caa934b2762434fb604d20), uint256(0x073774e95f3fd7a7a6a4822820a4d2e6ca217555bba4156c1cc1e2811969acd4));
        vk.gammaABC[20] = Pairing.G1Point(uint256(0x265353822a5f7959bba0b52b972992ab1a8499f0df06487749e2721d72f29ff8), uint256(0x11a8f6ee437a2036fc26f25d0e8ac038294e1e4085e2dddb05c2171335cb7c10));
        vk.gammaABC[21] = Pairing.G1Point(uint256(0x0f3bbc660b6e4587b87eade7c4bf067df946a0363e7abd13bc59f0b8921ee66a), uint256(0x1f02265a913190b22653f696e2a2884f74c6ecd798724f3884e8a3b329359564));
        vk.gammaABC[22] = Pairing.G1Point(uint256(0x2a4175e9b8f6d40d7d6b8b72cf9cbc497eeebe729a7e34683f6ce143215754f4), uint256(0x05adb444afb6fbc8283e18fca9f5c42357243a4deb5068e95680b5b4b52238ba));
        vk.gammaABC[23] = Pairing.G1Point(uint256(0x144746e47644350f652088e926bc8ecbc732cc9ce2891664a45f1ce9e27560b3), uint256(0x200b7ece3d411d798a605849369c75fad525b3ef6f1b88fdc018d940a3a8f118));
        vk.gammaABC[24] = Pairing.G1Point(uint256(0x10ddd31e1ebfefca5bed0f360ef2270f047b9c66d304df6b1c82b684c0029125), uint256(0x282b953ccd16a9e584847c18edd7ea012197203b8e86b6738ebe8897640266a3));
        vk.gammaABC[25] = Pairing.G1Point(uint256(0x2b192c3058da9267b3cc545f0c35e202e7ef193bb85f7d4293aa016eb5c71f1c), uint256(0x23123f3b4d0917335a7a3d8ffcddb0bfa347d53c2c3b973c904a39803ba65105));
        vk.gammaABC[26] = Pairing.G1Point(uint256(0x012582d844097c9adce1075d079dd44d53a8efc7f37724ac5ef04b04516de139), uint256(0x2291851d2955cc3c9370f188d217bcd798deef483be2ccc48ecc06112db8745e));
        vk.gammaABC[27] = Pairing.G1Point(uint256(0x1d71d3675f4f5801c778df81e6e632b350af73a4aa690736962d568baa8ae99d), uint256(0x00fd6d1e6d1915719415043d1a1bcf2efd1341dd257a631548d31d2175471b2b));
        vk.gammaABC[28] = Pairing.G1Point(uint256(0x1b2e37853d0b5d8cdef06b423993ec4db3699d8158061a13b0dc05919328cad2), uint256(0x00f47b9b86d19dcbc4f8319d118676fb90d17b2bb803311b0a68c746b9926b23));
        vk.gammaABC[29] = Pairing.G1Point(uint256(0x2517a6876da9275856f6f51c027a8ace1d89eeec69792ba2ef255be7065ba12f), uint256(0x05774559aa6c946ef12c928feb677a12007a74320ab001a5db087772602941fd));
        vk.gammaABC[30] = Pairing.G1Point(uint256(0x1c16b85a3fc08ca82fe5880baa6841594113514219ffae225c4af672ceb2feaf), uint256(0x03fc72d460264512b0bdc56cca9efb18fa54a78dc57b0d726e96259e152b4370));
        vk.gammaABC[31] = Pairing.G1Point(uint256(0x2ae46025bb07fba2ef7236a54654b32cc50d246fe05dc9e44a8465fd26d46fa3), uint256(0x27b473316a02f7590edae1f553a49e95e010dae71b5bd2fe629ab7447740571f));
        vk.gammaABC[32] = Pairing.G1Point(uint256(0x251f2c329e8f6cb17ae2053fa81f4b74269de81c14911feff557576c16860120), uint256(0x167f79f3240a462974e14630585692be515544eea5e13f4e3a39b57768633a12));
        vk.gammaABC[33] = Pairing.G1Point(uint256(0x0985664b49000580a45d51f58fabcf476e3c6c515b839a30b7d8eb5f9834c59d), uint256(0x08e201aaa2b26b57d19063bb5bacdabc46e0319a034431618f799a8186a8dcca));
        vk.gammaABC[34] = Pairing.G1Point(uint256(0x05f01afc67433ec7282efe42c34b22b64f94116bbdde1d7ad5ef731090aa2d3a), uint256(0x0e04b7131773af80099f98be6c68c2747b08f59724eac8e6bfb52c13b7e5a75a));
        vk.gammaABC[35] = Pairing.G1Point(uint256(0x1cc98a5f506534819c622b857d5c9f88790217a08804ea7dfc39e642c765ed38), uint256(0x14ee43a530e609704fdc0f8cd932d107c88c739b4c9e18c58519ddd7ec18c9ab));
        vk.gammaABC[36] = Pairing.G1Point(uint256(0x10aec29870a7013da8094354b897895002da48b79923fe108a2eaf1e6cc4a7b3), uint256(0x2f14a2a5b7f9d8ab8414de9512d062659843156f8340596d53d3275974ab63f7));
        vk.gammaABC[37] = Pairing.G1Point(uint256(0x06c42cf14b3e2d6bf9577f3e0fe170476cbda0ede50272a659202bbae135e3b4), uint256(0x0cfa7f509d7550a2f211d11b4430f95bd1a65d98dd7aabf3c4a40fd86c587cef));
        vk.gammaABC[38] = Pairing.G1Point(uint256(0x207df349671927450e37ab79d71b151bd4111404267e04df62be3b6f656b492b), uint256(0x0c38e90911fe0c9aff1cf2b2ffa01eb93d79d654df5f586da53a9e367699bac6));
        vk.gammaABC[39] = Pairing.G1Point(uint256(0x10aef87eb08306f9050a4702c843ca087155c482a1be73b0a2df5182412382ea), uint256(0x28bee4ba02e562b972f461aa32d4fab472fbc6df2163bda1b6b058f86e1892d6));
        vk.gammaABC[40] = Pairing.G1Point(uint256(0x1dd73042cad07d303e5423f5362e52f14f883309c24d7432c15466e57ca522bc), uint256(0x0504c19dfd78e3f1b79f9d58cead0480a28a50caa9fe65f74dbbd159cd09bdc8));
        vk.gammaABC[41] = Pairing.G1Point(uint256(0x1056e559aa155b1dbd40260b17a689220b895c0b86316724e781d1e3754c1025), uint256(0x08e884c522b4afede82012a6ac8bc12157ba4eaa0d2f3bb45ff227e4e3c80f49));
        vk.gammaABC[42] = Pairing.G1Point(uint256(0x16521cce3186e5ec64e2a31120686b97ada79f2e6f3708667db9a48a2295e622), uint256(0x2ea12eda292cb88990f4138132ed84a24b095bb0aedc6776f82f625233ed3489));
        vk.gammaABC[43] = Pairing.G1Point(uint256(0x0c390b49a15c5044564f836cc9e2f41d3279416a34c7809c1bdb89795e1a35ac), uint256(0x1a40c9d14294dced34023422f6d091e389e735be87fe47f6aab5c39f61dffe92));
        vk.gammaABC[44] = Pairing.G1Point(uint256(0x2c9c04628dafe680ef0ef413fc7b9a4bf7022734f6ddd6af2776a1072f8ec66b), uint256(0x2469b546cef727e178512a888779b9cfd38fb81df5f18523bd5281ad55952b16));
        vk.gammaABC[45] = Pairing.G1Point(uint256(0x253b993a9f5a2977e4ec284d0823d69a4f91290fb038fbd23bbf4101e5dda0b9), uint256(0x072b06ca5c0349bfa86ecc0fe1694cc8f68f95d3b7996761623b2f451264a1ba));
        vk.gammaABC[46] = Pairing.G1Point(uint256(0x05a04112a094536d44f647b95cc67b7b3a960678af36cd559c33c49a2681c373), uint256(0x156fbcd12b86b0e57416ee9c4f15dc061b6447e64be69831bc1df1ba6ab84bdc));
        vk.gammaABC[47] = Pairing.G1Point(uint256(0x2f651203a850a649083f61e9d6ac0175b24f3a120c1b88ab8bbd16c2216f603f), uint256(0x2247d7f0af4f8ca2d67a275e4862db4bd0f3caaa0227a2223d5c814d7a22c4ce));
        vk.gammaABC[48] = Pairing.G1Point(uint256(0x0fdfb5f1299082c81dc42f79bf1c34547a36c514f64ab34c5c92ce13279ab271), uint256(0x22231e518780b8f7431a77b7b9c81895aa3f7605d265134ae51ae3c98d225cc6));
        vk.gammaABC[49] = Pairing.G1Point(uint256(0x2430f48cfb210ac9eb61ed468842efa6bdedb68a62c9df9f60ec8871cf6cfa02), uint256(0x0bd7991232cbaf384540ca46250a9ea617cf64235c3cbbcfed5c9ea21dd17e0f));
        vk.gammaABC[50] = Pairing.G1Point(uint256(0x2e6298490f01a107ff3cd84fd5e650020468aa672e588c78acd4a9084e76f5da), uint256(0x234ff75c5ca9a3dafb73c25e6c71c6c8b8fb2fb75f52915697510e2bd719dd2b));
        vk.gammaABC[51] = Pairing.G1Point(uint256(0x2c8b0bc33126565ed5200e63ef60bfe2a60ec44fd3be8d9fba0cc2bac669d9a7), uint256(0x260178e1cb016a864b441e9fd8e3f580f54a69dbe0b2c910b862287d8ff1a4e0));
        vk.gammaABC[52] = Pairing.G1Point(uint256(0x2ce7ecaffcb6e5cb05f07e39b0990773c2efe1777ab258476d60232ecf9c71bb), uint256(0x172e1049ad61530baa3de6af117ce8f1d3db3d9ebcd92344624d4dd2d8fe66e3));
        vk.gammaABC[53] = Pairing.G1Point(uint256(0x157be777e085e7838907bf48e73acaf352710d3284ec8dfaf01403db424d9e46), uint256(0x1f4a9c88bd6e1e922a274679c896b3dff4c5ba2eac6b1d9778585764ae4c8972));
        vk.gammaABC[54] = Pairing.G1Point(uint256(0x2404dca978510206841b3069888598c13c4a674e17e912e793cc60005282827c), uint256(0x00a06d266fae402822f3af888cf55def91cdef1d6698809d90aa666faacb3634));
        vk.gammaABC[55] = Pairing.G1Point(uint256(0x0dd7eae8e2291e082c6995cadf103acfad99dd6276f59cc8f58a3ceb96d1af7b), uint256(0x2f337b74c84c796b37e88b6eec311bb9d36f396b8515f45ec7e8f91b0d35b214));
        vk.gammaABC[56] = Pairing.G1Point(uint256(0x1c0f4aca224dd570bc76b854b95d5ef1d3da5513eb1b6e12eae7b9e6fb63ff80), uint256(0x21a56c7e0f8f96e15c1ec3524d98e9aceeea4effeaad143a60c9eceb74123b52));
        vk.gammaABC[57] = Pairing.G1Point(uint256(0x2ad0fdebf4ea6fc604ae3a37df654a4ee5aaf81e7975243545487f3baf1f4d18), uint256(0x093a2ba775d905ecf2148a4d369a0d6a76c2bcd79a3110a0b2976e20c1bd00ba));
        vk.gammaABC[58] = Pairing.G1Point(uint256(0x221bbc01dc0df7f7c5ee1e9e12cc9d6ae1cdae1314cc18417a55058c0428aae9), uint256(0x037e63cbd407fb037cdf926248c55d00dc4e63408de674d5dde2105a3ff259c8));
        vk.gammaABC[59] = Pairing.G1Point(uint256(0x039a2bf33db6c052fbc9ebb291b5844c3aabc68e82ae494984aabf6e439a85b1), uint256(0x0ff49ca7616fc5f59a45836da5185ec5191bb608f94f2fed8312369388c98255));
        vk.gammaABC[60] = Pairing.G1Point(uint256(0x04f12d1308915a6010480e2d33046cc39e15cf76b8110690ca0e416893efb996), uint256(0x1ecd6cac6d3fa1b86cc328ddf9bcfe56c89a0865b98fda466348401c7f16b77f));
        vk.gammaABC[61] = Pairing.G1Point(uint256(0x1712d3839ede03ec02cb804c3d81ecb8035988b2971b78f71bd633fe435190b8), uint256(0x0fb454886c7d16129b0a52630d0f42e58bd7620c14e13c6c251848344209b655));
        vk.gammaABC[62] = Pairing.G1Point(uint256(0x2aac082db8e94ddbe2162aff25c1623df614b0b1196a8b8c3376eecfc8302815), uint256(0x045e5b672e50ab3e49eaddfc6a62f25bd6716f6aa31f402c10ed06fb04bcb31f));
        vk.gammaABC[63] = Pairing.G1Point(uint256(0x301464c03e9c13957c6ffe21dbfb39102e7ff7a551fb7744fb248b3a0f3aa8fd), uint256(0x07902c0b413d11f1096a2a1feb913add6c19db5dd24f5fdedd422c74b1ed77b3));
        vk.gammaABC[64] = Pairing.G1Point(uint256(0x22aca085e9d9f5bfd775a6d926a6a97b2fc3fbfc1c2dd3551757aff144824627), uint256(0x004b21f6efcd6cad47e4e6a3566a4769cb5f7b3346f248fdda1ac72cf535bbe5));
        vk.gammaABC[65] = Pairing.G1Point(uint256(0x26eeb45982dd6532747887588cdf399c26c783332147213a4a63ad0fc6da9395), uint256(0x28a032973f9a34fb7a145071b60b470eec977a15d3532318900e86f49382c531));
        vk.gammaABC[66] = Pairing.G1Point(uint256(0x293eda534d7b851502bd35a044523f20eb1ae4f3a9c87c573b837b19c1e42aad), uint256(0x28255cae5ae487cf283c66afa37aa3596c6066936c98ab1d89c3e1f5641a595b));
        vk.gammaABC[67] = Pairing.G1Point(uint256(0x1ac49c0f53050f10fe96918085f334ed5ccaca973cbdb2edb98836f642dcb4f4), uint256(0x2befd323890b3ee077a5eae7896ba7ed07f186dd6351c0c6b14a4f4498ee736c));
        vk.gammaABC[68] = Pairing.G1Point(uint256(0x091ffe0658476ed3d7809810a4fe950f9105404654a6cfb81e032a0a2cc82518), uint256(0x0a9264c8810866a149476c0e2f32e420ea6fe4b0b7ee8da248aeaa7b37043a43));
        vk.gammaABC[69] = Pairing.G1Point(uint256(0x16f9a43879f9531238549d47b86619cbd31e16a4ab94c0613b7bb19f9a1cc64d), uint256(0x2b0f96ed7a44d035bd9e4bb56d44aad81d64e4de6905165f621befa936527211));
        vk.gammaABC[70] = Pairing.G1Point(uint256(0x09934f2f5bde93618680e3d11b9b0369bb0fd83d5787e38510c6cbf140404d07), uint256(0x2ca54897c04e85a1d29ed2f39e40edb67f30438e2dd21759ee5bd2feaaa5990b));
        vk.gammaABC[71] = Pairing.G1Point(uint256(0x046c1693176a7a4c2077ae46617e6fe0fd6133ff1467faa650bca616852af431), uint256(0x2bc0af9fcdf9b09e42b2b6f68a4bd3dda6acc842b43c793a5b7c292d7bc3c921));
        vk.gammaABC[72] = Pairing.G1Point(uint256(0x0c5c22458c4a0a28077ae0580ff2c08fb5454be974c121c636838360ff26e507), uint256(0x1e7a828146ae6388cd462be71d0c501a6207fede02c071c6b455816a56ef680b));
        vk.gammaABC[73] = Pairing.G1Point(uint256(0x068739813618979a8261b7383fb116e8331010717e35ebc6a117abcabaa16c1e), uint256(0x2f95639b1a1ef83be63339867ef98fdc4841f42002089cce6edfccdccc039c27));
        vk.gammaABC[74] = Pairing.G1Point(uint256(0x1ec34bea6fe96e2bec1482b6cbe7ac290cef00b601538fb52dc47a06ffbcee2b), uint256(0x2b5e243e072ac8705b261a7d59ab8f73ba5f6b28b25877db053c981b203f07d5));
        vk.gammaABC[75] = Pairing.G1Point(uint256(0x1cfcba36f5608009a9a3e403b70ea4862341ec60d06ac66ffc3be9ad6e80c093), uint256(0x10d98fec96fc3d1578e6962080c9e08e0ccb71a89687ca4c61243a9ed507ed2e));
        vk.gammaABC[76] = Pairing.G1Point(uint256(0x1c70f66d4b5b2179c725dbff7df5dfa21164d6bffad91d0d5393ffae2540020a), uint256(0x14d12e19c9ada10fa7d7ebe015e388c78bd5440d269a8c00752ee780527edf25));
        vk.gammaABC[77] = Pairing.G1Point(uint256(0x0b327258ba1de59e47752a4f80955dfeb8d0a123277ff496fd4c26b90922f010), uint256(0x1a5dc2e55a6cae8d5f2a63a7097fa34e866864ac8080e2044b8f84911bd093cb));
        vk.gammaABC[78] = Pairing.G1Point(uint256(0x265cfae9c4437bc39fdfe3f0885220c1a7ac9b84074fff926aed5e3606e8dc63), uint256(0x24a305f808610ebd847ac561e4a1fb4eee73bcacee1be9b4f824fc5b40ae483e));
        vk.gammaABC[79] = Pairing.G1Point(uint256(0x1bbeb6d2ca5b21638c358081fecd07a0f799570077acf89d56544b72bba1bc88), uint256(0x1b705090a929f3ce34633135846484a031d9e7d8f2f01119b27275336b8d419d));
        vk.gammaABC[80] = Pairing.G1Point(uint256(0x294f6f2a0088fe7a7c0a63c7df45323112ad78ce62349152f2307acda5c6bc52), uint256(0x1c2d9bd11bf00619500868b339def95464a6691769284a5d7dc3904cb404e273));
        vk.gammaABC[81] = Pairing.G1Point(uint256(0x0f45c777c9f35df085b8697185b8c615dc818c717b327e1bcdca579fed2ffa37), uint256(0x1d3aa58da4ffa55fc08436b4098977d994dbad3389d0c6f697b8c2dfaa829355));
        vk.gammaABC[82] = Pairing.G1Point(uint256(0x177bb7ff7acbbe93245d63294288c15709a37fe83597b75b686c1c7508042444), uint256(0x16801ec50b5ae6b2ac04c4e0027839c59a276a21c3b8a9729480f4fd86169836));
        vk.gammaABC[83] = Pairing.G1Point(uint256(0x012ce8037f20ba1df66a41da3521ca3cd222f48dccc3a0a579e0ceb6418b5a1e), uint256(0x2ab609549ba191d03e3f968d1a6aa57312d591daa23f6d1dd9996b4de90f920d));
        vk.gammaABC[84] = Pairing.G1Point(uint256(0x17ca9eb2591439be4b39cd511e10f4d059ea9c5f51ff41f448f4d0116c5ca2a5), uint256(0x0f61fa48eda7a5463f83018649808a2b5c1208575adc372cb47d7f804bd23a45));
        vk.gammaABC[85] = Pairing.G1Point(uint256(0x0abaa74d288f8116a3d9beaa05e521a0faf2dfe8cdb1e65db90ae728878521da), uint256(0x05a64012ff1a9c1f835fb6ae89aa89d4c3560b7b5bad19565f1f3401876f6f42));
        vk.gammaABC[86] = Pairing.G1Point(uint256(0x0b9602cb6ffb90582ea57bba2bdeee554565f867bb66b87f312ca65b302f4a13), uint256(0x1fde9a8e7545223f2c3a10b31ab5a945441c04013ebbda790bd6eb42385c5445));
        vk.gammaABC[87] = Pairing.G1Point(uint256(0x21d8fe2ac727a6b080b18ef03739b74be3d3e4c57b37d548eb62a14d11c1b433), uint256(0x28bd534a18be2895c2f43c18f3cb965ff3bf8324c5c7e8442779843cba412807));
        vk.gammaABC[88] = Pairing.G1Point(uint256(0x249cfc25b283964faecccbc11487e2827cc0fe4b91f79d535c6ba754d404de9d), uint256(0x24af6e84a13f2ee13732e02bdb3ffab1cf36024c35a27b0c32c1bfe9bbd5cbc3));
        vk.gammaABC[89] = Pairing.G1Point(uint256(0x1b891655d12ee63d307f8b4e096f2c9a343dd508e6c2c61239cf55c5c37ee2a7), uint256(0x159d163d57b0ea56d1793ba4d57d921439840e7012326c03299ff51574489d42));
        vk.gammaABC[90] = Pairing.G1Point(uint256(0x2993e47e44b0f511ca420133e1b2a927e5b56484c2709dc4be98b8a4e470a3dd), uint256(0x1b75e2d29c903b1449cc61862a00f50baf2bc2af9442cef6ba55858f9df1634b));
        vk.gammaABC[91] = Pairing.G1Point(uint256(0x1503d79050a410f2caf5df0dde783e496d7997f3a9e6c46bdcc9408a5c8aec1d), uint256(0x282454e17741318cfbf14ab46addd3c0659f58165fc29ef2bb1b5a72da9b94e2));
        vk.gammaABC[92] = Pairing.G1Point(uint256(0x243e8719871bd190e146ea8d46212e5fec579b28d9a3755e306cf3a7ceecb635), uint256(0x1a7c17b0a51404cf11e47a9fdd721e46faaa59b69dbad75b012d5719272523ab));
        vk.gammaABC[93] = Pairing.G1Point(uint256(0x22e05ca6dbbb787b833e3dce7ba1d55da24f944b22da06da6644238a3341bce6), uint256(0x105cee7f83843a808acb1bad881cd5b2c9390eecebbbac465cf718cf2063b8f2));
        vk.gammaABC[94] = Pairing.G1Point(uint256(0x0b25e545b51d52b2adfe26d995fcf0de86e2374e251a5869ad436d5ea1ddbaf7), uint256(0x02607349ad9c7814271884b050de2547e016f57d466e3d385b65b4c6671298d1));
        vk.gammaABC[95] = Pairing.G1Point(uint256(0x2dd4ff8501832346d4a2964f010fcc81dd63ddfe53eddd0893c8fd8deb6c18fb), uint256(0x16e8eaa8be142e8e91288cbb252d7312da99f92ad18dd7dd01ba9b2c24c8e415));
        vk.gammaABC[96] = Pairing.G1Point(uint256(0x10360ce6949b961af2a3c98ddcbafa79eb3f36b5e1274a06da8f571daa06d91e), uint256(0x014f79270a43acecab916d55fb2e12aeba01a4c54edb1aa9f2055a01a71a90b2));
        vk.gammaABC[97] = Pairing.G1Point(uint256(0x0d8519559bc9c02580949231ec4b519531ef0f230c2627e94adcaa5907bac32e), uint256(0x0b61bb94c5d59c91659c1ed73056405e0b1f6a392be09840850c0cb7e71d2c71));
        vk.gammaABC[98] = Pairing.G1Point(uint256(0x1287fc2448588c8a46be2dcda6001fd1fb02a65593325d34158ef3418a2481f7), uint256(0x0dffc419f264a824f9dd9ecda4bc5048bc31ac4cf4f6641631f40d37ad31c5c1));
        vk.gammaABC[99] = Pairing.G1Point(uint256(0x299e44dad0844a5b6bb6e3ab3f41d70f43084b90b42b29695fc346746ff13b23), uint256(0x15f8e18fc3f45d8777e1d1af3ca1811e169ea5ab1712ac69b8ac746e4bc3731f));
        vk.gammaABC[100] = Pairing.G1Point(uint256(0x2189256af6f8de3c2cbb4a4df1c916e42de6dbe655355d8f8493d60ca7148860), uint256(0x14bd9099e03e90e55b1e5a6847f43cde305e7b03d646ce6d98fb5ab78f61620b));
        vk.gammaABC[101] = Pairing.G1Point(uint256(0x23dc60e39274fa379ed76ade46b91059ce56d7d4c22c0041778bc31623fc0f66), uint256(0x0c4c51df56e1ac90b8a3b93864d8904f5d2c0410a3596770753cf2084d3b977f));
        vk.gammaABC[102] = Pairing.G1Point(uint256(0x11bd09e8c2a31aa3980d42124188e8ccd5cabbc7f5f431a4512ae43878b4d2d5), uint256(0x15bff6e9e5dc2fdb6e8806b92756f2bcbe2b0f013a814b1f586e291e5f6f214e));
        vk.gammaABC[103] = Pairing.G1Point(uint256(0x1f47565c181bddf2440e2c54b368286cc940e5a761b3b581137a36e6e08fa04c), uint256(0x0908cf1abb71212be540c1445ecbb15718d381568b465817959865ad976eefb9));
        vk.gammaABC[104] = Pairing.G1Point(uint256(0x089ecaae7321cf4bbc8f546abb687429b3e7e1d6796576a6aaa7e1e5d7ec29a1), uint256(0x2484dea9d4b87040042e13f68d409faba9636ad8dc76339fa1203e565ed86e00));
        vk.gammaABC[105] = Pairing.G1Point(uint256(0x0c7ae82b6e2bfac790c6fd0d2020d93250bbc45a776649da3ad078504282d1ac), uint256(0x0a03854b7c0637514e1ac692783a244b33cb8decefa572a5fc0da68b2674bcf0));
        vk.gammaABC[106] = Pairing.G1Point(uint256(0x0c0a1049ac1a6fe46ad27af1711378fa0491e0af3e1896afc081cd95166601ce), uint256(0x1811aa845c9eaea9e351d3f6e1e81c1c7335e7eb17602a5d9048d1e777d8c57b));
        vk.gammaABC[107] = Pairing.G1Point(uint256(0x21d02c606fb10aff6d9db9199ead9d03d9a0163b1b5b4ea1f339517192a47fe4), uint256(0x04965e5b560f3dbc221122810fb25724cda59a46b5ab37636c7596f63a67692a));
        vk.gammaABC[108] = Pairing.G1Point(uint256(0x050fddf083cd831b17ad7aedd5827447a9318a07e6af54bcc7af1a2d66d21c08), uint256(0x0d996985134737f67324d92fc9b543d964509660df073e4dd4d913cc27c0e71b));
        vk.gammaABC[109] = Pairing.G1Point(uint256(0x0af1efa722dbef1bebc95b792dca99cb31ce8b28cd5613e13789bfb916241e23), uint256(0x1674d9f8396be117fb649a466ac13dfd43523dc833e8f7ea77660ba4b194bacc));
        vk.gammaABC[110] = Pairing.G1Point(uint256(0x081deb13c3276847bed54748e16069e6c4877271e7b0adeb391e66b81aa26206), uint256(0x0e028b853dbbfafec49701a21d344537a828d9c3f66468d9a2fc7f53caea4a48));
        vk.gammaABC[111] = Pairing.G1Point(uint256(0x03347cc7f820e35b8cb84ec955ce7b2f0db7a68a8fa2be1ea0e516d33f9b8ece), uint256(0x281a1d29df216977893cf11cf716fdf144c201f479da5fe2ff4a8cd83b5cfb6a));
        vk.gammaABC[112] = Pairing.G1Point(uint256(0x0644fb9d9ea9e470f533ddee3def218a92a4ce42e207de058876742975e30680), uint256(0x1d45d7aa13eaddd7e163ec2f316462c859e0095deeaeca9ce1b5ddaf782c4323));
        vk.gammaABC[113] = Pairing.G1Point(uint256(0x239e267d123508b4967aa6cfa645a494347c519aca947039fe1954f9aa13a0bb), uint256(0x13c70e1ece8773a08bf315d91cf079cdc49d7471d339a6441a92609b0349bd64));
        vk.gammaABC[114] = Pairing.G1Point(uint256(0x1555c0191bee646828e698d4c4050be658b772a2e87977116bb85e685398f78a), uint256(0x197f999099a9ed5602afdadf753936d73601e4f2c8877ac56665d8b73450b1ee));
        vk.gammaABC[115] = Pairing.G1Point(uint256(0x177051b6bf35bd5f7d3e6ed4e0d498058aff41d36e2a50e975ed66bcaa806395), uint256(0x100f15c5c0b5b48e703fa8f45b6d60ae66790b420d5c1f57bea0087bfef7a400));
        vk.gammaABC[116] = Pairing.G1Point(uint256(0x0b72849a9129da6a56a810bf9008bd4f53b6f1e954b66b1c6c17c9d609fb0308), uint256(0x0c1b0d46beee397d1ad5203d3cc5bc05a23038fef3622c965f90bda05def1730));
        vk.gammaABC[117] = Pairing.G1Point(uint256(0x2a29b06072b35001aa52e4da4875f34f2879634acc9ffb312727e8f01d5801e7), uint256(0x0f810e40c47152bebbd1d5d4d94f88ad22417fa5000d0df3acb9469fc30998e1));
        vk.gammaABC[118] = Pairing.G1Point(uint256(0x0a4af1b9c3d48d1b0d6d56d9e8f2169d6d12cac821de2a4da2105697f6218d9b), uint256(0x056dbb20e8dfef98a4c8ffbe8bae14d9b49b44460878d312114b17d5140228b5));
        vk.gammaABC[119] = Pairing.G1Point(uint256(0x1e53805b8f231b3537d4eccf9fb205863f2c5d52f85cc66e98eebeee2d15d225), uint256(0x2214e2d2b75a46fb7f9b5b5da3d878d6960653ec2b2b87dac3066225efbc6e52));
        vk.gammaABC[120] = Pairing.G1Point(uint256(0x1de73e4a428eb2ae983020bb8810dcc111d67e5c08782831fe5402466194598f), uint256(0x2f4bb85bad09d2c95b3910573bd52291bf02261afcf206d8038c9b5af6290a4f));
        vk.gammaABC[121] = Pairing.G1Point(uint256(0x02b4c615079b63ba5eeee8f4cd910a18cd32cf026a8ef2d70bf7d02d0133b385), uint256(0x198a5ea2f871e470cf23a0e17be4cd7896f290312b58ed6ad6692777c3bd0cb4));
        vk.gammaABC[122] = Pairing.G1Point(uint256(0x2a9bcf4985ea7e921218cdb6b4ea9aed0e8a8744e6649eb517db5b45c8efd82c), uint256(0x2858095a9a95c712842f15004eb8a9a57c1eecf8890802e09d072054165ee134));
        vk.gammaABC[123] = Pairing.G1Point(uint256(0x0a070abb8c91e92393bd8d3a4c8ed526d1a21c5c82e18f011e49aafb24f46735), uint256(0x153b78b139218324f74c03c174176599d01d8c5ad56daba0438a2bcfcdf9aacb));
        vk.gammaABC[124] = Pairing.G1Point(uint256(0x2eb181a1a89986a45c0f690e88498e5dca82fbf34482aeb2825e92a1818baffa), uint256(0x1df3a7451a7fd5543fa906f05ab2319164403f56ec675474f493bb1cc2690e59));
        vk.gammaABC[125] = Pairing.G1Point(uint256(0x053780ea52ca0df3d90533bb1094cb8aa1a2f64d7c3f3b12796c9f4630dffb2d), uint256(0x2d5305cfdff4be4e60fc20c39ab5b0bb9c0e5adce6e63b5e699b66fa3932b39f));
        vk.gammaABC[126] = Pairing.G1Point(uint256(0x2feb8addd05a36fc4f13135ab30ace6b942add95702238b99d33211bbd7200a9), uint256(0x28add308285dfd6b7b43159a249d3c4015ade0813689923257b06bdaa8fb4130));
        vk.gammaABC[127] = Pairing.G1Point(uint256(0x09330b4c91b5055f36706696bb46bdb8f9557bdb78e277b25d4bd8270a7b40f1), uint256(0x2702e7369ea721e6655b176183c23276bb1727aefacd4db7f42146c3c61e8647));
        vk.gammaABC[128] = Pairing.G1Point(uint256(0x2aa11bf3646d1cb0334041f51b35ecb879f5b07f70f765e0fca6279186434555), uint256(0x033f4240e398a4e84f05085e582a7573dc2f5a6e56bc3faef711861e943da9cb));
        vk.gammaABC[129] = Pairing.G1Point(uint256(0x1f16a84f01c2953d127120384e676337195dffed9b6ad22357605074d0339abd), uint256(0x1d4999aa4f6b8ea6e6fa2523768098d4be54145286368b561e22145f13798e82));
        vk.gammaABC[130] = Pairing.G1Point(uint256(0x17e8e9c14c22bcbc63508f172f12c00c4b774fb4f13d74bccea989be71661762), uint256(0x047914e10ec0668acc84a799e68036028135fe62ffb49dd9d6444708b8d76e4d));
        vk.gammaABC[131] = Pairing.G1Point(uint256(0x0834f797ee0c7ecd77f1fc43f68537e607ca771c617e5da2bdf12b98beef321e), uint256(0x2cf81741d0f036d9087caab8c8d10d98ad15532bfee242cf229c7cc35f28df05));
        vk.gammaABC[132] = Pairing.G1Point(uint256(0x2fdf6137975a56fbd35e71fb221853b2883b1eee805dd6fca93a41d28ff7c11f), uint256(0x2d5580a8985ba241c102dd3ff2b4c47b69015d03d5c09b87d012a4fdc1d99521));
        vk.gammaABC[133] = Pairing.G1Point(uint256(0x0c0db13fbaeada9580d2938a0f3148fede788809517cf1333d8c1dbfeead654b), uint256(0x04f62ee8262951c7dab731354d8300842e708a8b0fcd12994834d19887d76ac0));
        vk.gammaABC[134] = Pairing.G1Point(uint256(0x12db0587436d5f8751f156cb0853c502c907ac6446859973d362faf6bfcb9aed), uint256(0x28288d32a3f89bd69cd9fecf0c0c3289a9d8de83c0253c1e4ccf62bc3649dc50));
        vk.gammaABC[135] = Pairing.G1Point(uint256(0x245d887945a146795149e4c13e0509b643ba7b030307e4d27c9bb9c2a09f6dc2), uint256(0x22c1eae79c28b4207bd649d26bbd6cfd7624eba649e800f8d2f8e5e4c7616ea8));
        vk.gammaABC[136] = Pairing.G1Point(uint256(0x1b7fffb84afcb612ae0d2e8b93f8277e6f10bf27ea039e9b3460e8e926356cb0), uint256(0x22a364713e36ed32be2cd4c8eff53774a337b9fa3b2b755e8655b0e3aa96d926));
        vk.gammaABC[137] = Pairing.G1Point(uint256(0x17713bad0ca31ec1e8e8e370257ca48ace78425425c068955fc74890aeaa533a), uint256(0x0c71202a4faf24469fe7a046eb84b3a69d2ea6fcdf688ef9d1d42de3ab78f47e));
        vk.gammaABC[138] = Pairing.G1Point(uint256(0x1aea62c9e4e8531db3bfad32a51d5e2fc84ffb437085381abbac5ca63e12acfc), uint256(0x11f080e809990a9d2dd96b9dd3eb38f54cbda27e67e04e746dd46a745aeb6f53));
        vk.gammaABC[139] = Pairing.G1Point(uint256(0x0e5c3245c670a85c63615cb4023e1b8012f3c4862476fe8e18788ed4ec713a0a), uint256(0x08f59924874c92081957f912b13c97f1336a6afc34a0b73d2c37b17bf539a472));
        vk.gammaABC[140] = Pairing.G1Point(uint256(0x0455c4f299f0219a4050cf541dcb1e76647dab2c1be73e02f3e4142548841b31), uint256(0x1e331ac9cf09ebe8fa28e5d603985dbf91113c12f612b40fb385bcc8c12f2ee6));
        vk.gammaABC[141] = Pairing.G1Point(uint256(0x1701e4048312622c23468d2bf08bcb4626dd2cdaf1c9355ca496cf9a99c4227c), uint256(0x2155d067104d0dede52227657aa811b7c1169adfae3e95a142b128e9ebc49895));
        vk.gammaABC[142] = Pairing.G1Point(uint256(0x2e7f4ec937ff9071663dd57e6d1041f83cd6a5f1ed9d1fd0e6813d7d7c5647e1), uint256(0x19bd01983736068d337f00598503a9e3fae6423e19a99df57b0c6bbb74e92440));
        vk.gammaABC[143] = Pairing.G1Point(uint256(0x2155f5bd5ed798d7a30cac6078a5ba1108767264e23f23785ccb630451cd8ccc), uint256(0x0c67cce34fda4a6766c1ab5b0417c56ec854ec3a37149b8c08bc58327b0c6bfe));
        vk.gammaABC[144] = Pairing.G1Point(uint256(0x26e022b1111eef2cf6d42f3c1ff155e79ddab1ef2a514a20bf0bb8d85d4b75b1), uint256(0x102462458de9d4b89f4b4e0a08de2066ec2daa37c66ed381ab73b3e66f315290));
        vk.gammaABC[145] = Pairing.G1Point(uint256(0x1f4da386157a24d894de6c57bc055f045193020fab3525e376f2e8bed37c7358), uint256(0x15845c4619b365b9b319c19f144a4b575f441a6f7eceb2b83a93c1f76af747ea));
        vk.gammaABC[146] = Pairing.G1Point(uint256(0x235d8f44a09bce8d3b9461f200b1649269a395d6ff664c64b2c49ea06c871075), uint256(0x143e15883ba817d247570f9ade024a6e75149cddcd355f219a7fcf207bfa5629));
        vk.gammaABC[147] = Pairing.G1Point(uint256(0x075696f360701e02f5b725edde3d495c9e2a69db5a08fff5e559dd7a29819570), uint256(0x2943ba5aaabcea3aa113a0f8ddb03121514623537d5c3be6a9d5802569702ae2));
        vk.gammaABC[148] = Pairing.G1Point(uint256(0x1a442d8af20ce79aafbc0d5a99180fc18d740cf1f283e3cb97354933e6eccf56), uint256(0x0381d6e543d69c1924ea2fdd5937accbecdf8e7a7ce761a5572744c55055b67f));
        vk.gammaABC[149] = Pairing.G1Point(uint256(0x3033be493f87df7e52591aaecea92bdc4539556d03ff0e377f4f1968069f2f03), uint256(0x0b36800fe089cb01d14532f039c8db9560fecf83099ce88e237f4daafd6b4be8));
        vk.gammaABC[150] = Pairing.G1Point(uint256(0x14a0d2ce7744204ea89637eea68b2e06b915620da290ca8bd698cb5bb9d0791a), uint256(0x24cf2629f7ff9cbaee9616ce153104d1cdb74680e366bf586d8805a653d7b577));
        vk.gammaABC[151] = Pairing.G1Point(uint256(0x09aaa398197980ef1e9a0ac8fb845b4285679cbc94fce77a02b17500b02dc936), uint256(0x26a22a22f70fa4440411504871ec9046fc6a256781fc2367768e9eed27001f56));
        vk.gammaABC[152] = Pairing.G1Point(uint256(0x1953100377833ff551af91794462c425ed0341946891436a8ea00ad68a6efc37), uint256(0x1f84138b7ce1dcb2ee0e330a5ab56c41929c9f8215c79a1a0dade1cf0635e8cb));
        vk.gammaABC[153] = Pairing.G1Point(uint256(0x0ed950423d126452035e93011a58daa46937eafe6fe8e1ea834087309ed9d5fe), uint256(0x2efe08b941bf878074c334b9e4fa27d2c0e833850014acf229952fe2d93474ac));
        vk.gammaABC[154] = Pairing.G1Point(uint256(0x1edf9a20ec3dda1449019d1c551f6ffcbc9c85f5c634cb5c4bb46fcbed57a916), uint256(0x154f8a500a41f0fd01c8ed17d5feaf28c5e0fa489e3b94ccb955acd44a4b0e43));
        vk.gammaABC[155] = Pairing.G1Point(uint256(0x1a6ef09a8ec66d310348f6342a90e53ee6675b62f0b733f609b3115e07f2eace), uint256(0x194749c95426fed7d39b25947faa3bd0aae6b150e918099be77a09a7996b704d));
        vk.gammaABC[156] = Pairing.G1Point(uint256(0x1c63acae9d1104cf1fecbf10169841c25d23c3c92a2b10722cb9d56acdc46661), uint256(0x2992e7686260ca4ef9a07e9d13a7ac174ac6bda5ccacf10f745de78a318d3375));
        vk.gammaABC[157] = Pairing.G1Point(uint256(0x084fb9ee6bea0f06d366c304f97630f3f54f979a8de8781512d0d27a791bcd6d), uint256(0x29c649b4918bf3983cd3c5c266ac41d3e753d5eabc49d53d092992da47f3b421));
        vk.gammaABC[158] = Pairing.G1Point(uint256(0x083f50e66c70ece2611b9c3ffa7ecfd232af2ddc88cdb5f0a60f252acb44794e), uint256(0x0a2a74647294a4d442ffd78764a3f63ec7dec617f43809f0bc6a8ebc764c3623));
        vk.gammaABC[159] = Pairing.G1Point(uint256(0x2d56176873d624797cece2ea382e8c3d9fd23f8834130516219938d97a2c6b3a), uint256(0x15dd1b39551652a39a38729e061051d7ce6dc5fb605f65925d1ad06639defd37));
        vk.gammaABC[160] = Pairing.G1Point(uint256(0x2bcfb0c5145f6ae40c7440e83018dda69d21a9c71e75208a443e00c47dda3abf), uint256(0x1da619ca58e3b28b5c6aa238cfe18c602fb8c84d150ba7b0eec3342cd50c06c0));
        vk.gammaABC[161] = Pairing.G1Point(uint256(0x16ca8000ab49b94c20b35941224722fed7c645bc4b8c32e83ad9dd86d59b6a1e), uint256(0x20b3418dde6eaf6216de3de6c284b265c97138ee31b86f2b64250aaa0c147d89));
        vk.gammaABC[162] = Pairing.G1Point(uint256(0x23ca3433107b69e102f583294f75376ec593136ec9c22af171d50438c801256d), uint256(0x1e0203506fb431a034605a4744f0112fa9f00674976ea38cbf12d2f9ed18ac01));
        vk.gammaABC[163] = Pairing.G1Point(uint256(0x1995e6d00a6b037056258fbd239efca0d6c46c8648f50fad1894a8c58d062620), uint256(0x11085a27537207d2c6d2a7d274fda0cee8de98dc6284f2fc81cef9c3bf6d18bb));
        vk.gammaABC[164] = Pairing.G1Point(uint256(0x13908a21c202cc7ec4f7c1aadf2721b981a373c6886b1e93d2387c501a520bc4), uint256(0x1728f087b432c8e727723195a6031373097072f783c59f41f4e6793dc17622f9));
        vk.gammaABC[165] = Pairing.G1Point(uint256(0x1be697ea0cc301b34dcfdaf44b632affdf79bdd2b15ae89810b0cf81fb86dc16), uint256(0x0475bb8af9cb294ba781b99d9e4c29de1aaa5a69d9d66a8ebbeadadcdf760dc4));
        vk.gammaABC[166] = Pairing.G1Point(uint256(0x0261e38fd747f8e723a8289d2ddbc0c551fa21e0dafc78977a56965dae92e40c), uint256(0x0e3bd84e29fce130acaf198405e98f0b5314554190ba2b714b43133711be7fd7));
        vk.gammaABC[167] = Pairing.G1Point(uint256(0x27e5b3b1a0a1b5fc346890416867393411afcb87dea845ec61bbca89a0f0c852), uint256(0x0ae70156d8cef66529e81de88f3f211d0f0fa98f2e6f4ea47c731b00142c5a76));
        vk.gammaABC[168] = Pairing.G1Point(uint256(0x25b2ba4df96b04ffa969d71e2b93688219b05e60f192a782702fafcf8e5d9532), uint256(0x0f1afc088e0415e5388d542ced3ebdbd93d3be0854e643d11f02fb52194ab62c));
        vk.gammaABC[169] = Pairing.G1Point(uint256(0x0e095271797244a66dd6c50130667032498e0b68ace8979468b17d57a42a0dfd), uint256(0x1f776f0d38a464db5f53f823c024b0bd8c38d23110955e33a843e77e81d38ddc));
        vk.gammaABC[170] = Pairing.G1Point(uint256(0x0580055206af6b3a0ba9161a5bb3644dfc4ff05bcbbeb0ca77823ecafbfa6ab8), uint256(0x2559a143d7122655738b0dd0e87c154f97ebdab6add99827564b117337d804cb));
        vk.gammaABC[171] = Pairing.G1Point(uint256(0x2127a9bbabb303f8438783270dea442ae6b09af297303582f3d735a29f69721c), uint256(0x027a8099133beee963b03280d8301af3776b224adde11b60371cd49514637186));
        vk.gammaABC[172] = Pairing.G1Point(uint256(0x25a58d46cae8b4f90b8ab9303876395f3701c0ce8356e5cce2fd9b827bbbb939), uint256(0x002af7bd572e403e57377ec641d19794911124476a87eaf627ca42c899190cb2));
        vk.gammaABC[173] = Pairing.G1Point(uint256(0x1c13a1786ebf9665a06b10434a2078d1418e08f13a0119a60ad02347208c46ef), uint256(0x1aa3bc3e9b980075aa2682224c380a22795262cb8b34ca21bd274a48be1527bf));
        vk.gammaABC[174] = Pairing.G1Point(uint256(0x022ff9a603b42703a516129360bbf11de3c900337b93fc8ea8a025592407ed6f), uint256(0x20389df666b8987670da451caa07baeea0296b466b374fcc27a252f0275286fc));
        vk.gammaABC[175] = Pairing.G1Point(uint256(0x2565fd6f24d6bea972fc401c1818c883818528dda00922ebed43fb58cccde456), uint256(0x247fa56bbf73eaecdeac9efd110427f877e42130eb361d650b3285cda07d4715));
        vk.gammaABC[176] = Pairing.G1Point(uint256(0x127d98400ca1b5012a87247ada8d4e87ce8e03a203f1a33d4074af941d9cddc6), uint256(0x1d2ae289a8856e1a3dfb27a270cee7b88f924ad8c916f4886eaee8215517d4e5));
        vk.gammaABC[177] = Pairing.G1Point(uint256(0x1be8b7aa0eb9ab236989833a0304f3db0ca879cb505ba635109f0d069092c83a), uint256(0x0aad7585a022a91c9076cf5358e7c5f55d78060515a1d3d6fb26767a579d19f2));
        vk.gammaABC[178] = Pairing.G1Point(uint256(0x098d38e89a549b92fdfd686359f497f5eb4430bbf561c01cd1f4f01c8143a667), uint256(0x0d81d7b2c9b03f00ad1d3800e4930e8ff15e4e1f79047c247c55562a321e2592));
        vk.gammaABC[179] = Pairing.G1Point(uint256(0x00ce8a3eb65f7d5265f8fb5b2aa21493069e575c4e0c52e265a84e00e35bea86), uint256(0x0b8e9f0eeb8704f96c3d3ed6449024bb012196e4cd5b6079eae6f11550a7630f));
        vk.gammaABC[180] = Pairing.G1Point(uint256(0x1d8a95c65a56fb6aca3fb3af4822379b5fc69932af05fdaa3679daca585ca369), uint256(0x21f71f9421b7b1ac4237edb1c00fbbb7218f32d62bcaf25ca1c6fbb223d75731));
        vk.gammaABC[181] = Pairing.G1Point(uint256(0x1f72e90633ba2196a7887a58b07f9a412be607c05020c2ed1c76a86fd665d48a), uint256(0x0aea99a3537a46cc352188073e19484d091c9efb1eb5642bd02208f9f001e809));
        vk.gammaABC[182] = Pairing.G1Point(uint256(0x29d13f7f84da97933f1790408713dc64ecd503988ac3ec2e92f6f37eca1f6a66), uint256(0x2e9717375fb526eb260498b598fd5a4d2889da2f1587a5827d88ec0b9be86d24));
        vk.gammaABC[183] = Pairing.G1Point(uint256(0x2e3cc5c18399a4adef70ecc4f822808ea5a084509777771d4d65a7bbdb2358a8), uint256(0x206adad59d4e2cd68393d2295dfd5277268972c436d2f7cfe244daf7634d2d9b));
        vk.gammaABC[184] = Pairing.G1Point(uint256(0x003da1d8bb69d1bd118e76e5525d4ed8d1373f3566e6d6a1eeb39ae46fdb3472), uint256(0x07a5f233436d57a8a494eebbcc5522fa7b87b8c6df2de8a4a3cd95ec3202e1cc));
        vk.gammaABC[185] = Pairing.G1Point(uint256(0x05dfc355c0bca2c48a9c89e62101ca294bfb4569fd725d3533d3fc8e2a00e20b), uint256(0x2e76188ba521a8c29627cebcc6d83b61c85a359506f31e19db18639cf517f3e1));
        vk.gammaABC[186] = Pairing.G1Point(uint256(0x1ed2f9c8e41cb3f22c59b0dbdfa33689fcf3a86c1506e75e85d93e636fd5640a), uint256(0x26588fde874fc766e3c0cd1de31cda290eddc926214b738dd62e7ccb69912b82));
        vk.gammaABC[187] = Pairing.G1Point(uint256(0x1837fe1d06a11f308528a9c4102001231851155bad4b86d0c739f9bdf0507bcc), uint256(0x1282902f0c1a0010b74937ccff232c229738131582b7700b1493b55148a80c5d));
        vk.gammaABC[188] = Pairing.G1Point(uint256(0x2bd08513dcfa101480f53cadc65c0adc82683f0900d0f1d2202c8e7fa13cf94a), uint256(0x14b2fef9eef43c6e904a1d8e505860daf1c3300a0722b76d57ef3e6c24f0134c));
        vk.gammaABC[189] = Pairing.G1Point(uint256(0x1b122eea85e798b90db924c9ac6dbd656c7e52f0fc62721a1213bfc1c80e0d70), uint256(0x248c9661a174b53e03ae0aff50ecd855c08f47049cb0071587149c8d669dbaa8));
        vk.gammaABC[190] = Pairing.G1Point(uint256(0x29fd72049839ede7f4802a78b0adfd70aeaa411cf95a9a88810f919d65656613), uint256(0x10fff55ac1443a7b4b726478a74a4592c8c39f3932150d517d1b16b77cfd0e69));
        vk.gammaABC[191] = Pairing.G1Point(uint256(0x26db1b888a09c29facf8e13c0142be2583b4ad49a5dfdd7eae02d463146977ee), uint256(0x039db1855b080fe6a85844be85c6934d461ab59bf63aa3a6ffeef8e8d31fa48d));
        vk.gammaABC[192] = Pairing.G1Point(uint256(0x060c919d8f3d7dded6d431b0fa0d68f9a3c0544bdc414d65edb76a8ab80f185c), uint256(0x13e4a39bb15b2049138a0deda4240689e3a5b6233d46ee8f17c7fa07c7a47758));
        vk.gammaABC[193] = Pairing.G1Point(uint256(0x061ada6a00277120c9298c35bb992ff90e228f52468f286b3491ea61f401d0ed), uint256(0x228e9261331af0e5108f424b7b4aa1ed37e4aac6dc69d9636a5cd4aeed62cfbc));
        vk.gammaABC[194] = Pairing.G1Point(uint256(0x2543fbcedc2aea8cc88cbd55f0cba6bf32918de083291240d7a5abf27a42c0b3), uint256(0x166274d59664fde0115b6c24747392ae52596719279353c6fe03b571ad5f8733));
        vk.gammaABC[195] = Pairing.G1Point(uint256(0x0ea1dfbf492dc0e4050e6a5120c027f4cf39b34854ba0d051e27c7bf0dced192), uint256(0x143f4a5f450dd73a76b69681f1b9969dc775ad70170be47301d5ab0b0633454b));
        vk.gammaABC[196] = Pairing.G1Point(uint256(0x1f8005d0659a6386b75908d12329cd143f90271c18ebda56118d9d4d38ff1bea), uint256(0x2ed1c6e3e284f495672992c466b93a00534b84c4b59ddcbf454c0dfb6f64a61d));
        vk.gammaABC[197] = Pairing.G1Point(uint256(0x2f2e5dc6f1ff9e609a7f658e9c30f42c5878882258b87d7eaea766cfa39cc27d), uint256(0x28fe67d7b69ec46332797ffeacaab3480987165f6ea919cc6a58db2c75da0024));
        vk.gammaABC[198] = Pairing.G1Point(uint256(0x13c04ad7c64cdb94d34418a72cae4686ba9e06d7a3c1c69f5dc0bccedfe0bb86), uint256(0x10176d96c836f58d8e5656f5aae3b9aab70ca0ab0287a65a7423a410d63fd4c1));
        vk.gammaABC[199] = Pairing.G1Point(uint256(0x288eeabd04d30339352fd08343e280462c969e706763682bbc399249eda537cf), uint256(0x28d98c0e6485e2cb323088c6925e9a509157494e3ebd966155b9baf7ff60225e));
        vk.gammaABC[200] = Pairing.G1Point(uint256(0x14abf4eca0f34fc4df33c77f2c92f7d908b2f51ba6b7db249593401e97a96c26), uint256(0x242c3a52d4cbe06fe0063a8ed13db9a61e918e89da210bc7737bb891d928a2de));
        vk.gammaABC[201] = Pairing.G1Point(uint256(0x2bf98b7c8f3c4d3635fa2421dd0480290992d0f32b28f5692bbd10ce2c4dc9bd), uint256(0x1a370ac18d5aa4611751215d5e55539da7a4814d9d62b39f91901adaf076fab1));
        vk.gammaABC[202] = Pairing.G1Point(uint256(0x2b8e940e2d3fddce37c2ee77ce58a733a52d67c0c94973436cdb1a609f9963ee), uint256(0x029b71cb2eb5d25fefbf3a49df3c7e05d30ca0a6acec44b9e5b7bfd33406ab4b));
        vk.gammaABC[203] = Pairing.G1Point(uint256(0x12792a9a8e342bd75bf8cd7021511b424dad6d4a37908942571b4cd36cfa3daf), uint256(0x1461c7f2d161e9fcea4257faee0a6b4cef3f438a2ac3eafa0064afd6fc5e8e03));
        vk.gammaABC[204] = Pairing.G1Point(uint256(0x122f3ff82eda11f3a8e3b0b8fb8f2aeafac0b681d4277dd33488f3210968de00), uint256(0x096ac1f90ff9f158cdc694a83abc1c8f79b6fee10bb3ba231638f2b5dac1ea0a));
        vk.gammaABC[205] = Pairing.G1Point(uint256(0x1571bbc1b82af1121d2b2c6b69b59caf4de25a652efc8001e947929060eb08f7), uint256(0x0b0a9156a5879f463eeca98aea3c61305a835e89b259a7172ad36d66acf5e099));
        vk.gammaABC[206] = Pairing.G1Point(uint256(0x2a6c1e84e5eca2d6ebdd6a9c0612c1dd98c46c58f1bc2aae23703bf436cc3e6e), uint256(0x0233ded15e453af09119bf5ef6fb74e9eebc481a4063332f1d7e3a3b28aff31c));
        vk.gammaABC[207] = Pairing.G1Point(uint256(0x266bb82717e009533fe13a0a0f7125e017f3ed43a2a291558cbe5c755555a49e), uint256(0x22ea2bd44826ac58fb4c39d43043148fec343882fc1c6ec13b738e9f46b85fa2));
        vk.gammaABC[208] = Pairing.G1Point(uint256(0x0a295000363eb9e1e77e8f79056670096debe1442a84243eec80b04886778b14), uint256(0x2c355040ca19b5c5eb2a7033d16ebe3efa46c214800b7be9403aeb489cd2729f));
        vk.gammaABC[209] = Pairing.G1Point(uint256(0x277af93bb9eb40b085cdfea567aba10fa95bb6c5715aad3c84f4be482a9b78a4), uint256(0x2ef2382aedfdbdfcb9ad1cca152d856c7cac1baa0b393ef94682f2eb9bf8f957));
        vk.gammaABC[210] = Pairing.G1Point(uint256(0x2a308a3dd9e4dd37e0e4d8606a1d071e5beb73d609f04aec222478b8f91c9590), uint256(0x2d7feada366909678362f88c1a92c578d69c9fe14eb72e741100823d692817a5));
        vk.gammaABC[211] = Pairing.G1Point(uint256(0x0132c64e0d23c8d22660ee09540b5eee0836d4324a990cb57411e389d228a05f), uint256(0x0c4338856fcf14f551f4e05b70a18b8e26a06dd9c412ce4c159412331c997518));
        vk.gammaABC[212] = Pairing.G1Point(uint256(0x14790902ffd734c1df1a84f144b17ffd9dd332bef89b09764cb2ebc595ce3380), uint256(0x1c44b29e8d7507e7d1651b2f09d78fd0d97f651071e0f9be3d609dd8bf48f533));
        vk.gammaABC[213] = Pairing.G1Point(uint256(0x1da991017f2a2a381f8de708dddfb0b007203d809dbb2b7324154b70d62d8ebd), uint256(0x0c72dc86caf1092bfb8a076fb1d2688f6b3528c7244e6c49f9fce806461f3c34));
        vk.gammaABC[214] = Pairing.G1Point(uint256(0x12fbdeb1bc775d430ba195680b349bd1a842268220ec0bf60ed42d295adc5372), uint256(0x1563d465366a519ff81fc94b28f8c516120cb4526888f78189a3802e1bdcbb36));
        vk.gammaABC[215] = Pairing.G1Point(uint256(0x2ef01d457d83930fb54e53b07a8b162bee80dfe132c0fd1016725a484d2f13e7), uint256(0x0a54691ccf86a83c300d8412110b0b8b20fc2e2eed37ff99119180ff1b3b9723));
        vk.gammaABC[216] = Pairing.G1Point(uint256(0x0a5cec4a63a86da0176626942b13ed87827195452d9b8d818fca23e9853f5a93), uint256(0x1719ca46b75d612da90e3b5e36016384c8ae34d73286ece6eab18502d6d2eee5));
        vk.gammaABC[217] = Pairing.G1Point(uint256(0x18a8080ced35c02e47b4bffadd6ec3c068273e8fc3581a4034e67f6b0e456b72), uint256(0x0c9db3731be52aade33b3d1fb9605b38863dc096f68dd7d31369b70c5eabded6));
        vk.gammaABC[218] = Pairing.G1Point(uint256(0x298e72ef5faa48ce4464e62b478520059a19ef44d1210cdec1a6d8514ab236f6), uint256(0x2340f953fa214c4712541853c3c0e8efa186873e9e743b2b6e43079b775fd030));
        vk.gammaABC[219] = Pairing.G1Point(uint256(0x01a26b374048f6359c894092868a867d8a4c1a49dd9fff4b1fda439d33f9de16), uint256(0x1059052f2db3e6597db5c12707f4c87d568a6e1ac74d2af18db84f936065d9dc));
        vk.gammaABC[220] = Pairing.G1Point(uint256(0x28d8f22685541e69b8f6b27615ed06848f5541bb5d754be1213cf21dc91b3b67), uint256(0x0e6b4ea3f8f25e55a809b88120d1d3b77455a4541be77cfb7ddcb663888dd448));
        vk.gammaABC[221] = Pairing.G1Point(uint256(0x2bef076451fe6588d4623c1fdb913569356cdd4c8502c5eb50afc586a971713c), uint256(0x170a75db87e89e1e8ddac4818f245cc2cdb14ef2f332f31c49052c53254c8434));
        vk.gammaABC[222] = Pairing.G1Point(uint256(0x00b2267ede5cc10c03470d540953ffdea9185a629ec1214a675eecdea3998d2a), uint256(0x0e76db502a8b0b266d2bdae9d5dc85e7140d33c5f656890a46062117e5a54a1b));
        vk.gammaABC[223] = Pairing.G1Point(uint256(0x2d8aed4798711f4dfee2bfdf78c3ffad48bf12f56e63e91ef50af79f39be6570), uint256(0x018dae331970ee92b5e3f78f37476d1c2bc19671519c35db792f3dd4be39b91c));
        vk.gammaABC[224] = Pairing.G1Point(uint256(0x2486b0a7070b4f57acc6b13dba636d5b36968a44cd05a1e7f955153a176d64f8), uint256(0x13891f1b70706d31e6828cb8e11c684c186eb30624040ce5c9e9d4a8975b422e));
        vk.gammaABC[225] = Pairing.G1Point(uint256(0x2a82020fc5e765f93594b963cc98919478d9c2fe9e5094b5b0791725bd5b9894), uint256(0x0aa74767169eda1814fc4f6ca8f9fb6fa0777f7de3eb89ae667f2598549f4b3e));
        vk.gammaABC[226] = Pairing.G1Point(uint256(0x2df8cde8f881fa59966f992ba820532f7ed603da08bc02aada01204c31c4ab9e), uint256(0x0516b7776e6011c7db2ea62c4131969cedaf8000dca664f46bf010d823e90421));
        vk.gammaABC[227] = Pairing.G1Point(uint256(0x01c5c35efe204a4a757ce2c76736b3373fe97693558d5f1ad8c86baa6817ed6f), uint256(0x186809c97eded26bed8d8c38c4721a51c3c87058f317d712c7dd0cfe23f7e526));
        vk.gammaABC[228] = Pairing.G1Point(uint256(0x222b5cb588bb7ddf0be768f2882c5c5d5141d61f18b62a352393643ac00ce9e1), uint256(0x2a40c68942a79b7751bc2d323a5912d4d2cc7ad24b758e9302d9439120a06b9d));
        vk.gammaABC[229] = Pairing.G1Point(uint256(0x084825d68a985a0bb47ff483785b00d59ce05fe405257e68ac4f3458041e4ddd), uint256(0x0151e125bc2b66409e5fb352918bfafb7f28c63961b2f98f89bd59af29fcd38c));
        vk.gammaABC[230] = Pairing.G1Point(uint256(0x05637a5e4f12fa58833e1732ac3218aec5d5c7d2de00cab1061b0c13d66dbdad), uint256(0x15ede44b2a924024da67ed319dd9fb8b1f4bad421b74b06a5c96a242c821d2b0));
        vk.gammaABC[231] = Pairing.G1Point(uint256(0x1a0ca11d81d8d4f9f66aa93b57cfbcfc2dc41bbada8aa80b418308ce6c1875a4), uint256(0x060d4c3f002c04eae1702ddee61ad965957fdf8eb0c657b3243eae325c9ff4e1));
        vk.gammaABC[232] = Pairing.G1Point(uint256(0x29da69aa5e9de77e178b6a0b32e07240ea2373fcbf337c96de2c1b4236b06095), uint256(0x1c0ac14c02163dc3330f8be3e274c7d318ab97dff1b70fca6e9cdf6196233f98));
        vk.gammaABC[233] = Pairing.G1Point(uint256(0x0d8603390033d9f32e0beb6bff49b277a461bfaa37462fc95751a1d7d50db764), uint256(0x0af971a0e1d9788ea3f5586269263406f1523280a7c27aa7aa3e9c0ebfff6a4d));
        vk.gammaABC[234] = Pairing.G1Point(uint256(0x2cdeec136b275e9357be1467cf3cf93c320d5d99bd0d928edd3e48b2226f280f), uint256(0x224a0e719a5f039cb3493601f38261f91faae1bb1e2751764c0fd78b1b979029));
        vk.gammaABC[235] = Pairing.G1Point(uint256(0x0feef0e40c75e4ee004098f481831f7d71964166988e394367b3a328a0c7ab7e), uint256(0x26f26a064f007464794f107384f113d39eaf7181a105b1088b5889f46d2ce654));
        vk.gammaABC[236] = Pairing.G1Point(uint256(0x1f313901e2488e2ceaabe94cf0a2083c04d9b216feb414db8773ec5e25ab04ce), uint256(0x21e571e458786ce8f81b45f96203ca5c39330a9f9961f59bde9794cfde27fcf1));
        vk.gammaABC[237] = Pairing.G1Point(uint256(0x18493db1a03ef7cd4e14c497fea81582ad5c04df64aafc104a187030c756cb27), uint256(0x2880a0bfc3f0180e66b7adcb1b97c89b1fbcb5f28227a9635ec1ac81771f6de5));
        vk.gammaABC[238] = Pairing.G1Point(uint256(0x0d7e7ac682744611860947387bfe84b6610cee1a89fe46a70332b98ce65e1045), uint256(0x21492ed6ae808778bcdf2e8271954f458cbb973c031c580eccff702842deaea5));
        vk.gammaABC[239] = Pairing.G1Point(uint256(0x110d16d9378c3f5248f3151f3f35333cf38300112a0e5b50aade49d99833ea37), uint256(0x2c436a3db5233b311e35f68d10281b4a87c0c993ff8d9d7c6e15ac98f68aaa11));
        vk.gammaABC[240] = Pairing.G1Point(uint256(0x05db59e95e13f33829d49bd7b4af08e83ea25f37cbfb724f78a9cb316cd14ebf), uint256(0x1811b70836741620f1e48879702b65448d7346f0f39bd9c8a14f41e6bfe29636));
        vk.gammaABC[241] = Pairing.G1Point(uint256(0x18ef9e51ded8c649ede3165c4872ee6b57e37e86da2f48d5018b7a485f968546), uint256(0x2bfd240d77b833b15dcda1e622619e84e013d221ec0043961e43dc8602af11a9));
        vk.gammaABC[242] = Pairing.G1Point(uint256(0x2e8fdf68e2f346c5817a7e4ea33e5f6a2be7f6e05b6de62d33b381547eea230d), uint256(0x0443bd0825a850932c48561bd3be51e2f831c3f3be7f53740e9b323ef78d3e43));
        vk.gammaABC[243] = Pairing.G1Point(uint256(0x1f3e67a7ea3be3b63150c48bf2234ddcf52e6046a67240d8c0420756502f0279), uint256(0x142c18c4dacd22e05d32771148325bd696e97c2d57e070661152b7bdcb48996c));
        vk.gammaABC[244] = Pairing.G1Point(uint256(0x1e828746319c24553fa56f466d7a78049c19873dffb4372f07e8d88be9502c0a), uint256(0x2739e48793a87b2118cd8c5dd548f08df8d25cce550628bb22c49495729de10e));
        vk.gammaABC[245] = Pairing.G1Point(uint256(0x05ab75bda5435fd26300558c50a1ca8969cf7d103c42eb30ead5c2a2114dbdbf), uint256(0x0e109e1bbc0831f9d1a676f08391d20a073c5818e8116c21695ddd25340703c3));
        vk.gammaABC[246] = Pairing.G1Point(uint256(0x0ef51b66db7269c4635bbe66ae7ed45727ff0f395e0725819f564ff1c2cdc810), uint256(0x1045c9078e13970131c9373eb5a41bc2e7b837458ee806f9c39dfd505be3471e));
        vk.gammaABC[247] = Pairing.G1Point(uint256(0x0ab9ded9f7fbf555cfa10aeb99dc57698ac3eff601bdcfc2b7c23b056cd3dba2), uint256(0x0fbe9eab8a994a6a338b6b68031ba883d93821e56abd80d84c231d90a3ac90fa));
        vk.gammaABC[248] = Pairing.G1Point(uint256(0x0c70e4bd68c3c4e33a0a4e184b618af81ea320f5719eecc5caf5b41ee0d44a75), uint256(0x027c0ef13697be079a624cd152fef7a4d0c52b3beca5eefb775dcd836251ed2f));
        vk.gammaABC[249] = Pairing.G1Point(uint256(0x1b9d5904ed14fa31f42d65b67f1bda943c0a19e403459397450397e9344bf507), uint256(0x28d34b2545d582384bb66bb3a046184aa770ca774398944c0586a3793a556b47));
        vk.gammaABC[250] = Pairing.G1Point(uint256(0x1799e9872dfa7a96f39742c7895634019f2eb0da24df5af28499f525c703472a), uint256(0x241a7cc1088013341819bce8244c85c76d9279f9cfb61ae6b44287c4653e21fe));
        vk.gammaABC[251] = Pairing.G1Point(uint256(0x0f1b3198a75cd92dbe8a9da50d17f15f639132a2c4f83b1b44ab78978bdcafca), uint256(0x2ce5b59e6d6cdf1be662f291eab5196f1c1d318a13bbb17128f351fb5c12b99f));
        vk.gammaABC[252] = Pairing.G1Point(uint256(0x144c4b58e7aacff260b0ffbcc3c6c02972d402f711131ac21bd393122293b400), uint256(0x022a506f243d974c7c224e3a454f11520b72ccf667c8bfdd28bf1ac10afd0f01));
        vk.gammaABC[253] = Pairing.G1Point(uint256(0x1fb8be44d9c1dd9dcc24a28c1d3ec306920e435bd1b401d564cb40ea0aa8e226), uint256(0x2bf7bd5a72624a7df4477b224a6a015bc2e7f0bde9af573ba5188c0997f44fb9));
        vk.gammaABC[254] = Pairing.G1Point(uint256(0x26fce16835a0a52afedde20d12d0b28ec65aebfe1f15eae333503fef4413df5e), uint256(0x115d573590abf963e3df0f75ead63101a7ab6329afce6f5fcc06a34c9451d860));
        vk.gammaABC[255] = Pairing.G1Point(uint256(0x162370ee6b0558ce0cee287d596f3932fd742bc68be42f0394cde4ead0c4e1fb), uint256(0x287824c6f6d7d9532948113935b3637971541342de5a165f796970bd2114e446));
        vk.gammaABC[256] = Pairing.G1Point(uint256(0x0b8783d1cb3bea8a66271c88bada7f7b58fbbbf175868f1a145e20c3e319bb80), uint256(0x032ac75f48b3cc9648970dab4c5b1e01a0c0f5bb133f078b24c748645184a00a));
        vk.gammaABC[257] = Pairing.G1Point(uint256(0x1a32a0abcac325c1816441f19e98c960888703d2f367b61378f9c28b9aad70f3), uint256(0x167c56c4c341d46d082c7696ce98423354be9057ef1bd047243c579ff8b8577c));
        vk.gammaABC[258] = Pairing.G1Point(uint256(0x2b9f1bac05b640ee1010a8f6c41a5f9a2dfdf6d8a03d282e0d8f21f635718faf), uint256(0x2d42af2a3ff10d6bc3dc343fcd9196cc607d5a4885c4b33c81ceeeb52d1bdc7f));
        vk.gammaABC[259] = Pairing.G1Point(uint256(0x16821a84e0cb7293530e615d0fc51274270e7d48613e399b102ad7da0e7ba0fb), uint256(0x11dd3a155c39a96fade094d463a26eeee6508491fdcbbf4e722fe7f41092f571));
        vk.gammaABC[260] = Pairing.G1Point(uint256(0x0a94bbc5d59da9cf2ad55e4d7a642177ddefba584fc6bc4726bba5ae1e2534f7), uint256(0x027ad1d0ffffc11f2adeaa65d22a018b95406ea826b2183cfa6685214532e517));
        vk.gammaABC[261] = Pairing.G1Point(uint256(0x2f911f698ec30ceff415a578d2a680bbd856433036e985e5e3992c6bd2c64c8b), uint256(0x26e323df1184013bc9176b80f9806f6905b081a1362e8bf68f0d6a06c6772ab0));
        vk.gammaABC[262] = Pairing.G1Point(uint256(0x2c3544ba2a6ee7ca1e8806aa2562e312150d5dbf4febb905a4101cbcebaab831), uint256(0x1f2c29c71adef44b4d3e8c05d746e4bedf0b1142efc94a203e40e2cc3320ab88));
        vk.gammaABC[263] = Pairing.G1Point(uint256(0x2fd2cc2ad9feb2204ac893ae86eab59ed5d9356c45cd64bf63292549d4de2871), uint256(0x01396efc5eec77c0f4e84895cd1c18ec25ae2d73d76a20e91eb86091db9b44ca));
        vk.gammaABC[264] = Pairing.G1Point(uint256(0x01de1a3967582840bffc5e7fecd2a4eae1be0a8adc0df2e7db72c8780688c7fc), uint256(0x2d6ba9e339b00528004b6982b6777f47a569a2b53117a3a8474e36a91553cd52));
        vk.gammaABC[265] = Pairing.G1Point(uint256(0x062f5bdeb79b9445efdff8aea57d417b49a3005b852896b17488e84a3e457bbf), uint256(0x226e84de082514a749cd128f924f4c2cf168b3d26a9c5e29ba97fa949bb1936c));
        vk.gammaABC[266] = Pairing.G1Point(uint256(0x00426ddc1086d7b3f48f08ac48c701663837dc196caae5cc8bc0e723e8f9a747), uint256(0x1b05c5e7b5c0e8ca74f95c71b251c90b915f8c20008f6eb4d0a16a9b1bd26f90));
        vk.gammaABC[267] = Pairing.G1Point(uint256(0x07566a8afa1cecb6dc1d916007656345121fc9eebbe9234635daf34d732f73de), uint256(0x009ed2cdda2317e30b34c02f556e554765044132220ad0145e1b60f91282c9ca));
        vk.gammaABC[268] = Pairing.G1Point(uint256(0x12b21221f9a1804678346018d735609fd82ee0e844a92e5a8798f603b860f4d6), uint256(0x0362e47357400045baf39cd700ad7d8dc35010973f6d0dd0d87ceeb8acc8a385));
        vk.gammaABC[269] = Pairing.G1Point(uint256(0x09ceea973cd4b344fb5d148d35f5b3af90fcc8eec66916a598c74532d7f3ac50), uint256(0x1544f4f78b7e8dbf9e4a65a69fdbd9c3f94cace0de099ca7cce37bfd16e71e31));
        vk.gammaABC[270] = Pairing.G1Point(uint256(0x1623b7a80ccd746e79d9a8108c6bdf1dbc1dec4c1450b90d167da793d8cf8788), uint256(0x1abf886961e12a455ad12e8ffe4a8b4a24a0e05d4fe9e56f913b426310a8be2f));
        vk.gammaABC[271] = Pairing.G1Point(uint256(0x215fdd0be7547938bf059d1cde0cf9defeb85c881f9b500feadbfeaba891ef0a), uint256(0x024492d3ab5bbb9011a5c2629b25a8cc7380d8679fe430bd59699723132e2761));
        vk.gammaABC[272] = Pairing.G1Point(uint256(0x1f7384391e0181d0dc574e55096a59634a11f629e2646612c263c6c58436e0f8), uint256(0x01e6089753c37a5da33cdd62ab19d0a4f4f68cf5cff2a69629ebac86b86735c6));
        vk.gammaABC[273] = Pairing.G1Point(uint256(0x1987515c4e722a71d0c8f9510670c106fd996f60e985dc17c95cf0e92cf57333), uint256(0x2c835eeacfb358d1206f15d6e16d9377c5e2484cfb9bb1dadc9e070218732901));
        vk.gammaABC[274] = Pairing.G1Point(uint256(0x2b8695727ab57346d831660f768a03b9a8e1593524d4aeabf3c5187a97c7ea00), uint256(0x1b9613f709088aee8e9a474408a146217a0898934ecbff7788c81ce095e55482));
        vk.gammaABC[275] = Pairing.G1Point(uint256(0x25f30200066e058e11d5ce76c864d5ddfe4cb44eba17ba3a3231ebafc597a238), uint256(0x07a37e65b4b053bf41d555e7466d40162c1a26fc2a58142c5d85bb024d85620d));
        vk.gammaABC[276] = Pairing.G1Point(uint256(0x0250d557baae22dfea5a87d67d71adae184d3011a163dc293441f7e8cb227a8e), uint256(0x09ac5e8c3d9275e9acdf2dd5de3cf2aa4c216e542cf8837b893650f0c55cd4ea));
        vk.gammaABC[277] = Pairing.G1Point(uint256(0x0dc3ab19c3f5c0852ecbb0575e210743445af4dd4ebbfca18b60b5a34c558013), uint256(0x0dc5b3fee0b9f754322e0bcbcfd3cbad329b4c1063d37f77acabe2996bc85e70));
        vk.gammaABC[278] = Pairing.G1Point(uint256(0x115ed3fe1c597f9a9001d5a071b2ef63058340f06bac8544d08fc4ab3cfc92c5), uint256(0x1258c9f81c5a7442393fde1a6428385355f20af81fdac4516ce45ee2f984e46f));
        vk.gammaABC[279] = Pairing.G1Point(uint256(0x0029b4357b76c4fe7f6a87e636b38b0da02f3d1aa9ee3f8b19eab8c92b0620d0), uint256(0x2f9529793d30f2a47b0751adce28378a82dd7342ee0f4452a177f1633496680e));
        vk.gammaABC[280] = Pairing.G1Point(uint256(0x0a5bede584d32c14a5cb4d6e53511515ce5412be6468463029b8bed70aff4e25), uint256(0x1e4993ff4300f156fb1cfedd4e0a6bf9a208859c4642ff4e0e2345bae15ebf95));
        vk.gammaABC[281] = Pairing.G1Point(uint256(0x2c1face850e2b6c8d56f5c3848b6c38eaf6ef06edf239cfc85f2b2a0b41a186b), uint256(0x2e1ef71fb2bba4ed48c5a95db99e0388891fa6a3cb66b6722924609f2edc67e0));
        vk.gammaABC[282] = Pairing.G1Point(uint256(0x133fcbdde88715588d7f07bc9cc8f8e4c0829a7e416bd0b665cd4f293cd1af9e), uint256(0x1cf79d06669af45e06a99f63360825d95cb8fc8d92c527c65c55da2b1f3232ca));
        vk.gammaABC[283] = Pairing.G1Point(uint256(0x1c5f41a43520f86530952a4060b97df7b7e3d9557c32a27d9c926a7038728c4d), uint256(0x3014647303859dd50bb50a36b03be0a821603a6f795ee02d8b4a255ca9d7e4dd));
        vk.gammaABC[284] = Pairing.G1Point(uint256(0x2e841a96412ca4c7313603b34d0f023cd8bdaf61da7842da23fd58c73ed73068), uint256(0x29480b34bc81a02d7d4e9428f7d8fbaa1d7710af4fa01525038ac9a2eae88a8a));
        vk.gammaABC[285] = Pairing.G1Point(uint256(0x1fde5dc99954642fd9697d0c02aa5f4fcfce3c767a052f2e504ecdd71b977ece), uint256(0x2d8c454a1bed1aed68118f0742f7018324a6f25ecbd4111b0837133aacda6920));
        vk.gammaABC[286] = Pairing.G1Point(uint256(0x26b6fcedff7fa9c5290a5e4e7a81ea76234f23840a821541efaf3013d17cf70d), uint256(0x28a02740977577f10631ca26647477f75a6245c034cf1b3e4ca238f1820b65b6));
        vk.gammaABC[287] = Pairing.G1Point(uint256(0x1c51b0f714d774031e6f01d5e9830bdad74b72f5338e1016313c40e340c4fca5), uint256(0x0c9aac788071f1bdd23c350bd1bd0542c5c480ff6ac51ee21da7a08a54f0dbcd));
        vk.gammaABC[288] = Pairing.G1Point(uint256(0x04b37843b0d28ca6a47da17fbb2452b1f65342e22653e6580925554494c8913c), uint256(0x043ba726b39d84ed288b29029ed1987d7f6a09facfb024052018c8a53d0eea81));
        vk.gammaABC[289] = Pairing.G1Point(uint256(0x154174f02fac1c4fd8a056ca4ac94bc36e720834917a347f1e0c5ec0a2f406f9), uint256(0x2a96adf55a94c389f32608dcb24e076d55fb57ddab185180033949e9c41c6976));
        vk.gammaABC[290] = Pairing.G1Point(uint256(0x1535f60c9b0c1d815f7ce9280b739570907fef28f1354d6328e478580584d93b), uint256(0x168c7a129d53b6d418a6d0f5d998dff350eb362a26ed7aa89082624b072d42d7));
        vk.gammaABC[291] = Pairing.G1Point(uint256(0x247ba6803489398485f0bd55ff44b6b136d54e355ecb21a1afadc091ee0faafe), uint256(0x236fad03e3a2c78c5f2c2a1d8c9b85d3d771eed251759380734a933674c1c173));
        vk.gammaABC[292] = Pairing.G1Point(uint256(0x2b5e81da6d77ec53751d39f4befc5d5074a4303e0359b28012eb99f915ee3c59), uint256(0x0d8b8ff0701fbc6d1d4093a65101fcb576a4fba146941977f4c4c3849a59bd3e));
        vk.gammaABC[293] = Pairing.G1Point(uint256(0x19cd09ced87283a89659359c6ea5510ca6ef07abc51ab0fef08b9c91d8a55c87), uint256(0x10c894ca70ff5ea9a43fa1933c4780c5ab22e822934e45c6ca2660c65bcd38a9));
        vk.gammaABC[294] = Pairing.G1Point(uint256(0x1169d899219a29c312d5521dca8a281d2ce5b76dba14f58f0610673772cc385e), uint256(0x03bc3b2bd779bf5f1a33e702843f5f99289ca21300b277072f2a03a5f61bb9ec));
        vk.gammaABC[295] = Pairing.G1Point(uint256(0x0245a9df559ae0994f19f5883868d7002faf37feb4781eb487b74d71b8211a9e), uint256(0x2d39f06a6ee96d3c3766508e342eb1b8b33e9eb143f81174a0779d4e8aa213e0));
        vk.gammaABC[296] = Pairing.G1Point(uint256(0x2c2a75a6c82c11b0ef25a50582e591aa54c8db05b2053e97643c449114a0ec22), uint256(0x1952f4032876ca5bb43f87f62e8b1ca436f185f3864c56e71ec50c009fb14c91));
        vk.gammaABC[297] = Pairing.G1Point(uint256(0x266088a90680bc6177794843d2618c44cd93f3e107718d1c67ecc915df3389ae), uint256(0x2ddaefb734383b2e7dd9b9f02c9ca17b5b1e4f9379248a7e37c0cb2a7b7c659a));
        vk.gammaABC[298] = Pairing.G1Point(uint256(0x0fa1f80e77de57d82f50ddbb4b6a11e49702ceb1e785ed96e7e23e349622bb46), uint256(0x10bc8199e01aaf0f8f6801654c7408ee2ec9397645e5322d6ef1fd3962ba4693));
        vk.gammaABC[299] = Pairing.G1Point(uint256(0x1846e0319649edbe28f6684ceabbede083fc46a80afd4990ce582baf231e1c22), uint256(0x18aafd25a7197f1be6811c32eaf786ee1b4ddf858699ea98b60fbdfe426bd789));
        vk.gammaABC[300] = Pairing.G1Point(uint256(0x10879975c6d353cf4031fa2649ab272f88c1e5c92f9d8eaa35c95038e6be7fa4), uint256(0x008956fc788c0774c60a02eea45243eea3b381e8c90b403f5fe3645e38a65182));
        vk.gammaABC[301] = Pairing.G1Point(uint256(0x0fffc4c902ae504dc19f28dafb3ae1ec362ff2d820475354039510259fbea952), uint256(0x2c7247cb7c74aefec402f8b460eb1008e246508d35b81d919b0a629f88c2e31d));
        vk.gammaABC[302] = Pairing.G1Point(uint256(0x2781ef9e0a79646e5babe78e9da100c90be0e2cd2c2c7cc57ff7f4f1daf2dfe8), uint256(0x03e534e65e198e448576ac5331a0cc05b39626397b722aaa2b7adf123ffa5401));
        vk.gammaABC[303] = Pairing.G1Point(uint256(0x0a03f5659e5e907d2969071f5c4176ed68f7a7eb26bc1cffb4ac9d54f6bd9ec1), uint256(0x0c21e5fb3ee4deaaeb24d8ddb3b25f12e32b4f1d219d0d540de25dc3fe0c4736));
        vk.gammaABC[304] = Pairing.G1Point(uint256(0x23810bfe4b8b22bdc4c579a0102425b15b045c17df9d0401df22f370ef70d7c0), uint256(0x05b69236d3c3f5cc68d9d3f8394b6a112bb66cb84f8b5eacefcddfd6ea4d6e83));
        vk.gammaABC[305] = Pairing.G1Point(uint256(0x25395e8ad8c8bb5b07a0d7e2afba0fa47150e90b5f5915d9a4abfa7b302a2199), uint256(0x27ab1d18b412e0260af1dad70620743f1c8240fadf7159419eaad9dfa965f1bf));
        vk.gammaABC[306] = Pairing.G1Point(uint256(0x1b8f14dc9a9ac740753483002c54edce7c9ff5943d240bd94f9b132fd4f0d4af), uint256(0x0ed37d3d677dacae61126e128c92dfbe6b514c06c0503104915a2811efdc60a7));
        vk.gammaABC[307] = Pairing.G1Point(uint256(0x0b063e1a43ef91c15a9b48ecb9b82e16cf78ebc75a3b297ac28a0bd996983655), uint256(0x2c28d341227638c14fed1a1eb4419dff27a5f4b0c7b382f88fa6e938eab74317));
        vk.gammaABC[308] = Pairing.G1Point(uint256(0x02a77cdc010e7b82e86ae55f118f7ab3c9d96b5330a1b4483ff8980ef8aa8a8d), uint256(0x2110c8872221777b73c89597071a10d50650be64449caec63689493ccd1b2645));
        vk.gammaABC[309] = Pairing.G1Point(uint256(0x1234b5ba2bcacdbfc8c24a024359a77832cd1c361eedf0d755bbf9a855315d12), uint256(0x28af06b5f9b0d08cd0d696cd668dafbeaa907e406cc1b20f0793f89f7106b827));
        vk.gammaABC[310] = Pairing.G1Point(uint256(0x2ad7d3d3da2489cbad37c2d2fd6862fdc54d00116367da8296a49c949d853eb0), uint256(0x16c828dceeb0fe26aea1375e363193a39704a1be59fcd16ac4d7ab04cd013dcd));
        vk.gammaABC[311] = Pairing.G1Point(uint256(0x17a1c79346e4c5eb94d7811fb9a2231d5f22f87ef4adad4cf38cc4d3c60cc93f), uint256(0x25ab611b2dbb160cfc4cb15d07d4fae2235389c2b075cc2bd95816b9e1993ad7));
        vk.gammaABC[312] = Pairing.G1Point(uint256(0x2823a065c22395b7542c5932c15be0ac749b39af7c073e99be13b1635f583c2d), uint256(0x1ce9deb1f561c9e577317992f35395e13a1b2e156cc371bdba7d4c3c2a5613e1));
        vk.gammaABC[313] = Pairing.G1Point(uint256(0x11a3e4fcb6a9ae0de5f829c4f6f70950399dc72796792e067f6fa88a3e44aed5), uint256(0x2d86c8ab77042b5576095732f5676d8eac878442cb23efe44e0ef6fdce919e9e));
        vk.gammaABC[314] = Pairing.G1Point(uint256(0x11fd5d8ffb4e9597d0f547ec26304a09d91648a1e3875103130ea72a46598d4f), uint256(0x11e761d6972dccc9f8c13d328860ba9fa54f41b08adadb2fd9a00d09a6cdfcae));
        vk.gammaABC[315] = Pairing.G1Point(uint256(0x2c980edda7d56744c4a5e40c5a62144ce364ca055eb06f5a892aeb5933d89fce), uint256(0x1544bf44ca3b60efae17d7fb7d3246e579d57d5580a3e38e3d4bc1d59fab3b24));
        vk.gammaABC[316] = Pairing.G1Point(uint256(0x018d4ea42ef2008937229439974bf20bab1e230cdb615a3ad5469fa591263a8a), uint256(0x149f1da4a8ac1309d454e0fdc1a42a8200a8ed79110626804722ef1ef38aa93e));
        vk.gammaABC[317] = Pairing.G1Point(uint256(0x2bc7c3767d5bb8a66b7a4c8bd44cc5af14386989df3ee3870001eb9c5516e044), uint256(0x22f42c091ae020e04281681a26ac4ffc3adc2f1647be82c94e988896fedded3c));
        vk.gammaABC[318] = Pairing.G1Point(uint256(0x080da3a09cc545723fd165127a7a82bb57bfee031743c0e0f68be1bf8278cd45), uint256(0x2e62d03047bdf84f47b87abf116c5a321998a4f58f0f05dced44532687253786));
        vk.gammaABC[319] = Pairing.G1Point(uint256(0x097d639465f488d32d1e0676562f52cab465dda0a8ca4fd47a71411abf14232f), uint256(0x00dac8da6d82127975ac56ab8182352fdbe634323ff05bba376a3b580025b49e));
        vk.gammaABC[320] = Pairing.G1Point(uint256(0x1f1374a247c9d6da9aa802aa1a96a6b96450bd17811c3744b4a241fe7c266486), uint256(0x0d7d2a94eda721d2d068a0864295b8ee7afa4d2e0227b35a2c82259a16b5a19b));
        vk.gammaABC[321] = Pairing.G1Point(uint256(0x1e60b2a0ae2908a5a6abeb0f312b0fc9ca26b1d3c301edc77ed13452c4684688), uint256(0x0cc77a464721551b068007ece3d59d91fbab01c8374ffcfece50aca4e1159dc3));
        vk.gammaABC[322] = Pairing.G1Point(uint256(0x1ac8a2238afaaa1f920ee6a6582c52fa09c62fc644a1b3951bd8ed85350ecf00), uint256(0x183ca82b57d872719a20aba7254121cdaaf5adba8d28f2a8f1b8954a2cc7e0fc));
        vk.gammaABC[323] = Pairing.G1Point(uint256(0x2860a544a9acb80c9057658a1c183a0f8c91b117cf1d6a8f459f96fd4296f317), uint256(0x24c9c02d8cd8eeaf18cdd910723264b7dec1e9cf4aae1f8f1a04fc7fe7f60411));
        vk.gammaABC[324] = Pairing.G1Point(uint256(0x16a498455dfc55a7f722483d268b64bd0a119635c46f1e73abb2b452cbb450cf), uint256(0x074336f03412c8b52a77abbdac92096ab849f0ecfcb911296c0750759bd77b36));
        vk.gammaABC[325] = Pairing.G1Point(uint256(0x06e4fdb270b033edf396a5dfb5c24684cfe63a80267ac786abb0af3c7051d18c), uint256(0x19b9ac2a66a56b64309d8cf830938f118df51bb3be0b910caa7aadc3fa50ee89));
        vk.gammaABC[326] = Pairing.G1Point(uint256(0x2c94b0e2072c23f9f819eca0b80b4fab46b68cc1351950c2826c4cd292a9f42c), uint256(0x038bf87b3d3f0b5ddf839512b99d1121e805b6ed7780e90383a5f254987ce0b3));
        vk.gammaABC[327] = Pairing.G1Point(uint256(0x1b8cacfb09c0acb01e2f6640741ad125a10a25aee0fd68fc09014f4e604b4907), uint256(0x06e885641f84dd50d75ce7b4c54346ee55cc547c624815e76fb5c5c55f50105f));
        vk.gammaABC[328] = Pairing.G1Point(uint256(0x2f5d96e43a921ae0b918852b66e3355f7c7b419d5c4213acce8a72e5083cb61a), uint256(0x1e94f56482712962b87d8d76c40b69038bc3fcd6d22b2a340cd236781e737be1));
        vk.gammaABC[329] = Pairing.G1Point(uint256(0x2650459c708129729006f20b9d3b572b4694a77456479e748af1a778247240cd), uint256(0x27c95d6756d7d8d88eacec086eabbe1560a07e4ad24f0020268f401374784b9a));
        vk.gammaABC[330] = Pairing.G1Point(uint256(0x0d45853a6fd104a2df82c5832bfe45584c29c4f9f37a76737e53551a0885a935), uint256(0x25ee8b178ca9a874923b19dcbe3a6a55e9043c2bbdb98330c2436e810138b526));
        vk.gammaABC[331] = Pairing.G1Point(uint256(0x1f4b9266e3a46b56b7c4089faee15ef291bc943612aaff917f314a6f0deb9a2a), uint256(0x0e420b1fb7232d2e09d670107a0c6815ce74188b6268d4c5dd7b2f7f8809cf56));
        vk.gammaABC[332] = Pairing.G1Point(uint256(0x1a270629becfab8bd6913299f26f14078406b14c386adbd104e33eaa5597687c), uint256(0x2f0d2ca735efa636790ca31981a8b3085f16923dd748467a0c1f93591432962a));
        vk.gammaABC[333] = Pairing.G1Point(uint256(0x0c80ee0ea0e15ea3e927e7d60490d70bf914f07b701985feff850bf5f0e6c5a1), uint256(0x0b41ba9a391493f556e4cd896dcc4720e342a54522fdcdb7def4fa9cb045ba74));
        vk.gammaABC[334] = Pairing.G1Point(uint256(0x1b748f416b86b101d6de2f86ede5d9341be66b929be4cfd4297f337154cd4349), uint256(0x10891ea2e29be2c2e5bd390fe04463750a5ef06f18561bae6dd8d28d8f72ef82));
        vk.gammaABC[335] = Pairing.G1Point(uint256(0x1539a2a446d1f7b51c507bea84798fb2c3b91d8e8339e25171b5c07813c67408), uint256(0x2c5f903c165645d842d31fb9e02c0f04827b6ae72ea3d0f4c7807916381d05f9));
        vk.gammaABC[336] = Pairing.G1Point(uint256(0x01eb5864cef7bd334cb447f844f9f81e7e4aa100d7c928d357583381977995b5), uint256(0x14293bd387e178f8b4df540842c7fd713352e0830cbd85c01c3940d307cded00));
        vk.gammaABC[337] = Pairing.G1Point(uint256(0x0dac62f3c245ca1bf286042bf24939c001f8a5fa091c253df127ae7d895b5972), uint256(0x144b1bb28b4bd23657c354b2a5f40b026f347b5c45597e88224a810dcb0bead4));
        vk.gammaABC[338] = Pairing.G1Point(uint256(0x2d5a1f0ada99053b5b053c7c25fd2f0b3b3c6e78c7d46639a75f0ea0c5542997), uint256(0x27929e3d0d26a84be9ea01ea6e51903257e89fc4f698997f8100a2557de37080));
        vk.gammaABC[339] = Pairing.G1Point(uint256(0x28df7ee16cce36a384bc21dfe89234a7731f22d5c9c568d5077d8339ee229a9e), uint256(0x23368c2f7ce24416a19dacae610bcbf5b35fba947b4a99279e2f078b440690de));
        vk.gammaABC[340] = Pairing.G1Point(uint256(0x0c26bcfebd5b895115e0a969d9d4bd64aeedb89222764a3efd8a1bdb4e763c6c), uint256(0x04fdf556b9f73fa4643848feb5dd7db0e2574f7605735ac7e7db3188a120c602));
        vk.gammaABC[341] = Pairing.G1Point(uint256(0x0f3d72a653a8daacc5ea5ecd519cbab3389898419e08b5e5e38f86e9601d8e5b), uint256(0x0d5b3ecd8d1c8e31fe7782a8a0c4daf0e1d189b5d4c6cf9dc5f3e86515a6cc48));
        vk.gammaABC[342] = Pairing.G1Point(uint256(0x26cbc0c7547f766154e6e4625c397f2b9668f19400c06c2bf2ce0e4672d10068), uint256(0x24445c8945b3e656b1b26ace52dee8e2b812b8645008ba3dc884a1dc8e23387d));
        vk.gammaABC[343] = Pairing.G1Point(uint256(0x255c9e956105857af94649c9166af7c82f9597f4a5ad6d5921f6ca105d571d2e), uint256(0x145c1f463d44e6057d8e219a81cc0cb9ce3d289aab999bb032527cc21d7d76b2));
        vk.gammaABC[344] = Pairing.G1Point(uint256(0x09a63becf15538bd7695188abedb7b8b7c83b8288fb427df0fcfa982aacf7b1a), uint256(0x1e724d10fdf3c4f4442013d1c5d7e06e294bdd1e9d80340b9cff7c997b2395b9));
        vk.gammaABC[345] = Pairing.G1Point(uint256(0x0c4089e0d0467b5d31285063a923a56463af2a0a6b2fa5a21e458725d3c2305a), uint256(0x23f54d55257c66fdb21b37ee6f0b1c3f893144ddc83e5ff640ec26e99e9d1926));
        vk.gammaABC[346] = Pairing.G1Point(uint256(0x29f0b9f792c38f06d3647a1b7f58e25539324961016e36eb501d5925dc6e4148), uint256(0x16538fee71907a75d5ef6e9213aa07e6392cb075fb4f2c0f55519afbcc04e13e));
        vk.gammaABC[347] = Pairing.G1Point(uint256(0x214a3cadfaef0337072eb9f56abc4f738dfcba24cde8f28319387304228bd3a2), uint256(0x187d21e619092787bc6e63abf0f2c7026f700d033280397d535034c93845c436));
        vk.gammaABC[348] = Pairing.G1Point(uint256(0x0468fc664ba7a393dbc5adb2b8eab3f1647e7f13d7a90180677a5455a0d6f59b), uint256(0x1e9735131d83e96ed42ba130cf4ff2a220931c771010bf26e57416221b40a312));
        vk.gammaABC[349] = Pairing.G1Point(uint256(0x04c0041a6e40c41d5ffa9eaa25d08f1fcd043ad9aa60a436a38c3b62efd0212c), uint256(0x245ca33622b78beeaacb06063b85390b4b24446f3f486b2a2c025d5f6aad6ad2));
        vk.gammaABC[350] = Pairing.G1Point(uint256(0x15be4b3ac8ee7eea8ef6cee293515ac7f1f6e28db2b9f299f1391c978eca8f69), uint256(0x18ba3f6a92fc12c33c0b66cd6fb651976c0ab1b2d985644873147b227fafe50a));
        vk.gammaABC[351] = Pairing.G1Point(uint256(0x2a768e738e3303d879acb95e08965bb55d63078e32d6082e60f8ad74df083675), uint256(0x00332c3fecbb0f7a27a33efa9ccddad1fc04f909f36cffcc503039fd6b0d7508));
        vk.gammaABC[352] = Pairing.G1Point(uint256(0x24e60c031da496ed4fd77888edd691e568653d156ec1e6de52c0b1088e3b185a), uint256(0x026f64683675c1d554eb25b489077cf91b0bd8696788dfbf31446ea2772f5238));
        vk.gammaABC[353] = Pairing.G1Point(uint256(0x1f7a992d07da8628a5bf1fdc6c109a1d3f6cacaa3d6d10dc8ba822dd8e7e19a1), uint256(0x07cef944062ea7473471bb8f1e48cd613e1e83e8f94177d44c1cb26d062c3d3b));
        vk.gammaABC[354] = Pairing.G1Point(uint256(0x251675c25448952789f2354e701973de2680b4fdde4a5559f8c2245f50c75ed6), uint256(0x2ae44e86e74abb78d3e85ac5baca9c1a1f14f6c543d94434b98834ecf5c89155));
        vk.gammaABC[355] = Pairing.G1Point(uint256(0x01749e189fa7e10a1d8684ac5a7029122768a70912c0459b2a9afd153e799f5d), uint256(0x2bd070b5c12d82e36fd9647a27db84b5c05660e6e8bb1db40650544adfe092c5));
        vk.gammaABC[356] = Pairing.G1Point(uint256(0x2f3fef38a76fa952bceade6228b539ad75931ccd000d0bf274b8515828192243), uint256(0x0f6bbfc827c85b57a50310594f9bf3902aa647b14d8c40e05089170f1e8198a9));
        vk.gammaABC[357] = Pairing.G1Point(uint256(0x0a3b5e0e5b6f85cf91bf98cffad10e526da81ae18dfb4fad0590f4e3c2da1144), uint256(0x19fd4abe86359d3e1a11a7385ad8a8f6022fbd9db0a031db427ca6c69eace5c2));
        vk.gammaABC[358] = Pairing.G1Point(uint256(0x0cf8127aa76fb7d54d63be831580093ac715665bd991f11e5a9670c5902c27c6), uint256(0x20cd3c2239ec0dc5b7b4ebab59676b26b57cd681b983134d286f0be0f79769be));
        vk.gammaABC[359] = Pairing.G1Point(uint256(0x036493392450192d7547b1441f59df8901ec9bfd702c9a655c1b3a5dcb4a46e0), uint256(0x2c908551a29b36e7b5e1ebfebe4d687917b4b4d25ac0d877233a33fff7b60911));
        vk.gammaABC[360] = Pairing.G1Point(uint256(0x0888377cf3a966d80fa3a73c5f8cfc0cf763c3b5d1804b15c426325b4214337b), uint256(0x2109ecf12ab485f966108f0caed1c1d4ca4d078c5122bf342bf99e1877765a79));
        vk.gammaABC[361] = Pairing.G1Point(uint256(0x1a25ed36f1733a3406bbfedadc55afe59fa988ac093717917c77ad4074d7434d), uint256(0x0d3eb484bd1ce9921cd842cff7942e2cbf15f5e83ff2592002488698f16dcd39));
        vk.gammaABC[362] = Pairing.G1Point(uint256(0x1c502ded320a9cf684ec6cf921857cbd24b03eb37d0eed1e7985d0249f5c4f91), uint256(0x0cd9413c2c8dde59ca4243d892c864b5e31ecca0a45f279368f933b79a69b698));
        vk.gammaABC[363] = Pairing.G1Point(uint256(0x1c43c5a610d535cc38dfb0ee2db1f532a369dd199b8ff60c81444230480ed909), uint256(0x26f5e242b92dece2ce8f8841e68c2cf53d23c93f45738aabdf22f5055476d08a));
        vk.gammaABC[364] = Pairing.G1Point(uint256(0x24c8316f69bda6c9474ba6146d4e9d75ddc513c8f0b8efd36222b3c612d58f16), uint256(0x19ce9555c7df3020bb1d1695be0c82c5280b7ce1fb18a1d2faaff78b26908e48));
        vk.gammaABC[365] = Pairing.G1Point(uint256(0x0d7f0cfd7ddab65fbb806ab2a1cc5690e36430eda77695882d3d69de2bdcde7f), uint256(0x1b63f5b9e8a8ffc71764bd36594ea0e8e6ac1c76eb7d4ff8ac9d1a06a40d2638));
        vk.gammaABC[366] = Pairing.G1Point(uint256(0x2139d736059f281be3709a7eb6373ff7858f6509acd605d361e981d9b5277d59), uint256(0x000677478be634319a68704f7d78cec608c8b5fca37d6050f62107b072003446));
        vk.gammaABC[367] = Pairing.G1Point(uint256(0x0af386e9381978c5dbcb506f75c89664033ed2ef69ea44a4c13637f14fd40e08), uint256(0x022a655039e4fe568ede32bf54957d668b3a59b830613b7299ad2db24fb69959));
        vk.gammaABC[368] = Pairing.G1Point(uint256(0x03507b9a848f88f5a09c56796d46547c6f5cca87674b0928e8a19ae624aba1e5), uint256(0x13179a631ab2a1ea74fd49349a7c560c5919ed61c095fb5c48a9cb25a596c3ff));
        vk.gammaABC[369] = Pairing.G1Point(uint256(0x1e0dc00b9f48f811307c15d1c04ca37ca28d9d1372d20b3c4cb5f734fcaffd0c), uint256(0x2d153b0e969ee84464bb1349ea6e16038584ee8c97639c7a946e247469a3c4b0));
        vk.gammaABC[370] = Pairing.G1Point(uint256(0x2133ae559aff89764462cde3043633b5af7b787c7ffda67aabdbd8c57740d315), uint256(0x0be55290955b92ea4a2164b96ab535a8522ba847903a6734060f9c10e881de6d));
        vk.gammaABC[371] = Pairing.G1Point(uint256(0x17d0f8acabf320251d19c0f5a35f6d205f25596a1b2c1d20297a0faec6995a86), uint256(0x1101b5c7ac3d193b900ca3c2c99443d73429e4a2d024ef8b6f62d268e6f8c52f));
        vk.gammaABC[372] = Pairing.G1Point(uint256(0x1aaf706a252a9cf4dad2aedba36a99cc66e7c83fee0de602025d889fc9ce3c9b), uint256(0x1f88927aad62078e399b7c410fd97fee0fd073732c238508b9e23d3bc9b6cde4));
        vk.gammaABC[373] = Pairing.G1Point(uint256(0x2fda807a226a99d67743d9e04891aedb7dbc87936a387e55cc8040930764500f), uint256(0x0aec2f639837c709d947530f3789a593d74086e361e665d82be984c09708ea86));
        vk.gammaABC[374] = Pairing.G1Point(uint256(0x11b5555bbf2a3c8bd8ba0f062cd7805047150fe80c0d6cb33f025e4ac7c3ca07), uint256(0x2fd1bb48098701c228b02ce877ef2fd0c474a168352eb8d1305d91f6b29ad040));
        vk.gammaABC[375] = Pairing.G1Point(uint256(0x12f9859db6fbb112f917f8a3f9e9d166ae3628d7d6f7e4a6394aa2616ce34b87), uint256(0x17f81e875534a6e23c1146f830af3472d88bb1aa039e879f4ec2b1ba5148333c));
        vk.gammaABC[376] = Pairing.G1Point(uint256(0x25ee849d8a2f6f8d5269227c5e165724e27b19942c200ee72ce3d8b7f8675896), uint256(0x02f94434f9d6c8d9e324a841bf7a109193cb2c2df8039c37596857396d04c7c2));
        vk.gammaABC[377] = Pairing.G1Point(uint256(0x2419acbe0b8b13a98b6e872f9c117712a454288afe7072f64874088b63e96c26), uint256(0x027d5919cc8854600cc8a66ca27ca3bbb78e5520f46facb4ac051854d84cb591));
        vk.gammaABC[378] = Pairing.G1Point(uint256(0x22e44c594be8532a0071b4c5f592161a11d34a3f0a67322e278a23b36c0a0720), uint256(0x1bcae7b65956297d09cd6de2cb9bcd03aa177de47644f6b5ed2475b027a109e9));
        vk.gammaABC[379] = Pairing.G1Point(uint256(0x2187c3d15aaaae5ddc653b46a7aaac09c6939170f87a1f890327b4ab351a1027), uint256(0x00d6b4171c4b2eef17f0cfaa41d909c74f9b1831b08c82f01f99f65d4073fc56));
        vk.gammaABC[380] = Pairing.G1Point(uint256(0x059b191be0860720a08394fb159552c46c612aebd7163d2631058d696f1df8e9), uint256(0x23b6970b48dc0c8bbf56d06a94b8957af2c90f105cff104106acf6bc66dd170a));
        vk.gammaABC[381] = Pairing.G1Point(uint256(0x080c36a1689ab1901ff0c2147bac644eb928b041bdf7d64adb643cd9b275c644), uint256(0x2f66a19d5eda3812531e5e74a5cebc8ea45846027f33a12a13d40d342dbda66f));
        vk.gammaABC[382] = Pairing.G1Point(uint256(0x158aa70b75af3ce6ce0db990f708d5df8d3d6ea83b238d0f3dfc502d217dd4c2), uint256(0x203a725aeea419139b580f33ccf00e6a938e5736270f6d958dd86dfba5db543d));
        vk.gammaABC[383] = Pairing.G1Point(uint256(0x1b43dfb955064397dbb3f3a25e82df38f03daf716bff9d4662fac07e90e417e0), uint256(0x25d7b59033e536f8f4f1fa87c40c7dd100364250f64b4c125a53b8fdbdc6c3cb));
        vk.gammaABC[384] = Pairing.G1Point(uint256(0x2f71dfccf0f5f1141c1392b0a9b9d9d9da13e22ebf3c3587b13b48b00c7491d4), uint256(0x0241f277cafaa3189d4cf3ae32152aae62828e4d7888582725ff7be85cd19f9c));
        vk.gammaABC[385] = Pairing.G1Point(uint256(0x15aa8a6aa587fa8ab9e5bdb5d8ff16460e30594942e0fefd89fcd5710989b65a), uint256(0x08a627aa57ce06b3253f995b1d2b7a4a9059a1053c08f469ea9b23f2466ae4da));
        vk.gammaABC[386] = Pairing.G1Point(uint256(0x1d597a08ac999fe5a855e864515cf94f8763d4b289e85c9da418dca86bcdded0), uint256(0x1be15ef1483ec8cdd0ead0b47b2a94c09b71a1f18f71bf4285a4818f4dde07a9));
        vk.gammaABC[387] = Pairing.G1Point(uint256(0x2835569827693245644c3c4c3a467a6decd367dbc446d9016978aa2d631e17cc), uint256(0x031ff699f7ef6bb77d18e3f5b8b07cd1d1ea5d9edd17332fef7f500fe270c15a));
        vk.gammaABC[388] = Pairing.G1Point(uint256(0x0036f9e9d6a5c5eabba83864b4d3c84c60c53476d432071db5192a5d5e4df4a3), uint256(0x18f5d918329d8a379c32509c534e426b925c26c0650bf9c329161dca3ebf7c6d));
        vk.gammaABC[389] = Pairing.G1Point(uint256(0x28da3005d55c97da6cd7b4d437e6cbaa7f70b96e43bca03f3fd96150f0ce4ce4), uint256(0x080c30dee44fa49eea3c7d1bf5206befab03a41359b88cf196423408ebdeea63));
        vk.gammaABC[390] = Pairing.G1Point(uint256(0x136694c30d871ae49de0f140c97c8da9535bd69bbe1e816d6f38d0825311f560), uint256(0x2a4b57bc54f771e2aeed688a3407530503cb54280c2e180f87a36ca4339c181b));
        vk.gammaABC[391] = Pairing.G1Point(uint256(0x28e5b93e5fe57cafd0946d0419f7d376bee1f7b22025df1eede0b89409e0aec0), uint256(0x0f53e61be5a9813937d528d0d92873081fe80c075a078d77550b587b6937d786));
        vk.gammaABC[392] = Pairing.G1Point(uint256(0x1897277dbdf4599c4c021ca1ee8a9a6a14b00081b727835ab89614d6f046fcc8), uint256(0x1f8d8405238a0390230850a3739d8940b094cac6822ec98a24ce45e7c7292317));
        vk.gammaABC[393] = Pairing.G1Point(uint256(0x0fa9e372742ca5c1981e37b3e28a39196d85539463bba04a92e8859e879e7095), uint256(0x16c484f19fdf17e9574eb2f18e820be9c94773c6340abf057162bb4b7cea0225));
        vk.gammaABC[394] = Pairing.G1Point(uint256(0x0fc06cf5fd86e266cb01740b3d5c973ec15ef6128f1931cf5b176f1629a2f17a), uint256(0x1ca528e569d67bf45b901c979f642defd2335fa16f8a9f168ec8a49e8bb12bf5));
        vk.gammaABC[395] = Pairing.G1Point(uint256(0x1037686d260618b2daf0376d66762b6d36982050ca13625fa4dba3f9e57afcd8), uint256(0x2720c328c49166bffd598ee8ad8b346d7e36b74753c9a255d3b57c17aa2b053f));
        vk.gammaABC[396] = Pairing.G1Point(uint256(0x1011a24095ba66af304efe638f5441423dcedad6c16fde18442abac0c9179d12), uint256(0x1f9fbdfdfb0d3e25b61f0b8635d88ceeef511ffc6cdc0b28b3fa5da1ed904718));
        vk.gammaABC[397] = Pairing.G1Point(uint256(0x305d746537ad7bc6e0c49369812ff21e7ec5197fe01e8f4c05412eaa02c81796), uint256(0x2198c35b5d322257a1dfde7e8105a6f7f7a64b549ccd156bbb5ea3bed97a0cc1));
        vk.gammaABC[398] = Pairing.G1Point(uint256(0x2ba4185eb174612fa289e67f24e27b2536867cb4bb7f7973b8e50c14b31cb1ba), uint256(0x1e3949406b45b78ccceef5ab23a8c60a425861fcf9757a9a40c7141c35315d24));
        vk.gammaABC[399] = Pairing.G1Point(uint256(0x1615f2fb9fb5e8d64e8551198c425eb2b449f00ff85fdbcc6bb499e6f2e6a21d), uint256(0x2302240b12963b29a7d25a090754c52a4dc5e5eadfbf10c239e057246ee68685));
        vk.gammaABC[400] = Pairing.G1Point(uint256(0x21649e28fd3929094aedb74f9420145ab6376701191817636992a73f06408b0d), uint256(0x06cc1497e8f21875bad42005c03ee521e6040ead9d740ab0d2beef84038065e5));
        vk.gammaABC[401] = Pairing.G1Point(uint256(0x2f42ca44102181097c1d93f85f4ce8eacd19dcf49ebe4dcb34236db2926719cf), uint256(0x04ab4c39b04c6d6101e0f6ea6df0c482bd88e14a31be5cb71114bab90f5c3e38));
        vk.gammaABC[402] = Pairing.G1Point(uint256(0x2f15427bd08c2d8ff3ebbc8a6a5275435aff3320dc7169325204fd160960b99d), uint256(0x1b4c42a5179feeb0d985334afe4db3e3e4a9f31328449794bcb667caf4df50f9));
        vk.gammaABC[403] = Pairing.G1Point(uint256(0x2b1e31521641b9e4e9f3aa92db049403c697af6a6d501813b1247e98161101f4), uint256(0x15bbee81ed61805c92e0898479a8be60d31057a474631bf2fe2eb32569f913d8));
        vk.gammaABC[404] = Pairing.G1Point(uint256(0x0526afeadcb6aad99053cf02478ad1a47c605f7b6fc42021bc8232a36d46c841), uint256(0x1e8779ef2515fdcb07ca6115a0c64ea4de3d19364067e2eae9bb2ebe8b158943));
        vk.gammaABC[405] = Pairing.G1Point(uint256(0x13b47cb2847bae6edc0083e5fb8b5e7e455672752d46f6ae7b90fbacbc4ea4de), uint256(0x1722d78363e9e3996e0776eafd9731d012ceb98609574a5caa54b2b22b49bc9f));
        vk.gammaABC[406] = Pairing.G1Point(uint256(0x019ba1a9e40e2632952bf0f83f280f3324ee7d13b3a95d9bd9c97ee9bf3ed8a9), uint256(0x2754cb9838593fcc4c34c22bb26ea3129e381dcc67515546d87d3438bdb31ca0));
        vk.gammaABC[407] = Pairing.G1Point(uint256(0x1ba76b13a6d1a05675778afd5f9af1ff6c27436794d6a95be366cb939d3031c8), uint256(0x09abd3a93f6921b8cc266eec47349e910c090485b9abf58f6fd5f47f267d8602));
        vk.gammaABC[408] = Pairing.G1Point(uint256(0x1e3975f2200a14604c83fa083aa06b300ec9eb32833279afd547df92fc507b62), uint256(0x09e3b987a21f566b5c14fa5013dd93027af1ca2df426f146e5df0402c1e7b36c));
        vk.gammaABC[409] = Pairing.G1Point(uint256(0x25ba694cbc56e015289f6c93d9dc109de5beb85102ea42bf0b0feda8ba81c4cd), uint256(0x0c55b3865b22730b8a6e7e66d1ced341b482a5519969acd63c63e7908d58421e));
        vk.gammaABC[410] = Pairing.G1Point(uint256(0x231aa6372437c318f0cbd7b774253eb6c4d983d83e346c134d3be4c06f1e07b7), uint256(0x11f5d682dd117c5cdc9366aeba0b685342747dc3b31703ae537ff49d2d6845b2));
        vk.gammaABC[411] = Pairing.G1Point(uint256(0x0f337a8797ec48e0bc66241713ae4e2fc546ff980962927b9a580ae2c617d0f2), uint256(0x022dbdd7be8113cd59c56774fd893e0e27dda6c6a1eba80347e8344bbbe0e61e));
        vk.gammaABC[412] = Pairing.G1Point(uint256(0x0702c706ac21e9c32fc16c8a7c059b2d9f4c1077523434317bd152cd3f1bef62), uint256(0x1e93030c31eb66de2e574d6ad90cedfe612b196bd25063d4c90c962406dc6592));
        vk.gammaABC[413] = Pairing.G1Point(uint256(0x117ca9559288abcb5047eac7000212b8fc6969254d7b124794527ae3b061dd9c), uint256(0x1001d6949ea5023cc1fcd159071e32624df45f0669af162b39e4142435830709));
        vk.gammaABC[414] = Pairing.G1Point(uint256(0x076bfae8a708f453ce1d34929713de8ae62572d232b10f81d7a554754049f09e), uint256(0x12517d861a30f686e0a4258b1a9439a19000aa873cbb0bee5eb91fd39920a98f));
        vk.gammaABC[415] = Pairing.G1Point(uint256(0x05c0d8926ca3739b56b75e8657262d6f260a2103a33a0d031b30b5629e5a3a4a), uint256(0x205b12c3d45f0ff38d7004a9209dc292c67df81e3f37a7482f9420823ca17547));
        vk.gammaABC[416] = Pairing.G1Point(uint256(0x24858120889b07c1ed0241b2d4e76ad2af66fce57dd39f05c5f2ff2df4029211), uint256(0x0a0391986072205da4ab24df664676c5f121b9d33ce137d8943b8245062a7eda));
        vk.gammaABC[417] = Pairing.G1Point(uint256(0x2a042c5c65abefea9c59b56c44614da0457cf2f494e83c3e39d0ac08e0f6643e), uint256(0x0686f969ba8d5d8a869529e324b42ae5deb7ba88b09fdb46b62d4a48e8d07dc1));
        vk.gammaABC[418] = Pairing.G1Point(uint256(0x1b2215654f83d5722b06daa6ae1030c9e3d6a14c9e751a628434f937b211848d), uint256(0x1285336629dd3cefbf56a78578a234625b110b4da7867f5ff802b014b2ca832a));
        vk.gammaABC[419] = Pairing.G1Point(uint256(0x0d13151f5fb6a46a483785a47ad2af4fb830ab0eb81298874a4270724d0d5329), uint256(0x2b92464d77ed871a406ada9ffb497bc2c34ff560723ea45f25e878c40f8f988f));
        vk.gammaABC[420] = Pairing.G1Point(uint256(0x02d266e77f6b77ebbb02e3281a006379d1058605c55913220a8428aac2c0c6c7), uint256(0x2b91c8f8a58edbf597592afc103617c54a25ec78fd53f3dd354038ef94923490));
        vk.gammaABC[421] = Pairing.G1Point(uint256(0x0f020d0743b9200d151602aa537b2d0ba8f24f61ed794d9e3f96c60687cbb241), uint256(0x1608b5757078cf1782913924ad5393f1725856600bcf83cc667434c94ab56c24));
        vk.gammaABC[422] = Pairing.G1Point(uint256(0x221634fdc0451d4a694bf4635801804e4732066d1510c7a1f2e03856488213d0), uint256(0x11e68111cac329696e0c52b9c3f4053aa8b44ca07fc85f14a512d9a2376db786));
        vk.gammaABC[423] = Pairing.G1Point(uint256(0x071ccc02ae5a7a72c2cc0e4ee59201d128604007e5051569d0645e516723fdef), uint256(0x22c83395a2a5a566f9805744bf9ffeda331fdb24b8394b5957cd7184d8daebfd));
        vk.gammaABC[424] = Pairing.G1Point(uint256(0x095c09ed7ca5f413f9301cb69005eb886ada4e599357d07c35f90b130b8fc108), uint256(0x0666bd659e3be2c619c53464e8b110bda2b61ce7db373dc491ccd2aa22e06445));
        vk.gammaABC[425] = Pairing.G1Point(uint256(0x0802b346b130e930d875b84f67bfd1ef1d7b2c04f3996bb8eeb237c382fe19ea), uint256(0x2368adb1da1de680e9cf7478a3b025ac93a3a9a3cba63c3e2c2e1e85bfb57908));
        vk.gammaABC[426] = Pairing.G1Point(uint256(0x1128b803d1cca3fead743120262aba47f6a02b801d3328ee1de5e7127a2d6692), uint256(0x215dfa593a72aafb72c65c5cf707436e99e50dedb63251f8bf8d3dfdbd5bd30b));
        vk.gammaABC[427] = Pairing.G1Point(uint256(0x2d0839ef711f5f140f0bfc86ab84c3f24a8de9a00e8824e1d910939b64e03639), uint256(0x1398d28164434a9c6e9db430bfcc1c6a220b113072fc26b5d62bccef48c2e784));
        vk.gammaABC[428] = Pairing.G1Point(uint256(0x26ffade18e7b787a880b6340ea141f6e4355ac4e395c64707271f24e0daf22d5), uint256(0x191449cde431c1b8e3b110148b652b7e09e4c520cf0e1539756da4ae9ab378d2));
        vk.gammaABC[429] = Pairing.G1Point(uint256(0x0d4d04d4e98f3877d96e309fd976979752ae9e65ecc089ca9f4f0062ef91fd4a), uint256(0x2d89d929f918c3c9e6eff03fcae96753caf7690181468e083ea8740d1662837e));
        vk.gammaABC[430] = Pairing.G1Point(uint256(0x290e8090ec1a19584e97e0cf247975a9b09dcfb94752d0d8ed14f54639a197b8), uint256(0x0543f29d40e9f929b4778969a8f8be7c50bcaa17996cba51709e4ee1c05b1247));
        vk.gammaABC[431] = Pairing.G1Point(uint256(0x0d490a98f7ab93f79f1e6136aa0d9d3659d46c7cc02cade208aba54933b610eb), uint256(0x2ef3e074506e2ab9c5779edb63aa384ca83b7777a67ae3d9416077a0470746ee));
        vk.gammaABC[432] = Pairing.G1Point(uint256(0x1b53dd487395513e348acbfe61ccc1c6f31051d4c4cdf7df9f2aeb217456f62a), uint256(0x173a689c7aa0394114047aa832b0c11edc0d440d6380cc8864a28f8f53c3e84f));
        vk.gammaABC[433] = Pairing.G1Point(uint256(0x0c3ad9357b4dc2fcc0bbb3582f9b56fbc54b3cccd539982da0889f1697a74b96), uint256(0x1d07a7a9146732d65ebcf845a7a38ae681d9047021e88e6453232a7e0fd66ab6));
        vk.gammaABC[434] = Pairing.G1Point(uint256(0x1ea1965e8c6c156c7a9f493de78b090e2d4072ce93f55fc8dbc9a6f3da766d9c), uint256(0x06a014a414f4de12b09791c93019ea330634b29243dfa37a1759e0ac419ba1f7));
        vk.gammaABC[435] = Pairing.G1Point(uint256(0x0df8200948ca24145e60ac57ba498567135505a5c96ee6de865ec2ee2e729ccb), uint256(0x2d0fcc4f9e48d2788a20704ae401c4021ffb03011655ea957dce94e44e15684b));
        vk.gammaABC[436] = Pairing.G1Point(uint256(0x02e48921bbf6fd67711173e47b84fc09c90c1e09c239272e6cb6620ece85057f), uint256(0x30629dde0bbedd1cb99511d56d781fe2bb5472bbeb1308a9ddfb749242a44047));
        vk.gammaABC[437] = Pairing.G1Point(uint256(0x1df3af6b07a5939be6c5461d3397d8f0a3621f1694995274304fdd82e34e1b6a), uint256(0x1785f64743004a67bd3c552b27208392efd09697e06419f02632e20eaa90fd17));
        vk.gammaABC[438] = Pairing.G1Point(uint256(0x1d7ac347bef7a5d89336c02180ac807d0a36b1d8c02338e278961924a6e0d6b7), uint256(0x1bd164727555c50a78c05616dfd4dfdaaaad91897225db2bb3d6447bd0d285f8));
        vk.gammaABC[439] = Pairing.G1Point(uint256(0x182cab016a622b785b1b164851691475ade8907744bb7abe5987b71fa3e90e27), uint256(0x07f88b7a49ddd1fdb4b764a25e036337c5baa58a117afb1525c5faccf73c761b));
        vk.gammaABC[440] = Pairing.G1Point(uint256(0x10312134f95e6769cc6f2a52ef1a59bedf27f50c0df0dc91860d94eec9178912), uint256(0x1e9046d3868f588d36ec8a7056d29354d5f97fd1cc50616befaf0381cc29cf57));
        vk.gammaABC[441] = Pairing.G1Point(uint256(0x292d59665c77c474beb038fe68e830c9b1148ffb848d86a94f0e01b733e90178), uint256(0x2d66afd1ecc193d75b65334913127208a8396fa91708fcbd8005f96d4b6357bd));
        vk.gammaABC[442] = Pairing.G1Point(uint256(0x2a9f53a593e36890eee09ebcf62e05905d3f3c2c52c75a20e47701727610ba59), uint256(0x2d687aa0ac35221e467a7c5f56f1d5385dd9eeb5d4ae2073919a8b47f96b6cf6));
        vk.gammaABC[443] = Pairing.G1Point(uint256(0x1dcfcb43527988d7c8c0749748e5116f7c9d354907dfc009d09df86f7db46c8b), uint256(0x10a1e9bed34b6bed349479404915e92ceaf4c3cbfc2c2637adee31cc1acdab5d));
        vk.gammaABC[444] = Pairing.G1Point(uint256(0x2bf2833d658e7b37ad0e25ded4052fb427f329f96420b1ca1f11f5227b8c9677), uint256(0x2015c9e425cb5d29085a1a3ea9d0e7ce5bd1513b9b414cdadbb7fc3f25b85773));
        vk.gammaABC[445] = Pairing.G1Point(uint256(0x01f4e0fca58141e0f80191607ddb0242c5514abed997f9a6561ac7bde2973ac2), uint256(0x14abf57f7df025e12c7d7ece58a284f768c0307b177910cf227cd0dd4394d57a));
        vk.gammaABC[446] = Pairing.G1Point(uint256(0x101131cf8f123e13da157a4e5ed6d7a3c945a996a65d0250957b703b2634f6dd), uint256(0x10028a4373c887ab9427cc5fd95d72771930e6c19d2c29372630878b61c46c21));
        vk.gammaABC[447] = Pairing.G1Point(uint256(0x256e08d5139b868398f4a81d60d4818d5e8ff6d21272cd0a8de7e787b113a845), uint256(0x1e75799fa3fd5bfa9c9bd55ffb414f15c2d4a8f83ac278f0ccf682b3c8751591));
        vk.gammaABC[448] = Pairing.G1Point(uint256(0x0336a088a5d6331b6e2b2f683a2e48d269aafbfa363bce89311b755dba86590c), uint256(0x14dd789ef31cbb5b093b761187fae01a86435c8b696fce8c7618d666d5e2a494));
        vk.gammaABC[449] = Pairing.G1Point(uint256(0x0f9cf001e93a335fc115943db17c4e06d30812c40bd63eff97a8c17a04143159), uint256(0x201e41e3d4a67175faf503931383228f110ac515b2a0464452d0782d0368e30d));
        vk.gammaABC[450] = Pairing.G1Point(uint256(0x1e6c72cb5386c3398d34ee39c57aa4e68c7fc64ba6eaf9192aff67e2949fca70), uint256(0x0c86abb5c21651823fd7b429615a0338e176eb77c0a236ab2eb53b3bf877f138));
        vk.gammaABC[451] = Pairing.G1Point(uint256(0x1d377d8030d4068e6e476a3fefc93a5c910b46947a3406fbed8d67251091a89e), uint256(0x2ba74db404f9d7e7ad4ae6f54773125bab9538b708041cb0bd4ccc1443c75aa8));
        vk.gammaABC[452] = Pairing.G1Point(uint256(0x2b82717b4b4bbb86ee28daa44bd87285998fe7926e447bd687508c2e2231cf35), uint256(0x1bf8d064d450e5a813ddf88287ab4859eb626d298bede88850c1d9c14e7e0695));
        vk.gammaABC[453] = Pairing.G1Point(uint256(0x12c2ccdf95ef48c242818c009f561e8223a7521307728d597fd4df6bb29272cc), uint256(0x1342b404647eb67e4f046a440264132e3391043e321bb142f012aee56217c688));
        vk.gammaABC[454] = Pairing.G1Point(uint256(0x03aa359904790ede83c57f79e1979056d62f0c3a3c4ec19ab12729f1c896ae5d), uint256(0x00aafe776f7502e9bbed1ec064aa119663a1d88873da773c76cb3422b4834f1f));
        vk.gammaABC[455] = Pairing.G1Point(uint256(0x1cd4ff73f54700e3dd5b26e29751134481c5251fd619f8aca7d731e9f9dafe98), uint256(0x1ff192cf79c89c71001d77ad8143046dff2b96e2027789e069f6f831706cbcf1));
        vk.gammaABC[456] = Pairing.G1Point(uint256(0x106b696402e23e9da05957f2a66f87bcba446edcbc5d0a18870a8621bdcfe75c), uint256(0x0eec32979819facfd3d22da25291872ff3ceda9a6c1b1df1acec91c9414bac52));
        vk.gammaABC[457] = Pairing.G1Point(uint256(0x3041ba2bfbe85f1e6ed2f44b3a6a2e987ed134954b36a4266156afc94c195a68), uint256(0x05bfb23fb14bc0ecfb65662ff2263fedbabd90d13a7c3aa908ec5a2b5f1ce9ef));
        vk.gammaABC[458] = Pairing.G1Point(uint256(0x190d5930cf391b6a46b97590b44be906cb7cc15f4108e2a806e3263dc1033ca6), uint256(0x281136b28ec288153fd272a9f3fced2cf213dcdac92a27a5d7e8aa67476f485b));
        vk.gammaABC[459] = Pairing.G1Point(uint256(0x2ace19e2e5d19c9c1520173c64763860c9b0f52ba7a3fa3bb70714e373d43d1d), uint256(0x1730298cdebbd4e943cbd6991c93960dcb8d35aed31a3af1bfe81f63dd8aa2c2));
        vk.gammaABC[460] = Pairing.G1Point(uint256(0x1a2453d56e3a299c17693e44130e1916a3611511916f5e3727946483acadd896), uint256(0x2c0df4a4117b9d81b8b127bab5ebc1a2cddf69007058449741e3754e181c9188));
        vk.gammaABC[461] = Pairing.G1Point(uint256(0x12f4ce50eabf0e8954f52d5e0c961bde0f4ee9af967e6aa891efdf33e34782ee), uint256(0x2871102d422ee99d9efd74cdcccd0299ecf947374a4a7763c5988b298c84835e));
        vk.gammaABC[462] = Pairing.G1Point(uint256(0x1287443d3b6e2b877facffdb55d0890b142e8fbed79910707c848b547c9f9465), uint256(0x00ca6e1ae258e4b2aa11d8048d2420e15ef2d19ebe80b6675144b3a43756665b));
        vk.gammaABC[463] = Pairing.G1Point(uint256(0x0e2b9a5bcbd40b61d460bb1a194b44435fc508d23ab24fe44dda4a1098ffc5aa), uint256(0x1c1b1bf39bab63c95640fa88cd8853499aa439dc3988583394fb8c9321200806));
        vk.gammaABC[464] = Pairing.G1Point(uint256(0x296e81d736216b4a66863311d76214ba2ca7be1dec6ca972e966285abf00261f), uint256(0x100bbe7948d85b673b6c6a03301a3bda90227d3e79af24129916d0c4b34d0a0d));
        vk.gammaABC[465] = Pairing.G1Point(uint256(0x1ee0feab53f88b7d8c2e77b5c31e31f011ea39a73915f04af190cc9a611b4ff3), uint256(0x24de6220b9eab54d8e5013dd0193ece5e939c6355f3d057ecef1408dec90aeff));
        vk.gammaABC[466] = Pairing.G1Point(uint256(0x127495107c1fbe28e6004f657d530386c683ade0de7f5116353a7e7a362b7fc1), uint256(0x1b48e89f86cc2cad24096290cd6f2948889ac49866a9e9d5dfb09193a19b6930));
        vk.gammaABC[467] = Pairing.G1Point(uint256(0x112bf9170786cb1d42cd115ec985795fdb5c1202bff03d16f8c6f70cb37e28c0), uint256(0x1ec5aa5c93b5009473841b1378bd5f0cd0ef46ae1846a04b06cbeee036ab0547));
        vk.gammaABC[468] = Pairing.G1Point(uint256(0x2e9ebd06d72627312d47cc8fc4ccbd67c2f4145f7912f389005e012ad49baa4c), uint256(0x16aea13ef3feeb59c276e0a1f1d5c57eb0e0c80fd87d190b82bfbe76a978dd5a));
        vk.gammaABC[469] = Pairing.G1Point(uint256(0x0b3479b918f785dbc15d7a9dde31eeb1ba3a3980269834c024d63cf659ade147), uint256(0x03cb44784e1c2d5a5bffe2ba526bc42d629be37b5c56409f840c6418814ec9dd));
        vk.gammaABC[470] = Pairing.G1Point(uint256(0x1fe3b26af696380192f98578f1b2b55f52a3ab955df57d4e63b4f0f68415e034), uint256(0x06656e56322d25229462e32be3bc58898945012b1c96490dc9771ab7314a2c30));
        vk.gammaABC[471] = Pairing.G1Point(uint256(0x1600444f7ec4ca85c2e7d598ce2c4bc744c0c683591549e0f98abd986cb2daec), uint256(0x258d946d727ab0718dec256ffeed27b14005244af631bf6ebdd4c49df8e0d542));
        vk.gammaABC[472] = Pairing.G1Point(uint256(0x14b15516317560f804475b385bcb497445a295a521ebe2ba1bafeceb6ef79ac7), uint256(0x0586c483544712ec60dd5c55b653e67870e4a20eadb2958bf9bb3bd992c8c55a));
        vk.gammaABC[473] = Pairing.G1Point(uint256(0x0a616eefa2bbd2ea47b99d153a09da5f9054fb163d9048a533ddfc7335027048), uint256(0x1713048df81c04924d231469b3e818163cf7ebe94712706bc1202e6e3f47c6f6));
        vk.gammaABC[474] = Pairing.G1Point(uint256(0x02766d0da8c1f4d0795dc1ddb2768b49b26658ffa57f4c06a72dbcc38f2a1829), uint256(0x29184470459216556f4eb6c761fc3297357eb82ca9e6692440fda6c451e3fb5f));
        vk.gammaABC[475] = Pairing.G1Point(uint256(0x24878311d87ad2c120a2b23982878ec382a68a74e8d2c9b93000ba4629deec21), uint256(0x2a051603b7a587215cec8a48ffc4f7ef5437a6ab2a341522b34fe4e501b1db2b));
        vk.gammaABC[476] = Pairing.G1Point(uint256(0x17a4e472677c0a2d01a199f53de1d8b822a86cf9b988decf314d2f6e01cd0e88), uint256(0x146ed4cd37d71cb1805b68bd5c4a7d0a5f3f1a17dcc5dd97711c00e6b1c6b5e2));
        vk.gammaABC[477] = Pairing.G1Point(uint256(0x143988b9a1dc99c6a5ba38a3e2a86835b7c87c77603ba2dfcab969210706da78), uint256(0x1a768472865957aa343d0edc5cd58d7e4cefd0e310af404159bbad069d6111aa));
        vk.gammaABC[478] = Pairing.G1Point(uint256(0x1115f11ed842e70e5272b5c4adfccd6ac0322fed7b80524698aba5018c9cec1b), uint256(0x1c607e6f00ec671db8e4f3547dde84a0a1dbe292270b25e2de3ead97f81d4e93));
        vk.gammaABC[479] = Pairing.G1Point(uint256(0x1dcab27e1bee4200da0fee3b135f674d778dbc8f874b182d0962f855862b3ac1), uint256(0x0fe24ab160eb67f412a66d5813a28255b6e10ed57bae10149585f8a14cd16903));
        vk.gammaABC[480] = Pairing.G1Point(uint256(0x08d054bb01e41062641b209c1851e13d9085eeb757fabfffd6a0db69f23345fc), uint256(0x04d9992da4f9573e7fdf589999099bbf3f96ff310e9bb9e4136273529ff66f3c));
        vk.gammaABC[481] = Pairing.G1Point(uint256(0x1c09338483bc5f35880a4757a34119f10dd37155797757742ff1eeaf479f8a83), uint256(0x2910d64db851a9a42120db7e01c1b1474eb5c9060946094d47a4c4c74eafe848));
        vk.gammaABC[482] = Pairing.G1Point(uint256(0x08ed25745f0940303a66ec6861cdd4821735cfeee066027ae975e9d4d9203117), uint256(0x11d8233ad26985ef8792474176c771af49363a96d7d44f8701fdd11d69649bf7));
        vk.gammaABC[483] = Pairing.G1Point(uint256(0x23383c94ee2dabc7ab2d0264d4ac4579b64eb9b176d50f44b2de374d860e3c50), uint256(0x00e5c937230e99ad053c6c689363994a2001d994c1178eebb720fbe1342fd39b));
        vk.gammaABC[484] = Pairing.G1Point(uint256(0x2c54eae8337ccc483f2383b073837382e23ebb3d013d5a706b27a8f2d82e87aa), uint256(0x0ebf9a8c7623bd493a662d6a978c69b9d738e65b52d73418da753cfbd8cabd82));
        vk.gammaABC[485] = Pairing.G1Point(uint256(0x02a15643299a86a28ce5a81ce7b03aac50e077ade05e2429bf07aea294e31324), uint256(0x18d09000bffd71227b7175abffce592823ed69707a3bca693a16c1f12e4ba834));
        vk.gammaABC[486] = Pairing.G1Point(uint256(0x051566c4bfd86ebfa73eb9ee4a845136f559fcfd60c1e480437dfcc6b6dda379), uint256(0x1e7917e38d50d5db3fc9137d61fb8ee007bf2ebbf3dfc81ab8e3f88348b5ef57));
        vk.gammaABC[487] = Pairing.G1Point(uint256(0x1b106fc1db5ef2c933db2fabbbbab62b46270b792f08195318830afc61a41c20), uint256(0x13d24a4686423055fd5eac14032de014a9c122852f160feaba7f9f43a3ab90cf));
        vk.gammaABC[488] = Pairing.G1Point(uint256(0x26211bf3b407235e504a9ef339508b9be0e8e5abc7ab52baaebebe97347aa8e0), uint256(0x1f6cd01a6314c563d68fda0ce8cd625e04b67d0f7f27d1bd63bb0b5f38866c2b));
        vk.gammaABC[489] = Pairing.G1Point(uint256(0x1434043103e45bf83fcb1902d33adf98d815f9d27fbcfad0299d65def2083b8a), uint256(0x065563fc22eb915dbbf69df6ad86a620273ee2c52b2efbe16a433a910794a1e5));
        vk.gammaABC[490] = Pairing.G1Point(uint256(0x0d3be4766e6633061a87b264450f368f6eda7e541c95cafd631a9b8ea4da1504), uint256(0x214733c4014ca5ce00d70f3e87b1354c493477b9be9bc09b14890865a1920a3a));
        vk.gammaABC[491] = Pairing.G1Point(uint256(0x2e468c1cbcf11458a349dba0df72f1473859790e3c1522c76473aa68ff51edf2), uint256(0x0e29dfa3220802b4db1f1824783981da4f3a64daac1c79fd1f2a0a7151c560a8));
        vk.gammaABC[492] = Pairing.G1Point(uint256(0x1a8168efd6776c39fadf455ec1b4317443a33539f45429c274571bc5ee06cb62), uint256(0x15014b77d1628e4d97edc2ebd66968bb2841e56bee97630645c3e8a9df95b7a8));
        vk.gammaABC[493] = Pairing.G1Point(uint256(0x2afa3c16e1594ac50eb159c709ee635efac97f10863280f5471bfd601483f1d0), uint256(0x1ebb6a1107f91c470da70d3f46ce258b67e6cf81f7f0441ef422de200c1db6cb));
        vk.gammaABC[494] = Pairing.G1Point(uint256(0x0000a309225d12740e947ba0c9ace51207febc6afd65446021b3a4d44b9800a5), uint256(0x250111b53345dacc0bc8ca49ce5e39dd7a2a3aeef568f7a3e914d7cac8b1c29d));
        vk.gammaABC[495] = Pairing.G1Point(uint256(0x13fe0fa4a7e9b3b9884f2bab75dcba9f24258573470171776720ec9fe0c0e9d9), uint256(0x1ada0cb2d020fb54467adb356ba4fd60227d9916d2a6592d1796b3ee0f6d03c2));
        vk.gammaABC[496] = Pairing.G1Point(uint256(0x19fdf93e480026986b39a51d81dc0005d380423f1cde207eb0fb0bc3bf9c049b), uint256(0x119b2eea5d3f5a9475e5f7838bb47e4626d19ba06a97c7067cdbafdf1213820d));
        vk.gammaABC[497] = Pairing.G1Point(uint256(0x2b5a70d9f941d9e85beaee558a80b2355ab4e0660c8666453093522bfc943f7d), uint256(0x13dbdd5d6ec3ace704d53ae8244e45cdbf99b93d1c2961a0bd68e24d784e8b5f));
        vk.gammaABC[498] = Pairing.G1Point(uint256(0x095e4dbaee27759585b53bfb72915f3e587c4e3ec7352b931e3b3ce01a20f539), uint256(0x1dbd7f28cd218c9ce66e3a871acc623f8e7395e4e62fcb2cd386d15c471561ea));
        vk.gammaABC[499] = Pairing.G1Point(uint256(0x0e05a0949d9d8b197ba410e713fe8252dae7f32a016c5750b1e72283b1963a30), uint256(0x0882ae0b850391e7e6df63b03933f46f874936f3bb369b4147d7c2801b7b408d));
        vk.gammaABC[500] = Pairing.G1Point(uint256(0x075f8d95e07b093be57d0488d72fb5aad72d8b22480477cca10c00fd6e3a588a), uint256(0x044e0e50862be90984939f331c4903db5e8362226aa20b22cf155db9c79458f9));
        vk.gammaABC[501] = Pairing.G1Point(uint256(0x2de593dc3cd0f0afa35304e5d52cffd10d4d8c443aa221805bc6e2a88f848976), uint256(0x15f06e5f0461bbd666393effd86d72ecc3f82a87c94f3dc3671fada1d6e0522f));
        vk.gammaABC[502] = Pairing.G1Point(uint256(0x2c1150d77b350fba2d4d3a17dce0e0bf87c57c6fbc718009d98864d4cbddb56b), uint256(0x129f230ccfaf4ad726ba8145d9536b3cdf5a1fd168802d66103c09179e1d5837));
        vk.gammaABC[503] = Pairing.G1Point(uint256(0x042e97761e19cfb6008d86684eb534cade42b9374f05ddf783a3c19bc5139de9), uint256(0x05d65f5ac13d0638f587c65a1698903cb44535bc63e40367cb344a69735fff35));
        vk.gammaABC[504] = Pairing.G1Point(uint256(0x1e5be9245a3b4932310eb081c53520c48b6c68baed98b12d748622898864c9c1), uint256(0x2fcd02e33a342704844d3607f6b552245d95db6a8a3a577417d4db8c1adf0b1f));
        vk.gammaABC[505] = Pairing.G1Point(uint256(0x101972f488153f330575e03dd1cdc103467475003ad4591ba78a1067e31addf6), uint256(0x09b7036825194a3a01cb74bead7179084032edfd9816a4fec16c18edbb7866ed));
        vk.gammaABC[506] = Pairing.G1Point(uint256(0x0d6c800bc63bdd8b8c4f9037b7ca102af5c5c569747110f3cc744ebd4df14b46), uint256(0x2f618c0b85268a6b095a287580c1e428272b9ddf4622c5b68a359ebb30136de0));
        vk.gammaABC[507] = Pairing.G1Point(uint256(0x00ca58120fb9235c16a5feb14e1bb453fabf2d40e06fb29878bf652dffe309c0), uint256(0x23942fa40031378c47b2b7766791e09baf22a07c6fd30d7e9d6f754b5d1d7ba3));
        vk.gammaABC[508] = Pairing.G1Point(uint256(0x2e9867fe126accbefd50ae5cfb5efcf920309b34942651cd0e30a108228ed848), uint256(0x00ea4573290d35b7d9f5ec38b2a0444a23bff42419284dfa4e878e96f44c2f35));
        vk.gammaABC[509] = Pairing.G1Point(uint256(0x18301ec9276c52caca2674d099086d39cffd5aa94623ff304757c71d2f0f07ad), uint256(0x244ebbed0d008e4298a2dcccd4dc44b9020e5cc829f61c746a56901b26875213));
        vk.gammaABC[510] = Pairing.G1Point(uint256(0x2a9d4e361219187a94097f74d319264dc67c9de8cf9cbced5165f390ae7ad5ab), uint256(0x07f783b7b630c127c832f45d696fff5863b70127b947c7e96730c1185d05dc0a));
        vk.gammaABC[511] = Pairing.G1Point(uint256(0x0ab36292baae3ac20cd8b67d6b624f95350c7ca6627f7e01bee3ee9897433d69), uint256(0x284de295e48b7b947d22a8cdb81faf9b52753d3eff2cf6f49140f039b3965bbb));
        vk.gammaABC[512] = Pairing.G1Point(uint256(0x28831f4fba16e70758e736b093793127e42d008c66c6be3efab9d52c7d5e9e4e), uint256(0x1a719cf1001a05c583e73d78a8c06f8197c360565392349d523e5a30ee856d30));
        vk.gammaABC[513] = Pairing.G1Point(uint256(0x15e1c9d411e28ebc0784394c4ab91408ce5ee8f3e1292dcd790c80523d3e2733), uint256(0x185ae85732afe6312d2424977b15c4494c481450c0357bfffcd039a14e7f65aa));
        vk.gammaABC[514] = Pairing.G1Point(uint256(0x0caa5bdcbff499d53fb7826e413903b2caf26611f65b7e9865d2a81cc5552fea), uint256(0x25f67f051982f9a031dff4f767052461959f2a5c4101f57a733d7b10825695f9));
        vk.gammaABC[515] = Pairing.G1Point(uint256(0x1448049ef42de58ee0a1486d50dd0a901c8dda105b3a3c56b5d7ef1f435bab7f), uint256(0x30022c6406fde265cdc0c79a858f434c69691d076d4ff98166ae8277bdb26fe4));
        vk.gammaABC[516] = Pairing.G1Point(uint256(0x04ecdce782d5912795691f3dad9f7d9d4167c638de7381b4419f5c65d0ec0ade), uint256(0x0aa240208de4e4b94bfa004d575edd6dd9b9c4b47a63653356eb23a26c8a9df4));
        vk.gammaABC[517] = Pairing.G1Point(uint256(0x02fb233235013fafafac76664635e82a2e9dfbe1459e3a7f2c87222396355598), uint256(0x1683edd9974687db3dc45814ad54ef547d0835d6cc9eb5973de42cc1cef73686));
        vk.gammaABC[518] = Pairing.G1Point(uint256(0x1e08f2c835e72f5e8e52e1f4b50ce5bc02cb108201e4be442ad725e4b6692447), uint256(0x25956786144f120151ea5b9ec8011617b5139243f67e68c7a2de7c3d2bd6b0db));
        vk.gammaABC[519] = Pairing.G1Point(uint256(0x0bbc8d0b69c7d4ef778f3df8aa6ce3d3c3cbad2b6dbd19fc454f44ed00837ce6), uint256(0x103ec87da04acd1e29727add35a7ac14049d76186b073d46349ef5b74a4890b8));
        vk.gammaABC[520] = Pairing.G1Point(uint256(0x231a87ac13b065aa16cf70596ec3c535233cd80379189895563cda1b83acf499), uint256(0x297627f7dbd8d30b5e0e23773533519e5f3575d0533507d72837433ff81ec251));
        vk.gammaABC[521] = Pairing.G1Point(uint256(0x1c6288604211c09b9c37b8a2cc81ac606b18bb2831b7484df99cb12e1a537c0b), uint256(0x13b7a3bc1e5f97d16f27252cdb76b2667882f6afd15cfc1c26e16d9e68d047e7));
    }
    function verify(uint[] memory input, Proof memory proof) internal returns (uint) {
        VerifyingKey memory vk = verifyingKey();
        require(input.length + 1 == vk.gammaABC.length);
        // Compute the linear combination vk_x
        Pairing.G1Point memory vk_x = Pairing.G1Point(0, 0);
        for (uint i = 0; i < input.length; i++)
            vk_x = Pairing.addition(vk_x, Pairing.scalar_mul(vk.gammaABC[i + 1], input[i]));
        vk_x = Pairing.addition(vk_x, vk.gammaABC[0]);
        if(!Pairing.pairingProd4(
             proof.A, proof.B,
             Pairing.negate(vk_x), vk.gamma,
             Pairing.negate(proof.C), vk.delta,
             Pairing.negate(vk.a), vk.b)) return 1;
        return 0;
    }
    event Verified(string s);
    function verifyTx(
            uint[2] memory a,
            uint[2][2] memory b,
            uint[2] memory c,
            uint[521] memory input
        ) public returns (bool r) {
        Proof memory proof;
        proof.A = Pairing.G1Point(a[0], a[1]);
        proof.B = Pairing.G2Point([b[0][0], b[0][1]], [b[1][0], b[1][1]]);
        proof.C = Pairing.G1Point(c[0], c[1]);
        uint[] memory inputValues = new uint[](input.length);
        for(uint i = 0; i < input.length; i++){
            inputValues[i] = input[i];
        }
        if (verify(inputValues, proof) == 0) {
            emit Verified("Transaction successfully verified.");
            return true;
        } else {
            return false;
        }
    }
}
