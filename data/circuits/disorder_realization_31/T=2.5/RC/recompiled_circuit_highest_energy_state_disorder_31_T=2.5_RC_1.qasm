OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.5950522) q[0];
sx q[0];
rz(-0.15931436) q[0];
sx q[0];
rz(1.2567318) q[0];
rz(-2.9873084) q[1];
sx q[1];
rz(-1.0839387) q[1];
sx q[1];
rz(1.7854324) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0935554) q[0];
sx q[0];
rz(-1.1021656) q[0];
sx q[0];
rz(-2.8891536) q[0];
rz(-pi) q[1];
rz(0.21535318) q[2];
sx q[2];
rz(-0.15654187) q[2];
sx q[2];
rz(-0.48605824) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.4843085) q[1];
sx q[1];
rz(-2.2806232) q[1];
sx q[1];
rz(1.4906989) q[1];
rz(-pi) q[2];
rz(0.47812652) q[3];
sx q[3];
rz(-1.9476796) q[3];
sx q[3];
rz(1.865975) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.1821182) q[2];
sx q[2];
rz(-2.0017767) q[2];
sx q[2];
rz(0.095890447) q[2];
rz(-2.828756) q[3];
sx q[3];
rz(-1.0597119) q[3];
sx q[3];
rz(-0.50725168) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4023912) q[0];
sx q[0];
rz(-1.6283789) q[0];
sx q[0];
rz(1.8550523) q[0];
rz(-1.8402428) q[1];
sx q[1];
rz(-1.8845314) q[1];
sx q[1];
rz(-1.7110862) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55670083) q[0];
sx q[0];
rz(-1.624361) q[0];
sx q[0];
rz(-1.057583) q[0];
rz(-2.019577) q[2];
sx q[2];
rz(-1.1626662) q[2];
sx q[2];
rz(-0.47779065) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.9895091) q[1];
sx q[1];
rz(-2.9674917) q[1];
sx q[1];
rz(2.4071481) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4336171) q[3];
sx q[3];
rz(-1.3841419) q[3];
sx q[3];
rz(-2.942323) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.8282738) q[2];
sx q[2];
rz(-0.88130772) q[2];
sx q[2];
rz(-0.71448294) q[2];
rz(2.1099527) q[3];
sx q[3];
rz(-2.1561626) q[3];
sx q[3];
rz(2.6825582) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3752876) q[0];
sx q[0];
rz(-0.77852386) q[0];
sx q[0];
rz(-0.2226204) q[0];
rz(0.34046945) q[1];
sx q[1];
rz(-1.4687931) q[1];
sx q[1];
rz(-2.8276843) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4231611) q[0];
sx q[0];
rz(-0.88215798) q[0];
sx q[0];
rz(-2.5325943) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3000175) q[2];
sx q[2];
rz(-1.944659) q[2];
sx q[2];
rz(-0.29619994) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.43526486) q[1];
sx q[1];
rz(-1.2641331) q[1];
sx q[1];
rz(-2.4198662) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7438358) q[3];
sx q[3];
rz(-0.53103775) q[3];
sx q[3];
rz(0.49285313) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.24806222) q[2];
sx q[2];
rz(-1.5241104) q[2];
sx q[2];
rz(2.1171872) q[2];
rz(-1.3192568) q[3];
sx q[3];
rz(-0.28194591) q[3];
sx q[3];
rz(2.8446741) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1612741) q[0];
sx q[0];
rz(-1.7151105) q[0];
sx q[0];
rz(-1.0262161) q[0];
rz(-0.16464344) q[1];
sx q[1];
rz(-2.7963729) q[1];
sx q[1];
rz(2.3688597) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.39557657) q[0];
sx q[0];
rz(-0.74953929) q[0];
sx q[0];
rz(-0.54732449) q[0];
rz(-1.5748789) q[2];
sx q[2];
rz(-1.2130623) q[2];
sx q[2];
rz(-0.28752354) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.1387921) q[1];
sx q[1];
rz(-1.7874663) q[1];
sx q[1];
rz(2.2598221) q[1];
x q[2];
rz(-1.5192612) q[3];
sx q[3];
rz(-1.5750575) q[3];
sx q[3];
rz(-0.86458998) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.905895) q[2];
sx q[2];
rz(-2.5696281) q[2];
sx q[2];
rz(-0.31615654) q[2];
rz(0.21137992) q[3];
sx q[3];
rz(-1.4485161) q[3];
sx q[3];
rz(-2.0785418) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.137977) q[0];
sx q[0];
rz(-1.5253541) q[0];
sx q[0];
rz(0.77950087) q[0];
rz(1.7100854) q[1];
sx q[1];
rz(-1.6110438) q[1];
sx q[1];
rz(0.12758189) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4460642) q[0];
sx q[0];
rz(-1.4565153) q[0];
sx q[0];
rz(0.034868783) q[0];
x q[1];
rz(-0.95037912) q[2];
sx q[2];
rz(-2.1135672) q[2];
sx q[2];
rz(-1.448878) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.7292916) q[1];
sx q[1];
rz(-2.6874098) q[1];
sx q[1];
rz(-0.95328625) q[1];
rz(-pi) q[2];
rz(-2.3977816) q[3];
sx q[3];
rz(-2.3530745) q[3];
sx q[3];
rz(0.34495993) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.7369507) q[2];
sx q[2];
rz(-1.5664132) q[2];
sx q[2];
rz(2.289782) q[2];
rz(2.8972054) q[3];
sx q[3];
rz(-0.71072018) q[3];
sx q[3];
rz(-1.3062668) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7459625) q[0];
sx q[0];
rz(-1.4816477) q[0];
sx q[0];
rz(0.10598824) q[0];
rz(-1.4251739) q[1];
sx q[1];
rz(-1.9916078) q[1];
sx q[1];
rz(1.0894159) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.16817753) q[0];
sx q[0];
rz(-1.3333048) q[0];
sx q[0];
rz(-0.071478162) q[0];
rz(-0.52417119) q[2];
sx q[2];
rz(-2.3830915) q[2];
sx q[2];
rz(2.3189409) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.32192595) q[1];
sx q[1];
rz(-1.2533256) q[1];
sx q[1];
rz(-0.13369932) q[1];
rz(-pi) q[2];
rz(0.15481101) q[3];
sx q[3];
rz(-2.1795978) q[3];
sx q[3];
rz(-1.4325292) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.0561698) q[2];
sx q[2];
rz(-1.9499754) q[2];
sx q[2];
rz(2.5540111) q[2];
rz(-1.9116631) q[3];
sx q[3];
rz(-2.7432224) q[3];
sx q[3];
rz(-2.5689094) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9375482) q[0];
sx q[0];
rz(-1.2296822) q[0];
sx q[0];
rz(0.44573927) q[0];
rz(-0.84272376) q[1];
sx q[1];
rz(-2.3363967) q[1];
sx q[1];
rz(-0.79949784) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.56109062) q[0];
sx q[0];
rz(-1.9724625) q[0];
sx q[0];
rz(-0.40720823) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4183774) q[2];
sx q[2];
rz(-2.2951381) q[2];
sx q[2];
rz(2.3220313) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-3.1285727) q[1];
sx q[1];
rz(-1.1718084) q[1];
sx q[1];
rz(-0.47069957) q[1];
rz(-pi) q[2];
rz(1.6715253) q[3];
sx q[3];
rz(-0.58222773) q[3];
sx q[3];
rz(-1.738172) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.2004956) q[2];
sx q[2];
rz(-2.7636187) q[2];
sx q[2];
rz(2.359158) q[2];
rz(-2.2109219) q[3];
sx q[3];
rz(-1.8698255) q[3];
sx q[3];
rz(-2.8455287) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.74872577) q[0];
sx q[0];
rz(-2.1935232) q[0];
sx q[0];
rz(-2.0178846) q[0];
rz(0.8404845) q[1];
sx q[1];
rz(-1.9978465) q[1];
sx q[1];
rz(-1.7488272) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26396093) q[0];
sx q[0];
rz(-3.0288806) q[0];
sx q[0];
rz(1.1219704) q[0];
x q[1];
rz(-2.8181067) q[2];
sx q[2];
rz(-1.3703823) q[2];
sx q[2];
rz(0.85497626) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.1906793) q[1];
sx q[1];
rz(-1.2251405) q[1];
sx q[1];
rz(0.18195621) q[1];
rz(-3.0245549) q[3];
sx q[3];
rz(-2.5212332) q[3];
sx q[3];
rz(0.24347603) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.8380518) q[2];
sx q[2];
rz(-2.3736062) q[2];
sx q[2];
rz(1.6221907) q[2];
rz(1.0649902) q[3];
sx q[3];
rz(-1.0737373) q[3];
sx q[3];
rz(0.85552335) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72724718) q[0];
sx q[0];
rz(-2.1942744) q[0];
sx q[0];
rz(1.1668246) q[0];
rz(-0.30544063) q[1];
sx q[1];
rz(-1.069953) q[1];
sx q[1];
rz(-2.6397612) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0286642) q[0];
sx q[0];
rz(-1.6643066) q[0];
sx q[0];
rz(-2.9458955) q[0];
x q[1];
rz(0.73955543) q[2];
sx q[2];
rz(-0.90074343) q[2];
sx q[2];
rz(0.53436992) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.8055685) q[1];
sx q[1];
rz(-0.37266392) q[1];
sx q[1];
rz(1.7737081) q[1];
rz(-pi) q[2];
x q[2];
rz(0.54399854) q[3];
sx q[3];
rz(-2.0680313) q[3];
sx q[3];
rz(1.0656644) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.66494232) q[2];
sx q[2];
rz(-2.3822337) q[2];
sx q[2];
rz(0.13258983) q[2];
rz(-2.8442123) q[3];
sx q[3];
rz(-1.5404276) q[3];
sx q[3];
rz(-0.69445777) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7731758) q[0];
sx q[0];
rz(-1.6650124) q[0];
sx q[0];
rz(-1.6126527) q[0];
rz(-1.7932786) q[1];
sx q[1];
rz(-1.7650083) q[1];
sx q[1];
rz(2.7682159) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2397176) q[0];
sx q[0];
rz(-1.7267701) q[0];
sx q[0];
rz(-1.2084991) q[0];
x q[1];
rz(1.2676722) q[2];
sx q[2];
rz(-2.0344007) q[2];
sx q[2];
rz(1.3410717) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.11168682) q[1];
sx q[1];
rz(-0.77757793) q[1];
sx q[1];
rz(-1.0295246) q[1];
rz(-pi) q[2];
rz(0.67383234) q[3];
sx q[3];
rz(-2.8008409) q[3];
sx q[3];
rz(1.0125404) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.0384963) q[2];
sx q[2];
rz(-0.080048397) q[2];
sx q[2];
rz(2.6888964) q[2];
rz(0.90940851) q[3];
sx q[3];
rz(-1.012864) q[3];
sx q[3];
rz(-0.67210853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7047003) q[0];
sx q[0];
rz(-1.1672651) q[0];
sx q[0];
rz(1.1080909) q[0];
rz(-2.3469901) q[1];
sx q[1];
rz(-1.4789076) q[1];
sx q[1];
rz(2.2750003) q[1];
rz(0.022877175) q[2];
sx q[2];
rz(-2.2616017) q[2];
sx q[2];
rz(0.2837413) q[2];
rz(-2.361479) q[3];
sx q[3];
rz(-1.2419392) q[3];
sx q[3];
rz(1.8942647) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
