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
rz(-2.6211205) q[0];
sx q[0];
rz(-1.0360798) q[0];
sx q[0];
rz(0.64089027) q[0];
rz(-0.23735292) q[1];
sx q[1];
rz(-2.4999455) q[1];
sx q[1];
rz(3.0566888) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2332524) q[0];
sx q[0];
rz(-1.4734992) q[0];
sx q[0];
rz(-2.2761274) q[0];
x q[1];
rz(0.48975421) q[2];
sx q[2];
rz(-0.066250525) q[2];
sx q[2];
rz(-0.79023933) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.3255439) q[1];
sx q[1];
rz(-2.8286165) q[1];
sx q[1];
rz(-2.7626687) q[1];
x q[2];
rz(-2.747211) q[3];
sx q[3];
rz(-2.1254) q[3];
sx q[3];
rz(1.7930052) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.889633) q[2];
sx q[2];
rz(-1.9388988) q[2];
sx q[2];
rz(1.5084051) q[2];
rz(-2.3097322) q[3];
sx q[3];
rz(-2.8155477) q[3];
sx q[3];
rz(1.3445725) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2107596) q[0];
sx q[0];
rz(-0.63527125) q[0];
sx q[0];
rz(-0.26500901) q[0];
rz(2.3459072) q[1];
sx q[1];
rz(-0.34514752) q[1];
sx q[1];
rz(-1.4211242) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4583295) q[0];
sx q[0];
rz(-2.4340672) q[0];
sx q[0];
rz(-1.1573155) q[0];
x q[1];
rz(2.1876343) q[2];
sx q[2];
rz(-1.9618615) q[2];
sx q[2];
rz(0.80634889) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.15884033) q[1];
sx q[1];
rz(-1.5271865) q[1];
sx q[1];
rz(-1.3189454) q[1];
x q[2];
rz(1.7461997) q[3];
sx q[3];
rz(-2.9713745) q[3];
sx q[3];
rz(2.0147918) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.77218324) q[2];
sx q[2];
rz(-2.1464244) q[2];
sx q[2];
rz(0.65917242) q[2];
rz(-2.2199421) q[3];
sx q[3];
rz(-2.0989336) q[3];
sx q[3];
rz(1.5590182) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7121861) q[0];
sx q[0];
rz(-0.61539188) q[0];
sx q[0];
rz(1.2834826) q[0];
rz(-2.7146924) q[1];
sx q[1];
rz(-0.59395298) q[1];
sx q[1];
rz(1.1675507) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99770633) q[0];
sx q[0];
rz(-0.50181118) q[0];
sx q[0];
rz(-2.3104179) q[0];
rz(-pi) q[1];
rz(-1.0658355) q[2];
sx q[2];
rz(-1.0965818) q[2];
sx q[2];
rz(1.1467288) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.1192273) q[1];
sx q[1];
rz(-1.4767547) q[1];
sx q[1];
rz(-0.29081197) q[1];
rz(-pi) q[2];
rz(0.57691388) q[3];
sx q[3];
rz(-2.8261999) q[3];
sx q[3];
rz(-1.2030935) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.9350962) q[2];
sx q[2];
rz(-1.4859716) q[2];
sx q[2];
rz(-0.11779724) q[2];
rz(1.3055034) q[3];
sx q[3];
rz(-2.2820303) q[3];
sx q[3];
rz(-0.99087244) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.30462581) q[0];
sx q[0];
rz(-1.0517629) q[0];
sx q[0];
rz(-3.1090609) q[0];
rz(-2.1578728) q[1];
sx q[1];
rz(-0.81147057) q[1];
sx q[1];
rz(-0.54289114) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6257831) q[0];
sx q[0];
rz(-2.2383732) q[0];
sx q[0];
rz(-1.2207009) q[0];
x q[1];
rz(3.0383598) q[2];
sx q[2];
rz(-1.8816173) q[2];
sx q[2];
rz(-0.23676591) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.11527744) q[1];
sx q[1];
rz(-1.8247073) q[1];
sx q[1];
rz(2.3350969) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.832295) q[3];
sx q[3];
rz(-1.4595539) q[3];
sx q[3];
rz(2.5865366) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.5618318) q[2];
sx q[2];
rz(-1.1299955) q[2];
sx q[2];
rz(-2.180991) q[2];
rz(-2.1660755) q[3];
sx q[3];
rz(-0.87149182) q[3];
sx q[3];
rz(-0.52198854) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8303216) q[0];
sx q[0];
rz(-2.4595342) q[0];
sx q[0];
rz(-0.3325381) q[0];
rz(-0.63811103) q[1];
sx q[1];
rz(-0.6256012) q[1];
sx q[1];
rz(0.45613751) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.79955937) q[0];
sx q[0];
rz(-1.4385975) q[0];
sx q[0];
rz(-3.0782736) q[0];
rz(-pi) q[1];
rz(0.47686968) q[2];
sx q[2];
rz(-2.2945291) q[2];
sx q[2];
rz(-0.70666955) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.8673381) q[1];
sx q[1];
rz(-0.98877866) q[1];
sx q[1];
rz(-0.10669218) q[1];
rz(-pi) q[2];
x q[2];
rz(0.92335574) q[3];
sx q[3];
rz(-1.8643987) q[3];
sx q[3];
rz(1.3577611) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.0629603) q[2];
sx q[2];
rz(-0.76138622) q[2];
sx q[2];
rz(2.6893993) q[2];
rz(0.26442987) q[3];
sx q[3];
rz(-2.9954438) q[3];
sx q[3];
rz(-1.2368081) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0753157) q[0];
sx q[0];
rz(-0.83106581) q[0];
sx q[0];
rz(-2.4237295) q[0];
rz(1.588899) q[1];
sx q[1];
rz(-1.540686) q[1];
sx q[1];
rz(0.20725651) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8533181) q[0];
sx q[0];
rz(-2.69876) q[0];
sx q[0];
rz(0.063837039) q[0];
rz(-pi) q[1];
rz(3.0704801) q[2];
sx q[2];
rz(-1.9926535) q[2];
sx q[2];
rz(2.6049169) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.9178324) q[1];
sx q[1];
rz(-2.4759935) q[1];
sx q[1];
rz(0.83719801) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8317811) q[3];
sx q[3];
rz(-1.4647604) q[3];
sx q[3];
rz(0.98766726) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.76102105) q[2];
sx q[2];
rz(-2.0863159) q[2];
sx q[2];
rz(-2.5804248) q[2];
rz(1.0593972) q[3];
sx q[3];
rz(-2.0249517) q[3];
sx q[3];
rz(-1.4580956) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.92784944) q[0];
sx q[0];
rz(-0.97074592) q[0];
sx q[0];
rz(2.2947776) q[0];
rz(-0.98833409) q[1];
sx q[1];
rz(-1.3761995) q[1];
sx q[1];
rz(0.83289897) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7306379) q[0];
sx q[0];
rz(-2.9175076) q[0];
sx q[0];
rz(-0.48416324) q[0];
rz(-pi) q[1];
rz(1.49176) q[2];
sx q[2];
rz(-1.7124868) q[2];
sx q[2];
rz(0.92565432) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.33796577) q[1];
sx q[1];
rz(-2.4110893) q[1];
sx q[1];
rz(-2.1419163) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4571413) q[3];
sx q[3];
rz(-1.6196005) q[3];
sx q[3];
rz(1.5391853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.6758468) q[2];
sx q[2];
rz(-0.70285672) q[2];
sx q[2];
rz(0.77009002) q[2];
rz(0.13104023) q[3];
sx q[3];
rz(-1.4821056) q[3];
sx q[3];
rz(1.8946764) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.34762621) q[0];
sx q[0];
rz(-2.2063231) q[0];
sx q[0];
rz(-2.6499709) q[0];
rz(2.8288016) q[1];
sx q[1];
rz(-2.8520695) q[1];
sx q[1];
rz(3.0591931) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1954049) q[0];
sx q[0];
rz(-3.051149) q[0];
sx q[0];
rz(1.7121332) q[0];
rz(-2.5748207) q[2];
sx q[2];
rz(-2.1852036) q[2];
sx q[2];
rz(-2.4943309) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.8696599) q[1];
sx q[1];
rz(-0.23889506) q[1];
sx q[1];
rz(-0.82456692) q[1];
rz(-pi) q[2];
rz(1.8833722) q[3];
sx q[3];
rz(-1.4736611) q[3];
sx q[3];
rz(-2.7526906) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.4440492) q[2];
sx q[2];
rz(-1.2080344) q[2];
sx q[2];
rz(-3.1235798) q[2];
rz(2.2703914) q[3];
sx q[3];
rz(-0.33659354) q[3];
sx q[3];
rz(-0.66756162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.95053259) q[0];
sx q[0];
rz(-0.55452269) q[0];
sx q[0];
rz(-0.46448034) q[0];
rz(-0.49098268) q[1];
sx q[1];
rz(-1.4717439) q[1];
sx q[1];
rz(1.4674662) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.53290908) q[0];
sx q[0];
rz(-2.1957046) q[0];
sx q[0];
rz(-1.5708357) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2651132) q[2];
sx q[2];
rz(-2.0569498) q[2];
sx q[2];
rz(-0.4403688) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.0027568) q[1];
sx q[1];
rz(-1.8187539) q[1];
sx q[1];
rz(-3.0791705) q[1];
rz(-2.4201066) q[3];
sx q[3];
rz(-2.2389447) q[3];
sx q[3];
rz(0.20340445) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(3.1269425) q[2];
sx q[2];
rz(-1.4913968) q[2];
sx q[2];
rz(-1.3313741) q[2];
rz(2.2018382) q[3];
sx q[3];
rz(-2.0846114) q[3];
sx q[3];
rz(1.9542998) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2347655) q[0];
sx q[0];
rz(-2.6563788) q[0];
sx q[0];
rz(-2.1038798) q[0];
rz(-1.0543793) q[1];
sx q[1];
rz(-1.8260006) q[1];
sx q[1];
rz(1.0494999) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9572311) q[0];
sx q[0];
rz(-1.5210694) q[0];
sx q[0];
rz(0.00018361365) q[0];
rz(-pi) q[1];
x q[1];
rz(2.439211) q[2];
sx q[2];
rz(-2.4535013) q[2];
sx q[2];
rz(-1.976672) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.29925525) q[1];
sx q[1];
rz(-1.6016869) q[1];
sx q[1];
rz(0.96986356) q[1];
rz(-pi) q[2];
rz(-2.1621733) q[3];
sx q[3];
rz(-2.0775177) q[3];
sx q[3];
rz(3.0924071) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.99027571) q[2];
sx q[2];
rz(-1.8345202) q[2];
sx q[2];
rz(-1.6952197) q[2];
rz(-0.15898786) q[3];
sx q[3];
rz(-1.9828601) q[3];
sx q[3];
rz(-2.3915763) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1405519) q[0];
sx q[0];
rz(-2.1646071) q[0];
sx q[0];
rz(0.80143308) q[0];
rz(1.4871545) q[1];
sx q[1];
rz(-2.9947037) q[1];
sx q[1];
rz(-0.53302232) q[1];
rz(0.85114807) q[2];
sx q[2];
rz(-1.0590886) q[2];
sx q[2];
rz(1.5589489) q[2];
rz(-0.728353) q[3];
sx q[3];
rz(-2.3621758) q[3];
sx q[3];
rz(-2.4324535) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
