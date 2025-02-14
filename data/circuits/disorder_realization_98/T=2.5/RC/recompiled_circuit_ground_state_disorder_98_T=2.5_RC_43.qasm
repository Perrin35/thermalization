OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.15220517) q[0];
sx q[0];
rz(2.9289991) q[0];
sx q[0];
rz(8.6906035) q[0];
rz(-1.9250159) q[1];
sx q[1];
rz(2.7956378) q[1];
sx q[1];
rz(12.591127) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9182565) q[0];
sx q[0];
rz(-1.6331216) q[0];
sx q[0];
rz(1.5274836) q[0];
rz(0.11344929) q[2];
sx q[2];
rz(-1.3384124) q[2];
sx q[2];
rz(-2.8040545) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.5879257) q[1];
sx q[1];
rz(-0.30372639) q[1];
sx q[1];
rz(0.17596735) q[1];
rz(1.3583899) q[3];
sx q[3];
rz(-1.6254566) q[3];
sx q[3];
rz(-0.89059356) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.88432246) q[2];
sx q[2];
rz(-1.4929644) q[2];
sx q[2];
rz(2.8327668) q[2];
rz(2.3158) q[3];
sx q[3];
rz(-2.1335996) q[3];
sx q[3];
rz(-2.3776313) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9708213) q[0];
sx q[0];
rz(-1.229137) q[0];
sx q[0];
rz(2.1139076) q[0];
rz(2.3254501) q[1];
sx q[1];
rz(-2.0232537) q[1];
sx q[1];
rz(-0.81248409) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.9922794) q[0];
sx q[0];
rz(-1.9198787) q[0];
sx q[0];
rz(-2.9724253) q[0];
rz(-pi) q[1];
rz(0.036244321) q[2];
sx q[2];
rz(-0.37359056) q[2];
sx q[2];
rz(-0.61818733) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.634944) q[1];
sx q[1];
rz(-2.6104984) q[1];
sx q[1];
rz(1.2960766) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2151617) q[3];
sx q[3];
rz(-1.8994678) q[3];
sx q[3];
rz(1.6021014) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.5673148) q[2];
sx q[2];
rz(-1.4894166) q[2];
sx q[2];
rz(0.91892773) q[2];
rz(2.610176) q[3];
sx q[3];
rz(-1.0454949) q[3];
sx q[3];
rz(3.1004068) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4384005) q[0];
sx q[0];
rz(-1.0696609) q[0];
sx q[0];
rz(1.137314) q[0];
rz(-1.7510341) q[1];
sx q[1];
rz(-2.5847692) q[1];
sx q[1];
rz(0.73296076) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3175388) q[0];
sx q[0];
rz(-1.8630872) q[0];
sx q[0];
rz(2.7089416) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.175019) q[2];
sx q[2];
rz(-1.086965) q[2];
sx q[2];
rz(-1.0072198) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.7480173) q[1];
sx q[1];
rz(-2.0055341) q[1];
sx q[1];
rz(-1.0487021) q[1];
rz(2.0195229) q[3];
sx q[3];
rz(-1.5767234) q[3];
sx q[3];
rz(-2.0172265) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.0040032337) q[2];
sx q[2];
rz(-1.4121476) q[2];
sx q[2];
rz(-2.1027193) q[2];
rz(1.002257) q[3];
sx q[3];
rz(-1.301731) q[3];
sx q[3];
rz(-0.29016289) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9349979) q[0];
sx q[0];
rz(-2.0552141) q[0];
sx q[0];
rz(-2.2344053) q[0];
rz(2.9838003) q[1];
sx q[1];
rz(-1.0148427) q[1];
sx q[1];
rz(1.2812322) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8195551) q[0];
sx q[0];
rz(-2.0896308) q[0];
sx q[0];
rz(0.71625336) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0220535) q[2];
sx q[2];
rz(-2.3733778) q[2];
sx q[2];
rz(2.8295585) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.80824796) q[1];
sx q[1];
rz(-0.95511694) q[1];
sx q[1];
rz(-2.0991825) q[1];
rz(-pi) q[2];
x q[2];
rz(0.90003711) q[3];
sx q[3];
rz(-2.0931912) q[3];
sx q[3];
rz(-0.99985048) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.8944051) q[2];
sx q[2];
rz(-2.2053568) q[2];
sx q[2];
rz(1.2897162) q[2];
rz(-1.7973409) q[3];
sx q[3];
rz(-1.684609) q[3];
sx q[3];
rz(0.70014203) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46071389) q[0];
sx q[0];
rz(-2.6432156) q[0];
sx q[0];
rz(1.0182678) q[0];
rz(2.0385888) q[1];
sx q[1];
rz(-2.4296727) q[1];
sx q[1];
rz(1.3404554) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0606275) q[0];
sx q[0];
rz(-2.5239391) q[0];
sx q[0];
rz(-0.65184848) q[0];
rz(-pi) q[1];
rz(-1.6685772) q[2];
sx q[2];
rz(-1.9488397) q[2];
sx q[2];
rz(-2.5387272) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.75234883) q[1];
sx q[1];
rz(-1.6977786) q[1];
sx q[1];
rz(-1.4402706) q[1];
rz(-1.6667073) q[3];
sx q[3];
rz(-0.54236327) q[3];
sx q[3];
rz(2.9090442) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.5420142) q[2];
sx q[2];
rz(-0.88853637) q[2];
sx q[2];
rz(1.0020024) q[2];
rz(-2.4501948) q[3];
sx q[3];
rz(-1.1804429) q[3];
sx q[3];
rz(1.6208167) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6028676) q[0];
sx q[0];
rz(-2.0724917) q[0];
sx q[0];
rz(-2.5163203) q[0];
rz(-1.193115) q[1];
sx q[1];
rz(-2.4517877) q[1];
sx q[1];
rz(2.5742721) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0085908169) q[0];
sx q[0];
rz(-0.80256217) q[0];
sx q[0];
rz(2.8961097) q[0];
rz(-pi) q[1];
rz(0.0074947798) q[2];
sx q[2];
rz(-0.82641331) q[2];
sx q[2];
rz(1.9454582) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.3750227) q[1];
sx q[1];
rz(-0.61900292) q[1];
sx q[1];
rz(1.5260215) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9653699) q[3];
sx q[3];
rz(-1.3941951) q[3];
sx q[3];
rz(-3.0295102) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.304473) q[2];
sx q[2];
rz(-1.8581055) q[2];
sx q[2];
rz(-1.5817969) q[2];
rz(2.5290153) q[3];
sx q[3];
rz(-1.9809096) q[3];
sx q[3];
rz(-0.15577236) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0953858) q[0];
sx q[0];
rz(-1.950773) q[0];
sx q[0];
rz(-2.0147391) q[0];
rz(-1.2696179) q[1];
sx q[1];
rz(-1.1879299) q[1];
sx q[1];
rz(0.52106214) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.71741912) q[0];
sx q[0];
rz(-1.9195286) q[0];
sx q[0];
rz(-0.86535378) q[0];
rz(1.4748739) q[2];
sx q[2];
rz(-2.811199) q[2];
sx q[2];
rz(-1.0164193) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.4849629) q[1];
sx q[1];
rz(-0.99492517) q[1];
sx q[1];
rz(-2.5099975) q[1];
rz(-pi) q[2];
rz(0.54102202) q[3];
sx q[3];
rz(-1.8086026) q[3];
sx q[3];
rz(-2.1062807) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.0236987) q[2];
sx q[2];
rz(-1.7774899) q[2];
sx q[2];
rz(1.914631) q[2];
rz(-1.6795233) q[3];
sx q[3];
rz(-1.5173802) q[3];
sx q[3];
rz(2.7000361) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0284477) q[0];
sx q[0];
rz(-2.1118836) q[0];
sx q[0];
rz(2.5944769) q[0];
rz(3.0534577) q[1];
sx q[1];
rz(-1.793975) q[1];
sx q[1];
rz(-2.7190582) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8478017) q[0];
sx q[0];
rz(-1.5489086) q[0];
sx q[0];
rz(1.8518015) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9401378) q[2];
sx q[2];
rz(-1.0779194) q[2];
sx q[2];
rz(-2.3702459) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.40228981) q[1];
sx q[1];
rz(-1.5147494) q[1];
sx q[1];
rz(1.5016374) q[1];
rz(1.3874153) q[3];
sx q[3];
rz(-1.7088582) q[3];
sx q[3];
rz(-2.3473397) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.2253458) q[2];
sx q[2];
rz(-1.6371181) q[2];
sx q[2];
rz(-2.8543191) q[2];
rz(2.5987127) q[3];
sx q[3];
rz(-0.77016872) q[3];
sx q[3];
rz(-0.058102593) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3995689) q[0];
sx q[0];
rz(-2.4574807) q[0];
sx q[0];
rz(2.3642819) q[0];
rz(2.2151392) q[1];
sx q[1];
rz(-1.9994206) q[1];
sx q[1];
rz(-1.3949589) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.762601) q[0];
sx q[0];
rz(-1.1830353) q[0];
sx q[0];
rz(0.34404018) q[0];
x q[1];
rz(2.5424273) q[2];
sx q[2];
rz(-1.3627421) q[2];
sx q[2];
rz(-2.4924202) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.7662188) q[1];
sx q[1];
rz(-1.2412984) q[1];
sx q[1];
rz(1.0887515) q[1];
rz(-pi) q[2];
rz(-2.5245669) q[3];
sx q[3];
rz(-2.113638) q[3];
sx q[3];
rz(1.9844685) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.7159783) q[2];
sx q[2];
rz(-1.7691111) q[2];
sx q[2];
rz(2.051029) q[2];
rz(0.96450949) q[3];
sx q[3];
rz(-1.0114074) q[3];
sx q[3];
rz(-1.9264268) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5665117) q[0];
sx q[0];
rz(-0.33841857) q[0];
sx q[0];
rz(-1.6315208) q[0];
rz(3.1072726) q[1];
sx q[1];
rz(-1.3395373) q[1];
sx q[1];
rz(0.43201772) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5199909) q[0];
sx q[0];
rz(-1.5120602) q[0];
sx q[0];
rz(-2.1099173) q[0];
rz(-pi) q[1];
rz(1.2197184) q[2];
sx q[2];
rz(-1.4939927) q[2];
sx q[2];
rz(1.6258282) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.2434531) q[1];
sx q[1];
rz(-1.6695078) q[1];
sx q[1];
rz(2.9961186) q[1];
rz(-pi) q[2];
x q[2];
rz(0.26507399) q[3];
sx q[3];
rz(-2.1360096) q[3];
sx q[3];
rz(0.044818002) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.67046514) q[2];
sx q[2];
rz(-1.1521143) q[2];
sx q[2];
rz(2.067789) q[2];
rz(-2.9527169) q[3];
sx q[3];
rz(-0.60868588) q[3];
sx q[3];
rz(2.7591738) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.043924532) q[0];
sx q[0];
rz(-0.32079874) q[0];
sx q[0];
rz(-0.70377845) q[0];
rz(-1.3532928) q[1];
sx q[1];
rz(-2.4663993) q[1];
sx q[1];
rz(-0.83723062) q[1];
rz(-0.7424217) q[2];
sx q[2];
rz(-1.0992285) q[2];
sx q[2];
rz(0.64115094) q[2];
rz(-1.6394284) q[3];
sx q[3];
rz(-1.1123688) q[3];
sx q[3];
rz(0.49751626) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
