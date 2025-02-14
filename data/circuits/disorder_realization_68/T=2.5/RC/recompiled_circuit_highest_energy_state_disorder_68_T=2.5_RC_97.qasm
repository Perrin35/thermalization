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
rz(-1.4425059) q[0];
sx q[0];
rz(-1.2966172) q[0];
sx q[0];
rz(1.698864) q[0];
rz(2.7891085) q[1];
sx q[1];
rz(3.6678996) q[1];
sx q[1];
rz(7.6511135) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5697081) q[0];
sx q[0];
rz(-1.3131327) q[0];
sx q[0];
rz(-0.37057693) q[0];
rz(1.5330731) q[2];
sx q[2];
rz(-1.9262656) q[2];
sx q[2];
rz(-1.0743574) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.753073) q[1];
sx q[1];
rz(-2.6021829) q[1];
sx q[1];
rz(-2.8864157) q[1];
rz(2.3867704) q[3];
sx q[3];
rz(-2.0275381) q[3];
sx q[3];
rz(2.7025102) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0736531) q[2];
sx q[2];
rz(-1.0509793) q[2];
sx q[2];
rz(1.2037207) q[2];
rz(0.97483557) q[3];
sx q[3];
rz(-2.1741512) q[3];
sx q[3];
rz(-1.2459285) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.75342733) q[0];
sx q[0];
rz(-2.0797256) q[0];
sx q[0];
rz(0.88428307) q[0];
rz(2.4321804) q[1];
sx q[1];
rz(-1.8953036) q[1];
sx q[1];
rz(0.90219227) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1195148) q[0];
sx q[0];
rz(-1.365322) q[0];
sx q[0];
rz(-0.46066649) q[0];
x q[1];
rz(-1.7184751) q[2];
sx q[2];
rz(-2.5294249) q[2];
sx q[2];
rz(-1.8107121) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.883256) q[1];
sx q[1];
rz(-2.751103) q[1];
sx q[1];
rz(-2.4526333) q[1];
rz(-1.7092) q[3];
sx q[3];
rz(-1.4908199) q[3];
sx q[3];
rz(3.0113335) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.117729) q[2];
sx q[2];
rz(-2.8459097) q[2];
sx q[2];
rz(-1.4884865) q[2];
rz(-1.9272517) q[3];
sx q[3];
rz(-1.4597273) q[3];
sx q[3];
rz(1.2795718) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0073256) q[0];
sx q[0];
rz(-1.3995582) q[0];
sx q[0];
rz(-0.16635995) q[0];
rz(2.4436489) q[1];
sx q[1];
rz(-2.3614387) q[1];
sx q[1];
rz(1.4770329) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1082527) q[0];
sx q[0];
rz(-2.0128002) q[0];
sx q[0];
rz(3.1382794) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8700902) q[2];
sx q[2];
rz(-1.2997441) q[2];
sx q[2];
rz(1.0752016) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.055226) q[1];
sx q[1];
rz(-0.7387195) q[1];
sx q[1];
rz(-2.2768873) q[1];
rz(-1.5904347) q[3];
sx q[3];
rz(-2.0691039) q[3];
sx q[3];
rz(-2.7036959) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.6537689) q[2];
sx q[2];
rz(-1.1710125) q[2];
sx q[2];
rz(2.5269395) q[2];
rz(1.4608308) q[3];
sx q[3];
rz(-1.5771644) q[3];
sx q[3];
rz(-3.0998668) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.6093269) q[0];
sx q[0];
rz(-1.3294687) q[0];
sx q[0];
rz(-1.5186658) q[0];
rz(-2.3329349) q[1];
sx q[1];
rz(-1.2963908) q[1];
sx q[1];
rz(-0.048695806) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8762465) q[0];
sx q[0];
rz(-0.9101724) q[0];
sx q[0];
rz(1.9991849) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.28566177) q[2];
sx q[2];
rz(-0.39644074) q[2];
sx q[2];
rz(0.91914058) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.073507) q[1];
sx q[1];
rz(-2.0754015) q[1];
sx q[1];
rz(1.0166974) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1761769) q[3];
sx q[3];
rz(-1.2413238) q[3];
sx q[3];
rz(-0.62240619) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.101717) q[2];
sx q[2];
rz(-2.3287435) q[2];
sx q[2];
rz(-1.2561049) q[2];
rz(3.0806372) q[3];
sx q[3];
rz(-0.90164369) q[3];
sx q[3];
rz(0.92756334) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1561279) q[0];
sx q[0];
rz(-2.0138854) q[0];
sx q[0];
rz(-0.10966478) q[0];
rz(-2.1127286) q[1];
sx q[1];
rz(-1.8548465) q[1];
sx q[1];
rz(-1.5493772) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3774535) q[0];
sx q[0];
rz(-1.8528588) q[0];
sx q[0];
rz(-0.46746032) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.36269168) q[2];
sx q[2];
rz(-2.0433934) q[2];
sx q[2];
rz(0.23841274) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.2031496) q[1];
sx q[1];
rz(-2.0324273) q[1];
sx q[1];
rz(2.0600956) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3098833) q[3];
sx q[3];
rz(-3.0906537) q[3];
sx q[3];
rz(0.65492958) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.5900383) q[2];
sx q[2];
rz(-2.8503214) q[2];
sx q[2];
rz(0.48198286) q[2];
rz(0.41989741) q[3];
sx q[3];
rz(-2.0764543) q[3];
sx q[3];
rz(-0.70986748) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9919306) q[0];
sx q[0];
rz(-0.88352942) q[0];
sx q[0];
rz(-0.71257198) q[0];
rz(0.86992162) q[1];
sx q[1];
rz(-0.37418071) q[1];
sx q[1];
rz(-3.1178927) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7837304) q[0];
sx q[0];
rz(-2.0867769) q[0];
sx q[0];
rz(-1.7586772) q[0];
rz(-pi) q[1];
rz(2.4871102) q[2];
sx q[2];
rz(-2.3238782) q[2];
sx q[2];
rz(-0.47374642) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.5620251) q[1];
sx q[1];
rz(-1.4050086) q[1];
sx q[1];
rz(0.001494249) q[1];
rz(2.6610801) q[3];
sx q[3];
rz(-1.7578205) q[3];
sx q[3];
rz(2.9399237) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.80453834) q[2];
sx q[2];
rz(-0.64547515) q[2];
sx q[2];
rz(-2.0013334) q[2];
rz(-1.987847) q[3];
sx q[3];
rz(-1.9269582) q[3];
sx q[3];
rz(1.8203189) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.425728) q[0];
sx q[0];
rz(-1.3874929) q[0];
sx q[0];
rz(-2.7737889) q[0];
rz(1.6416719) q[1];
sx q[1];
rz(-0.92253128) q[1];
sx q[1];
rz(2.0466764) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.33312265) q[0];
sx q[0];
rz(-0.54352353) q[0];
sx q[0];
rz(2.031206) q[0];
rz(-pi) q[1];
rz(1.7000654) q[2];
sx q[2];
rz(-2.1174927) q[2];
sx q[2];
rz(1.136029) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.3125632) q[1];
sx q[1];
rz(-0.2364279) q[1];
sx q[1];
rz(-3.0655686) q[1];
rz(-pi) q[2];
x q[2];
rz(3.1203116) q[3];
sx q[3];
rz(-1.9742734) q[3];
sx q[3];
rz(1.9627987) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.0423923) q[2];
sx q[2];
rz(-1.3513869) q[2];
sx q[2];
rz(3.0858827) q[2];
rz(-2.4646711) q[3];
sx q[3];
rz(-2.5261295) q[3];
sx q[3];
rz(2.2630579) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1894576) q[0];
sx q[0];
rz(-2.5791233) q[0];
sx q[0];
rz(2.1761555) q[0];
rz(1.2646487) q[1];
sx q[1];
rz(-0.74600428) q[1];
sx q[1];
rz(-2.8146578) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.0069939216) q[0];
sx q[0];
rz(-0.36834684) q[0];
sx q[0];
rz(-1.2587496) q[0];
rz(-1.5599653) q[2];
sx q[2];
rz(-2.0167266) q[2];
sx q[2];
rz(-1.8919093) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.2911243) q[1];
sx q[1];
rz(-1.7479436) q[1];
sx q[1];
rz(1.3590823) q[1];
rz(-1.502874) q[3];
sx q[3];
rz(-1.1227896) q[3];
sx q[3];
rz(-0.94326708) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.38291976) q[2];
sx q[2];
rz(-0.70432538) q[2];
sx q[2];
rz(-2.3455589) q[2];
rz(0.087513611) q[3];
sx q[3];
rz(-2.5329866) q[3];
sx q[3];
rz(2.2039738) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7079119) q[0];
sx q[0];
rz(-1.3734564) q[0];
sx q[0];
rz(0.32980907) q[0];
rz(-3.0682796) q[1];
sx q[1];
rz(-0.70711702) q[1];
sx q[1];
rz(2.0475533) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1739376) q[0];
sx q[0];
rz(-2.6693845) q[0];
sx q[0];
rz(-3.0595991) q[0];
rz(-1.173063) q[2];
sx q[2];
rz(-1.6327883) q[2];
sx q[2];
rz(-0.50794377) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.0713404) q[1];
sx q[1];
rz(-2.1161181) q[1];
sx q[1];
rz(1.7252183) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.002203) q[3];
sx q[3];
rz(-1.6069222) q[3];
sx q[3];
rz(-2.2808035) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.0988934) q[2];
sx q[2];
rz(-2.3697479) q[2];
sx q[2];
rz(1.137255) q[2];
rz(-2.5134261) q[3];
sx q[3];
rz(-2.7558432) q[3];
sx q[3];
rz(1.0073193) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6657418) q[0];
sx q[0];
rz(-0.51097149) q[0];
sx q[0];
rz(-1.092859) q[0];
rz(1.8765556) q[1];
sx q[1];
rz(-1.8362412) q[1];
sx q[1];
rz(-2.1078033) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2726934) q[0];
sx q[0];
rz(-2.0420235) q[0];
sx q[0];
rz(-0.97288218) q[0];
rz(0.7598147) q[2];
sx q[2];
rz(-1.7489479) q[2];
sx q[2];
rz(2.9896328) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.8333697) q[1];
sx q[1];
rz(-0.78082456) q[1];
sx q[1];
rz(-2.5379347) q[1];
x q[2];
rz(-0.035461144) q[3];
sx q[3];
rz(-1.7093813) q[3];
sx q[3];
rz(-1.1690804) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.25896245) q[2];
sx q[2];
rz(-2.4179103) q[2];
sx q[2];
rz(0.21931973) q[2];
rz(-2.3389881) q[3];
sx q[3];
rz(-1.31253) q[3];
sx q[3];
rz(2.820211) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6312859) q[0];
sx q[0];
rz(-1.8230556) q[0];
sx q[0];
rz(-2.0429116) q[0];
rz(0.17768271) q[1];
sx q[1];
rz(-2.1251353) q[1];
sx q[1];
rz(2.9716117) q[1];
rz(-0.86831696) q[2];
sx q[2];
rz(-1.5219896) q[2];
sx q[2];
rz(0.513214) q[2];
rz(-3.0232676) q[3];
sx q[3];
rz(-1.3814303) q[3];
sx q[3];
rz(-0.43683972) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
