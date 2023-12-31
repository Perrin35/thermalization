OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.9053469) q[0];
sx q[0];
rz(-0.72609225) q[0];
sx q[0];
rz(-0.2015764) q[0];
rz(0.4959313) q[1];
sx q[1];
rz(-0.5402686) q[1];
sx q[1];
rz(-0.93710605) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0159338) q[0];
sx q[0];
rz(-2.4056245) q[0];
sx q[0];
rz(-0.65471162) q[0];
x q[1];
rz(2.4891698) q[2];
sx q[2];
rz(-2.4180275) q[2];
sx q[2];
rz(3.0384118) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.01552445) q[1];
sx q[1];
rz(-1.4926732) q[1];
sx q[1];
rz(-1.3874345) q[1];
rz(-3.044371) q[3];
sx q[3];
rz(-1.3029649) q[3];
sx q[3];
rz(2.3107446) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.15930882) q[2];
sx q[2];
rz(-2.8757877) q[2];
sx q[2];
rz(3.0554331) q[2];
rz(-2.384095) q[3];
sx q[3];
rz(-0.75452724) q[3];
sx q[3];
rz(1.0936201) q[3];
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
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5646097) q[0];
sx q[0];
rz(-1.6596721) q[0];
sx q[0];
rz(-1.8151059) q[0];
rz(1.2558698) q[1];
sx q[1];
rz(-1.565226) q[1];
sx q[1];
rz(-2.870141) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4685681) q[0];
sx q[0];
rz(-0.66267555) q[0];
sx q[0];
rz(-0.89232348) q[0];
rz(-pi) q[1];
rz(-2.6255252) q[2];
sx q[2];
rz(-1.8415673) q[2];
sx q[2];
rz(-2.5004636) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.506362) q[1];
sx q[1];
rz(-0.87345424) q[1];
sx q[1];
rz(0.49763775) q[1];
rz(-pi) q[2];
rz(0.80941697) q[3];
sx q[3];
rz(-1.3076233) q[3];
sx q[3];
rz(-2.0078878) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.0469971) q[2];
sx q[2];
rz(-1.3826028) q[2];
sx q[2];
rz(-2.5644152) q[2];
rz(-0.92352891) q[3];
sx q[3];
rz(-2.2369592) q[3];
sx q[3];
rz(-1.7318168) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8975824) q[0];
sx q[0];
rz(-1.0616466) q[0];
sx q[0];
rz(-2.999021) q[0];
rz(1.7890731) q[1];
sx q[1];
rz(-2.1058154) q[1];
sx q[1];
rz(-0.20908633) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7861048) q[0];
sx q[0];
rz(-2.3007563) q[0];
sx q[0];
rz(-0.79746042) q[0];
rz(2.694391) q[2];
sx q[2];
rz(-0.7913835) q[2];
sx q[2];
rz(-2.7176822) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.28029728) q[1];
sx q[1];
rz(-1.1251083) q[1];
sx q[1];
rz(-0.46511005) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.61061065) q[3];
sx q[3];
rz(-0.36452499) q[3];
sx q[3];
rz(-3.1162457) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.7769988) q[2];
sx q[2];
rz(-2.2475188) q[2];
sx q[2];
rz(-2.7374632) q[2];
rz(1.2858307) q[3];
sx q[3];
rz(-1.1288246) q[3];
sx q[3];
rz(-2.9742441) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4362815) q[0];
sx q[0];
rz(-3.0349338) q[0];
sx q[0];
rz(0.96281111) q[0];
rz(-0.46936938) q[1];
sx q[1];
rz(-2.5517187) q[1];
sx q[1];
rz(-0.00096360047) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9379291) q[0];
sx q[0];
rz(-2.2712049) q[0];
sx q[0];
rz(-0.80313375) q[0];
rz(-pi) q[1];
x q[1];
rz(3.1374627) q[2];
sx q[2];
rz(-0.12400707) q[2];
sx q[2];
rz(-2.4069402) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.5818134) q[1];
sx q[1];
rz(-1.0043136) q[1];
sx q[1];
rz(-2.986746) q[1];
rz(1.4812874) q[3];
sx q[3];
rz(-1.0032017) q[3];
sx q[3];
rz(-1.9247418) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.0659539) q[2];
sx q[2];
rz(-1.5116189) q[2];
sx q[2];
rz(-3.0299419) q[2];
rz(-0.81104898) q[3];
sx q[3];
rz(-2.6871197) q[3];
sx q[3];
rz(-3.1276935) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1902996) q[0];
sx q[0];
rz(-2.2409029) q[0];
sx q[0];
rz(-0.35650373) q[0];
rz(-0.50645343) q[1];
sx q[1];
rz(-1.3566147) q[1];
sx q[1];
rz(-0.26062632) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.168419) q[0];
sx q[0];
rz(-1.6126313) q[0];
sx q[0];
rz(2.9996458) q[0];
rz(-pi) q[1];
rz(1.8243276) q[2];
sx q[2];
rz(-1.8142482) q[2];
sx q[2];
rz(1.7340811) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.4958447) q[1];
sx q[1];
rz(-2.145088) q[1];
sx q[1];
rz(1.5917642) q[1];
x q[2];
rz(2.3652606) q[3];
sx q[3];
rz(-0.20271248) q[3];
sx q[3];
rz(-0.67938882) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.9115209) q[2];
sx q[2];
rz(-1.6610961) q[2];
sx q[2];
rz(-2.9122706) q[2];
rz(-0.54245943) q[3];
sx q[3];
rz(-0.31033236) q[3];
sx q[3];
rz(-2.1054161) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3689573) q[0];
sx q[0];
rz(-0.90894037) q[0];
sx q[0];
rz(1.4916346) q[0];
rz(-2.1024599) q[1];
sx q[1];
rz(-1.2960478) q[1];
sx q[1];
rz(-1.414149) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3154253) q[0];
sx q[0];
rz(-2.8624479) q[0];
sx q[0];
rz(0.17255737) q[0];
rz(0.51775198) q[2];
sx q[2];
rz(-1.3377829) q[2];
sx q[2];
rz(-1.053996) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.8687739) q[1];
sx q[1];
rz(-1.0827853) q[1];
sx q[1];
rz(-0.52656071) q[1];
rz(-pi) q[2];
rz(1.9023444) q[3];
sx q[3];
rz(-0.22833951) q[3];
sx q[3];
rz(2.6339298) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.1386537) q[2];
sx q[2];
rz(-2.5230375) q[2];
sx q[2];
rz(0.027739851) q[2];
rz(-2.6489143) q[3];
sx q[3];
rz(-1.6511107) q[3];
sx q[3];
rz(1.7475351) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7572927) q[0];
sx q[0];
rz(-1.5889656) q[0];
sx q[0];
rz(0.12167715) q[0];
rz(1.9901468) q[1];
sx q[1];
rz(-2.6897488) q[1];
sx q[1];
rz(-0.40245232) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0284751) q[0];
sx q[0];
rz(-2.3987781) q[0];
sx q[0];
rz(-1.9755367) q[0];
rz(-pi) q[1];
rz(-1.5624814) q[2];
sx q[2];
rz(-1.5247716) q[2];
sx q[2];
rz(-0.41551057) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.3869785) q[1];
sx q[1];
rz(-1.7474035) q[1];
sx q[1];
rz(0.28680027) q[1];
rz(-pi) q[2];
rz(1.9190556) q[3];
sx q[3];
rz(-1.5969443) q[3];
sx q[3];
rz(0.98046434) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.014331269) q[2];
sx q[2];
rz(-1.1703337) q[2];
sx q[2];
rz(2.2793615) q[2];
rz(-0.47752738) q[3];
sx q[3];
rz(-1.1815485) q[3];
sx q[3];
rz(-2.2235218) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35571337) q[0];
sx q[0];
rz(-2.9861351) q[0];
sx q[0];
rz(0.73295897) q[0];
rz(-3.0015302) q[1];
sx q[1];
rz(-2.1439794) q[1];
sx q[1];
rz(-1.0345116) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2037172) q[0];
sx q[0];
rz(-2.8563742) q[0];
sx q[0];
rz(2.2144149) q[0];
x q[1];
rz(0.10266281) q[2];
sx q[2];
rz(-1.6409988) q[2];
sx q[2];
rz(-1.7662802) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.4701177) q[1];
sx q[1];
rz(-2.5744994) q[1];
sx q[1];
rz(-1.5017897) q[1];
rz(1.6789867) q[3];
sx q[3];
rz(-1.1547609) q[3];
sx q[3];
rz(-0.070852208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.2404279) q[2];
sx q[2];
rz(-1.9463836) q[2];
sx q[2];
rz(0.77587664) q[2];
rz(-2.4173229) q[3];
sx q[3];
rz(-2.3359719) q[3];
sx q[3];
rz(1.4340713) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4416606) q[0];
sx q[0];
rz(-2.3775546) q[0];
sx q[0];
rz(0.016816703) q[0];
rz(3.1230714) q[1];
sx q[1];
rz(-1.8073945) q[1];
sx q[1];
rz(-2.3628078) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.22395615) q[0];
sx q[0];
rz(-2.7178239) q[0];
sx q[0];
rz(0.64026041) q[0];
rz(-pi) q[1];
rz(-0.2741371) q[2];
sx q[2];
rz(-2.0275653) q[2];
sx q[2];
rz(-2.0704839) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.76965895) q[1];
sx q[1];
rz(-1.2049335) q[1];
sx q[1];
rz(-2.8918173) q[1];
rz(1.7917463) q[3];
sx q[3];
rz(-1.6337992) q[3];
sx q[3];
rz(0.27730478) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.6918216) q[2];
sx q[2];
rz(-1.8372953) q[2];
sx q[2];
rz(0.71643913) q[2];
rz(1.6843494) q[3];
sx q[3];
rz(-1.4102035) q[3];
sx q[3];
rz(1.289207) q[3];
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
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0338106) q[0];
sx q[0];
rz(-2.6066715) q[0];
sx q[0];
rz(0.15144908) q[0];
rz(-0.17627136) q[1];
sx q[1];
rz(-1.9468032) q[1];
sx q[1];
rz(0.72296468) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.82397126) q[0];
sx q[0];
rz(-2.9569607) q[0];
sx q[0];
rz(2.1872107) q[0];
x q[1];
rz(-1.9480013) q[2];
sx q[2];
rz(-1.4133487) q[2];
sx q[2];
rz(2.8709656) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.1180229) q[1];
sx q[1];
rz(-1.6260864) q[1];
sx q[1];
rz(1.0667332) q[1];
rz(-1.1687752) q[3];
sx q[3];
rz(-0.82377454) q[3];
sx q[3];
rz(2.0603501) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.4621949) q[2];
sx q[2];
rz(-1.2827736) q[2];
sx q[2];
rz(2.6160713) q[2];
rz(-0.28371352) q[3];
sx q[3];
rz(-1.107639) q[3];
sx q[3];
rz(2.7450558) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59239607) q[0];
sx q[0];
rz(-2.017673) q[0];
sx q[0];
rz(2.8155433) q[0];
rz(-1.0992959) q[1];
sx q[1];
rz(-2.9922843) q[1];
sx q[1];
rz(-0.86984787) q[1];
rz(2.2407871) q[2];
sx q[2];
rz(-1.1849891) q[2];
sx q[2];
rz(1.6580788) q[2];
rz(1.7952193) q[3];
sx q[3];
rz(-1.2524458) q[3];
sx q[3];
rz(2.8129775) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
