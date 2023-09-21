OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.9010889) q[0];
sx q[0];
rz(-0.17740372) q[0];
sx q[0];
rz(-2.0071964) q[0];
rz(1.1881243) q[1];
sx q[1];
rz(-2.1048574) q[1];
sx q[1];
rz(-0.66361767) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3463319) q[0];
sx q[0];
rz(-1.594461) q[0];
sx q[0];
rz(-2.6803826) q[0];
rz(-pi) q[1];
rz(-0.3590091) q[2];
sx q[2];
rz(-0.81072545) q[2];
sx q[2];
rz(0.63149482) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.9259778) q[1];
sx q[1];
rz(-1.3091334) q[1];
sx q[1];
rz(1.2222626) q[1];
x q[2];
rz(0.10856467) q[3];
sx q[3];
rz(-1.3285471) q[3];
sx q[3];
rz(-3.0755175) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.2628281) q[2];
sx q[2];
rz(-0.44439134) q[2];
sx q[2];
rz(-0.051068548) q[2];
rz(-2.5845394) q[3];
sx q[3];
rz(-0.80018187) q[3];
sx q[3];
rz(-1.5867656) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5550845) q[0];
sx q[0];
rz(-2.3095135) q[0];
sx q[0];
rz(2.5449975) q[0];
rz(-2.3157628) q[1];
sx q[1];
rz(-1.4412216) q[1];
sx q[1];
rz(1.2260431) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81235028) q[0];
sx q[0];
rz(-1.3717522) q[0];
sx q[0];
rz(0.45254405) q[0];
rz(-pi) q[1];
rz(-2.3184899) q[2];
sx q[2];
rz(-2.1777993) q[2];
sx q[2];
rz(2.5603103) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.023149816) q[1];
sx q[1];
rz(-2.9217302) q[1];
sx q[1];
rz(1.5688194) q[1];
x q[2];
rz(0.0094718178) q[3];
sx q[3];
rz(-2.1433899) q[3];
sx q[3];
rz(-0.026281683) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.3423959) q[2];
sx q[2];
rz(-1.1725972) q[2];
sx q[2];
rz(-0.33102316) q[2];
rz(-2.3349169) q[3];
sx q[3];
rz(-1.2253864) q[3];
sx q[3];
rz(-1.7139009) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1598635) q[0];
sx q[0];
rz(-1.230343) q[0];
sx q[0];
rz(-1.8925517) q[0];
rz(-3.0535835) q[1];
sx q[1];
rz(-1.1227337) q[1];
sx q[1];
rz(-1.0294611) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2217076) q[0];
sx q[0];
rz(-1.3516597) q[0];
sx q[0];
rz(2.1973781) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8573895) q[2];
sx q[2];
rz(-0.9409875) q[2];
sx q[2];
rz(1.0857925) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.0440812) q[1];
sx q[1];
rz(-1.4315726) q[1];
sx q[1];
rz(-0.15686762) q[1];
rz(-pi) q[2];
rz(-3.0089278) q[3];
sx q[3];
rz(-2.3972315) q[3];
sx q[3];
rz(2.6877407) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(3.133698) q[2];
sx q[2];
rz(-1.4174613) q[2];
sx q[2];
rz(0.50764817) q[2];
rz(-1.7525904) q[3];
sx q[3];
rz(-2.8116083) q[3];
sx q[3];
rz(1.1631789) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65748173) q[0];
sx q[0];
rz(-0.27814516) q[0];
sx q[0];
rz(1.5959651) q[0];
rz(1.0428628) q[1];
sx q[1];
rz(-1.9680126) q[1];
sx q[1];
rz(-1.625659) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31947485) q[0];
sx q[0];
rz(-1.428077) q[0];
sx q[0];
rz(0.6832173) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4342986) q[2];
sx q[2];
rz(-2.578306) q[2];
sx q[2];
rz(-1.0922645) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-3.1277395) q[1];
sx q[1];
rz(-1.8426367) q[1];
sx q[1];
rz(0.66687648) q[1];
x q[2];
rz(0.32228542) q[3];
sx q[3];
rz(-2.5282113) q[3];
sx q[3];
rz(1.8005288) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.3778014) q[2];
sx q[2];
rz(-2.2968473) q[2];
sx q[2];
rz(-1.2949004) q[2];
rz(-3.0002248) q[3];
sx q[3];
rz(-0.54261345) q[3];
sx q[3];
rz(2.1550089) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.35448733) q[0];
sx q[0];
rz(-1.1802477) q[0];
sx q[0];
rz(1.5198583) q[0];
rz(-0.63201085) q[1];
sx q[1];
rz(-2.4193587) q[1];
sx q[1];
rz(-2.246726) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2762404) q[0];
sx q[0];
rz(-1.1603174) q[0];
sx q[0];
rz(2.4549237) q[0];
rz(0.61770265) q[2];
sx q[2];
rz(-0.82759826) q[2];
sx q[2];
rz(-0.4862116) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.85256165) q[1];
sx q[1];
rz(-0.82419306) q[1];
sx q[1];
rz(2.186071) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5785061) q[3];
sx q[3];
rz(-1.2843411) q[3];
sx q[3];
rz(-2.8849998) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.52508369) q[2];
sx q[2];
rz(-0.36642763) q[2];
sx q[2];
rz(1.1425225) q[2];
rz(2.3948495) q[3];
sx q[3];
rz(-1.8871504) q[3];
sx q[3];
rz(1.07553) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.58364761) q[0];
sx q[0];
rz(-0.32145158) q[0];
sx q[0];
rz(1.3775795) q[0];
rz(2.6691943) q[1];
sx q[1];
rz(-2.6230085) q[1];
sx q[1];
rz(-0.46498743) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2705921) q[0];
sx q[0];
rz(-2.5809079) q[0];
sx q[0];
rz(-2.2339348) q[0];
rz(-pi) q[1];
x q[1];
rz(0.34157413) q[2];
sx q[2];
rz(-1.5014868) q[2];
sx q[2];
rz(-0.14788936) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.38575129) q[1];
sx q[1];
rz(-2.5045536) q[1];
sx q[1];
rz(1.6511276) q[1];
x q[2];
rz(-1.3607849) q[3];
sx q[3];
rz(-2.2689153) q[3];
sx q[3];
rz(-0.39892808) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.49089367) q[2];
sx q[2];
rz(-1.8785672) q[2];
sx q[2];
rz(-2.6965551) q[2];
rz(0.93368357) q[3];
sx q[3];
rz(-1.7088339) q[3];
sx q[3];
rz(-2.8745108) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.12748195) q[0];
sx q[0];
rz(-1.575379) q[0];
sx q[0];
rz(1.4200462) q[0];
rz(-0.02380112) q[1];
sx q[1];
rz(-2.5271466) q[1];
sx q[1];
rz(-0.15596095) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.50990803) q[0];
sx q[0];
rz(-1.8242867) q[0];
sx q[0];
rz(1.7476728) q[0];
rz(2.4620352) q[2];
sx q[2];
rz(-0.62629269) q[2];
sx q[2];
rz(-1.9192413) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.8487726) q[1];
sx q[1];
rz(-2.9584868) q[1];
sx q[1];
rz(-1.8715026) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6050623) q[3];
sx q[3];
rz(-1.1324258) q[3];
sx q[3];
rz(0.28552548) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-3.0722787) q[2];
sx q[2];
rz(-0.57702714) q[2];
sx q[2];
rz(1.0726661) q[2];
rz(-2.8159451) q[3];
sx q[3];
rz(-1.1762534) q[3];
sx q[3];
rz(-1.6252888) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9496562) q[0];
sx q[0];
rz(-0.40238109) q[0];
sx q[0];
rz(0.33777133) q[0];
rz(-1.0900963) q[1];
sx q[1];
rz(-2.1665159) q[1];
sx q[1];
rz(-2.8930194) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.08911207) q[0];
sx q[0];
rz(-2.5626474) q[0];
sx q[0];
rz(0.46778932) q[0];
rz(-pi) q[1];
x q[1];
rz(0.48034251) q[2];
sx q[2];
rz(-2.0311653) q[2];
sx q[2];
rz(0.72887052) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.8057115) q[1];
sx q[1];
rz(-0.98166775) q[1];
sx q[1];
rz(-2.5475403) q[1];
rz(-1.1484654) q[3];
sx q[3];
rz(-0.18581192) q[3];
sx q[3];
rz(-1.8728135) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.8481855) q[2];
sx q[2];
rz(-0.53331393) q[2];
sx q[2];
rz(-0.08671134) q[2];
rz(2.6596206) q[3];
sx q[3];
rz(-0.50589365) q[3];
sx q[3];
rz(1.1530676) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6329353) q[0];
sx q[0];
rz(-1.0077227) q[0];
sx q[0];
rz(-2.8588262) q[0];
rz(-0.70156082) q[1];
sx q[1];
rz(-0.82156721) q[1];
sx q[1];
rz(1.3185906) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.628172) q[0];
sx q[0];
rz(-1.4903729) q[0];
sx q[0];
rz(-2.0128987) q[0];
rz(-pi) q[1];
x q[1];
rz(0.96005) q[2];
sx q[2];
rz(-2.1337482) q[2];
sx q[2];
rz(2.8658531) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.5732167) q[1];
sx q[1];
rz(-1.425256) q[1];
sx q[1];
rz(-2.2578866) q[1];
x q[2];
rz(-0.75026476) q[3];
sx q[3];
rz(-1.8599469) q[3];
sx q[3];
rz(0.17444785) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.41137722) q[2];
sx q[2];
rz(-2.8432379) q[2];
sx q[2];
rz(2.7424157) q[2];
rz(-2.2579851) q[3];
sx q[3];
rz(-1.6688321) q[3];
sx q[3];
rz(1.9201027) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76319641) q[0];
sx q[0];
rz(-2.3662687) q[0];
sx q[0];
rz(-1.5989074) q[0];
rz(1.0653161) q[1];
sx q[1];
rz(-2.1677446) q[1];
sx q[1];
rz(-1.261196) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1051837) q[0];
sx q[0];
rz(-1.1622218) q[0];
sx q[0];
rz(1.4559792) q[0];
x q[1];
rz(0.032632685) q[2];
sx q[2];
rz(-0.59851461) q[2];
sx q[2];
rz(-2.0805217) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.65834261) q[1];
sx q[1];
rz(-1.5066506) q[1];
sx q[1];
rz(0.49883962) q[1];
rz(-1.4114755) q[3];
sx q[3];
rz(-0.89309249) q[3];
sx q[3];
rz(0.89980984) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.5579055) q[2];
sx q[2];
rz(-2.5186899) q[2];
sx q[2];
rz(-0.56979257) q[2];
rz(-1.2184881) q[3];
sx q[3];
rz(-0.64703882) q[3];
sx q[3];
rz(2.5789554) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31310836) q[0];
sx q[0];
rz(-0.88043558) q[0];
sx q[0];
rz(1.2783929) q[0];
rz(0.60824153) q[1];
sx q[1];
rz(-2.6651762) q[1];
sx q[1];
rz(-2.6574635) q[1];
rz(-2.1424978) q[2];
sx q[2];
rz(-1.5928762) q[2];
sx q[2];
rz(1.5230509) q[2];
rz(-1.431987) q[3];
sx q[3];
rz(-0.98416735) q[3];
sx q[3];
rz(-2.9744801) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];