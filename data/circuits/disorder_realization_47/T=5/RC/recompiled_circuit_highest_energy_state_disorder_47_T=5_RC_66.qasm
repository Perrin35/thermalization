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
rz(1.6506305) q[0];
sx q[0];
rz(-1.4181674) q[0];
sx q[0];
rz(1.6562847) q[0];
rz(-2.0909042) q[1];
sx q[1];
rz(-1.3991855) q[1];
sx q[1];
rz(-2.4032226) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6577539) q[0];
sx q[0];
rz(-1.8138348) q[0];
sx q[0];
rz(1.7263822) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9887669) q[2];
sx q[2];
rz(-1.9597561) q[2];
sx q[2];
rz(-2.9940384) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.3787427) q[1];
sx q[1];
rz(-2.0113011) q[1];
sx q[1];
rz(1.4164914) q[1];
x q[2];
rz(1.2655696) q[3];
sx q[3];
rz(-2.1681227) q[3];
sx q[3];
rz(2.0862153) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.38306132) q[2];
sx q[2];
rz(-2.8585275) q[2];
sx q[2];
rz(-2.7281249) q[2];
rz(2.5773898) q[3];
sx q[3];
rz(-1.6744303) q[3];
sx q[3];
rz(1.9570501) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41459945) q[0];
sx q[0];
rz(-2.0818384) q[0];
sx q[0];
rz(-0.10398908) q[0];
rz(3.0005786) q[1];
sx q[1];
rz(-2.4555989) q[1];
sx q[1];
rz(-2.7242421) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39825059) q[0];
sx q[0];
rz(-2.5210327) q[0];
sx q[0];
rz(2.1055566) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6029458) q[2];
sx q[2];
rz(-0.36395437) q[2];
sx q[2];
rz(3.0917496) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.1097393) q[1];
sx q[1];
rz(-1.268034) q[1];
sx q[1];
rz(-0.043655386) q[1];
rz(-pi) q[2];
rz(0.52882282) q[3];
sx q[3];
rz(-1.0695056) q[3];
sx q[3];
rz(-2.5426585) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.61791164) q[2];
sx q[2];
rz(-1.0289501) q[2];
sx q[2];
rz(2.9376612) q[2];
rz(1.5150874) q[3];
sx q[3];
rz(-2.6650059) q[3];
sx q[3];
rz(1.509607) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5390891) q[0];
sx q[0];
rz(-1.6747549) q[0];
sx q[0];
rz(-0.3983101) q[0];
rz(-2.0643945) q[1];
sx q[1];
rz(-0.64882433) q[1];
sx q[1];
rz(-2.5881252) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3368126) q[0];
sx q[0];
rz(-2.3051728) q[0];
sx q[0];
rz(-1.5848716) q[0];
x q[1];
rz(1.7856423) q[2];
sx q[2];
rz(-2.2707377) q[2];
sx q[2];
rz(-0.46589771) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.0981507) q[1];
sx q[1];
rz(-0.69313184) q[1];
sx q[1];
rz(-2.5767703) q[1];
rz(-0.74456498) q[3];
sx q[3];
rz(-0.25601632) q[3];
sx q[3];
rz(-0.42041966) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(3.0674949) q[2];
sx q[2];
rz(-0.40253887) q[2];
sx q[2];
rz(-1.2158166) q[2];
rz(2.1408234) q[3];
sx q[3];
rz(-1.7583022) q[3];
sx q[3];
rz(-1.8892939) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(0.85553402) q[0];
sx q[0];
rz(-2.0834041) q[0];
sx q[0];
rz(1.4372987) q[0];
rz(1.3358491) q[1];
sx q[1];
rz(-0.62594405) q[1];
sx q[1];
rz(2.6864973) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6949161) q[0];
sx q[0];
rz(-1.4908264) q[0];
sx q[0];
rz(2.5237501) q[0];
rz(-pi) q[1];
rz(2.9326523) q[2];
sx q[2];
rz(-1.9358228) q[2];
sx q[2];
rz(1.2301333) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-3.0008853) q[1];
sx q[1];
rz(-1.326084) q[1];
sx q[1];
rz(2.7974456) q[1];
rz(0.42174815) q[3];
sx q[3];
rz(-0.22796002) q[3];
sx q[3];
rz(-2.3123119) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.8670292) q[2];
sx q[2];
rz(-0.44963351) q[2];
sx q[2];
rz(-0.8320128) q[2];
rz(0.653382) q[3];
sx q[3];
rz(-0.91975206) q[3];
sx q[3];
rz(2.466989) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(1.5117383) q[0];
sx q[0];
rz(-1.8269704) q[0];
sx q[0];
rz(2.4526556) q[0];
rz(2.9298933) q[1];
sx q[1];
rz(-1.4392821) q[1];
sx q[1];
rz(0.45417085) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8261738) q[0];
sx q[0];
rz(-2.4103904) q[0];
sx q[0];
rz(0.8765799) q[0];
rz(-pi) q[1];
rz(2.0549704) q[2];
sx q[2];
rz(-0.6956898) q[2];
sx q[2];
rz(1.5351968) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.7579753) q[1];
sx q[1];
rz(-1.0838008) q[1];
sx q[1];
rz(2.3112554) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8332043) q[3];
sx q[3];
rz(-1.4091564) q[3];
sx q[3];
rz(1.9899881) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.6058558) q[2];
sx q[2];
rz(-1.3593295) q[2];
sx q[2];
rz(2.6306756) q[2];
rz(1.9262675) q[3];
sx q[3];
rz(-1.5769703) q[3];
sx q[3];
rz(-0.63466614) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.11584347) q[0];
sx q[0];
rz(-0.31084335) q[0];
sx q[0];
rz(2.6281443) q[0];
rz(2.2250941) q[1];
sx q[1];
rz(-1.1767574) q[1];
sx q[1];
rz(1.7880012) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3037864) q[0];
sx q[0];
rz(-0.45188658) q[0];
sx q[0];
rz(-0.19048283) q[0];
x q[1];
rz(-0.36717461) q[2];
sx q[2];
rz(-1.9556139) q[2];
sx q[2];
rz(-0.095331017) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.533355) q[1];
sx q[1];
rz(-1.5974733) q[1];
sx q[1];
rz(-0.4877301) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.94025375) q[3];
sx q[3];
rz(-3.0315786) q[3];
sx q[3];
rz(2.2785254) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.4765656) q[2];
sx q[2];
rz(-2.5883784) q[2];
sx q[2];
rz(3.108007) q[2];
rz(0.50470662) q[3];
sx q[3];
rz(-1.4387771) q[3];
sx q[3];
rz(-2.3869042) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.018709239) q[0];
sx q[0];
rz(-1.2241192) q[0];
sx q[0];
rz(-1.1705742) q[0];
rz(2.9508044) q[1];
sx q[1];
rz(-0.46563322) q[1];
sx q[1];
rz(-0.64291397) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.87827728) q[0];
sx q[0];
rz(-1.8366792) q[0];
sx q[0];
rz(0.10590597) q[0];
x q[1];
rz(-2.229081) q[2];
sx q[2];
rz(-2.7325919) q[2];
sx q[2];
rz(0.3106948) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.2819124) q[1];
sx q[1];
rz(-1.2676867) q[1];
sx q[1];
rz(-1.3202841) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4516109) q[3];
sx q[3];
rz(-1.9728807) q[3];
sx q[3];
rz(-1.6329671) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.1428895) q[2];
sx q[2];
rz(-0.10991749) q[2];
sx q[2];
rz(1.9996803) q[2];
rz(-3.111908) q[3];
sx q[3];
rz(-1.5342865) q[3];
sx q[3];
rz(0.77965411) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90758816) q[0];
sx q[0];
rz(-1.1827844) q[0];
sx q[0];
rz(1.8335861) q[0];
rz(-1.1169149) q[1];
sx q[1];
rz(-1.6447379) q[1];
sx q[1];
rz(0.21211472) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.708688) q[0];
sx q[0];
rz(-1.2554662) q[0];
sx q[0];
rz(-0.84724119) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8050952) q[2];
sx q[2];
rz(-1.8673225) q[2];
sx q[2];
rz(-1.4594452) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.6740283) q[1];
sx q[1];
rz(-1.0564359) q[1];
sx q[1];
rz(2.6326724) q[1];
rz(-pi) q[2];
rz(2.2905825) q[3];
sx q[3];
rz(-2.27423) q[3];
sx q[3];
rz(1.2860166) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.1954631) q[2];
sx q[2];
rz(-1.3529494) q[2];
sx q[2];
rz(-1.8606261) q[2];
rz(-1.7383176) q[3];
sx q[3];
rz(-2.2303228) q[3];
sx q[3];
rz(0.72428552) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.69287777) q[0];
sx q[0];
rz(-0.71664482) q[0];
sx q[0];
rz(0.88248673) q[0];
rz(0.29300434) q[1];
sx q[1];
rz(-0.92331433) q[1];
sx q[1];
rz(-2.015347) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4871656) q[0];
sx q[0];
rz(-1.9357271) q[0];
sx q[0];
rz(2.4449744) q[0];
rz(-pi) q[1];
x q[1];
rz(3.1358997) q[2];
sx q[2];
rz(-2.5416592) q[2];
sx q[2];
rz(1.987074) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.3224015) q[1];
sx q[1];
rz(-1.5171851) q[1];
sx q[1];
rz(0.80162572) q[1];
rz(-pi) q[2];
rz(-2.7357941) q[3];
sx q[3];
rz(-1.8631336) q[3];
sx q[3];
rz(0.97989156) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.61665159) q[2];
sx q[2];
rz(-2.5280819) q[2];
sx q[2];
rz(2.9743527) q[2];
rz(3.1104769) q[3];
sx q[3];
rz(-1.9626706) q[3];
sx q[3];
rz(-0.64605609) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6328218) q[0];
sx q[0];
rz(-1.333586) q[0];
sx q[0];
rz(2.3311145) q[0];
rz(1.3088016) q[1];
sx q[1];
rz(-2.2056613) q[1];
sx q[1];
rz(-3.0442309) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0705436) q[0];
sx q[0];
rz(-2.3412421) q[0];
sx q[0];
rz(-2.9092824) q[0];
x q[1];
rz(1.951959) q[2];
sx q[2];
rz(-1.544017) q[2];
sx q[2];
rz(-1.8130515) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.78079911) q[1];
sx q[1];
rz(-2.153984) q[1];
sx q[1];
rz(-1.6287644) q[1];
rz(-1.3033569) q[3];
sx q[3];
rz(-1.4576266) q[3];
sx q[3];
rz(-1.3133996) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.8436766) q[2];
sx q[2];
rz(-1.0138136) q[2];
sx q[2];
rz(-1.8664912) q[2];
rz(0.44025931) q[3];
sx q[3];
rz(-0.85846725) q[3];
sx q[3];
rz(2.8662203) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.020582747) q[0];
sx q[0];
rz(-2.5204211) q[0];
sx q[0];
rz(-0.68516635) q[0];
rz(2.5195925) q[1];
sx q[1];
rz(-1.2983464) q[1];
sx q[1];
rz(1.3141528) q[1];
rz(1.4595819) q[2];
sx q[2];
rz(-1.6485452) q[2];
sx q[2];
rz(-2.9948276) q[2];
rz(1.3986599) q[3];
sx q[3];
rz(-0.93440104) q[3];
sx q[3];
rz(1.222119) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
