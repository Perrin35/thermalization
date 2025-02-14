OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.8134107) q[0];
sx q[0];
rz(-2.1017177) q[0];
sx q[0];
rz(-1.9591969) q[0];
rz(0.27529588) q[1];
sx q[1];
rz(7.2841865) q[1];
sx q[1];
rz(10.871973) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6436734) q[0];
sx q[0];
rz(-2.5324973) q[0];
sx q[0];
rz(2.4793371) q[0];
rz(-pi) q[1];
x q[1];
rz(0.30302466) q[2];
sx q[2];
rz(-0.87867295) q[2];
sx q[2];
rz(-0.36811583) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.027982546) q[1];
sx q[1];
rz(-2.01568) q[1];
sx q[1];
rz(1.0335733) q[1];
rz(-pi) q[2];
rz(2.0979904) q[3];
sx q[3];
rz(-1.8917682) q[3];
sx q[3];
rz(-2.2466602) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.64624661) q[2];
sx q[2];
rz(-2.4301811) q[2];
sx q[2];
rz(1.4900253) q[2];
rz(1.1680158) q[3];
sx q[3];
rz(-1.370859) q[3];
sx q[3];
rz(-1.542154) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1646351) q[0];
sx q[0];
rz(-0.63666207) q[0];
sx q[0];
rz(2.013999) q[0];
rz(-0.98059869) q[1];
sx q[1];
rz(-2.6041928) q[1];
sx q[1];
rz(-1.5375536) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6036768) q[0];
sx q[0];
rz(-1.5620908) q[0];
sx q[0];
rz(1.1786657) q[0];
x q[1];
rz(1.4784548) q[2];
sx q[2];
rz(-1.5545111) q[2];
sx q[2];
rz(-3.1402325) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.1852946) q[1];
sx q[1];
rz(-1.0377656) q[1];
sx q[1];
rz(0.1406488) q[1];
x q[2];
rz(-1.8530773) q[3];
sx q[3];
rz(-1.936463) q[3];
sx q[3];
rz(2.7136841) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.8424707) q[2];
sx q[2];
rz(-0.8490347) q[2];
sx q[2];
rz(-0.15057286) q[2];
rz(-1.8922837) q[3];
sx q[3];
rz(-0.95730296) q[3];
sx q[3];
rz(-1.5540436) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7097178) q[0];
sx q[0];
rz(-1.4282325) q[0];
sx q[0];
rz(-0.10794434) q[0];
rz(2.3297294) q[1];
sx q[1];
rz(-0.64473647) q[1];
sx q[1];
rz(0.20400861) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1402748) q[0];
sx q[0];
rz(-1.55758) q[0];
sx q[0];
rz(-2.2030438) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1930253) q[2];
sx q[2];
rz(-1.1576414) q[2];
sx q[2];
rz(0.042095351) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.2141238) q[1];
sx q[1];
rz(-1.0368616) q[1];
sx q[1];
rz(-3.0005891) q[1];
rz(0.35514851) q[3];
sx q[3];
rz(-2.2844076) q[3];
sx q[3];
rz(2.5220152) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.418499) q[2];
sx q[2];
rz(-0.30947026) q[2];
sx q[2];
rz(1.2014368) q[2];
rz(2.5712683) q[3];
sx q[3];
rz(-1.0446965) q[3];
sx q[3];
rz(0.23049155) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0975321) q[0];
sx q[0];
rz(-0.965913) q[0];
sx q[0];
rz(-3.0934546) q[0];
rz(-2.8037181) q[1];
sx q[1];
rz(-2.477555) q[1];
sx q[1];
rz(-2.7545676) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.52464763) q[0];
sx q[0];
rz(-2.1957186) q[0];
sx q[0];
rz(2.8760103) q[0];
rz(-pi) q[1];
rz(-2.0104644) q[2];
sx q[2];
rz(-0.8789458) q[2];
sx q[2];
rz(0.95702167) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.6070575) q[1];
sx q[1];
rz(-1.1967772) q[1];
sx q[1];
rz(2.1675571) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0122672) q[3];
sx q[3];
rz(-2.8366907) q[3];
sx q[3];
rz(-2.4559742) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7016697) q[2];
sx q[2];
rz(-2.5950626) q[2];
sx q[2];
rz(-2.2129391) q[2];
rz(-0.52740151) q[3];
sx q[3];
rz(-1.4832067) q[3];
sx q[3];
rz(1.5206913) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.034390282) q[0];
sx q[0];
rz(-0.90029383) q[0];
sx q[0];
rz(0.1874371) q[0];
rz(-1.772359) q[1];
sx q[1];
rz(-2.1010459) q[1];
sx q[1];
rz(-1.7231411) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47309067) q[0];
sx q[0];
rz(-0.62431931) q[0];
sx q[0];
rz(-1.6416277) q[0];
x q[1];
rz(2.4985053) q[2];
sx q[2];
rz(-1.4355112) q[2];
sx q[2];
rz(-0.39150086) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.95067731) q[1];
sx q[1];
rz(-0.29914221) q[1];
sx q[1];
rz(-1.0043112) q[1];
x q[2];
rz(0.62331919) q[3];
sx q[3];
rz(-1.6409564) q[3];
sx q[3];
rz(-2.1753785) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.15293846) q[2];
sx q[2];
rz(-2.8230437) q[2];
sx q[2];
rz(0.62206507) q[2];
rz(2.7234744) q[3];
sx q[3];
rz(-1.6303635) q[3];
sx q[3];
rz(2.6869584) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.37586) q[0];
sx q[0];
rz(-1.3272165) q[0];
sx q[0];
rz(-0.8263998) q[0];
rz(1.3580258) q[1];
sx q[1];
rz(-0.91550255) q[1];
sx q[1];
rz(0.079039097) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5306471) q[0];
sx q[0];
rz(-1.4904596) q[0];
sx q[0];
rz(-1.4974844) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0108931) q[2];
sx q[2];
rz(-0.93564763) q[2];
sx q[2];
rz(-3.1237313) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.0896645) q[1];
sx q[1];
rz(-0.66510495) q[1];
sx q[1];
rz(2.0716578) q[1];
rz(2.6812234) q[3];
sx q[3];
rz(-1.7475024) q[3];
sx q[3];
rz(-1.9702156) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.11321154) q[2];
sx q[2];
rz(-1.6190642) q[2];
sx q[2];
rz(-1.1979423) q[2];
rz(-2.1146663) q[3];
sx q[3];
rz(-2.9668861) q[3];
sx q[3];
rz(-3.0026156) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.86676627) q[0];
sx q[0];
rz(-0.7678031) q[0];
sx q[0];
rz(0.36198947) q[0];
rz(-2.844574) q[1];
sx q[1];
rz(-1.8880918) q[1];
sx q[1];
rz(-0.46989918) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36951652) q[0];
sx q[0];
rz(-0.84699357) q[0];
sx q[0];
rz(2.2053201) q[0];
rz(-pi) q[1];
rz(-3.1387822) q[2];
sx q[2];
rz(-0.96544525) q[2];
sx q[2];
rz(-2.6058818) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-3.1220144) q[1];
sx q[1];
rz(-1.2570001) q[1];
sx q[1];
rz(-0.0079082735) q[1];
x q[2];
rz(-1.6401902) q[3];
sx q[3];
rz(-2.5186144) q[3];
sx q[3];
rz(-1.5322073) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.517259) q[2];
sx q[2];
rz(-0.91906157) q[2];
sx q[2];
rz(2.2993235) q[2];
rz(1.9159348) q[3];
sx q[3];
rz(-1.3282447) q[3];
sx q[3];
rz(-1.3262424) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
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
rz(0.4351742) q[0];
sx q[0];
rz(-1.5442727) q[0];
sx q[0];
rz(2.401974) q[0];
rz(2.836272) q[1];
sx q[1];
rz(-1.0605597) q[1];
sx q[1];
rz(1.7230497) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5180018) q[0];
sx q[0];
rz(-0.64096795) q[0];
sx q[0];
rz(-2.8056754) q[0];
rz(-pi) q[1];
rz(0.16295675) q[2];
sx q[2];
rz(-0.96471404) q[2];
sx q[2];
rz(-1.9418756) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.7463136) q[1];
sx q[1];
rz(-2.6398217) q[1];
sx q[1];
rz(-0.65478874) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7867602) q[3];
sx q[3];
rz(-1.7056173) q[3];
sx q[3];
rz(-2.1444837) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.50489396) q[2];
sx q[2];
rz(-0.86504522) q[2];
sx q[2];
rz(-1.0587586) q[2];
rz(-2.5508358) q[3];
sx q[3];
rz(-0.96610132) q[3];
sx q[3];
rz(3.0920715) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(0.010794086) q[0];
sx q[0];
rz(-0.32242355) q[0];
sx q[0];
rz(-1.4149433) q[0];
rz(-2.2825799) q[1];
sx q[1];
rz(-2.450727) q[1];
sx q[1];
rz(2.1376999) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.067779439) q[0];
sx q[0];
rz(-1.3155903) q[0];
sx q[0];
rz(1.7408446) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2357227) q[2];
sx q[2];
rz(-0.94099578) q[2];
sx q[2];
rz(0.12459151) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.5145961) q[1];
sx q[1];
rz(-0.8452943) q[1];
sx q[1];
rz(1.4519361) q[1];
rz(-pi) q[2];
rz(1.6565768) q[3];
sx q[3];
rz(-2.3137233) q[3];
sx q[3];
rz(2.0606526) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.88177195) q[2];
sx q[2];
rz(-1.5871781) q[2];
sx q[2];
rz(-0.47015321) q[2];
rz(1.5035037) q[3];
sx q[3];
rz(-0.94109002) q[3];
sx q[3];
rz(2.9269384) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.692602) q[0];
sx q[0];
rz(-0.34994352) q[0];
sx q[0];
rz(1.0119447) q[0];
rz(-0.28255209) q[1];
sx q[1];
rz(-0.91767445) q[1];
sx q[1];
rz(0.73175398) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5846281) q[0];
sx q[0];
rz(-3.0935253) q[0];
sx q[0];
rz(-2.9352504) q[0];
rz(-pi) q[1];
rz(-2.9596161) q[2];
sx q[2];
rz(-1.739758) q[2];
sx q[2];
rz(2.6017058) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.080345237) q[1];
sx q[1];
rz(-0.75704439) q[1];
sx q[1];
rz(-2.5233242) q[1];
rz(-1.8153207) q[3];
sx q[3];
rz(-2.7774307) q[3];
sx q[3];
rz(-2.9579666) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.67873502) q[2];
sx q[2];
rz(-1.5827468) q[2];
sx q[2];
rz(-0.25657121) q[2];
rz(0.76777846) q[3];
sx q[3];
rz(-2.0753746) q[3];
sx q[3];
rz(-1.0919051) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.11138042) q[0];
sx q[0];
rz(-2.31507) q[0];
sx q[0];
rz(2.4872381) q[0];
rz(2.4334473) q[1];
sx q[1];
rz(-1.000052) q[1];
sx q[1];
rz(2.4194385) q[1];
rz(-1.234645) q[2];
sx q[2];
rz(-0.79094255) q[2];
sx q[2];
rz(-3.0002181) q[2];
rz(3.1307316) q[3];
sx q[3];
rz(-1.9739816) q[3];
sx q[3];
rz(-1.0021742) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
