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
rz(1.3923378) q[0];
sx q[0];
rz(-2.7633986) q[0];
sx q[0];
rz(-0.35051546) q[0];
rz(0.89138436) q[1];
sx q[1];
rz(3.9337629) q[1];
sx q[1];
rz(13.73929) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.354686) q[0];
sx q[0];
rz(-1.3808492) q[0];
sx q[0];
rz(2.9482916) q[0];
x q[1];
rz(0.15526659) q[2];
sx q[2];
rz(-1.7874996) q[2];
sx q[2];
rz(-2.7086176) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.4306372) q[1];
sx q[1];
rz(-1.733755) q[1];
sx q[1];
rz(-0.32117543) q[1];
rz(-pi) q[2];
x q[2];
rz(0.47178206) q[3];
sx q[3];
rz(-1.7241447) q[3];
sx q[3];
rz(2.9073496) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.9709836) q[2];
sx q[2];
rz(-1.8656518) q[2];
sx q[2];
rz(-0.21273908) q[2];
rz(-3.0563266) q[3];
sx q[3];
rz(-1.250896) q[3];
sx q[3];
rz(-1.2712449) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.79412115) q[0];
sx q[0];
rz(-2.329282) q[0];
sx q[0];
rz(1.043327) q[0];
rz(3.0087545) q[1];
sx q[1];
rz(-2.9265407) q[1];
sx q[1];
rz(-1.1588233) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2272233) q[0];
sx q[0];
rz(-2.717319) q[0];
sx q[0];
rz(2.7844564) q[0];
x q[1];
rz(0.01387502) q[2];
sx q[2];
rz(-1.4094947) q[2];
sx q[2];
rz(-0.32729766) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.093250153) q[1];
sx q[1];
rz(-2.5777317) q[1];
sx q[1];
rz(-2.6287931) q[1];
x q[2];
rz(2.8125127) q[3];
sx q[3];
rz(-1.1606154) q[3];
sx q[3];
rz(3.0966334) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.1966689) q[2];
sx q[2];
rz(-0.66201869) q[2];
sx q[2];
rz(-0.72466737) q[2];
rz(2.479018) q[3];
sx q[3];
rz(-2.1734838) q[3];
sx q[3];
rz(2.841943) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9579983) q[0];
sx q[0];
rz(-1.8738926) q[0];
sx q[0];
rz(2.2268353) q[0];
rz(0.58685189) q[1];
sx q[1];
rz(-1.4205168) q[1];
sx q[1];
rz(0.14952001) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.11598524) q[0];
sx q[0];
rz(-1.7113026) q[0];
sx q[0];
rz(-2.0951659) q[0];
rz(-pi) q[1];
rz(-1.2039886) q[2];
sx q[2];
rz(-1.2292394) q[2];
sx q[2];
rz(1.2836518) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.55184522) q[1];
sx q[1];
rz(-2.0974297) q[1];
sx q[1];
rz(-0.92411516) q[1];
x q[2];
rz(-2.2633576) q[3];
sx q[3];
rz(-2.1562169) q[3];
sx q[3];
rz(-0.4957605) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.84825039) q[2];
sx q[2];
rz(-1.205227) q[2];
sx q[2];
rz(-2.3411574) q[2];
rz(-2.6750001) q[3];
sx q[3];
rz(-2.0355909) q[3];
sx q[3];
rz(-2.485937) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7993497) q[0];
sx q[0];
rz(-1.8514587) q[0];
sx q[0];
rz(2.2829862) q[0];
rz(-0.41060064) q[1];
sx q[1];
rz(-2.938439) q[1];
sx q[1];
rz(2.1536749) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3666375) q[0];
sx q[0];
rz(-1.2420085) q[0];
sx q[0];
rz(0.15952296) q[0];
rz(0.12148492) q[2];
sx q[2];
rz(-1.0787691) q[2];
sx q[2];
rz(0.58771261) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.6776442) q[1];
sx q[1];
rz(-1.8419187) q[1];
sx q[1];
rz(-1.2605002) q[1];
x q[2];
rz(-0.32295042) q[3];
sx q[3];
rz(-1.3846372) q[3];
sx q[3];
rz(1.8017839) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.9504488) q[2];
sx q[2];
rz(-1.7473651) q[2];
sx q[2];
rz(-2.6960755) q[2];
rz(2.4281003) q[3];
sx q[3];
rz(-2.4092509) q[3];
sx q[3];
rz(-2.2215686) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5934061) q[0];
sx q[0];
rz(-2.0596518) q[0];
sx q[0];
rz(-1.5997546) q[0];
rz(-1.2708739) q[1];
sx q[1];
rz(-1.4022695) q[1];
sx q[1];
rz(2.612203) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4146136) q[0];
sx q[0];
rz(-1.8185934) q[0];
sx q[0];
rz(1.9290061) q[0];
rz(-pi) q[1];
rz(0.61754333) q[2];
sx q[2];
rz(-1.4055335) q[2];
sx q[2];
rz(2.3338855) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.8258851) q[1];
sx q[1];
rz(-1.8216019) q[1];
sx q[1];
rz(-2.6205089) q[1];
rz(-1.9560341) q[3];
sx q[3];
rz(-2.7636313) q[3];
sx q[3];
rz(0.015008275) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.9993837) q[2];
sx q[2];
rz(-1.1915519) q[2];
sx q[2];
rz(-1.3231529) q[2];
rz(0.38149825) q[3];
sx q[3];
rz(-2.0254717) q[3];
sx q[3];
rz(2.748446) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8489654) q[0];
sx q[0];
rz(-0.75858527) q[0];
sx q[0];
rz(-2.1204156) q[0];
rz(1.609833) q[1];
sx q[1];
rz(-0.94667089) q[1];
sx q[1];
rz(-0.80345947) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.96413104) q[0];
sx q[0];
rz(-1.2636856) q[0];
sx q[0];
rz(1.403128) q[0];
x q[1];
rz(2.95058) q[2];
sx q[2];
rz(-0.41897853) q[2];
sx q[2];
rz(1.6167757) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.4157214) q[1];
sx q[1];
rz(-0.5827924) q[1];
sx q[1];
rz(0.21907138) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.62688503) q[3];
sx q[3];
rz(-1.4895559) q[3];
sx q[3];
rz(-2.8089942) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.6959186) q[2];
sx q[2];
rz(-2.5542104) q[2];
sx q[2];
rz(2.7109801) q[2];
rz(0.63498354) q[3];
sx q[3];
rz(-0.018714232) q[3];
sx q[3];
rz(1.8733321) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9950614) q[0];
sx q[0];
rz(-0.54556161) q[0];
sx q[0];
rz(1.5437641) q[0];
rz(0.68430463) q[1];
sx q[1];
rz(-1.6938208) q[1];
sx q[1];
rz(-2.3421471) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.862683) q[0];
sx q[0];
rz(-0.076795243) q[0];
sx q[0];
rz(3.0712295) q[0];
x q[1];
rz(0.85549037) q[2];
sx q[2];
rz(-1.5134619) q[2];
sx q[2];
rz(2.3279026) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.23161665) q[1];
sx q[1];
rz(-1.5607395) q[1];
sx q[1];
rz(-0.46044989) q[1];
x q[2];
rz(-0.83752172) q[3];
sx q[3];
rz(-1.234741) q[3];
sx q[3];
rz(2.1934862) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.53576175) q[2];
sx q[2];
rz(-2.3039218) q[2];
sx q[2];
rz(-2.3568995) q[2];
rz(-0.11624087) q[3];
sx q[3];
rz(-1.616547) q[3];
sx q[3];
rz(-1.2522662) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
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
rz(-0.5410974) q[0];
sx q[0];
rz(-1.9116115) q[0];
sx q[0];
rz(-2.1121209) q[0];
rz(-0.12044278) q[1];
sx q[1];
rz(-1.8414626) q[1];
sx q[1];
rz(-2.5708503) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.014387696) q[0];
sx q[0];
rz(-0.61349166) q[0];
sx q[0];
rz(0.86742371) q[0];
rz(-pi) q[1];
rz(2.4230401) q[2];
sx q[2];
rz(-1.4863803) q[2];
sx q[2];
rz(0.29044232) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.9214258) q[1];
sx q[1];
rz(-1.1296185) q[1];
sx q[1];
rz(0.36018546) q[1];
rz(-1.086497) q[3];
sx q[3];
rz(-2.7653793) q[3];
sx q[3];
rz(1.0602601) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.268078) q[2];
sx q[2];
rz(-2.2669078) q[2];
sx q[2];
rz(-0.52379215) q[2];
rz(1.1526456) q[3];
sx q[3];
rz(-0.58061424) q[3];
sx q[3];
rz(-1.7136278) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48083392) q[0];
sx q[0];
rz(-1.2013712) q[0];
sx q[0];
rz(0.42385605) q[0];
rz(-0.86680523) q[1];
sx q[1];
rz(-1.024217) q[1];
sx q[1];
rz(0.20326916) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0829859) q[0];
sx q[0];
rz(-1.6174249) q[0];
sx q[0];
rz(-0.65445047) q[0];
x q[1];
rz(1.5166984) q[2];
sx q[2];
rz(-1.8803758) q[2];
sx q[2];
rz(-0.54540173) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.31873676) q[1];
sx q[1];
rz(-1.7137495) q[1];
sx q[1];
rz(0.28495423) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8762693) q[3];
sx q[3];
rz(-1.8169699) q[3];
sx q[3];
rz(0.15428972) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.0554492) q[2];
sx q[2];
rz(-1.2148427) q[2];
sx q[2];
rz(0.27935585) q[2];
rz(-1.8481567) q[3];
sx q[3];
rz(-1.3118298) q[3];
sx q[3];
rz(-1.2156585) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9781037) q[0];
sx q[0];
rz(-2.3390529) q[0];
sx q[0];
rz(-1.4168903) q[0];
rz(-2.2143927) q[1];
sx q[1];
rz(-1.5617153) q[1];
sx q[1];
rz(-1.7318116) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5418842) q[0];
sx q[0];
rz(-2.4711631) q[0];
sx q[0];
rz(1.7267358) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0538834) q[2];
sx q[2];
rz(-1.8894686) q[2];
sx q[2];
rz(-2.0234194) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.3626707) q[1];
sx q[1];
rz(-0.92423981) q[1];
sx q[1];
rz(1.0151267) q[1];
x q[2];
rz(1.7396443) q[3];
sx q[3];
rz(-2.8956684) q[3];
sx q[3];
rz(-0.84154522) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.37567821) q[2];
sx q[2];
rz(-1.8199074) q[2];
sx q[2];
rz(-2.1709757) q[2];
rz(-2.4961903) q[3];
sx q[3];
rz(-0.97964764) q[3];
sx q[3];
rz(-1.0055044) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.29393016) q[0];
sx q[0];
rz(-1.6456589) q[0];
sx q[0];
rz(1.8346067) q[0];
rz(2.7597799) q[1];
sx q[1];
rz(-1.3419071) q[1];
sx q[1];
rz(0.76795427) q[1];
rz(-0.51580372) q[2];
sx q[2];
rz(-2.2602409) q[2];
sx q[2];
rz(0.088574499) q[2];
rz(1.1102724) q[3];
sx q[3];
rz(-1.9721748) q[3];
sx q[3];
rz(1.0655793) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
