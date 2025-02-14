OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.72803175) q[0];
sx q[0];
rz(-0.95946884) q[0];
sx q[0];
rz(2.5897107) q[0];
rz(-0.38504398) q[1];
sx q[1];
rz(-1.3621962) q[1];
sx q[1];
rz(0.41866067) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6543579) q[0];
sx q[0];
rz(-2.4208768) q[0];
sx q[0];
rz(2.0134175) q[0];
rz(-pi) q[1];
rz(-0.44960449) q[2];
sx q[2];
rz(-1.4105964) q[2];
sx q[2];
rz(-0.9148324) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.8238099) q[1];
sx q[1];
rz(-1.6595073) q[1];
sx q[1];
rz(2.937708) q[1];
rz(-pi) q[2];
rz(0.017114279) q[3];
sx q[3];
rz(-1.5403707) q[3];
sx q[3];
rz(2.3516114) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.60575134) q[2];
sx q[2];
rz(-1.8391106) q[2];
sx q[2];
rz(-0.38899404) q[2];
rz(-0.89748663) q[3];
sx q[3];
rz(-0.56098452) q[3];
sx q[3];
rz(-1.2787904) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
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
rz(1.2314583) q[0];
sx q[0];
rz(-0.95160216) q[0];
sx q[0];
rz(-2.8125473) q[0];
rz(1.344205) q[1];
sx q[1];
rz(-0.75733328) q[1];
sx q[1];
rz(-0.12408852) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24606516) q[0];
sx q[0];
rz(-2.82282) q[0];
sx q[0];
rz(-2.2970339) q[0];
rz(-pi) q[1];
rz(2.8555238) q[2];
sx q[2];
rz(-1.7033615) q[2];
sx q[2];
rz(-0.90434597) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.263843) q[1];
sx q[1];
rz(-1.9274492) q[1];
sx q[1];
rz(-2.043173) q[1];
x q[2];
rz(2.2749008) q[3];
sx q[3];
rz(-1.7287489) q[3];
sx q[3];
rz(2.7622836) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.8590392) q[2];
sx q[2];
rz(-0.48290792) q[2];
sx q[2];
rz(0.60749751) q[2];
rz(-0.88827682) q[3];
sx q[3];
rz(-1.5253303) q[3];
sx q[3];
rz(-1.4265149) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68179503) q[0];
sx q[0];
rz(-1.8495704) q[0];
sx q[0];
rz(2.9070396) q[0];
rz(1.0598496) q[1];
sx q[1];
rz(-1.9748961) q[1];
sx q[1];
rz(2.2134728) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7674196) q[0];
sx q[0];
rz(-2.1819127) q[0];
sx q[0];
rz(-2.0815297) q[0];
rz(1.9869204) q[2];
sx q[2];
rz(-1.1319931) q[2];
sx q[2];
rz(-3.0896387) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.8660714) q[1];
sx q[1];
rz(-1.8981985) q[1];
sx q[1];
rz(-1.7129461) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.15367963) q[3];
sx q[3];
rz(-1.2046332) q[3];
sx q[3];
rz(-1.2602214) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.5782535) q[2];
sx q[2];
rz(-1.4207062) q[2];
sx q[2];
rz(-1.2467747) q[2];
rz(-2.9747544) q[3];
sx q[3];
rz(-2.2127547) q[3];
sx q[3];
rz(-1.5290574) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4406776) q[0];
sx q[0];
rz(-3.0755141) q[0];
sx q[0];
rz(-0.75827688) q[0];
rz(1.5658763) q[1];
sx q[1];
rz(-0.58454746) q[1];
sx q[1];
rz(-0.87055269) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1567932) q[0];
sx q[0];
rz(-2.198303) q[0];
sx q[0];
rz(-2.1294562) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4192018) q[2];
sx q[2];
rz(-2.4263627) q[2];
sx q[2];
rz(1.4211637) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.9782171) q[1];
sx q[1];
rz(-1.723147) q[1];
sx q[1];
rz(-1.9780897) q[1];
x q[2];
rz(-1.4985256) q[3];
sx q[3];
rz(-2.1974583) q[3];
sx q[3];
rz(2.5570803) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.7190651) q[2];
sx q[2];
rz(-1.0079404) q[2];
sx q[2];
rz(2.3818805) q[2];
rz(0.64940137) q[3];
sx q[3];
rz(-1.6067182) q[3];
sx q[3];
rz(1.5788797) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0287057) q[0];
sx q[0];
rz(-2.4876471) q[0];
sx q[0];
rz(1.1829859) q[0];
rz(0.68823367) q[1];
sx q[1];
rz(-1.3166683) q[1];
sx q[1];
rz(0.36929718) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0645341) q[0];
sx q[0];
rz(-2.4774533) q[0];
sx q[0];
rz(0.044483552) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.30226548) q[2];
sx q[2];
rz(-2.4068953) q[2];
sx q[2];
rz(-0.96978984) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-3.0329513) q[1];
sx q[1];
rz(-0.88283112) q[1];
sx q[1];
rz(-2.9017519) q[1];
rz(-2.3595722) q[3];
sx q[3];
rz(-0.45392515) q[3];
sx q[3];
rz(2.5097414) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.94283048) q[2];
sx q[2];
rz(-1.3254712) q[2];
sx q[2];
rz(0.99786264) q[2];
rz(0.61740795) q[3];
sx q[3];
rz(-2.182775) q[3];
sx q[3];
rz(-2.3203885) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.038641039) q[0];
sx q[0];
rz(-1.2741673) q[0];
sx q[0];
rz(-1.8983023) q[0];
rz(-0.78168166) q[1];
sx q[1];
rz(-1.8676753) q[1];
sx q[1];
rz(-2.8807358) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2170625) q[0];
sx q[0];
rz(-0.53981656) q[0];
sx q[0];
rz(-2.7227719) q[0];
rz(-pi) q[1];
rz(-0.86791188) q[2];
sx q[2];
rz(-0.54722584) q[2];
sx q[2];
rz(-2.1457714) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.1697537) q[1];
sx q[1];
rz(-1.1128281) q[1];
sx q[1];
rz(-2.9724246) q[1];
rz(-pi) q[2];
rz(-0.93918899) q[3];
sx q[3];
rz(-1.7429461) q[3];
sx q[3];
rz(0.84014713) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.47238749) q[2];
sx q[2];
rz(-2.278625) q[2];
sx q[2];
rz(2.6589987) q[2];
rz(0.11416301) q[3];
sx q[3];
rz(-2.0518905) q[3];
sx q[3];
rz(-2.193006) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.362185) q[0];
sx q[0];
rz(-3.0857093) q[0];
sx q[0];
rz(0.26595297) q[0];
rz(-0.42568046) q[1];
sx q[1];
rz(-1.2454147) q[1];
sx q[1];
rz(-2.6735305) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.45586203) q[0];
sx q[0];
rz(-1.0294559) q[0];
sx q[0];
rz(-3.1336407) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0487814) q[2];
sx q[2];
rz(-1.2315237) q[2];
sx q[2];
rz(-2.3606127) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.1359766) q[1];
sx q[1];
rz(-1.3141201) q[1];
sx q[1];
rz(2.0784432) q[1];
rz(0.031074957) q[3];
sx q[3];
rz(-1.0508176) q[3];
sx q[3];
rz(-1.7620373) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.257306) q[2];
sx q[2];
rz(-1.851119) q[2];
sx q[2];
rz(-0.5536983) q[2];
rz(0.92343679) q[3];
sx q[3];
rz(-0.90673509) q[3];
sx q[3];
rz(3.0807909) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4629352) q[0];
sx q[0];
rz(-0.23660062) q[0];
sx q[0];
rz(2.4635354) q[0];
rz(-2.1642115) q[1];
sx q[1];
rz(-1.9934318) q[1];
sx q[1];
rz(-1.4378907) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17409731) q[0];
sx q[0];
rz(-2.6405442) q[0];
sx q[0];
rz(-1.1606085) q[0];
rz(0.0043785574) q[2];
sx q[2];
rz(-2.5776049) q[2];
sx q[2];
rz(2.0019238) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.4297144) q[1];
sx q[1];
rz(-0.3466045) q[1];
sx q[1];
rz(2.4019353) q[1];
rz(2.3947761) q[3];
sx q[3];
rz(-1.1594611) q[3];
sx q[3];
rz(2.8432027) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.73072726) q[2];
sx q[2];
rz(-0.56034708) q[2];
sx q[2];
rz(2.4746258) q[2];
rz(0.11735958) q[3];
sx q[3];
rz(-1.2753692) q[3];
sx q[3];
rz(-1.7608775) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9475107) q[0];
sx q[0];
rz(-1.5503333) q[0];
sx q[0];
rz(0.34341735) q[0];
rz(-1.4944448) q[1];
sx q[1];
rz(-2.1518555) q[1];
sx q[1];
rz(0.48430482) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0665047) q[0];
sx q[0];
rz(-2.4874788) q[0];
sx q[0];
rz(0.44371407) q[0];
rz(1.1555919) q[2];
sx q[2];
rz(-1.1744497) q[2];
sx q[2];
rz(-1.1434778) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.6086413) q[1];
sx q[1];
rz(-1.7923454) q[1];
sx q[1];
rz(2.5480342) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8922594) q[3];
sx q[3];
rz(-1.7847381) q[3];
sx q[3];
rz(2.6765424) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.9453498) q[2];
sx q[2];
rz(-2.0413155) q[2];
sx q[2];
rz(-2.2958882) q[2];
rz(1.2196994) q[3];
sx q[3];
rz(-1.2518576) q[3];
sx q[3];
rz(2.9367101) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6824816) q[0];
sx q[0];
rz(-0.24523188) q[0];
sx q[0];
rz(2.5526175) q[0];
rz(-0.67539769) q[1];
sx q[1];
rz(-2.1928936) q[1];
sx q[1];
rz(1.4640456) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1896473) q[0];
sx q[0];
rz(-0.69301096) q[0];
sx q[0];
rz(-2.243744) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4824994) q[2];
sx q[2];
rz(-0.73511926) q[2];
sx q[2];
rz(2.5500848) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.78641075) q[1];
sx q[1];
rz(-0.9847329) q[1];
sx q[1];
rz(-1.5152452) q[1];
rz(-pi) q[2];
rz(-0.69461639) q[3];
sx q[3];
rz(-1.2053849) q[3];
sx q[3];
rz(1.7323187) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.7318763) q[2];
sx q[2];
rz(-1.7475374) q[2];
sx q[2];
rz(1.124292) q[2];
rz(2.2935947) q[3];
sx q[3];
rz(-0.8046937) q[3];
sx q[3];
rz(0.030357411) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.547434) q[0];
sx q[0];
rz(-1.1347329) q[0];
sx q[0];
rz(-1.5190079) q[0];
rz(2.2183954) q[1];
sx q[1];
rz(-2.1122439) q[1];
sx q[1];
rz(-1.5069638) q[1];
rz(-0.30994762) q[2];
sx q[2];
rz(-1.7480231) q[2];
sx q[2];
rz(2.9668948) q[2];
rz(-2.2062929) q[3];
sx q[3];
rz(-2.7705396) q[3];
sx q[3];
rz(-2.0075575) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
