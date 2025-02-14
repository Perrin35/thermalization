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
rz(1.0518987) q[0];
sx q[0];
rz(2.9419152) q[0];
sx q[0];
rz(10.022104) q[0];
rz(2.6019793) q[1];
sx q[1];
rz(-2.5909254) q[1];
sx q[1];
rz(-0.82672969) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0866643) q[0];
sx q[0];
rz(-2.5205527) q[0];
sx q[0];
rz(-0.6426471) q[0];
rz(-pi) q[1];
rz(2.8887822) q[2];
sx q[2];
rz(-0.86517715) q[2];
sx q[2];
rz(-2.9575612) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.93365) q[1];
sx q[1];
rz(-0.68653715) q[1];
sx q[1];
rz(1.0892434) q[1];
rz(-pi) q[2];
rz(2.9197745) q[3];
sx q[3];
rz(-2.6008587) q[3];
sx q[3];
rz(-2.3774035) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.4248767) q[2];
sx q[2];
rz(-1.1307319) q[2];
sx q[2];
rz(-2.7190599) q[2];
rz(-1.2381964) q[3];
sx q[3];
rz(-2.7178552) q[3];
sx q[3];
rz(2.197926) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0331405) q[0];
sx q[0];
rz(-2.5863681) q[0];
sx q[0];
rz(-0.92414498) q[0];
rz(0.67584258) q[1];
sx q[1];
rz(-1.6249514) q[1];
sx q[1];
rz(1.0327551) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7387787) q[0];
sx q[0];
rz(-2.1354851) q[0];
sx q[0];
rz(-1.1346456) q[0];
rz(-pi) q[1];
x q[1];
rz(0.17272213) q[2];
sx q[2];
rz(-1.5859436) q[2];
sx q[2];
rz(-2.8081389) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.101866) q[1];
sx q[1];
rz(-1.1168774) q[1];
sx q[1];
rz(-0.32923876) q[1];
x q[2];
rz(-2.5268484) q[3];
sx q[3];
rz(-0.82332506) q[3];
sx q[3];
rz(0.73574443) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.5702901) q[2];
sx q[2];
rz(-1.4455659) q[2];
sx q[2];
rz(-2.9147713) q[2];
rz(1.4443385) q[3];
sx q[3];
rz(-0.11059136) q[3];
sx q[3];
rz(0.9816696) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.084376) q[0];
sx q[0];
rz(-2.9380874) q[0];
sx q[0];
rz(-0.89619005) q[0];
rz(2.8016688) q[1];
sx q[1];
rz(-1.8819921) q[1];
sx q[1];
rz(-0.47421727) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8086169) q[0];
sx q[0];
rz(-1.364437) q[0];
sx q[0];
rz(1.1550857) q[0];
rz(0.034696636) q[2];
sx q[2];
rz(-1.3857174) q[2];
sx q[2];
rz(0.68962956) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.6354628) q[1];
sx q[1];
rz(-2.3610507) q[1];
sx q[1];
rz(2.5012963) q[1];
rz(1.7960707) q[3];
sx q[3];
rz(-0.74639136) q[3];
sx q[3];
rz(2.3994544) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.077896) q[2];
sx q[2];
rz(-1.6113969) q[2];
sx q[2];
rz(-1.873675) q[2];
rz(0.39342132) q[3];
sx q[3];
rz(-0.370341) q[3];
sx q[3];
rz(2.2441277) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9758519) q[0];
sx q[0];
rz(-2.6065047) q[0];
sx q[0];
rz(0.4210327) q[0];
rz(-0.16972217) q[1];
sx q[1];
rz(-1.38849) q[1];
sx q[1];
rz(1.1441182) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.73070684) q[0];
sx q[0];
rz(-3.0907486) q[0];
sx q[0];
rz(2.6579141) q[0];
rz(-0.11387352) q[2];
sx q[2];
rz(-1.9673229) q[2];
sx q[2];
rz(0.20908326) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.3402945) q[1];
sx q[1];
rz(-2.3173875) q[1];
sx q[1];
rz(2.8293239) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.74419501) q[3];
sx q[3];
rz(-2.5928516) q[3];
sx q[3];
rz(2.5732267) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.052553136) q[2];
sx q[2];
rz(-1.7567987) q[2];
sx q[2];
rz(2.4844737) q[2];
rz(-2.0402015) q[3];
sx q[3];
rz(-1.907932) q[3];
sx q[3];
rz(-1.5737981) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9927486) q[0];
sx q[0];
rz(-2.086047) q[0];
sx q[0];
rz(1.8788991) q[0];
rz(2.8008723) q[1];
sx q[1];
rz(-1.4422528) q[1];
sx q[1];
rz(-2.5943894) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.094484821) q[0];
sx q[0];
rz(-1.5105627) q[0];
sx q[0];
rz(1.747986) q[0];
rz(-pi) q[1];
rz(-1.5261602) q[2];
sx q[2];
rz(-1.4603764) q[2];
sx q[2];
rz(2.8912889) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.2325523) q[1];
sx q[1];
rz(-2.4989933) q[1];
sx q[1];
rz(1.065083) q[1];
rz(2.2464348) q[3];
sx q[3];
rz(-0.46592281) q[3];
sx q[3];
rz(1.2873481) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.0146694) q[2];
sx q[2];
rz(-1.1875443) q[2];
sx q[2];
rz(0.32153258) q[2];
rz(-0.95162904) q[3];
sx q[3];
rz(-1.245433) q[3];
sx q[3];
rz(1.8754225) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5132009) q[0];
sx q[0];
rz(-2.7096847) q[0];
sx q[0];
rz(2.5950854) q[0];
rz(0.47738099) q[1];
sx q[1];
rz(-0.42196208) q[1];
sx q[1];
rz(1.2786) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1882478) q[0];
sx q[0];
rz(-0.45397568) q[0];
sx q[0];
rz(2.2843642) q[0];
rz(-pi) q[1];
rz(-1.5029656) q[2];
sx q[2];
rz(-2.1178341) q[2];
sx q[2];
rz(3.1233968) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.0562699) q[1];
sx q[1];
rz(-0.58093089) q[1];
sx q[1];
rz(-2.3694498) q[1];
x q[2];
rz(-1.8493722) q[3];
sx q[3];
rz(-2.4083353) q[3];
sx q[3];
rz(-2.5391308) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.64855856) q[2];
sx q[2];
rz(-0.15668046) q[2];
sx q[2];
rz(-0.9542166) q[2];
rz(-0.75718015) q[3];
sx q[3];
rz(-1.6482407) q[3];
sx q[3];
rz(-2.2455588) q[3];
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
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2467982) q[0];
sx q[0];
rz(-0.049012683) q[0];
sx q[0];
rz(3.0931296) q[0];
rz(-2.5080644) q[1];
sx q[1];
rz(-0.79094473) q[1];
sx q[1];
rz(1.6652426) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.033166151) q[0];
sx q[0];
rz(-1.3241741) q[0];
sx q[0];
rz(0.89065243) q[0];
rz(-pi) q[1];
x q[1];
rz(0.76870512) q[2];
sx q[2];
rz(-2.8787347) q[2];
sx q[2];
rz(2.4117087) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.4004835) q[1];
sx q[1];
rz(-1.2874075) q[1];
sx q[1];
rz(-2.7251935) q[1];
rz(-pi) q[2];
rz(-1.3562629) q[3];
sx q[3];
rz(-1.3140142) q[3];
sx q[3];
rz(-2.7329426) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.1518636) q[2];
sx q[2];
rz(-0.96111861) q[2];
sx q[2];
rz(2.9748919) q[2];
rz(1.9549595) q[3];
sx q[3];
rz(-1.9125166) q[3];
sx q[3];
rz(-2.1738539) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1879021) q[0];
sx q[0];
rz(-0.1583651) q[0];
sx q[0];
rz(2.4358391) q[0];
rz(1.3450274) q[1];
sx q[1];
rz(-1.2860362) q[1];
sx q[1];
rz(0.3261303) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6767384) q[0];
sx q[0];
rz(-0.87552545) q[0];
sx q[0];
rz(2.763163) q[0];
rz(-pi) q[1];
rz(3.1295849) q[2];
sx q[2];
rz(-2.084765) q[2];
sx q[2];
rz(-1.3552841) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.7520515) q[1];
sx q[1];
rz(-1.737575) q[1];
sx q[1];
rz(-0.15845815) q[1];
rz(-pi) q[2];
rz(1.2318065) q[3];
sx q[3];
rz(-2.7010446) q[3];
sx q[3];
rz(-0.46432981) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.2976133) q[2];
sx q[2];
rz(-1.16301) q[2];
sx q[2];
rz(-0.45912099) q[2];
rz(-1.8504359) q[3];
sx q[3];
rz(-1.2192817) q[3];
sx q[3];
rz(-3.0572157) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5463663) q[0];
sx q[0];
rz(-1.8084753) q[0];
sx q[0];
rz(3.1072733) q[0];
rz(-1.8345087) q[1];
sx q[1];
rz(-0.75051761) q[1];
sx q[1];
rz(-1.4141356) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3668418) q[0];
sx q[0];
rz(-2.2829608) q[0];
sx q[0];
rz(-1.2690498) q[0];
rz(2.6488536) q[2];
sx q[2];
rz(-0.73046267) q[2];
sx q[2];
rz(-1.3606537) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-3.0486949) q[1];
sx q[1];
rz(-2.1259397) q[1];
sx q[1];
rz(2.0616966) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8933252) q[3];
sx q[3];
rz(-1.5756996) q[3];
sx q[3];
rz(2.4941187) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.2109334) q[2];
sx q[2];
rz(-2.2948269) q[2];
sx q[2];
rz(0.47908121) q[2];
rz(0.72862285) q[3];
sx q[3];
rz(-1.9727547) q[3];
sx q[3];
rz(-0.59520477) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0625192) q[0];
sx q[0];
rz(-2.4010824) q[0];
sx q[0];
rz(0.24653521) q[0];
rz(-1.1277699) q[1];
sx q[1];
rz(-2.0083387) q[1];
sx q[1];
rz(0.22334982) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6205231) q[0];
sx q[0];
rz(-2.6436673) q[0];
sx q[0];
rz(-3.0715406) q[0];
rz(-pi) q[1];
rz(-2.6010314) q[2];
sx q[2];
rz(-1.7258919) q[2];
sx q[2];
rz(-2.0412091) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.57695147) q[1];
sx q[1];
rz(-1.511876) q[1];
sx q[1];
rz(1.8693187) q[1];
x q[2];
rz(-2.559313) q[3];
sx q[3];
rz(-0.49685541) q[3];
sx q[3];
rz(-0.41801807) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.8862306) q[2];
sx q[2];
rz(-0.88408154) q[2];
sx q[2];
rz(-0.13038005) q[2];
rz(-2.5731795) q[3];
sx q[3];
rz(-1.7805028) q[3];
sx q[3];
rz(0.4140678) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2100621) q[0];
sx q[0];
rz(-0.93183403) q[0];
sx q[0];
rz(0.43781042) q[0];
rz(1.1424278) q[1];
sx q[1];
rz(-1.2429968) q[1];
sx q[1];
rz(-1.4243781) q[1];
rz(2.3571711) q[2];
sx q[2];
rz(-2.5696345) q[2];
sx q[2];
rz(-0.77789236) q[2];
rz(-1.2503446) q[3];
sx q[3];
rz(-1.1966101) q[3];
sx q[3];
rz(0.93991652) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
