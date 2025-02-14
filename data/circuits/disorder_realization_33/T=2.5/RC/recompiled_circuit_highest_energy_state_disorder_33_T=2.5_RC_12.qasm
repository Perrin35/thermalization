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
rz(0.1641195) q[0];
sx q[0];
rz(-2.1174705) q[0];
sx q[0];
rz(0.48409387) q[0];
rz(-0.76264277) q[1];
sx q[1];
rz(-1.5781382) q[1];
sx q[1];
rz(0.054952316) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8380175) q[0];
sx q[0];
rz(-0.9009046) q[0];
sx q[0];
rz(-0.039867) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0829551) q[2];
sx q[2];
rz(-1.8203041) q[2];
sx q[2];
rz(2.3326258) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.7701227) q[1];
sx q[1];
rz(-2.9767163) q[1];
sx q[1];
rz(0.6368963) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1922429) q[3];
sx q[3];
rz(-1.8481405) q[3];
sx q[3];
rz(-2.8327033) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.310828) q[2];
sx q[2];
rz(-2.1711633) q[2];
sx q[2];
rz(-0.27466276) q[2];
rz(-0.87485391) q[3];
sx q[3];
rz(-2.0503876) q[3];
sx q[3];
rz(-2.8680475) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(pi/2) q[3];
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
rz(0.65983588) q[0];
sx q[0];
rz(-0.6898703) q[0];
sx q[0];
rz(0.39875317) q[0];
rz(-2.8975471) q[1];
sx q[1];
rz(-1.9052541) q[1];
sx q[1];
rz(-1.1057378) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.53351346) q[0];
sx q[0];
rz(-1.7987804) q[0];
sx q[0];
rz(0.45324676) q[0];
rz(-pi) q[1];
rz(-0.57184394) q[2];
sx q[2];
rz(-0.99858054) q[2];
sx q[2];
rz(-0.68293152) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.0827786) q[1];
sx q[1];
rz(-1.5734488) q[1];
sx q[1];
rz(3.1383731) q[1];
x q[2];
rz(1.1532743) q[3];
sx q[3];
rz(-2.8510749) q[3];
sx q[3];
rz(-0.49714303) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.4411321) q[2];
sx q[2];
rz(-1.9596142) q[2];
sx q[2];
rz(2.0665118) q[2];
rz(-0.69617802) q[3];
sx q[3];
rz(-1.8680365) q[3];
sx q[3];
rz(2.9785494) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8657846) q[0];
sx q[0];
rz(-1.2514665) q[0];
sx q[0];
rz(-2.9248917) q[0];
rz(1.7614583) q[1];
sx q[1];
rz(-2.7237027) q[1];
sx q[1];
rz(-2.5221141) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.86176819) q[0];
sx q[0];
rz(-3.0935043) q[0];
sx q[0];
rz(-2.2347941) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.046929788) q[2];
sx q[2];
rz(-1.4478193) q[2];
sx q[2];
rz(2.8356981) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.8113831) q[1];
sx q[1];
rz(-0.37322497) q[1];
sx q[1];
rz(-2.1887145) q[1];
rz(-0.9914753) q[3];
sx q[3];
rz(-1.0414755) q[3];
sx q[3];
rz(-0.73562276) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.2565903) q[2];
sx q[2];
rz(-1.8424748) q[2];
sx q[2];
rz(3.0510862) q[2];
rz(-0.45804405) q[3];
sx q[3];
rz(-2.6543255) q[3];
sx q[3];
rz(1.1080144) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
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
rz(1.3312382) q[0];
sx q[0];
rz(-2.2585223) q[0];
sx q[0];
rz(-0.93165398) q[0];
rz(-2.9725507) q[1];
sx q[1];
rz(-0.5642429) q[1];
sx q[1];
rz(1.0394675) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0213892) q[0];
sx q[0];
rz(-1.3303555) q[0];
sx q[0];
rz(-0.73657764) q[0];
rz(2.9157964) q[2];
sx q[2];
rz(-1.3232627) q[2];
sx q[2];
rz(0.087866656) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.4770866) q[1];
sx q[1];
rz(-1.5110104) q[1];
sx q[1];
rz(1.3633481) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5040197) q[3];
sx q[3];
rz(-0.3488144) q[3];
sx q[3];
rz(-1.3138527) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.52183759) q[2];
sx q[2];
rz(-1.1014742) q[2];
sx q[2];
rz(2.6452765) q[2];
rz(-0.79658341) q[3];
sx q[3];
rz(-2.8556672) q[3];
sx q[3];
rz(-2.0485785) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26688823) q[0];
sx q[0];
rz(-2.5570091) q[0];
sx q[0];
rz(0.97187483) q[0];
rz(-0.12174363) q[1];
sx q[1];
rz(-1.5309445) q[1];
sx q[1];
rz(-1.8121388) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15518256) q[0];
sx q[0];
rz(-1.2647795) q[0];
sx q[0];
rz(-3.0772692) q[0];
rz(-pi) q[1];
rz(2.0636286) q[2];
sx q[2];
rz(-1.2000286) q[2];
sx q[2];
rz(0.28647067) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.5962474) q[1];
sx q[1];
rz(-1.1584499) q[1];
sx q[1];
rz(-0.72280563) q[1];
rz(1.6498783) q[3];
sx q[3];
rz(-1.5615765) q[3];
sx q[3];
rz(-0.10445933) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.10151265) q[2];
sx q[2];
rz(-1.9209346) q[2];
sx q[2];
rz(0.15138781) q[2];
rz(-2.3769412) q[3];
sx q[3];
rz(-0.83573666) q[3];
sx q[3];
rz(-1.5966655) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(2.0695892) q[0];
sx q[0];
rz(-1.4729426) q[0];
sx q[0];
rz(-1.0714916) q[0];
rz(-2.6103861) q[1];
sx q[1];
rz(-1.6516282) q[1];
sx q[1];
rz(-2.5154617) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3817978) q[0];
sx q[0];
rz(-3.1239428) q[0];
sx q[0];
rz(1.8438898) q[0];
rz(1.8161723) q[2];
sx q[2];
rz(-0.98773709) q[2];
sx q[2];
rz(-0.32223216) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.8529359) q[1];
sx q[1];
rz(-0.96785883) q[1];
sx q[1];
rz(-1.1974242) q[1];
rz(-pi) q[2];
rz(-2.0141861) q[3];
sx q[3];
rz(-2.0282241) q[3];
sx q[3];
rz(1.2113038) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.7834187) q[2];
sx q[2];
rz(-0.61966115) q[2];
sx q[2];
rz(1.0142856) q[2];
rz(-1.2554393) q[3];
sx q[3];
rz(-1.7858601) q[3];
sx q[3];
rz(2.0881418) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1805873) q[0];
sx q[0];
rz(-2.2659681) q[0];
sx q[0];
rz(0.92051202) q[0];
rz(-0.11407425) q[1];
sx q[1];
rz(-1.736085) q[1];
sx q[1];
rz(2.0106409) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3823272) q[0];
sx q[0];
rz(-2.5438519) q[0];
sx q[0];
rz(3.1315342) q[0];
x q[1];
rz(0.99920706) q[2];
sx q[2];
rz(-2.0519369) q[2];
sx q[2];
rz(0.66758093) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.29207001) q[1];
sx q[1];
rz(-2.2300868) q[1];
sx q[1];
rz(-0.16772049) q[1];
rz(-pi) q[2];
x q[2];
rz(0.5163631) q[3];
sx q[3];
rz(-1.6404021) q[3];
sx q[3];
rz(0.35610896) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.34746927) q[2];
sx q[2];
rz(-2.4994714) q[2];
sx q[2];
rz(0.40880173) q[2];
rz(-0.18320228) q[3];
sx q[3];
rz(-1.5513523) q[3];
sx q[3];
rz(1.1387775) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.113753) q[0];
sx q[0];
rz(-3.0945859) q[0];
sx q[0];
rz(-0.29079944) q[0];
rz(-2.193702) q[1];
sx q[1];
rz(-0.46698505) q[1];
sx q[1];
rz(1.9416169) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24058293) q[0];
sx q[0];
rz(-1.173061) q[0];
sx q[0];
rz(0.57698864) q[0];
x q[1];
rz(-2.4418996) q[2];
sx q[2];
rz(-2.202575) q[2];
sx q[2];
rz(-0.80365411) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.7941798) q[1];
sx q[1];
rz(-0.29736667) q[1];
sx q[1];
rz(-1.7296687) q[1];
x q[2];
rz(1.0710217) q[3];
sx q[3];
rz(-2.933549) q[3];
sx q[3];
rz(2.2676454) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.7802508) q[2];
sx q[2];
rz(-1.6929408) q[2];
sx q[2];
rz(-1.7928436) q[2];
rz(-1.6382943) q[3];
sx q[3];
rz(-2.2595451) q[3];
sx q[3];
rz(2.9564814) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0272442) q[0];
sx q[0];
rz(-2.769727) q[0];
sx q[0];
rz(-2.0674904) q[0];
rz(0.4298003) q[1];
sx q[1];
rz(-1.0440412) q[1];
sx q[1];
rz(-0.46357402) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.060424711) q[0];
sx q[0];
rz(-1.3827818) q[0];
sx q[0];
rz(3.0449165) q[0];
rz(-pi) q[1];
rz(-1.7724636) q[2];
sx q[2];
rz(-1.8275765) q[2];
sx q[2];
rz(0.078291206) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-3.1164828) q[1];
sx q[1];
rz(-1.1949202) q[1];
sx q[1];
rz(-2.2459338) q[1];
rz(-pi) q[2];
rz(1.4809247) q[3];
sx q[3];
rz(-2.1276132) q[3];
sx q[3];
rz(-1.1695605) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.8672436) q[2];
sx q[2];
rz(-0.63259071) q[2];
sx q[2];
rz(0.80844936) q[2];
rz(2.6570053) q[3];
sx q[3];
rz(-2.7633568) q[3];
sx q[3];
rz(-2.7073879) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3633858) q[0];
sx q[0];
rz(-0.45553842) q[0];
sx q[0];
rz(-1.1325915) q[0];
rz(3.1032108) q[1];
sx q[1];
rz(-1.8033586) q[1];
sx q[1];
rz(-0.29676944) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1365741) q[0];
sx q[0];
rz(-1.2941735) q[0];
sx q[0];
rz(-2.121595) q[0];
rz(-pi) q[1];
rz(2.0149258) q[2];
sx q[2];
rz(-1.160356) q[2];
sx q[2];
rz(-0.12516147) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.3546875) q[1];
sx q[1];
rz(-0.2345492) q[1];
sx q[1];
rz(0.3292747) q[1];
rz(0.88480437) q[3];
sx q[3];
rz(-1.2594885) q[3];
sx q[3];
rz(-1.6303568) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.177629) q[2];
sx q[2];
rz(-1.12135) q[2];
sx q[2];
rz(-2.5753944) q[2];
rz(-1.1244134) q[3];
sx q[3];
rz(-3.0163613) q[3];
sx q[3];
rz(-1.9521149) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.88846702) q[0];
sx q[0];
rz(-2.4430226) q[0];
sx q[0];
rz(-3.0342614) q[0];
rz(-0.26168564) q[1];
sx q[1];
rz(-1.9410004) q[1];
sx q[1];
rz(0.95071361) q[1];
rz(-3.059584) q[2];
sx q[2];
rz(-1.5463943) q[2];
sx q[2];
rz(0.80873185) q[2];
rz(0.85121831) q[3];
sx q[3];
rz(-2.2891582) q[3];
sx q[3];
rz(-1.4386406) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
