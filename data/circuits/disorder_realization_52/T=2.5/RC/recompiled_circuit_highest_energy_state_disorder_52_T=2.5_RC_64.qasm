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
rz(1.1531416) q[0];
sx q[0];
rz(-0.81557953) q[0];
sx q[0];
rz(2.3834035) q[0];
rz(-0.012501333) q[1];
sx q[1];
rz(4.9795436) q[1];
sx q[1];
rz(10.995168) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.16906315) q[0];
sx q[0];
rz(-0.95989812) q[0];
sx q[0];
rz(-1.1283895) q[0];
rz(-pi) q[1];
rz(-1.1711898) q[2];
sx q[2];
rz(-1.6328703) q[2];
sx q[2];
rz(0.8555748) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.6504892) q[1];
sx q[1];
rz(-1.2105618) q[1];
sx q[1];
rz(-0.061042518) q[1];
x q[2];
rz(-1.8800354) q[3];
sx q[3];
rz(-1.6282363) q[3];
sx q[3];
rz(-0.32306898) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.6115173) q[2];
sx q[2];
rz(-0.0081743058) q[2];
sx q[2];
rz(2.7008936) q[2];
rz(-0.079744451) q[3];
sx q[3];
rz(-3.1414882) q[3];
sx q[3];
rz(2.0174446) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1752862) q[0];
sx q[0];
rz(-0.056862406) q[0];
sx q[0];
rz(2.9735907) q[0];
rz(-3.1213144) q[1];
sx q[1];
rz(-0.30826491) q[1];
sx q[1];
rz(1.6049989) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.56242053) q[0];
sx q[0];
rz(-1.354166) q[0];
sx q[0];
rz(2.9013293) q[0];
rz(-pi) q[1];
rz(-3.1390983) q[2];
sx q[2];
rz(-1.5606784) q[2];
sx q[2];
rz(-1.5455139) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.3602996) q[1];
sx q[1];
rz(-1.5754682) q[1];
sx q[1];
rz(-0.0010990573) q[1];
x q[2];
rz(-1.5797938) q[3];
sx q[3];
rz(-1.5285042) q[3];
sx q[3];
rz(-2.472714) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.34540471) q[2];
sx q[2];
rz(-0.91190839) q[2];
sx q[2];
rz(-1.4076642) q[2];
rz(-2.0965072) q[3];
sx q[3];
rz(-3.0920691) q[3];
sx q[3];
rz(2.869587) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3097836) q[0];
sx q[0];
rz(-0.97667664) q[0];
sx q[0];
rz(-0.56104863) q[0];
rz(0.27684119) q[1];
sx q[1];
rz(-0.012877348) q[1];
sx q[1];
rz(-1.3078088) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9914068) q[0];
sx q[0];
rz(-1.6117818) q[0];
sx q[0];
rz(1.2738373) q[0];
rz(0.0073892825) q[2];
sx q[2];
rz(-1.5708013) q[2];
sx q[2];
rz(-0.50732343) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.0289291) q[1];
sx q[1];
rz(-1.512292) q[1];
sx q[1];
rz(-2.1444291) q[1];
x q[2];
rz(-2.6827319) q[3];
sx q[3];
rz(-2.1956177) q[3];
sx q[3];
rz(0.030908728) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.3448559) q[2];
sx q[2];
rz(-3.1414746) q[2];
sx q[2];
rz(-0.5893839) q[2];
rz(-1.2423337) q[3];
sx q[3];
rz(-3.1292249) q[3];
sx q[3];
rz(1.3602268) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1347443) q[0];
sx q[0];
rz(-2.6297748) q[0];
sx q[0];
rz(1.3486598) q[0];
rz(0.0066561247) q[1];
sx q[1];
rz(-1.321188) q[1];
sx q[1];
rz(-3.1080918) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2804256) q[0];
sx q[0];
rz(-0.67641947) q[0];
sx q[0];
rz(-1.2926213) q[0];
rz(-pi) q[1];
rz(1.5662996) q[2];
sx q[2];
rz(-1.6904313) q[2];
sx q[2];
rz(0.41882354) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.12043145) q[1];
sx q[1];
rz(-0.26909262) q[1];
sx q[1];
rz(-1.5325559) q[1];
rz(-1.7351522) q[3];
sx q[3];
rz(-1.4504315) q[3];
sx q[3];
rz(0.09935483) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.032430705) q[2];
sx q[2];
rz(-0.0065294821) q[2];
sx q[2];
rz(-0.29864857) q[2];
rz(-1.8715035) q[3];
sx q[3];
rz(-0.01615571) q[3];
sx q[3];
rz(-3.0884009) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.771516) q[0];
sx q[0];
rz(-1.5144441) q[0];
sx q[0];
rz(-2.694743) q[0];
rz(-2.9554548) q[1];
sx q[1];
rz(-0.061807241) q[1];
sx q[1];
rz(-1.4131379) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9139149) q[0];
sx q[0];
rz(-2.9083038) q[0];
sx q[0];
rz(-1.7167709) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1485708) q[2];
sx q[2];
rz(-1.3729501) q[2];
sx q[2];
rz(2.5289218) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.4889989) q[1];
sx q[1];
rz(-1.5212544) q[1];
sx q[1];
rz(-3.0964353) q[1];
rz(-pi) q[2];
rz(-2.2373393) q[3];
sx q[3];
rz(-2.7891141) q[3];
sx q[3];
rz(2.3942687) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.3073005) q[2];
sx q[2];
rz(-1.5895546) q[2];
sx q[2];
rz(-0.49868047) q[2];
rz(2.5636766) q[3];
sx q[3];
rz(-2.6579865) q[3];
sx q[3];
rz(-0.59670603) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96108288) q[0];
sx q[0];
rz(-1.120765) q[0];
sx q[0];
rz(0.34641308) q[0];
rz(-0.60180426) q[1];
sx q[1];
rz(-1.5608984) q[1];
sx q[1];
rz(2.3880889) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0929209) q[0];
sx q[0];
rz(-1.762383) q[0];
sx q[0];
rz(1.3956153) q[0];
rz(-pi) q[1];
rz(2.1624203) q[2];
sx q[2];
rz(-3.0074928) q[2];
sx q[2];
rz(0.025452415) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.9484472) q[1];
sx q[1];
rz(-2.3420534) q[1];
sx q[1];
rz(-2.1667797) q[1];
x q[2];
rz(-3.0493204) q[3];
sx q[3];
rz(-1.408395) q[3];
sx q[3];
rz(-2.5889461) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.57009131) q[2];
sx q[2];
rz(-3.1381021) q[2];
sx q[2];
rz(1.5241148) q[2];
rz(0.12024719) q[3];
sx q[3];
rz(-0.0032987981) q[3];
sx q[3];
rz(-0.53774589) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26871249) q[0];
sx q[0];
rz(-2.217642) q[0];
sx q[0];
rz(0.22126108) q[0];
rz(1.4606754) q[1];
sx q[1];
rz(-0.9333846) q[1];
sx q[1];
rz(3.0642919) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.59931) q[0];
sx q[0];
rz(-3.0998383) q[0];
sx q[0];
rz(1.5357273) q[0];
rz(-pi) q[1];
x q[1];
rz(3.1368106) q[2];
sx q[2];
rz(-1.5619229) q[2];
sx q[2];
rz(2.9298669) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.0822466) q[1];
sx q[1];
rz(-1.4993854) q[1];
sx q[1];
rz(2.9661428) q[1];
rz(3.0139913) q[3];
sx q[3];
rz(-1.3584879) q[3];
sx q[3];
rz(1.9636167) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.3495425) q[2];
sx q[2];
rz(-0.011186102) q[2];
sx q[2];
rz(-0.95996094) q[2];
rz(2.8121484) q[3];
sx q[3];
rz(-3.1335148) q[3];
sx q[3];
rz(2.2913057) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
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
rz(-2.195381) q[0];
sx q[0];
rz(-0.61310261) q[0];
sx q[0];
rz(3.0323113) q[0];
rz(0.37336135) q[1];
sx q[1];
rz(-2.3318718) q[1];
sx q[1];
rz(1.2304617) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4407935) q[0];
sx q[0];
rz(-2.1820118) q[0];
sx q[0];
rz(-0.81328765) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9443342) q[2];
sx q[2];
rz(-1.7575118) q[2];
sx q[2];
rz(-3.124539) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.8671678) q[1];
sx q[1];
rz(-1.503957) q[1];
sx q[1];
rz(3.1301366) q[1];
rz(-pi) q[2];
x q[2];
rz(0.91992141) q[3];
sx q[3];
rz(-1.8200785) q[3];
sx q[3];
rz(-1.9712308) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.5665148) q[2];
sx q[2];
rz(-1.9058303) q[2];
sx q[2];
rz(1.8129978) q[2];
rz(-1.7447507) q[3];
sx q[3];
rz(-0.003740398) q[3];
sx q[3];
rz(2.1197135) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0777271) q[0];
sx q[0];
rz(-1.6985748) q[0];
sx q[0];
rz(2.5685837) q[0];
rz(2.833448) q[1];
sx q[1];
rz(-2.7317218) q[1];
sx q[1];
rz(-1.0073957) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0125748) q[0];
sx q[0];
rz(-1.5652302) q[0];
sx q[0];
rz(1.6773423) q[0];
rz(1.8780872) q[2];
sx q[2];
rz(-2.9973534) q[2];
sx q[2];
rz(2.7968614) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.7511661) q[1];
sx q[1];
rz(-1.5301643) q[1];
sx q[1];
rz(1.6540113) q[1];
rz(-pi) q[2];
rz(1.5406403) q[3];
sx q[3];
rz(-1.5573676) q[3];
sx q[3];
rz(-1.3659988) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.31743) q[2];
sx q[2];
rz(-2.514826) q[2];
sx q[2];
rz(2.7516348) q[2];
rz(0.069084875) q[3];
sx q[3];
rz(-3.1324813) q[3];
sx q[3];
rz(-2.8100584) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.020141715) q[0];
sx q[0];
rz(-2.3905601) q[0];
sx q[0];
rz(-0.48594117) q[0];
rz(2.2700229) q[1];
sx q[1];
rz(-1.8337245) q[1];
sx q[1];
rz(-1.4942687) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0594306) q[0];
sx q[0];
rz(-0.63157394) q[0];
sx q[0];
rz(-2.6078014) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5266085) q[2];
sx q[2];
rz(-0.95643294) q[2];
sx q[2];
rz(0.0690661) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.68077129) q[1];
sx q[1];
rz(-2.6893565) q[1];
sx q[1];
rz(2.3994563) q[1];
rz(-pi) q[2];
rz(1.5288197) q[3];
sx q[3];
rz(-1.4577565) q[3];
sx q[3];
rz(2.6338995) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.5668874) q[2];
sx q[2];
rz(-3.0988099) q[2];
sx q[2];
rz(0.031878397) q[2];
rz(2.3632862) q[3];
sx q[3];
rz(-0.0068155546) q[3];
sx q[3];
rz(0.2955029) q[3];
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
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42302172) q[0];
sx q[0];
rz(-1.5324677) q[0];
sx q[0];
rz(1.8146429) q[0];
rz(3.0146535) q[1];
sx q[1];
rz(-2.9025684) q[1];
sx q[1];
rz(-2.92166) q[1];
rz(1.4311287) q[2];
sx q[2];
rz(-1.5652547) q[2];
sx q[2];
rz(1.8143285) q[2];
rz(-1.4999785) q[3];
sx q[3];
rz(-0.67588617) q[3];
sx q[3];
rz(-2.5839154) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
