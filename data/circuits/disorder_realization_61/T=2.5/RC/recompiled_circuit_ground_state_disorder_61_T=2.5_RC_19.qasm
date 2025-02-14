OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.9899848) q[0];
sx q[0];
rz(4.0816981) q[0];
sx q[0];
rz(8.8844086) q[0];
rz(-2.6481533) q[1];
sx q[1];
rz(-0.72055888) q[1];
sx q[1];
rz(-2.5277353) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2342398) q[0];
sx q[0];
rz(-2.1932903) q[0];
sx q[0];
rz(1.7982152) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5905321) q[2];
sx q[2];
rz(-2.8497549) q[2];
sx q[2];
rz(0.29668754) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.39345523) q[1];
sx q[1];
rz(-2.2917213) q[1];
sx q[1];
rz(-1.1925405) q[1];
rz(3.0560232) q[3];
sx q[3];
rz(-0.24805476) q[3];
sx q[3];
rz(2.5906467) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.2962239) q[2];
sx q[2];
rz(-1.2426528) q[2];
sx q[2];
rz(-2.5073012) q[2];
rz(0.93572179) q[3];
sx q[3];
rz(-2.8959385) q[3];
sx q[3];
rz(-2.3989357) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3212386) q[0];
sx q[0];
rz(-1.6004434) q[0];
sx q[0];
rz(2.6974005) q[0];
rz(-2.3815637) q[1];
sx q[1];
rz(-2.0030231) q[1];
sx q[1];
rz(0.98145032) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.32641706) q[0];
sx q[0];
rz(-1.278864) q[0];
sx q[0];
rz(0.53262226) q[0];
rz(-pi) q[1];
rz(0.75656105) q[2];
sx q[2];
rz(-2.8943099) q[2];
sx q[2];
rz(0.73723388) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.3722575) q[1];
sx q[1];
rz(-0.44722873) q[1];
sx q[1];
rz(-0.20453899) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.4325474) q[3];
sx q[3];
rz(-1.3290231) q[3];
sx q[3];
rz(1.6530619) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.11357073) q[2];
sx q[2];
rz(-0.61442033) q[2];
sx q[2];
rz(-1.9363972) q[2];
rz(2.6049854) q[3];
sx q[3];
rz(-1.9034932) q[3];
sx q[3];
rz(1.7369778) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9179012) q[0];
sx q[0];
rz(-1.8603928) q[0];
sx q[0];
rz(-0.83876383) q[0];
rz(-0.63610786) q[1];
sx q[1];
rz(-1.6517703) q[1];
sx q[1];
rz(0.020523358) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8021401) q[0];
sx q[0];
rz(-2.1198556) q[0];
sx q[0];
rz(-2.2744176) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0383561) q[2];
sx q[2];
rz(-0.9976495) q[2];
sx q[2];
rz(-1.8351042) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.518259) q[1];
sx q[1];
rz(-1.9837399) q[1];
sx q[1];
rz(-1.5790126) q[1];
rz(-2.9841524) q[3];
sx q[3];
rz(-1.8529112) q[3];
sx q[3];
rz(2.0705786) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.3402349) q[2];
sx q[2];
rz(-1.6131718) q[2];
sx q[2];
rz(-2.197263) q[2];
rz(-2.1186192) q[3];
sx q[3];
rz(-1.2228271) q[3];
sx q[3];
rz(-2.003722) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
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
rz(2.2706547) q[0];
sx q[0];
rz(-1.9671257) q[0];
sx q[0];
rz(-1.9842072) q[0];
rz(1.4986787) q[1];
sx q[1];
rz(-1.4721556) q[1];
sx q[1];
rz(-1.1882943) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0152215) q[0];
sx q[0];
rz(-1.0980716) q[0];
sx q[0];
rz(1.0593828) q[0];
rz(0.073578667) q[2];
sx q[2];
rz(-0.54983222) q[2];
sx q[2];
rz(-2.4940707) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.3664361) q[1];
sx q[1];
rz(-0.85298733) q[1];
sx q[1];
rz(-2.3597673) q[1];
x q[2];
rz(-1.1340302) q[3];
sx q[3];
rz(-2.7705857) q[3];
sx q[3];
rz(-1.8414258) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.0231861) q[2];
sx q[2];
rz(-1.964317) q[2];
sx q[2];
rz(-2.4510621) q[2];
rz(-1.1150507) q[3];
sx q[3];
rz(-0.77459049) q[3];
sx q[3];
rz(1.5007277) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0951776) q[0];
sx q[0];
rz(-1.7061808) q[0];
sx q[0];
rz(-1.0272367) q[0];
rz(-1.3014303) q[1];
sx q[1];
rz(-2.4899028) q[1];
sx q[1];
rz(-2.8177736) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1808609) q[0];
sx q[0];
rz(-1.9692076) q[0];
sx q[0];
rz(0.84503998) q[0];
rz(1.0179881) q[2];
sx q[2];
rz(-1.5587285) q[2];
sx q[2];
rz(0.34876212) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.65460881) q[1];
sx q[1];
rz(-1.8909847) q[1];
sx q[1];
rz(2.5694808) q[1];
rz(-0.15669723) q[3];
sx q[3];
rz(-2.0771871) q[3];
sx q[3];
rz(1.8739669) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.7258437) q[2];
sx q[2];
rz(-2.9314633) q[2];
sx q[2];
rz(-2.6591163) q[2];
rz(-1.9735362) q[3];
sx q[3];
rz(-1.9246293) q[3];
sx q[3];
rz(2.4647958) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2060858) q[0];
sx q[0];
rz(-0.85010234) q[0];
sx q[0];
rz(0.8859984) q[0];
rz(0.63367263) q[1];
sx q[1];
rz(-2.251667) q[1];
sx q[1];
rz(-0.69127965) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6768178) q[0];
sx q[0];
rz(-0.64660836) q[0];
sx q[0];
rz(3.0809513) q[0];
rz(-pi) q[1];
rz(-0.69181594) q[2];
sx q[2];
rz(-1.5216646) q[2];
sx q[2];
rz(0.092158801) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-3.1289802) q[1];
sx q[1];
rz(-0.75769934) q[1];
sx q[1];
rz(1.4441667) q[1];
x q[2];
rz(-0.25005682) q[3];
sx q[3];
rz(-2.3365006) q[3];
sx q[3];
rz(0.4522194) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.0014235) q[2];
sx q[2];
rz(-2.3052577) q[2];
sx q[2];
rz(-0.26958618) q[2];
rz(0.94830281) q[3];
sx q[3];
rz(-1.6092665) q[3];
sx q[3];
rz(-1.74291) q[3];
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
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7798994) q[0];
sx q[0];
rz(-2.7044856) q[0];
sx q[0];
rz(-2.1345188) q[0];
rz(-2.7359447) q[1];
sx q[1];
rz(-2.5463153) q[1];
sx q[1];
rz(-0.55353037) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6104975) q[0];
sx q[0];
rz(-0.26108867) q[0];
sx q[0];
rz(1.7946662) q[0];
x q[1];
rz(2.5895025) q[2];
sx q[2];
rz(-2.0497397) q[2];
sx q[2];
rz(-0.29651422) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.036052536) q[1];
sx q[1];
rz(-0.80784384) q[1];
sx q[1];
rz(-0.64942645) q[1];
rz(-0.062640142) q[3];
sx q[3];
rz(-2.7609112) q[3];
sx q[3];
rz(2.9306987) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.3125399) q[2];
sx q[2];
rz(-1.092814) q[2];
sx q[2];
rz(-1.2139758) q[2];
rz(2.35516) q[3];
sx q[3];
rz(-1.7460456) q[3];
sx q[3];
rz(3.1184375) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.19602747) q[0];
sx q[0];
rz(-2.8493311) q[0];
sx q[0];
rz(1.3319525) q[0];
rz(2.5281483) q[1];
sx q[1];
rz(-2.1211801) q[1];
sx q[1];
rz(2.6920998) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4411104) q[0];
sx q[0];
rz(-1.4876517) q[0];
sx q[0];
rz(-2.0594199) q[0];
rz(-pi) q[1];
rz(1.23486) q[2];
sx q[2];
rz(-1.4892231) q[2];
sx q[2];
rz(-2.5334266) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.39169381) q[1];
sx q[1];
rz(-2.3428681) q[1];
sx q[1];
rz(2.9500089) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9233016) q[3];
sx q[3];
rz(-2.8779753) q[3];
sx q[3];
rz(-1.6363615) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.8248262) q[2];
sx q[2];
rz(-0.68244857) q[2];
sx q[2];
rz(0.60834926) q[2];
rz(-1.9781205) q[3];
sx q[3];
rz(-1.7721662) q[3];
sx q[3];
rz(2.5782862) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1083531) q[0];
sx q[0];
rz(-2.2305363) q[0];
sx q[0];
rz(0.60229993) q[0];
rz(-2.1741518) q[1];
sx q[1];
rz(-0.93092218) q[1];
sx q[1];
rz(2.0379351) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8062144) q[0];
sx q[0];
rz(-2.4455482) q[0];
sx q[0];
rz(0.93675128) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5561364) q[2];
sx q[2];
rz(-2.2073064) q[2];
sx q[2];
rz(2.9755862) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.97776088) q[1];
sx q[1];
rz(-1.1196616) q[1];
sx q[1];
rz(0.0013063858) q[1];
rz(-pi) q[2];
rz(2.1488701) q[3];
sx q[3];
rz(-1.6274656) q[3];
sx q[3];
rz(-0.68777675) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.2063107) q[2];
sx q[2];
rz(-1.3056359) q[2];
sx q[2];
rz(-0.3375816) q[2];
rz(0.28389367) q[3];
sx q[3];
rz(-0.68113911) q[3];
sx q[3];
rz(-2.9547227) q[3];
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
x q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5481446) q[0];
sx q[0];
rz(-1.633506) q[0];
sx q[0];
rz(3.0737851) q[0];
rz(-1.1495122) q[1];
sx q[1];
rz(-1.5510473) q[1];
sx q[1];
rz(-0.4745208) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7260925) q[0];
sx q[0];
rz(-1.584238) q[0];
sx q[0];
rz(1.8095762) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7065918) q[2];
sx q[2];
rz(-1.8095922) q[2];
sx q[2];
rz(-0.39993024) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.89868916) q[1];
sx q[1];
rz(-1.6575282) q[1];
sx q[1];
rz(-0.0068706339) q[1];
rz(-pi) q[2];
x q[2];
rz(0.81409295) q[3];
sx q[3];
rz(-2.11907) q[3];
sx q[3];
rz(-2.2473638) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.48156753) q[2];
sx q[2];
rz(-2.2109172) q[2];
sx q[2];
rz(-1.0732667) q[2];
rz(0.16452161) q[3];
sx q[3];
rz(-1.8142895) q[3];
sx q[3];
rz(-2.8209414) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
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
rz(0.58707033) q[0];
sx q[0];
rz(-1.5656492) q[0];
sx q[0];
rz(1.5026305) q[0];
rz(-1.5837689) q[1];
sx q[1];
rz(-2.0638034) q[1];
sx q[1];
rz(2.6149909) q[1];
rz(2.4921992) q[2];
sx q[2];
rz(-0.55772256) q[2];
sx q[2];
rz(-1.100308) q[2];
rz(-1.4889553) q[3];
sx q[3];
rz(-1.428953) q[3];
sx q[3];
rz(-0.088464213) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
