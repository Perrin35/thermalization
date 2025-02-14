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
rz(4.1657148) q[0];
sx q[0];
rz(9.9088718) q[0];
rz(2.3789499) q[1];
sx q[1];
rz(-1.5634544) q[1];
sx q[1];
rz(-0.054952316) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.30357519) q[0];
sx q[0];
rz(-2.2406881) q[0];
sx q[0];
rz(3.1017257) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.820716) q[2];
sx q[2];
rz(-1.5139765) q[2];
sx q[2];
rz(-2.394258) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.311885) q[1];
sx q[1];
rz(-1.6685608) q[1];
sx q[1];
rz(3.0086161) q[1];
x q[2];
rz(0.94934978) q[3];
sx q[3];
rz(-1.8481405) q[3];
sx q[3];
rz(0.30888939) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.83076465) q[2];
sx q[2];
rz(-0.97042933) q[2];
sx q[2];
rz(-2.8669299) q[2];
rz(2.2667387) q[3];
sx q[3];
rz(-2.0503876) q[3];
sx q[3];
rz(0.27354512) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4817568) q[0];
sx q[0];
rz(-2.4517224) q[0];
sx q[0];
rz(-2.7428395) q[0];
rz(0.2440456) q[1];
sx q[1];
rz(-1.9052541) q[1];
sx q[1];
rz(2.0358548) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.53351346) q[0];
sx q[0];
rz(-1.7987804) q[0];
sx q[0];
rz(-2.6883459) q[0];
x q[1];
rz(0.57184394) q[2];
sx q[2];
rz(-0.99858054) q[2];
sx q[2];
rz(0.68293152) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.9404906) q[1];
sx q[1];
rz(-3.1374212) q[1];
sx q[1];
rz(-0.68912403) q[1];
x q[2];
rz(-3.0209474) q[3];
sx q[3];
rz(-1.3058666) q[3];
sx q[3];
rz(-2.210833) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.4411321) q[2];
sx q[2];
rz(-1.1819785) q[2];
sx q[2];
rz(-2.0665118) q[2];
rz(-2.4454146) q[3];
sx q[3];
rz(-1.2735561) q[3];
sx q[3];
rz(-0.16304326) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8657846) q[0];
sx q[0];
rz(-1.2514665) q[0];
sx q[0];
rz(-2.9248917) q[0];
rz(1.3801344) q[1];
sx q[1];
rz(-2.7237027) q[1];
sx q[1];
rz(-0.61947852) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.769128) q[0];
sx q[0];
rz(-1.6004246) q[0];
sx q[0];
rz(1.6086786) q[0];
x q[1];
rz(1.9335494) q[2];
sx q[2];
rz(-0.13158509) q[2];
sx q[2];
rz(-3.0818444) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.6782376) q[1];
sx q[1];
rz(-1.2690407) q[1];
sx q[1];
rz(-0.22308087) q[1];
rz(-2.5314919) q[3];
sx q[3];
rz(-2.0629632) q[3];
sx q[3];
rz(0.51612332) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.8850024) q[2];
sx q[2];
rz(-1.2991178) q[2];
sx q[2];
rz(0.090506434) q[2];
rz(2.6835486) q[3];
sx q[3];
rz(-0.48726714) q[3];
sx q[3];
rz(-1.1080144) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8103545) q[0];
sx q[0];
rz(-2.2585223) q[0];
sx q[0];
rz(2.2099387) q[0];
rz(-2.9725507) q[1];
sx q[1];
rz(-0.5642429) q[1];
sx q[1];
rz(-2.1021252) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66327205) q[0];
sx q[0];
rz(-0.86000681) q[0];
sx q[0];
rz(1.8904314) q[0];
rz(-0.8458237) q[2];
sx q[2];
rz(-2.8081144) q[2];
sx q[2];
rz(-0.84104702) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.66450602) q[1];
sx q[1];
rz(-1.5110104) q[1];
sx q[1];
rz(1.7782446) q[1];
rz(-pi) q[2];
rz(-2.8572731) q[3];
sx q[3];
rz(-1.7756724) q[3];
sx q[3];
rz(-2.7903872) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.52183759) q[2];
sx q[2];
rz(-1.1014742) q[2];
sx q[2];
rz(0.49631611) q[2];
rz(-2.3450092) q[3];
sx q[3];
rz(-0.28592548) q[3];
sx q[3];
rz(-2.0485785) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8747044) q[0];
sx q[0];
rz(-2.5570091) q[0];
sx q[0];
rz(-0.97187483) q[0];
rz(-3.019849) q[1];
sx q[1];
rz(-1.6106482) q[1];
sx q[1];
rz(-1.8121388) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.055453528) q[0];
sx q[0];
rz(-2.8290966) q[0];
sx q[0];
rz(1.3700831) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2586063) q[2];
sx q[2];
rz(-2.5341883) q[2];
sx q[2];
rz(0.6907874) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.5426104) q[1];
sx q[1];
rz(-2.3282781) q[1];
sx q[1];
rz(-2.5573362) q[1];
rz(-pi) q[2];
rz(1.4546118) q[3];
sx q[3];
rz(-0.079616485) q[3];
sx q[3];
rz(1.5594359) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-3.04008) q[2];
sx q[2];
rz(-1.2206581) q[2];
sx q[2];
rz(-0.15138781) q[2];
rz(-0.76465145) q[3];
sx q[3];
rz(-0.83573666) q[3];
sx q[3];
rz(-1.5449272) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0720035) q[0];
sx q[0];
rz(-1.4729426) q[0];
sx q[0];
rz(1.0714916) q[0];
rz(2.6103861) q[1];
sx q[1];
rz(-1.6516282) q[1];
sx q[1];
rz(2.5154617) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3817978) q[0];
sx q[0];
rz(-0.017649895) q[0];
sx q[0];
rz(-1.8438898) q[0];
x q[1];
rz(1.8161723) q[2];
sx q[2];
rz(-2.1538556) q[2];
sx q[2];
rz(0.32223216) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.5007361) q[1];
sx q[1];
rz(-1.2656414) q[1];
sx q[1];
rz(-2.5049097) q[1];
rz(-pi) q[2];
x q[2];
rz(0.49895309) q[3];
sx q[3];
rz(-1.1756983) q[3];
sx q[3];
rz(-2.5753491) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.7834187) q[2];
sx q[2];
rz(-2.5219315) q[2];
sx q[2];
rz(1.0142856) q[2];
rz(-1.8861534) q[3];
sx q[3];
rz(-1.3557326) q[3];
sx q[3];
rz(-1.0534508) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1805873) q[0];
sx q[0];
rz(-0.87562457) q[0];
sx q[0];
rz(-0.92051202) q[0];
rz(-0.11407425) q[1];
sx q[1];
rz(-1.4055077) q[1];
sx q[1];
rz(-2.0106409) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.74709725) q[0];
sx q[0];
rz(-2.1685026) q[0];
sx q[0];
rz(1.5776442) q[0];
x q[1];
rz(-2.5860687) q[2];
sx q[2];
rz(-2.0709403) q[2];
sx q[2];
rz(1.9490567) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.29207001) q[1];
sx q[1];
rz(-0.91150586) q[1];
sx q[1];
rz(-0.16772049) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4907964) q[3];
sx q[3];
rz(-1.0558075) q[3];
sx q[3];
rz(-1.1752216) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.7941234) q[2];
sx q[2];
rz(-2.4994714) q[2];
sx q[2];
rz(0.40880173) q[2];
rz(2.9583904) q[3];
sx q[3];
rz(-1.5513523) q[3];
sx q[3];
rz(-2.0028152) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.113753) q[0];
sx q[0];
rz(-3.0945859) q[0];
sx q[0];
rz(0.29079944) q[0];
rz(0.94789061) q[1];
sx q[1];
rz(-0.46698505) q[1];
sx q[1];
rz(1.9416169) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5644218) q[0];
sx q[0];
rz(-2.0977712) q[0];
sx q[0];
rz(1.1061263) q[0];
rz(-pi) q[1];
rz(2.3340204) q[2];
sx q[2];
rz(-1.0243197) q[2];
sx q[2];
rz(1.228491) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.7941798) q[1];
sx q[1];
rz(-2.844226) q[1];
sx q[1];
rz(1.411924) q[1];
rz(-pi) q[2];
rz(0.10082106) q[3];
sx q[3];
rz(-1.3885048) q[3];
sx q[3];
rz(1.3829447) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.7802508) q[2];
sx q[2];
rz(-1.6929408) q[2];
sx q[2];
rz(-1.7928436) q[2];
rz(-1.5032984) q[3];
sx q[3];
rz(-2.2595451) q[3];
sx q[3];
rz(0.1851113) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1143484) q[0];
sx q[0];
rz(-2.769727) q[0];
sx q[0];
rz(2.0674904) q[0];
rz(-0.4298003) q[1];
sx q[1];
rz(-2.0975515) q[1];
sx q[1];
rz(2.6780186) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.060424711) q[0];
sx q[0];
rz(-1.3827818) q[0];
sx q[0];
rz(-3.0449165) q[0];
rz(-pi) q[1];
rz(-0.65166574) q[2];
sx q[2];
rz(-0.32512384) q[2];
sx q[2];
rz(-0.59949694) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.2598272) q[1];
sx q[1];
rz(-2.1912088) q[1];
sx q[1];
rz(-0.46807351) q[1];
x q[2];
rz(0.55863278) q[3];
sx q[3];
rz(-1.6470634) q[3];
sx q[3];
rz(-2.7879451) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.27434906) q[2];
sx q[2];
rz(-0.63259071) q[2];
sx q[2];
rz(-0.80844936) q[2];
rz(-0.48458734) q[3];
sx q[3];
rz(-0.37823585) q[3];
sx q[3];
rz(-0.43420473) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3633858) q[0];
sx q[0];
rz(-0.45553842) q[0];
sx q[0];
rz(2.0090012) q[0];
rz(3.1032108) q[1];
sx q[1];
rz(-1.8033586) q[1];
sx q[1];
rz(-0.29676944) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0050186) q[0];
sx q[0];
rz(-1.2941735) q[0];
sx q[0];
rz(-2.121595) q[0];
rz(-1.1266668) q[2];
sx q[2];
rz(-1.160356) q[2];
sx q[2];
rz(3.0164312) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.78690517) q[1];
sx q[1];
rz(-2.9070435) q[1];
sx q[1];
rz(2.812318) q[1];
rz(0.88480437) q[3];
sx q[3];
rz(-1.8821041) q[3];
sx q[3];
rz(1.6303568) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.177629) q[2];
sx q[2];
rz(-1.12135) q[2];
sx q[2];
rz(2.5753944) q[2];
rz(-2.0171793) q[3];
sx q[3];
rz(-0.12523139) q[3];
sx q[3];
rz(-1.9521149) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.88846702) q[0];
sx q[0];
rz(-0.69857004) q[0];
sx q[0];
rz(0.10733124) q[0];
rz(-2.879907) q[1];
sx q[1];
rz(-1.2005922) q[1];
sx q[1];
rz(-2.190879) q[1];
rz(-1.5952806) q[2];
sx q[2];
rz(-1.4888121) q[2];
sx q[2];
rz(-0.76406995) q[2];
rz(2.2903743) q[3];
sx q[3];
rz(-0.85243445) q[3];
sx q[3];
rz(1.7029521) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
