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
rz(-2.9774732) q[0];
sx q[0];
rz(-1.0241221) q[0];
sx q[0];
rz(2.6574988) q[0];
rz(2.3789499) q[1];
sx q[1];
rz(4.7197309) q[1];
sx q[1];
rz(9.3698256) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8380175) q[0];
sx q[0];
rz(-2.2406881) q[0];
sx q[0];
rz(-3.1017257) q[0];
x q[1];
rz(-0.058637549) q[2];
sx q[2];
rz(-1.3212886) q[2];
sx q[2];
rz(0.80896689) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.311885) q[1];
sx q[1];
rz(-1.4730318) q[1];
sx q[1];
rz(-0.13297653) q[1];
rz(-1.1160158) q[3];
sx q[3];
rz(-2.4686128) q[3];
sx q[3];
rz(-2.2448886) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.310828) q[2];
sx q[2];
rz(-2.1711633) q[2];
sx q[2];
rz(0.27466276) q[2];
rz(0.87485391) q[3];
sx q[3];
rz(-2.0503876) q[3];
sx q[3];
rz(2.8680475) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65983588) q[0];
sx q[0];
rz(-0.6898703) q[0];
sx q[0];
rz(-2.7428395) q[0];
rz(2.8975471) q[1];
sx q[1];
rz(-1.9052541) q[1];
sx q[1];
rz(-2.0358548) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2139521) q[0];
sx q[0];
rz(-2.0114779) q[0];
sx q[0];
rz(1.8233612) q[0];
rz(-pi) q[1];
rz(-2.5697487) q[2];
sx q[2];
rz(-2.1430121) q[2];
sx q[2];
rz(2.4586611) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.2011021) q[1];
sx q[1];
rz(-0.0041714287) q[1];
sx q[1];
rz(-0.68912403) q[1];
rz(-pi) q[2];
rz(3.0209474) q[3];
sx q[3];
rz(-1.3058666) q[3];
sx q[3];
rz(2.210833) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.70046052) q[2];
sx q[2];
rz(-1.9596142) q[2];
sx q[2];
rz(-2.0665118) q[2];
rz(2.4454146) q[3];
sx q[3];
rz(-1.8680365) q[3];
sx q[3];
rz(2.9785494) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8657846) q[0];
sx q[0];
rz(-1.2514665) q[0];
sx q[0];
rz(2.9248917) q[0];
rz(-1.7614583) q[1];
sx q[1];
rz(-0.41788995) q[1];
sx q[1];
rz(0.61947852) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.769128) q[0];
sx q[0];
rz(-1.541168) q[0];
sx q[0];
rz(-1.6086786) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2080433) q[2];
sx q[2];
rz(-0.13158509) q[2];
sx q[2];
rz(-0.059748273) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.463355) q[1];
sx q[1];
rz(-1.872552) q[1];
sx q[1];
rz(-2.9185118) q[1];
rz(2.3893395) q[3];
sx q[3];
rz(-2.3779388) q[3];
sx q[3];
rz(1.6490761) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.8850024) q[2];
sx q[2];
rz(-1.2991178) q[2];
sx q[2];
rz(-3.0510862) q[2];
rz(-2.6835486) q[3];
sx q[3];
rz(-0.48726714) q[3];
sx q[3];
rz(-2.0335782) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8103545) q[0];
sx q[0];
rz(-0.88307035) q[0];
sx q[0];
rz(-2.2099387) q[0];
rz(2.9725507) q[1];
sx q[1];
rz(-0.5642429) q[1];
sx q[1];
rz(2.1021252) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66327205) q[0];
sx q[0];
rz(-0.86000681) q[0];
sx q[0];
rz(1.8904314) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9157964) q[2];
sx q[2];
rz(-1.3232627) q[2];
sx q[2];
rz(-3.053726) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.1828945) q[1];
sx q[1];
rz(-2.9258203) q[1];
sx q[1];
rz(1.853626) q[1];
rz(-pi) q[2];
rz(1.7839892) q[3];
sx q[3];
rz(-1.8490095) q[3];
sx q[3];
rz(-1.9813862) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.52183759) q[2];
sx q[2];
rz(-1.1014742) q[2];
sx q[2];
rz(0.49631611) q[2];
rz(2.3450092) q[3];
sx q[3];
rz(-2.8556672) q[3];
sx q[3];
rz(-2.0485785) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8747044) q[0];
sx q[0];
rz(-2.5570091) q[0];
sx q[0];
rz(-2.1697178) q[0];
rz(-0.12174363) q[1];
sx q[1];
rz(-1.6106482) q[1];
sx q[1];
rz(1.8121388) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.055453528) q[0];
sx q[0];
rz(-2.8290966) q[0];
sx q[0];
rz(1.3700831) q[0];
rz(-pi) q[1];
rz(-2.2586063) q[2];
sx q[2];
rz(-0.60740439) q[2];
sx q[2];
rz(0.6907874) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.36523) q[1];
sx q[1];
rz(-0.91971469) q[1];
sx q[1];
rz(2.0988223) q[1];
rz(-pi) q[2];
rz(0.0092486898) q[3];
sx q[3];
rz(-1.4917177) q[3];
sx q[3];
rz(1.6759863) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-3.04008) q[2];
sx q[2];
rz(-1.2206581) q[2];
sx q[2];
rz(2.9902048) q[2];
rz(2.3769412) q[3];
sx q[3];
rz(-0.83573666) q[3];
sx q[3];
rz(-1.5449272) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0720035) q[0];
sx q[0];
rz(-1.66865) q[0];
sx q[0];
rz(-1.0714916) q[0];
rz(-2.6103861) q[1];
sx q[1];
rz(-1.4899645) q[1];
sx q[1];
rz(-0.62613097) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.75979489) q[0];
sx q[0];
rz(-0.017649895) q[0];
sx q[0];
rz(1.8438898) q[0];
rz(-2.5444736) q[2];
sx q[2];
rz(-1.3665939) q[2];
sx q[2];
rz(2.0300421) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.6408566) q[1];
sx q[1];
rz(-1.2656414) q[1];
sx q[1];
rz(2.5049097) q[1];
x q[2];
rz(0.71685426) q[3];
sx q[3];
rz(-0.62590137) q[3];
sx q[3];
rz(2.7519873) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7834187) q[2];
sx q[2];
rz(-2.5219315) q[2];
sx q[2];
rz(1.0142856) q[2];
rz(-1.8861534) q[3];
sx q[3];
rz(-1.7858601) q[3];
sx q[3];
rz(1.0534508) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(2.1805873) q[0];
sx q[0];
rz(-2.2659681) q[0];
sx q[0];
rz(0.92051202) q[0];
rz(-0.11407425) q[1];
sx q[1];
rz(-1.4055077) q[1];
sx q[1];
rz(1.1309518) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3823272) q[0];
sx q[0];
rz(-0.59774071) q[0];
sx q[0];
rz(3.1315342) q[0];
rz(-2.1423856) q[2];
sx q[2];
rz(-1.0896557) q[2];
sx q[2];
rz(2.4740117) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.9662094) q[1];
sx q[1];
rz(-1.7031324) q[1];
sx q[1];
rz(0.90465178) q[1];
rz(-pi) q[2];
rz(-1.4907964) q[3];
sx q[3];
rz(-2.0857852) q[3];
sx q[3];
rz(1.1752216) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.34746927) q[2];
sx q[2];
rz(-0.64212126) q[2];
sx q[2];
rz(0.40880173) q[2];
rz(0.18320228) q[3];
sx q[3];
rz(-1.5902404) q[3];
sx q[3];
rz(-2.0028152) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.113753) q[0];
sx q[0];
rz(-0.04700679) q[0];
sx q[0];
rz(2.8507932) q[0];
rz(2.193702) q[1];
sx q[1];
rz(-2.6746076) q[1];
sx q[1];
rz(1.9416169) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24058293) q[0];
sx q[0];
rz(-1.173061) q[0];
sx q[0];
rz(0.57698864) q[0];
x q[1];
rz(2.3340204) q[2];
sx q[2];
rz(-2.117273) q[2];
sx q[2];
rz(1.9131017) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.0702182) q[1];
sx q[1];
rz(-1.6171675) q[1];
sx q[1];
rz(1.2769615) q[1];
rz(-pi) q[2];
rz(2.070571) q[3];
sx q[3];
rz(-2.933549) q[3];
sx q[3];
rz(0.87394729) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.36134186) q[2];
sx q[2];
rz(-1.6929408) q[2];
sx q[2];
rz(1.7928436) q[2];
rz(1.6382943) q[3];
sx q[3];
rz(-0.88204757) q[3];
sx q[3];
rz(2.9564814) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
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
rz(-0.46357402) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.060424711) q[0];
sx q[0];
rz(-1.3827818) q[0];
sx q[0];
rz(0.096676143) q[0];
rz(-pi) q[1];
rz(2.4899269) q[2];
sx q[2];
rz(-0.32512384) q[2];
sx q[2];
rz(-0.59949694) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.9755755) q[1];
sx q[1];
rz(-2.3834627) q[1];
sx q[1];
rz(-1.0075955) q[1];
x q[2];
rz(0.55863278) q[3];
sx q[3];
rz(-1.6470634) q[3];
sx q[3];
rz(0.35364756) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.8672436) q[2];
sx q[2];
rz(-0.63259071) q[2];
sx q[2];
rz(2.3331433) q[2];
rz(2.6570053) q[3];
sx q[3];
rz(-2.7633568) q[3];
sx q[3];
rz(-2.7073879) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
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
rz(-1.3633858) q[0];
sx q[0];
rz(-2.6860542) q[0];
sx q[0];
rz(2.0090012) q[0];
rz(3.1032108) q[1];
sx q[1];
rz(-1.3382341) q[1];
sx q[1];
rz(-2.8448232) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39958056) q[0];
sx q[0];
rz(-2.098408) q[0];
sx q[0];
rz(0.3216089) q[0];
rz(-pi) q[1];
x q[1];
rz(0.77905853) q[2];
sx q[2];
rz(-0.59528661) q[2];
sx q[2];
rz(-2.1434458) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.78690517) q[1];
sx q[1];
rz(-2.9070435) q[1];
sx q[1];
rz(0.3292747) q[1];
rz(-pi) q[2];
rz(2.2567883) q[3];
sx q[3];
rz(-1.8821041) q[3];
sx q[3];
rz(1.5112359) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.177629) q[2];
sx q[2];
rz(-1.12135) q[2];
sx q[2];
rz(-0.56619823) q[2];
rz(2.0171793) q[3];
sx q[3];
rz(-3.0163613) q[3];
sx q[3];
rz(-1.9521149) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.88846702) q[0];
sx q[0];
rz(-2.4430226) q[0];
sx q[0];
rz(-3.0342614) q[0];
rz(-2.879907) q[1];
sx q[1];
rz(-1.2005922) q[1];
sx q[1];
rz(-2.190879) q[1];
rz(-1.546312) q[2];
sx q[2];
rz(-1.6527805) q[2];
sx q[2];
rz(2.3775227) q[2];
rz(-0.85121831) q[3];
sx q[3];
rz(-0.85243445) q[3];
sx q[3];
rz(1.7029521) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
