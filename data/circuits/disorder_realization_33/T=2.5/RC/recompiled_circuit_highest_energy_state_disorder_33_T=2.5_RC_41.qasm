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
rz(-0.76264277) q[1];
sx q[1];
rz(-1.5781382) q[1];
sx q[1];
rz(-3.0866403) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
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
rz(-3.0829551) q[2];
sx q[2];
rz(-1.8203041) q[2];
sx q[2];
rz(-2.3326258) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.82970769) q[1];
sx q[1];
rz(-1.6685608) q[1];
sx q[1];
rz(-0.13297653) q[1];
rz(-0.94934978) q[3];
sx q[3];
rz(-1.8481405) q[3];
sx q[3];
rz(2.8327033) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.83076465) q[2];
sx q[2];
rz(-2.1711633) q[2];
sx q[2];
rz(0.27466276) q[2];
rz(0.87485391) q[3];
sx q[3];
rz(-2.0503876) q[3];
sx q[3];
rz(-0.27354512) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65983588) q[0];
sx q[0];
rz(-0.6898703) q[0];
sx q[0];
rz(2.7428395) q[0];
rz(0.2440456) q[1];
sx q[1];
rz(-1.9052541) q[1];
sx q[1];
rz(-1.1057378) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2139521) q[0];
sx q[0];
rz(-1.1301148) q[0];
sx q[0];
rz(1.3182314) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.87200882) q[2];
sx q[2];
rz(-0.78561831) q[2];
sx q[2];
rz(-2.9532022) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.2011021) q[1];
sx q[1];
rz(-3.1374212) q[1];
sx q[1];
rz(0.68912403) q[1];
rz(-0.12064528) q[3];
sx q[3];
rz(-1.835726) q[3];
sx q[3];
rz(0.93075965) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.70046052) q[2];
sx q[2];
rz(-1.1819785) q[2];
sx q[2];
rz(1.0750809) q[2];
rz(-2.4454146) q[3];
sx q[3];
rz(-1.8680365) q[3];
sx q[3];
rz(-2.9785494) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8657846) q[0];
sx q[0];
rz(-1.8901261) q[0];
sx q[0];
rz(-2.9248917) q[0];
rz(1.7614583) q[1];
sx q[1];
rz(-0.41788995) q[1];
sx q[1];
rz(2.5221141) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3724646) q[0];
sx q[0];
rz(-1.6004246) q[0];
sx q[0];
rz(1.6086786) q[0];
rz(-pi) q[1];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.9668321) q[1];
sx q[1];
rz(-1.7836387) q[1];
sx q[1];
rz(-1.2618466) q[1];
rz(2.3893395) q[3];
sx q[3];
rz(-2.3779388) q[3];
sx q[3];
rz(-1.4925166) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.2565903) q[2];
sx q[2];
rz(-1.2991178) q[2];
sx q[2];
rz(-0.090506434) q[2];
rz(0.45804405) q[3];
sx q[3];
rz(-0.48726714) q[3];
sx q[3];
rz(1.1080144) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8103545) q[0];
sx q[0];
rz(-2.2585223) q[0];
sx q[0];
rz(-2.2099387) q[0];
rz(-2.9725507) q[1];
sx q[1];
rz(-2.5773498) q[1];
sx q[1];
rz(2.1021252) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0213892) q[0];
sx q[0];
rz(-1.3303555) q[0];
sx q[0];
rz(-0.73657764) q[0];
rz(-pi) q[1];
rz(2.2957689) q[2];
sx q[2];
rz(-2.8081144) q[2];
sx q[2];
rz(-0.84104702) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.89371496) q[1];
sx q[1];
rz(-1.3637241) q[1];
sx q[1];
rz(-3.0805001) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.63757293) q[3];
sx q[3];
rz(-0.3488144) q[3];
sx q[3];
rz(-1.82774) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.6197551) q[2];
sx q[2];
rz(-1.1014742) q[2];
sx q[2];
rz(-0.49631611) q[2];
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
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26688823) q[0];
sx q[0];
rz(-2.5570091) q[0];
sx q[0];
rz(0.97187483) q[0];
rz(0.12174363) q[1];
sx q[1];
rz(-1.5309445) q[1];
sx q[1];
rz(1.8121388) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9864101) q[0];
sx q[0];
rz(-1.2647795) q[0];
sx q[0];
rz(-0.06432342) q[0];
x q[1];
rz(-2.2586063) q[2];
sx q[2];
rz(-2.5341883) q[2];
sx q[2];
rz(2.4508053) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.5426104) q[1];
sx q[1];
rz(-2.3282781) q[1];
sx q[1];
rz(-2.5573362) q[1];
rz(1.4917144) q[3];
sx q[3];
rz(-1.5800161) q[3];
sx q[3];
rz(-0.10445933) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-3.04008) q[2];
sx q[2];
rz(-1.9209346) q[2];
sx q[2];
rz(-2.9902048) q[2];
rz(0.76465145) q[3];
sx q[3];
rz(-2.305856) q[3];
sx q[3];
rz(1.5966655) q[3];
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
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0720035) q[0];
sx q[0];
rz(-1.4729426) q[0];
sx q[0];
rz(-2.0701011) q[0];
rz(-2.6103861) q[1];
sx q[1];
rz(-1.6516282) q[1];
sx q[1];
rz(-2.5154617) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.75979489) q[0];
sx q[0];
rz(-3.1239428) q[0];
sx q[0];
rz(1.2977029) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5444736) q[2];
sx q[2];
rz(-1.7749987) q[2];
sx q[2];
rz(-2.0300421) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.31589139) q[1];
sx q[1];
rz(-0.69676149) q[1];
sx q[1];
rz(0.48721643) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.49895309) q[3];
sx q[3];
rz(-1.1756983) q[3];
sx q[3];
rz(-0.56624352) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.35817394) q[2];
sx q[2];
rz(-0.61966115) q[2];
sx q[2];
rz(-1.0142856) q[2];
rz(-1.2554393) q[3];
sx q[3];
rz(-1.3557326) q[3];
sx q[3];
rz(1.0534508) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1805873) q[0];
sx q[0];
rz(-2.2659681) q[0];
sx q[0];
rz(2.2210806) q[0];
rz(-3.0275184) q[1];
sx q[1];
rz(-1.4055077) q[1];
sx q[1];
rz(-1.1309518) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3944954) q[0];
sx q[0];
rz(-2.1685026) q[0];
sx q[0];
rz(1.5639485) q[0];
x q[1];
rz(-0.99920706) q[2];
sx q[2];
rz(-2.0519369) q[2];
sx q[2];
rz(2.4740117) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.1753833) q[1];
sx q[1];
rz(-1.7031324) q[1];
sx q[1];
rz(-2.2369409) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6507963) q[3];
sx q[3];
rz(-1.0558075) q[3];
sx q[3];
rz(-1.966371) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.34746927) q[2];
sx q[2];
rz(-2.4994714) q[2];
sx q[2];
rz(-2.7327909) q[2];
rz(0.18320228) q[3];
sx q[3];
rz(-1.5902404) q[3];
sx q[3];
rz(-2.0028152) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0278397) q[0];
sx q[0];
rz(-0.04700679) q[0];
sx q[0];
rz(0.29079944) q[0];
rz(-0.94789061) q[1];
sx q[1];
rz(-0.46698505) q[1];
sx q[1];
rz(-1.9416169) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9010097) q[0];
sx q[0];
rz(-1.9685317) q[0];
sx q[0];
rz(-0.57698864) q[0];
rz(-pi) q[1];
rz(0.80757226) q[2];
sx q[2];
rz(-1.0243197) q[2];
sx q[2];
rz(-1.228491) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.7941798) q[1];
sx q[1];
rz(-2.844226) q[1];
sx q[1];
rz(1.411924) q[1];
rz(-pi) q[2];
x q[2];
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
rz(-0.36134186) q[2];
sx q[2];
rz(-1.6929408) q[2];
sx q[2];
rz(-1.7928436) q[2];
rz(1.6382943) q[3];
sx q[3];
rz(-2.2595451) q[3];
sx q[3];
rz(0.1851113) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0272442) q[0];
sx q[0];
rz(-0.37186563) q[0];
sx q[0];
rz(1.0741023) q[0];
rz(-2.7117924) q[1];
sx q[1];
rz(-1.0440412) q[1];
sx q[1];
rz(-0.46357402) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6493452) q[0];
sx q[0];
rz(-1.4758291) q[0];
sx q[0];
rz(-1.3819206) q[0];
x q[1];
rz(2.8797382) q[2];
sx q[2];
rz(-1.3758278) q[2];
sx q[2];
rz(-1.5443813) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.2598272) q[1];
sx q[1];
rz(-0.9503839) q[1];
sx q[1];
rz(2.6735191) q[1];
rz(-2.5829599) q[3];
sx q[3];
rz(-1.6470634) q[3];
sx q[3];
rz(0.35364756) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.27434906) q[2];
sx q[2];
rz(-0.63259071) q[2];
sx q[2];
rz(0.80844936) q[2];
rz(-0.48458734) q[3];
sx q[3];
rz(-0.37823585) q[3];
sx q[3];
rz(2.7073879) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3633858) q[0];
sx q[0];
rz(-2.6860542) q[0];
sx q[0];
rz(-2.0090012) q[0];
rz(0.038381902) q[1];
sx q[1];
rz(-1.3382341) q[1];
sx q[1];
rz(2.8448232) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7420121) q[0];
sx q[0];
rz(-2.098408) q[0];
sx q[0];
rz(0.3216089) q[0];
rz(-pi) q[1];
rz(2.0149258) q[2];
sx q[2];
rz(-1.9812366) q[2];
sx q[2];
rz(0.12516147) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.0368288) q[1];
sx q[1];
rz(-1.4955758) q[1];
sx q[1];
rz(-0.22237088) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.39408306) q[3];
sx q[3];
rz(-2.218045) q[3];
sx q[3];
rz(0.18612063) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.177629) q[2];
sx q[2];
rz(-1.12135) q[2];
sx q[2];
rz(0.56619823) q[2];
rz(2.0171793) q[3];
sx q[3];
rz(-3.0163613) q[3];
sx q[3];
rz(-1.9521149) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2531256) q[0];
sx q[0];
rz(-0.69857004) q[0];
sx q[0];
rz(0.10733124) q[0];
rz(-2.879907) q[1];
sx q[1];
rz(-1.2005922) q[1];
sx q[1];
rz(-2.190879) q[1];
rz(2.8520201) q[2];
sx q[2];
rz(-3.0560383) q[2];
sx q[2];
rz(-0.4734931) q[2];
rz(0.86033173) q[3];
sx q[3];
rz(-2.0900149) q[3];
sx q[3];
rz(2.7505977) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
