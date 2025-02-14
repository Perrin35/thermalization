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
rz(0.64463717) q[0];
sx q[0];
rz(-2.5519389) q[0];
sx q[0];
rz(1.0876422) q[0];
rz(-0.015406869) q[1];
sx q[1];
rz(3.5313731) q[1];
sx q[1];
rz(11.402147) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.919214) q[0];
sx q[0];
rz(-1.5660677) q[0];
sx q[0];
rz(-1.7008002) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9109601) q[2];
sx q[2];
rz(-0.90848604) q[2];
sx q[2];
rz(1.2746122) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.1401745) q[1];
sx q[1];
rz(-1.3571897) q[1];
sx q[1];
rz(1.0402388) q[1];
rz(-pi) q[2];
rz(-1.6887929) q[3];
sx q[3];
rz(-1.6462277) q[3];
sx q[3];
rz(0.82619595) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.2414134) q[2];
sx q[2];
rz(-2.5799077) q[2];
sx q[2];
rz(0.64603311) q[2];
rz(2.8863886) q[3];
sx q[3];
rz(-0.45707688) q[3];
sx q[3];
rz(-1.3946165) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0274886) q[0];
sx q[0];
rz(-0.33559594) q[0];
sx q[0];
rz(-2.3345729) q[0];
rz(-1.457816) q[1];
sx q[1];
rz(-0.57676637) q[1];
sx q[1];
rz(1.3438276) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3407077) q[0];
sx q[0];
rz(-1.721764) q[0];
sx q[0];
rz(-2.8503391) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0033002) q[2];
sx q[2];
rz(-1.9229182) q[2];
sx q[2];
rz(-2.9825236) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.67664424) q[1];
sx q[1];
rz(-1.0847155) q[1];
sx q[1];
rz(2.2653511) q[1];
rz(-pi) q[2];
rz(-0.044780894) q[3];
sx q[3];
rz(-0.84242994) q[3];
sx q[3];
rz(-1.7573259) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.1809711) q[2];
sx q[2];
rz(-1.0701067) q[2];
sx q[2];
rz(-0.19372678) q[2];
rz(-1.6173897) q[3];
sx q[3];
rz(-2.3834855) q[3];
sx q[3];
rz(0.1787506) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5209565) q[0];
sx q[0];
rz(-2.9036324) q[0];
sx q[0];
rz(2.1654907) q[0];
rz(2.2528265) q[1];
sx q[1];
rz(-0.74405324) q[1];
sx q[1];
rz(-2.9685453) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7905549) q[0];
sx q[0];
rz(-1.5332744) q[0];
sx q[0];
rz(-0.0028126082) q[0];
rz(2.6209774) q[2];
sx q[2];
rz(-1.6797425) q[2];
sx q[2];
rz(-1.1586939) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.0806607) q[1];
sx q[1];
rz(-1.0679886) q[1];
sx q[1];
rz(-1.7601556) q[1];
x q[2];
rz(1.6794231) q[3];
sx q[3];
rz(-1.2501688) q[3];
sx q[3];
rz(-2.4214793) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.139107) q[2];
sx q[2];
rz(-1.3287013) q[2];
sx q[2];
rz(-0.71387449) q[2];
rz(0.21126963) q[3];
sx q[3];
rz(-1.0372838) q[3];
sx q[3];
rz(-0.56536388) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7738889) q[0];
sx q[0];
rz(-1.9581032) q[0];
sx q[0];
rz(1.1327889) q[0];
rz(-2.6702113) q[1];
sx q[1];
rz(-1.2939204) q[1];
sx q[1];
rz(2.7105791) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7119766) q[0];
sx q[0];
rz(-1.4578447) q[0];
sx q[0];
rz(-2.9289401) q[0];
rz(-pi) q[1];
rz(2.3296146) q[2];
sx q[2];
rz(-0.50305191) q[2];
sx q[2];
rz(2.4131218) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.9543332) q[1];
sx q[1];
rz(-1.0403087) q[1];
sx q[1];
rz(-1.9217291) q[1];
x q[2];
rz(2.1924344) q[3];
sx q[3];
rz(-1.3484869) q[3];
sx q[3];
rz(0.75321001) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.6212578) q[2];
sx q[2];
rz(-2.4020577) q[2];
sx q[2];
rz(-0.4062824) q[2];
rz(1.2064365) q[3];
sx q[3];
rz(-2.0483569) q[3];
sx q[3];
rz(-0.58667293) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
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
rz(-2.5093812) q[0];
sx q[0];
rz(-2.2659232) q[0];
sx q[0];
rz(-2.8152554) q[0];
rz(2.1874766) q[1];
sx q[1];
rz(-2.4931144) q[1];
sx q[1];
rz(-0.023524806) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1237121) q[0];
sx q[0];
rz(-2.5766226) q[0];
sx q[0];
rz(2.3199757) q[0];
x q[1];
rz(0.65803501) q[2];
sx q[2];
rz(-1.3172704) q[2];
sx q[2];
rz(-0.70632284) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.6885029) q[1];
sx q[1];
rz(-1.9132691) q[1];
sx q[1];
rz(-0.30330412) q[1];
rz(-1.7106423) q[3];
sx q[3];
rz(-1.8636216) q[3];
sx q[3];
rz(1.7738284) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.894459) q[2];
sx q[2];
rz(-0.45866141) q[2];
sx q[2];
rz(-0.66391724) q[2];
rz(0.89020056) q[3];
sx q[3];
rz(-1.592344) q[3];
sx q[3];
rz(-0.012705407) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3506055) q[0];
sx q[0];
rz(-2.7101639) q[0];
sx q[0];
rz(2.8107693) q[0];
rz(-2.779003) q[1];
sx q[1];
rz(-1.9622842) q[1];
sx q[1];
rz(-2.6258452) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2543751) q[0];
sx q[0];
rz(-1.8083113) q[0];
sx q[0];
rz(-0.22254469) q[0];
x q[1];
rz(-0.2771122) q[2];
sx q[2];
rz(-2.5293969) q[2];
sx q[2];
rz(3.1316568) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.679299) q[1];
sx q[1];
rz(-0.88074917) q[1];
sx q[1];
rz(3.0412004) q[1];
x q[2];
rz(-3.0119275) q[3];
sx q[3];
rz(-0.80012459) q[3];
sx q[3];
rz(-2.7050381) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.7911239) q[2];
sx q[2];
rz(-1.6914657) q[2];
sx q[2];
rz(2.3340732) q[2];
rz(3.0430072) q[3];
sx q[3];
rz(-2.2435296) q[3];
sx q[3];
rz(-2.0487823) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.43807855) q[0];
sx q[0];
rz(-2.5683537) q[0];
sx q[0];
rz(-2.692063) q[0];
rz(-0.685177) q[1];
sx q[1];
rz(-0.40760577) q[1];
sx q[1];
rz(2.5221241) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44815608) q[0];
sx q[0];
rz(-1.5991028) q[0];
sx q[0];
rz(2.8319915) q[0];
x q[1];
rz(-1.6190104) q[2];
sx q[2];
rz(-2.0985262) q[2];
sx q[2];
rz(-0.22942782) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.072314315) q[1];
sx q[1];
rz(-0.89767805) q[1];
sx q[1];
rz(2.0050383) q[1];
rz(-pi) q[2];
rz(-2.9569217) q[3];
sx q[3];
rz(-1.369056) q[3];
sx q[3];
rz(-1.4278442) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.5770136) q[2];
sx q[2];
rz(-1.4489633) q[2];
sx q[2];
rz(2.9435834) q[2];
rz(-1.7791344) q[3];
sx q[3];
rz(-0.30006108) q[3];
sx q[3];
rz(-2.4411966) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9971767) q[0];
sx q[0];
rz(-0.62189019) q[0];
sx q[0];
rz(-1.2855592) q[0];
rz(2.8649435) q[1];
sx q[1];
rz(-1.7729019) q[1];
sx q[1];
rz(-0.10861529) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7347057) q[0];
sx q[0];
rz(-0.50841516) q[0];
sx q[0];
rz(-1.7395354) q[0];
rz(-pi) q[1];
rz(-1.4608624) q[2];
sx q[2];
rz(-2.514719) q[2];
sx q[2];
rz(-3.0981321) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-3.120879) q[1];
sx q[1];
rz(-2.2710137) q[1];
sx q[1];
rz(0.28182272) q[1];
x q[2];
rz(-0.1257841) q[3];
sx q[3];
rz(-1.9544807) q[3];
sx q[3];
rz(-0.48641274) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.8331929) q[2];
sx q[2];
rz(-1.1966285) q[2];
sx q[2];
rz(1.2742554) q[2];
rz(-1.7878923) q[3];
sx q[3];
rz(-0.43930587) q[3];
sx q[3];
rz(-2.275009) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(0.15247791) q[0];
sx q[0];
rz(-2.415933) q[0];
sx q[0];
rz(3.0210378) q[0];
rz(2.4001135) q[1];
sx q[1];
rz(-1.9375786) q[1];
sx q[1];
rz(0.075686879) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8814319) q[0];
sx q[0];
rz(-2.7921225) q[0];
sx q[0];
rz(-0.49728877) q[0];
x q[1];
rz(1.7255177) q[2];
sx q[2];
rz(-2.675229) q[2];
sx q[2];
rz(-1.5049962) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.8119252) q[1];
sx q[1];
rz(-0.65663785) q[1];
sx q[1];
rz(-3.1001904) q[1];
rz(2.5406557) q[3];
sx q[3];
rz(-1.4190201) q[3];
sx q[3];
rz(2.5545718) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.072448298) q[2];
sx q[2];
rz(-0.88280237) q[2];
sx q[2];
rz(-0.34633386) q[2];
rz(-2.196178) q[3];
sx q[3];
rz(-1.8190705) q[3];
sx q[3];
rz(0.94714981) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.866975) q[0];
sx q[0];
rz(-1.9434384) q[0];
sx q[0];
rz(0.44788885) q[0];
rz(2.2121494) q[1];
sx q[1];
rz(-1.7421236) q[1];
sx q[1];
rz(-2.4670752) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4607118) q[0];
sx q[0];
rz(-1.7708017) q[0];
sx q[0];
rz(-2.3671211) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.10402502) q[2];
sx q[2];
rz(-1.3859704) q[2];
sx q[2];
rz(1.3074753) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.0793987) q[1];
sx q[1];
rz(-0.17339686) q[1];
sx q[1];
rz(1.5773415) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.075427051) q[3];
sx q[3];
rz(-1.9884681) q[3];
sx q[3];
rz(0.088069629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.931539) q[2];
sx q[2];
rz(-1.209582) q[2];
sx q[2];
rz(-2.8741969) q[2];
rz(-1.1072655) q[3];
sx q[3];
rz(-2.3591154) q[3];
sx q[3];
rz(-0.84899181) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4103107) q[0];
sx q[0];
rz(-1.5189497) q[0];
sx q[0];
rz(2.0547163) q[0];
rz(2.3375753) q[1];
sx q[1];
rz(-1.4878648) q[1];
sx q[1];
rz(1.0555242) q[1];
rz(1.4236535) q[2];
sx q[2];
rz(-1.6416807) q[2];
sx q[2];
rz(0.15105187) q[2];
rz(-0.08392423) q[3];
sx q[3];
rz(-2.2369583) q[3];
sx q[3];
rz(0.53085622) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
