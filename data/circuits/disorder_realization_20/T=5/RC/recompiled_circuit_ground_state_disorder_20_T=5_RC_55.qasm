OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.5372758) q[0];
sx q[0];
rz(-0.24157) q[0];
sx q[0];
rz(0.33302745) q[0];
rz(2.0060519) q[1];
sx q[1];
rz(-0.82692868) q[1];
sx q[1];
rz(-0.64396042) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49227958) q[0];
sx q[0];
rz(-2.4631073) q[0];
sx q[0];
rz(2.0114824) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8415235) q[2];
sx q[2];
rz(-1.7207533) q[2];
sx q[2];
rz(-1.4407002) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.3463979) q[1];
sx q[1];
rz(-1.193622) q[1];
sx q[1];
rz(-1.8257797) q[1];
rz(0.87857492) q[3];
sx q[3];
rz(-2.6248616) q[3];
sx q[3];
rz(2.0661092) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.48646271) q[2];
sx q[2];
rz(-0.4643521) q[2];
sx q[2];
rz(0.74074024) q[2];
rz(-1.5517392) q[3];
sx q[3];
rz(-0.71151763) q[3];
sx q[3];
rz(2.5488502) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
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
rz(-2.7084259) q[0];
sx q[0];
rz(-1.1335224) q[0];
sx q[0];
rz(0.11696996) q[0];
rz(0.50432694) q[1];
sx q[1];
rz(-1.802899) q[1];
sx q[1];
rz(2.8574944) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5586276) q[0];
sx q[0];
rz(-2.5694048) q[0];
sx q[0];
rz(1.4165322) q[0];
rz(-pi) q[1];
rz(-2.5710377) q[2];
sx q[2];
rz(-2.2584256) q[2];
sx q[2];
rz(-0.90435435) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.5631905) q[1];
sx q[1];
rz(-3.0147073) q[1];
sx q[1];
rz(1.8475501) q[1];
rz(-1.8485214) q[3];
sx q[3];
rz(-1.3284995) q[3];
sx q[3];
rz(1.5838069) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.32610193) q[2];
sx q[2];
rz(-2.2559866) q[2];
sx q[2];
rz(-2.3808114) q[2];
rz(-2.271999) q[3];
sx q[3];
rz(-1.875501) q[3];
sx q[3];
rz(-2.6884955) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3608383) q[0];
sx q[0];
rz(-2.0212845) q[0];
sx q[0];
rz(-0.25217062) q[0];
rz(1.4422656) q[1];
sx q[1];
rz(-0.95528722) q[1];
sx q[1];
rz(-2.7853277) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8667135) q[0];
sx q[0];
rz(-1.5769813) q[0];
sx q[0];
rz(1.5112108) q[0];
rz(0.61932694) q[2];
sx q[2];
rz(-2.4973923) q[2];
sx q[2];
rz(-1.8875811) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.85689304) q[1];
sx q[1];
rz(-0.93888043) q[1];
sx q[1];
rz(-0.17669682) q[1];
rz(-pi) q[2];
rz(3.1396418) q[3];
sx q[3];
rz(-1.8271128) q[3];
sx q[3];
rz(-2.0511829) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.1189271) q[2];
sx q[2];
rz(-1.7513195) q[2];
sx q[2];
rz(-1.6290132) q[2];
rz(-3.1318943) q[3];
sx q[3];
rz(-1.9111948) q[3];
sx q[3];
rz(-0.62140083) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(1.6464722) q[0];
sx q[0];
rz(-2.5992114) q[0];
sx q[0];
rz(0.25076184) q[0];
rz(1.3465025) q[1];
sx q[1];
rz(-2.612412) q[1];
sx q[1];
rz(-1.4432602) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0899857) q[0];
sx q[0];
rz(-2.3653226) q[0];
sx q[0];
rz(1.8818186) q[0];
rz(-pi) q[1];
rz(-0.94130959) q[2];
sx q[2];
rz(-2.513859) q[2];
sx q[2];
rz(0.24294397) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.52407284) q[1];
sx q[1];
rz(-2.1012602) q[1];
sx q[1];
rz(2.4071715) q[1];
rz(-pi) q[2];
rz(-1.8399393) q[3];
sx q[3];
rz(-1.3175822) q[3];
sx q[3];
rz(-1.9631752) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.5859588) q[2];
sx q[2];
rz(-1.1052174) q[2];
sx q[2];
rz(0.24331681) q[2];
rz(-0.58468741) q[3];
sx q[3];
rz(-0.57716113) q[3];
sx q[3];
rz(-0.028133597) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.993416) q[0];
sx q[0];
rz(-2.7758444) q[0];
sx q[0];
rz(0.17955968) q[0];
rz(-1.9163632) q[1];
sx q[1];
rz(-1.4045249) q[1];
sx q[1];
rz(-0.45825759) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1711463) q[0];
sx q[0];
rz(-1.866328) q[0];
sx q[0];
rz(1.4895205) q[0];
rz(-0.29693895) q[2];
sx q[2];
rz(-0.7504979) q[2];
sx q[2];
rz(0.82210449) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.51650713) q[1];
sx q[1];
rz(-1.436215) q[1];
sx q[1];
rz(0.075304042) q[1];
rz(-0.36009501) q[3];
sx q[3];
rz(-2.2392139) q[3];
sx q[3];
rz(-0.51262142) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.3852343) q[2];
sx q[2];
rz(-0.42432722) q[2];
sx q[2];
rz(0.9217841) q[2];
rz(-1.2515757) q[3];
sx q[3];
rz(-1.4082963) q[3];
sx q[3];
rz(-3.0799227) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.09403041) q[0];
sx q[0];
rz(-2.229409) q[0];
sx q[0];
rz(-2.9929274) q[0];
rz(-2.8246763) q[1];
sx q[1];
rz(-0.47873679) q[1];
sx q[1];
rz(-2.5247578) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.29051) q[0];
sx q[0];
rz(-1.6491873) q[0];
sx q[0];
rz(-3.0799887) q[0];
rz(-0.89646879) q[2];
sx q[2];
rz(-2.5988262) q[2];
sx q[2];
rz(1.8405869) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.3909797) q[1];
sx q[1];
rz(-1.8209848) q[1];
sx q[1];
rz(0.98586086) q[1];
rz(2.8657593) q[3];
sx q[3];
rz(-2.0778928) q[3];
sx q[3];
rz(-2.3292993) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.2769015) q[2];
sx q[2];
rz(-1.428705) q[2];
sx q[2];
rz(-0.0083262715) q[2];
rz(2.5152123) q[3];
sx q[3];
rz(-2.4255987) q[3];
sx q[3];
rz(3.1410419) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3061227) q[0];
sx q[0];
rz(-1.1518814) q[0];
sx q[0];
rz(-0.48208153) q[0];
rz(0.029190633) q[1];
sx q[1];
rz(-1.8455285) q[1];
sx q[1];
rz(2.3775502) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4493443) q[0];
sx q[0];
rz(-1.6381761) q[0];
sx q[0];
rz(-1.199556) q[0];
rz(0.025310658) q[2];
sx q[2];
rz(-0.9111852) q[2];
sx q[2];
rz(-0.74756223) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.5611539) q[1];
sx q[1];
rz(-0.62404666) q[1];
sx q[1];
rz(-1.8638205) q[1];
rz(3.1336083) q[3];
sx q[3];
rz(-1.1698876) q[3];
sx q[3];
rz(-2.8454091) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.46334106) q[2];
sx q[2];
rz(-1.8222787) q[2];
sx q[2];
rz(2.6574262) q[2];
rz(-0.68228996) q[3];
sx q[3];
rz(-2.4874918) q[3];
sx q[3];
rz(3.0068523) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
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
rz(3.057137) q[0];
sx q[0];
rz(-0.77780044) q[0];
sx q[0];
rz(1.8498259) q[0];
rz(-0.46503398) q[1];
sx q[1];
rz(-2.6227622) q[1];
sx q[1];
rz(-2.8344287) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.43873337) q[0];
sx q[0];
rz(-0.68638681) q[0];
sx q[0];
rz(0.86875654) q[0];
rz(-pi) q[1];
rz(-0.95057733) q[2];
sx q[2];
rz(-2.3831316) q[2];
sx q[2];
rz(1.72067) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.7502082) q[1];
sx q[1];
rz(-2.5473875) q[1];
sx q[1];
rz(-1.2886402) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.028178111) q[3];
sx q[3];
rz(-0.62528505) q[3];
sx q[3];
rz(-2.6373088) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.18053599) q[2];
sx q[2];
rz(-0.44027105) q[2];
sx q[2];
rz(-0.72009909) q[2];
rz(-2.2911206) q[3];
sx q[3];
rz(-1.1668147) q[3];
sx q[3];
rz(-1.4115964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9367323) q[0];
sx q[0];
rz(-0.16848773) q[0];
sx q[0];
rz(2.4334461) q[0];
rz(-1.6234966) q[1];
sx q[1];
rz(-0.93815362) q[1];
sx q[1];
rz(-2.785397) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4771381) q[0];
sx q[0];
rz(-1.9915446) q[0];
sx q[0];
rz(-0.53945213) q[0];
rz(-pi) q[1];
rz(-2.5237066) q[2];
sx q[2];
rz(-1.5466585) q[2];
sx q[2];
rz(-1.9113505) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.5070008) q[1];
sx q[1];
rz(-0.92220491) q[1];
sx q[1];
rz(1.6084987) q[1];
x q[2];
rz(0.40608866) q[3];
sx q[3];
rz(-1.1392987) q[3];
sx q[3];
rz(-0.0270947) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.0922962) q[2];
sx q[2];
rz(-2.3342817) q[2];
sx q[2];
rz(-0.88579196) q[2];
rz(0.93585912) q[3];
sx q[3];
rz(-1.1805308) q[3];
sx q[3];
rz(0.22127557) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7037999) q[0];
sx q[0];
rz(-0.047310345) q[0];
sx q[0];
rz(1.6920775) q[0];
rz(-1.0225147) q[1];
sx q[1];
rz(-0.48373628) q[1];
sx q[1];
rz(-1.6922916) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.81240772) q[0];
sx q[0];
rz(-1.158876) q[0];
sx q[0];
rz(-0.054188577) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3751276) q[2];
sx q[2];
rz(-2.577707) q[2];
sx q[2];
rz(-1.6414798) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.8359747) q[1];
sx q[1];
rz(-0.59594369) q[1];
sx q[1];
rz(-1.9553595) q[1];
rz(-pi) q[2];
rz(-1.3053042) q[3];
sx q[3];
rz(-1.7137626) q[3];
sx q[3];
rz(-2.6022823) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.25541043) q[2];
sx q[2];
rz(-1.1899199) q[2];
sx q[2];
rz(2.4667242) q[2];
rz(2.1700962) q[3];
sx q[3];
rz(-0.70236218) q[3];
sx q[3];
rz(1.491588) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7991199) q[0];
sx q[0];
rz(-1.5457038) q[0];
sx q[0];
rz(-0.85734838) q[0];
rz(-2.9091861) q[1];
sx q[1];
rz(-2.1288165) q[1];
sx q[1];
rz(-1.7850599) q[1];
rz(-2.2131481) q[2];
sx q[2];
rz(-0.63277638) q[2];
sx q[2];
rz(0.28296726) q[2];
rz(-0.078217004) q[3];
sx q[3];
rz(-0.84368869) q[3];
sx q[3];
rz(0.61144184) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
