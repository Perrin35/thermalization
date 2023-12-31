OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.47473946) q[0];
sx q[0];
rz(-0.82959509) q[0];
sx q[0];
rz(-2.9876246) q[0];
rz(0.83377588) q[1];
sx q[1];
rz(-2.1492465) q[1];
sx q[1];
rz(-0.33831236) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6688639) q[0];
sx q[0];
rz(-2.7358486) q[0];
sx q[0];
rz(2.2292482) q[0];
rz(-pi) q[1];
x q[1];
rz(0.89262427) q[2];
sx q[2];
rz(-1.8632338) q[2];
sx q[2];
rz(-2.6543648) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(3.0677883) q[1];
sx q[1];
rz(-1.5177625) q[1];
sx q[1];
rz(0.28633134) q[1];
x q[2];
rz(-0.015720856) q[3];
sx q[3];
rz(-1.058488) q[3];
sx q[3];
rz(-0.44954625) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.14264318) q[2];
sx q[2];
rz(-2.8012186) q[2];
sx q[2];
rz(-1.9677229) q[2];
rz(3.0657892) q[3];
sx q[3];
rz(-1.9971763) q[3];
sx q[3];
rz(0.092806667) q[3];
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
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8006111) q[0];
sx q[0];
rz(-2.0759463) q[0];
sx q[0];
rz(-0.064963438) q[0];
rz(-2.5669572) q[1];
sx q[1];
rz(-0.42962933) q[1];
sx q[1];
rz(-1.8992791) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.500538) q[0];
sx q[0];
rz(-0.28028742) q[0];
sx q[0];
rz(2.6602402) q[0];
rz(-2.20349) q[2];
sx q[2];
rz(-1.4422851) q[2];
sx q[2];
rz(3.1046257) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.9638483) q[1];
sx q[1];
rz(-1.7301534) q[1];
sx q[1];
rz(2.6049155) q[1];
rz(-pi) q[2];
x q[2];
rz(0.70988016) q[3];
sx q[3];
rz(-1.1871561) q[3];
sx q[3];
rz(0.82304728) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.3339281) q[2];
sx q[2];
rz(-1.0753205) q[2];
sx q[2];
rz(2.5578965) q[2];
rz(2.5675473) q[3];
sx q[3];
rz(-2.0160926) q[3];
sx q[3];
rz(-0.13124245) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7211001) q[0];
sx q[0];
rz(-2.2476966) q[0];
sx q[0];
rz(2.4131391) q[0];
rz(-1.4942253) q[1];
sx q[1];
rz(-0.39847001) q[1];
sx q[1];
rz(-1.0167936) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25123888) q[0];
sx q[0];
rz(-1.4892502) q[0];
sx q[0];
rz(-1.5348877) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0605893) q[2];
sx q[2];
rz(-2.042633) q[2];
sx q[2];
rz(2.8420198) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.2336725) q[1];
sx q[1];
rz(-0.83041149) q[1];
sx q[1];
rz(-0.18632142) q[1];
x q[2];
rz(-1.8524283) q[3];
sx q[3];
rz(-1.6179807) q[3];
sx q[3];
rz(-0.55164528) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.8016522) q[2];
sx q[2];
rz(-1.6093971) q[2];
sx q[2];
rz(-2.0920848) q[2];
rz(0.63878757) q[3];
sx q[3];
rz(-0.63126606) q[3];
sx q[3];
rz(1.9558186) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3291572) q[0];
sx q[0];
rz(-1.2673459) q[0];
sx q[0];
rz(1.6695492) q[0];
rz(0.73515785) q[1];
sx q[1];
rz(-0.77886326) q[1];
sx q[1];
rz(2.8947815) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4287764) q[0];
sx q[0];
rz(-1.4220337) q[0];
sx q[0];
rz(-2.2674198) q[0];
rz(-pi) q[1];
x q[1];
rz(0.47841448) q[2];
sx q[2];
rz(-1.7583499) q[2];
sx q[2];
rz(2.0699376) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.4775131) q[1];
sx q[1];
rz(-1.941136) q[1];
sx q[1];
rz(-1.5143637) q[1];
rz(-1.1770583) q[3];
sx q[3];
rz(-2.4871832) q[3];
sx q[3];
rz(0.27458336) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.6639158) q[2];
sx q[2];
rz(-2.0229979) q[2];
sx q[2];
rz(-1.6003312) q[2];
rz(-0.70704308) q[3];
sx q[3];
rz(-1.0975081) q[3];
sx q[3];
rz(-2.2533806) q[3];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.33070579) q[0];
sx q[0];
rz(-2.4208477) q[0];
sx q[0];
rz(1.8141618) q[0];
rz(1.5785626) q[1];
sx q[1];
rz(-0.47416082) q[1];
sx q[1];
rz(-0.24838233) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36149704) q[0];
sx q[0];
rz(-2.1108315) q[0];
sx q[0];
rz(1.6105152) q[0];
rz(-pi) q[1];
rz(2.4400473) q[2];
sx q[2];
rz(-2.8994459) q[2];
sx q[2];
rz(-2.0602351) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.0866962) q[1];
sx q[1];
rz(-2.5270562) q[1];
sx q[1];
rz(-0.30026786) q[1];
x q[2];
rz(-2.8304843) q[3];
sx q[3];
rz(-2.4690383) q[3];
sx q[3];
rz(0.5141408) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.7053232) q[2];
sx q[2];
rz(-0.98781172) q[2];
sx q[2];
rz(-0.79745897) q[2];
rz(0.38875368) q[3];
sx q[3];
rz(-2.5377486) q[3];
sx q[3];
rz(0.50271547) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0080863) q[0];
sx q[0];
rz(-3.0651423) q[0];
sx q[0];
rz(-1.7957934) q[0];
rz(-2.0603518) q[1];
sx q[1];
rz(-1.9045647) q[1];
sx q[1];
rz(3.016901) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.161675) q[0];
sx q[0];
rz(-1.3945762) q[0];
sx q[0];
rz(1.6582703) q[0];
rz(-2.5325534) q[2];
sx q[2];
rz(-1.8134724) q[2];
sx q[2];
rz(-2.2720624) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.2248968) q[1];
sx q[1];
rz(-1.5142913) q[1];
sx q[1];
rz(0.7373666) q[1];
x q[2];
rz(0.74495875) q[3];
sx q[3];
rz(-2.8248441) q[3];
sx q[3];
rz(-1.484364) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.6283915) q[2];
sx q[2];
rz(-1.1867563) q[2];
sx q[2];
rz(-0.61895269) q[2];
rz(-1.0533054) q[3];
sx q[3];
rz(-2.9635933) q[3];
sx q[3];
rz(2.2119904) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.58182794) q[0];
sx q[0];
rz(-1.7511837) q[0];
sx q[0];
rz(-1.0429617) q[0];
rz(0.46328059) q[1];
sx q[1];
rz(-1.1136585) q[1];
sx q[1];
rz(-1.0707062) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2830414) q[0];
sx q[0];
rz(-1.0867449) q[0];
sx q[0];
rz(-1.5701243) q[0];
rz(-pi) q[1];
rz(0.35371874) q[2];
sx q[2];
rz(-2.7732447) q[2];
sx q[2];
rz(3.0596717) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-3.0106525) q[1];
sx q[1];
rz(-2.3619235) q[1];
sx q[1];
rz(-0.14203771) q[1];
x q[2];
rz(-0.33220746) q[3];
sx q[3];
rz(-1.5899961) q[3];
sx q[3];
rz(1.7393877) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.36859194) q[2];
sx q[2];
rz(-2.4020782) q[2];
sx q[2];
rz(0.32361844) q[2];
rz(2.1598024) q[3];
sx q[3];
rz(-2.2798645) q[3];
sx q[3];
rz(2.0543082) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48802808) q[0];
sx q[0];
rz(-1.4166778) q[0];
sx q[0];
rz(1.2063684) q[0];
rz(-1.2127097) q[1];
sx q[1];
rz(-2.2884463) q[1];
sx q[1];
rz(0.94747296) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.065814106) q[0];
sx q[0];
rz(-1.6718719) q[0];
sx q[0];
rz(-1.0321192) q[0];
rz(-pi) q[1];
x q[1];
rz(2.410789) q[2];
sx q[2];
rz(-0.23554221) q[2];
sx q[2];
rz(-1.3256324) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.2944813) q[1];
sx q[1];
rz(-2.1335019) q[1];
sx q[1];
rz(0.28114762) q[1];
rz(1.9318337) q[3];
sx q[3];
rz(-0.79663888) q[3];
sx q[3];
rz(-2.3494997) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.56132135) q[2];
sx q[2];
rz(-1.7611971) q[2];
sx q[2];
rz(-0.46978152) q[2];
rz(-1.8404768) q[3];
sx q[3];
rz(-1.7062635) q[3];
sx q[3];
rz(-0.27967134) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8895421) q[0];
sx q[0];
rz(-2.7520576) q[0];
sx q[0];
rz(1.3289733) q[0];
rz(0.7912311) q[1];
sx q[1];
rz(-0.33114094) q[1];
sx q[1];
rz(2.9387617) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.77712599) q[0];
sx q[0];
rz(-1.6098032) q[0];
sx q[0];
rz(1.5968621) q[0];
rz(-2.3896396) q[2];
sx q[2];
rz(-0.30138902) q[2];
sx q[2];
rz(-1.4434659) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.5962332) q[1];
sx q[1];
rz(-1.4701478) q[1];
sx q[1];
rz(3.0920995) q[1];
rz(-pi) q[2];
rz(2.8392302) q[3];
sx q[3];
rz(-2.6609169) q[3];
sx q[3];
rz(-2.9722948) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.41708502) q[2];
sx q[2];
rz(-0.26792002) q[2];
sx q[2];
rz(0.99651304) q[2];
rz(0.35342446) q[3];
sx q[3];
rz(-0.74520183) q[3];
sx q[3];
rz(2.3341808) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0614232) q[0];
sx q[0];
rz(-0.81289476) q[0];
sx q[0];
rz(-2.9598575) q[0];
rz(-3.0985447) q[1];
sx q[1];
rz(-2.4964066) q[1];
sx q[1];
rz(0.28082401) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.94175324) q[0];
sx q[0];
rz(-2.0240677) q[0];
sx q[0];
rz(-2.0487294) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9279465) q[2];
sx q[2];
rz(-1.3823969) q[2];
sx q[2];
rz(-1.633916) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.909543) q[1];
sx q[1];
rz(-3.0735525) q[1];
sx q[1];
rz(2.1310991) q[1];
x q[2];
rz(0.19631581) q[3];
sx q[3];
rz(-1.8972978) q[3];
sx q[3];
rz(1.6895837) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.8250371) q[2];
sx q[2];
rz(-1.8871769) q[2];
sx q[2];
rz(2.5184856) q[2];
rz(1.0021707) q[3];
sx q[3];
rz(-1.3391756) q[3];
sx q[3];
rz(2.5785057) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6476718) q[0];
sx q[0];
rz(-1.5681842) q[0];
sx q[0];
rz(1.6012123) q[0];
rz(-2.2676246) q[1];
sx q[1];
rz(-1.0653492) q[1];
sx q[1];
rz(-3.031562) q[1];
rz(0.80550823) q[2];
sx q[2];
rz(-2.0279573) q[2];
sx q[2];
rz(1.3062994) q[2];
rz(0.86250967) q[3];
sx q[3];
rz(-2.6211092) q[3];
sx q[3];
rz(0.72343788) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
