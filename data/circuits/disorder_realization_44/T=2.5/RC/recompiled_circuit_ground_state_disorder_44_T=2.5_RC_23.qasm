OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.7437416) q[0];
sx q[0];
rz(-0.97111094) q[0];
sx q[0];
rz(2.4308393) q[0];
rz(0.96762586) q[1];
sx q[1];
rz(-3.1275446) q[1];
sx q[1];
rz(0.79298055) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4913038) q[0];
sx q[0];
rz(-2.2772335) q[0];
sx q[0];
rz(0.44235787) q[0];
x q[1];
rz(-1.8618199) q[2];
sx q[2];
rz(-3.0918192) q[2];
sx q[2];
rz(2.7437944) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.087634511) q[1];
sx q[1];
rz(-2.2184508) q[1];
sx q[1];
rz(2.1235076) q[1];
x q[2];
rz(2.3550911) q[3];
sx q[3];
rz(-2.3701027) q[3];
sx q[3];
rz(-2.2685693) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.8636318) q[2];
sx q[2];
rz(-0.42676723) q[2];
sx q[2];
rz(-0.58941907) q[2];
rz(0.74728084) q[3];
sx q[3];
rz(-3.0487479) q[3];
sx q[3];
rz(-2.5417627) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0624307) q[0];
sx q[0];
rz(-0.23167647) q[0];
sx q[0];
rz(0.9011426) q[0];
rz(0.98183739) q[1];
sx q[1];
rz(-0.58126175) q[1];
sx q[1];
rz(2.1493886) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.85800303) q[0];
sx q[0];
rz(-2.186728) q[0];
sx q[0];
rz(-2.9735656) q[0];
x q[1];
rz(-2.9770933) q[2];
sx q[2];
rz(-0.82969159) q[2];
sx q[2];
rz(-2.6890597) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.7957669) q[1];
sx q[1];
rz(-0.92928934) q[1];
sx q[1];
rz(-2.262425) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3808591) q[3];
sx q[3];
rz(-0.84322107) q[3];
sx q[3];
rz(0.47970495) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.046752669) q[2];
sx q[2];
rz(-2.6111626) q[2];
sx q[2];
rz(-0.025010427) q[2];
rz(-2.181634) q[3];
sx q[3];
rz(-3.0955866) q[3];
sx q[3];
rz(-3.0733601) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.18441021) q[0];
sx q[0];
rz(-3.130641) q[0];
sx q[0];
rz(2.4175194) q[0];
rz(-0.11101668) q[1];
sx q[1];
rz(-2.5357775) q[1];
sx q[1];
rz(3.1287126) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6136821) q[0];
sx q[0];
rz(-2.9146195) q[0];
sx q[0];
rz(-1.2140973) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.78978844) q[2];
sx q[2];
rz(-0.26726535) q[2];
sx q[2];
rz(-0.59556669) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.6659505) q[1];
sx q[1];
rz(-0.29120177) q[1];
sx q[1];
rz(2.5282574) q[1];
rz(-pi) q[2];
rz(0.88408459) q[3];
sx q[3];
rz(-1.8662765) q[3];
sx q[3];
rz(-2.0990163) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.24590242) q[2];
sx q[2];
rz(-0.7220214) q[2];
sx q[2];
rz(-0.32430172) q[2];
rz(-0.4970099) q[3];
sx q[3];
rz(-2.9091166) q[3];
sx q[3];
rz(2.9089109) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4683891) q[0];
sx q[0];
rz(-2.9717746) q[0];
sx q[0];
rz(-2.66535) q[0];
rz(-2.6561123) q[1];
sx q[1];
rz(-0.49518934) q[1];
sx q[1];
rz(2.9229497) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6784558) q[0];
sx q[0];
rz(-1.3443724) q[0];
sx q[0];
rz(-1.1562267) q[0];
rz(-pi) q[1];
rz(1.5831469) q[2];
sx q[2];
rz(-2.2943779) q[2];
sx q[2];
rz(-2.6820094) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.2656897) q[1];
sx q[1];
rz(-0.69702083) q[1];
sx q[1];
rz(-0.69209309) q[1];
rz(-pi) q[2];
rz(0.63204445) q[3];
sx q[3];
rz(-0.38292745) q[3];
sx q[3];
rz(-2.6296089) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.789088) q[2];
sx q[2];
rz(-0.75864351) q[2];
sx q[2];
rz(2.5701806) q[2];
rz(0.93829751) q[3];
sx q[3];
rz(-0.58656991) q[3];
sx q[3];
rz(2.9686484) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.57662302) q[0];
sx q[0];
rz(-2.7374856) q[0];
sx q[0];
rz(2.670766) q[0];
rz(2.2305523) q[1];
sx q[1];
rz(-2.7016579) q[1];
sx q[1];
rz(-0.43112531) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3456988) q[0];
sx q[0];
rz(-0.78471334) q[0];
sx q[0];
rz(0.037952947) q[0];
rz(0.25283739) q[2];
sx q[2];
rz(-2.188854) q[2];
sx q[2];
rz(2.035066) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.8226128) q[1];
sx q[1];
rz(-3.0676004) q[1];
sx q[1];
rz(-1.1896108) q[1];
rz(2.9399559) q[3];
sx q[3];
rz(-1.4343702) q[3];
sx q[3];
rz(2.0755943) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.33991459) q[2];
sx q[2];
rz(-0.0047923294) q[2];
sx q[2];
rz(2.8002296) q[2];
rz(2.7692128) q[3];
sx q[3];
rz(-2.5058993) q[3];
sx q[3];
rz(2.671833) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85207283) q[0];
sx q[0];
rz(-2.2884123) q[0];
sx q[0];
rz(0.12292718) q[0];
rz(-2.7561103) q[1];
sx q[1];
rz(-0.79817927) q[1];
sx q[1];
rz(2.4179359) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.052304) q[0];
sx q[0];
rz(-1.6160242) q[0];
sx q[0];
rz(-1.7448698) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.18823024) q[2];
sx q[2];
rz(-1.2410795) q[2];
sx q[2];
rz(-0.70590574) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.78711975) q[1];
sx q[1];
rz(-0.90098375) q[1];
sx q[1];
rz(-0.77490999) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9272363) q[3];
sx q[3];
rz(-1.949614) q[3];
sx q[3];
rz(0.54730584) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.75189292) q[2];
sx q[2];
rz(-0.60201001) q[2];
sx q[2];
rz(-0.67328084) q[2];
rz(0.26345396) q[3];
sx q[3];
rz(-0.47502381) q[3];
sx q[3];
rz(2.8365005) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.70169705) q[0];
sx q[0];
rz(-2.3450527) q[0];
sx q[0];
rz(0.51629603) q[0];
rz(-2.3856178) q[1];
sx q[1];
rz(-2.1831903) q[1];
sx q[1];
rz(0.36852536) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4936016) q[0];
sx q[0];
rz(-0.89009919) q[0];
sx q[0];
rz(-2.2590911) q[0];
rz(-pi) q[1];
rz(-0.045727878) q[2];
sx q[2];
rz(-0.70529443) q[2];
sx q[2];
rz(-0.7208342) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.8766899) q[1];
sx q[1];
rz(-1.6769092) q[1];
sx q[1];
rz(1.4011032) q[1];
rz(2.1462899) q[3];
sx q[3];
rz(-0.80484521) q[3];
sx q[3];
rz(2.9965049) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.4964909) q[2];
sx q[2];
rz(-3.0848905) q[2];
sx q[2];
rz(-2.6594095) q[2];
rz(2.9218946) q[3];
sx q[3];
rz(-0.69742656) q[3];
sx q[3];
rz(-2.4612332) q[3];
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
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41736233) q[0];
sx q[0];
rz(-0.14886947) q[0];
sx q[0];
rz(-2.505488) q[0];
rz(-0.96027374) q[1];
sx q[1];
rz(-0.93820131) q[1];
sx q[1];
rz(-2.2669534) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8281404) q[0];
sx q[0];
rz(-2.9813672) q[0];
sx q[0];
rz(1.942722) q[0];
rz(-3.0205171) q[2];
sx q[2];
rz(-1.6796475) q[2];
sx q[2];
rz(-1.4360652) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.88747178) q[1];
sx q[1];
rz(-1.9560818) q[1];
sx q[1];
rz(-2.0068568) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2611924) q[3];
sx q[3];
rz(-1.0171618) q[3];
sx q[3];
rz(-2.2601489) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.8735698) q[2];
sx q[2];
rz(-2.7251254) q[2];
sx q[2];
rz(2.2250309) q[2];
rz(-0.77740866) q[3];
sx q[3];
rz(-2.58367) q[3];
sx q[3];
rz(2.6072445) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.022631835) q[0];
sx q[0];
rz(-2.4038778) q[0];
sx q[0];
rz(-2.90888) q[0];
rz(-2.4932056) q[1];
sx q[1];
rz(-2.5189923) q[1];
sx q[1];
rz(0.56786215) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2875339) q[0];
sx q[0];
rz(-1.9224478) q[0];
sx q[0];
rz(-1.132347) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8428734) q[2];
sx q[2];
rz(-1.3634184) q[2];
sx q[2];
rz(-0.033471154) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.372124) q[1];
sx q[1];
rz(-0.25126496) q[1];
sx q[1];
rz(-2.1970046) q[1];
x q[2];
rz(-2.4106823) q[3];
sx q[3];
rz(-0.95940351) q[3];
sx q[3];
rz(2.9711593) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.76967543) q[2];
sx q[2];
rz(-2.9480675) q[2];
sx q[2];
rz(-2.589321) q[2];
rz(0.28198379) q[3];
sx q[3];
rz(-2.7115287) q[3];
sx q[3];
rz(-3.0388487) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2118537) q[0];
sx q[0];
rz(-0.064082853) q[0];
sx q[0];
rz(0.1317568) q[0];
rz(0.53648221) q[1];
sx q[1];
rz(-2.9832612) q[1];
sx q[1];
rz(0.90824711) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2836766) q[0];
sx q[0];
rz(-2.2649797) q[0];
sx q[0];
rz(-0.038900872) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.26226173) q[2];
sx q[2];
rz(-2.024767) q[2];
sx q[2];
rz(-1.639546) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.40501198) q[1];
sx q[1];
rz(-0.93376505) q[1];
sx q[1];
rz(-0.65959658) q[1];
rz(-pi) q[2];
rz(-0.0018168505) q[3];
sx q[3];
rz(-1.7967136) q[3];
sx q[3];
rz(0.29669138) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.6324255) q[2];
sx q[2];
rz(-0.91370344) q[2];
sx q[2];
rz(2.8636279) q[2];
rz(0.5667423) q[3];
sx q[3];
rz(-2.9394579) q[3];
sx q[3];
rz(-2.22866) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81000281) q[0];
sx q[0];
rz(-1.4053874) q[0];
sx q[0];
rz(-0.94595861) q[0];
rz(-0.48238659) q[1];
sx q[1];
rz(-1.9079897) q[1];
sx q[1];
rz(2.447396) q[1];
rz(-2.161088) q[2];
sx q[2];
rz(-1.7750778) q[2];
sx q[2];
rz(1.9951174) q[2];
rz(1.514074) q[3];
sx q[3];
rz(-2.6346907) q[3];
sx q[3];
rz(-2.9148867) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
