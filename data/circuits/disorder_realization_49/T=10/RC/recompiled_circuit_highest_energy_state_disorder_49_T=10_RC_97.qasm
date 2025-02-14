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
rz(0.3857412) q[0];
sx q[0];
rz(2.1585611) q[0];
sx q[0];
rz(9.9821363) q[0];
rz(-1.9269257) q[1];
sx q[1];
rz(7.1375224) q[1];
sx q[1];
rz(8.4827276) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0165812) q[0];
sx q[0];
rz(-1.6894369) q[0];
sx q[0];
rz(-1.3604547) q[0];
rz(-pi) q[1];
rz(-0.73150191) q[2];
sx q[2];
rz(-0.78353751) q[2];
sx q[2];
rz(-2.6611626) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.57191932) q[1];
sx q[1];
rz(-1.633857) q[1];
sx q[1];
rz(-0.38625269) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.60987925) q[3];
sx q[3];
rz(-1.0945265) q[3];
sx q[3];
rz(-2.7000138) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.4152834) q[2];
sx q[2];
rz(-1.3741263) q[2];
sx q[2];
rz(-1.3709925) q[2];
rz(2.3224984) q[3];
sx q[3];
rz(-0.84688014) q[3];
sx q[3];
rz(3.0313361) q[3];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.46721989) q[0];
sx q[0];
rz(-1.9085566) q[0];
sx q[0];
rz(2.9357173) q[0];
rz(1.8241833) q[1];
sx q[1];
rz(-1.859788) q[1];
sx q[1];
rz(2.6726216) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9896133) q[0];
sx q[0];
rz(-1.5474404) q[0];
sx q[0];
rz(-0.65520136) q[0];
rz(2.9250336) q[2];
sx q[2];
rz(-1.1149841) q[2];
sx q[2];
rz(0.45879455) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.4552942) q[1];
sx q[1];
rz(-1.0166234) q[1];
sx q[1];
rz(0.55880736) q[1];
rz(-pi) q[2];
rz(-2.5303115) q[3];
sx q[3];
rz(-2.287902) q[3];
sx q[3];
rz(-2.1016013) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.59490243) q[2];
sx q[2];
rz(-0.83260834) q[2];
sx q[2];
rz(2.5140095) q[2];
rz(1.0707431) q[3];
sx q[3];
rz(-2.6379733) q[3];
sx q[3];
rz(0.32881769) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
rz(1.4048432) q[0];
sx q[0];
rz(-0.81031814) q[0];
sx q[0];
rz(2.3531083) q[0];
rz(-2.5704747) q[1];
sx q[1];
rz(-1.8126789) q[1];
sx q[1];
rz(-1.6956537) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.80966893) q[0];
sx q[0];
rz(-1.8394711) q[0];
sx q[0];
rz(0.013703811) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6035516) q[2];
sx q[2];
rz(-1.5376411) q[2];
sx q[2];
rz(-2.3743001) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.4169374) q[1];
sx q[1];
rz(-1.8677399) q[1];
sx q[1];
rz(-1.4663803) q[1];
rz(-pi) q[2];
rz(-0.51579185) q[3];
sx q[3];
rz(-0.63871562) q[3];
sx q[3];
rz(2.0206919) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.0995522) q[2];
sx q[2];
rz(-2.6991762) q[2];
sx q[2];
rz(0.37910795) q[2];
rz(1.7673309) q[3];
sx q[3];
rz(-0.97135025) q[3];
sx q[3];
rz(-0.9837392) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9336201) q[0];
sx q[0];
rz(-0.75478983) q[0];
sx q[0];
rz(-0.91243139) q[0];
rz(-1.747793) q[1];
sx q[1];
rz(-0.50519609) q[1];
sx q[1];
rz(2.8392653) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.63589225) q[0];
sx q[0];
rz(-1.5109343) q[0];
sx q[0];
rz(-0.12714956) q[0];
rz(-pi) q[1];
rz(0.85643391) q[2];
sx q[2];
rz(-1.8765479) q[2];
sx q[2];
rz(-0.46692525) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.6136377) q[1];
sx q[1];
rz(-0.21371811) q[1];
sx q[1];
rz(-1.3655013) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4563645) q[3];
sx q[3];
rz(-0.99894612) q[3];
sx q[3];
rz(0.88129191) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.48431524) q[2];
sx q[2];
rz(-2.4399098) q[2];
sx q[2];
rz(-0.62937984) q[2];
rz(-1.4194007) q[3];
sx q[3];
rz(-1.281176) q[3];
sx q[3];
rz(-2.4331376) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0486384) q[0];
sx q[0];
rz(-1.0372256) q[0];
sx q[0];
rz(-1.3193489) q[0];
rz(-1.0880148) q[1];
sx q[1];
rz(-1.271778) q[1];
sx q[1];
rz(-1.6768657) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.78878337) q[0];
sx q[0];
rz(-1.6029198) q[0];
sx q[0];
rz(-2.2706768) q[0];
rz(1.6032463) q[2];
sx q[2];
rz(-2.2032149) q[2];
sx q[2];
rz(-2.4995952) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.32680997) q[1];
sx q[1];
rz(-1.9091354) q[1];
sx q[1];
rz(-1.3698575) q[1];
rz(2.5358236) q[3];
sx q[3];
rz(-1.7585227) q[3];
sx q[3];
rz(2.5180205) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.69270837) q[2];
sx q[2];
rz(-0.11881891) q[2];
sx q[2];
rz(-0.23692712) q[2];
rz(0.98053011) q[3];
sx q[3];
rz(-1.4292932) q[3];
sx q[3];
rz(-0.94394365) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15810814) q[0];
sx q[0];
rz(-0.10949245) q[0];
sx q[0];
rz(-0.14516251) q[0];
rz(-1.6204087) q[1];
sx q[1];
rz(-2.3416134) q[1];
sx q[1];
rz(-3.1343585) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.67381672) q[0];
sx q[0];
rz(-2.4002808) q[0];
sx q[0];
rz(1.2178161) q[0];
x q[1];
rz(-2.5356338) q[2];
sx q[2];
rz(-2.3176386) q[2];
sx q[2];
rz(0.29203019) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.045102677) q[1];
sx q[1];
rz(-1.4670925) q[1];
sx q[1];
rz(1.8058877) q[1];
rz(-pi) q[2];
rz(0.16902828) q[3];
sx q[3];
rz(-0.16847408) q[3];
sx q[3];
rz(3.1210312) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.89473692) q[2];
sx q[2];
rz(-1.9523018) q[2];
sx q[2];
rz(-0.57752937) q[2];
rz(-1.9891116) q[3];
sx q[3];
rz(-2.5566176) q[3];
sx q[3];
rz(1.729689) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5685527) q[0];
sx q[0];
rz(-0.19565208) q[0];
sx q[0];
rz(2.0797119) q[0];
rz(1.3878239) q[1];
sx q[1];
rz(-1.9111218) q[1];
sx q[1];
rz(-1.498163) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.8138431) q[0];
sx q[0];
rz(-0.44506207) q[0];
sx q[0];
rz(0.31223865) q[0];
rz(0.64548026) q[2];
sx q[2];
rz(-1.3041851) q[2];
sx q[2];
rz(1.6763761) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.64387262) q[1];
sx q[1];
rz(-1.7959698) q[1];
sx q[1];
rz(2.1399645) q[1];
x q[2];
rz(-2.5779547) q[3];
sx q[3];
rz(-2.5803714) q[3];
sx q[3];
rz(-1.3432251) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.7428703) q[2];
sx q[2];
rz(-0.12464945) q[2];
sx q[2];
rz(-2.2275662) q[2];
rz(-2.3877609) q[3];
sx q[3];
rz(-1.9099312) q[3];
sx q[3];
rz(2.6851173) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1330426) q[0];
sx q[0];
rz(-0.32519105) q[0];
sx q[0];
rz(0.010566674) q[0];
rz(1.0039302) q[1];
sx q[1];
rz(-1.7251451) q[1];
sx q[1];
rz(-3.1026057) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9462601) q[0];
sx q[0];
rz(-1.8882013) q[0];
sx q[0];
rz(2.9575149) q[0];
rz(-0.029185295) q[2];
sx q[2];
rz(-0.31692908) q[2];
sx q[2];
rz(-2.2259797) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.85248831) q[1];
sx q[1];
rz(-0.90910289) q[1];
sx q[1];
rz(2.0579335) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4812864) q[3];
sx q[3];
rz(-2.2319372) q[3];
sx q[3];
rz(1.6474219) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.5616592) q[2];
sx q[2];
rz(-1.9114405) q[2];
sx q[2];
rz(-0.099543355) q[2];
rz(1.2933412) q[3];
sx q[3];
rz(-1.2852531) q[3];
sx q[3];
rz(-2.7660363) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.67000166) q[0];
sx q[0];
rz(-1.7741859) q[0];
sx q[0];
rz(1.8294096) q[0];
rz(0.52681628) q[1];
sx q[1];
rz(-2.3878038) q[1];
sx q[1];
rz(2.3048293) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6569516) q[0];
sx q[0];
rz(-1.27702) q[0];
sx q[0];
rz(2.8696069) q[0];
x q[1];
rz(0.86076059) q[2];
sx q[2];
rz(-0.75131159) q[2];
sx q[2];
rz(1.6565135) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.93406127) q[1];
sx q[1];
rz(-2.4467) q[1];
sx q[1];
rz(0.41568218) q[1];
x q[2];
rz(-2.760254) q[3];
sx q[3];
rz(-1.1313038) q[3];
sx q[3];
rz(0.56909305) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.3706563) q[2];
sx q[2];
rz(-2.5281361) q[2];
sx q[2];
rz(1.9343617) q[2];
rz(1.8218482) q[3];
sx q[3];
rz(-1.5650257) q[3];
sx q[3];
rz(2.3453662) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
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
rz(-3.0470444) q[0];
sx q[0];
rz(-0.47484174) q[0];
sx q[0];
rz(0.71665254) q[0];
rz(1.5776177) q[1];
sx q[1];
rz(-1.3037222) q[1];
sx q[1];
rz(2.5242453) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8113239) q[0];
sx q[0];
rz(-1.399938) q[0];
sx q[0];
rz(0.85612423) q[0];
rz(2.0390338) q[2];
sx q[2];
rz(-1.4371294) q[2];
sx q[2];
rz(-2.432761) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.3315369) q[1];
sx q[1];
rz(-2.202534) q[1];
sx q[1];
rz(-2.1239566) q[1];
x q[2];
rz(2.8598911) q[3];
sx q[3];
rz(-2.1920077) q[3];
sx q[3];
rz(3.1007183) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.8120332) q[2];
sx q[2];
rz(-0.38745189) q[2];
sx q[2];
rz(-2.6623902) q[2];
rz(1.6241578) q[3];
sx q[3];
rz(-1.4045709) q[3];
sx q[3];
rz(1.4994538) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0064405) q[0];
sx q[0];
rz(-1.7870446) q[0];
sx q[0];
rz(2.5743299) q[0];
rz(-2.3334423) q[1];
sx q[1];
rz(-1.5222526) q[1];
sx q[1];
rz(1.6076988) q[1];
rz(-1.7980747) q[2];
sx q[2];
rz(-2.6103316) q[2];
sx q[2];
rz(-1.1634315) q[2];
rz(0.2740673) q[3];
sx q[3];
rz(-0.82809877) q[3];
sx q[3];
rz(-0.86622151) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
