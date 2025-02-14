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
rz(-2.1844644) q[0];
sx q[0];
rz(-1.6532093) q[0];
sx q[0];
rz(-1.1487577) q[0];
rz(0.59612885) q[1];
sx q[1];
rz(-0.40607536) q[1];
sx q[1];
rz(-1.4546855) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5243149) q[0];
sx q[0];
rz(-1.9348659) q[0];
sx q[0];
rz(0.84366701) q[0];
rz(-pi) q[1];
rz(0.11496719) q[2];
sx q[2];
rz(-0.76412725) q[2];
sx q[2];
rz(-0.72944356) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.3439181) q[1];
sx q[1];
rz(-0.82268381) q[1];
sx q[1];
rz(1.3901346) q[1];
x q[2];
rz(-0.65806453) q[3];
sx q[3];
rz(-1.6474373) q[3];
sx q[3];
rz(-1.5358401) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.9601606) q[2];
sx q[2];
rz(-0.97603193) q[2];
sx q[2];
rz(1.9350447) q[2];
rz(1.9801961) q[3];
sx q[3];
rz(-0.56777081) q[3];
sx q[3];
rz(-1.5706221) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1777765) q[0];
sx q[0];
rz(-2.0160567) q[0];
sx q[0];
rz(0.97208446) q[0];
rz(2.8731335) q[1];
sx q[1];
rz(-1.700289) q[1];
sx q[1];
rz(-2.6867902) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9869315) q[0];
sx q[0];
rz(-1.6973746) q[0];
sx q[0];
rz(-1.8693699) q[0];
x q[1];
rz(-1.182453) q[2];
sx q[2];
rz(-1.4253311) q[2];
sx q[2];
rz(2.5487473) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.6356335) q[1];
sx q[1];
rz(-1.1065496) q[1];
sx q[1];
rz(1.5134769) q[1];
x q[2];
rz(1.0271026) q[3];
sx q[3];
rz(-1.0732526) q[3];
sx q[3];
rz(3.080201) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.0316169) q[2];
sx q[2];
rz(-0.72828186) q[2];
sx q[2];
rz(2.1419683) q[2];
rz(0.085974606) q[3];
sx q[3];
rz(-1.5558259) q[3];
sx q[3];
rz(-0.011367817) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40083945) q[0];
sx q[0];
rz(-1.4680306) q[0];
sx q[0];
rz(0.22415796) q[0];
rz(1.594918) q[1];
sx q[1];
rz(-1.5983351) q[1];
sx q[1];
rz(-1.3067783) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49840333) q[0];
sx q[0];
rz(-2.1191594) q[0];
sx q[0];
rz(-1.7140618) q[0];
x q[1];
rz(2.4073657) q[2];
sx q[2];
rz(-2.1598744) q[2];
sx q[2];
rz(-2.7393419) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.1968888) q[1];
sx q[1];
rz(-2.534229) q[1];
sx q[1];
rz(-0.783226) q[1];
x q[2];
rz(-1.8238092) q[3];
sx q[3];
rz(-0.91716498) q[3];
sx q[3];
rz(-2.831649) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.811502) q[2];
sx q[2];
rz(-1.6197438) q[2];
sx q[2];
rz(-2.1211993) q[2];
rz(-1.329782) q[3];
sx q[3];
rz(-1.1251757) q[3];
sx q[3];
rz(1.5248732) q[3];
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
rz(-pi) q[3];
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
rz(-0.40600768) q[0];
sx q[0];
rz(-1.0574874) q[0];
sx q[0];
rz(1.1757346) q[0];
rz(-0.27578393) q[1];
sx q[1];
rz(-0.83668721) q[1];
sx q[1];
rz(2.5036459) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36786554) q[0];
sx q[0];
rz(-1.6035032) q[0];
sx q[0];
rz(-1.5814533) q[0];
rz(-pi) q[1];
x q[1];
rz(0.3023382) q[2];
sx q[2];
rz(-2.4983642) q[2];
sx q[2];
rz(2.0143721) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.13220683) q[1];
sx q[1];
rz(-1.8957545) q[1];
sx q[1];
rz(2.4924949) q[1];
x q[2];
rz(0.31973394) q[3];
sx q[3];
rz(-1.481062) q[3];
sx q[3];
rz(3.0872726) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.63567579) q[2];
sx q[2];
rz(-1.5779147) q[2];
sx q[2];
rz(-1.4471311) q[2];
rz(-0.21008374) q[3];
sx q[3];
rz(-0.45322067) q[3];
sx q[3];
rz(-0.92857462) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2337445) q[0];
sx q[0];
rz(-0.15234983) q[0];
sx q[0];
rz(-2.0727169) q[0];
rz(1.2999889) q[1];
sx q[1];
rz(-1.4062107) q[1];
sx q[1];
rz(-1.9357505) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5652454) q[0];
sx q[0];
rz(-1.6505213) q[0];
sx q[0];
rz(-1.1673156) q[0];
rz(-pi) q[1];
x q[1];
rz(0.0095179518) q[2];
sx q[2];
rz(-2.2235907) q[2];
sx q[2];
rz(1.643484) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.5609392) q[1];
sx q[1];
rz(-1.4672797) q[1];
sx q[1];
rz(-0.85080457) q[1];
x q[2];
rz(-0.48128328) q[3];
sx q[3];
rz(-1.4247155) q[3];
sx q[3];
rz(-1.2969601) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.892889) q[2];
sx q[2];
rz(-0.37962571) q[2];
sx q[2];
rz(-3.0730263) q[2];
rz(0.93962234) q[3];
sx q[3];
rz(-1.7328123) q[3];
sx q[3];
rz(-1.2287963) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.762961) q[0];
sx q[0];
rz(-1.7447504) q[0];
sx q[0];
rz(-0.5214386) q[0];
rz(1.5154845) q[1];
sx q[1];
rz(-1.3896959) q[1];
sx q[1];
rz(0.62561402) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8541295) q[0];
sx q[0];
rz(-1.03106) q[0];
sx q[0];
rz(3.1053393) q[0];
rz(1.6997153) q[2];
sx q[2];
rz(-1.9471418) q[2];
sx q[2];
rz(-2.3835045) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.8486937) q[1];
sx q[1];
rz(-1.1176511) q[1];
sx q[1];
rz(-2.921656) q[1];
rz(-3.0066177) q[3];
sx q[3];
rz(-0.42536456) q[3];
sx q[3];
rz(-0.36799001) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.0079415) q[2];
sx q[2];
rz(-0.59591728) q[2];
sx q[2];
rz(-0.1055183) q[2];
rz(-2.1775235) q[3];
sx q[3];
rz(-1.9897507) q[3];
sx q[3];
rz(0.9849557) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0725726) q[0];
sx q[0];
rz(-2.9264086) q[0];
sx q[0];
rz(-1.4129114) q[0];
rz(0.37881306) q[1];
sx q[1];
rz(-1.2544931) q[1];
sx q[1];
rz(-2.5884195) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.39554292) q[0];
sx q[0];
rz(-0.93727124) q[0];
sx q[0];
rz(-2.6468011) q[0];
rz(2.5752284) q[2];
sx q[2];
rz(-0.16728044) q[2];
sx q[2];
rz(-2.9110514) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.70771101) q[1];
sx q[1];
rz(-1.6110782) q[1];
sx q[1];
rz(-0.28357419) q[1];
rz(-pi) q[2];
rz(-1.3101226) q[3];
sx q[3];
rz(-1.3743152) q[3];
sx q[3];
rz(0.5082265) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.6806014) q[2];
sx q[2];
rz(-2.1059771) q[2];
sx q[2];
rz(-0.41213948) q[2];
rz(-2.4814217) q[3];
sx q[3];
rz(-2.6001402) q[3];
sx q[3];
rz(-0.49303833) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.93765813) q[0];
sx q[0];
rz(-1.0451319) q[0];
sx q[0];
rz(2.9397021) q[0];
rz(-2.0837325) q[1];
sx q[1];
rz(-0.36483929) q[1];
sx q[1];
rz(-0.26022628) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5903918) q[0];
sx q[0];
rz(-1.0559715) q[0];
sx q[0];
rz(-1.1333254) q[0];
x q[1];
rz(2.0864395) q[2];
sx q[2];
rz(-1.1844352) q[2];
sx q[2];
rz(2.6576415) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.59108666) q[1];
sx q[1];
rz(-0.45190736) q[1];
sx q[1];
rz(0.68303789) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.21589888) q[3];
sx q[3];
rz(-2.4697018) q[3];
sx q[3];
rz(0.011530576) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.95320931) q[2];
sx q[2];
rz(-1.7069495) q[2];
sx q[2];
rz(2.5592213) q[2];
rz(-2.6992056) q[3];
sx q[3];
rz(-1.9059537) q[3];
sx q[3];
rz(-0.83627397) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9198832) q[0];
sx q[0];
rz(-1.5279122) q[0];
sx q[0];
rz(1.7145994) q[0];
rz(-1.3555948) q[1];
sx q[1];
rz(-1.4878863) q[1];
sx q[1];
rz(-1.4367163) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6528931) q[0];
sx q[0];
rz(-2.821229) q[0];
sx q[0];
rz(-1.5071177) q[0];
rz(-0.53578761) q[2];
sx q[2];
rz(-2.5965207) q[2];
sx q[2];
rz(1.6824012) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.88594) q[1];
sx q[1];
rz(-2.8999834) q[1];
sx q[1];
rz(1.5109946) q[1];
rz(-1.0441699) q[3];
sx q[3];
rz(-1.7071402) q[3];
sx q[3];
rz(0.52140331) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.23472486) q[2];
sx q[2];
rz(-1.2805254) q[2];
sx q[2];
rz(-1.5391763) q[2];
rz(2.5336044) q[3];
sx q[3];
rz(-1.4092849) q[3];
sx q[3];
rz(-0.68979818) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7749087) q[0];
sx q[0];
rz(-2.4803949) q[0];
sx q[0];
rz(2.3181424) q[0];
rz(-0.89556328) q[1];
sx q[1];
rz(-1.7915553) q[1];
sx q[1];
rz(-2.7604738) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6247647) q[0];
sx q[0];
rz(-0.73113686) q[0];
sx q[0];
rz(0.69964377) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7627583) q[2];
sx q[2];
rz(-0.4288579) q[2];
sx q[2];
rz(-0.90012401) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.1824519) q[1];
sx q[1];
rz(-0.63196665) q[1];
sx q[1];
rz(-2.0672634) q[1];
rz(-pi) q[2];
x q[2];
rz(0.62750979) q[3];
sx q[3];
rz(-2.0352279) q[3];
sx q[3];
rz(0.04963683) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.33619189) q[2];
sx q[2];
rz(-2.7298584) q[2];
sx q[2];
rz(2.1545048) q[2];
rz(1.6030715) q[3];
sx q[3];
rz(-0.58731949) q[3];
sx q[3];
rz(0.2400329) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9170452) q[0];
sx q[0];
rz(-1.9018835) q[0];
sx q[0];
rz(-1.6214669) q[0];
rz(0.70980258) q[1];
sx q[1];
rz(-2.5820422) q[1];
sx q[1];
rz(-1.6834264) q[1];
rz(-1.456121) q[2];
sx q[2];
rz(-1.2085452) q[2];
sx q[2];
rz(2.1160407) q[2];
rz(2.5153574) q[3];
sx q[3];
rz(-2.2953485) q[3];
sx q[3];
rz(2.9957537) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
