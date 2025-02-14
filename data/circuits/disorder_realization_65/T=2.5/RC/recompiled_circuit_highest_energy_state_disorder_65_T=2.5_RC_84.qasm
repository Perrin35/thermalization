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
rz(0.12198099) q[0];
sx q[0];
rz(3.0966336) q[0];
sx q[0];
rz(9.3293204) q[0];
rz(1.8010315) q[1];
sx q[1];
rz(-0.065923318) q[1];
sx q[1];
rz(0.59360582) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.18842342) q[0];
sx q[0];
rz(-1.0894408) q[0];
sx q[0];
rz(-1.166422) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6907211) q[2];
sx q[2];
rz(-1.7210135) q[2];
sx q[2];
rz(0.37301979) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.2954457) q[1];
sx q[1];
rz(-2.4762003) q[1];
sx q[1];
rz(0.34122463) q[1];
x q[2];
rz(-2.096522) q[3];
sx q[3];
rz(-1.8107256) q[3];
sx q[3];
rz(0.6434427) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(3.0815214) q[2];
sx q[2];
rz(-1.2520484) q[2];
sx q[2];
rz(-0.56212765) q[2];
rz(2.8020322) q[3];
sx q[3];
rz(-1.5226676) q[3];
sx q[3];
rz(1.08574) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.958441) q[0];
sx q[0];
rz(-0.14587942) q[0];
sx q[0];
rz(-1.1613783) q[0];
rz(-2.9005652) q[1];
sx q[1];
rz(-2.3190505) q[1];
sx q[1];
rz(-2.8384812) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68837092) q[0];
sx q[0];
rz(-1.3728598) q[0];
sx q[0];
rz(3.0730547) q[0];
rz(-pi) q[1];
rz(2.436736) q[2];
sx q[2];
rz(-0.19071707) q[2];
sx q[2];
rz(2.5097367) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.0063340291) q[1];
sx q[1];
rz(-0.63992441) q[1];
sx q[1];
rz(2.695822) q[1];
rz(-pi) q[2];
rz(2.0419028) q[3];
sx q[3];
rz(-2.0031824) q[3];
sx q[3];
rz(1.8023091) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.8384398) q[2];
sx q[2];
rz(-1.2829605) q[2];
sx q[2];
rz(-0.34070936) q[2];
rz(-0.43863145) q[3];
sx q[3];
rz(-0.80629587) q[3];
sx q[3];
rz(0.058569245) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4230147) q[0];
sx q[0];
rz(-2.5231762) q[0];
sx q[0];
rz(-2.3037236) q[0];
rz(-2.2184929) q[1];
sx q[1];
rz(-1.9337312) q[1];
sx q[1];
rz(-0.34742483) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9956857) q[0];
sx q[0];
rz(-1.6646228) q[0];
sx q[0];
rz(0.037414649) q[0];
rz(-pi) q[1];
rz(1.8272039) q[2];
sx q[2];
rz(-1.5359582) q[2];
sx q[2];
rz(-2.9098791) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.4975884) q[1];
sx q[1];
rz(-1.0217182) q[1];
sx q[1];
rz(2.1357029) q[1];
x q[2];
rz(0.33073552) q[3];
sx q[3];
rz(-1.8285255) q[3];
sx q[3];
rz(0.39604076) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.0190987) q[2];
sx q[2];
rz(-0.118003) q[2];
sx q[2];
rz(0.72354358) q[2];
rz(-1.2875693) q[3];
sx q[3];
rz(-2.5836594) q[3];
sx q[3];
rz(-0.045469835) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(-0.36640722) q[0];
sx q[0];
rz(-0.43908304) q[0];
sx q[0];
rz(-2.0570237) q[0];
rz(1.1868125) q[1];
sx q[1];
rz(-1.908952) q[1];
sx q[1];
rz(-1.67217) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.050768269) q[0];
sx q[0];
rz(-1.5406797) q[0];
sx q[0];
rz(1.5539021) q[0];
x q[1];
rz(1.9590366) q[2];
sx q[2];
rz(-1.2080384) q[2];
sx q[2];
rz(1.658939) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.332843) q[1];
sx q[1];
rz(-1.9522138) q[1];
sx q[1];
rz(-2.0705968) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3923548) q[3];
sx q[3];
rz(-0.47296528) q[3];
sx q[3];
rz(-1.7974896) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.12545086) q[2];
sx q[2];
rz(-2.2405388) q[2];
sx q[2];
rz(0.93552843) q[2];
rz(0.15141307) q[3];
sx q[3];
rz(-0.97984034) q[3];
sx q[3];
rz(2.6305731) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.22684111) q[0];
sx q[0];
rz(-1.5837357) q[0];
sx q[0];
rz(-0.75411183) q[0];
rz(-2.5304645) q[1];
sx q[1];
rz(-1.1437623) q[1];
sx q[1];
rz(-2.4653844) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0026086) q[0];
sx q[0];
rz(-1.7327285) q[0];
sx q[0];
rz(1.2404299) q[0];
x q[1];
rz(-1.0173747) q[2];
sx q[2];
rz(-0.62741919) q[2];
sx q[2];
rz(1.2469893) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.3067627) q[1];
sx q[1];
rz(-1.3993629) q[1];
sx q[1];
rz(-1.7333561) q[1];
rz(-pi) q[2];
rz(-0.12376484) q[3];
sx q[3];
rz(-0.58933697) q[3];
sx q[3];
rz(-2.2791354) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.0040697441) q[2];
sx q[2];
rz(-0.71228945) q[2];
sx q[2];
rz(0.75312692) q[2];
rz(2.5774041) q[3];
sx q[3];
rz(-1.1078395) q[3];
sx q[3];
rz(-2.3851725) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(0.18993987) q[0];
sx q[0];
rz(-0.48007444) q[0];
sx q[0];
rz(2.9204364) q[0];
rz(-3.0033374) q[1];
sx q[1];
rz(-1.6860551) q[1];
sx q[1];
rz(0.55603212) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.69029278) q[0];
sx q[0];
rz(-0.31421071) q[0];
sx q[0];
rz(0.67361535) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.37555917) q[2];
sx q[2];
rz(-1.2068159) q[2];
sx q[2];
rz(-0.057372626) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.1631524) q[1];
sx q[1];
rz(-1.5915055) q[1];
sx q[1];
rz(-2.9881575) q[1];
rz(-pi) q[2];
rz(-1.8121427) q[3];
sx q[3];
rz(-0.46306073) q[3];
sx q[3];
rz(1.1487414) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.14744645) q[2];
sx q[2];
rz(-2.3799956) q[2];
sx q[2];
rz(2.2046294) q[2];
rz(0.41038904) q[3];
sx q[3];
rz(-0.97987163) q[3];
sx q[3];
rz(0.65569896) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.74528247) q[0];
sx q[0];
rz(-0.87438011) q[0];
sx q[0];
rz(-2.1790047) q[0];
rz(0.67086041) q[1];
sx q[1];
rz(-0.52898359) q[1];
sx q[1];
rz(2.3765391) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6214692) q[0];
sx q[0];
rz(-1.4111821) q[0];
sx q[0];
rz(1.1667211) q[0];
rz(-0.069326055) q[2];
sx q[2];
rz(-1.6717807) q[2];
sx q[2];
rz(2.947383) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.9817762) q[1];
sx q[1];
rz(-0.66513956) q[1];
sx q[1];
rz(-2.239091) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.12419923) q[3];
sx q[3];
rz(-1.936125) q[3];
sx q[3];
rz(-0.78094966) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.44136167) q[2];
sx q[2];
rz(-0.92741489) q[2];
sx q[2];
rz(-2.3564763) q[2];
rz(-2.650812) q[3];
sx q[3];
rz(-1.1622585) q[3];
sx q[3];
rz(-2.6665915) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9908776) q[0];
sx q[0];
rz(-0.10902037) q[0];
sx q[0];
rz(0.14608598) q[0];
rz(0.27169216) q[1];
sx q[1];
rz(-0.79333317) q[1];
sx q[1];
rz(-0.21154107) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1871227) q[0];
sx q[0];
rz(-0.7899219) q[0];
sx q[0];
rz(-0.23440897) q[0];
rz(3.0258133) q[2];
sx q[2];
rz(-1.9405138) q[2];
sx q[2];
rz(-0.055495128) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.52741963) q[1];
sx q[1];
rz(-2.895311) q[1];
sx q[1];
rz(0.28913943) q[1];
rz(0.11084307) q[3];
sx q[3];
rz(-2.1680084) q[3];
sx q[3];
rz(-2.4750714) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.81296706) q[2];
sx q[2];
rz(-0.93572891) q[2];
sx q[2];
rz(0.82175559) q[2];
rz(1.3385319) q[3];
sx q[3];
rz(-1.0889564) q[3];
sx q[3];
rz(-0.72430044) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.61447918) q[0];
sx q[0];
rz(-0.69632691) q[0];
sx q[0];
rz(1.3592199) q[0];
rz(2.1837153) q[1];
sx q[1];
rz(-2.8045037) q[1];
sx q[1];
rz(-0.27761308) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2569297) q[0];
sx q[0];
rz(-0.96159928) q[0];
sx q[0];
rz(-2.4692187) q[0];
rz(-2.8238966) q[2];
sx q[2];
rz(-0.80149273) q[2];
sx q[2];
rz(1.7048175) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.902144) q[1];
sx q[1];
rz(-0.43593299) q[1];
sx q[1];
rz(-0.3050584) q[1];
rz(-2.7630799) q[3];
sx q[3];
rz(-1.7309824) q[3];
sx q[3];
rz(-1.2293564) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.9684888) q[2];
sx q[2];
rz(-2.5471881) q[2];
sx q[2];
rz(1.8572726) q[2];
rz(0.77749085) q[3];
sx q[3];
rz(-2.5992664) q[3];
sx q[3];
rz(2.8418181) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1182627) q[0];
sx q[0];
rz(-1.4813923) q[0];
sx q[0];
rz(-3.135664) q[0];
rz(1.802035) q[1];
sx q[1];
rz(-2.1052994) q[1];
sx q[1];
rz(-0.48066995) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4673432) q[0];
sx q[0];
rz(-3.1200069) q[0];
sx q[0];
rz(-2.4186446) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0239915) q[2];
sx q[2];
rz(-0.725774) q[2];
sx q[2];
rz(-0.99601907) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.0152032) q[1];
sx q[1];
rz(-1.9815832) q[1];
sx q[1];
rz(-2.4900764) q[1];
x q[2];
rz(1.3044429) q[3];
sx q[3];
rz(-1.7031809) q[3];
sx q[3];
rz(-1.3416106) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.9504451) q[2];
sx q[2];
rz(-0.55926776) q[2];
sx q[2];
rz(2.9648103) q[2];
rz(0.92987972) q[3];
sx q[3];
rz(-1.9527438) q[3];
sx q[3];
rz(-2.1117579) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.69242351) q[0];
sx q[0];
rz(-1.4207191) q[0];
sx q[0];
rz(-1.0906013) q[0];
rz(-1.0724267) q[1];
sx q[1];
rz(-1.6630395) q[1];
sx q[1];
rz(1.6987775) q[1];
rz(-2.6497447) q[2];
sx q[2];
rz(-1.5651697) q[2];
sx q[2];
rz(-1.6357505) q[2];
rz(-3.0092893) q[3];
sx q[3];
rz(-2.5322661) q[3];
sx q[3];
rz(0.54064396) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
