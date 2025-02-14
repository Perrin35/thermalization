OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.39785102) q[0];
sx q[0];
rz(4.1127036) q[0];
sx q[0];
rz(10.135531) q[0];
rz(-2.1739668) q[1];
sx q[1];
rz(-0.014048014) q[1];
sx q[1];
rz(2.3486121) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.377787) q[0];
sx q[0];
rz(-1.9024693) q[0];
sx q[0];
rz(0.8140696) q[0];
rz(1.5231126) q[2];
sx q[2];
rz(-1.5565201) q[2];
sx q[2];
rz(-0.88231495) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.087634511) q[1];
sx q[1];
rz(-2.2184508) q[1];
sx q[1];
rz(-2.1235076) q[1];
rz(-pi) q[2];
rz(-0.96785424) q[3];
sx q[3];
rz(-2.0857028) q[3];
sx q[3];
rz(3.0647354) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.8636318) q[2];
sx q[2];
rz(-0.42676723) q[2];
sx q[2];
rz(-0.58941907) q[2];
rz(-2.3943118) q[3];
sx q[3];
rz(-3.0487479) q[3];
sx q[3];
rz(0.59982991) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0624307) q[0];
sx q[0];
rz(-0.23167647) q[0];
sx q[0];
rz(-0.9011426) q[0];
rz(-0.98183739) q[1];
sx q[1];
rz(-2.5603309) q[1];
sx q[1];
rz(2.1493886) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.57241523) q[0];
sx q[0];
rz(-0.63555866) q[0];
sx q[0];
rz(-1.8028238) q[0];
rz(-pi) q[1];
rz(1.7478701) q[2];
sx q[2];
rz(-0.75573993) q[2];
sx q[2];
rz(-0.21142879) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.3766831) q[1];
sx q[1];
rz(-1.0344532) q[1];
sx q[1];
rz(-2.3715644) q[1];
rz(-pi) q[2];
rz(-2.3808591) q[3];
sx q[3];
rz(-0.84322107) q[3];
sx q[3];
rz(-2.6618877) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.046752669) q[2];
sx q[2];
rz(-2.6111626) q[2];
sx q[2];
rz(-0.025010427) q[2];
rz(0.95995861) q[3];
sx q[3];
rz(-3.0955866) q[3];
sx q[3];
rz(-3.0733601) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.18441021) q[0];
sx q[0];
rz(-3.130641) q[0];
sx q[0];
rz(-2.4175194) q[0];
rz(3.030576) q[1];
sx q[1];
rz(-0.60581517) q[1];
sx q[1];
rz(-3.1287126) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4469876) q[0];
sx q[0];
rz(-1.6494538) q[0];
sx q[0];
rz(-1.3576677) q[0];
x q[1];
rz(-0.19042966) q[2];
sx q[2];
rz(-1.7594751) q[2];
sx q[2];
rz(1.3946472) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.9835123) q[1];
sx q[1];
rz(-1.8077824) q[1];
sx q[1];
rz(-1.7416341) q[1];
rz(2.7666088) q[3];
sx q[3];
rz(-2.2224748) q[3];
sx q[3];
rz(0.76256547) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.24590242) q[2];
sx q[2];
rz(-0.7220214) q[2];
sx q[2];
rz(-0.32430172) q[2];
rz(-2.6445828) q[3];
sx q[3];
rz(-2.9091166) q[3];
sx q[3];
rz(0.23268172) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
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
rz(-2.6464033) q[1];
sx q[1];
rz(0.21864299) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2061342) q[0];
sx q[0];
rz(-1.1674321) q[0];
sx q[0];
rz(-2.8950188) q[0];
rz(-1.5584458) q[2];
sx q[2];
rz(-2.2943779) q[2];
sx q[2];
rz(0.45958322) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-3.0899892) q[1];
sx q[1];
rz(-1.0538424) q[1];
sx q[1];
rz(-2.0614784) q[1];
rz(-pi) q[2];
x q[2];
rz(0.63204445) q[3];
sx q[3];
rz(-2.7586652) q[3];
sx q[3];
rz(2.6296089) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.3525047) q[2];
sx q[2];
rz(-0.75864351) q[2];
sx q[2];
rz(2.5701806) q[2];
rz(2.2032951) q[3];
sx q[3];
rz(-2.5550227) q[3];
sx q[3];
rz(-0.17294426) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
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
rz(-0.47082666) q[0];
rz(2.2305523) q[1];
sx q[1];
rz(-2.7016579) q[1];
sx q[1];
rz(2.7104673) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.74228263) q[0];
sx q[0];
rz(-2.3547908) q[0];
sx q[0];
rz(1.6086701) q[0];
rz(-pi) q[1];
x q[1];
rz(0.93740119) q[2];
sx q[2];
rz(-1.7761111) q[2];
sx q[2];
rz(-2.5287153) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.8226128) q[1];
sx q[1];
rz(-3.0676004) q[1];
sx q[1];
rz(-1.9519819) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5406932) q[3];
sx q[3];
rz(-2.89866) q[3];
sx q[3];
rz(-3.0593135) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.33991459) q[2];
sx q[2];
rz(-0.0047923294) q[2];
sx q[2];
rz(2.8002296) q[2];
rz(-2.7692128) q[3];
sx q[3];
rz(-2.5058993) q[3];
sx q[3];
rz(0.46975964) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2895198) q[0];
sx q[0];
rz(-0.85318035) q[0];
sx q[0];
rz(-3.0186655) q[0];
rz(2.7561103) q[1];
sx q[1];
rz(-2.3434134) q[1];
sx q[1];
rz(2.4179359) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0892886) q[0];
sx q[0];
rz(-1.5255685) q[0];
sx q[0];
rz(1.7448698) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9060109) q[2];
sx q[2];
rz(-1.7487757) q[2];
sx q[2];
rz(2.2151057) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.8116709) q[1];
sx q[1];
rz(-2.151346) q[1];
sx q[1];
rz(2.4075721) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0615929) q[3];
sx q[3];
rz(-0.43269581) q[3];
sx q[3];
rz(3.1263292) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.75189292) q[2];
sx q[2];
rz(-0.60201001) q[2];
sx q[2];
rz(2.4683118) q[2];
rz(2.8781387) q[3];
sx q[3];
rz(-0.47502381) q[3];
sx q[3];
rz(0.30509216) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.70169705) q[0];
sx q[0];
rz(-0.79653996) q[0];
sx q[0];
rz(-2.6252966) q[0];
rz(-0.75597489) q[1];
sx q[1];
rz(-2.1831903) q[1];
sx q[1];
rz(-0.36852536) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8719292) q[0];
sx q[0];
rz(-0.9270398) q[0];
sx q[0];
rz(2.4764349) q[0];
x q[1];
rz(3.0958648) q[2];
sx q[2];
rz(-0.70529443) q[2];
sx q[2];
rz(-0.7208342) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.8175537) q[1];
sx q[1];
rz(-1.7395258) q[1];
sx q[1];
rz(0.1076474) q[1];
x q[2];
rz(-2.6266699) q[3];
sx q[3];
rz(-0.92149261) q[3];
sx q[3];
rz(-2.2442371) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4964909) q[2];
sx q[2];
rz(-3.0848905) q[2];
sx q[2];
rz(-2.6594095) q[2];
rz(-2.9218946) q[3];
sx q[3];
rz(-2.4441661) q[3];
sx q[3];
rz(0.68035948) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41736233) q[0];
sx q[0];
rz(-0.14886947) q[0];
sx q[0];
rz(0.6361047) q[0];
rz(-0.96027374) q[1];
sx q[1];
rz(-0.93820131) q[1];
sx q[1];
rz(-2.2669534) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4518294) q[0];
sx q[0];
rz(-1.7199819) q[0];
sx q[0];
rz(-0.058663603) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4611488) q[2];
sx q[2];
rz(-1.4504408) q[2];
sx q[2];
rz(-0.14794825) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.284901) q[1];
sx q[1];
rz(-1.1686348) q[1];
sx q[1];
rz(-0.42070893) q[1];
rz(-pi) q[2];
rz(-0.80031826) q[3];
sx q[3];
rz(-2.2860675) q[3];
sx q[3];
rz(-3.0190938) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.26802289) q[2];
sx q[2];
rz(-2.7251254) q[2];
sx q[2];
rz(0.91656172) q[2];
rz(-2.364184) q[3];
sx q[3];
rz(-0.55792266) q[3];
sx q[3];
rz(2.6072445) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1189608) q[0];
sx q[0];
rz(-0.73771483) q[0];
sx q[0];
rz(0.23271261) q[0];
rz(-2.4932056) q[1];
sx q[1];
rz(-2.5189923) q[1];
sx q[1];
rz(0.56786215) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0584326) q[0];
sx q[0];
rz(-0.55480236) q[0];
sx q[0];
rz(-2.2834999) q[0];
x q[1];
rz(-0.90654208) q[2];
sx q[2];
rz(-2.8010606) q[2];
sx q[2];
rz(-2.1730925) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.1901924) q[1];
sx q[1];
rz(-1.4245598) q[1];
sx q[1];
rz(1.7758572) q[1];
x q[2];
rz(0.73091032) q[3];
sx q[3];
rz(-0.95940351) q[3];
sx q[3];
rz(2.9711593) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.76967543) q[2];
sx q[2];
rz(-2.9480675) q[2];
sx q[2];
rz(2.589321) q[2];
rz(0.28198379) q[3];
sx q[3];
rz(-0.43006399) q[3];
sx q[3];
rz(3.0388487) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2118537) q[0];
sx q[0];
rz(-3.0775098) q[0];
sx q[0];
rz(-0.1317568) q[0];
rz(0.53648221) q[1];
sx q[1];
rz(-0.15833144) q[1];
sx q[1];
rz(-0.90824711) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85791608) q[0];
sx q[0];
rz(-0.87661298) q[0];
sx q[0];
rz(-3.1026918) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0824049) q[2];
sx q[2];
rz(-0.51967144) q[2];
sx q[2];
rz(-2.1888417) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.40501198) q[1];
sx q[1];
rz(-2.2078276) q[1];
sx q[1];
rz(2.4819961) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.0018168505) q[3];
sx q[3];
rz(-1.7967136) q[3];
sx q[3];
rz(0.29669138) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.6324255) q[2];
sx q[2];
rz(-2.2278892) q[2];
sx q[2];
rz(2.8636279) q[2];
rz(0.5667423) q[3];
sx q[3];
rz(-2.9394579) q[3];
sx q[3];
rz(0.91293269) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3315898) q[0];
sx q[0];
rz(-1.7362052) q[0];
sx q[0];
rz(2.195634) q[0];
rz(0.48238659) q[1];
sx q[1];
rz(-1.2336029) q[1];
sx q[1];
rz(-0.69419669) q[1];
rz(-0.24438582) q[2];
sx q[2];
rz(-0.99437154) q[2];
sx q[2];
rz(-2.852358) q[2];
rz(-2.0770155) q[3];
sx q[3];
rz(-1.5983221) q[3];
sx q[3];
rz(1.7479001) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
