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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7638057) q[0];
sx q[0];
rz(-1.2391234) q[0];
sx q[0];
rz(-0.8140696) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5231126) q[2];
sx q[2];
rz(-1.5850726) q[2];
sx q[2];
rz(2.2592777) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-3.0539581) q[1];
sx q[1];
rz(-0.9231418) q[1];
sx q[1];
rz(-2.1235076) q[1];
rz(-pi) q[2];
rz(-0.7865016) q[3];
sx q[3];
rz(-2.3701027) q[3];
sx q[3];
rz(-2.2685693) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.8636318) q[2];
sx q[2];
rz(-2.7148254) q[2];
sx q[2];
rz(-2.5521736) q[2];
rz(-2.3943118) q[3];
sx q[3];
rz(-3.0487479) q[3];
sx q[3];
rz(0.59982991) q[3];
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
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0624307) q[0];
sx q[0];
rz(-2.9099162) q[0];
sx q[0];
rz(2.2404501) q[0];
rz(2.1597553) q[1];
sx q[1];
rz(-0.58126175) q[1];
sx q[1];
rz(-2.1493886) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3311148) q[0];
sx q[0];
rz(-1.7077291) q[0];
sx q[0];
rz(0.94816072) q[0];
rz(0.16449931) q[2];
sx q[2];
rz(-2.3119011) q[2];
sx q[2];
rz(2.6890597) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.8507311) q[1];
sx q[1];
rz(-0.90576101) q[1];
sx q[1];
rz(-0.70677251) q[1];
x q[2];
rz(-0.76073356) q[3];
sx q[3];
rz(-2.2983716) q[3];
sx q[3];
rz(0.47970495) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.09484) q[2];
sx q[2];
rz(-0.53043008) q[2];
sx q[2];
rz(-0.025010427) q[2];
rz(0.95995861) q[3];
sx q[3];
rz(-0.046006087) q[3];
sx q[3];
rz(3.0733601) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9571824) q[0];
sx q[0];
rz(-3.130641) q[0];
sx q[0];
rz(0.72407323) q[0];
rz(0.11101668) q[1];
sx q[1];
rz(-0.60581517) q[1];
sx q[1];
rz(-0.012880005) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89319481) q[0];
sx q[0];
rz(-1.3583368) q[0];
sx q[0];
rz(-0.080470632) q[0];
x q[1];
rz(-2.951163) q[2];
sx q[2];
rz(-1.7594751) q[2];
sx q[2];
rz(-1.3946472) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.15808039) q[1];
sx q[1];
rz(-1.8077824) q[1];
sx q[1];
rz(-1.7416341) q[1];
rz(-pi) q[2];
rz(2.0184085) q[3];
sx q[3];
rz(-2.4035998) q[3];
sx q[3];
rz(0.18692218) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.8956902) q[2];
sx q[2];
rz(-2.4195713) q[2];
sx q[2];
rz(-0.32430172) q[2];
rz(-2.6445828) q[3];
sx q[3];
rz(-2.9091166) q[3];
sx q[3];
rz(-2.9089109) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.67320353) q[0];
sx q[0];
rz(-2.9717746) q[0];
sx q[0];
rz(2.66535) q[0];
rz(0.4854804) q[1];
sx q[1];
rz(-0.49518934) q[1];
sx q[1];
rz(2.9229497) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9354585) q[0];
sx q[0];
rz(-1.1674321) q[0];
sx q[0];
rz(-2.8950188) q[0];
rz(1.5584458) q[2];
sx q[2];
rz(-0.84721476) q[2];
sx q[2];
rz(0.45958322) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.051603457) q[1];
sx q[1];
rz(-1.0538424) q[1];
sx q[1];
rz(2.0614784) q[1];
rz(-0.63204445) q[3];
sx q[3];
rz(-0.38292745) q[3];
sx q[3];
rz(-0.51198375) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.3525047) q[2];
sx q[2];
rz(-2.3829491) q[2];
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
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
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
rz(-2.5649696) q[0];
sx q[0];
rz(-2.7374856) q[0];
sx q[0];
rz(0.47082666) q[0];
rz(0.91104031) q[1];
sx q[1];
rz(-2.7016579) q[1];
sx q[1];
rz(0.43112531) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.39931) q[0];
sx q[0];
rz(-2.3547908) q[0];
sx q[0];
rz(1.6086701) q[0];
x q[1];
rz(2.2041915) q[2];
sx q[2];
rz(-1.3654815) q[2];
sx q[2];
rz(0.61287731) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.4404802) q[1];
sx q[1];
rz(-1.639469) q[1];
sx q[1];
rz(3.1140226) q[1];
rz(-0.60089941) q[3];
sx q[3];
rz(-0.24293262) q[3];
sx q[3];
rz(3.0593135) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.8016781) q[2];
sx q[2];
rz(-0.0047923294) q[2];
sx q[2];
rz(2.8002296) q[2];
rz(-0.3723799) q[3];
sx q[3];
rz(-0.63569331) q[3];
sx q[3];
rz(-2.671833) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
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
rz(-2.3434134) q[1];
sx q[1];
rz(0.72365671) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.22985499) q[0];
sx q[0];
rz(-0.17979547) q[0];
sx q[0];
rz(1.3151965) q[0];
x q[1];
rz(-1.0704109) q[2];
sx q[2];
rz(-2.7636409) q[2];
sx q[2];
rz(0.17411451) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.78711975) q[1];
sx q[1];
rz(-0.90098375) q[1];
sx q[1];
rz(-0.77490999) q[1];
x q[2];
rz(2.0615929) q[3];
sx q[3];
rz(-0.43269581) q[3];
sx q[3];
rz(-0.015263488) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3896997) q[2];
sx q[2];
rz(-2.5395826) q[2];
sx q[2];
rz(-2.4683118) q[2];
rz(2.8781387) q[3];
sx q[3];
rz(-2.6665688) q[3];
sx q[3];
rz(2.8365005) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.70169705) q[0];
sx q[0];
rz(-2.3450527) q[0];
sx q[0];
rz(-0.51629603) q[0];
rz(-2.3856178) q[1];
sx q[1];
rz(-0.95840234) q[1];
sx q[1];
rz(2.7730673) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4936016) q[0];
sx q[0];
rz(-2.2514935) q[0];
sx q[0];
rz(2.2590911) q[0];
rz(-pi) q[1];
rz(1.5318977) q[2];
sx q[2];
rz(-2.2752011) q[2];
sx q[2];
rz(-0.6608085) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.8766899) q[1];
sx q[1];
rz(-1.4646834) q[1];
sx q[1];
rz(-1.7404895) q[1];
x q[2];
rz(2.1462899) q[3];
sx q[3];
rz(-0.80484521) q[3];
sx q[3];
rz(-0.14508776) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.64510173) q[2];
sx q[2];
rz(-3.0848905) q[2];
sx q[2];
rz(-0.48218316) q[2];
rz(2.9218946) q[3];
sx q[3];
rz(-2.4441661) q[3];
sx q[3];
rz(-0.68035948) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7242303) q[0];
sx q[0];
rz(-0.14886947) q[0];
sx q[0];
rz(-2.505488) q[0];
rz(-0.96027374) q[1];
sx q[1];
rz(-0.93820131) q[1];
sx q[1];
rz(0.87463921) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3134523) q[0];
sx q[0];
rz(-0.16022542) q[0];
sx q[0];
rz(-1.942722) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6804439) q[2];
sx q[2];
rz(-1.6911518) q[2];
sx q[2];
rz(0.14794825) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.0046152701) q[1];
sx q[1];
rz(-2.5680313) q[1];
sx q[1];
rz(-0.80570813) q[1];
rz(-pi) q[2];
rz(-2.465807) q[3];
sx q[3];
rz(-0.99832557) q[3];
sx q[3];
rz(1.0990717) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.26802289) q[2];
sx q[2];
rz(-0.41646725) q[2];
sx q[2];
rz(0.91656172) q[2];
rz(-2.364184) q[3];
sx q[3];
rz(-2.58367) q[3];
sx q[3];
rz(-2.6072445) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1189608) q[0];
sx q[0];
rz(-0.73771483) q[0];
sx q[0];
rz(-2.90888) q[0];
rz(-0.64838707) q[1];
sx q[1];
rz(-2.5189923) q[1];
sx q[1];
rz(2.5737305) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.083160087) q[0];
sx q[0];
rz(-2.5867903) q[0];
sx q[0];
rz(-2.2834999) q[0];
rz(2.2350506) q[2];
sx q[2];
rz(-2.8010606) q[2];
sx q[2];
rz(-2.1730925) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.7306914) q[1];
sx q[1];
rz(-1.3679549) q[1];
sx q[1];
rz(-0.14932015) q[1];
x q[2];
rz(0.73091032) q[3];
sx q[3];
rz(-0.95940351) q[3];
sx q[3];
rz(2.9711593) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.76967543) q[2];
sx q[2];
rz(-0.19352517) q[2];
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
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(0.929739) q[0];
sx q[0];
rz(-3.0775098) q[0];
sx q[0];
rz(0.1317568) q[0];
rz(0.53648221) q[1];
sx q[1];
rz(-0.15833144) q[1];
sx q[1];
rz(-0.90824711) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4038178) q[0];
sx q[0];
rz(-1.6006915) q[0];
sx q[0];
rz(-2.2653518) q[0];
x q[1];
rz(1.0824049) q[2];
sx q[2];
rz(-0.51967144) q[2];
sx q[2];
rz(-2.1888417) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.51146347) q[1];
sx q[1];
rz(-0.88246934) q[1];
sx q[1];
rz(0.87911112) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5628917) q[3];
sx q[3];
rz(-0.22592446) q[3];
sx q[3];
rz(-2.853012) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.50916719) q[2];
sx q[2];
rz(-0.91370344) q[2];
sx q[2];
rz(-2.8636279) q[2];
rz(0.5667423) q[3];
sx q[3];
rz(-0.20213474) q[3];
sx q[3];
rz(2.22866) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
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
rz(0.81000281) q[0];
sx q[0];
rz(-1.4053874) q[0];
sx q[0];
rz(-0.94595861) q[0];
rz(2.6592061) q[1];
sx q[1];
rz(-1.9079897) q[1];
sx q[1];
rz(2.447396) q[1];
rz(0.24438582) q[2];
sx q[2];
rz(-2.1472211) q[2];
sx q[2];
rz(0.28923464) q[2];
rz(0.031470555) q[3];
sx q[3];
rz(-2.0768055) q[3];
sx q[3];
rz(0.16184645) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
