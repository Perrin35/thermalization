OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.79787624) q[0];
sx q[0];
rz(2.1581082) q[0];
sx q[0];
rz(9.6165514) q[0];
rz(-2.6311488) q[1];
sx q[1];
rz(-1.3671083) q[1];
sx q[1];
rz(1.3956611) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6073942) q[0];
sx q[0];
rz(-1.5572131) q[0];
sx q[0];
rz(-3.109038) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.67740719) q[2];
sx q[2];
rz(-0.29149326) q[2];
sx q[2];
rz(-0.10805932) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.6916445) q[1];
sx q[1];
rz(-0.80300036) q[1];
sx q[1];
rz(1.2951485) q[1];
rz(1.1825425) q[3];
sx q[3];
rz(-0.40273977) q[3];
sx q[3];
rz(0.54408636) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.50769606) q[2];
sx q[2];
rz(-0.32749978) q[2];
sx q[2];
rz(0.92602777) q[2];
rz(1.5121459) q[3];
sx q[3];
rz(-0.67073268) q[3];
sx q[3];
rz(1.9984455) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.73285) q[0];
sx q[0];
rz(-2.4091305) q[0];
sx q[0];
rz(-0.23250411) q[0];
rz(2.0022424) q[1];
sx q[1];
rz(-1.573223) q[1];
sx q[1];
rz(0.92612902) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6059604) q[0];
sx q[0];
rz(-1.7754007) q[0];
sx q[0];
rz(1.3406517) q[0];
x q[1];
rz(-1.4732292) q[2];
sx q[2];
rz(-1.2446523) q[2];
sx q[2];
rz(-0.17352428) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.88484832) q[1];
sx q[1];
rz(-1.7769281) q[1];
sx q[1];
rz(0.54289674) q[1];
rz(-2.0235463) q[3];
sx q[3];
rz(-2.0914075) q[3];
sx q[3];
rz(-1.758616) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.53655475) q[2];
sx q[2];
rz(-1.6190517) q[2];
sx q[2];
rz(0.79263318) q[2];
rz(2.7301181) q[3];
sx q[3];
rz(-0.94235197) q[3];
sx q[3];
rz(0.80504942) q[3];
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
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8234392) q[0];
sx q[0];
rz(-1.0018438) q[0];
sx q[0];
rz(0.26082984) q[0];
rz(-2.8584495) q[1];
sx q[1];
rz(-1.4500376) q[1];
sx q[1];
rz(2.0859437) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5427534) q[0];
sx q[0];
rz(-0.93523798) q[0];
sx q[0];
rz(2.928614) q[0];
rz(0.0051202444) q[2];
sx q[2];
rz(-2.5084506) q[2];
sx q[2];
rz(2.1767165) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.187589) q[1];
sx q[1];
rz(-2.2936645) q[1];
sx q[1];
rz(-1.4299117) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6840213) q[3];
sx q[3];
rz(-2.6432496) q[3];
sx q[3];
rz(-1.6087974) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.61567125) q[2];
sx q[2];
rz(-1.0893704) q[2];
sx q[2];
rz(-1.2582568) q[2];
rz(2.1028171) q[3];
sx q[3];
rz(-0.98826161) q[3];
sx q[3];
rz(-1.725215) q[3];
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
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0106169) q[0];
sx q[0];
rz(-1.5357786) q[0];
sx q[0];
rz(-3.0220939) q[0];
rz(-0.15826982) q[1];
sx q[1];
rz(-2.6893078) q[1];
sx q[1];
rz(-1.4403042) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7410204) q[0];
sx q[0];
rz(-1.2109247) q[0];
sx q[0];
rz(0.82525702) q[0];
rz(0.079147804) q[2];
sx q[2];
rz(-1.2368349) q[2];
sx q[2];
rz(-0.755366) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.56934565) q[1];
sx q[1];
rz(-1.3075446) q[1];
sx q[1];
rz(-2.6514228) q[1];
rz(2.7915032) q[3];
sx q[3];
rz(-1.105827) q[3];
sx q[3];
rz(1.657682) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.4629472) q[2];
sx q[2];
rz(-0.1354278) q[2];
sx q[2];
rz(-2.7079008) q[2];
rz(2.4667451) q[3];
sx q[3];
rz(-1.0899455) q[3];
sx q[3];
rz(-1.1564144) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40815142) q[0];
sx q[0];
rz(-1.0557405) q[0];
sx q[0];
rz(-0.59930402) q[0];
rz(-2.377548) q[1];
sx q[1];
rz(-2.1115477) q[1];
sx q[1];
rz(-2.5291671) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9117994) q[0];
sx q[0];
rz(-0.62055885) q[0];
sx q[0];
rz(-1.3156652) q[0];
rz(-0.71686042) q[2];
sx q[2];
rz(-1.5713501) q[2];
sx q[2];
rz(3.1231511) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.45578411) q[1];
sx q[1];
rz(-2.1791885) q[1];
sx q[1];
rz(-1.7753106) q[1];
rz(-pi) q[2];
rz(-2.6153485) q[3];
sx q[3];
rz(-1.11107) q[3];
sx q[3];
rz(-2.3731874) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.4446438) q[2];
sx q[2];
rz(-0.79978839) q[2];
sx q[2];
rz(-2.6714228) q[2];
rz(2.7759975) q[3];
sx q[3];
rz(-1.9496893) q[3];
sx q[3];
rz(2.5308334) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.325901) q[0];
sx q[0];
rz(-0.79623643) q[0];
sx q[0];
rz(2.4796487) q[0];
rz(0.081261948) q[1];
sx q[1];
rz(-2.6632402) q[1];
sx q[1];
rz(0.14019664) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.88410035) q[0];
sx q[0];
rz(-2.0292768) q[0];
sx q[0];
rz(1.8152587) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.95023167) q[2];
sx q[2];
rz(-1.8355614) q[2];
sx q[2];
rz(-2.1193412) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.6684456) q[1];
sx q[1];
rz(-2.6700017) q[1];
sx q[1];
rz(-1.5525329) q[1];
x q[2];
rz(1.3883038) q[3];
sx q[3];
rz(-1.914822) q[3];
sx q[3];
rz(-1.9545912) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.7774272) q[2];
sx q[2];
rz(-1.825793) q[2];
sx q[2];
rz(2.8167456) q[2];
rz(-0.15803629) q[3];
sx q[3];
rz(-0.41301781) q[3];
sx q[3];
rz(-1.0154999) q[3];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.82191104) q[0];
sx q[0];
rz(-0.076229036) q[0];
sx q[0];
rz(2.2538189) q[0];
rz(2.7582788) q[1];
sx q[1];
rz(-0.8822459) q[1];
sx q[1];
rz(0.15957889) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4604546) q[0];
sx q[0];
rz(-0.41550203) q[0];
sx q[0];
rz(-2.2592708) q[0];
rz(-pi) q[1];
rz(2.1453484) q[2];
sx q[2];
rz(-2.0483082) q[2];
sx q[2];
rz(1.4487131) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.3899052) q[1];
sx q[1];
rz(-1.4426822) q[1];
sx q[1];
rz(0.37813487) q[1];
rz(0.28954808) q[3];
sx q[3];
rz(-1.411045) q[3];
sx q[3];
rz(-1.8600132) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.5923656) q[2];
sx q[2];
rz(-1.9840019) q[2];
sx q[2];
rz(1.7249736) q[2];
rz(1.8727411) q[3];
sx q[3];
rz(-0.78059355) q[3];
sx q[3];
rz(-0.95808539) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.45605993) q[0];
sx q[0];
rz(-2.0262418) q[0];
sx q[0];
rz(-2.2400895) q[0];
rz(-2.447336) q[1];
sx q[1];
rz(-1.6777104) q[1];
sx q[1];
rz(-1.136397) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0471056) q[0];
sx q[0];
rz(-1.4102954) q[0];
sx q[0];
rz(-3.1094527) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.55113422) q[2];
sx q[2];
rz(-1.2493842) q[2];
sx q[2];
rz(1.8773735) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.9688837) q[1];
sx q[1];
rz(-0.69689059) q[1];
sx q[1];
rz(1.642176) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4122167) q[3];
sx q[3];
rz(-1.1664026) q[3];
sx q[3];
rz(2.803363) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(3.0006492) q[2];
sx q[2];
rz(-1.2691754) q[2];
sx q[2];
rz(-2.3285274) q[2];
rz(-2.5874169) q[3];
sx q[3];
rz(-1.7289836) q[3];
sx q[3];
rz(0.36177844) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5799705) q[0];
sx q[0];
rz(-2.0273835) q[0];
sx q[0];
rz(0.17350523) q[0];
rz(-0.42501998) q[1];
sx q[1];
rz(-1.5360906) q[1];
sx q[1];
rz(0.82957155) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.94151992) q[0];
sx q[0];
rz(-2.3931112) q[0];
sx q[0];
rz(0.21416382) q[0];
rz(-0.7220975) q[2];
sx q[2];
rz(-3.0797662) q[2];
sx q[2];
rz(2.823148) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.049343892) q[1];
sx q[1];
rz(-2.4929765) q[1];
sx q[1];
rz(0.34347024) q[1];
rz(2.8052727) q[3];
sx q[3];
rz(-1.5097396) q[3];
sx q[3];
rz(1.2741736) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.35036626) q[2];
sx q[2];
rz(-1.841505) q[2];
sx q[2];
rz(1.680797) q[2];
rz(2.1469877) q[3];
sx q[3];
rz(-0.9959144) q[3];
sx q[3];
rz(-2.2554956) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9643672) q[0];
sx q[0];
rz(-2.565964) q[0];
sx q[0];
rz(3.0395569) q[0];
rz(-2.521934) q[1];
sx q[1];
rz(-1.6435868) q[1];
sx q[1];
rz(-2.5443351) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5383947) q[0];
sx q[0];
rz(-1.704633) q[0];
sx q[0];
rz(0.60026818) q[0];
rz(-pi) q[1];
rz(-2.7623457) q[2];
sx q[2];
rz(-0.37988099) q[2];
sx q[2];
rz(-2.3178562) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.44587943) q[1];
sx q[1];
rz(-1.5339601) q[1];
sx q[1];
rz(-2.1365154) q[1];
rz(-pi) q[2];
rz(0.47145505) q[3];
sx q[3];
rz(-1.5416036) q[3];
sx q[3];
rz(-1.4271133) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.9336046) q[2];
sx q[2];
rz(-0.36078578) q[2];
sx q[2];
rz(2.2130845) q[2];
rz(-1.1601296) q[3];
sx q[3];
rz(-1.4312294) q[3];
sx q[3];
rz(2.3742356) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.74850294) q[0];
sx q[0];
rz(-0.7001644) q[0];
sx q[0];
rz(-2.9099332) q[0];
rz(2.0139991) q[1];
sx q[1];
rz(-0.53378202) q[1];
sx q[1];
rz(3.1376874) q[1];
rz(0.071795287) q[2];
sx q[2];
rz(-0.89569246) q[2];
sx q[2];
rz(2.8099974) q[2];
rz(-2.2670408) q[3];
sx q[3];
rz(-1.8378432) q[3];
sx q[3];
rz(1.7624553) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
