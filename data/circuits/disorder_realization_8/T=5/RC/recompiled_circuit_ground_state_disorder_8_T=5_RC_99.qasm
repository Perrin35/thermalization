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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.64144999) q[0];
sx q[0];
rz(-0.035273835) q[0];
sx q[0];
rz(-2.7462237) q[0];
x q[1];
rz(0.22966603) q[2];
sx q[2];
rz(-1.3896835) q[2];
sx q[2];
rz(-2.119198) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.073133999) q[1];
sx q[1];
rz(-1.3737118) q[1];
sx q[1];
rz(2.3545594) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9816876) q[3];
sx q[3];
rz(-1.9420338) q[3];
sx q[3];
rz(-2.1791636) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.6338966) q[2];
sx q[2];
rz(-0.32749978) q[2];
sx q[2];
rz(0.92602777) q[2];
rz(1.5121459) q[3];
sx q[3];
rz(-0.67073268) q[3];
sx q[3];
rz(-1.1431471) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.73285) q[0];
sx q[0];
rz(-2.4091305) q[0];
sx q[0];
rz(0.23250411) q[0];
rz(1.1393503) q[1];
sx q[1];
rz(-1.5683697) q[1];
sx q[1];
rz(0.92612902) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8208198) q[0];
sx q[0];
rz(-0.30673393) q[0];
sx q[0];
rz(-2.308918) q[0];
rz(-pi) q[1];
rz(2.813999) q[2];
sx q[2];
rz(-1.4783876) q[2];
sx q[2];
rz(-1.7756697) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.7826816) q[1];
sx q[1];
rz(-2.5645683) q[1];
sx q[1];
rz(-0.38459528) q[1];
x q[2];
rz(0.56761928) q[3];
sx q[3];
rz(-1.960037) q[3];
sx q[3];
rz(-0.4252227) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.6050379) q[2];
sx q[2];
rz(-1.522541) q[2];
sx q[2];
rz(-2.3489595) q[2];
rz(-0.41147453) q[3];
sx q[3];
rz(-2.1992407) q[3];
sx q[3];
rz(2.3365432) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3181535) q[0];
sx q[0];
rz(-2.1397488) q[0];
sx q[0];
rz(-0.26082984) q[0];
rz(-0.28314319) q[1];
sx q[1];
rz(-1.6915551) q[1];
sx q[1];
rz(-1.0556489) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59883927) q[0];
sx q[0];
rz(-0.93523798) q[0];
sx q[0];
rz(0.21297867) q[0];
rz(-0.63313578) q[2];
sx q[2];
rz(-1.5738259) q[2];
sx q[2];
rz(2.5315447) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.3987777) q[1];
sx q[1];
rz(-2.4075634) q[1];
sx q[1];
rz(0.15785288) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0664472) q[3];
sx q[3];
rz(-1.6248253) q[3];
sx q[3];
rz(-3.0040405) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.61567125) q[2];
sx q[2];
rz(-2.0522223) q[2];
sx q[2];
rz(1.2582568) q[2];
rz(-2.1028171) q[3];
sx q[3];
rz(-2.153331) q[3];
sx q[3];
rz(-1.725215) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
sx q[3];
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
rz(-0.1309758) q[0];
sx q[0];
rz(-1.6058141) q[0];
sx q[0];
rz(-3.0220939) q[0];
rz(0.15826982) q[1];
sx q[1];
rz(-0.4522849) q[1];
sx q[1];
rz(1.7012885) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.8058469) q[0];
sx q[0];
rz(-2.3290538) q[0];
sx q[0];
rz(-1.0643908) q[0];
rz(-pi) q[1];
x q[1];
rz(0.079147804) q[2];
sx q[2];
rz(-1.9047577) q[2];
sx q[2];
rz(-2.3862267) q[2];
x q[3];
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
rz(-pi) q[2];
x q[2];
rz(-2.7915032) q[3];
sx q[3];
rz(-1.105827) q[3];
sx q[3];
rz(1.4839107) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.67864546) q[2];
sx q[2];
rz(-0.1354278) q[2];
sx q[2];
rz(-0.43369183) q[2];
rz(-2.4667451) q[3];
sx q[3];
rz(-2.0516472) q[3];
sx q[3];
rz(1.9851782) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7334412) q[0];
sx q[0];
rz(-2.0858522) q[0];
sx q[0];
rz(-2.5422886) q[0];
rz(2.377548) q[1];
sx q[1];
rz(-2.1115477) q[1];
sx q[1];
rz(2.5291671) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9117994) q[0];
sx q[0];
rz(-2.5210338) q[0];
sx q[0];
rz(1.3156652) q[0];
rz(-pi) q[1];
rz(0.00084288518) q[2];
sx q[2];
rz(-0.71686059) q[2];
sx q[2];
rz(1.5529902) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.6858085) q[1];
sx q[1];
rz(-0.9624042) q[1];
sx q[1];
rz(-1.7753106) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.52624412) q[3];
sx q[3];
rz(-2.0305227) q[3];
sx q[3];
rz(0.7684052) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.4446438) q[2];
sx q[2];
rz(-0.79978839) q[2];
sx q[2];
rz(-0.4701699) q[2];
rz(2.7759975) q[3];
sx q[3];
rz(-1.1919034) q[3];
sx q[3];
rz(-2.5308334) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-2.325901) q[0];
sx q[0];
rz(-2.3453562) q[0];
sx q[0];
rz(0.66194397) q[0];
rz(-3.0603307) q[1];
sx q[1];
rz(-0.47835246) q[1];
sx q[1];
rz(-0.14019664) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37085846) q[0];
sx q[0];
rz(-2.6261289) q[0];
sx q[0];
rz(2.6856642) q[0];
x q[1];
rz(2.0070932) q[2];
sx q[2];
rz(-0.66777523) q[2];
sx q[2];
rz(-2.2420924) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.4936476) q[1];
sx q[1];
rz(-2.0423023) q[1];
sx q[1];
rz(3.1322797) q[1];
rz(1.3883038) q[3];
sx q[3];
rz(-1.914822) q[3];
sx q[3];
rz(1.1870015) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.36416546) q[2];
sx q[2];
rz(-1.3157996) q[2];
sx q[2];
rz(0.3248471) q[2];
rz(2.9835564) q[3];
sx q[3];
rz(-0.41301781) q[3];
sx q[3];
rz(-1.0154999) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.82191104) q[0];
sx q[0];
rz(-0.076229036) q[0];
sx q[0];
rz(-2.2538189) q[0];
rz(0.38331389) q[1];
sx q[1];
rz(-2.2593468) q[1];
sx q[1];
rz(-2.9820138) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.051285714) q[0];
sx q[0];
rz(-1.2538099) q[0];
sx q[0];
rz(-0.27329926) q[0];
rz(-pi) q[1];
rz(0.99624421) q[2];
sx q[2];
rz(-1.0932845) q[2];
sx q[2];
rz(-1.6928796) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.9099906) q[1];
sx q[1];
rz(-1.195915) q[1];
sx q[1];
rz(1.433062) q[1];
rz(-pi) q[2];
rz(2.6278213) q[3];
sx q[3];
rz(-2.8119866) q[3];
sx q[3];
rz(0.20121516) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.5923656) q[2];
sx q[2];
rz(-1.9840019) q[2];
sx q[2];
rz(-1.4166191) q[2];
rz(-1.2688515) q[3];
sx q[3];
rz(-0.78059355) q[3];
sx q[3];
rz(2.1835073) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6855327) q[0];
sx q[0];
rz(-2.0262418) q[0];
sx q[0];
rz(-2.2400895) q[0];
rz(2.447336) q[1];
sx q[1];
rz(-1.6777104) q[1];
sx q[1];
rz(1.136397) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6601453) q[0];
sx q[0];
rz(-1.602523) q[0];
sx q[0];
rz(1.4102139) q[0];
x q[1];
rz(-1.9433697) q[2];
sx q[2];
rz(-2.0907474) q[2];
sx q[2];
rz(0.49836788) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.7982805) q[1];
sx q[1];
rz(-1.5250051) q[1];
sx q[1];
rz(-2.2664323) q[1];
rz(-0.57089048) q[3];
sx q[3];
rz(-0.81557214) q[3];
sx q[3];
rz(-2.3237128) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.14094341) q[2];
sx q[2];
rz(-1.8724172) q[2];
sx q[2];
rz(0.81306523) q[2];
rz(0.55417577) q[3];
sx q[3];
rz(-1.412609) q[3];
sx q[3];
rz(-0.36177844) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5799705) q[0];
sx q[0];
rz(-1.1142092) q[0];
sx q[0];
rz(-2.9680874) q[0];
rz(0.42501998) q[1];
sx q[1];
rz(-1.5360906) q[1];
sx q[1];
rz(-0.82957155) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2000727) q[0];
sx q[0];
rz(-2.3931112) q[0];
sx q[0];
rz(-0.21416382) q[0];
rz(-pi) q[1];
rz(-1.5299022) q[2];
sx q[2];
rz(-1.5244134) q[2];
sx q[2];
rz(-2.1001018) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.47118716) q[1];
sx q[1];
rz(-0.96574942) q[1];
sx q[1];
rz(-1.3208645) q[1];
rz(-2.8052727) q[3];
sx q[3];
rz(-1.631853) q[3];
sx q[3];
rz(-1.8674191) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.7912264) q[2];
sx q[2];
rz(-1.3000877) q[2];
sx q[2];
rz(1.680797) q[2];
rz(2.1469877) q[3];
sx q[3];
rz(-2.1456783) q[3];
sx q[3];
rz(2.2554956) q[3];
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
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9643672) q[0];
sx q[0];
rz(-2.565964) q[0];
sx q[0];
rz(3.0395569) q[0];
rz(-0.61965865) q[1];
sx q[1];
rz(-1.6435868) q[1];
sx q[1];
rz(2.5443351) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0181018) q[0];
sx q[0];
rz(-2.1649556) q[0];
sx q[0];
rz(-1.4090562) q[0];
x q[1];
rz(0.3551746) q[2];
sx q[2];
rz(-1.7085129) q[2];
sx q[2];
rz(-2.0400408) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.44587943) q[1];
sx q[1];
rz(-1.6076325) q[1];
sx q[1];
rz(-2.1365154) q[1];
rz(-pi) q[2];
x q[2];
rz(0.47145505) q[3];
sx q[3];
rz(-1.5416036) q[3];
sx q[3];
rz(-1.4271133) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.20798802) q[2];
sx q[2];
rz(-2.7808069) q[2];
sx q[2];
rz(0.92850816) q[2];
rz(1.9814631) q[3];
sx q[3];
rz(-1.7103633) q[3];
sx q[3];
rz(-2.3742356) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.74850294) q[0];
sx q[0];
rz(-2.4414283) q[0];
sx q[0];
rz(0.23165942) q[0];
rz(-1.1275935) q[1];
sx q[1];
rz(-0.53378202) q[1];
sx q[1];
rz(3.1376874) q[1];
rz(-1.660158) q[2];
sx q[2];
rz(-0.67831525) q[2];
sx q[2];
rz(-0.44616551) q[2];
rz(-1.1675904) q[3];
sx q[3];
rz(-0.73763631) q[3];
sx q[3];
rz(0.49745001) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
