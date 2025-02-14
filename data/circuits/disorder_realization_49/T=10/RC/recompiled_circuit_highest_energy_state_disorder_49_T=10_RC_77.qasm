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
rz(-2.7558514) q[0];
sx q[0];
rz(-2.1585611) q[0];
sx q[0];
rz(-2.5842343) q[0];
rz(-1.9269257) q[1];
sx q[1];
rz(-2.2872556) q[1];
sx q[1];
rz(-2.1995423) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.048006213) q[0];
sx q[0];
rz(-0.24106179) q[0];
sx q[0];
rz(2.0895435) q[0];
x q[1];
rz(0.6379794) q[2];
sx q[2];
rz(-2.061741) q[2];
sx q[2];
rz(-0.52396905) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.57191932) q[1];
sx q[1];
rz(-1.5077356) q[1];
sx q[1];
rz(-2.75534) q[1];
x q[2];
rz(-0.73318847) q[3];
sx q[3];
rz(-0.7546784) q[3];
sx q[3];
rz(1.4316991) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.72630924) q[2];
sx q[2];
rz(-1.3741263) q[2];
sx q[2];
rz(1.7706002) q[2];
rz(2.3224984) q[3];
sx q[3];
rz(-2.2947125) q[3];
sx q[3];
rz(-3.0313361) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6743728) q[0];
sx q[0];
rz(-1.9085566) q[0];
sx q[0];
rz(-0.20587532) q[0];
rz(1.8241833) q[1];
sx q[1];
rz(-1.859788) q[1];
sx q[1];
rz(2.6726216) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.38843537) q[0];
sx q[0];
rz(-2.4860365) q[0];
sx q[0];
rz(0.038319328) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0360175) q[2];
sx q[2];
rz(-1.7649467) q[2];
sx q[2];
rz(2.126136) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.81509619) q[1];
sx q[1];
rz(-2.3760258) q[1];
sx q[1];
rz(0.86242843) q[1];
x q[2];
rz(-0.75403611) q[3];
sx q[3];
rz(-2.0181351) q[3];
sx q[3];
rz(-3.0424117) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.59490243) q[2];
sx q[2];
rz(-0.83260834) q[2];
sx q[2];
rz(-0.62758315) q[2];
rz(-1.0707431) q[3];
sx q[3];
rz(-0.50361931) q[3];
sx q[3];
rz(0.32881769) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-1.4048432) q[0];
sx q[0];
rz(-0.81031814) q[0];
sx q[0];
rz(0.78848439) q[0];
rz(-2.5704747) q[1];
sx q[1];
rz(-1.8126789) q[1];
sx q[1];
rz(1.4459389) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3841032) q[0];
sx q[0];
rz(-1.5575842) q[0];
sx q[0];
rz(1.3020975) q[0];
rz(-pi) q[1];
x q[1];
rz(0.064632434) q[2];
sx q[2];
rz(-0.53896133) q[2];
sx q[2];
rz(0.85899437) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.4169374) q[1];
sx q[1];
rz(-1.8677399) q[1];
sx q[1];
rz(1.4663803) q[1];
x q[2];
rz(0.57351968) q[3];
sx q[3];
rz(-1.2723426) q[3];
sx q[3];
rz(2.2645166) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.0420405) q[2];
sx q[2];
rz(-2.6991762) q[2];
sx q[2];
rz(-0.37910795) q[2];
rz(1.3742617) q[3];
sx q[3];
rz(-2.1702424) q[3];
sx q[3];
rz(2.1578535) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(-2.9336201) q[0];
sx q[0];
rz(-0.75478983) q[0];
sx q[0];
rz(2.2291613) q[0];
rz(-1.3937996) q[1];
sx q[1];
rz(-0.50519609) q[1];
sx q[1];
rz(0.30232731) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.94255208) q[0];
sx q[0];
rz(-1.4438757) q[0];
sx q[0];
rz(1.5104483) q[0];
rz(-pi) q[1];
rz(2.7458397) q[2];
sx q[2];
rz(-2.2456093) q[2];
sx q[2];
rz(-0.84854919) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.15790882) q[1];
sx q[1];
rz(-1.6140466) q[1];
sx q[1];
rz(-1.7801579) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4563645) q[3];
sx q[3];
rz(-0.99894612) q[3];
sx q[3];
rz(-2.2603007) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.6572774) q[2];
sx q[2];
rz(-0.70168287) q[2];
sx q[2];
rz(-0.62937984) q[2];
rz(-1.4194007) q[3];
sx q[3];
rz(-1.281176) q[3];
sx q[3];
rz(0.70845503) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.092954271) q[0];
sx q[0];
rz(-2.104367) q[0];
sx q[0];
rz(-1.3193489) q[0];
rz(-2.0535779) q[1];
sx q[1];
rz(-1.8698147) q[1];
sx q[1];
rz(-1.6768657) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7438904) q[0];
sx q[0];
rz(-2.4410998) q[0];
sx q[0];
rz(-1.6206436) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6032463) q[2];
sx q[2];
rz(-2.2032149) q[2];
sx q[2];
rz(2.4995952) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.3114918) q[1];
sx q[1];
rz(-1.3813918) q[1];
sx q[1];
rz(0.34475355) q[1];
rz(-pi) q[2];
x q[2];
rz(0.3220114) q[3];
sx q[3];
rz(-0.63068855) q[3];
sx q[3];
rz(1.2104152) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.4488843) q[2];
sx q[2];
rz(-0.11881891) q[2];
sx q[2];
rz(2.9046655) q[2];
rz(-2.1610625) q[3];
sx q[3];
rz(-1.4292932) q[3];
sx q[3];
rz(2.197649) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15810814) q[0];
sx q[0];
rz(-3.0321002) q[0];
sx q[0];
rz(-2.9964301) q[0];
rz(-1.521184) q[1];
sx q[1];
rz(-0.79997921) q[1];
sx q[1];
rz(0.0072341166) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1370498) q[0];
sx q[0];
rz(-0.88456735) q[0];
sx q[0];
rz(0.30651019) q[0];
x q[1];
rz(2.4154948) q[2];
sx q[2];
rz(-2.0019922) q[2];
sx q[2];
rz(-2.3027755) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.045102677) q[1];
sx q[1];
rz(-1.4670925) q[1];
sx q[1];
rz(1.3357049) q[1];
rz(-2.9754753) q[3];
sx q[3];
rz(-1.542585) q[3];
sx q[3];
rz(1.3835554) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.89473692) q[2];
sx q[2];
rz(-1.9523018) q[2];
sx q[2];
rz(-2.5640633) q[2];
rz(-1.9891116) q[3];
sx q[3];
rz(-0.58497506) q[3];
sx q[3];
rz(-1.729689) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5685527) q[0];
sx q[0];
rz(-2.9459406) q[0];
sx q[0];
rz(2.0797119) q[0];
rz(1.3878239) q[1];
sx q[1];
rz(-1.2304708) q[1];
sx q[1];
rz(-1.6434297) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.8138431) q[0];
sx q[0];
rz(-0.44506207) q[0];
sx q[0];
rz(2.829354) q[0];
rz(-2.4961124) q[2];
sx q[2];
rz(-1.8374075) q[2];
sx q[2];
rz(1.4652166) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.0688144) q[1];
sx q[1];
rz(-2.1238951) q[1];
sx q[1];
rz(2.8760853) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.56363799) q[3];
sx q[3];
rz(-0.56122128) q[3];
sx q[3];
rz(-1.3432251) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.7428703) q[2];
sx q[2];
rz(-3.0169432) q[2];
sx q[2];
rz(-2.2275662) q[2];
rz(0.75383178) q[3];
sx q[3];
rz(-1.9099312) q[3];
sx q[3];
rz(-0.45647538) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1330426) q[0];
sx q[0];
rz(-2.8164016) q[0];
sx q[0];
rz(-0.010566674) q[0];
rz(1.0039302) q[1];
sx q[1];
rz(-1.4164475) q[1];
sx q[1];
rz(-0.038987003) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4083791) q[0];
sx q[0];
rz(-0.36535242) q[0];
sx q[0];
rz(1.062458) q[0];
rz(-pi) q[1];
rz(3.1124074) q[2];
sx q[2];
rz(-0.31692908) q[2];
sx q[2];
rz(0.91561299) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.2891043) q[1];
sx q[1];
rz(-0.90910289) q[1];
sx q[1];
rz(1.0836592) q[1];
rz(-pi) q[2];
rz(-1.6603062) q[3];
sx q[3];
rz(-0.90965547) q[3];
sx q[3];
rz(-1.6474219) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.5799334) q[2];
sx q[2];
rz(-1.9114405) q[2];
sx q[2];
rz(-0.099543355) q[2];
rz(1.8482515) q[3];
sx q[3];
rz(-1.8563396) q[3];
sx q[3];
rz(0.37555638) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.67000166) q[0];
sx q[0];
rz(-1.3674068) q[0];
sx q[0];
rz(-1.8294096) q[0];
rz(0.52681628) q[1];
sx q[1];
rz(-0.75378886) q[1];
sx q[1];
rz(0.83676338) q[1];
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
rz(-pi) q[1];
rz(2.5946765) q[2];
sx q[2];
rz(-1.0267057) q[2];
sx q[2];
rz(0.61887276) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.30984658) q[1];
sx q[1];
rz(-1.832331) q[1];
sx q[1];
rz(2.4900718) q[1];
x q[2];
rz(-0.9012192) q[3];
sx q[3];
rz(-0.57358426) q[3];
sx q[3];
rz(0.18665126) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.7709363) q[2];
sx q[2];
rz(-2.5281361) q[2];
sx q[2];
rz(1.9343617) q[2];
rz(1.3197445) q[3];
sx q[3];
rz(-1.5765669) q[3];
sx q[3];
rz(2.3453662) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0470444) q[0];
sx q[0];
rz(-2.6667509) q[0];
sx q[0];
rz(-2.4249401) q[0];
rz(1.5776177) q[1];
sx q[1];
rz(-1.3037222) q[1];
sx q[1];
rz(-0.61734739) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0475395) q[0];
sx q[0];
rz(-2.2729024) q[0];
sx q[0];
rz(-0.2245837) q[0];
x q[1];
rz(1.1025589) q[2];
sx q[2];
rz(-1.7044633) q[2];
sx q[2];
rz(0.70883162) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.52386284) q[1];
sx q[1];
rz(-0.81392127) q[1];
sx q[1];
rz(0.62266962) q[1];
rz(-pi) q[2];
rz(-1.9412597) q[3];
sx q[3];
rz(-2.4672567) q[3];
sx q[3];
rz(0.5023027) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.3295595) q[2];
sx q[2];
rz(-2.7541408) q[2];
sx q[2];
rz(2.6623902) q[2];
rz(-1.6241578) q[3];
sx q[3];
rz(-1.4045709) q[3];
sx q[3];
rz(1.6421389) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1351521) q[0];
sx q[0];
rz(-1.3545481) q[0];
sx q[0];
rz(-0.56726278) q[0];
rz(0.80815036) q[1];
sx q[1];
rz(-1.5222526) q[1];
sx q[1];
rz(1.6076988) q[1];
rz(-3.0099536) q[2];
sx q[2];
rz(-2.0870118) q[2];
sx q[2];
rz(1.7160838) q[2];
rz(-0.2740673) q[3];
sx q[3];
rz(-2.3134939) q[3];
sx q[3];
rz(2.2753711) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
