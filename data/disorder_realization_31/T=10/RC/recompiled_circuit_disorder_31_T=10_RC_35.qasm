OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.6150317) q[0];
sx q[0];
rz(-0.57305133) q[0];
sx q[0];
rz(-2.2990062) q[0];
rz(2.1057582) q[1];
sx q[1];
rz(8.3254568) q[1];
sx q[1];
rz(7.96666) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6095088) q[0];
sx q[0];
rz(-2.3713787) q[0];
sx q[0];
rz(0.30755933) q[0];
rz(0.027878472) q[2];
sx q[2];
rz(-2.5275702) q[2];
sx q[2];
rz(2.8369454) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.7588501) q[1];
sx q[1];
rz(-1.1945801) q[1];
sx q[1];
rz(-1.7099027) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9620142) q[3];
sx q[3];
rz(-0.60185963) q[3];
sx q[3];
rz(-1.7367401) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.6136916) q[2];
sx q[2];
rz(-1.0062904) q[2];
sx q[2];
rz(2.9620985) q[2];
rz(-1.2256631) q[3];
sx q[3];
rz(-1.3464728) q[3];
sx q[3];
rz(-2.3195482) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3935788) q[0];
sx q[0];
rz(-0.8809692) q[0];
sx q[0];
rz(2.8161312) q[0];
rz(-1.7851967) q[1];
sx q[1];
rz(-1.0486832) q[1];
sx q[1];
rz(-1.1546086) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3852859) q[0];
sx q[0];
rz(-0.025408832) q[0];
sx q[0];
rz(0.73861648) q[0];
rz(-pi) q[1];
rz(-2.1968368) q[2];
sx q[2];
rz(-1.895004) q[2];
sx q[2];
rz(-2.7446483) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.3719912) q[1];
sx q[1];
rz(-2.3725315) q[1];
sx q[1];
rz(-0.12188697) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1784337) q[3];
sx q[3];
rz(-0.31664407) q[3];
sx q[3];
rz(-2.7505927) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.6894199) q[2];
sx q[2];
rz(-1.2499115) q[2];
sx q[2];
rz(-0.88341218) q[2];
rz(2.6702821) q[3];
sx q[3];
rz(-1.703197) q[3];
sx q[3];
rz(2.3538891) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
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
rz(2.8283591) q[0];
sx q[0];
rz(-1.6468843) q[0];
sx q[0];
rz(-1.5154243) q[0];
rz(2.5405163) q[1];
sx q[1];
rz(-0.54769146) q[1];
sx q[1];
rz(1.0916969) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1900345) q[0];
sx q[0];
rz(-2.1265731) q[0];
sx q[0];
rz(-2.4843198) q[0];
rz(-pi) q[1];
rz(1.0513564) q[2];
sx q[2];
rz(-2.4161077) q[2];
sx q[2];
rz(1.5067593) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.8852383) q[1];
sx q[1];
rz(-1.9372205) q[1];
sx q[1];
rz(-2.3430235) q[1];
x q[2];
rz(-3.0136209) q[3];
sx q[3];
rz(-1.5068753) q[3];
sx q[3];
rz(1.909006) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.8213356) q[2];
sx q[2];
rz(-0.50575033) q[2];
sx q[2];
rz(0.88095218) q[2];
rz(1.3736003) q[3];
sx q[3];
rz(-1.526984) q[3];
sx q[3];
rz(-2.1239471) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3110733) q[0];
sx q[0];
rz(-1.3922465) q[0];
sx q[0];
rz(2.7048892) q[0];
rz(2.9084335) q[1];
sx q[1];
rz(-1.8893087) q[1];
sx q[1];
rz(-0.31035796) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.719602) q[0];
sx q[0];
rz(-1.4388226) q[0];
sx q[0];
rz(-0.3814401) q[0];
x q[1];
rz(2.9878346) q[2];
sx q[2];
rz(-2.4506844) q[2];
sx q[2];
rz(1.5916057) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.6314108) q[1];
sx q[1];
rz(-1.5893755) q[1];
sx q[1];
rz(1.2127962) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6944828) q[3];
sx q[3];
rz(-1.5354772) q[3];
sx q[3];
rz(0.24231054) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.13005304) q[2];
sx q[2];
rz(-0.71802846) q[2];
sx q[2];
rz(-2.0641573) q[2];
rz(-3.0854026) q[3];
sx q[3];
rz(-2.5037933) q[3];
sx q[3];
rz(-1.594054) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.189165) q[0];
sx q[0];
rz(-1.0452894) q[0];
sx q[0];
rz(-0.24965832) q[0];
rz(-1.5646308) q[1];
sx q[1];
rz(-0.77762929) q[1];
sx q[1];
rz(0.87019428) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10660431) q[0];
sx q[0];
rz(-0.37811324) q[0];
sx q[0];
rz(-2.5614221) q[0];
x q[1];
rz(0.69663163) q[2];
sx q[2];
rz(-1.7156892) q[2];
sx q[2];
rz(-2.2640413) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.89228499) q[1];
sx q[1];
rz(-1.882949) q[1];
sx q[1];
rz(-2.409163) q[1];
x q[2];
rz(1.8364041) q[3];
sx q[3];
rz(-0.57816539) q[3];
sx q[3];
rz(-2.2418914) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.27328086) q[2];
sx q[2];
rz(-1.8153278) q[2];
sx q[2];
rz(-2.4678521) q[2];
rz(2.8379748) q[3];
sx q[3];
rz(-1.9165336) q[3];
sx q[3];
rz(-1.822086) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7917787) q[0];
sx q[0];
rz(-0.93739167) q[0];
sx q[0];
rz(-2.8836024) q[0];
rz(0.42516431) q[1];
sx q[1];
rz(-0.95562569) q[1];
sx q[1];
rz(-1.649883) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8269577) q[0];
sx q[0];
rz(-1.1420113) q[0];
sx q[0];
rz(1.6786806) q[0];
rz(-0.50478023) q[2];
sx q[2];
rz(-0.74540388) q[2];
sx q[2];
rz(-2.8842852) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.7548435) q[1];
sx q[1];
rz(-1.0068839) q[1];
sx q[1];
rz(0.90653231) q[1];
rz(1.2061938) q[3];
sx q[3];
rz(-1.4758665) q[3];
sx q[3];
rz(-1.7942384) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.1288746) q[2];
sx q[2];
rz(-0.96367633) q[2];
sx q[2];
rz(-3.0498665) q[2];
rz(-0.84364676) q[3];
sx q[3];
rz(-2.1618312) q[3];
sx q[3];
rz(0.89404026) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3180852) q[0];
sx q[0];
rz(-1.894269) q[0];
sx q[0];
rz(-2.7303625) q[0];
rz(0.86589083) q[1];
sx q[1];
rz(-0.31232467) q[1];
sx q[1];
rz(-0.033989865) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.8695553) q[0];
sx q[0];
rz(-2.3453658) q[0];
sx q[0];
rz(1.7982593) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2074239) q[2];
sx q[2];
rz(-2.0206497) q[2];
sx q[2];
rz(2.5069782) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.46376343) q[1];
sx q[1];
rz(-2.4619108) q[1];
sx q[1];
rz(-1.4903279) q[1];
rz(-pi) q[2];
rz(-1.7185228) q[3];
sx q[3];
rz(-0.96102321) q[3];
sx q[3];
rz(-0.33965286) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.6035446) q[2];
sx q[2];
rz(-0.59843439) q[2];
sx q[2];
rz(2.2650488) q[2];
rz(-2.792568) q[3];
sx q[3];
rz(-1.2004431) q[3];
sx q[3];
rz(-2.9984737) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2440764) q[0];
sx q[0];
rz(-1.4325457) q[0];
sx q[0];
rz(-2.7476655) q[0];
rz(-2.774033) q[1];
sx q[1];
rz(-1.7575248) q[1];
sx q[1];
rz(1.6961018) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.60364265) q[0];
sx q[0];
rz(-2.0041487) q[0];
sx q[0];
rz(0.005714697) q[0];
rz(0.58480279) q[2];
sx q[2];
rz(-2.1091166) q[2];
sx q[2];
rz(-0.98758299) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.1576924) q[1];
sx q[1];
rz(-2.1196113) q[1];
sx q[1];
rz(-1.3859315) q[1];
x q[2];
rz(-0.66283488) q[3];
sx q[3];
rz(-2.088306) q[3];
sx q[3];
rz(2.585632) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.8470856) q[2];
sx q[2];
rz(-0.89035788) q[2];
sx q[2];
rz(-0.40714804) q[2];
rz(-1.6242705) q[3];
sx q[3];
rz(-1.1573236) q[3];
sx q[3];
rz(0.24967641) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.80609926) q[0];
sx q[0];
rz(-0.51500106) q[0];
sx q[0];
rz(1.2517713) q[0];
rz(0.66954008) q[1];
sx q[1];
rz(-1.1839097) q[1];
sx q[1];
rz(0.30977419) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1447434) q[0];
sx q[0];
rz(-1.0486756) q[0];
sx q[0];
rz(1.4247308) q[0];
rz(-pi) q[1];
rz(-1.6236213) q[2];
sx q[2];
rz(-2.9644358) q[2];
sx q[2];
rz(-1.0747386) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.5619547) q[1];
sx q[1];
rz(-1.6288174) q[1];
sx q[1];
rz(-1.7768363) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8272607) q[3];
sx q[3];
rz(-2.1520352) q[3];
sx q[3];
rz(-2.6687711) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.4391675) q[2];
sx q[2];
rz(-0.71321407) q[2];
sx q[2];
rz(-1.9343728) q[2];
rz(2.1045945) q[3];
sx q[3];
rz(-1.2456649) q[3];
sx q[3];
rz(2.4859378) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.050215125) q[0];
sx q[0];
rz(-1.3239048) q[0];
sx q[0];
rz(-1.9357095) q[0];
rz(-0.58569113) q[1];
sx q[1];
rz(-1.0810477) q[1];
sx q[1];
rz(-1.6419798) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8382032) q[0];
sx q[0];
rz(-1.2801542) q[0];
sx q[0];
rz(-0.10586664) q[0];
rz(1.502938) q[2];
sx q[2];
rz(-2.3336764) q[2];
sx q[2];
rz(-0.18005904) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.3824532) q[1];
sx q[1];
rz(-1.1083974) q[1];
sx q[1];
rz(-2.9049302) q[1];
rz(-pi) q[2];
rz(-3.1390926) q[3];
sx q[3];
rz(-1.7435939) q[3];
sx q[3];
rz(-2.1713184) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.2575834) q[2];
sx q[2];
rz(-1.7929701) q[2];
sx q[2];
rz(-0.94669) q[2];
rz(-0.36869129) q[3];
sx q[3];
rz(-1.576141) q[3];
sx q[3];
rz(0.45599109) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35836999) q[0];
sx q[0];
rz(-1.9932278) q[0];
sx q[0];
rz(2.7182462) q[0];
rz(-0.070925698) q[1];
sx q[1];
rz(-1.4535041) q[1];
sx q[1];
rz(2.8765875) q[1];
rz(1.0735738) q[2];
sx q[2];
rz(-1.1560658) q[2];
sx q[2];
rz(1.8552468) q[2];
rz(1.0740888) q[3];
sx q[3];
rz(-2.2278193) q[3];
sx q[3];
rz(-0.57886119) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
