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
rz(5.7558007) q[0];
sx q[0];
rz(5.7100073) q[0];
sx q[0];
rz(10.041458) q[0];
rz(2.1332027) q[1];
sx q[1];
rz(-0.73928666) q[1];
sx q[1];
rz(-2.8032805) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36114281) q[0];
sx q[0];
rz(-2.8964213) q[0];
sx q[0];
rz(0.14921363) q[0];
rz(-pi) q[1];
rz(1.1461444) q[2];
sx q[2];
rz(-2.3156347) q[2];
sx q[2];
rz(-2.7934157) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.1627312) q[1];
sx q[1];
rz(-0.56038364) q[1];
sx q[1];
rz(-0.25516971) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9659392) q[3];
sx q[3];
rz(-0.82161108) q[3];
sx q[3];
rz(-2.780811) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.7484635) q[2];
sx q[2];
rz(-1.4145565) q[2];
sx q[2];
rz(-2.4836922) q[2];
rz(2.6864478) q[3];
sx q[3];
rz(-2.9908266) q[3];
sx q[3];
rz(1.8337102) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8667792) q[0];
sx q[0];
rz(-1.4115189) q[0];
sx q[0];
rz(-2.859512) q[0];
rz(-2.8981949) q[1];
sx q[1];
rz(-1.2671821) q[1];
sx q[1];
rz(0.74877053) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.80303148) q[0];
sx q[0];
rz(-2.411532) q[0];
sx q[0];
rz(-2.8499313) q[0];
x q[1];
rz(-1.44913) q[2];
sx q[2];
rz(-1.5632544) q[2];
sx q[2];
rz(-2.6040524) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.0091331) q[1];
sx q[1];
rz(-2.0722425) q[1];
sx q[1];
rz(3.0363068) q[1];
rz(-0.33269675) q[3];
sx q[3];
rz(-2.930958) q[3];
sx q[3];
rz(-0.11158195) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.8613209) q[2];
sx q[2];
rz(-1.5254285) q[2];
sx q[2];
rz(1.8005499) q[2];
rz(1.1567814) q[3];
sx q[3];
rz(-2.3564434) q[3];
sx q[3];
rz(-2.3780499) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.84083104) q[0];
sx q[0];
rz(-2.9756727) q[0];
sx q[0];
rz(-0.65735835) q[0];
rz(1.9961458) q[1];
sx q[1];
rz(-2.1871388) q[1];
sx q[1];
rz(2.3232536) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9377146) q[0];
sx q[0];
rz(-2.7174207) q[0];
sx q[0];
rz(-0.61221497) q[0];
rz(-1.4253699) q[2];
sx q[2];
rz(-2.8231695) q[2];
sx q[2];
rz(-0.91998902) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.0249306) q[1];
sx q[1];
rz(-2.0567472) q[1];
sx q[1];
rz(2.3930753) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7210427) q[3];
sx q[3];
rz(-2.683542) q[3];
sx q[3];
rz(2.2786105) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.98693097) q[2];
sx q[2];
rz(-1.8795452) q[2];
sx q[2];
rz(-1.3531125) q[2];
rz(2.1766369) q[3];
sx q[3];
rz(-1.302364) q[3];
sx q[3];
rz(-2.7195948) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6956536) q[0];
sx q[0];
rz(-2.4762479) q[0];
sx q[0];
rz(1.9406142) q[0];
rz(-3.1383842) q[1];
sx q[1];
rz(-1.4326347) q[1];
sx q[1];
rz(-2.9046955) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1334522) q[0];
sx q[0];
rz(-2.4137375) q[0];
sx q[0];
rz(0.86489622) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4400993) q[2];
sx q[2];
rz(-0.937619) q[2];
sx q[2];
rz(-2.6038632) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.57509267) q[1];
sx q[1];
rz(-1.3761576) q[1];
sx q[1];
rz(2.1823078) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3194189) q[3];
sx q[3];
rz(-0.94566761) q[3];
sx q[3];
rz(1.353598) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.0995471) q[2];
sx q[2];
rz(-1.2914265) q[2];
sx q[2];
rz(2.715204) q[2];
rz(2.1060627) q[3];
sx q[3];
rz(-1.4617045) q[3];
sx q[3];
rz(2.5199913) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7972888) q[0];
sx q[0];
rz(-0.042003691) q[0];
sx q[0];
rz(2.7916743) q[0];
rz(-1.2087076) q[1];
sx q[1];
rz(-1.2827001) q[1];
sx q[1];
rz(-0.0016317687) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2372555) q[0];
sx q[0];
rz(-1.7371103) q[0];
sx q[0];
rz(-0.24760274) q[0];
rz(-pi) q[1];
rz(0.12589215) q[2];
sx q[2];
rz(-0.29968867) q[2];
sx q[2];
rz(0.82848779) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.220158) q[1];
sx q[1];
rz(-1.618127) q[1];
sx q[1];
rz(1.3937852) q[1];
rz(-2.9309421) q[3];
sx q[3];
rz(-1.8351166) q[3];
sx q[3];
rz(-1.5791062) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.95167595) q[2];
sx q[2];
rz(-2.2168171) q[2];
sx q[2];
rz(2.9803989) q[2];
rz(-0.4392043) q[3];
sx q[3];
rz(-0.74857155) q[3];
sx q[3];
rz(2.2883033) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31140232) q[0];
sx q[0];
rz(-1.7364194) q[0];
sx q[0];
rz(-2.3468974) q[0];
rz(1.3658124) q[1];
sx q[1];
rz(-2.3410485) q[1];
sx q[1];
rz(0.47439233) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.77664033) q[0];
sx q[0];
rz(-2.4863003) q[0];
sx q[0];
rz(-2.220015) q[0];
rz(-0.55616711) q[2];
sx q[2];
rz(-2.5092297) q[2];
sx q[2];
rz(1.4390505) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.5833334) q[1];
sx q[1];
rz(-2.4129902) q[1];
sx q[1];
rz(-0.63459756) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1544564) q[3];
sx q[3];
rz(-1.2452092) q[3];
sx q[3];
rz(-1.9501571) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.5279493) q[2];
sx q[2];
rz(-1.2443292) q[2];
sx q[2];
rz(-0.5274241) q[2];
rz(-0.6066277) q[3];
sx q[3];
rz(-0.18692034) q[3];
sx q[3];
rz(2.0691779) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2622751) q[0];
sx q[0];
rz(-0.19392218) q[0];
sx q[0];
rz(-2.344017) q[0];
rz(0.89353117) q[1];
sx q[1];
rz(-2.306566) q[1];
sx q[1];
rz(2.8614047) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2036616) q[0];
sx q[0];
rz(-1.9660006) q[0];
sx q[0];
rz(2.9523938) q[0];
rz(-pi) q[1];
x q[1];
rz(0.75887078) q[2];
sx q[2];
rz(-2.3312097) q[2];
sx q[2];
rz(0.23960613) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.97238648) q[1];
sx q[1];
rz(-1.5833588) q[1];
sx q[1];
rz(2.0533086) q[1];
rz(-pi) q[2];
rz(2.0810764) q[3];
sx q[3];
rz(-2.4852537) q[3];
sx q[3];
rz(2.9440232) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.7402652) q[2];
sx q[2];
rz(-1.836931) q[2];
sx q[2];
rz(-1.9291482) q[2];
rz(2.9017743) q[3];
sx q[3];
rz(-2.4739154) q[3];
sx q[3];
rz(2.5734606) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3577393) q[0];
sx q[0];
rz(-1.724406) q[0];
sx q[0];
rz(1.3215815) q[0];
rz(-0.18121885) q[1];
sx q[1];
rz(-1.8099338) q[1];
sx q[1];
rz(3.0785353) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6609333) q[0];
sx q[0];
rz(-1.4278605) q[0];
sx q[0];
rz(2.3602135) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8562154) q[2];
sx q[2];
rz(-1.4865424) q[2];
sx q[2];
rz(-0.38049305) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.588187) q[1];
sx q[1];
rz(-1.3253731) q[1];
sx q[1];
rz(-2.7977944) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.1386915) q[3];
sx q[3];
rz(-1.7333247) q[3];
sx q[3];
rz(-1.8123466) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.41398373) q[2];
sx q[2];
rz(-1.0980282) q[2];
sx q[2];
rz(1.6723527) q[2];
rz(1.3937048) q[3];
sx q[3];
rz(-1.4398451) q[3];
sx q[3];
rz(-0.20974717) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.53836981) q[0];
sx q[0];
rz(-2.5312238) q[0];
sx q[0];
rz(0.99697733) q[0];
rz(-0.793055) q[1];
sx q[1];
rz(-1.4184364) q[1];
sx q[1];
rz(2.0319895) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.0080896688) q[0];
sx q[0];
rz(-0.64787302) q[0];
sx q[0];
rz(0.4956719) q[0];
x q[1];
rz(1.9969606) q[2];
sx q[2];
rz(-0.89489102) q[2];
sx q[2];
rz(3.0578095) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.0834137) q[1];
sx q[1];
rz(-2.1201057) q[1];
sx q[1];
rz(2.5843402) q[1];
rz(-pi) q[2];
x q[2];
rz(0.63639464) q[3];
sx q[3];
rz(-1.4348467) q[3];
sx q[3];
rz(-3.1042447) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.7290466) q[2];
sx q[2];
rz(-1.6763326) q[2];
sx q[2];
rz(-1.2726146) q[2];
rz(2.0160969) q[3];
sx q[3];
rz(-2.9477305) q[3];
sx q[3];
rz(-0.079843609) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
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
rz(2.7703055) q[0];
sx q[0];
rz(-1.9885539) q[0];
sx q[0];
rz(-0.09819296) q[0];
rz(-1.498361) q[1];
sx q[1];
rz(-1.1474835) q[1];
sx q[1];
rz(-0.8383382) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.22238092) q[0];
sx q[0];
rz(-1.3940689) q[0];
sx q[0];
rz(-0.76120283) q[0];
x q[1];
rz(1.5151843) q[2];
sx q[2];
rz(-1.3190373) q[2];
sx q[2];
rz(2.1625569) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.7994218) q[1];
sx q[1];
rz(-1.8525181) q[1];
sx q[1];
rz(-1.1200302) q[1];
rz(-pi) q[2];
rz(-1.3499979) q[3];
sx q[3];
rz(-1.9933369) q[3];
sx q[3];
rz(-2.411946) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.7958293) q[2];
sx q[2];
rz(-0.886262) q[2];
sx q[2];
rz(-2.8161827) q[2];
rz(0.77643967) q[3];
sx q[3];
rz(-1.0151981) q[3];
sx q[3];
rz(2.5344892) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0008739) q[0];
sx q[0];
rz(-1.9372531) q[0];
sx q[0];
rz(-1.0304864) q[0];
rz(-2.8554032) q[1];
sx q[1];
rz(-1.20594) q[1];
sx q[1];
rz(1.1084569) q[1];
rz(-2.7306225) q[2];
sx q[2];
rz(-2.3373418) q[2];
sx q[2];
rz(1.892754) q[2];
rz(0.48578942) q[3];
sx q[3];
rz(-0.64668568) q[3];
sx q[3];
rz(0.034737094) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
