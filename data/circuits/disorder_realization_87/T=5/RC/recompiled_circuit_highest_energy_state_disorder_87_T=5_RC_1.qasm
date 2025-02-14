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
rz(2.6142081) q[0];
sx q[0];
rz(-2.5684147) q[0];
sx q[0];
rz(2.5249124) q[0];
rz(-1.0083899) q[1];
sx q[1];
rz(-2.402306) q[1];
sx q[1];
rz(-0.33831212) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2074006) q[0];
sx q[0];
rz(-1.3284042) q[0];
sx q[0];
rz(-1.5336179) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3502928) q[2];
sx q[2];
rz(-1.2630579) q[2];
sx q[2];
rz(-1.6215633) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.1627312) q[1];
sx q[1];
rz(-0.56038364) q[1];
sx q[1];
rz(0.25516971) q[1];
x q[2];
rz(-0.1756535) q[3];
sx q[3];
rz(-0.82161108) q[3];
sx q[3];
rz(2.780811) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.3931291) q[2];
sx q[2];
rz(-1.7270361) q[2];
sx q[2];
rz(2.4836922) q[2];
rz(-2.6864478) q[3];
sx q[3];
rz(-2.9908266) q[3];
sx q[3];
rz(1.3078825) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2748134) q[0];
sx q[0];
rz(-1.4115189) q[0];
sx q[0];
rz(0.28208062) q[0];
rz(2.8981949) q[1];
sx q[1];
rz(-1.8744105) q[1];
sx q[1];
rz(-2.3928221) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3385612) q[0];
sx q[0];
rz(-2.411532) q[0];
sx q[0];
rz(0.29166136) q[0];
rz(-3.1339945) q[2];
sx q[2];
rz(-1.4491334) q[2];
sx q[2];
rz(1.0341782) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.91606027) q[1];
sx q[1];
rz(-0.5114564) q[1];
sx q[1];
rz(-1.7602073) q[1];
x q[2];
rz(-1.5010819) q[3];
sx q[3];
rz(-1.3718714) q[3];
sx q[3];
rz(-2.6903542) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.28027174) q[2];
sx q[2];
rz(-1.6161641) q[2];
sx q[2];
rz(-1.3410428) q[2];
rz(1.9848112) q[3];
sx q[3];
rz(-0.78514922) q[3];
sx q[3];
rz(-2.3780499) q[3];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.84083104) q[0];
sx q[0];
rz(-2.9756727) q[0];
sx q[0];
rz(0.65735835) q[0];
rz(1.1454469) q[1];
sx q[1];
rz(-2.1871388) q[1];
sx q[1];
rz(-2.3232536) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86037246) q[0];
sx q[0];
rz(-1.2272626) q[0];
sx q[0];
rz(-1.8247129) q[0];
x q[1];
rz(-3.0938593) q[2];
sx q[2];
rz(-1.8857419) q[2];
sx q[2];
rz(0.76698179) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.2215449) q[1];
sx q[1];
rz(-0.86600477) q[1];
sx q[1];
rz(2.4815766) q[1];
rz(3.0679296) q[3];
sx q[3];
rz(-1.1182925) q[3];
sx q[3];
rz(2.1114001) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.98693097) q[2];
sx q[2];
rz(-1.2620474) q[2];
sx q[2];
rz(-1.7884802) q[2];
rz(2.1766369) q[3];
sx q[3];
rz(-1.8392287) q[3];
sx q[3];
rz(2.7195948) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44593909) q[0];
sx q[0];
rz(-2.4762479) q[0];
sx q[0];
rz(-1.2009784) q[0];
rz(3.1383842) q[1];
sx q[1];
rz(-1.708958) q[1];
sx q[1];
rz(-2.9046955) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.0081404765) q[0];
sx q[0];
rz(-2.4137375) q[0];
sx q[0];
rz(0.86489622) q[0];
rz(-pi) q[1];
rz(0.17572524) q[2];
sx q[2];
rz(-0.64470664) q[2];
sx q[2];
rz(0.75632655) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.5665) q[1];
sx q[1];
rz(-1.3761576) q[1];
sx q[1];
rz(0.95928488) q[1];
rz(2.3194189) q[3];
sx q[3];
rz(-0.94566761) q[3];
sx q[3];
rz(-1.353598) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.0995471) q[2];
sx q[2];
rz(-1.8501661) q[2];
sx q[2];
rz(-2.715204) q[2];
rz(1.03553) q[3];
sx q[3];
rz(-1.4617045) q[3];
sx q[3];
rz(-2.5199913) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.34430382) q[0];
sx q[0];
rz(-3.099589) q[0];
sx q[0];
rz(-0.34991831) q[0];
rz(1.2087076) q[1];
sx q[1];
rz(-1.2827001) q[1];
sx q[1];
rz(0.0016317687) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2372555) q[0];
sx q[0];
rz(-1.4044824) q[0];
sx q[0];
rz(2.8939899) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5320184) q[2];
sx q[2];
rz(-1.2735521) q[2];
sx q[2];
rz(2.4448038) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.7824911) q[1];
sx q[1];
rz(-1.3939855) q[1];
sx q[1];
rz(0.048080877) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9309421) q[3];
sx q[3];
rz(-1.306476) q[3];
sx q[3];
rz(-1.5624865) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.1899167) q[2];
sx q[2];
rz(-0.92477551) q[2];
sx q[2];
rz(-0.16119371) q[2];
rz(0.4392043) q[3];
sx q[3];
rz(-2.3930211) q[3];
sx q[3];
rz(-0.85328931) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31140232) q[0];
sx q[0];
rz(-1.7364194) q[0];
sx q[0];
rz(2.3468974) q[0];
rz(-1.3658124) q[1];
sx q[1];
rz(-0.8005442) q[1];
sx q[1];
rz(-2.6672003) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.77664033) q[0];
sx q[0];
rz(-2.4863003) q[0];
sx q[0];
rz(2.220015) q[0];
rz(-pi) q[1];
rz(2.5849336) q[2];
sx q[2];
rz(-1.2534598) q[2];
sx q[2];
rz(0.59653004) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.5582592) q[1];
sx q[1];
rz(-2.4129902) q[1];
sx q[1];
rz(0.63459756) q[1];
x q[2];
rz(2.1544564) q[3];
sx q[3];
rz(-1.8963834) q[3];
sx q[3];
rz(1.9501571) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.5279493) q[2];
sx q[2];
rz(-1.8972634) q[2];
sx q[2];
rz(-2.6141686) q[2];
rz(0.6066277) q[3];
sx q[3];
rz(-0.18692034) q[3];
sx q[3];
rz(1.0724148) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8793176) q[0];
sx q[0];
rz(-0.19392218) q[0];
sx q[0];
rz(-2.344017) q[0];
rz(-2.2480615) q[1];
sx q[1];
rz(-0.83502665) q[1];
sx q[1];
rz(-2.8614047) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.93793106) q[0];
sx q[0];
rz(-1.9660006) q[0];
sx q[0];
rz(0.18919887) q[0];
rz(0.75887078) q[2];
sx q[2];
rz(-0.810383) q[2];
sx q[2];
rz(-0.23960613) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.60498991) q[1];
sx q[1];
rz(-2.0532673) q[1];
sx q[1];
rz(-3.1274113) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0810764) q[3];
sx q[3];
rz(-2.4852537) q[3];
sx q[3];
rz(-2.9440232) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.4013275) q[2];
sx q[2];
rz(-1.3046616) q[2];
sx q[2];
rz(1.2124445) q[2];
rz(0.23981833) q[3];
sx q[3];
rz(-2.4739154) q[3];
sx q[3];
rz(-2.5734606) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(2.3577393) q[0];
sx q[0];
rz(-1.4171866) q[0];
sx q[0];
rz(-1.3215815) q[0];
rz(-0.18121885) q[1];
sx q[1];
rz(-1.3316589) q[1];
sx q[1];
rz(0.063057335) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6609333) q[0];
sx q[0];
rz(-1.4278605) q[0];
sx q[0];
rz(2.3602135) q[0];
rz(-1.2793853) q[2];
sx q[2];
rz(-0.29726899) q[2];
sx q[2];
rz(-0.91100541) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.5279505) q[1];
sx q[1];
rz(-2.7220352) q[1];
sx q[1];
rz(0.63906868) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5531055) q[3];
sx q[3];
rz(-2.9790386) q[3];
sx q[3];
rz(1.3471732) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.41398373) q[2];
sx q[2];
rz(-2.0435645) q[2];
sx q[2];
rz(-1.46924) q[2];
rz(1.3937048) q[3];
sx q[3];
rz(-1.4398451) q[3];
sx q[3];
rz(2.9318455) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.53836981) q[0];
sx q[0];
rz(-0.61036888) q[0];
sx q[0];
rz(-0.99697733) q[0];
rz(-2.3485377) q[1];
sx q[1];
rz(-1.7231562) q[1];
sx q[1];
rz(2.0319895) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.133503) q[0];
sx q[0];
rz(-0.64787302) q[0];
sx q[0];
rz(2.6459207) q[0];
rz(-pi) q[1];
x q[1];
rz(0.72203861) q[2];
sx q[2];
rz(-1.2424316) q[2];
sx q[2];
rz(1.9313081) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.80191699) q[1];
sx q[1];
rz(-1.1028506) q[1];
sx q[1];
rz(-0.94593327) q[1];
x q[2];
rz(0.22623541) q[3];
sx q[3];
rz(-0.64877766) q[3];
sx q[3];
rz(1.7895376) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.7290466) q[2];
sx q[2];
rz(-1.6763326) q[2];
sx q[2];
rz(-1.868978) q[2];
rz(2.0160969) q[3];
sx q[3];
rz(-2.9477305) q[3];
sx q[3];
rz(-0.079843609) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.37128714) q[0];
sx q[0];
rz(-1.1530387) q[0];
sx q[0];
rz(3.0433997) q[0];
rz(1.498361) q[1];
sx q[1];
rz(-1.1474835) q[1];
sx q[1];
rz(0.8383382) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.22238092) q[0];
sx q[0];
rz(-1.7475238) q[0];
sx q[0];
rz(2.3803898) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6264084) q[2];
sx q[2];
rz(-1.3190373) q[2];
sx q[2];
rz(-0.97903573) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.34217087) q[1];
sx q[1];
rz(-1.8525181) q[1];
sx q[1];
rz(1.1200302) q[1];
x q[2];
rz(1.7915947) q[3];
sx q[3];
rz(-1.1482557) q[3];
sx q[3];
rz(-0.72964668) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.3457634) q[2];
sx q[2];
rz(-2.2553307) q[2];
sx q[2];
rz(2.8161827) q[2];
rz(2.365153) q[3];
sx q[3];
rz(-2.1263945) q[3];
sx q[3];
rz(-0.6071035) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14071874) q[0];
sx q[0];
rz(-1.9372531) q[0];
sx q[0];
rz(-1.0304864) q[0];
rz(2.8554032) q[1];
sx q[1];
rz(-1.9356526) q[1];
sx q[1];
rz(-2.0331358) q[1];
rz(-2.3807965) q[2];
sx q[2];
rz(-1.8626871) q[2];
sx q[2];
rz(-2.5260851) q[2];
rz(1.9097044) q[3];
sx q[3];
rz(-2.1327426) q[3];
sx q[3];
rz(-0.54973092) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
