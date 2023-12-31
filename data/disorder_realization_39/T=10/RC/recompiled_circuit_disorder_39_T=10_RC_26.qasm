OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.0467779) q[0];
sx q[0];
rz(-1.0682286) q[0];
sx q[0];
rz(-0.46407035) q[0];
rz(1.9595454) q[1];
sx q[1];
rz(6.2160677) q[1];
sx q[1];
rz(10.709229) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0425134) q[0];
sx q[0];
rz(-0.49603841) q[0];
sx q[0];
rz(-0.92143671) q[0];
x q[1];
rz(2.651865) q[2];
sx q[2];
rz(-1.5999319) q[2];
sx q[2];
rz(-2.4907128) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.7317036) q[1];
sx q[1];
rz(-2.0239502) q[1];
sx q[1];
rz(-1.975592) q[1];
rz(0.15070446) q[3];
sx q[3];
rz(-1.5871443) q[3];
sx q[3];
rz(2.5022262) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.7444732) q[2];
sx q[2];
rz(-0.93966547) q[2];
sx q[2];
rz(2.822067) q[2];
rz(0.5784936) q[3];
sx q[3];
rz(-0.47839034) q[3];
sx q[3];
rz(2.4676676) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7330866) q[0];
sx q[0];
rz(-1.5698313) q[0];
sx q[0];
rz(2.615036) q[0];
rz(-0.56354848) q[1];
sx q[1];
rz(-1.6105517) q[1];
sx q[1];
rz(0.79663509) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54290402) q[0];
sx q[0];
rz(-1.2447378) q[0];
sx q[0];
rz(-2.23404) q[0];
x q[1];
rz(-0.30226207) q[2];
sx q[2];
rz(-0.59603359) q[2];
sx q[2];
rz(2.9994534) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.7318774) q[1];
sx q[1];
rz(-2.5341923) q[1];
sx q[1];
rz(3.1190447) q[1];
rz(-pi) q[2];
rz(-2.1809686) q[3];
sx q[3];
rz(-1.989813) q[3];
sx q[3];
rz(-0.39715365) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.9600296) q[2];
sx q[2];
rz(-1.7627565) q[2];
sx q[2];
rz(2.3550418) q[2];
rz(-2.6484047) q[3];
sx q[3];
rz(-1.9599873) q[3];
sx q[3];
rz(0.44737679) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2717473) q[0];
sx q[0];
rz(-2.2646876) q[0];
sx q[0];
rz(-1.440381) q[0];
rz(0.72021833) q[1];
sx q[1];
rz(-2.5455988) q[1];
sx q[1];
rz(2.4386141) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5876578) q[0];
sx q[0];
rz(-2.8371187) q[0];
sx q[0];
rz(-1.9703883) q[0];
rz(-pi) q[1];
rz(-2.5419652) q[2];
sx q[2];
rz(-2.5587974) q[2];
sx q[2];
rz(1.8011013) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.86537251) q[1];
sx q[1];
rz(-0.11905383) q[1];
sx q[1];
rz(-0.40465506) q[1];
x q[2];
rz(-2.1941575) q[3];
sx q[3];
rz(-1.4294335) q[3];
sx q[3];
rz(-2.3879104) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.8748223) q[2];
sx q[2];
rz(-1.4845994) q[2];
sx q[2];
rz(2.2300569) q[2];
rz(-0.91397816) q[3];
sx q[3];
rz(-1.711288) q[3];
sx q[3];
rz(-2.3142557) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5754958) q[0];
sx q[0];
rz(-0.63593447) q[0];
sx q[0];
rz(2.3655868) q[0];
rz(-1.874118) q[1];
sx q[1];
rz(-1.1088561) q[1];
sx q[1];
rz(-0.56328303) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.253809) q[0];
sx q[0];
rz(-0.94622181) q[0];
sx q[0];
rz(0.33872351) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.69022501) q[2];
sx q[2];
rz(-1.8340655) q[2];
sx q[2];
rz(0.0036247591) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.5967484) q[1];
sx q[1];
rz(-1.8425643) q[1];
sx q[1];
rz(1.5775058) q[1];
x q[2];
rz(-2.8068845) q[3];
sx q[3];
rz(-0.74581205) q[3];
sx q[3];
rz(1.0583744) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.6528066) q[2];
sx q[2];
rz(-2.2061009) q[2];
sx q[2];
rz(0.164786) q[2];
rz(2.9131043) q[3];
sx q[3];
rz(-2.861664) q[3];
sx q[3];
rz(2.1951108) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8673458) q[0];
sx q[0];
rz(-0.86425257) q[0];
sx q[0];
rz(-2.2891323) q[0];
rz(0.35119855) q[1];
sx q[1];
rz(-2.5636667) q[1];
sx q[1];
rz(1.9794827) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4441372) q[0];
sx q[0];
rz(-1.1778957) q[0];
sx q[0];
rz(0.59209728) q[0];
rz(-pi) q[1];
rz(-2.0931307) q[2];
sx q[2];
rz(-1.221721) q[2];
sx q[2];
rz(-2.366684) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.4967242) q[1];
sx q[1];
rz(-0.86360303) q[1];
sx q[1];
rz(1.0134539) q[1];
rz(-pi) q[2];
rz(1.5557628) q[3];
sx q[3];
rz(-0.84926499) q[3];
sx q[3];
rz(-2.815747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.44624415) q[2];
sx q[2];
rz(-2.7480795) q[2];
sx q[2];
rz(-1.2472786) q[2];
rz(-0.013109664) q[3];
sx q[3];
rz(-1.40991) q[3];
sx q[3];
rz(2.9124027) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2728249) q[0];
sx q[0];
rz(-1.2478963) q[0];
sx q[0];
rz(0.75063467) q[0];
rz(-1.1095095) q[1];
sx q[1];
rz(-2.7710319) q[1];
sx q[1];
rz(1.3060588) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40697843) q[0];
sx q[0];
rz(-2.4771871) q[0];
sx q[0];
rz(2.2218496) q[0];
rz(-pi) q[1];
rz(-0.46361228) q[2];
sx q[2];
rz(-1.7525275) q[2];
sx q[2];
rz(2.6169427) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.8667824) q[1];
sx q[1];
rz(-0.85046235) q[1];
sx q[1];
rz(1.5809098) q[1];
rz(1.4218876) q[3];
sx q[3];
rz(-2.1544666) q[3];
sx q[3];
rz(-0.35051171) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.3874454) q[2];
sx q[2];
rz(-1.0070609) q[2];
sx q[2];
rz(0.60069096) q[2];
rz(-2.1389652) q[3];
sx q[3];
rz(-2.6032344) q[3];
sx q[3];
rz(0.35282648) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3951185) q[0];
sx q[0];
rz(-1.4557319) q[0];
sx q[0];
rz(2.0948998) q[0];
rz(-1.5294317) q[1];
sx q[1];
rz(-1.4271095) q[1];
sx q[1];
rz(0.41710645) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1244303) q[0];
sx q[0];
rz(-2.9201047) q[0];
sx q[0];
rz(-1.4593967) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8344526) q[2];
sx q[2];
rz(-1.8516314) q[2];
sx q[2];
rz(1.3020696) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.6542146) q[1];
sx q[1];
rz(-2.0301135) q[1];
sx q[1];
rz(-0.76141255) q[1];
x q[2];
rz(0.65317513) q[3];
sx q[3];
rz(-1.3305059) q[3];
sx q[3];
rz(-1.0550635) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(3.1356915) q[2];
sx q[2];
rz(-0.4726755) q[2];
sx q[2];
rz(-1.3674412) q[2];
rz(-2.588429) q[3];
sx q[3];
rz(-1.4033214) q[3];
sx q[3];
rz(-0.52545351) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-2.816514) q[0];
sx q[0];
rz(-1.8186318) q[0];
sx q[0];
rz(1.1897855) q[0];
rz(-1.4272383) q[1];
sx q[1];
rz(-2.863527) q[1];
sx q[1];
rz(-0.11238012) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1988797) q[0];
sx q[0];
rz(-0.50857022) q[0];
sx q[0];
rz(1.2099427) q[0];
rz(-pi) q[1];
rz(1.6582279) q[2];
sx q[2];
rz(-1.2382675) q[2];
sx q[2];
rz(2.33193) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.1689414) q[1];
sx q[1];
rz(-0.87247889) q[1];
sx q[1];
rz(-2.6886743) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3120679) q[3];
sx q[3];
rz(-1.6885763) q[3];
sx q[3];
rz(-3.0974914) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.0671976) q[2];
sx q[2];
rz(-1.9125568) q[2];
sx q[2];
rz(-1.7929662) q[2];
rz(1.9366692) q[3];
sx q[3];
rz(-1.4068853) q[3];
sx q[3];
rz(-2.9437734) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7444721) q[0];
sx q[0];
rz(-1.4646126) q[0];
sx q[0];
rz(-0.034974139) q[0];
rz(-2.2947625) q[1];
sx q[1];
rz(-1.7085107) q[1];
sx q[1];
rz(-2.2299178) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.002279) q[0];
sx q[0];
rz(-0.3542491) q[0];
sx q[0];
rz(0.87689633) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6976835) q[2];
sx q[2];
rz(-1.03089) q[2];
sx q[2];
rz(-1.1586231) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.4166491) q[1];
sx q[1];
rz(-1.1133725) q[1];
sx q[1];
rz(-1.7979421) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0802286) q[3];
sx q[3];
rz(-1.0883696) q[3];
sx q[3];
rz(-0.8182984) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.4336865) q[2];
sx q[2];
rz(-2.6277379) q[2];
sx q[2];
rz(-2.8923477) q[2];
rz(0.76672673) q[3];
sx q[3];
rz(-1.8079575) q[3];
sx q[3];
rz(-2.7419817) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7837759) q[0];
sx q[0];
rz(-1.0483085) q[0];
sx q[0];
rz(0.37208474) q[0];
rz(2.5601939) q[1];
sx q[1];
rz(-2.364368) q[1];
sx q[1];
rz(-1.7262329) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6045195) q[0];
sx q[0];
rz(-1.8329289) q[0];
sx q[0];
rz(1.7781236) q[0];
rz(-pi) q[1];
x q[1];
rz(2.084311) q[2];
sx q[2];
rz(-2.1272749) q[2];
sx q[2];
rz(1.8884115) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.33503867) q[1];
sx q[1];
rz(-1.5015258) q[1];
sx q[1];
rz(-0.075242234) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0269208) q[3];
sx q[3];
rz(-1.9034991) q[3];
sx q[3];
rz(-2.016891) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.8156478) q[2];
sx q[2];
rz(-0.14020136) q[2];
sx q[2];
rz(0.20467219) q[2];
rz(-1.7278016) q[3];
sx q[3];
rz(-1.9982343) q[3];
sx q[3];
rz(-1.0958825) q[3];
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
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.50080147) q[0];
sx q[0];
rz(-1.4688107) q[0];
sx q[0];
rz(-2.2039913) q[0];
rz(-1.5851371) q[1];
sx q[1];
rz(-1.382348) q[1];
sx q[1];
rz(1.4437645) q[1];
rz(1.3265058) q[2];
sx q[2];
rz(-2.6752224) q[2];
sx q[2];
rz(0.35717076) q[2];
rz(0.99863573) q[3];
sx q[3];
rz(-1.640366) q[3];
sx q[3];
rz(0.58983005) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
