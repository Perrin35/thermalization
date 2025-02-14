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
rz(-0.93936062) q[0];
sx q[0];
rz(4.2606104) q[0];
sx q[0];
rz(9.4584447) q[0];
rz(1.6839924) q[1];
sx q[1];
rz(-1.3352609) q[1];
sx q[1];
rz(-1.508498) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2062708) q[0];
sx q[0];
rz(-1.6007206) q[0];
sx q[0];
rz(-1.4564092) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2766507) q[2];
sx q[2];
rz(-1.642881) q[2];
sx q[2];
rz(0.26157899) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.37225393) q[1];
sx q[1];
rz(-2.4318691) q[1];
sx q[1];
rz(2.918887) q[1];
x q[2];
rz(-0.63920984) q[3];
sx q[3];
rz(-2.2230801) q[3];
sx q[3];
rz(-0.88310421) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.17026751) q[2];
sx q[2];
rz(-2.2496932) q[2];
sx q[2];
rz(0.83667052) q[2];
rz(0.80638805) q[3];
sx q[3];
rz(-0.92525768) q[3];
sx q[3];
rz(-2.9562318) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0760913) q[0];
sx q[0];
rz(-1.0668904) q[0];
sx q[0];
rz(-1.9245032) q[0];
rz(-2.941046) q[1];
sx q[1];
rz(-2.3417818) q[1];
sx q[1];
rz(0.54380551) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1233305) q[0];
sx q[0];
rz(-1.2549434) q[0];
sx q[0];
rz(-1.3297775) q[0];
rz(2.5466304) q[2];
sx q[2];
rz(-1.5820855) q[2];
sx q[2];
rz(-1.4608698) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.041808279) q[1];
sx q[1];
rz(-1.1893032) q[1];
sx q[1];
rz(0.31111335) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.614566) q[3];
sx q[3];
rz(-2.4659116) q[3];
sx q[3];
rz(2.0213493) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.66889888) q[2];
sx q[2];
rz(-2.3987179) q[2];
sx q[2];
rz(1.5354068) q[2];
rz(-1.499048) q[3];
sx q[3];
rz(-0.58321548) q[3];
sx q[3];
rz(-1.8892052) q[3];
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
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1010308) q[0];
sx q[0];
rz(-2.5767548) q[0];
sx q[0];
rz(-0.0095796883) q[0];
rz(-1.0293695) q[1];
sx q[1];
rz(-1.5734438) q[1];
sx q[1];
rz(1.3363438) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.97735092) q[0];
sx q[0];
rz(-1.0625496) q[0];
sx q[0];
rz(-2.4524816) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8200978) q[2];
sx q[2];
rz(-1.090637) q[2];
sx q[2];
rz(-2.8804422) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(3.134446) q[1];
sx q[1];
rz(-1.2255417) q[1];
sx q[1];
rz(-0.18642872) q[1];
rz(-pi) q[2];
rz(-1.9668749) q[3];
sx q[3];
rz(-1.4142591) q[3];
sx q[3];
rz(0.96352947) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.10793992) q[2];
sx q[2];
rz(-0.5113217) q[2];
sx q[2];
rz(1.8281724) q[2];
rz(3.1016453) q[3];
sx q[3];
rz(-2.2291456) q[3];
sx q[3];
rz(2.0136755) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(-0.1333328) q[0];
sx q[0];
rz(-1.2902211) q[0];
sx q[0];
rz(2.7449352) q[0];
rz(2.2465514) q[1];
sx q[1];
rz(-1.7170693) q[1];
sx q[1];
rz(2.3779714) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7437282) q[0];
sx q[0];
rz(-1.3258805) q[0];
sx q[0];
rz(-3.0586277) q[0];
rz(-pi) q[1];
rz(3.0210546) q[2];
sx q[2];
rz(-1.5539031) q[2];
sx q[2];
rz(2.8617045) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.7553923) q[1];
sx q[1];
rz(-1.5441276) q[1];
sx q[1];
rz(1.8083354) q[1];
rz(-pi) q[2];
x q[2];
rz(0.29080363) q[3];
sx q[3];
rz(-1.606308) q[3];
sx q[3];
rz(2.5574656) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.414173) q[2];
sx q[2];
rz(-1.6215934) q[2];
sx q[2];
rz(-2.011389) q[2];
rz(-1.1388418) q[3];
sx q[3];
rz(-1.1915221) q[3];
sx q[3];
rz(1.9738919) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.30911699) q[0];
sx q[0];
rz(-2.583355) q[0];
sx q[0];
rz(-2.6600237) q[0];
rz(-0.32749495) q[1];
sx q[1];
rz(-1.7227252) q[1];
sx q[1];
rz(-0.96424261) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6511232) q[0];
sx q[0];
rz(-1.3523605) q[0];
sx q[0];
rz(-1.1743702) q[0];
rz(-pi) q[1];
rz(-2.9498001) q[2];
sx q[2];
rz(-0.99770498) q[2];
sx q[2];
rz(-1.4258143) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.89470071) q[1];
sx q[1];
rz(-2.8671226) q[1];
sx q[1];
rz(2.0361221) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7751415) q[3];
sx q[3];
rz(-1.1875523) q[3];
sx q[3];
rz(0.97571532) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.1714736) q[2];
sx q[2];
rz(-1.118266) q[2];
sx q[2];
rz(1.4946651) q[2];
rz(-2.1286879) q[3];
sx q[3];
rz(-2.2154112) q[3];
sx q[3];
rz(-0.51500285) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9377015) q[0];
sx q[0];
rz(-2.1370115) q[0];
sx q[0];
rz(-1.0585693) q[0];
rz(-0.92264289) q[1];
sx q[1];
rz(-2.1244996) q[1];
sx q[1];
rz(-0.10291084) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.12891586) q[0];
sx q[0];
rz(-0.8230477) q[0];
sx q[0];
rz(0.91632592) q[0];
rz(-1.2936959) q[2];
sx q[2];
rz(-2.1418946) q[2];
sx q[2];
rz(-1.7383611) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.7988777) q[1];
sx q[1];
rz(-0.56523501) q[1];
sx q[1];
rz(-1.2223627) q[1];
rz(-pi) q[2];
x q[2];
rz(0.0062557022) q[3];
sx q[3];
rz(-0.85438529) q[3];
sx q[3];
rz(2.7992918) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.70204488) q[2];
sx q[2];
rz(-2.6921258) q[2];
sx q[2];
rz(-2.8257545) q[2];
rz(2.1465007) q[3];
sx q[3];
rz(-1.8683878) q[3];
sx q[3];
rz(0.72527138) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.134326) q[0];
sx q[0];
rz(-0.26679978) q[0];
sx q[0];
rz(-2.5079978) q[0];
rz(0.38465706) q[1];
sx q[1];
rz(-1.9622012) q[1];
sx q[1];
rz(2.7124009) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10776183) q[0];
sx q[0];
rz(-1.1035447) q[0];
sx q[0];
rz(2.856887) q[0];
rz(-0.68558461) q[2];
sx q[2];
rz(-1.597817) q[2];
sx q[2];
rz(2.9179058) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.8365578) q[1];
sx q[1];
rz(-1.9227131) q[1];
sx q[1];
rz(2.5266527) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6771992) q[3];
sx q[3];
rz(-1.0611609) q[3];
sx q[3];
rz(-1.2023752) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(3.0523494) q[2];
sx q[2];
rz(-2.3640335) q[2];
sx q[2];
rz(0.26697978) q[2];
rz(-3.0698245) q[3];
sx q[3];
rz(-1.0177344) q[3];
sx q[3];
rz(-1.4192386) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4778022) q[0];
sx q[0];
rz(-0.1013805) q[0];
sx q[0];
rz(-1.4005533) q[0];
rz(1.1446704) q[1];
sx q[1];
rz(-1.0619699) q[1];
sx q[1];
rz(2.5544419) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5349197) q[0];
sx q[0];
rz(-0.79100709) q[0];
sx q[0];
rz(-0.61137166) q[0];
rz(-pi) q[1];
rz(-1.4540423) q[2];
sx q[2];
rz(-1.1551305) q[2];
sx q[2];
rz(-1.0041725) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.352601) q[1];
sx q[1];
rz(-1.3459715) q[1];
sx q[1];
rz(-2.4487752) q[1];
rz(-2.1903992) q[3];
sx q[3];
rz(-1.3096598) q[3];
sx q[3];
rz(-1.2984683) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.8553541) q[2];
sx q[2];
rz(-1.9619532) q[2];
sx q[2];
rz(2.722495) q[2];
rz(1.6176443) q[3];
sx q[3];
rz(-2.2231299) q[3];
sx q[3];
rz(-1.6379448) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5986901) q[0];
sx q[0];
rz(-1.221523) q[0];
sx q[0];
rz(0.30938095) q[0];
rz(-1.7912553) q[1];
sx q[1];
rz(-0.58809909) q[1];
sx q[1];
rz(2.5877171) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0623035) q[0];
sx q[0];
rz(-1.252018) q[0];
sx q[0];
rz(-0.43894025) q[0];
x q[1];
rz(-1.1872841) q[2];
sx q[2];
rz(-1.4330136) q[2];
sx q[2];
rz(-0.99064186) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.9123037) q[1];
sx q[1];
rz(-0.69432753) q[1];
sx q[1];
rz(2.0635384) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7782544) q[3];
sx q[3];
rz(-0.73550341) q[3];
sx q[3];
rz(2.9309762) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.89173335) q[2];
sx q[2];
rz(-1.5313818) q[2];
sx q[2];
rz(1.3472793) q[2];
rz(-0.76752082) q[3];
sx q[3];
rz(-1.1859505) q[3];
sx q[3];
rz(2.2931113) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7193741) q[0];
sx q[0];
rz(-0.30000559) q[0];
sx q[0];
rz(-1.9268217) q[0];
rz(-1.0019373) q[1];
sx q[1];
rz(-1.5129713) q[1];
sx q[1];
rz(-0.92438662) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0262605) q[0];
sx q[0];
rz(-2.1450274) q[0];
sx q[0];
rz(2.6747245) q[0];
rz(2.1154249) q[2];
sx q[2];
rz(-1.5942425) q[2];
sx q[2];
rz(0.12870994) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.9634642) q[1];
sx q[1];
rz(-0.22451065) q[1];
sx q[1];
rz(3.0369395) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2460292) q[3];
sx q[3];
rz(-2.319284) q[3];
sx q[3];
rz(-3.0117717) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.34710106) q[2];
sx q[2];
rz(-2.2091986) q[2];
sx q[2];
rz(-0.36716983) q[2];
rz(1.1651039) q[3];
sx q[3];
rz(-2.1737183) q[3];
sx q[3];
rz(0.86618209) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9544871) q[0];
sx q[0];
rz(-0.73910537) q[0];
sx q[0];
rz(0.11866971) q[0];
rz(0.26609303) q[1];
sx q[1];
rz(-0.91231822) q[1];
sx q[1];
rz(-1.4516861) q[1];
rz(-2.3824206) q[2];
sx q[2];
rz(-2.7552032) q[2];
sx q[2];
rz(-0.0061782171) q[2];
rz(-2.9333276) q[3];
sx q[3];
rz(-0.97934813) q[3];
sx q[3];
rz(-1.4772268) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
