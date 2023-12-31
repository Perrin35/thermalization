OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.4296071) q[0];
sx q[0];
rz(-0.40661943) q[0];
sx q[0];
rz(-2.8924195) q[0];
rz(-0.063440032) q[1];
sx q[1];
rz(4.1133147) q[1];
sx q[1];
rz(9.9749554) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7101947) q[0];
sx q[0];
rz(-2.9203127) q[0];
sx q[0];
rz(-1.4472603) q[0];
rz(-pi) q[1];
rz(2.1044188) q[2];
sx q[2];
rz(-1.954477) q[2];
sx q[2];
rz(1.2984315) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.5397415) q[1];
sx q[1];
rz(-1.1907693) q[1];
sx q[1];
rz(-2.3439581) q[1];
rz(-pi) q[2];
rz(-0.95991858) q[3];
sx q[3];
rz(-1.5023408) q[3];
sx q[3];
rz(-1.5271036) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.7258519) q[2];
sx q[2];
rz(-0.44885138) q[2];
sx q[2];
rz(2.501781) q[2];
rz(2.2885627) q[3];
sx q[3];
rz(-2.5363688) q[3];
sx q[3];
rz(-2.7602592) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2663651) q[0];
sx q[0];
rz(-1.1667644) q[0];
sx q[0];
rz(2.8711328) q[0];
rz(-0.71331435) q[1];
sx q[1];
rz(-1.0353054) q[1];
sx q[1];
rz(1.5126022) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0727901) q[0];
sx q[0];
rz(-1.2501984) q[0];
sx q[0];
rz(-2.1945303) q[0];
rz(-1.5741882) q[2];
sx q[2];
rz(-1.8545824) q[2];
sx q[2];
rz(2.0397182) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.1158893) q[1];
sx q[1];
rz(-1.70417) q[1];
sx q[1];
rz(2.2504836) q[1];
rz(-pi) q[2];
rz(0.87725957) q[3];
sx q[3];
rz(-2.7765397) q[3];
sx q[3];
rz(3.0960494) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.9397395) q[2];
sx q[2];
rz(-0.24213232) q[2];
sx q[2];
rz(2.3382323) q[2];
rz(1.057829) q[3];
sx q[3];
rz(-1.6488896) q[3];
sx q[3];
rz(-3.1159744) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
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
rz(1.8818883) q[0];
sx q[0];
rz(-1.0972247) q[0];
sx q[0];
rz(0.29552466) q[0];
rz(2.9064536) q[1];
sx q[1];
rz(-1.7087015) q[1];
sx q[1];
rz(-0.74584109) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5999818) q[0];
sx q[0];
rz(-0.2395425) q[0];
sx q[0];
rz(1.177686) q[0];
rz(-pi) q[1];
rz(-1.4065811) q[2];
sx q[2];
rz(-1.0501554) q[2];
sx q[2];
rz(-2.8525713) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.44887603) q[1];
sx q[1];
rz(-1.8489031) q[1];
sx q[1];
rz(0.49726185) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6607473) q[3];
sx q[3];
rz(-2.6391609) q[3];
sx q[3];
rz(2.0744827) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.088034257) q[2];
sx q[2];
rz(-3.1159248) q[2];
sx q[2];
rz(0.68874613) q[2];
rz(-3.0912494) q[3];
sx q[3];
rz(-0.91209948) q[3];
sx q[3];
rz(1.5215727) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.922309) q[0];
sx q[0];
rz(-2.121121) q[0];
sx q[0];
rz(3.0396089) q[0];
rz(0.12022262) q[1];
sx q[1];
rz(-0.52821237) q[1];
sx q[1];
rz(-0.27329683) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.8768339) q[0];
sx q[0];
rz(-0.33565531) q[0];
sx q[0];
rz(1.8434974) q[0];
rz(0.64812135) q[2];
sx q[2];
rz(-1.1037877) q[2];
sx q[2];
rz(1.8304706) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.4218688) q[1];
sx q[1];
rz(-2.8641652) q[1];
sx q[1];
rz(-1.006702) q[1];
rz(-0.64951879) q[3];
sx q[3];
rz(-1.4276541) q[3];
sx q[3];
rz(-0.79917819) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.79167241) q[2];
sx q[2];
rz(-2.4035954) q[2];
sx q[2];
rz(-2.6376574) q[2];
rz(-0.079581633) q[3];
sx q[3];
rz(-1.1527529) q[3];
sx q[3];
rz(0.3751522) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2824317) q[0];
sx q[0];
rz(-2.0895884) q[0];
sx q[0];
rz(0.082745634) q[0];
rz(0.67963183) q[1];
sx q[1];
rz(-1.4928763) q[1];
sx q[1];
rz(2.1544429) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8866681) q[0];
sx q[0];
rz(-2.3912171) q[0];
sx q[0];
rz(-0.95959856) q[0];
rz(-pi) q[1];
rz(-1.5802025) q[2];
sx q[2];
rz(-2.214589) q[2];
sx q[2];
rz(-2.8807092) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.6189177) q[1];
sx q[1];
rz(-0.72026157) q[1];
sx q[1];
rz(-2.656225) q[1];
rz(-1.6721252) q[3];
sx q[3];
rz(-2.7469027) q[3];
sx q[3];
rz(2.8916388) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.8905028) q[2];
sx q[2];
rz(-0.7901929) q[2];
sx q[2];
rz(-0.21128543) q[2];
rz(-0.42090297) q[3];
sx q[3];
rz(-2.588721) q[3];
sx q[3];
rz(3.0781854) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6565276) q[0];
sx q[0];
rz(-1.8873029) q[0];
sx q[0];
rz(2.3497537) q[0];
rz(-0.99545288) q[1];
sx q[1];
rz(-2.1738926) q[1];
sx q[1];
rz(1.2794367) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6364481) q[0];
sx q[0];
rz(-2.8594058) q[0];
sx q[0];
rz(-1.5178773) q[0];
rz(-pi) q[1];
x q[1];
rz(0.40462599) q[2];
sx q[2];
rz(-0.64794108) q[2];
sx q[2];
rz(0.65587703) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.8705604) q[1];
sx q[1];
rz(-0.522627) q[1];
sx q[1];
rz(1.51103) q[1];
rz(-pi) q[2];
rz(1.7518696) q[3];
sx q[3];
rz(-1.7098134) q[3];
sx q[3];
rz(-2.3982323) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.8906158) q[2];
sx q[2];
rz(-1.0303409) q[2];
sx q[2];
rz(2.3708169) q[2];
rz(-1.4701014) q[3];
sx q[3];
rz(-2.7272868) q[3];
sx q[3];
rz(-0.18338403) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8188748) q[0];
sx q[0];
rz(-0.6495496) q[0];
sx q[0];
rz(0.055667002) q[0];
rz(2.9260013) q[1];
sx q[1];
rz(-0.76342738) q[1];
sx q[1];
rz(-3.1380222) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5647033) q[0];
sx q[0];
rz(-1.2029552) q[0];
sx q[0];
rz(-2.6951615) q[0];
rz(0.75242075) q[2];
sx q[2];
rz(-1.8634897) q[2];
sx q[2];
rz(2.9257286) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(3.0842972) q[1];
sx q[1];
rz(-2.1665386) q[1];
sx q[1];
rz(0.36235313) q[1];
rz(-3.0655541) q[3];
sx q[3];
rz(-1.7685316) q[3];
sx q[3];
rz(-0.34464371) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.5443762) q[2];
sx q[2];
rz(-2.1901013) q[2];
sx q[2];
rz(-0.92010951) q[2];
rz(-0.19872935) q[3];
sx q[3];
rz(-1.2343497) q[3];
sx q[3];
rz(2.7238817) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6253117) q[0];
sx q[0];
rz(-1.5100864) q[0];
sx q[0];
rz(-1.7238808) q[0];
rz(-2.7334546) q[1];
sx q[1];
rz(-2.1022271) q[1];
sx q[1];
rz(-0.67869854) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6060651) q[0];
sx q[0];
rz(-1.8341656) q[0];
sx q[0];
rz(0.33959629) q[0];
rz(2.6762814) q[2];
sx q[2];
rz(-2.2530472) q[2];
sx q[2];
rz(2.0899783) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.855367) q[1];
sx q[1];
rz(-1.6337506) q[1];
sx q[1];
rz(1.0463868) q[1];
rz(-pi) q[2];
rz(-2.802556) q[3];
sx q[3];
rz(-2.6936274) q[3];
sx q[3];
rz(-2.8807004) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.17710182) q[2];
sx q[2];
rz(-1.0937546) q[2];
sx q[2];
rz(2.2237681) q[2];
rz(1.142189) q[3];
sx q[3];
rz(-2.2873788) q[3];
sx q[3];
rz(-2.2043622) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1552102) q[0];
sx q[0];
rz(-2.4676403) q[0];
sx q[0];
rz(0.25979364) q[0];
rz(-2.4329176) q[1];
sx q[1];
rz(-0.27888137) q[1];
sx q[1];
rz(2.6146467) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5567112) q[0];
sx q[0];
rz(-1.3538133) q[0];
sx q[0];
rz(2.7816539) q[0];
rz(-pi) q[1];
rz(0.91295816) q[2];
sx q[2];
rz(-1.7857988) q[2];
sx q[2];
rz(-1.022033) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.867013) q[1];
sx q[1];
rz(-1.8362987) q[1];
sx q[1];
rz(2.8200779) q[1];
rz(-pi) q[2];
rz(-2.6960877) q[3];
sx q[3];
rz(-1.3078491) q[3];
sx q[3];
rz(-1.477369) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.32593411) q[2];
sx q[2];
rz(-0.84647536) q[2];
sx q[2];
rz(0.79088598) q[2];
rz(-2.7549426) q[3];
sx q[3];
rz(-1.0145885) q[3];
sx q[3];
rz(2.8044243) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40795046) q[0];
sx q[0];
rz(-2.9691417) q[0];
sx q[0];
rz(0.98544425) q[0];
rz(-0.5685637) q[1];
sx q[1];
rz(-2.0303346) q[1];
sx q[1];
rz(-0.3607761) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8737313) q[0];
sx q[0];
rz(-0.16927969) q[0];
sx q[0];
rz(1.7196359) q[0];
rz(0.12014328) q[2];
sx q[2];
rz(-1.7842245) q[2];
sx q[2];
rz(-3.0394768) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.68041486) q[1];
sx q[1];
rz(-0.83592452) q[1];
sx q[1];
rz(-1.7635036) q[1];
x q[2];
rz(-2.4264614) q[3];
sx q[3];
rz(-0.22278856) q[3];
sx q[3];
rz(0.95105329) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.0976022) q[2];
sx q[2];
rz(-2.4520935) q[2];
sx q[2];
rz(-2.1981751) q[2];
rz(0.57389456) q[3];
sx q[3];
rz(-2.6608163) q[3];
sx q[3];
rz(1.3963612) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5671134) q[0];
sx q[0];
rz(-1.6945101) q[0];
sx q[0];
rz(2.2816818) q[0];
rz(-1.7815331) q[1];
sx q[1];
rz(-1.0018476) q[1];
sx q[1];
rz(-0.76029653) q[1];
rz(-1.7760989) q[2];
sx q[2];
rz(-1.6663972) q[2];
sx q[2];
rz(-1.321928) q[2];
rz(0.73137024) q[3];
sx q[3];
rz(-0.89805713) q[3];
sx q[3];
rz(-3.0039207) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
