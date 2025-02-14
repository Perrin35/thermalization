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
rz(-2.8415866) q[0];
sx q[0];
rz(-2.0553148) q[0];
sx q[0];
rz(0.80817428) q[0];
rz(1.9975245) q[1];
sx q[1];
rz(-2.3699528) q[1];
sx q[1];
rz(1.5789403) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90955665) q[0];
sx q[0];
rz(-2.252451) q[0];
sx q[0];
rz(-0.85083346) q[0];
x q[1];
rz(1.0536532) q[2];
sx q[2];
rz(-2.8364193) q[2];
sx q[2];
rz(0.15811731) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.3951125) q[1];
sx q[1];
rz(-1.0653138) q[1];
sx q[1];
rz(-0.12684568) q[1];
rz(1.6049625) q[3];
sx q[3];
rz(-1.6164268) q[3];
sx q[3];
rz(-2.6229317) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.943104) q[2];
sx q[2];
rz(-2.9297332) q[2];
sx q[2];
rz(1.0342106) q[2];
rz(-0.69283038) q[3];
sx q[3];
rz(-2.0878017) q[3];
sx q[3];
rz(-2.0809765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8387872) q[0];
sx q[0];
rz(-3.1133339) q[0];
sx q[0];
rz(-2.5651108) q[0];
rz(-0.020596404) q[1];
sx q[1];
rz(-2.7437904) q[1];
sx q[1];
rz(1.0284665) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.27273681) q[0];
sx q[0];
rz(-1.8276365) q[0];
sx q[0];
rz(1.3480074) q[0];
rz(1.9443545) q[2];
sx q[2];
rz(-1.8578863) q[2];
sx q[2];
rz(2.8127828) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.8081144) q[1];
sx q[1];
rz(-1.9030182) q[1];
sx q[1];
rz(1.9971041) q[1];
rz(-pi) q[2];
rz(-2.4941958) q[3];
sx q[3];
rz(-1.8893554) q[3];
sx q[3];
rz(2.7872374) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.2241406) q[2];
sx q[2];
rz(-2.1805306) q[2];
sx q[2];
rz(-1.2646328) q[2];
rz(-2.1680016) q[3];
sx q[3];
rz(-1.4507989) q[3];
sx q[3];
rz(0.26315954) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
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
rz(-0.50315404) q[0];
sx q[0];
rz(-1.3042903) q[0];
sx q[0];
rz(-0.85357443) q[0];
rz(-1.0517993) q[1];
sx q[1];
rz(-2.167326) q[1];
sx q[1];
rz(3.1265756) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.029808345) q[0];
sx q[0];
rz(-1.5795603) q[0];
sx q[0];
rz(-0.86017227) q[0];
x q[1];
rz(-3.0454841) q[2];
sx q[2];
rz(-1.5920361) q[2];
sx q[2];
rz(2.6526895) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.4045118) q[1];
sx q[1];
rz(-0.85857262) q[1];
sx q[1];
rz(-1.668307) q[1];
rz(-2.851753) q[3];
sx q[3];
rz(-1.6320758) q[3];
sx q[3];
rz(1.0672399) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.1198279) q[2];
sx q[2];
rz(-1.652635) q[2];
sx q[2];
rz(-2.8601698) q[2];
rz(-1.4540539) q[3];
sx q[3];
rz(-1.931793) q[3];
sx q[3];
rz(-0.76930261) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5467095) q[0];
sx q[0];
rz(-1.895772) q[0];
sx q[0];
rz(0.17387867) q[0];
rz(1.8100544) q[1];
sx q[1];
rz(-1.4215819) q[1];
sx q[1];
rz(-1.4283659) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15584942) q[0];
sx q[0];
rz(-1.5850889) q[0];
sx q[0];
rz(1.3784268) q[0];
rz(-0.32555737) q[2];
sx q[2];
rz(-2.7635305) q[2];
sx q[2];
rz(-2.4488673) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.62828583) q[1];
sx q[1];
rz(-2.2277955) q[1];
sx q[1];
rz(-1.5221859) q[1];
x q[2];
rz(2.3263116) q[3];
sx q[3];
rz(-0.51183701) q[3];
sx q[3];
rz(-2.3533604) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.4398769) q[2];
sx q[2];
rz(-0.82650799) q[2];
sx q[2];
rz(2.7244549) q[2];
rz(-0.16962984) q[3];
sx q[3];
rz(-3.0483584) q[3];
sx q[3];
rz(-0.72181845) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7222662) q[0];
sx q[0];
rz(-0.37501431) q[0];
sx q[0];
rz(0.30935031) q[0];
rz(-2.3720692) q[1];
sx q[1];
rz(-2.0517495) q[1];
sx q[1];
rz(-1.6228898) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.076974) q[0];
sx q[0];
rz(-0.91320052) q[0];
sx q[0];
rz(0.31975759) q[0];
x q[1];
rz(0.74414545) q[2];
sx q[2];
rz(-0.89140201) q[2];
sx q[2];
rz(1.8216004) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.1479234) q[1];
sx q[1];
rz(-1.5412314) q[1];
sx q[1];
rz(-1.7600585) q[1];
rz(-pi) q[2];
rz(1.7047911) q[3];
sx q[3];
rz(-1.4729285) q[3];
sx q[3];
rz(-2.0201473) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.3994483) q[2];
sx q[2];
rz(-0.98230201) q[2];
sx q[2];
rz(-2.1461416) q[2];
rz(0.71850145) q[3];
sx q[3];
rz(-1.813443) q[3];
sx q[3];
rz(2.3464581) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28498483) q[0];
sx q[0];
rz(-1.6866848) q[0];
sx q[0];
rz(1.6108151) q[0];
rz(1.6948304) q[1];
sx q[1];
rz(-1.4434573) q[1];
sx q[1];
rz(1.4379427) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.7890678) q[0];
sx q[0];
rz(-1.1630327) q[0];
sx q[0];
rz(-3.1192794) q[0];
x q[1];
rz(1.9447359) q[2];
sx q[2];
rz(-1.0122006) q[2];
sx q[2];
rz(2.2064759) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.6561778) q[1];
sx q[1];
rz(-2.6256769) q[1];
sx q[1];
rz(2.9713216) q[1];
rz(0.0030869129) q[3];
sx q[3];
rz(-2.3820138) q[3];
sx q[3];
rz(0.39955968) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.1630254) q[2];
sx q[2];
rz(-1.359553) q[2];
sx q[2];
rz(-1.9095518) q[2];
rz(1.9949404) q[3];
sx q[3];
rz(-1.3403284) q[3];
sx q[3];
rz(2.8310217) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6335886) q[0];
sx q[0];
rz(-2.3242943) q[0];
sx q[0];
rz(2.001413) q[0];
rz(-0.89705244) q[1];
sx q[1];
rz(-1.3373809) q[1];
sx q[1];
rz(0.030489771) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.545118) q[0];
sx q[0];
rz(-1.6870572) q[0];
sx q[0];
rz(-1.7932939) q[0];
x q[1];
rz(2.4989481) q[2];
sx q[2];
rz(-1.3054928) q[2];
sx q[2];
rz(0.20390262) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.5329689) q[1];
sx q[1];
rz(-0.86809671) q[1];
sx q[1];
rz(-2.7547902) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.072149) q[3];
sx q[3];
rz(-0.79550084) q[3];
sx q[3];
rz(1.7444057) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.19411479) q[2];
sx q[2];
rz(-1.0741445) q[2];
sx q[2];
rz(-1.2987632) q[2];
rz(-1.4878368) q[3];
sx q[3];
rz(-1.0349118) q[3];
sx q[3];
rz(-1.3207159) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3461935) q[0];
sx q[0];
rz(-2.8454056) q[0];
sx q[0];
rz(1.3854223) q[0];
rz(-1.9253383) q[1];
sx q[1];
rz(-1.6938035) q[1];
sx q[1];
rz(-0.48929712) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3922244) q[0];
sx q[0];
rz(-0.67182344) q[0];
sx q[0];
rz(3.1098614) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9070315) q[2];
sx q[2];
rz(-1.6052941) q[2];
sx q[2];
rz(0.86755841) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.9605105) q[1];
sx q[1];
rz(-1.1000634) q[1];
sx q[1];
rz(-0.1058284) q[1];
x q[2];
rz(1.9527444) q[3];
sx q[3];
rz(-0.26774613) q[3];
sx q[3];
rz(-1.6588039) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.74799246) q[2];
sx q[2];
rz(-1.1425428) q[2];
sx q[2];
rz(0.46621123) q[2];
rz(1.6097869) q[3];
sx q[3];
rz(-1.5038871) q[3];
sx q[3];
rz(2.8209749) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6521249) q[0];
sx q[0];
rz(-0.89972275) q[0];
sx q[0];
rz(0.049064431) q[0];
rz(-0.24066726) q[1];
sx q[1];
rz(-2.1235695) q[1];
sx q[1];
rz(-3.1191471) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3291859) q[0];
sx q[0];
rz(-2.5617449) q[0];
sx q[0];
rz(-1.5270698) q[0];
x q[1];
rz(2.7540665) q[2];
sx q[2];
rz(-1.8948613) q[2];
sx q[2];
rz(-1.9112196) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.53246236) q[1];
sx q[1];
rz(-2.3098619) q[1];
sx q[1];
rz(-2.4799281) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0231699) q[3];
sx q[3];
rz(-0.9370155) q[3];
sx q[3];
rz(2.3364802) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.13751328) q[2];
sx q[2];
rz(-1.235032) q[2];
sx q[2];
rz(-0.9683041) q[2];
rz(2.3051895) q[3];
sx q[3];
rz(-0.29763779) q[3];
sx q[3];
rz(1.3692726) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4150998) q[0];
sx q[0];
rz(-1.2905755) q[0];
sx q[0];
rz(-2.7628164) q[0];
rz(-3.1355766) q[1];
sx q[1];
rz(-0.52979398) q[1];
sx q[1];
rz(2.5078497) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7597719) q[0];
sx q[0];
rz(-2.2359747) q[0];
sx q[0];
rz(2.0048672) q[0];
rz(-pi) q[1];
rz(1.3258377) q[2];
sx q[2];
rz(-2.5542521) q[2];
sx q[2];
rz(2.7224685) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.72530109) q[1];
sx q[1];
rz(-2.0547199) q[1];
sx q[1];
rz(1.0199706) q[1];
x q[2];
rz(-2.2489088) q[3];
sx q[3];
rz(-1.0859981) q[3];
sx q[3];
rz(0.78692737) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.4358431) q[2];
sx q[2];
rz(-1.6195932) q[2];
sx q[2];
rz(2.5283234) q[2];
rz(2.4663726) q[3];
sx q[3];
rz(-0.37874159) q[3];
sx q[3];
rz(2.3238382) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2262065) q[0];
sx q[0];
rz(-1.362726) q[0];
sx q[0];
rz(1.2363731) q[0];
rz(-0.18648237) q[1];
sx q[1];
rz(-1.7541371) q[1];
sx q[1];
rz(-1.6288155) q[1];
rz(-0.45966799) q[2];
sx q[2];
rz(-1.5528233) q[2];
sx q[2];
rz(-0.96832392) q[2];
rz(-0.9486089) q[3];
sx q[3];
rz(-1.4962248) q[3];
sx q[3];
rz(0.20635508) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
