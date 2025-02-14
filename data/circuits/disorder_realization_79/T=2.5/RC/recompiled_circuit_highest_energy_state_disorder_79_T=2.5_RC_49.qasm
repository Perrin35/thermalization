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
rz(0.6495629) q[0];
sx q[0];
rz(-0.49512884) q[0];
sx q[0];
rz(-2.8861217) q[0];
rz(-0.35222346) q[1];
sx q[1];
rz(4.5402266) q[1];
sx q[1];
rz(11.160688) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.79069239) q[0];
sx q[0];
rz(-3.1381099) q[0];
sx q[0];
rz(-3.0538959) q[0];
rz(-2.6256526) q[2];
sx q[2];
rz(-1.6637515) q[2];
sx q[2];
rz(-0.0001212349) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.97326255) q[1];
sx q[1];
rz(-1.574206) q[1];
sx q[1];
rz(1.1876606) q[1];
rz(-0.56804232) q[3];
sx q[3];
rz(-1.5431817) q[3];
sx q[3];
rz(1.9909007) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.5070255) q[2];
sx q[2];
rz(-0.031012379) q[2];
sx q[2];
rz(0.8325038) q[2];
rz(-0.80820525) q[3];
sx q[3];
rz(-0.015232239) q[3];
sx q[3];
rz(0.25850779) q[3];
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
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7611564) q[0];
sx q[0];
rz(-0.34270898) q[0];
sx q[0];
rz(-0.1880745) q[0];
rz(-3.0694118) q[1];
sx q[1];
rz(-1.0344104) q[1];
sx q[1];
rz(-1.5123051) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24969026) q[0];
sx q[0];
rz(-1.238126) q[0];
sx q[0];
rz(1.4494677) q[0];
x q[1];
rz(-0.030628344) q[2];
sx q[2];
rz(-1.5248486) q[2];
sx q[2];
rz(2.6776938) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.9058414) q[1];
sx q[1];
rz(-3.0260147) q[1];
sx q[1];
rz(0.53821941) q[1];
rz(-pi) q[2];
rz(2.333669) q[3];
sx q[3];
rz(-1.9821385) q[3];
sx q[3];
rz(-1.5802671) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.9709116) q[2];
sx q[2];
rz(-2.4965681) q[2];
sx q[2];
rz(-1.2930653) q[2];
rz(-0.34717789) q[3];
sx q[3];
rz(-2.8721589) q[3];
sx q[3];
rz(-0.080304317) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3229356) q[0];
sx q[0];
rz(-1.8523536) q[0];
sx q[0];
rz(2.6710508) q[0];
rz(-1.200354) q[1];
sx q[1];
rz(-0.73089868) q[1];
sx q[1];
rz(1.8847195) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3124533) q[0];
sx q[0];
rz(-1.0395673) q[0];
sx q[0];
rz(-0.4108528) q[0];
rz(-pi) q[1];
rz(2.1684709) q[2];
sx q[2];
rz(-0.16301708) q[2];
sx q[2];
rz(-0.68372852) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.4670406) q[1];
sx q[1];
rz(-1.5853154) q[1];
sx q[1];
rz(0.0013703636) q[1];
x q[2];
rz(-2.8506365) q[3];
sx q[3];
rz(-0.97948631) q[3];
sx q[3];
rz(1.3475498) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.54258004) q[2];
sx q[2];
rz(-2.161669) q[2];
sx q[2];
rz(2.2194594) q[2];
rz(1.3433836) q[3];
sx q[3];
rz(-0.94815367) q[3];
sx q[3];
rz(1.5141727) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2682997) q[0];
sx q[0];
rz(-2.101185) q[0];
sx q[0];
rz(-2.701395) q[0];
rz(-1.531155) q[1];
sx q[1];
rz(-1.4819769) q[1];
sx q[1];
rz(0.24756113) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.089624) q[0];
sx q[0];
rz(-2.0735093) q[0];
sx q[0];
rz(2.6460656) q[0];
x q[1];
rz(-1.4804968) q[2];
sx q[2];
rz(-1.4019499) q[2];
sx q[2];
rz(2.894998) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.9369071) q[1];
sx q[1];
rz(-1.2720894) q[1];
sx q[1];
rz(-3.0644817) q[1];
rz(-pi) q[2];
rz(0.64856802) q[3];
sx q[3];
rz(-1.5463136) q[3];
sx q[3];
rz(-0.45059965) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.89746785) q[2];
sx q[2];
rz(-1.030587) q[2];
sx q[2];
rz(-1.4424651) q[2];
rz(-0.20081946) q[3];
sx q[3];
rz(-1.9229869) q[3];
sx q[3];
rz(2.0012746) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
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
rz(-0.15631256) q[0];
sx q[0];
rz(-0.15287481) q[0];
sx q[0];
rz(-2.5961764) q[0];
rz(3.0975869) q[1];
sx q[1];
rz(-0.017887201) q[1];
sx q[1];
rz(2.6652179) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2891083) q[0];
sx q[0];
rz(-1.8140287) q[0];
sx q[0];
rz(1.6575251) q[0];
x q[1];
rz(0.40619855) q[2];
sx q[2];
rz(-1.9425689) q[2];
sx q[2];
rz(0.68454725) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.093845) q[1];
sx q[1];
rz(-0.72870164) q[1];
sx q[1];
rz(-1.8456259) q[1];
rz(3.0073193) q[3];
sx q[3];
rz(-2.0760787) q[3];
sx q[3];
rz(-3.0485632) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.5215317) q[2];
sx q[2];
rz(-0.69053495) q[2];
sx q[2];
rz(2.7812092) q[2];
rz(-2.9195869) q[3];
sx q[3];
rz(-1.6111776) q[3];
sx q[3];
rz(1.41159) q[3];
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
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5905269) q[0];
sx q[0];
rz(-1.5883625) q[0];
sx q[0];
rz(1.5850413) q[0];
rz(0.12697728) q[1];
sx q[1];
rz(-1.304909) q[1];
sx q[1];
rz(3.051905) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4596649) q[0];
sx q[0];
rz(-2.0832422) q[0];
sx q[0];
rz(0.95109032) q[0];
rz(-pi) q[1];
rz(2.3562535) q[2];
sx q[2];
rz(-1.1113104) q[2];
sx q[2];
rz(1.6096514) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.9759102) q[1];
sx q[1];
rz(-1.0780562) q[1];
sx q[1];
rz(-1.8330281) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5406214) q[3];
sx q[3];
rz(-1.8716164) q[3];
sx q[3];
rz(1.6114637) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.0330607) q[2];
sx q[2];
rz(-0.33691275) q[2];
sx q[2];
rz(-2.8899371) q[2];
rz(3.1015977) q[3];
sx q[3];
rz(-2.7997041) q[3];
sx q[3];
rz(2.6778636) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7104178) q[0];
sx q[0];
rz(-0.091910563) q[0];
sx q[0];
rz(0.41496667) q[0];
rz(1.6522495) q[1];
sx q[1];
rz(-3.0840315) q[1];
sx q[1];
rz(-0.32589486) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1292757) q[0];
sx q[0];
rz(-1.1291478) q[0];
sx q[0];
rz(-0.4037767) q[0];
rz(-2.0891248) q[2];
sx q[2];
rz(-1.4926693) q[2];
sx q[2];
rz(0.6860439) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.4517026) q[1];
sx q[1];
rz(-0.57194114) q[1];
sx q[1];
rz(-2.6533935) q[1];
rz(-pi) q[2];
rz(3.065686) q[3];
sx q[3];
rz(-1.8733288) q[3];
sx q[3];
rz(-0.74750604) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.2275647) q[2];
sx q[2];
rz(-1.3332557) q[2];
sx q[2];
rz(-2.0951159) q[2];
rz(-2.922831) q[3];
sx q[3];
rz(-2.1121787) q[3];
sx q[3];
rz(-1.7441162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6662812) q[0];
sx q[0];
rz(-0.0175716) q[0];
sx q[0];
rz(-0.46324357) q[0];
rz(-0.26495588) q[1];
sx q[1];
rz(-3.1400561) q[1];
sx q[1];
rz(-1.6418246) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4508598) q[0];
sx q[0];
rz(-1.5520649) q[0];
sx q[0];
rz(-2.0035726) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6959571) q[2];
sx q[2];
rz(-1.5210954) q[2];
sx q[2];
rz(1.7222863) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.8275576) q[1];
sx q[1];
rz(-1.5998987) q[1];
sx q[1];
rz(1.3489789) q[1];
rz(-1.8345233) q[3];
sx q[3];
rz(-0.26104195) q[3];
sx q[3];
rz(-1.3021951) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.2396592) q[2];
sx q[2];
rz(-2.689211) q[2];
sx q[2];
rz(-1.8233914) q[2];
rz(3.1384595) q[3];
sx q[3];
rz(-1.9414709) q[3];
sx q[3];
rz(-1.9759294) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5779293) q[0];
sx q[0];
rz(-2.2888373) q[0];
sx q[0];
rz(0.71459115) q[0];
rz(2.928012) q[1];
sx q[1];
rz(-3.112308) q[1];
sx q[1];
rz(-1.2108796) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0830363) q[0];
sx q[0];
rz(-2.8175111) q[0];
sx q[0];
rz(1.9084318) q[0];
rz(-pi) q[1];
rz(2.4461652) q[2];
sx q[2];
rz(-1.8423242) q[2];
sx q[2];
rz(2.5144387) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.80187243) q[1];
sx q[1];
rz(-2.0959217) q[1];
sx q[1];
rz(0.32466356) q[1];
x q[2];
rz(-0.87745046) q[3];
sx q[3];
rz(-0.2444548) q[3];
sx q[3];
rz(-2.2742184) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.96656281) q[2];
sx q[2];
rz(-0.021447072) q[2];
sx q[2];
rz(0.19801298) q[2];
rz(-0.35061947) q[3];
sx q[3];
rz(-1.5712761) q[3];
sx q[3];
rz(1.9982136) q[3];
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
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38458934) q[0];
sx q[0];
rz(-0.92246711) q[0];
sx q[0];
rz(-2.5272227) q[0];
rz(0.059582926) q[1];
sx q[1];
rz(-3.0501084) q[1];
sx q[1];
rz(1.4245859) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2587268) q[0];
sx q[0];
rz(-0.75260163) q[0];
sx q[0];
rz(-1.2105788) q[0];
x q[1];
rz(-0.015083357) q[2];
sx q[2];
rz(-3.0102979) q[2];
sx q[2];
rz(2.7867336) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.9944133) q[1];
sx q[1];
rz(-0.1754027) q[1];
sx q[1];
rz(-2.0270623) q[1];
x q[2];
rz(-0.21496694) q[3];
sx q[3];
rz(-1.1720814) q[3];
sx q[3];
rz(0.48290184) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.1018579) q[2];
sx q[2];
rz(-2.8701344) q[2];
sx q[2];
rz(-0.52379918) q[2];
rz(-1.0456746) q[3];
sx q[3];
rz(-1.8664675) q[3];
sx q[3];
rz(0.4642134) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47281784) q[0];
sx q[0];
rz(-1.5325118) q[0];
sx q[0];
rz(-1.4788628) q[0];
rz(-1.4115903) q[1];
sx q[1];
rz(-0.23163207) q[1];
sx q[1];
rz(-3.0805265) q[1];
rz(-1.2712485) q[2];
sx q[2];
rz(-0.94956492) q[2];
sx q[2];
rz(1.9506394) q[2];
rz(-1.7674592) q[3];
sx q[3];
rz(-2.0296367) q[3];
sx q[3];
rz(2.2214132) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
