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
rz(0.76437104) q[0];
sx q[0];
rz(1.8005014) q[0];
sx q[0];
rz(10.330248) q[0];
rz(-1.634693) q[1];
sx q[1];
rz(-0.96087471) q[1];
sx q[1];
rz(-1.5009872) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.95661058) q[0];
sx q[0];
rz(-2.2652103) q[0];
sx q[0];
rz(2.882769) q[0];
rz(-pi) q[1];
rz(-1.1314279) q[2];
sx q[2];
rz(-1.7305534) q[2];
sx q[2];
rz(1.8607832) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.562362) q[1];
sx q[1];
rz(-0.048431245) q[1];
sx q[1];
rz(2.6391451) q[1];
x q[2];
rz(1.0085868) q[3];
sx q[3];
rz(-1.2744181) q[3];
sx q[3];
rz(-2.5135771) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.91244873) q[2];
sx q[2];
rz(-1.7488166) q[2];
sx q[2];
rz(2.8196715) q[2];
rz(0.95494444) q[3];
sx q[3];
rz(-2.3117282) q[3];
sx q[3];
rz(1.1436852) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(-0.8973273) q[0];
sx q[0];
rz(-1.0538415) q[0];
sx q[0];
rz(0.51189297) q[0];
rz(-1.5882675) q[1];
sx q[1];
rz(-0.48738185) q[1];
sx q[1];
rz(2.0194676) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.666709) q[0];
sx q[0];
rz(-2.8132952) q[0];
sx q[0];
rz(-0.31995456) q[0];
rz(-pi) q[1];
x q[1];
rz(2.000359) q[2];
sx q[2];
rz(-1.5341752) q[2];
sx q[2];
rz(0.0096732339) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.8583688) q[1];
sx q[1];
rz(-2.5852381) q[1];
sx q[1];
rz(2.2409641) q[1];
x q[2];
rz(-2.0418704) q[3];
sx q[3];
rz(-2.2212914) q[3];
sx q[3];
rz(2.4471403) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.7684218) q[2];
sx q[2];
rz(-1.827652) q[2];
sx q[2];
rz(2.6785417) q[2];
rz(2.566346) q[3];
sx q[3];
rz(-1.7762215) q[3];
sx q[3];
rz(-0.099253207) q[3];
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
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7333882) q[0];
sx q[0];
rz(-2.2938804) q[0];
sx q[0];
rz(-2.7114765) q[0];
rz(-0.46562132) q[1];
sx q[1];
rz(-0.72048134) q[1];
sx q[1];
rz(-0.63708416) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89403215) q[0];
sx q[0];
rz(-1.5567509) q[0];
sx q[0];
rz(1.5298903) q[0];
rz(1.5445721) q[2];
sx q[2];
rz(-2.3560212) q[2];
sx q[2];
rz(2.3674813) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.3780669) q[1];
sx q[1];
rz(-2.3664306) q[1];
sx q[1];
rz(-2.7833392) q[1];
rz(-1.9608843) q[3];
sx q[3];
rz(-2.0566787) q[3];
sx q[3];
rz(2.7220243) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.35885262) q[2];
sx q[2];
rz(-1.0461067) q[2];
sx q[2];
rz(-2.4165238) q[2];
rz(1.8105761) q[3];
sx q[3];
rz(-1.015181) q[3];
sx q[3];
rz(-2.5031808) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26938874) q[0];
sx q[0];
rz(-1.4545472) q[0];
sx q[0];
rz(1.4440906) q[0];
rz(3.0929502) q[1];
sx q[1];
rz(-1.8154058) q[1];
sx q[1];
rz(-2.8526502) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.671593) q[0];
sx q[0];
rz(-0.44737383) q[0];
sx q[0];
rz(0.98487052) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7213351) q[2];
sx q[2];
rz(-1.5427914) q[2];
sx q[2];
rz(-1.752252) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.7270941) q[1];
sx q[1];
rz(-1.7448404) q[1];
sx q[1];
rz(3.0309309) q[1];
rz(-pi) q[2];
rz(-1.6252609) q[3];
sx q[3];
rz(-0.90622444) q[3];
sx q[3];
rz(1.2913454) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.284953) q[2];
sx q[2];
rz(-0.53264561) q[2];
sx q[2];
rz(0.3802158) q[2];
rz(-0.63816655) q[3];
sx q[3];
rz(-1.2297945) q[3];
sx q[3];
rz(1.1285454) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3087092) q[0];
sx q[0];
rz(-0.55279151) q[0];
sx q[0];
rz(-1.0943476) q[0];
rz(0.87567466) q[1];
sx q[1];
rz(-2.5119669) q[1];
sx q[1];
rz(-2.0424776) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7849784) q[0];
sx q[0];
rz(-1.7081424) q[0];
sx q[0];
rz(-2.7524968) q[0];
rz(-2.7988222) q[2];
sx q[2];
rz(-1.2452599) q[2];
sx q[2];
rz(-2.8936762) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.588394) q[1];
sx q[1];
rz(-2.1171682) q[1];
sx q[1];
rz(-2.9562034) q[1];
rz(-pi) q[2];
rz(-1.7775675) q[3];
sx q[3];
rz(-2.4768157) q[3];
sx q[3];
rz(2.1867276) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.8984453) q[2];
sx q[2];
rz(-1.3546319) q[2];
sx q[2];
rz(-0.97664991) q[2];
rz(0.53168932) q[3];
sx q[3];
rz(-0.44638005) q[3];
sx q[3];
rz(-0.71438742) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0584745) q[0];
sx q[0];
rz(-2.0396905) q[0];
sx q[0];
rz(1.8699159) q[0];
rz(-2.7187128) q[1];
sx q[1];
rz(-1.6115178) q[1];
sx q[1];
rz(1.6082825) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6278224) q[0];
sx q[0];
rz(-1.5830399) q[0];
sx q[0];
rz(-1.6374541) q[0];
rz(-pi) q[1];
rz(1.6663867) q[2];
sx q[2];
rz(-3.0431192) q[2];
sx q[2];
rz(-2.5661039) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.3376006) q[1];
sx q[1];
rz(-1.2678267) q[1];
sx q[1];
rz(1.5125809) q[1];
x q[2];
rz(1.815237) q[3];
sx q[3];
rz(-0.77271739) q[3];
sx q[3];
rz(2.6561007) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.61234683) q[2];
sx q[2];
rz(-1.8491448) q[2];
sx q[2];
rz(0.38522729) q[2];
rz(3.0569844) q[3];
sx q[3];
rz(-0.43360964) q[3];
sx q[3];
rz(-0.87578526) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2200634) q[0];
sx q[0];
rz(-1.055287) q[0];
sx q[0];
rz(0.39749843) q[0];
rz(-2.6037604) q[1];
sx q[1];
rz(-0.42306867) q[1];
sx q[1];
rz(0.0040815512) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5631112) q[0];
sx q[0];
rz(-2.3627776) q[0];
sx q[0];
rz(-1.388873) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2572558) q[2];
sx q[2];
rz(-1.8914273) q[2];
sx q[2];
rz(0.31242958) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.7023125) q[1];
sx q[1];
rz(-1.951362) q[1];
sx q[1];
rz(0.61813942) q[1];
rz(2.3621906) q[3];
sx q[3];
rz(-2.2891697) q[3];
sx q[3];
rz(1.8441594) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.5092545) q[2];
sx q[2];
rz(-1.1804487) q[2];
sx q[2];
rz(-2.9285367) q[2];
rz(0.80900711) q[3];
sx q[3];
rz(-3.1000948) q[3];
sx q[3];
rz(1.8607148) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0179366) q[0];
sx q[0];
rz(-1.9391215) q[0];
sx q[0];
rz(-1.1630455) q[0];
rz(-2.842438) q[1];
sx q[1];
rz(-1.9367633) q[1];
sx q[1];
rz(0.97602731) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7928182) q[0];
sx q[0];
rz(-0.26991699) q[0];
sx q[0];
rz(-2.7151373) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7485745) q[2];
sx q[2];
rz(-1.7209098) q[2];
sx q[2];
rz(0.77265384) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.4055172) q[1];
sx q[1];
rz(-1.8092522) q[1];
sx q[1];
rz(1.2082491) q[1];
x q[2];
rz(0.42404948) q[3];
sx q[3];
rz(-1.7982139) q[3];
sx q[3];
rz(-3.1229179) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.8396847) q[2];
sx q[2];
rz(-0.32543108) q[2];
sx q[2];
rz(-3.0111266) q[2];
rz(1.8179551) q[3];
sx q[3];
rz(-1.8925083) q[3];
sx q[3];
rz(0.29449335) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.94006222) q[0];
sx q[0];
rz(-0.37684965) q[0];
sx q[0];
rz(1.4991722) q[0];
rz(0.54010737) q[1];
sx q[1];
rz(-2.3441548) q[1];
sx q[1];
rz(-1.7971719) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3279874) q[0];
sx q[0];
rz(-1.9623956) q[0];
sx q[0];
rz(2.8285976) q[0];
rz(-2.1189519) q[2];
sx q[2];
rz(-1.022199) q[2];
sx q[2];
rz(0.61778574) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.1229748) q[1];
sx q[1];
rz(-2.5135165) q[1];
sx q[1];
rz(2.0043892) q[1];
rz(-pi) q[2];
rz(-2.2910883) q[3];
sx q[3];
rz(-2.0652186) q[3];
sx q[3];
rz(2.0317584) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.53238791) q[2];
sx q[2];
rz(-2.8049073) q[2];
sx q[2];
rz(1.4538291) q[2];
rz(-0.42803556) q[3];
sx q[3];
rz(-1.5902218) q[3];
sx q[3];
rz(-1.154703) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55353272) q[0];
sx q[0];
rz(-2.689671) q[0];
sx q[0];
rz(1.4916627) q[0];
rz(0.55117575) q[1];
sx q[1];
rz(-1.0124413) q[1];
sx q[1];
rz(-0.59250441) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.036464036) q[0];
sx q[0];
rz(-1.8499287) q[0];
sx q[0];
rz(-1.3238086) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5985322) q[2];
sx q[2];
rz(-2.0351699) q[2];
sx q[2];
rz(0.65242243) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.4514102) q[1];
sx q[1];
rz(-1.8333149) q[1];
sx q[1];
rz(-0.39509829) q[1];
rz(-0.49919923) q[3];
sx q[3];
rz(-1.5268402) q[3];
sx q[3];
rz(1.4756965) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.7117915) q[2];
sx q[2];
rz(-2.1000803) q[2];
sx q[2];
rz(3.0832624) q[2];
rz(-0.37060261) q[3];
sx q[3];
rz(-2.8648418) q[3];
sx q[3];
rz(-3.1378194) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8534828) q[0];
sx q[0];
rz(-1.7848889) q[0];
sx q[0];
rz(-3.0369192) q[0];
rz(1.7755605) q[1];
sx q[1];
rz(-0.80148253) q[1];
sx q[1];
rz(-2.186224) q[1];
rz(-2.7544207) q[2];
sx q[2];
rz(-2.5324814) q[2];
sx q[2];
rz(-1.0003288) q[2];
rz(-2.1128863) q[3];
sx q[3];
rz(-1.1114612) q[3];
sx q[3];
rz(-2.3041861) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
