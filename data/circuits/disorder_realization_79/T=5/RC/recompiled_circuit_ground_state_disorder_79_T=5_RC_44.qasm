OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.3747342) q[0];
sx q[0];
rz(-2.2280966) q[0];
sx q[0];
rz(-1.7537533) q[0];
rz(-2.9517382) q[1];
sx q[1];
rz(-1.318765) q[1];
sx q[1];
rz(2.8554822) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0510584) q[0];
sx q[0];
rz(-1.5947475) q[0];
sx q[0];
rz(0.46047237) q[0];
rz(-pi) q[1];
rz(1.6880116) q[2];
sx q[2];
rz(-2.5273537) q[2];
sx q[2];
rz(0.87488824) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.3158947) q[1];
sx q[1];
rz(-1.311439) q[1];
sx q[1];
rz(0.56866912) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6872365) q[3];
sx q[3];
rz(-0.44346657) q[3];
sx q[3];
rz(0.93772353) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.046595786) q[2];
sx q[2];
rz(-1.4904212) q[2];
sx q[2];
rz(0.68520927) q[2];
rz(-1.4677216) q[3];
sx q[3];
rz(-2.0592212) q[3];
sx q[3];
rz(2.8732324) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14107038) q[0];
sx q[0];
rz(-1.5474316) q[0];
sx q[0];
rz(-0.94456124) q[0];
rz(-0.073307723) q[1];
sx q[1];
rz(-0.83275515) q[1];
sx q[1];
rz(-1.8836969) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2908173) q[0];
sx q[0];
rz(-0.75277872) q[0];
sx q[0];
rz(2.2673009) q[0];
rz(-pi) q[1];
rz(1.3665359) q[2];
sx q[2];
rz(-2.0957207) q[2];
sx q[2];
rz(-1.7634942) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.98698178) q[1];
sx q[1];
rz(-1.4263337) q[1];
sx q[1];
rz(-0.39244907) q[1];
x q[2];
rz(0.03713921) q[3];
sx q[3];
rz(-1.5243153) q[3];
sx q[3];
rz(0.4872191) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.8153317) q[2];
sx q[2];
rz(-1.7504642) q[2];
sx q[2];
rz(-2.1168671) q[2];
rz(2.4552086) q[3];
sx q[3];
rz(-2.4095583) q[3];
sx q[3];
rz(-2.2288403) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.82655418) q[0];
sx q[0];
rz(-1.2737561) q[0];
sx q[0];
rz(-2.3231373) q[0];
rz(-1.2089027) q[1];
sx q[1];
rz(-1.6345638) q[1];
sx q[1];
rz(-0.85982972) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8743961) q[0];
sx q[0];
rz(-0.98496297) q[0];
sx q[0];
rz(-0.39975014) q[0];
rz(-0.53162127) q[2];
sx q[2];
rz(-0.79382703) q[2];
sx q[2];
rz(0.17228157) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.7581186) q[1];
sx q[1];
rz(-2.2920597) q[1];
sx q[1];
rz(-0.63754941) q[1];
rz(-pi) q[2];
rz(1.1016261) q[3];
sx q[3];
rz(-0.2519603) q[3];
sx q[3];
rz(-2.2445238) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.7439421) q[2];
sx q[2];
rz(-2.8068779) q[2];
sx q[2];
rz(-3.0261377) q[2];
rz(0.3913106) q[3];
sx q[3];
rz(-1.0262998) q[3];
sx q[3];
rz(1.1439884) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7855969) q[0];
sx q[0];
rz(-1.1859256) q[0];
sx q[0];
rz(2.9010229) q[0];
rz(-1.7754414) q[1];
sx q[1];
rz(-1.6484478) q[1];
sx q[1];
rz(-1.8919401) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25729733) q[0];
sx q[0];
rz(-0.69778555) q[0];
sx q[0];
rz(-1.1691514) q[0];
rz(-pi) q[1];
rz(-0.95486887) q[2];
sx q[2];
rz(-1.2754585) q[2];
sx q[2];
rz(1.5463285) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.27423672) q[1];
sx q[1];
rz(-0.70046591) q[1];
sx q[1];
rz(-2.7718616) q[1];
rz(-pi) q[2];
rz(-2.4576538) q[3];
sx q[3];
rz(-1.9762282) q[3];
sx q[3];
rz(0.33567521) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.8852641) q[2];
sx q[2];
rz(-0.54577959) q[2];
sx q[2];
rz(-1.3679999) q[2];
rz(0.3041501) q[3];
sx q[3];
rz(-1.5658028) q[3];
sx q[3];
rz(2.4450541) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.91365415) q[0];
sx q[0];
rz(-2.9253687) q[0];
sx q[0];
rz(-2.6044593) q[0];
rz(-2.3917603) q[1];
sx q[1];
rz(-2.0593819) q[1];
sx q[1];
rz(0.73712635) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0219624) q[0];
sx q[0];
rz(-3.0263794) q[0];
sx q[0];
rz(2.553279) q[0];
rz(3.0546435) q[2];
sx q[2];
rz(-1.7450404) q[2];
sx q[2];
rz(-2.6373495) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.8512042) q[1];
sx q[1];
rz(-1.8294639) q[1];
sx q[1];
rz(1.438739) q[1];
x q[2];
rz(1.4786105) q[3];
sx q[3];
rz(-1.2488197) q[3];
sx q[3];
rz(0.76065242) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.4445112) q[2];
sx q[2];
rz(-0.35327521) q[2];
sx q[2];
rz(-2.9949761) q[2];
rz(-1.692449) q[3];
sx q[3];
rz(-1.861898) q[3];
sx q[3];
rz(2.6176738) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.523282) q[0];
sx q[0];
rz(-2.1250516) q[0];
sx q[0];
rz(-1.4215533) q[0];
rz(1.2591741) q[1];
sx q[1];
rz(-2.4651395) q[1];
sx q[1];
rz(-2.7377985) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.97966563) q[0];
sx q[0];
rz(-0.84507361) q[0];
sx q[0];
rz(0.37958522) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.25713276) q[2];
sx q[2];
rz(-1.3916573) q[2];
sx q[2];
rz(0.15878294) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.660259) q[1];
sx q[1];
rz(-1.0254045) q[1];
sx q[1];
rz(-0.86472269) q[1];
rz(-2.850612) q[3];
sx q[3];
rz(-2.432219) q[3];
sx q[3];
rz(-3.0090699) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.28973618) q[2];
sx q[2];
rz(-1.4137665) q[2];
sx q[2];
rz(0.16656052) q[2];
rz(-1.5038495) q[3];
sx q[3];
rz(-0.47347355) q[3];
sx q[3];
rz(3.04305) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2378167) q[0];
sx q[0];
rz(-0.38919583) q[0];
sx q[0];
rz(-0.87752262) q[0];
rz(0.96218836) q[1];
sx q[1];
rz(-1.1233556) q[1];
sx q[1];
rz(-0.35433623) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41397944) q[0];
sx q[0];
rz(-0.2383735) q[0];
sx q[0];
rz(0.075043126) q[0];
rz(-pi) q[1];
rz(1.4260068) q[2];
sx q[2];
rz(-1.5476523) q[2];
sx q[2];
rz(2.1320081) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.42592624) q[1];
sx q[1];
rz(-0.86326438) q[1];
sx q[1];
rz(2.3051579) q[1];
x q[2];
rz(-0.50994344) q[3];
sx q[3];
rz(-1.5048319) q[3];
sx q[3];
rz(-3.0086781) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.3366036) q[2];
sx q[2];
rz(-2.4038834) q[2];
sx q[2];
rz(-2.0965516) q[2];
rz(-0.40237829) q[3];
sx q[3];
rz(-1.9741524) q[3];
sx q[3];
rz(1.6654525) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7104915) q[0];
sx q[0];
rz(-3.0248108) q[0];
sx q[0];
rz(-2.2905599) q[0];
rz(-0.35375133) q[1];
sx q[1];
rz(-1.7314792) q[1];
sx q[1];
rz(0.85465777) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2609278) q[0];
sx q[0];
rz(-1.3890084) q[0];
sx q[0];
rz(2.8559235) q[0];
rz(0.53422569) q[2];
sx q[2];
rz(-0.28817982) q[2];
sx q[2];
rz(-2.2792918) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.51616299) q[1];
sx q[1];
rz(-0.99402797) q[1];
sx q[1];
rz(-1.7712084) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7386916) q[3];
sx q[3];
rz(-1.6454927) q[3];
sx q[3];
rz(2.3011329) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.9510368) q[2];
sx q[2];
rz(-0.64728105) q[2];
sx q[2];
rz(2.5212042) q[2];
rz(0.22647151) q[3];
sx q[3];
rz(-1.7445931) q[3];
sx q[3];
rz(2.5605719) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9613551) q[0];
sx q[0];
rz(-1.2018452) q[0];
sx q[0];
rz(-3.1264547) q[0];
rz(2.119428) q[1];
sx q[1];
rz(-2.731555) q[1];
sx q[1];
rz(-1.2000363) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.43754345) q[0];
sx q[0];
rz(-1.3945197) q[0];
sx q[0];
rz(-2.2793819) q[0];
x q[1];
rz(-0.59320088) q[2];
sx q[2];
rz(-1.6199379) q[2];
sx q[2];
rz(1.0047439) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.2851499) q[1];
sx q[1];
rz(-2.0822444) q[1];
sx q[1];
rz(-3.1072381) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3288378) q[3];
sx q[3];
rz(-2.4852607) q[3];
sx q[3];
rz(1.1876653) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.52593652) q[2];
sx q[2];
rz(-0.91292149) q[2];
sx q[2];
rz(-2.9021662) q[2];
rz(3.0827403) q[3];
sx q[3];
rz(-2.3299496) q[3];
sx q[3];
rz(-2.8656901) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9489768) q[0];
sx q[0];
rz(-1.9276351) q[0];
sx q[0];
rz(2.1666727) q[0];
rz(0.67543593) q[1];
sx q[1];
rz(-0.91921872) q[1];
sx q[1];
rz(2.1598037) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4083207) q[0];
sx q[0];
rz(-1.3106723) q[0];
sx q[0];
rz(0.97113804) q[0];
rz(-pi) q[1];
rz(-2.9862829) q[2];
sx q[2];
rz(-2.7255645) q[2];
sx q[2];
rz(1.95259) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.1002994) q[1];
sx q[1];
rz(-0.86242341) q[1];
sx q[1];
rz(-2.7977976) q[1];
rz(1.6493787) q[3];
sx q[3];
rz(-2.693697) q[3];
sx q[3];
rz(-3.0495897) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.5741817) q[2];
sx q[2];
rz(-2.1705706) q[2];
sx q[2];
rz(1.2249472) q[2];
rz(-2.775906) q[3];
sx q[3];
rz(-2.1920125) q[3];
sx q[3];
rz(2.001781) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.097261978) q[0];
sx q[0];
rz(-1.8987569) q[0];
sx q[0];
rz(1.4415997) q[0];
rz(0.080009566) q[1];
sx q[1];
rz(-1.3792104) q[1];
sx q[1];
rz(3.1389799) q[1];
rz(-1.710113) q[2];
sx q[2];
rz(-2.0124544) q[2];
sx q[2];
rz(-2.6228946) q[2];
rz(-3.0946685) q[3];
sx q[3];
rz(-0.45294807) q[3];
sx q[3];
rz(3.1035085) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
