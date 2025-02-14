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
rz(1.7914766) q[0];
sx q[0];
rz(-0.67576367) q[0];
sx q[0];
rz(-0.0078553353) q[0];
rz(0.76771843) q[1];
sx q[1];
rz(-1.5239198) q[1];
sx q[1];
rz(2.89892) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4123101) q[0];
sx q[0];
rz(-0.72107102) q[0];
sx q[0];
rz(1.7626053) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9353634) q[2];
sx q[2];
rz(-2.0040214) q[2];
sx q[2];
rz(-1.3490167) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.3625177) q[1];
sx q[1];
rz(-1.9240161) q[1];
sx q[1];
rz(-1.5959969) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4759329) q[3];
sx q[3];
rz(-1.7659617) q[3];
sx q[3];
rz(-1.631537) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0365888) q[2];
sx q[2];
rz(-1.236311) q[2];
sx q[2];
rz(-2.429602) q[2];
rz(-2.8365734) q[3];
sx q[3];
rz(-0.22235338) q[3];
sx q[3];
rz(-3.0960848) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5104093) q[0];
sx q[0];
rz(-1.864186) q[0];
sx q[0];
rz(-1.7741868) q[0];
rz(-0.039904682) q[1];
sx q[1];
rz(-2.3343562) q[1];
sx q[1];
rz(2.5423999) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0586226) q[0];
sx q[0];
rz(-0.48958594) q[0];
sx q[0];
rz(0.24607308) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.97936432) q[2];
sx q[2];
rz(-2.2970846) q[2];
sx q[2];
rz(-1.3853816) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.6015905) q[1];
sx q[1];
rz(-2.4036602) q[1];
sx q[1];
rz(-0.5088536) q[1];
rz(-pi) q[2];
rz(2.8017524) q[3];
sx q[3];
rz(-0.71469864) q[3];
sx q[3];
rz(0.13091892) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.082531884) q[2];
sx q[2];
rz(-0.56190562) q[2];
sx q[2];
rz(1.833029) q[2];
rz(-0.023905309) q[3];
sx q[3];
rz(-1.6126361) q[3];
sx q[3];
rz(-1.5786952) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6858653) q[0];
sx q[0];
rz(-2.4754334) q[0];
sx q[0];
rz(-2.300793) q[0];
rz(1.0288382) q[1];
sx q[1];
rz(-2.7327171) q[1];
sx q[1];
rz(-1.4604481) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59305916) q[0];
sx q[0];
rz(-1.7139627) q[0];
sx q[0];
rz(-0.4645223) q[0];
rz(-pi) q[1];
rz(1.9424136) q[2];
sx q[2];
rz(-1.1254246) q[2];
sx q[2];
rz(-2.1567313) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.3884133) q[1];
sx q[1];
rz(-2.3824661) q[1];
sx q[1];
rz(0.49465804) q[1];
x q[2];
rz(0.60734235) q[3];
sx q[3];
rz(-1.6173956) q[3];
sx q[3];
rz(2.6452015) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.8847522) q[2];
sx q[2];
rz(-1.2105056) q[2];
sx q[2];
rz(-2.9849226) q[2];
rz(-1.6414075) q[3];
sx q[3];
rz(-1.2431966) q[3];
sx q[3];
rz(-2.959804) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0750065) q[0];
sx q[0];
rz(-1.8445419) q[0];
sx q[0];
rz(-0.080667607) q[0];
rz(1.1219885) q[1];
sx q[1];
rz(-2.5009218) q[1];
sx q[1];
rz(1.5871619) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4951404) q[0];
sx q[0];
rz(-2.3123154) q[0];
sx q[0];
rz(-0.089228169) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8055259) q[2];
sx q[2];
rz(-3.0185351) q[2];
sx q[2];
rz(-0.50257909) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.746628) q[1];
sx q[1];
rz(-1.0735128) q[1];
sx q[1];
rz(-0.73594603) q[1];
x q[2];
rz(-2.5347363) q[3];
sx q[3];
rz(-1.8491866) q[3];
sx q[3];
rz(-0.75137072) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.88177219) q[2];
sx q[2];
rz(-2.0505003) q[2];
sx q[2];
rz(2.1176977) q[2];
rz(2.3648868) q[3];
sx q[3];
rz(-1.1506162) q[3];
sx q[3];
rz(-1.0030494) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7316932) q[0];
sx q[0];
rz(-2.2875146) q[0];
sx q[0];
rz(2.5140629) q[0];
rz(-0.69951406) q[1];
sx q[1];
rz(-2.1414089) q[1];
sx q[1];
rz(-0.90739179) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9357244) q[0];
sx q[0];
rz(-2.3875934) q[0];
sx q[0];
rz(3.0831017) q[0];
x q[1];
rz(-2.9512915) q[2];
sx q[2];
rz(-1.7277328) q[2];
sx q[2];
rz(-2.4630594) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.40014711) q[1];
sx q[1];
rz(-0.65175599) q[1];
sx q[1];
rz(0.94237997) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0346342) q[3];
sx q[3];
rz(-1.5201293) q[3];
sx q[3];
rz(-1.9127653) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.13713914) q[2];
sx q[2];
rz(-1.2848102) q[2];
sx q[2];
rz(-0.84490204) q[2];
rz(0.6984624) q[3];
sx q[3];
rz(-2.4013077) q[3];
sx q[3];
rz(-2.3499878) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.38248211) q[0];
sx q[0];
rz(-1.6805205) q[0];
sx q[0];
rz(1.8735029) q[0];
rz(1.5269439) q[1];
sx q[1];
rz(-2.0497649) q[1];
sx q[1];
rz(-0.39438927) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1067692) q[0];
sx q[0];
rz(-0.64578694) q[0];
sx q[0];
rz(1.4950947) q[0];
rz(-pi) q[1];
rz(-0.092124513) q[2];
sx q[2];
rz(-0.73706223) q[2];
sx q[2];
rz(1.2110405) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.7463267) q[1];
sx q[1];
rz(-1.6155287) q[1];
sx q[1];
rz(-0.6529863) q[1];
rz(2.9378618) q[3];
sx q[3];
rz(-2.7739848) q[3];
sx q[3];
rz(2.550761) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.4322728) q[2];
sx q[2];
rz(-3.0797112) q[2];
sx q[2];
rz(0.56149948) q[2];
rz(1.1310486) q[3];
sx q[3];
rz(-0.80315042) q[3];
sx q[3];
rz(2.6530182) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5442218) q[0];
sx q[0];
rz(-2.9591296) q[0];
sx q[0];
rz(-2.9083948) q[0];
rz(3.0274262) q[1];
sx q[1];
rz(-2.3515067) q[1];
sx q[1];
rz(-0.17359576) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.9812935) q[0];
sx q[0];
rz(-2.3620785) q[0];
sx q[0];
rz(2.1008089) q[0];
rz(2.7596682) q[2];
sx q[2];
rz(-1.9878824) q[2];
sx q[2];
rz(-2.797319) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.0376589) q[1];
sx q[1];
rz(-1.7915495) q[1];
sx q[1];
rz(-2.9751865) q[1];
x q[2];
rz(1.5585445) q[3];
sx q[3];
rz(-1.3821332) q[3];
sx q[3];
rz(2.0607299) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.0357828) q[2];
sx q[2];
rz(-2.1828987) q[2];
sx q[2];
rz(-0.85476533) q[2];
rz(-0.0828951) q[3];
sx q[3];
rz(-1.1707183) q[3];
sx q[3];
rz(0.054072592) q[3];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6524803) q[0];
sx q[0];
rz(-1.880045) q[0];
sx q[0];
rz(1.1789119) q[0];
rz(2.1401801) q[1];
sx q[1];
rz(-1.8779571) q[1];
sx q[1];
rz(2.9875535) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4074041) q[0];
sx q[0];
rz(-1.2449236) q[0];
sx q[0];
rz(-2.4209321) q[0];
rz(0.78275795) q[2];
sx q[2];
rz(-2.2680757) q[2];
sx q[2];
rz(0.00072602206) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.5370085) q[1];
sx q[1];
rz(-1.8698543) q[1];
sx q[1];
rz(1.0826151) q[1];
rz(2.8518744) q[3];
sx q[3];
rz(-1.8236092) q[3];
sx q[3];
rz(-0.022965206) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.4607294) q[2];
sx q[2];
rz(-1.0793842) q[2];
sx q[2];
rz(2.4737849) q[2];
rz(1.7364511) q[3];
sx q[3];
rz(-1.3677771) q[3];
sx q[3];
rz(2.5051129) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3619096) q[0];
sx q[0];
rz(-1.6661665) q[0];
sx q[0];
rz(-2.5478126) q[0];
rz(1.8608015) q[1];
sx q[1];
rz(-2.180876) q[1];
sx q[1];
rz(1.2978172) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0827328) q[0];
sx q[0];
rz(-1.642859) q[0];
sx q[0];
rz(2.5565992) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.62784451) q[2];
sx q[2];
rz(-1.825001) q[2];
sx q[2];
rz(2.0275627) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.1480963) q[1];
sx q[1];
rz(-1.9283174) q[1];
sx q[1];
rz(1.7054249) q[1];
rz(-0.16420096) q[3];
sx q[3];
rz(-0.81402011) q[3];
sx q[3];
rz(-2.9661953) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.7456776) q[2];
sx q[2];
rz(-3.119645) q[2];
sx q[2];
rz(1.5094666) q[2];
rz(0.11219003) q[3];
sx q[3];
rz(-1.0271007) q[3];
sx q[3];
rz(-1.7621015) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6701732) q[0];
sx q[0];
rz(-2.1355974) q[0];
sx q[0];
rz(2.8759586) q[0];
rz(0.98948014) q[1];
sx q[1];
rz(-1.2629291) q[1];
sx q[1];
rz(-0.33214733) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.51458329) q[0];
sx q[0];
rz(-0.3334612) q[0];
sx q[0];
rz(-1.5824806) q[0];
x q[1];
rz(0.99224706) q[2];
sx q[2];
rz(-1.3253115) q[2];
sx q[2];
rz(3.0089889) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.2492318) q[1];
sx q[1];
rz(-0.3438102) q[1];
sx q[1];
rz(1.4268141) q[1];
rz(-pi) q[2];
rz(-3.0933558) q[3];
sx q[3];
rz(-1.3833481) q[3];
sx q[3];
rz(-1.6316044) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.0418479) q[2];
sx q[2];
rz(-1.0724649) q[2];
sx q[2];
rz(-0.15038807) q[2];
rz(0.8485052) q[3];
sx q[3];
rz(-2.7917807) q[3];
sx q[3];
rz(1.8457886) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.083241845) q[0];
sx q[0];
rz(-1.8235089) q[0];
sx q[0];
rz(1.9705809) q[0];
rz(1.8104443) q[1];
sx q[1];
rz(-2.34453) q[1];
sx q[1];
rz(-1.1042368) q[1];
rz(-0.83126478) q[2];
sx q[2];
rz(-1.6192042) q[2];
sx q[2];
rz(-1.529196) q[2];
rz(-1.5965309) q[3];
sx q[3];
rz(-0.86418695) q[3];
sx q[3];
rz(1.0623111) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
