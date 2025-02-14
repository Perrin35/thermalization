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
rz(-2.6391368) q[0];
sx q[0];
rz(-1.1075736) q[0];
sx q[0];
rz(-1.3753608) q[0];
rz(-3.4604685) q[1];
sx q[1];
rz(4.7825216) q[1];
sx q[1];
rz(8.6934269) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36424822) q[0];
sx q[0];
rz(-0.70667446) q[0];
sx q[0];
rz(-2.0752491) q[0];
rz(-pi) q[1];
rz(0.66197239) q[2];
sx q[2];
rz(-1.5675002) q[2];
sx q[2];
rz(2.4087083) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.3870942) q[1];
sx q[1];
rz(-1.525314) q[1];
sx q[1];
rz(-0.78630739) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.5238138) q[3];
sx q[3];
rz(-2.9546545) q[3];
sx q[3];
rz(2.4863249) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.88087624) q[2];
sx q[2];
rz(-1.2141576) q[2];
sx q[2];
rz(-0.94240776) q[2];
rz(-0.85855329) q[3];
sx q[3];
rz(-0.61045727) q[3];
sx q[3];
rz(2.4220991) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.079916) q[0];
sx q[0];
rz(-2.583857) q[0];
sx q[0];
rz(-3.0829561) q[0];
rz(-3.0851641) q[1];
sx q[1];
rz(-1.3899048) q[1];
sx q[1];
rz(1.1172392) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4898713) q[0];
sx q[0];
rz(-0.84017032) q[0];
sx q[0];
rz(-2.1108759) q[0];
rz(1.9675206) q[2];
sx q[2];
rz(-2.8182321) q[2];
sx q[2];
rz(-0.67229313) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.62395019) q[1];
sx q[1];
rz(-0.25282598) q[1];
sx q[1];
rz(-1.235591) q[1];
x q[2];
rz(1.6916709) q[3];
sx q[3];
rz(-2.0434922) q[3];
sx q[3];
rz(-2.8877088) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.3650018) q[2];
sx q[2];
rz(-2.8114909) q[2];
sx q[2];
rz(-2.7351725) q[2];
rz(-0.31600076) q[3];
sx q[3];
rz(-1.4603115) q[3];
sx q[3];
rz(-2.3587904) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86421788) q[0];
sx q[0];
rz(-0.45806956) q[0];
sx q[0];
rz(-1.8055441) q[0];
rz(2.863073) q[1];
sx q[1];
rz(-1.7428935) q[1];
sx q[1];
rz(2.1727402) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.028038) q[0];
sx q[0];
rz(-0.73532205) q[0];
sx q[0];
rz(0.57082973) q[0];
x q[1];
rz(2.5843262) q[2];
sx q[2];
rz(-2.2695978) q[2];
sx q[2];
rz(3.0521309) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.90889764) q[1];
sx q[1];
rz(-1.2746342) q[1];
sx q[1];
rz(-3.0975268) q[1];
x q[2];
rz(-0.28964759) q[3];
sx q[3];
rz(-0.74637369) q[3];
sx q[3];
rz(-2.1044855) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.1442673) q[2];
sx q[2];
rz(-2.6582025) q[2];
sx q[2];
rz(2.9664795) q[2];
rz(-1.8363606) q[3];
sx q[3];
rz(-2.1988018) q[3];
sx q[3];
rz(-1.85873) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1433732) q[0];
sx q[0];
rz(-2.0046736) q[0];
sx q[0];
rz(2.112222) q[0];
rz(-0.81946212) q[1];
sx q[1];
rz(-1.9792604) q[1];
sx q[1];
rz(-2.3447461) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5982847) q[0];
sx q[0];
rz(-0.48439041) q[0];
sx q[0];
rz(-0.422307) q[0];
rz(-pi) q[1];
x q[1];
rz(0.0060002319) q[2];
sx q[2];
rz(-2.0987392) q[2];
sx q[2];
rz(2.1354298) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.0174344) q[1];
sx q[1];
rz(-2.7982111) q[1];
sx q[1];
rz(0.92059532) q[1];
x q[2];
rz(2.4791777) q[3];
sx q[3];
rz(-1.2235421) q[3];
sx q[3];
rz(1.9062689) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.13021079) q[2];
sx q[2];
rz(-0.036616651) q[2];
sx q[2];
rz(1.9127362) q[2];
rz(0.43904385) q[3];
sx q[3];
rz(-1.5758347) q[3];
sx q[3];
rz(-1.2307833) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48290408) q[0];
sx q[0];
rz(-1.6084325) q[0];
sx q[0];
rz(-2.3332692) q[0];
rz(-1.7841548) q[1];
sx q[1];
rz(-2.3517377) q[1];
sx q[1];
rz(-2.435991) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7019136) q[0];
sx q[0];
rz(-1.8712303) q[0];
sx q[0];
rz(1.6070532) q[0];
x q[1];
rz(0.23106261) q[2];
sx q[2];
rz(-2.2602644) q[2];
sx q[2];
rz(-1.5514411) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.8452705) q[1];
sx q[1];
rz(-2.8883838) q[1];
sx q[1];
rz(-1.5814387) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3427395) q[3];
sx q[3];
rz(-1.0980513) q[3];
sx q[3];
rz(-0.96351876) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.6610873) q[2];
sx q[2];
rz(-1.8578119) q[2];
sx q[2];
rz(-0.60574469) q[2];
rz(-0.12039603) q[3];
sx q[3];
rz(-1.8431289) q[3];
sx q[3];
rz(0.81407434) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26009387) q[0];
sx q[0];
rz(-2.3277178) q[0];
sx q[0];
rz(-0.0023728097) q[0];
rz(0.35898769) q[1];
sx q[1];
rz(-1.9406043) q[1];
sx q[1];
rz(-0.40828362) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7255032) q[0];
sx q[0];
rz(-1.444284) q[0];
sx q[0];
rz(-0.036065412) q[0];
x q[1];
rz(-1.194656) q[2];
sx q[2];
rz(-2.500211) q[2];
sx q[2];
rz(-0.94571404) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.5376512) q[1];
sx q[1];
rz(-1.3408325) q[1];
sx q[1];
rz(1.5020788) q[1];
x q[2];
rz(-2.5465132) q[3];
sx q[3];
rz(-2.236955) q[3];
sx q[3];
rz(-0.72197589) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.6458873) q[2];
sx q[2];
rz(-0.79621035) q[2];
sx q[2];
rz(-3.0118946) q[2];
rz(2.7221223) q[3];
sx q[3];
rz(-1.9286112) q[3];
sx q[3];
rz(0.81370846) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7523338) q[0];
sx q[0];
rz(-0.77521721) q[0];
sx q[0];
rz(-2.4580521) q[0];
rz(-0.93841249) q[1];
sx q[1];
rz(-1.7529731) q[1];
sx q[1];
rz(-1.2260812) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48213595) q[0];
sx q[0];
rz(-0.99471751) q[0];
sx q[0];
rz(1.129219) q[0];
rz(-2.5415594) q[2];
sx q[2];
rz(-1.4156315) q[2];
sx q[2];
rz(0.70916353) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.9270806) q[1];
sx q[1];
rz(-1.1225554) q[1];
sx q[1];
rz(2.8238676) q[1];
rz(2.2565186) q[3];
sx q[3];
rz(-0.7739017) q[3];
sx q[3];
rz(-2.2507071) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.73416) q[2];
sx q[2];
rz(-1.4340883) q[2];
sx q[2];
rz(2.4343991) q[2];
rz(-1.4736942) q[3];
sx q[3];
rz(-2.4652822) q[3];
sx q[3];
rz(-1.7128806) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4555175) q[0];
sx q[0];
rz(-2.7954743) q[0];
sx q[0];
rz(0.92380512) q[0];
rz(2.0763981) q[1];
sx q[1];
rz(-1.1895475) q[1];
sx q[1];
rz(2.0069897) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3649639) q[0];
sx q[0];
rz(-0.37742773) q[0];
sx q[0];
rz(-1.4616184) q[0];
rz(-pi) q[1];
rz(0.85757014) q[2];
sx q[2];
rz(-0.826322) q[2];
sx q[2];
rz(0.37304953) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.90701404) q[1];
sx q[1];
rz(-1.624214) q[1];
sx q[1];
rz(2.1018886) q[1];
x q[2];
rz(-2.5328358) q[3];
sx q[3];
rz(-1.2521825) q[3];
sx q[3];
rz(-0.99005885) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.137546) q[2];
sx q[2];
rz(-1.4038439) q[2];
sx q[2];
rz(2.498632) q[2];
rz(2.9709587) q[3];
sx q[3];
rz(-2.4891487) q[3];
sx q[3];
rz(-1.0954866) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3951025) q[0];
sx q[0];
rz(-3.009142) q[0];
sx q[0];
rz(-2.3353031) q[0];
rz(3.131648) q[1];
sx q[1];
rz(-2.5237623) q[1];
sx q[1];
rz(0.26201216) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4777139) q[0];
sx q[0];
rz(-2.0933525) q[0];
sx q[0];
rz(-2.9502003) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4100894) q[2];
sx q[2];
rz(-2.2101058) q[2];
sx q[2];
rz(-2.8455545) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.9741874) q[1];
sx q[1];
rz(-1.7767118) q[1];
sx q[1];
rz(1.1307101) q[1];
x q[2];
rz(2.1266009) q[3];
sx q[3];
rz(-2.5683346) q[3];
sx q[3];
rz(-2.3044293) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.3917196) q[2];
sx q[2];
rz(-0.55730692) q[2];
sx q[2];
rz(-2.9987175) q[2];
rz(-2.2376132) q[3];
sx q[3];
rz(-1.4986135) q[3];
sx q[3];
rz(-2.2601295) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3751462) q[0];
sx q[0];
rz(-1.2297933) q[0];
sx q[0];
rz(2.4850856) q[0];
rz(-2.2625066) q[1];
sx q[1];
rz(-1.9706985) q[1];
sx q[1];
rz(0.62969977) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7555576) q[0];
sx q[0];
rz(-0.81869253) q[0];
sx q[0];
rz(-0.56171239) q[0];
x q[1];
rz(1.0378605) q[2];
sx q[2];
rz(-1.8031954) q[2];
sx q[2];
rz(0.17017066) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.9835397) q[1];
sx q[1];
rz(-3.088542) q[1];
sx q[1];
rz(-2.1237462) q[1];
rz(2.0721758) q[3];
sx q[3];
rz(-2.5540941) q[3];
sx q[3];
rz(2.0794433) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.5188344) q[2];
sx q[2];
rz(-0.82509416) q[2];
sx q[2];
rz(-0.68010124) q[2];
rz(0.71637362) q[3];
sx q[3];
rz(-0.10321897) q[3];
sx q[3];
rz(-2.0228001) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5730561) q[0];
sx q[0];
rz(-1.2328883) q[0];
sx q[0];
rz(0.20986025) q[0];
rz(-0.14280351) q[1];
sx q[1];
rz(-2.0695984) q[1];
sx q[1];
rz(-2.4048068) q[1];
rz(-0.057351107) q[2];
sx q[2];
rz(-1.3290559) q[2];
sx q[2];
rz(0.17964687) q[2];
rz(0.96137267) q[3];
sx q[3];
rz(-1.537804) q[3];
sx q[3];
rz(0.15450879) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
