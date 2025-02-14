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
rz(0.50245589) q[0];
sx q[0];
rz(-2.034019) q[0];
sx q[0];
rz(-1.7662319) q[0];
rz(-3.4604685) q[1];
sx q[1];
rz(4.7825216) q[1];
sx q[1];
rz(8.6934269) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.80901805) q[0];
sx q[0];
rz(-1.8900196) q[0];
sx q[0];
rz(0.92895023) q[0];
rz(-pi) q[1];
rz(3.1362304) q[2];
sx q[2];
rz(-0.66197936) q[2];
sx q[2];
rz(2.299451) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.3870942) q[1];
sx q[1];
rz(-1.525314) q[1];
sx q[1];
rz(0.78630739) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9792487) q[3];
sx q[3];
rz(-1.663891) q[3];
sx q[3];
rz(-2.7422991) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.88087624) q[2];
sx q[2];
rz(-1.9274351) q[2];
sx q[2];
rz(-0.94240776) q[2];
rz(2.2830394) q[3];
sx q[3];
rz(-0.61045727) q[3];
sx q[3];
rz(-0.71949351) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.061676625) q[0];
sx q[0];
rz(-2.583857) q[0];
sx q[0];
rz(3.0829561) q[0];
rz(0.05642852) q[1];
sx q[1];
rz(-1.7516878) q[1];
sx q[1];
rz(-1.1172392) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3836544) q[0];
sx q[0];
rz(-2.2636739) q[0];
sx q[0];
rz(-0.52097676) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1740721) q[2];
sx q[2];
rz(-0.32336059) q[2];
sx q[2];
rz(-2.4692995) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.5176425) q[1];
sx q[1];
rz(-0.25282598) q[1];
sx q[1];
rz(-1.9060017) q[1];
x q[2];
rz(-0.23157236) q[3];
sx q[3];
rz(-0.48677126) q[3];
sx q[3];
rz(2.6269905) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.3650018) q[2];
sx q[2];
rz(-0.33010179) q[2];
sx q[2];
rz(-0.4064202) q[2];
rz(-2.8255919) q[3];
sx q[3];
rz(-1.4603115) q[3];
sx q[3];
rz(2.3587904) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86421788) q[0];
sx q[0];
rz(-2.6835231) q[0];
sx q[0];
rz(1.8055441) q[0];
rz(2.863073) q[1];
sx q[1];
rz(-1.7428935) q[1];
sx q[1];
rz(-0.96885243) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8271885) q[0];
sx q[0];
rz(-2.1705856) q[0];
sx q[0];
rz(1.1161854) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1325705) q[2];
sx q[2];
rz(-2.2780905) q[2];
sx q[2];
rz(-2.2826441) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.90889764) q[1];
sx q[1];
rz(-1.8669584) q[1];
sx q[1];
rz(-3.0975268) q[1];
rz(2.4163867) q[3];
sx q[3];
rz(-1.375633) q[3];
sx q[3];
rz(-0.31828398) q[3];
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
rz(1.3052321) q[3];
sx q[3];
rz(-0.94279083) q[3];
sx q[3];
rz(1.85873) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99821943) q[0];
sx q[0];
rz(-2.0046736) q[0];
sx q[0];
rz(2.112222) q[0];
rz(-2.3221305) q[1];
sx q[1];
rz(-1.9792604) q[1];
sx q[1];
rz(-0.79684657) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5982847) q[0];
sx q[0];
rz(-2.6572022) q[0];
sx q[0];
rz(-0.422307) q[0];
rz(-pi) q[1];
rz(-3.1355924) q[2];
sx q[2];
rz(-2.0987392) q[2];
sx q[2];
rz(2.1354298) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.0734892) q[1];
sx q[1];
rz(-1.7760381) q[1];
sx q[1];
rz(-1.8480592) q[1];
rz(-pi) q[2];
x q[2];
rz(0.66241499) q[3];
sx q[3];
rz(-1.9180505) q[3];
sx q[3];
rz(1.9062689) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-3.0113819) q[2];
sx q[2];
rz(-0.036616651) q[2];
sx q[2];
rz(1.9127362) q[2];
rz(-2.7025488) q[3];
sx q[3];
rz(-1.565758) q[3];
sx q[3];
rz(1.2307833) q[3];
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
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48290408) q[0];
sx q[0];
rz(-1.6084325) q[0];
sx q[0];
rz(2.3332692) q[0];
rz(1.3574379) q[1];
sx q[1];
rz(-2.3517377) q[1];
sx q[1];
rz(0.70560169) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3177174) q[0];
sx q[0];
rz(-2.8390445) q[0];
sx q[0];
rz(-3.0251193) q[0];
x q[1];
rz(2.91053) q[2];
sx q[2];
rz(-2.2602644) q[2];
sx q[2];
rz(1.5514411) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.8452705) q[1];
sx q[1];
rz(-0.25320881) q[1];
sx q[1];
rz(1.5601539) q[1];
rz(0.6198778) q[3];
sx q[3];
rz(-2.2408463) q[3];
sx q[3];
rz(2.1170962) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.48050532) q[2];
sx q[2];
rz(-1.8578119) q[2];
sx q[2];
rz(-2.535848) q[2];
rz(0.12039603) q[3];
sx q[3];
rz(-1.8431289) q[3];
sx q[3];
rz(-0.81407434) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26009387) q[0];
sx q[0];
rz(-0.81387481) q[0];
sx q[0];
rz(0.0023728097) q[0];
rz(2.782605) q[1];
sx q[1];
rz(-1.9406043) q[1];
sx q[1];
rz(-2.733309) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6946163) q[0];
sx q[0];
rz(-0.13152619) q[0];
sx q[0];
rz(-1.294554) q[0];
rz(-pi) q[1];
rz(2.177816) q[2];
sx q[2];
rz(-1.7923819) q[2];
sx q[2];
rz(0.9315679) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.6039415) q[1];
sx q[1];
rz(-1.3408325) q[1];
sx q[1];
rz(-1.6395139) q[1];
rz(-pi) q[2];
rz(0.8114641) q[3];
sx q[3];
rz(-2.0272019) q[3];
sx q[3];
rz(-1.8965669) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.4957054) q[2];
sx q[2];
rz(-0.79621035) q[2];
sx q[2];
rz(-0.1296981) q[2];
rz(2.7221223) q[3];
sx q[3];
rz(-1.2129815) q[3];
sx q[3];
rz(-0.81370846) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7523338) q[0];
sx q[0];
rz(-2.3663754) q[0];
sx q[0];
rz(2.4580521) q[0];
rz(2.2031802) q[1];
sx q[1];
rz(-1.3886195) q[1];
sx q[1];
rz(-1.9155115) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9447362) q[0];
sx q[0];
rz(-0.71030197) q[0];
sx q[0];
rz(-2.5596749) q[0];
x q[1];
rz(-1.7581045) q[2];
sx q[2];
rz(-0.97895998) q[2];
sx q[2];
rz(-0.96697742) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.56537762) q[1];
sx q[1];
rz(-2.5984796) q[1];
sx q[1];
rz(0.99467038) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5874373) q[3];
sx q[3];
rz(-2.1423577) q[3];
sx q[3];
rz(-3.1031648) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.40743264) q[2];
sx q[2];
rz(-1.4340883) q[2];
sx q[2];
rz(-2.4343991) q[2];
rz(1.6678984) q[3];
sx q[3];
rz(-0.67631045) q[3];
sx q[3];
rz(1.7128806) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4555175) q[0];
sx q[0];
rz(-0.34611836) q[0];
sx q[0];
rz(-0.92380512) q[0];
rz(-2.0763981) q[1];
sx q[1];
rz(-1.1895475) q[1];
sx q[1];
rz(-2.0069897) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3649639) q[0];
sx q[0];
rz(-0.37742773) q[0];
sx q[0];
rz(-1.4616184) q[0];
rz(-2.2840225) q[2];
sx q[2];
rz(-2.3152707) q[2];
sx q[2];
rz(2.7685431) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.90701404) q[1];
sx q[1];
rz(-1.624214) q[1];
sx q[1];
rz(-1.0397041) q[1];
rz(1.9530961) q[3];
sx q[3];
rz(-2.144882) q[3];
sx q[3];
rz(-2.7758383) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.137546) q[2];
sx q[2];
rz(-1.7377487) q[2];
sx q[2];
rz(-2.498632) q[2];
rz(2.9709587) q[3];
sx q[3];
rz(-2.4891487) q[3];
sx q[3];
rz(-1.0954866) q[3];
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
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.74649015) q[0];
sx q[0];
rz(-0.13245067) q[0];
sx q[0];
rz(-2.3353031) q[0];
rz(-0.0099446615) q[1];
sx q[1];
rz(-2.5237623) q[1];
sx q[1];
rz(-2.8795805) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1073955) q[0];
sx q[0];
rz(-2.5881564) q[0];
sx q[0];
rz(1.2518) q[0];
x q[1];
rz(0.21199439) q[2];
sx q[2];
rz(-0.65644453) q[2];
sx q[2];
rz(0.5613297) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.6422185) q[1];
sx q[1];
rz(-1.1406349) q[1];
sx q[1];
rz(-0.22689928) q[1];
rz(1.069182) q[3];
sx q[3];
rz(-1.2805689) q[3];
sx q[3];
rz(-0.25267664) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.749873) q[2];
sx q[2];
rz(-0.55730692) q[2];
sx q[2];
rz(2.9987175) q[2];
rz(0.90397942) q[3];
sx q[3];
rz(-1.4986135) q[3];
sx q[3];
rz(-2.2601295) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7664465) q[0];
sx q[0];
rz(-1.2297933) q[0];
sx q[0];
rz(-2.4850856) q[0];
rz(2.2625066) q[1];
sx q[1];
rz(-1.1708941) q[1];
sx q[1];
rz(-2.5118929) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3860351) q[0];
sx q[0];
rz(-2.3229001) q[0];
sx q[0];
rz(-0.56171239) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1348502) q[2];
sx q[2];
rz(-0.576888) q[2];
sx q[2];
rz(2.1132121) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.1765285) q[1];
sx q[1];
rz(-1.5986491) q[1];
sx q[1];
rz(-1.615953) q[1];
rz(-pi) q[2];
rz(1.042243) q[3];
sx q[3];
rz(-1.8404598) q[3];
sx q[3];
rz(2.2049513) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.62275824) q[2];
sx q[2];
rz(-0.82509416) q[2];
sx q[2];
rz(0.68010124) q[2];
rz(2.425219) q[3];
sx q[3];
rz(-3.0383737) q[3];
sx q[3];
rz(-2.0228001) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
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
rz(1.5730561) q[0];
sx q[0];
rz(-1.9087044) q[0];
sx q[0];
rz(-2.9317324) q[0];
rz(-0.14280351) q[1];
sx q[1];
rz(-2.0695984) q[1];
sx q[1];
rz(-2.4048068) q[1];
rz(3.0842415) q[2];
sx q[2];
rz(-1.3290559) q[2];
sx q[2];
rz(0.17964687) q[2];
rz(-0.96137267) q[3];
sx q[3];
rz(-1.6037887) q[3];
sx q[3];
rz(-2.9870839) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
