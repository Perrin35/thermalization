OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.7251627) q[0];
sx q[0];
rz(-3.0017612) q[0];
sx q[0];
rz(-0.60959417) q[0];
rz(-5.6145515) q[1];
sx q[1];
rz(0.86548391) q[1];
sx q[1];
rz(15.794985) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5370109) q[0];
sx q[0];
rz(-2.4049979) q[0];
sx q[0];
rz(1.0010615) q[0];
x q[1];
rz(-3.0084228) q[2];
sx q[2];
rz(-2.3699017) q[2];
sx q[2];
rz(-2.0193677) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.7510208) q[1];
sx q[1];
rz(-2.6387072) q[1];
sx q[1];
rz(-1.8608171) q[1];
rz(1.617241) q[3];
sx q[3];
rz(-1.9825476) q[3];
sx q[3];
rz(-2.9573033) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.13880754) q[2];
sx q[2];
rz(-1.3802718) q[2];
sx q[2];
rz(0.37386093) q[2];
rz(-2.8047681) q[3];
sx q[3];
rz(-1.5461494) q[3];
sx q[3];
rz(-0.22836223) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8841298) q[0];
sx q[0];
rz(-2.7682436) q[0];
sx q[0];
rz(-1.9467547) q[0];
rz(0.082611235) q[1];
sx q[1];
rz(-1.1673085) q[1];
sx q[1];
rz(-0.00037489051) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9276792) q[0];
sx q[0];
rz(-1.9880901) q[0];
sx q[0];
rz(2.7620865) q[0];
x q[1];
rz(2.9412494) q[2];
sx q[2];
rz(-1.3140972) q[2];
sx q[2];
rz(-2.8368907) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-3.0034069) q[1];
sx q[1];
rz(-2.1669743) q[1];
sx q[1];
rz(-1.2852933) q[1];
rz(1.625929) q[3];
sx q[3];
rz(-1.1109567) q[3];
sx q[3];
rz(0.22051792) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.3327545) q[2];
sx q[2];
rz(-2.9465582) q[2];
sx q[2];
rz(2.9648798) q[2];
rz(-0.79408944) q[3];
sx q[3];
rz(-2.3637171) q[3];
sx q[3];
rz(-0.55364048) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90213838) q[0];
sx q[0];
rz(-2.4376526) q[0];
sx q[0];
rz(-2.7368271) q[0];
rz(1.2813214) q[1];
sx q[1];
rz(-0.57360137) q[1];
sx q[1];
rz(-1.8331029) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.20298) q[0];
sx q[0];
rz(-0.70883195) q[0];
sx q[0];
rz(2.3908486) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8070418) q[2];
sx q[2];
rz(-2.0958732) q[2];
sx q[2];
rz(0.46307785) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.1723459) q[1];
sx q[1];
rz(-1.8241276) q[1];
sx q[1];
rz(-2.2971056) q[1];
rz(-0.19034068) q[3];
sx q[3];
rz(-2.6958709) q[3];
sx q[3];
rz(1.0090855) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.219316) q[2];
sx q[2];
rz(-0.52307659) q[2];
sx q[2];
rz(3.1211839) q[2];
rz(-1.071788) q[3];
sx q[3];
rz(-2.0139549) q[3];
sx q[3];
rz(2.4782457) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9999009) q[0];
sx q[0];
rz(-2.2321556) q[0];
sx q[0];
rz(-2.2256057) q[0];
rz(2.6782716) q[1];
sx q[1];
rz(-1.0659734) q[1];
sx q[1];
rz(2.0844918) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.076544882) q[0];
sx q[0];
rz(-2.0079945) q[0];
sx q[0];
rz(-1.5989499) q[0];
rz(-pi) q[1];
x q[1];
rz(0.82614233) q[2];
sx q[2];
rz(-0.75781265) q[2];
sx q[2];
rz(3.0405424) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.4790198) q[1];
sx q[1];
rz(-2.9111007) q[1];
sx q[1];
rz(-1.6307994) q[1];
x q[2];
rz(-2.9017157) q[3];
sx q[3];
rz(-0.7719709) q[3];
sx q[3];
rz(-3.1023657) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.2477734) q[2];
sx q[2];
rz(-0.63988581) q[2];
sx q[2];
rz(-2.7988953) q[2];
rz(-1.6977067) q[3];
sx q[3];
rz(-2.2659437) q[3];
sx q[3];
rz(0.7152344) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3559568) q[0];
sx q[0];
rz(-2.5781093) q[0];
sx q[0];
rz(3.0551531) q[0];
rz(1.3899639) q[1];
sx q[1];
rz(-0.93170634) q[1];
sx q[1];
rz(2.6729029) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8353032) q[0];
sx q[0];
rz(-1.9625184) q[0];
sx q[0];
rz(2.3098371) q[0];
rz(-pi) q[1];
rz(1.3895949) q[2];
sx q[2];
rz(-1.2672408) q[2];
sx q[2];
rz(-0.62523491) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.17051711) q[1];
sx q[1];
rz(-1.7726388) q[1];
sx q[1];
rz(1.438373) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.080805578) q[3];
sx q[3];
rz(-0.81597933) q[3];
sx q[3];
rz(-0.24085837) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.1420574) q[2];
sx q[2];
rz(-2.2613328) q[2];
sx q[2];
rz(2.2772677) q[2];
rz(2.9243829) q[3];
sx q[3];
rz(-2.5253798) q[3];
sx q[3];
rz(2.8488081) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7140759) q[0];
sx q[0];
rz(-1.7345411) q[0];
sx q[0];
rz(-0.19700225) q[0];
rz(1.3621832) q[1];
sx q[1];
rz(-1.4809337) q[1];
sx q[1];
rz(-2.8053455) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7173548) q[0];
sx q[0];
rz(-0.80168085) q[0];
sx q[0];
rz(0.63071155) q[0];
rz(-2.5480812) q[2];
sx q[2];
rz(-1.0700018) q[2];
sx q[2];
rz(-1.3154495) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.044144883) q[1];
sx q[1];
rz(-2.1598585) q[1];
sx q[1];
rz(1.5138022) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7429966) q[3];
sx q[3];
rz(-2.2599054) q[3];
sx q[3];
rz(0.25758753) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.53529915) q[2];
sx q[2];
rz(-1.665411) q[2];
sx q[2];
rz(0.84632787) q[2];
rz(-1.9019295) q[3];
sx q[3];
rz(-1.2318434) q[3];
sx q[3];
rz(0.0011750778) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7476615) q[0];
sx q[0];
rz(-2.1037536) q[0];
sx q[0];
rz(-2.8826662) q[0];
rz(-1.3461643) q[1];
sx q[1];
rz(-1.3849473) q[1];
sx q[1];
rz(-1.0940201) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8261995) q[0];
sx q[0];
rz(-2.0059735) q[0];
sx q[0];
rz(0.3776334) q[0];
rz(-3.1292186) q[2];
sx q[2];
rz(-1.3482598) q[2];
sx q[2];
rz(2.6518133) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.104407) q[1];
sx q[1];
rz(-1.6652602) q[1];
sx q[1];
rz(3.1135998) q[1];
rz(-1.592698) q[3];
sx q[3];
rz(-1.4727482) q[3];
sx q[3];
rz(0.4050307) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.1743494) q[2];
sx q[2];
rz(-1.3568342) q[2];
sx q[2];
rz(2.8209177) q[2];
rz(-0.43618068) q[3];
sx q[3];
rz(-2.6685721) q[3];
sx q[3];
rz(2.6935553) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2748579) q[0];
sx q[0];
rz(-0.47645706) q[0];
sx q[0];
rz(1.0048237) q[0];
rz(-0.59016219) q[1];
sx q[1];
rz(-2.2361123) q[1];
sx q[1];
rz(2.337713) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.77792149) q[0];
sx q[0];
rz(-1.6102618) q[0];
sx q[0];
rz(-0.58593336) q[0];
rz(-2.5875823) q[2];
sx q[2];
rz(-1.4129352) q[2];
sx q[2];
rz(-2.6487034) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.56402928) q[1];
sx q[1];
rz(-1.659435) q[1];
sx q[1];
rz(-3.1313529) q[1];
rz(-pi) q[2];
rz(0.6635267) q[3];
sx q[3];
rz(-1.6420206) q[3];
sx q[3];
rz(-2.0957259) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.8089495) q[2];
sx q[2];
rz(-2.2968764) q[2];
sx q[2];
rz(2.1311029) q[2];
rz(-0.60339749) q[3];
sx q[3];
rz(-1.6167275) q[3];
sx q[3];
rz(2.5914014) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5378961) q[0];
sx q[0];
rz(-1.864707) q[0];
sx q[0];
rz(0.4075152) q[0];
rz(0.28911668) q[1];
sx q[1];
rz(-1.1228077) q[1];
sx q[1];
rz(-0.75072748) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4780316) q[0];
sx q[0];
rz(-2.7285517) q[0];
sx q[0];
rz(2.3703299) q[0];
x q[1];
rz(-0.27955555) q[2];
sx q[2];
rz(-1.4923555) q[2];
sx q[2];
rz(-2.560175) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.1474485) q[1];
sx q[1];
rz(-1.7732883) q[1];
sx q[1];
rz(1.8840428) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5969855) q[3];
sx q[3];
rz(-1.2695754) q[3];
sx q[3];
rz(-1.6645886) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.0687381) q[2];
sx q[2];
rz(-2.4170503) q[2];
sx q[2];
rz(-1.9753974) q[2];
rz(1.3646305) q[3];
sx q[3];
rz(-2.3659673) q[3];
sx q[3];
rz(3.1197746) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5831379) q[0];
sx q[0];
rz(-0.82644176) q[0];
sx q[0];
rz(-1.7512084) q[0];
rz(0.33069262) q[1];
sx q[1];
rz(-2.3762517) q[1];
sx q[1];
rz(1.4601382) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4936192) q[0];
sx q[0];
rz(-0.68035347) q[0];
sx q[0];
rz(-0.86662678) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5286469) q[2];
sx q[2];
rz(-0.50439207) q[2];
sx q[2];
rz(1.991589) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.9863661) q[1];
sx q[1];
rz(-0.44253293) q[1];
sx q[1];
rz(2.7000484) q[1];
rz(-pi) q[2];
x q[2];
rz(0.7695997) q[3];
sx q[3];
rz(-2.1160612) q[3];
sx q[3];
rz(-1.256497) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.34974393) q[2];
sx q[2];
rz(-1.3246374) q[2];
sx q[2];
rz(-1.6798518) q[2];
rz(-1.1200303) q[3];
sx q[3];
rz(-2.5199065) q[3];
sx q[3];
rz(-2.6859443) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.6256975) q[0];
sx q[0];
rz(-1.5213756) q[0];
sx q[0];
rz(-1.9949927) q[0];
rz(-1.760578) q[1];
sx q[1];
rz(-1.2881423) q[1];
sx q[1];
rz(-1.2013411) q[1];
rz(2.7881651) q[2];
sx q[2];
rz(-0.8197439) q[2];
sx q[2];
rz(2.843429) q[2];
rz(2.667726) q[3];
sx q[3];
rz(-0.73054536) q[3];
sx q[3];
rz(1.2252145) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];