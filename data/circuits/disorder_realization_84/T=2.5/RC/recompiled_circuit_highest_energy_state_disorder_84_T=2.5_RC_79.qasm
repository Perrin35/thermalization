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
rz(-0.62491971) q[0];
sx q[0];
rz(-1.3345557) q[0];
sx q[0];
rz(-0.053243756) q[0];
rz(-2.6553395) q[1];
sx q[1];
rz(-3.0142205) q[1];
sx q[1];
rz(-1.706634) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7409748) q[0];
sx q[0];
rz(-1.8549998) q[0];
sx q[0];
rz(1.4301197) q[0];
x q[1];
rz(-1.6900914) q[2];
sx q[2];
rz(-1.8362459) q[2];
sx q[2];
rz(-0.99670289) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.3327121) q[1];
sx q[1];
rz(-1.1774027) q[1];
sx q[1];
rz(1.5129526) q[1];
rz(-pi) q[2];
x q[2];
rz(0.94287431) q[3];
sx q[3];
rz(-0.74340313) q[3];
sx q[3];
rz(-1.2974031) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.80483738) q[2];
sx q[2];
rz(-1.3433604) q[2];
sx q[2];
rz(2.8678144) q[2];
rz(-1.2182419) q[3];
sx q[3];
rz(-2.4680586) q[3];
sx q[3];
rz(-0.28606733) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
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
rz(3.1186721) q[0];
sx q[0];
rz(-2.3781222) q[0];
sx q[0];
rz(-2.3619695) q[0];
rz(0.66462213) q[1];
sx q[1];
rz(-1.9326991) q[1];
sx q[1];
rz(1.2453311) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1561403) q[0];
sx q[0];
rz(-0.13053556) q[0];
sx q[0];
rz(-1.077561) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.092463569) q[2];
sx q[2];
rz(-2.2801054) q[2];
sx q[2];
rz(0.48395983) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.1511953) q[1];
sx q[1];
rz(-1.510069) q[1];
sx q[1];
rz(-1.7833738) q[1];
rz(-pi) q[2];
rz(0.20540463) q[3];
sx q[3];
rz(-2.1202135) q[3];
sx q[3];
rz(-0.37976625) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.7094884) q[2];
sx q[2];
rz(-2.4435142) q[2];
sx q[2];
rz(0.049987642) q[2];
rz(1.6677808) q[3];
sx q[3];
rz(-1.4155017) q[3];
sx q[3];
rz(-0.93562359) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1525928) q[0];
sx q[0];
rz(-0.86243668) q[0];
sx q[0];
rz(1.1784026) q[0];
rz(-1.1475457) q[1];
sx q[1];
rz(-0.18745628) q[1];
sx q[1];
rz(-0.60290927) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41571799) q[0];
sx q[0];
rz(-2.9511669) q[0];
sx q[0];
rz(-2.1384396) q[0];
rz(-pi) q[1];
rz(0.52003543) q[2];
sx q[2];
rz(-2.1923861) q[2];
sx q[2];
rz(3.1261895) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.6711839) q[1];
sx q[1];
rz(-2.3215331) q[1];
sx q[1];
rz(2.5633308) q[1];
rz(-0.52899482) q[3];
sx q[3];
rz(-0.96312614) q[3];
sx q[3];
rz(2.0207456) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.19043645) q[2];
sx q[2];
rz(-2.3556605) q[2];
sx q[2];
rz(0.3832761) q[2];
rz(1.3112618) q[3];
sx q[3];
rz(-1.0878599) q[3];
sx q[3];
rz(-0.12152984) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.7168147) q[0];
sx q[0];
rz(-0.99262339) q[0];
sx q[0];
rz(-2.2130261) q[0];
rz(0.83944744) q[1];
sx q[1];
rz(-1.3176368) q[1];
sx q[1];
rz(2.3323257) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2977375) q[0];
sx q[0];
rz(-2.1240196) q[0];
sx q[0];
rz(-2.7064302) q[0];
x q[1];
rz(-0.35688422) q[2];
sx q[2];
rz(-0.45512558) q[2];
sx q[2];
rz(-0.17562107) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.2553808) q[1];
sx q[1];
rz(-1.501923) q[1];
sx q[1];
rz(-1.8167472) q[1];
x q[2];
rz(-2.094481) q[3];
sx q[3];
rz(-2.4274656) q[3];
sx q[3];
rz(1.7585825) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.75260085) q[2];
sx q[2];
rz(-1.5316803) q[2];
sx q[2];
rz(-1.4468225) q[2];
rz(-1.1270479) q[3];
sx q[3];
rz(-2.4067252) q[3];
sx q[3];
rz(0.62895044) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9160354) q[0];
sx q[0];
rz(-1.9434403) q[0];
sx q[0];
rz(1.4300562) q[0];
rz(2.9043708) q[1];
sx q[1];
rz(-2.1784541) q[1];
sx q[1];
rz(-2.846834) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65564102) q[0];
sx q[0];
rz(-1.427934) q[0];
sx q[0];
rz(0.44012196) q[0];
rz(-pi) q[1];
rz(1.1072269) q[2];
sx q[2];
rz(-0.71706248) q[2];
sx q[2];
rz(2.0328558) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.7476866) q[1];
sx q[1];
rz(-1.9853885) q[1];
sx q[1];
rz(2.0646329) q[1];
rz(2.5676457) q[3];
sx q[3];
rz(-1.6960295) q[3];
sx q[3];
rz(-2.4698911) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.1416867) q[2];
sx q[2];
rz(-2.9408231) q[2];
sx q[2];
rz(-3.0086009) q[2];
rz(2.3717132) q[3];
sx q[3];
rz(-1.8498288) q[3];
sx q[3];
rz(2.8225115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28984508) q[0];
sx q[0];
rz(-0.3388437) q[0];
sx q[0];
rz(-1.7293365) q[0];
rz(0.051636592) q[1];
sx q[1];
rz(-2.5187571) q[1];
sx q[1];
rz(-0.19270611) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15899236) q[0];
sx q[0];
rz(-1.5518414) q[0];
sx q[0];
rz(2.6277831) q[0];
x q[1];
rz(0.95780722) q[2];
sx q[2];
rz(-0.93767208) q[2];
sx q[2];
rz(1.9547966) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.6289857) q[1];
sx q[1];
rz(-0.64388212) q[1];
sx q[1];
rz(0.75276796) q[1];
rz(-pi) q[2];
rz(-1.1901549) q[3];
sx q[3];
rz(-2.1771113) q[3];
sx q[3];
rz(-2.9553138) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.073033832) q[2];
sx q[2];
rz(-1.2898338) q[2];
sx q[2];
rz(1.3812836) q[2];
rz(2.6265788) q[3];
sx q[3];
rz(-0.88740715) q[3];
sx q[3];
rz(1.380434) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9800022) q[0];
sx q[0];
rz(-1.9965633) q[0];
sx q[0];
rz(-2.7951151) q[0];
rz(1.0193635) q[1];
sx q[1];
rz(-2.4788224) q[1];
sx q[1];
rz(1.6874708) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2812735) q[0];
sx q[0];
rz(-0.77206445) q[0];
sx q[0];
rz(-0.66769974) q[0];
rz(-pi) q[1];
rz(2.2242111) q[2];
sx q[2];
rz(-1.1251984) q[2];
sx q[2];
rz(3.0609956) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.2994308) q[1];
sx q[1];
rz(-0.60363942) q[1];
sx q[1];
rz(-2.0718715) q[1];
x q[2];
rz(0.99549528) q[3];
sx q[3];
rz(-0.18333215) q[3];
sx q[3];
rz(-0.14106476) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.5399897) q[2];
sx q[2];
rz(-1.3301962) q[2];
sx q[2];
rz(-2.9023564) q[2];
rz(2.7573977) q[3];
sx q[3];
rz(-0.68238634) q[3];
sx q[3];
rz(0.95170963) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.20188986) q[0];
sx q[0];
rz(-1.6738142) q[0];
sx q[0];
rz(1.3772759) q[0];
rz(-1.9210303) q[1];
sx q[1];
rz(-0.86253291) q[1];
sx q[1];
rz(0.16622226) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3275422) q[0];
sx q[0];
rz(-0.53314994) q[0];
sx q[0];
rz(2.0943805) q[0];
rz(1.2532089) q[2];
sx q[2];
rz(-0.86195213) q[2];
sx q[2];
rz(-2.6203757) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.6320914) q[1];
sx q[1];
rz(-2.9281514) q[1];
sx q[1];
rz(-2.9963993) q[1];
x q[2];
rz(-2.1020426) q[3];
sx q[3];
rz(-2.663718) q[3];
sx q[3];
rz(0.22051624) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.347747) q[2];
sx q[2];
rz(-0.85024992) q[2];
sx q[2];
rz(-1.0469077) q[2];
rz(1.2547803) q[3];
sx q[3];
rz(-2.051765) q[3];
sx q[3];
rz(1.8449529) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.60085249) q[0];
sx q[0];
rz(-0.37373251) q[0];
sx q[0];
rz(2.0284213) q[0];
rz(2.3241849) q[1];
sx q[1];
rz(-1.6956885) q[1];
sx q[1];
rz(-0.36849749) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2558852) q[0];
sx q[0];
rz(-1.2987381) q[0];
sx q[0];
rz(-1.9291124) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3992052) q[2];
sx q[2];
rz(-1.6642804) q[2];
sx q[2];
rz(2.7274318) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.2713743) q[1];
sx q[1];
rz(-1.8370359) q[1];
sx q[1];
rz(0.72721796) q[1];
rz(-pi) q[2];
rz(-0.59898563) q[3];
sx q[3];
rz(-2.3040207) q[3];
sx q[3];
rz(2.4323104) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.9937399) q[2];
sx q[2];
rz(-0.61823121) q[2];
sx q[2];
rz(-1.7842133) q[2];
rz(-0.74710685) q[3];
sx q[3];
rz(-0.96304572) q[3];
sx q[3];
rz(-0.67556206) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.52637446) q[0];
sx q[0];
rz(-0.45658657) q[0];
sx q[0];
rz(-1.0427465) q[0];
rz(-2.3710947) q[1];
sx q[1];
rz(-2.1310525) q[1];
sx q[1];
rz(-0.76046336) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7211518) q[0];
sx q[0];
rz(-1.8095224) q[0];
sx q[0];
rz(1.9107642) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3907688) q[2];
sx q[2];
rz(-1.5329201) q[2];
sx q[2];
rz(1.2779209) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.14899602) q[1];
sx q[1];
rz(-1.1482052) q[1];
sx q[1];
rz(-1.6572957) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6728739) q[3];
sx q[3];
rz(-0.7112452) q[3];
sx q[3];
rz(-2.6576692) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.20327917) q[2];
sx q[2];
rz(-2.5310897) q[2];
sx q[2];
rz(-0.40943134) q[2];
rz(-1.5583386) q[3];
sx q[3];
rz(-0.46054545) q[3];
sx q[3];
rz(-2.515365) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.64774491) q[0];
sx q[0];
rz(-2.0997601) q[0];
sx q[0];
rz(0.44152015) q[0];
rz(1.582675) q[1];
sx q[1];
rz(-1.1987004) q[1];
sx q[1];
rz(-0.62335062) q[1];
rz(0.10298038) q[2];
sx q[2];
rz(-1.5126577) q[2];
sx q[2];
rz(2.589227) q[2];
rz(0.24258976) q[3];
sx q[3];
rz(-1.6995097) q[3];
sx q[3];
rz(-2.6293737) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
