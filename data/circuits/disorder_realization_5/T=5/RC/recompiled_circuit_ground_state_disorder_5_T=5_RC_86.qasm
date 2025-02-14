OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.65052819) q[0];
sx q[0];
rz(2.129038) q[0];
sx q[0];
rz(11.644047) q[0];
rz(-0.84696472) q[1];
sx q[1];
rz(4.6159336) q[1];
sx q[1];
rz(9.6428975) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.51061741) q[0];
sx q[0];
rz(-1.6625064) q[0];
sx q[0];
rz(0.46268483) q[0];
x q[1];
rz(-3.0384205) q[2];
sx q[2];
rz(-0.37984797) q[2];
sx q[2];
rz(1.0641629) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.0099758) q[1];
sx q[1];
rz(-2.9180315) q[1];
sx q[1];
rz(-1.3790591) q[1];
rz(0.078212528) q[3];
sx q[3];
rz(-1.1421176) q[3];
sx q[3];
rz(-3.0826867) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(3.0685136) q[2];
sx q[2];
rz(-1.928669) q[2];
sx q[2];
rz(2.6050513) q[2];
rz(2.1327298) q[3];
sx q[3];
rz(-2.3705685) q[3];
sx q[3];
rz(-2.3276276) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.63758481) q[0];
sx q[0];
rz(-1.8300087) q[0];
sx q[0];
rz(-3.1245533) q[0];
rz(-2.0783966) q[1];
sx q[1];
rz(-0.76779643) q[1];
sx q[1];
rz(-0.95471901) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3957152) q[0];
sx q[0];
rz(-1.3416738) q[0];
sx q[0];
rz(-1.627501) q[0];
x q[1];
rz(1.2715019) q[2];
sx q[2];
rz(-1.4134644) q[2];
sx q[2];
rz(2.2512539) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.3054637) q[1];
sx q[1];
rz(-1.3333798) q[1];
sx q[1];
rz(-1.8029455) q[1];
rz(3.122284) q[3];
sx q[3];
rz(-1.5175765) q[3];
sx q[3];
rz(-1.2217219) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.9118328) q[2];
sx q[2];
rz(-2.2288897) q[2];
sx q[2];
rz(1.1141874) q[2];
rz(-1.1344502) q[3];
sx q[3];
rz(-2.9454234) q[3];
sx q[3];
rz(2.9451059) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
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
rz(1.5512307) q[0];
sx q[0];
rz(-0.95142618) q[0];
sx q[0];
rz(2.9587342) q[0];
rz(-3.0602449) q[1];
sx q[1];
rz(-0.75129879) q[1];
sx q[1];
rz(-0.61838165) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8326022) q[0];
sx q[0];
rz(-1.497735) q[0];
sx q[0];
rz(2.354524) q[0];
rz(-2.9419627) q[2];
sx q[2];
rz(-2.3702894) q[2];
sx q[2];
rz(2.9735801) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.74677361) q[1];
sx q[1];
rz(-1.9095943) q[1];
sx q[1];
rz(-0.14751409) q[1];
rz(-0.62348714) q[3];
sx q[3];
rz(-2.1716431) q[3];
sx q[3];
rz(-2.2740325) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.1387834) q[2];
sx q[2];
rz(-1.323779) q[2];
sx q[2];
rz(-0.36661026) q[2];
rz(-1.2159411) q[3];
sx q[3];
rz(-1.1953019) q[3];
sx q[3];
rz(-2.478157) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4902896) q[0];
sx q[0];
rz(-1.9090575) q[0];
sx q[0];
rz(-0.7269727) q[0];
rz(-0.90826774) q[1];
sx q[1];
rz(-0.5779225) q[1];
sx q[1];
rz(1.04331) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3684961) q[0];
sx q[0];
rz(-2.2319622) q[0];
sx q[0];
rz(-1.5776683) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.415148) q[2];
sx q[2];
rz(-0.37964941) q[2];
sx q[2];
rz(-1.8774892) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.88731447) q[1];
sx q[1];
rz(-1.320854) q[1];
sx q[1];
rz(1.1123841) q[1];
x q[2];
rz(-2.1328385) q[3];
sx q[3];
rz(-1.2444454) q[3];
sx q[3];
rz(-2.1570753) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.72994453) q[2];
sx q[2];
rz(-0.32117716) q[2];
sx q[2];
rz(1.8357065) q[2];
rz(-3.0236687) q[3];
sx q[3];
rz(-2.6730461) q[3];
sx q[3];
rz(1.0606631) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0333198) q[0];
sx q[0];
rz(-2.1554027) q[0];
sx q[0];
rz(1.5340075) q[0];
rz(2.8458505) q[1];
sx q[1];
rz(-1.5402126) q[1];
sx q[1];
rz(-1.6827513) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0223607) q[0];
sx q[0];
rz(-1.6408477) q[0];
sx q[0];
rz(-2.3168867) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8334728) q[2];
sx q[2];
rz(-1.9240555) q[2];
sx q[2];
rz(-2.0120914) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.0572059) q[1];
sx q[1];
rz(-2.4090892) q[1];
sx q[1];
rz(1.8671579) q[1];
rz(-pi) q[2];
rz(-0.23083516) q[3];
sx q[3];
rz(-2.1394205) q[3];
sx q[3];
rz(-0.057536803) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.7532588) q[2];
sx q[2];
rz(-2.7497141) q[2];
sx q[2];
rz(0.64588109) q[2];
rz(-1.9258457) q[3];
sx q[3];
rz(-1.682621) q[3];
sx q[3];
rz(1.7997883) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.011768613) q[0];
sx q[0];
rz(-2.3079066) q[0];
sx q[0];
rz(-2.5339793) q[0];
rz(1.1786849) q[1];
sx q[1];
rz(-1.8557529) q[1];
sx q[1];
rz(2.7778621) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0640489) q[0];
sx q[0];
rz(-2.0850967) q[0];
sx q[0];
rz(0.27405996) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9035788) q[2];
sx q[2];
rz(-1.6146891) q[2];
sx q[2];
rz(-1.3258758) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.0180925) q[1];
sx q[1];
rz(-2.315633) q[1];
sx q[1];
rz(2.9781746) q[1];
rz(-pi) q[2];
x q[2];
rz(2.167114) q[3];
sx q[3];
rz(-1.5502366) q[3];
sx q[3];
rz(1.009915) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.58934775) q[2];
sx q[2];
rz(-3.037368) q[2];
sx q[2];
rz(1.1927401) q[2];
rz(-0.73181152) q[3];
sx q[3];
rz(-1.7164427) q[3];
sx q[3];
rz(-1.4881136) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8530387) q[0];
sx q[0];
rz(-2.9712501) q[0];
sx q[0];
rz(1.5227675) q[0];
rz(1.6743926) q[1];
sx q[1];
rz(-1.8312981) q[1];
sx q[1];
rz(0.67749643) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2427432) q[0];
sx q[0];
rz(-1.6209319) q[0];
sx q[0];
rz(0.87445796) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4676553) q[2];
sx q[2];
rz(-0.29602414) q[2];
sx q[2];
rz(0.37002555) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.65645331) q[1];
sx q[1];
rz(-0.59018007) q[1];
sx q[1];
rz(0.82832576) q[1];
x q[2];
rz(-0.73245184) q[3];
sx q[3];
rz(-1.6897413) q[3];
sx q[3];
rz(0.67057395) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.39685321) q[2];
sx q[2];
rz(-0.91529673) q[2];
sx q[2];
rz(-1.8358561) q[2];
rz(2.8030677) q[3];
sx q[3];
rz(-2.0207696) q[3];
sx q[3];
rz(-0.6849851) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6029538) q[0];
sx q[0];
rz(-0.21553497) q[0];
sx q[0];
rz(2.9576874) q[0];
rz(1.6869847) q[1];
sx q[1];
rz(-2.589476) q[1];
sx q[1];
rz(-2.6457381) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3199098) q[0];
sx q[0];
rz(-1.3047072) q[0];
sx q[0];
rz(-2.0771107) q[0];
x q[1];
rz(0.77347254) q[2];
sx q[2];
rz(-0.91433734) q[2];
sx q[2];
rz(-0.17420775) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.715818) q[1];
sx q[1];
rz(-1.1905021) q[1];
sx q[1];
rz(2.7575995) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0227376) q[3];
sx q[3];
rz(-1.0215852) q[3];
sx q[3];
rz(2.1191747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.054606525) q[2];
sx q[2];
rz(-2.2973674) q[2];
sx q[2];
rz(-1.1620713) q[2];
rz(1.6880796) q[3];
sx q[3];
rz(-2.1789357) q[3];
sx q[3];
rz(1.2801722) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7206955) q[0];
sx q[0];
rz(-1.9891885) q[0];
sx q[0];
rz(-2.5757117) q[0];
rz(-0.3903009) q[1];
sx q[1];
rz(-1.5719599) q[1];
sx q[1];
rz(-1.8338667) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2824664) q[0];
sx q[0];
rz(-2.0604366) q[0];
sx q[0];
rz(0.35207502) q[0];
rz(1.0729703) q[2];
sx q[2];
rz(-2.7171869) q[2];
sx q[2];
rz(-2.5967732) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.28748419) q[1];
sx q[1];
rz(-2.4696147) q[1];
sx q[1];
rz(0.42852105) q[1];
rz(-pi) q[2];
rz(-2.2447137) q[3];
sx q[3];
rz(-0.19052902) q[3];
sx q[3];
rz(-0.49303699) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.2719443) q[2];
sx q[2];
rz(-0.41676909) q[2];
sx q[2];
rz(1.0375674) q[2];
rz(-1.8218254) q[3];
sx q[3];
rz(-0.82873738) q[3];
sx q[3];
rz(0.64798361) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8808402) q[0];
sx q[0];
rz(-2.1639731) q[0];
sx q[0];
rz(0.64224893) q[0];
rz(-1.3308659) q[1];
sx q[1];
rz(-2.2763177) q[1];
sx q[1];
rz(2.9790402) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0560682) q[0];
sx q[0];
rz(-2.3352288) q[0];
sx q[0];
rz(1.0170522) q[0];
x q[1];
rz(-2.3847975) q[2];
sx q[2];
rz(-0.50163236) q[2];
sx q[2];
rz(-1.0197786) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.4288669) q[1];
sx q[1];
rz(-1.3955294) q[1];
sx q[1];
rz(2.0574942) q[1];
rz(-3.0216072) q[3];
sx q[3];
rz(-2.5549742) q[3];
sx q[3];
rz(-2.2722759) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.6315397) q[2];
sx q[2];
rz(-1.2756462) q[2];
sx q[2];
rz(-2.1007382) q[2];
rz(-0.96796525) q[3];
sx q[3];
rz(-2.0716045) q[3];
sx q[3];
rz(-0.66106838) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.80617245) q[0];
sx q[0];
rz(-1.5652884) q[0];
sx q[0];
rz(-0.096927222) q[0];
rz(-0.23282911) q[1];
sx q[1];
rz(-2.1131344) q[1];
sx q[1];
rz(0.0075385787) q[1];
rz(0.74749584) q[2];
sx q[2];
rz(-1.1807673) q[2];
sx q[2];
rz(1.8058106) q[2];
rz(0.44648689) q[3];
sx q[3];
rz(-2.3904908) q[3];
sx q[3];
rz(-0.6294546) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
