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
rz(0.12662521) q[0];
sx q[0];
rz(1.5746483) q[0];
sx q[0];
rz(9.9404542) q[0];
rz(0.22817837) q[1];
sx q[1];
rz(-0.74695865) q[1];
sx q[1];
rz(0.42626122) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.92862064) q[0];
sx q[0];
rz(-0.80982319) q[0];
sx q[0];
rz(0.23913236) q[0];
rz(-pi) q[1];
rz(-0.50349109) q[2];
sx q[2];
rz(-1.8241183) q[2];
sx q[2];
rz(1.5905141) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.961024) q[1];
sx q[1];
rz(-2.5371309) q[1];
sx q[1];
rz(-3.040041) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.32483806) q[3];
sx q[3];
rz(-2.1305314) q[3];
sx q[3];
rz(-1.4046275) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.8969741) q[2];
sx q[2];
rz(-1.8316869) q[2];
sx q[2];
rz(-0.2429602) q[2];
rz(-2.7729559) q[3];
sx q[3];
rz(-0.60550767) q[3];
sx q[3];
rz(1.9093556) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31545562) q[0];
sx q[0];
rz(-2.9189126) q[0];
sx q[0];
rz(-2.7742703) q[0];
rz(0.79633725) q[1];
sx q[1];
rz(-2.0834736) q[1];
sx q[1];
rz(-2.499089) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.028051826) q[0];
sx q[0];
rz(-1.100825) q[0];
sx q[0];
rz(2.0451106) q[0];
rz(-pi) q[1];
rz(0.60566492) q[2];
sx q[2];
rz(-2.4993651) q[2];
sx q[2];
rz(0.58711827) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.6141908) q[1];
sx q[1];
rz(-2.3878015) q[1];
sx q[1];
rz(0.68282737) q[1];
rz(-pi) q[2];
rz(2.0525371) q[3];
sx q[3];
rz(-1.6828487) q[3];
sx q[3];
rz(-0.58603906) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.6392886) q[2];
sx q[2];
rz(-0.93561155) q[2];
sx q[2];
rz(1.3703692) q[2];
rz(2.9141407) q[3];
sx q[3];
rz(-1.8919614) q[3];
sx q[3];
rz(-1.1915709) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3962536) q[0];
sx q[0];
rz(-0.17530137) q[0];
sx q[0];
rz(-1.1981717) q[0];
rz(-2.0612969) q[1];
sx q[1];
rz(-2.9273169) q[1];
sx q[1];
rz(0.094873039) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.085233) q[0];
sx q[0];
rz(-1.0329536) q[0];
sx q[0];
rz(2.8406124) q[0];
rz(-pi) q[1];
x q[1];
rz(0.40887654) q[2];
sx q[2];
rz(-0.95252242) q[2];
sx q[2];
rz(2.8586819) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.66902924) q[1];
sx q[1];
rz(-1.5298784) q[1];
sx q[1];
rz(-2.1075691) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9547988) q[3];
sx q[3];
rz(-1.3130762) q[3];
sx q[3];
rz(-0.4674165) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.0692811) q[2];
sx q[2];
rz(-0.60439622) q[2];
sx q[2];
rz(0.4064694) q[2];
rz(1.9629924) q[3];
sx q[3];
rz(-1.4402163) q[3];
sx q[3];
rz(-1.9788205) q[3];
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
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.56226319) q[0];
sx q[0];
rz(-3.0070906) q[0];
sx q[0];
rz(2.80559) q[0];
rz(2.621189) q[1];
sx q[1];
rz(-0.86417472) q[1];
sx q[1];
rz(-3.0373108) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.18096237) q[0];
sx q[0];
rz(-1.2237566) q[0];
sx q[0];
rz(2.7913559) q[0];
rz(-pi) q[1];
rz(-1.5566795) q[2];
sx q[2];
rz(-1.7272564) q[2];
sx q[2];
rz(-0.9597646) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.6426864) q[1];
sx q[1];
rz(-0.79267529) q[1];
sx q[1];
rz(-2.098987) q[1];
rz(1.453425) q[3];
sx q[3];
rz(-1.1201829) q[3];
sx q[3];
rz(0.71244682) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.5229554) q[2];
sx q[2];
rz(-1.0982265) q[2];
sx q[2];
rz(1.1445507) q[2];
rz(2.5942904) q[3];
sx q[3];
rz(-1.2283044) q[3];
sx q[3];
rz(1.087629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.8892141) q[0];
sx q[0];
rz(-0.19110282) q[0];
sx q[0];
rz(-2.5872173) q[0];
rz(0.91122183) q[1];
sx q[1];
rz(-1.8781885) q[1];
sx q[1];
rz(0.30141452) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49839334) q[0];
sx q[0];
rz(-1.7991156) q[0];
sx q[0];
rz(0.96781815) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7784987) q[2];
sx q[2];
rz(-1.1900717) q[2];
sx q[2];
rz(-2.8247339) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.099286) q[1];
sx q[1];
rz(-1.1908682) q[1];
sx q[1];
rz(0.2289339) q[1];
rz(1.9456995) q[3];
sx q[3];
rz(-1.8640567) q[3];
sx q[3];
rz(2.9305262) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.021412795) q[2];
sx q[2];
rz(-2.9605949) q[2];
sx q[2];
rz(0.0069228355) q[2];
rz(-2.210468) q[3];
sx q[3];
rz(-2.0102863) q[3];
sx q[3];
rz(-2.0324223) q[3];
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
x q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6329426) q[0];
sx q[0];
rz(-2.3045492) q[0];
sx q[0];
rz(-1.6078) q[0];
rz(-2.4389229) q[1];
sx q[1];
rz(-0.53121316) q[1];
sx q[1];
rz(0.30219561) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9990197) q[0];
sx q[0];
rz(-2.4732865) q[0];
sx q[0];
rz(-1.7354986) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.40001656) q[2];
sx q[2];
rz(-2.3850394) q[2];
sx q[2];
rz(-0.053002593) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.8540878) q[1];
sx q[1];
rz(-0.64784986) q[1];
sx q[1];
rz(1.2811529) q[1];
rz(-1.4190361) q[3];
sx q[3];
rz(-1.6223063) q[3];
sx q[3];
rz(-0.53799483) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.2493784) q[2];
sx q[2];
rz(-2.0912781) q[2];
sx q[2];
rz(2.2155679) q[2];
rz(1.2146436) q[3];
sx q[3];
rz(-0.49321431) q[3];
sx q[3];
rz(0.27629575) q[3];
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
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72325426) q[0];
sx q[0];
rz(-2.3746018) q[0];
sx q[0];
rz(-1.1085229) q[0];
rz(-0.33379894) q[1];
sx q[1];
rz(-1.8148986) q[1];
sx q[1];
rz(2.0557859) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3876824) q[0];
sx q[0];
rz(-0.22015239) q[0];
sx q[0];
rz(2.7698344) q[0];
rz(0.18098197) q[2];
sx q[2];
rz(-1.1753193) q[2];
sx q[2];
rz(-0.68812319) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.2183509) q[1];
sx q[1];
rz(-0.78303799) q[1];
sx q[1];
rz(-0.79773517) q[1];
x q[2];
rz(-0.13768519) q[3];
sx q[3];
rz(-0.9387278) q[3];
sx q[3];
rz(1.1732303) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.6256025) q[2];
sx q[2];
rz(-0.36449271) q[2];
sx q[2];
rz(-1.1642574) q[2];
rz(-2.8734251) q[3];
sx q[3];
rz(-1.2428913) q[3];
sx q[3];
rz(0.75781649) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5371573) q[0];
sx q[0];
rz(-2.6217961) q[0];
sx q[0];
rz(-1.9926158) q[0];
rz(0.2941429) q[1];
sx q[1];
rz(-1.5584757) q[1];
sx q[1];
rz(-2.099096) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.54759083) q[0];
sx q[0];
rz(-0.85193513) q[0];
sx q[0];
rz(2.6963364) q[0];
rz(0.78977079) q[2];
sx q[2];
rz(-1.4068479) q[2];
sx q[2];
rz(-0.81610926) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.6175872) q[1];
sx q[1];
rz(-1.397314) q[1];
sx q[1];
rz(-2.7702727) q[1];
x q[2];
rz(-0.028826272) q[3];
sx q[3];
rz(-2.2346063) q[3];
sx q[3];
rz(0.038906038) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.5440172) q[2];
sx q[2];
rz(-0.67108265) q[2];
sx q[2];
rz(-1.7768804) q[2];
rz(-2.4890066) q[3];
sx q[3];
rz(-2.3502974) q[3];
sx q[3];
rz(-0.64835382) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
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
rz(-1.6259916) q[0];
sx q[0];
rz(-0.37726548) q[0];
sx q[0];
rz(3.1242477) q[0];
rz(-3.1248202) q[1];
sx q[1];
rz(-0.78712946) q[1];
sx q[1];
rz(-2.9877072) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.34129225) q[0];
sx q[0];
rz(-1.8878536) q[0];
sx q[0];
rz(-1.4912259) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1915386) q[2];
sx q[2];
rz(-2.1851106) q[2];
sx q[2];
rz(1.4221734) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.9388403) q[1];
sx q[1];
rz(-1.2766663) q[1];
sx q[1];
rz(2.7619656) q[1];
x q[2];
rz(3.1176223) q[3];
sx q[3];
rz(-1.4949189) q[3];
sx q[3];
rz(-0.45437231) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.953557) q[2];
sx q[2];
rz(-2.3253658) q[2];
sx q[2];
rz(-0.44450644) q[2];
rz(2.9514173) q[3];
sx q[3];
rz(-1.5835652) q[3];
sx q[3];
rz(1.7688513) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.001215) q[0];
sx q[0];
rz(-1.2495406) q[0];
sx q[0];
rz(1.0706527) q[0];
rz(3.095678) q[1];
sx q[1];
rz(-1.4713947) q[1];
sx q[1];
rz(-2.3505223) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96595461) q[0];
sx q[0];
rz(-1.6583867) q[0];
sx q[0];
rz(1.2344633) q[0];
rz(-pi) q[1];
rz(1.8395109) q[2];
sx q[2];
rz(-1.1860606) q[2];
sx q[2];
rz(-0.045298227) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.6562499) q[1];
sx q[1];
rz(-2.5963077) q[1];
sx q[1];
rz(-2.4052909) q[1];
rz(-pi) q[2];
rz(-1.0656625) q[3];
sx q[3];
rz(-1.6871243) q[3];
sx q[3];
rz(1.1992421) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.4967686) q[2];
sx q[2];
rz(-0.18459979) q[2];
sx q[2];
rz(1.4727288) q[2];
rz(1.6139) q[3];
sx q[3];
rz(-1.0563285) q[3];
sx q[3];
rz(-1.9014026) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
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
rz(1.2033757) q[0];
sx q[0];
rz(-1.4980409) q[0];
sx q[0];
rz(-0.71312755) q[0];
rz(-0.51364246) q[1];
sx q[1];
rz(-1.4331663) q[1];
sx q[1];
rz(-1.6930361) q[1];
rz(1.8088874) q[2];
sx q[2];
rz(-1.0793964) q[2];
sx q[2];
rz(2.4894077) q[2];
rz(0.80908262) q[3];
sx q[3];
rz(-1.2590903) q[3];
sx q[3];
rz(2.1303582) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
