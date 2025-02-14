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
rz(-3.0149674) q[0];
sx q[0];
rz(-1.5746483) q[0];
sx q[0];
rz(-2.6259165) q[0];
rz(-2.9134143) q[1];
sx q[1];
rz(-2.394634) q[1];
sx q[1];
rz(-0.42626122) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.58881271) q[0];
sx q[0];
rz(-2.3511887) q[0];
sx q[0];
rz(-1.8145723) q[0];
x q[1];
rz(0.50349109) q[2];
sx q[2];
rz(-1.8241183) q[2];
sx q[2];
rz(1.5510786) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.1805686) q[1];
sx q[1];
rz(-2.5371309) q[1];
sx q[1];
rz(-0.1015517) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.32483806) q[3];
sx q[3];
rz(-2.1305314) q[3];
sx q[3];
rz(1.7369651) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.2446186) q[2];
sx q[2];
rz(-1.8316869) q[2];
sx q[2];
rz(0.2429602) q[2];
rz(0.36863676) q[3];
sx q[3];
rz(-0.60550767) q[3];
sx q[3];
rz(1.9093556) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.826137) q[0];
sx q[0];
rz(-0.22268) q[0];
sx q[0];
rz(0.36732236) q[0];
rz(-2.3452554) q[1];
sx q[1];
rz(-2.0834736) q[1];
sx q[1];
rz(-2.499089) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.87602654) q[0];
sx q[0];
rz(-2.4867704) q[0];
sx q[0];
rz(0.73237082) q[0];
x q[1];
rz(1.9733866) q[2];
sx q[2];
rz(-1.0559096) q[2];
sx q[2];
rz(-1.8414611) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.5629329) q[1];
sx q[1];
rz(-2.0173434) q[1];
sx q[1];
rz(2.5121653) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.12629892) q[3];
sx q[3];
rz(-1.0923315) q[3];
sx q[3];
rz(2.0984405) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.502304) q[2];
sx q[2];
rz(-2.2059811) q[2];
sx q[2];
rz(-1.3703692) q[2];
rz(0.227452) q[3];
sx q[3];
rz(-1.2496313) q[3];
sx q[3];
rz(-1.1915709) q[3];
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
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.74533904) q[0];
sx q[0];
rz(-2.9662913) q[0];
sx q[0];
rz(-1.9434209) q[0];
rz(-1.0802957) q[1];
sx q[1];
rz(-2.9273169) q[1];
sx q[1];
rz(-0.094873039) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.085233) q[0];
sx q[0];
rz(-2.108639) q[0];
sx q[0];
rz(-2.8406124) q[0];
rz(1.0611141) q[2];
sx q[2];
rz(-2.4154) q[2];
sx q[2];
rz(0.92483556) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.87743041) q[1];
sx q[1];
rz(-2.1070711) q[1];
sx q[1];
rz(0.047604851) q[1];
x q[2];
rz(-2.1839147) q[3];
sx q[3];
rz(-2.6827178) q[3];
sx q[3];
rz(1.6661299) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.0723116) q[2];
sx q[2];
rz(-0.60439622) q[2];
sx q[2];
rz(2.7351232) q[2];
rz(-1.1786002) q[3];
sx q[3];
rz(-1.4402163) q[3];
sx q[3];
rz(-1.9788205) q[3];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5793295) q[0];
sx q[0];
rz(-0.13450204) q[0];
sx q[0];
rz(0.33600268) q[0];
rz(-0.52040368) q[1];
sx q[1];
rz(-0.86417472) q[1];
sx q[1];
rz(0.10428183) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1395499) q[0];
sx q[0];
rz(-0.487953) q[0];
sx q[0];
rz(0.81172184) q[0];
rz(3.0523446) q[2];
sx q[2];
rz(-0.1570905) q[2];
sx q[2];
rz(-2.2721827) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.6809685) q[1];
sx q[1];
rz(-1.2036585) q[1];
sx q[1];
rz(-0.8511833) q[1];
rz(-pi) q[2];
rz(0.23747344) q[3];
sx q[3];
rz(-0.4646315) q[3];
sx q[3];
rz(0.97685087) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.61863724) q[2];
sx q[2];
rz(-1.0982265) q[2];
sx q[2];
rz(1.9970419) q[2];
rz(2.5942904) q[3];
sx q[3];
rz(-1.9132883) q[3];
sx q[3];
rz(2.0539637) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.8892141) q[0];
sx q[0];
rz(-2.9504898) q[0];
sx q[0];
rz(0.55437535) q[0];
rz(2.2303708) q[1];
sx q[1];
rz(-1.2634042) q[1];
sx q[1];
rz(0.30141452) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.75479555) q[0];
sx q[0];
rz(-0.63969958) q[0];
sx q[0];
rz(-1.1819089) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7784987) q[2];
sx q[2];
rz(-1.951521) q[2];
sx q[2];
rz(2.8247339) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.6032981) q[1];
sx q[1];
rz(-0.4406826) q[1];
sx q[1];
rz(2.0875817) q[1];
rz(2.2603792) q[3];
sx q[3];
rz(-2.6699319) q[3];
sx q[3];
rz(0.72615964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(3.1201799) q[2];
sx q[2];
rz(-0.18099774) q[2];
sx q[2];
rz(0.0069228355) q[2];
rz(0.93112469) q[3];
sx q[3];
rz(-2.0102863) q[3];
sx q[3];
rz(-2.0324223) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.50865006) q[0];
sx q[0];
rz(-0.83704346) q[0];
sx q[0];
rz(1.5337926) q[0];
rz(0.7026698) q[1];
sx q[1];
rz(-0.53121316) q[1];
sx q[1];
rz(-2.839397) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.93389171) q[0];
sx q[0];
rz(-2.2284628) q[0];
sx q[0];
rz(0.12872966) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.40001656) q[2];
sx q[2];
rz(-2.3850394) q[2];
sx q[2];
rz(-0.053002593) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.8540878) q[1];
sx q[1];
rz(-2.4937428) q[1];
sx q[1];
rz(-1.2811529) q[1];
rz(-pi) q[2];
rz(1.2421397) q[3];
sx q[3];
rz(-2.9813926) q[3];
sx q[3];
rz(-0.70806187) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.2493784) q[2];
sx q[2];
rz(-1.0503146) q[2];
sx q[2];
rz(2.2155679) q[2];
rz(-1.2146436) q[3];
sx q[3];
rz(-2.6483783) q[3];
sx q[3];
rz(-2.8652969) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-2.4183384) q[0];
sx q[0];
rz(-0.76699081) q[0];
sx q[0];
rz(-1.1085229) q[0];
rz(0.33379894) q[1];
sx q[1];
rz(-1.326694) q[1];
sx q[1];
rz(2.0557859) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7539102) q[0];
sx q[0];
rz(-0.22015239) q[0];
sx q[0];
rz(-0.37175827) q[0];
rz(0.18098197) q[2];
sx q[2];
rz(-1.1753193) q[2];
sx q[2];
rz(-0.68812319) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.8657596) q[1];
sx q[1];
rz(-1.0414904) q[1];
sx q[1];
rz(2.5342026) q[1];
rz(-0.9341888) q[3];
sx q[3];
rz(-1.459834) q[3];
sx q[3];
rz(2.8257089) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.6256025) q[2];
sx q[2];
rz(-2.7770999) q[2];
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
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6044354) q[0];
sx q[0];
rz(-2.6217961) q[0];
sx q[0];
rz(1.9926158) q[0];
rz(-0.2941429) q[1];
sx q[1];
rz(-1.5831169) q[1];
sx q[1];
rz(-2.099096) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1746848) q[0];
sx q[0];
rz(-2.3174441) q[0];
sx q[0];
rz(2.0281726) q[0];
x q[1];
rz(1.3399959) q[2];
sx q[2];
rz(-0.79446213) q[2];
sx q[2];
rz(0.59150254) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.1139026) q[1];
sx q[1];
rz(-1.2053145) q[1];
sx q[1];
rz(1.7566856) q[1];
x q[2];
rz(1.5339666) q[3];
sx q[3];
rz(-2.4772518) q[3];
sx q[3];
rz(0.0078594154) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.5440172) q[2];
sx q[2];
rz(-0.67108265) q[2];
sx q[2];
rz(1.7768804) q[2];
rz(-0.65258604) q[3];
sx q[3];
rz(-0.79129523) q[3];
sx q[3];
rz(-0.64835382) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6259916) q[0];
sx q[0];
rz(-0.37726548) q[0];
sx q[0];
rz(-3.1242477) q[0];
rz(0.01677244) q[1];
sx q[1];
rz(-2.3544632) q[1];
sx q[1];
rz(-0.15388547) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8003004) q[0];
sx q[0];
rz(-1.8878536) q[0];
sx q[0];
rz(1.6503667) q[0];
rz(2.4272404) q[2];
sx q[2];
rz(-1.0754943) q[2];
sx q[2];
rz(2.60204) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.9964136) q[1];
sx q[1];
rz(-0.47582483) q[1];
sx q[1];
rz(-0.68527542) q[1];
rz(3.1176223) q[3];
sx q[3];
rz(-1.6466738) q[3];
sx q[3];
rz(0.45437231) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.1880356) q[2];
sx q[2];
rz(-2.3253658) q[2];
sx q[2];
rz(-0.44450644) q[2];
rz(0.19017531) q[3];
sx q[3];
rz(-1.5580274) q[3];
sx q[3];
rz(-1.3727413) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
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
rz(1.001215) q[0];
sx q[0];
rz(-1.2495406) q[0];
sx q[0];
rz(2.0709399) q[0];
rz(3.095678) q[1];
sx q[1];
rz(-1.4713947) q[1];
sx q[1];
rz(-2.3505223) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5673253) q[0];
sx q[0];
rz(-1.2358032) q[0];
sx q[0];
rz(-3.0488324) q[0];
rz(-pi) q[1];
rz(-2.7439762) q[2];
sx q[2];
rz(-1.8194345) q[2];
sx q[2];
rz(-1.719081) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.715383) q[1];
sx q[1];
rz(-1.2150303) q[1];
sx q[1];
rz(2.7191616) q[1];
rz(-pi) q[2];
x q[2];
rz(0.13274712) q[3];
sx q[3];
rz(-2.0721966) q[3];
sx q[3];
rz(-2.8341334) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.4967686) q[2];
sx q[2];
rz(-2.9569929) q[2];
sx q[2];
rz(-1.4727288) q[2];
rz(1.5276927) q[3];
sx q[3];
rz(-1.0563285) q[3];
sx q[3];
rz(-1.24019) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.938217) q[0];
sx q[0];
rz(-1.6435517) q[0];
sx q[0];
rz(2.4284651) q[0];
rz(-0.51364246) q[1];
sx q[1];
rz(-1.4331663) q[1];
sx q[1];
rz(-1.6930361) q[1];
rz(1.3327053) q[2];
sx q[2];
rz(-2.0621962) q[2];
sx q[2];
rz(-0.65218492) q[2];
rz(2.33251) q[3];
sx q[3];
rz(-1.8825023) q[3];
sx q[3];
rz(-1.0112345) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
