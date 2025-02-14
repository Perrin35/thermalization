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
rz(0.99211168) q[0];
sx q[0];
rz(3.8881128) q[0];
sx q[0];
rz(9.9915656) q[0];
rz(-1.8911288) q[1];
sx q[1];
rz(-1.1486147) q[1];
sx q[1];
rz(-2.221938) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.88094372) q[0];
sx q[0];
rz(-2.3365417) q[0];
sx q[0];
rz(-0.22329231) q[0];
rz(-2.5819725) q[2];
sx q[2];
rz(-0.80794789) q[2];
sx q[2];
rz(-1.9576278) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.9360322) q[1];
sx q[1];
rz(-1.6171997) q[1];
sx q[1];
rz(-2.499799) q[1];
rz(-pi) q[2];
rz(-2.9018974) q[3];
sx q[3];
rz(-0.54876941) q[3];
sx q[3];
rz(1.5645998) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.2016242) q[2];
sx q[2];
rz(-2.2087966) q[2];
sx q[2];
rz(2.0882108) q[2];
rz(2.4273704) q[3];
sx q[3];
rz(-2.768399) q[3];
sx q[3];
rz(-1.8478954) q[3];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5193704) q[0];
sx q[0];
rz(-2.433233) q[0];
sx q[0];
rz(-2.569662) q[0];
rz(-1.8427303) q[1];
sx q[1];
rz(-0.79890257) q[1];
sx q[1];
rz(-2.3518708) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.84452) q[0];
sx q[0];
rz(-1.6668586) q[0];
sx q[0];
rz(-1.8734832) q[0];
rz(-1.7310119) q[2];
sx q[2];
rz(-2.2188257) q[2];
sx q[2];
rz(-0.90479484) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.6691298) q[1];
sx q[1];
rz(-1.2910559) q[1];
sx q[1];
rz(-0.6196687) q[1];
rz(-0.59047191) q[3];
sx q[3];
rz(-2.3749224) q[3];
sx q[3];
rz(-1.1696512) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.85084891) q[2];
sx q[2];
rz(-2.3727543) q[2];
sx q[2];
rz(-1.7791465) q[2];
rz(0.29081523) q[3];
sx q[3];
rz(-1.7751866) q[3];
sx q[3];
rz(-3.0349558) q[3];
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
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0027851) q[0];
sx q[0];
rz(-1.0951575) q[0];
sx q[0];
rz(2.441067) q[0];
rz(0.48577148) q[1];
sx q[1];
rz(-2.3120717) q[1];
sx q[1];
rz(0.22451678) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6922263) q[0];
sx q[0];
rz(-2.7590519) q[0];
sx q[0];
rz(2.7367949) q[0];
x q[1];
rz(-1.4204558) q[2];
sx q[2];
rz(-0.58359658) q[2];
sx q[2];
rz(2.4395669) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.74925923) q[1];
sx q[1];
rz(-2.1527004) q[1];
sx q[1];
rz(3.0514293) q[1];
rz(-pi) q[2];
rz(3.0279474) q[3];
sx q[3];
rz(-0.84535852) q[3];
sx q[3];
rz(0.64021969) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.756788) q[2];
sx q[2];
rz(-0.89185682) q[2];
sx q[2];
rz(-0.049081238) q[2];
rz(-0.067042025) q[3];
sx q[3];
rz(-1.1578355) q[3];
sx q[3];
rz(-2.4750347) q[3];
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
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1753569) q[0];
sx q[0];
rz(-0.55713621) q[0];
sx q[0];
rz(0.019158451) q[0];
rz(-1.7823904) q[1];
sx q[1];
rz(-1.4298341) q[1];
sx q[1];
rz(3.137099) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2990103) q[0];
sx q[0];
rz(-1.0753462) q[0];
sx q[0];
rz(2.5240077) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6439444) q[2];
sx q[2];
rz(-1.9714173) q[2];
sx q[2];
rz(1.7965574) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.055306074) q[1];
sx q[1];
rz(-1.9290566) q[1];
sx q[1];
rz(0.18017027) q[1];
x q[2];
rz(2.9380083) q[3];
sx q[3];
rz(-1.4798375) q[3];
sx q[3];
rz(-1.998015) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.4517453) q[2];
sx q[2];
rz(-1.0901901) q[2];
sx q[2];
rz(1.6881662) q[2];
rz(1.8870185) q[3];
sx q[3];
rz(-1.6202319) q[3];
sx q[3];
rz(-2.5793251) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3314064) q[0];
sx q[0];
rz(-1.2368546) q[0];
sx q[0];
rz(-1.1619262) q[0];
rz(-2.3792073) q[1];
sx q[1];
rz(-0.98680174) q[1];
sx q[1];
rz(2.7120178) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9305161) q[0];
sx q[0];
rz(-1.1628224) q[0];
sx q[0];
rz(3.1019449) q[0];
rz(7/(8*pi)) q[2];
sx q[2];
rz(-2.3633011) q[2];
sx q[2];
rz(-0.24519224) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.9813198) q[1];
sx q[1];
rz(-1.9076875) q[1];
sx q[1];
rz(-2.261922) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2240673) q[3];
sx q[3];
rz(-1.297373) q[3];
sx q[3];
rz(1.3454252) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.21542159) q[2];
sx q[2];
rz(-1.8502219) q[2];
sx q[2];
rz(-1.311709) q[2];
rz(-2.6808776) q[3];
sx q[3];
rz(-1.0034794) q[3];
sx q[3];
rz(-0.13319143) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.32955125) q[0];
sx q[0];
rz(-0.1596182) q[0];
sx q[0];
rz(1.106369) q[0];
rz(3.0270992) q[1];
sx q[1];
rz(-1.0527) q[1];
sx q[1];
rz(-3.0879424) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6686443) q[0];
sx q[0];
rz(-1.8264174) q[0];
sx q[0];
rz(-0.0059519569) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.099951115) q[2];
sx q[2];
rz(-2.3315773) q[2];
sx q[2];
rz(-2.9423327) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.1785894) q[1];
sx q[1];
rz(-0.51437639) q[1];
sx q[1];
rz(1.6806755) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1838004) q[3];
sx q[3];
rz(-1.0292064) q[3];
sx q[3];
rz(1.047618) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.8657118) q[2];
sx q[2];
rz(-1.8837453) q[2];
sx q[2];
rz(-2.6344521) q[2];
rz(-1.3745314) q[3];
sx q[3];
rz(-1.5406698) q[3];
sx q[3];
rz(1.9929632) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62992612) q[0];
sx q[0];
rz(-2.0138795) q[0];
sx q[0];
rz(-0.041286904) q[0];
rz(0.73829007) q[1];
sx q[1];
rz(-1.2246882) q[1];
sx q[1];
rz(0.98947492) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3408728) q[0];
sx q[0];
rz(-2.0400908) q[0];
sx q[0];
rz(-1.0020026) q[0];
x q[1];
rz(-0.47669741) q[2];
sx q[2];
rz(-0.39604353) q[2];
sx q[2];
rz(-0.11644289) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.7724645) q[1];
sx q[1];
rz(-0.99812344) q[1];
sx q[1];
rz(2.2989208) q[1];
x q[2];
rz(2.7691844) q[3];
sx q[3];
rz(-2.1914327) q[3];
sx q[3];
rz(-2.2351976) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.5422633) q[2];
sx q[2];
rz(-1.0487391) q[2];
sx q[2];
rz(-0.85912022) q[2];
rz(-0.011693444) q[3];
sx q[3];
rz(-1.499736) q[3];
sx q[3];
rz(3.0288127) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.82666731) q[0];
sx q[0];
rz(-0.59190094) q[0];
sx q[0];
rz(0.23183204) q[0];
rz(0.58206093) q[1];
sx q[1];
rz(-1.3287909) q[1];
sx q[1];
rz(1.9225165) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9909331) q[0];
sx q[0];
rz(-2.5703589) q[0];
sx q[0];
rz(-0.21749638) q[0];
rz(-pi) q[1];
x q[1];
rz(0.14149551) q[2];
sx q[2];
rz(-2.7853051) q[2];
sx q[2];
rz(-2.3959999) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.4640567) q[1];
sx q[1];
rz(-1.7221863) q[1];
sx q[1];
rz(-1.7580126) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1120511) q[3];
sx q[3];
rz(-0.89139056) q[3];
sx q[3];
rz(-0.32493362) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.1845188) q[2];
sx q[2];
rz(-1.3513214) q[2];
sx q[2];
rz(-0.78682023) q[2];
rz(1.5015073) q[3];
sx q[3];
rz(-0.4314751) q[3];
sx q[3];
rz(1.7155581) q[3];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1805434) q[0];
sx q[0];
rz(-1.3732055) q[0];
sx q[0];
rz(-0.94026646) q[0];
rz(2.6967948) q[1];
sx q[1];
rz(-0.85591379) q[1];
sx q[1];
rz(-2.99517) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.525) q[0];
sx q[0];
rz(-1.4448651) q[0];
sx q[0];
rz(2.8406744) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8690879) q[2];
sx q[2];
rz(-1.8159869) q[2];
sx q[2];
rz(1.6134855) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.2487029) q[1];
sx q[1];
rz(-1.8161402) q[1];
sx q[1];
rz(-0.69458436) q[1];
rz(-1.4692471) q[3];
sx q[3];
rz(-2.3595105) q[3];
sx q[3];
rz(2.8080432) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.08708295) q[2];
sx q[2];
rz(-1.4887709) q[2];
sx q[2];
rz(-2.3040237) q[2];
rz(-0.93291035) q[3];
sx q[3];
rz(-0.98740238) q[3];
sx q[3];
rz(-2.3605409) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26302108) q[0];
sx q[0];
rz(-1.2282547) q[0];
sx q[0];
rz(-1.7735057) q[0];
rz(-0.076016501) q[1];
sx q[1];
rz(-2.5612505) q[1];
sx q[1];
rz(-2.0215633) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.66066) q[0];
sx q[0];
rz(-1.2613861) q[0];
sx q[0];
rz(2.8641228) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4330288) q[2];
sx q[2];
rz(-1.8052793) q[2];
sx q[2];
rz(-1.8013148) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.6516797) q[1];
sx q[1];
rz(-1.2800242) q[1];
sx q[1];
rz(-0.099822961) q[1];
x q[2];
rz(1.0088167) q[3];
sx q[3];
rz(-0.19246092) q[3];
sx q[3];
rz(0.61103067) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.40176216) q[2];
sx q[2];
rz(-1.4515452) q[2];
sx q[2];
rz(0.24701992) q[2];
rz(-2.5714827) q[3];
sx q[3];
rz(-0.48738042) q[3];
sx q[3];
rz(-2.0232239) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0433255) q[0];
sx q[0];
rz(-1.387431) q[0];
sx q[0];
rz(-0.36547216) q[0];
rz(-0.098516057) q[1];
sx q[1];
rz(-2.5010074) q[1];
sx q[1];
rz(0.88190257) q[1];
rz(1.1578887) q[2];
sx q[2];
rz(-1.558424) q[2];
sx q[2];
rz(1.3267714) q[2];
rz(-1.8311114) q[3];
sx q[3];
rz(-1.9394939) q[3];
sx q[3];
rz(0.60422411) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
