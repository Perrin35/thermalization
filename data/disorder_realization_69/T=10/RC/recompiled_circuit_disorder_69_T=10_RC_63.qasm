OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.34713137) q[0];
sx q[0];
rz(-1.0153271) q[0];
sx q[0];
rz(0.46749687) q[0];
rz(-0.52019083) q[1];
sx q[1];
rz(-1.3462892) q[1];
sx q[1];
rz(2.4612114) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2517393) q[0];
sx q[0];
rz(-1.7008874) q[0];
sx q[0];
rz(-0.84530172) q[0];
rz(-1.275185) q[2];
sx q[2];
rz(-1.3829074) q[2];
sx q[2];
rz(2.0762028) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.28847028) q[1];
sx q[1];
rz(-2.5783357) q[1];
sx q[1];
rz(2.017615) q[1];
x q[2];
rz(-2.8075571) q[3];
sx q[3];
rz(-1.738027) q[3];
sx q[3];
rz(-2.2807896) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.9900069) q[2];
sx q[2];
rz(-2.1534584) q[2];
sx q[2];
rz(3.0541259) q[2];
rz(-0.72922373) q[3];
sx q[3];
rz(-0.37344033) q[3];
sx q[3];
rz(-0.27174404) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7049578) q[0];
sx q[0];
rz(-0.74902642) q[0];
sx q[0];
rz(1.0013642) q[0];
rz(-2.9691866) q[1];
sx q[1];
rz(-2.0253851) q[1];
sx q[1];
rz(-2.6175245) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66602089) q[0];
sx q[0];
rz(-0.7374987) q[0];
sx q[0];
rz(-2.0185508) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2696813) q[2];
sx q[2];
rz(-1.2362091) q[2];
sx q[2];
rz(1.7988234) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.9429051) q[1];
sx q[1];
rz(-2.6932979) q[1];
sx q[1];
rz(-2.5793736) q[1];
x q[2];
rz(0.53701138) q[3];
sx q[3];
rz(-1.4376663) q[3];
sx q[3];
rz(-0.46768026) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.15741631) q[2];
sx q[2];
rz(-2.6817862) q[2];
sx q[2];
rz(1.3400419) q[2];
rz(-0.79483461) q[3];
sx q[3];
rz(-2.0017616) q[3];
sx q[3];
rz(-0.036858233) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3443417) q[0];
sx q[0];
rz(-0.69673711) q[0];
sx q[0];
rz(0.62477338) q[0];
rz(0.96427381) q[1];
sx q[1];
rz(-2.6565629) q[1];
sx q[1];
rz(-2.952081) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42650578) q[0];
sx q[0];
rz(-0.77930342) q[0];
sx q[0];
rz(-2.7515609) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.465221) q[2];
sx q[2];
rz(-1.9441248) q[2];
sx q[2];
rz(-0.8651274) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.5901958) q[1];
sx q[1];
rz(-1.8101748) q[1];
sx q[1];
rz(1.5134642) q[1];
x q[2];
rz(1.8241006) q[3];
sx q[3];
rz(-2.5251303) q[3];
sx q[3];
rz(2.9951819) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.1069964) q[2];
sx q[2];
rz(-2.2980289) q[2];
sx q[2];
rz(-1.3872046) q[2];
rz(-2.710279) q[3];
sx q[3];
rz(-1.286819) q[3];
sx q[3];
rz(-1.0881902) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65961924) q[0];
sx q[0];
rz(-2.1932333) q[0];
sx q[0];
rz(-1.5455998) q[0];
rz(1.1384456) q[1];
sx q[1];
rz(-0.7754511) q[1];
sx q[1];
rz(2.0514354) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.50463146) q[0];
sx q[0];
rz(-1.22389) q[0];
sx q[0];
rz(2.7142801) q[0];
x q[1];
rz(-0.50699373) q[2];
sx q[2];
rz(-0.88627964) q[2];
sx q[2];
rz(0.13912858) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.66149368) q[1];
sx q[1];
rz(-1.1575932) q[1];
sx q[1];
rz(-0.6306298) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7146043) q[3];
sx q[3];
rz(-1.3823576) q[3];
sx q[3];
rz(-2.6915336) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.99889341) q[2];
sx q[2];
rz(-2.6901851) q[2];
sx q[2];
rz(2.2857655) q[2];
rz(-1.1936197) q[3];
sx q[3];
rz(-1.6222298) q[3];
sx q[3];
rz(2.2284171) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6511433) q[0];
sx q[0];
rz(-0.97390807) q[0];
sx q[0];
rz(2.9918616) q[0];
rz(2.1504452) q[1];
sx q[1];
rz(-1.9350355) q[1];
sx q[1];
rz(-1.1700464) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.67384185) q[0];
sx q[0];
rz(-1.8918599) q[0];
sx q[0];
rz(-2.2348316) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2511475) q[2];
sx q[2];
rz(-1.3157985) q[2];
sx q[2];
rz(-1.6971708) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.7331446) q[1];
sx q[1];
rz(-1.1322347) q[1];
sx q[1];
rz(2.0088197) q[1];
x q[2];
rz(-0.31801362) q[3];
sx q[3];
rz(-2.5807096) q[3];
sx q[3];
rz(-0.47919264) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.871792) q[2];
sx q[2];
rz(-1.5912) q[2];
sx q[2];
rz(-3.0043547) q[2];
rz(1.8042701) q[3];
sx q[3];
rz(-2.5604355) q[3];
sx q[3];
rz(-1.8491245) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.18734922) q[0];
sx q[0];
rz(-1.7631148) q[0];
sx q[0];
rz(1.3504008) q[0];
rz(0.84287914) q[1];
sx q[1];
rz(-2.402669) q[1];
sx q[1];
rz(0.67289105) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.85149375) q[0];
sx q[0];
rz(-0.58304542) q[0];
sx q[0];
rz(-1.3821938) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.024725155) q[2];
sx q[2];
rz(-0.85660663) q[2];
sx q[2];
rz(1.4662454) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.3402965) q[1];
sx q[1];
rz(-2.166689) q[1];
sx q[1];
rz(1.3188386) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6478496) q[3];
sx q[3];
rz(-1.9814081) q[3];
sx q[3];
rz(-1.3239469) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.792753) q[2];
sx q[2];
rz(-1.2092084) q[2];
sx q[2];
rz(2.7065281) q[2];
rz(-1.7815636) q[3];
sx q[3];
rz(-0.74917787) q[3];
sx q[3];
rz(-0.24766651) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.172794) q[0];
sx q[0];
rz(-2.2491169) q[0];
sx q[0];
rz(-2.0794179) q[0];
rz(2.0299714) q[1];
sx q[1];
rz(-1.2373135) q[1];
sx q[1];
rz(-1.4020845) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3081449) q[0];
sx q[0];
rz(-1.2045367) q[0];
sx q[0];
rz(-1.2783947) q[0];
rz(0.76485302) q[2];
sx q[2];
rz(-1.3991038) q[2];
sx q[2];
rz(-2.9219251) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.40560383) q[1];
sx q[1];
rz(-0.39499184) q[1];
sx q[1];
rz(0.40785457) q[1];
x q[2];
rz(0.13774638) q[3];
sx q[3];
rz(-0.98689729) q[3];
sx q[3];
rz(3.0018842) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.4222251) q[2];
sx q[2];
rz(-0.56754595) q[2];
sx q[2];
rz(-2.4613703) q[2];
rz(-0.42823544) q[3];
sx q[3];
rz(-1.8869583) q[3];
sx q[3];
rz(1.359882) q[3];
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
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.632804) q[0];
sx q[0];
rz(-1.7254242) q[0];
sx q[0];
rz(-1.5493786) q[0];
rz(2.8920065) q[1];
sx q[1];
rz(-1.9892178) q[1];
sx q[1];
rz(-0.54135281) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68529785) q[0];
sx q[0];
rz(-1.0066427) q[0];
sx q[0];
rz(1.5575404) q[0];
rz(0.63033732) q[2];
sx q[2];
rz(-0.54140831) q[2];
sx q[2];
rz(-2.2556717) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.036877) q[1];
sx q[1];
rz(-2.2762205) q[1];
sx q[1];
rz(0.55286644) q[1];
rz(-3.1369282) q[3];
sx q[3];
rz(-2.7612673) q[3];
sx q[3];
rz(-0.62957803) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-3.0970739) q[2];
sx q[2];
rz(-2.7842583) q[2];
sx q[2];
rz(-2.3642335) q[2];
rz(-2.2903531) q[3];
sx q[3];
rz(-1.0612396) q[3];
sx q[3];
rz(1.8204934) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7580496) q[0];
sx q[0];
rz(-1.8109011) q[0];
sx q[0];
rz(3.1316277) q[0];
rz(-1.0154356) q[1];
sx q[1];
rz(-2.3773057) q[1];
sx q[1];
rz(-1.6962956) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9706668) q[0];
sx q[0];
rz(-2.2609841) q[0];
sx q[0];
rz(0.70113457) q[0];
rz(-2.2465835) q[2];
sx q[2];
rz(-1.0969321) q[2];
sx q[2];
rz(2.5788139) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.6362793) q[1];
sx q[1];
rz(-1.5552102) q[1];
sx q[1];
rz(-0.63025766) q[1];
rz(0.85083665) q[3];
sx q[3];
rz(-0.81504226) q[3];
sx q[3];
rz(-1.9581865) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.4426667) q[2];
sx q[2];
rz(-2.2103504) q[2];
sx q[2];
rz(-0.60738579) q[2];
rz(1.7025042) q[3];
sx q[3];
rz(-1.3834229) q[3];
sx q[3];
rz(2.3506892) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.33567515) q[0];
sx q[0];
rz(-1.6573925) q[0];
sx q[0];
rz(-0.6138531) q[0];
rz(-2.0954258) q[1];
sx q[1];
rz(-0.26509735) q[1];
sx q[1];
rz(-0.39224958) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0454355) q[0];
sx q[0];
rz(-2.3491612) q[0];
sx q[0];
rz(2.5655454) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.62551542) q[2];
sx q[2];
rz(-1.0872456) q[2];
sx q[2];
rz(0.12550437) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.2892462) q[1];
sx q[1];
rz(-2.1834071) q[1];
sx q[1];
rz(-2.5958519) q[1];
rz(-pi) q[2];
rz(-3.0030389) q[3];
sx q[3];
rz(-0.20692736) q[3];
sx q[3];
rz(3.0272527) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.0739416) q[2];
sx q[2];
rz(-1.757688) q[2];
sx q[2];
rz(-2.5411141) q[2];
rz(-1.0673808) q[3];
sx q[3];
rz(-1.3200656) q[3];
sx q[3];
rz(2.0480806) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4338715) q[0];
sx q[0];
rz(-0.43294551) q[0];
sx q[0];
rz(-1.659163) q[0];
rz(0.99647635) q[1];
sx q[1];
rz(-1.4947718) q[1];
sx q[1];
rz(1.6047118) q[1];
rz(2.773182) q[2];
sx q[2];
rz(-2.340292) q[2];
sx q[2];
rz(3.1030263) q[2];
rz(0.92580933) q[3];
sx q[3];
rz(-0.79173294) q[3];
sx q[3];
rz(-1.790495) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
