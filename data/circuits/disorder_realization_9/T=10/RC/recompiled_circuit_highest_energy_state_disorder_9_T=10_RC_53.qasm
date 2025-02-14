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
rz(0.23586805) q[0];
sx q[0];
rz(1.9960825) q[0];
sx q[0];
rz(9.3762015) q[0];
rz(-1.6254758) q[1];
sx q[1];
rz(-1.0928417) q[1];
sx q[1];
rz(-2.4258274) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5950105) q[0];
sx q[0];
rz(-1.6002185) q[0];
sx q[0];
rz(-1.741623) q[0];
rz(-pi) q[1];
rz(2.2015436) q[2];
sx q[2];
rz(-2.7995281) q[2];
sx q[2];
rz(-2.2243481) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.29808334) q[1];
sx q[1];
rz(-2.2560826) q[1];
sx q[1];
rz(1.0451026) q[1];
rz(1.1668901) q[3];
sx q[3];
rz(-2.0527186) q[3];
sx q[3];
rz(1.9306077) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.48892659) q[2];
sx q[2];
rz(-0.55997866) q[2];
sx q[2];
rz(-2.0860591) q[2];
rz(-2.7355898) q[3];
sx q[3];
rz(-1.3136274) q[3];
sx q[3];
rz(-2.3199911) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
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
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9838298) q[0];
sx q[0];
rz(-2.1021748) q[0];
sx q[0];
rz(-2.5854172) q[0];
rz(1.3293386) q[1];
sx q[1];
rz(-0.61873299) q[1];
sx q[1];
rz(3.0583196) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0389685) q[0];
sx q[0];
rz(-0.73154615) q[0];
sx q[0];
rz(2.6464173) q[0];
rz(-pi) q[1];
x q[1];
rz(0.57751285) q[2];
sx q[2];
rz(-1.0339875) q[2];
sx q[2];
rz(2.2978738) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.00035380414) q[1];
sx q[1];
rz(-1.7906931) q[1];
sx q[1];
rz(0.45334894) q[1];
x q[2];
rz(0.86822416) q[3];
sx q[3];
rz(-0.69528841) q[3];
sx q[3];
rz(0.27588683) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.9437774) q[2];
sx q[2];
rz(-1.752172) q[2];
sx q[2];
rz(-2.6488292) q[2];
rz(1.0265776) q[3];
sx q[3];
rz(-0.30288282) q[3];
sx q[3];
rz(-1.4125642) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76348412) q[0];
sx q[0];
rz(-0.56270993) q[0];
sx q[0];
rz(-0.43877959) q[0];
rz(1.6890866) q[1];
sx q[1];
rz(-0.32183281) q[1];
sx q[1];
rz(-2.7486393) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1673812) q[0];
sx q[0];
rz(-1.0797473) q[0];
sx q[0];
rz(-0.031662861) q[0];
rz(-1.6268756) q[2];
sx q[2];
rz(-2.6125557) q[2];
sx q[2];
rz(0.38944405) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.23506) q[1];
sx q[1];
rz(-1.9882012) q[1];
sx q[1];
rz(3.0639265) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0370029) q[3];
sx q[3];
rz(-2.4896693) q[3];
sx q[3];
rz(3.1273747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.68985933) q[2];
sx q[2];
rz(-0.84443513) q[2];
sx q[2];
rz(-0.27352697) q[2];
rz(0.61255974) q[3];
sx q[3];
rz(-1.7275683) q[3];
sx q[3];
rz(0.85649049) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3278811) q[0];
sx q[0];
rz(-0.59309816) q[0];
sx q[0];
rz(-2.3408422) q[0];
rz(2.4294991) q[1];
sx q[1];
rz(-1.1199896) q[1];
sx q[1];
rz(-2.6584279) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.33934911) q[0];
sx q[0];
rz(-0.74269811) q[0];
sx q[0];
rz(0.00089641103) q[0];
x q[1];
rz(2.0676548) q[2];
sx q[2];
rz(-2.0570847) q[2];
sx q[2];
rz(-2.4118347) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.8193244) q[1];
sx q[1];
rz(-1.5419648) q[1];
sx q[1];
rz(0.64614702) q[1];
x q[2];
rz(-2.7399446) q[3];
sx q[3];
rz(-0.73026087) q[3];
sx q[3];
rz(0.83056565) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.157865) q[2];
sx q[2];
rz(-0.71061504) q[2];
sx q[2];
rz(-1.857081) q[2];
rz(2.6537248) q[3];
sx q[3];
rz(-1.7043461) q[3];
sx q[3];
rz(1.454486) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3294285) q[0];
sx q[0];
rz(-2.282113) q[0];
sx q[0];
rz(-2.676945) q[0];
rz(0.94841415) q[1];
sx q[1];
rz(-2.3034838) q[1];
sx q[1];
rz(1.4533739) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8888055) q[0];
sx q[0];
rz(-1.545419) q[0];
sx q[0];
rz(-1.5359098) q[0];
x q[1];
rz(-2.862013) q[2];
sx q[2];
rz(-0.76020066) q[2];
sx q[2];
rz(0.81499824) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.075045295) q[1];
sx q[1];
rz(-0.61457026) q[1];
sx q[1];
rz(3.0250508) q[1];
rz(-pi) q[2];
rz(1.0539758) q[3];
sx q[3];
rz(-0.91435104) q[3];
sx q[3];
rz(2.2051728) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.98443085) q[2];
sx q[2];
rz(-1.0474019) q[2];
sx q[2];
rz(2.7471527) q[2];
rz(-1.1557584) q[3];
sx q[3];
rz(-2.211536) q[3];
sx q[3];
rz(1.3084779) q[3];
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
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8406469) q[0];
sx q[0];
rz(-2.7257958) q[0];
sx q[0];
rz(1.0619324) q[0];
rz(-2.1901219) q[1];
sx q[1];
rz(-2.2679592) q[1];
sx q[1];
rz(-1.5207312) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1298286) q[0];
sx q[0];
rz(-1.6239663) q[0];
sx q[0];
rz(-1.8480248) q[0];
x q[1];
rz(-2.4967147) q[2];
sx q[2];
rz(-1.0928003) q[2];
sx q[2];
rz(1.0005815) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.9545867) q[1];
sx q[1];
rz(-2.329064) q[1];
sx q[1];
rz(1.7962667) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3502832) q[3];
sx q[3];
rz(-1.1854544) q[3];
sx q[3];
rz(2.6366608) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.1157397) q[2];
sx q[2];
rz(-1.226475) q[2];
sx q[2];
rz(0.24570492) q[2];
rz(-1.4541516) q[3];
sx q[3];
rz(-1.629963) q[3];
sx q[3];
rz(-0.83034849) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65569735) q[0];
sx q[0];
rz(-0.73807722) q[0];
sx q[0];
rz(2.4468716) q[0];
rz(3.0409536) q[1];
sx q[1];
rz(-1.0357608) q[1];
sx q[1];
rz(0.86520854) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14510205) q[0];
sx q[0];
rz(-1.9887513) q[0];
sx q[0];
rz(1.5254024) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2312282) q[2];
sx q[2];
rz(-1.5558793) q[2];
sx q[2];
rz(-0.93176022) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.45564902) q[1];
sx q[1];
rz(-0.63243094) q[1];
sx q[1];
rz(-1.8258105) q[1];
rz(-pi) q[2];
rz(-0.67460028) q[3];
sx q[3];
rz(-1.1230506) q[3];
sx q[3];
rz(-2.646776) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.18409099) q[2];
sx q[2];
rz(-2.3881193) q[2];
sx q[2];
rz(2.2313879) q[2];
rz(-1.4522067) q[3];
sx q[3];
rz(-1.9096749) q[3];
sx q[3];
rz(-2.6773101) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24377395) q[0];
sx q[0];
rz(-1.047736) q[0];
sx q[0];
rz(0.62328231) q[0];
rz(-2.9858164) q[1];
sx q[1];
rz(-2.364295) q[1];
sx q[1];
rz(0.26330858) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0513228) q[0];
sx q[0];
rz(-1.3557845) q[0];
sx q[0];
rz(3.1306491) q[0];
rz(-pi) q[1];
rz(-1.6313577) q[2];
sx q[2];
rz(-1.4944296) q[2];
sx q[2];
rz(1.5431984) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.0075078) q[1];
sx q[1];
rz(-0.99564394) q[1];
sx q[1];
rz(-1.47912) q[1];
x q[2];
rz(-2.3562779) q[3];
sx q[3];
rz(-1.1668596) q[3];
sx q[3];
rz(2.8304826) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.89173633) q[2];
sx q[2];
rz(-0.24954924) q[2];
sx q[2];
rz(-2.2535394) q[2];
rz(2.213721) q[3];
sx q[3];
rz(-1.1279305) q[3];
sx q[3];
rz(-2.1685062) q[3];
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
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.71858281) q[0];
sx q[0];
rz(-2.6319478) q[0];
sx q[0];
rz(-2.6408559) q[0];
rz(1.1205193) q[1];
sx q[1];
rz(-0.78103939) q[1];
sx q[1];
rz(2.3401071) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5716176) q[0];
sx q[0];
rz(-1.5688251) q[0];
sx q[0];
rz(0.070974102) q[0];
x q[1];
rz(2.4013585) q[2];
sx q[2];
rz(-0.70986219) q[2];
sx q[2];
rz(0.20698337) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.7940143) q[1];
sx q[1];
rz(-1.1433351) q[1];
sx q[1];
rz(-0.68750896) q[1];
rz(-pi) q[2];
rz(2.2128747) q[3];
sx q[3];
rz(-0.31826543) q[3];
sx q[3];
rz(-0.74198328) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.72006172) q[2];
sx q[2];
rz(-1.3214107) q[2];
sx q[2];
rz(-0.23019543) q[2];
rz(0.0097533334) q[3];
sx q[3];
rz(-1.6246656) q[3];
sx q[3];
rz(-2.9856288) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4474051) q[0];
sx q[0];
rz(-1.367584) q[0];
sx q[0];
rz(-2.3719924) q[0];
rz(0.96254483) q[1];
sx q[1];
rz(-1.8237518) q[1];
sx q[1];
rz(-0.98709551) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1709508) q[0];
sx q[0];
rz(-1.0449396) q[0];
sx q[0];
rz(-3.0609511) q[0];
x q[1];
rz(-2.9883698) q[2];
sx q[2];
rz(-1.330687) q[2];
sx q[2];
rz(2.3599412) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.35008123) q[1];
sx q[1];
rz(-1.529585) q[1];
sx q[1];
rz(0.033153127) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0760097) q[3];
sx q[3];
rz(-0.85799137) q[3];
sx q[3];
rz(0.1103758) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.4276245) q[2];
sx q[2];
rz(-0.93728137) q[2];
sx q[2];
rz(-0.06812185) q[2];
rz(2.6918329) q[3];
sx q[3];
rz(-2.7263548) q[3];
sx q[3];
rz(1.4382039) q[3];
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
rz(1.149067) q[0];
sx q[0];
rz(-1.5328007) q[0];
sx q[0];
rz(-1.3971064) q[0];
rz(0.29009157) q[1];
sx q[1];
rz(-0.42872226) q[1];
sx q[1];
rz(0.52679481) q[1];
rz(-1.220096) q[2];
sx q[2];
rz(-1.097622) q[2];
sx q[2];
rz(-1.9672774) q[2];
rz(-0.3984821) q[3];
sx q[3];
rz(-1.553411) q[3];
sx q[3];
rz(-0.51391272) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
