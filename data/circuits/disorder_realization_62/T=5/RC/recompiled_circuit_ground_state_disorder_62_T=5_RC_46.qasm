OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.5247076) q[0];
sx q[0];
rz(-2.4029713) q[0];
sx q[0];
rz(-0.11697669) q[0];
rz(-2.1891131) q[1];
sx q[1];
rz(-1.4701966) q[1];
sx q[1];
rz(2.267946) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1348717) q[0];
sx q[0];
rz(-2.3804579) q[0];
sx q[0];
rz(1.2740244) q[0];
x q[1];
rz(0.9175542) q[2];
sx q[2];
rz(-0.54346701) q[2];
sx q[2];
rz(2.562491) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.3404939) q[1];
sx q[1];
rz(-1.5135161) q[1];
sx q[1];
rz(-2.1199473) q[1];
rz(-pi) q[2];
x q[2];
rz(0.98883137) q[3];
sx q[3];
rz(-0.47401325) q[3];
sx q[3];
rz(-2.1287166) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.9690507) q[2];
sx q[2];
rz(-1.3323063) q[2];
sx q[2];
rz(1.2782028) q[2];
rz(-1.605002) q[3];
sx q[3];
rz(-1.9032109) q[3];
sx q[3];
rz(-2.8240375) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4441967) q[0];
sx q[0];
rz(-2.9228656) q[0];
sx q[0];
rz(-2.2636407) q[0];
rz(-0.76389337) q[1];
sx q[1];
rz(-2.0593675) q[1];
sx q[1];
rz(-2.0108932) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5661086) q[0];
sx q[0];
rz(-0.16666767) q[0];
sx q[0];
rz(-1.9420293) q[0];
rz(-pi) q[1];
rz(-1.133059) q[2];
sx q[2];
rz(-1.4302876) q[2];
sx q[2];
rz(-0.92012355) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.89040989) q[1];
sx q[1];
rz(-1.5835539) q[1];
sx q[1];
rz(2.3240672) q[1];
rz(-pi) q[2];
x q[2];
rz(0.14826397) q[3];
sx q[3];
rz(-2.5103323) q[3];
sx q[3];
rz(-0.12693044) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.6905602) q[2];
sx q[2];
rz(-1.5105543) q[2];
sx q[2];
rz(-0.65563273) q[2];
rz(1.2304652) q[3];
sx q[3];
rz(-0.96223193) q[3];
sx q[3];
rz(-0.83278304) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3967628) q[0];
sx q[0];
rz(-1.9704882) q[0];
sx q[0];
rz(-2.6387264) q[0];
rz(0.14547959) q[1];
sx q[1];
rz(-1.5301887) q[1];
sx q[1];
rz(-1.9810289) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46782079) q[0];
sx q[0];
rz(-1.2611817) q[0];
sx q[0];
rz(-1.8873812) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7021178) q[2];
sx q[2];
rz(-1.9177556) q[2];
sx q[2];
rz(-2.0974713) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.33659157) q[1];
sx q[1];
rz(-0.9550001) q[1];
sx q[1];
rz(0.47290441) q[1];
rz(-pi) q[2];
rz(2.4907001) q[3];
sx q[3];
rz(-2.7136991) q[3];
sx q[3];
rz(2.9127757) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.9460556) q[2];
sx q[2];
rz(-2.1058407) q[2];
sx q[2];
rz(-0.76988402) q[2];
rz(0.59512538) q[3];
sx q[3];
rz(-2.6502521) q[3];
sx q[3];
rz(0.46014053) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6339517) q[0];
sx q[0];
rz(-2.0946298) q[0];
sx q[0];
rz(1.3557583) q[0];
rz(-0.62375623) q[1];
sx q[1];
rz(-1.535894) q[1];
sx q[1];
rz(2.6502868) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3771916) q[0];
sx q[0];
rz(-2.1024414) q[0];
sx q[0];
rz(1.3926943) q[0];
rz(1.3613961) q[2];
sx q[2];
rz(-1.2970957) q[2];
sx q[2];
rz(1.4865246) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.92190336) q[1];
sx q[1];
rz(-1.8796726) q[1];
sx q[1];
rz(2.3050453) q[1];
x q[2];
rz(-2.5120751) q[3];
sx q[3];
rz(-0.77198085) q[3];
sx q[3];
rz(-0.38956583) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.19646159) q[2];
sx q[2];
rz(-1.6170231) q[2];
sx q[2];
rz(2.5932236) q[2];
rz(1.3301001) q[3];
sx q[3];
rz(-0.29812223) q[3];
sx q[3];
rz(-2.8032081) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2746628) q[0];
sx q[0];
rz(-0.45329705) q[0];
sx q[0];
rz(-2.0931661) q[0];
rz(3.0323845) q[1];
sx q[1];
rz(-1.5501153) q[1];
sx q[1];
rz(2.0352667) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8620548) q[0];
sx q[0];
rz(-1.7906655) q[0];
sx q[0];
rz(0.1105652) q[0];
rz(-pi) q[1];
rz(-2.8064617) q[2];
sx q[2];
rz(-0.8742399) q[2];
sx q[2];
rz(-2.1256688) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.26004654) q[1];
sx q[1];
rz(-1.4235919) q[1];
sx q[1];
rz(1.24415) q[1];
rz(-pi) q[2];
rz(1.2731092) q[3];
sx q[3];
rz(-1.3411689) q[3];
sx q[3];
rz(-2.7252293) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.76935524) q[2];
sx q[2];
rz(-1.2662338) q[2];
sx q[2];
rz(-2.9388169) q[2];
rz(-2.3619704) q[3];
sx q[3];
rz(-1.0896261) q[3];
sx q[3];
rz(2.7705418) q[3];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62696537) q[0];
sx q[0];
rz(-3.096464) q[0];
sx q[0];
rz(0.87920642) q[0];
rz(-2.5458287) q[1];
sx q[1];
rz(-0.87409449) q[1];
sx q[1];
rz(1.8858058) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6539772) q[0];
sx q[0];
rz(-0.10627986) q[0];
sx q[0];
rz(0.52070658) q[0];
rz(2.6379616) q[2];
sx q[2];
rz(-0.63490311) q[2];
sx q[2];
rz(1.1213999) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.3869898) q[1];
sx q[1];
rz(-0.55596724) q[1];
sx q[1];
rz(-0.43773361) q[1];
rz(-0.61579951) q[3];
sx q[3];
rz(-1.7753505) q[3];
sx q[3];
rz(-1.8425065) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.6229728) q[2];
sx q[2];
rz(-0.33766782) q[2];
sx q[2];
rz(-1.4818304) q[2];
rz(1.4431813) q[3];
sx q[3];
rz(-1.3866164) q[3];
sx q[3];
rz(1.1186918) q[3];
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
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54565322) q[0];
sx q[0];
rz(-2.312199) q[0];
sx q[0];
rz(2.948092) q[0];
rz(3.0962443) q[1];
sx q[1];
rz(-1.9769316) q[1];
sx q[1];
rz(-0.68101105) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65359241) q[0];
sx q[0];
rz(-0.18463273) q[0];
sx q[0];
rz(0.49166481) q[0];
rz(0.25410507) q[2];
sx q[2];
rz(-2.0850344) q[2];
sx q[2];
rz(-2.6387173) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.4894708) q[1];
sx q[1];
rz(-0.014111405) q[1];
sx q[1];
rz(-1.997214) q[1];
rz(-pi) q[2];
rz(2.8425334) q[3];
sx q[3];
rz(-2.7388696) q[3];
sx q[3];
rz(-2.2880768) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.8650774) q[2];
sx q[2];
rz(-2.4962208) q[2];
sx q[2];
rz(1.3481677) q[2];
rz(-1.7776141) q[3];
sx q[3];
rz(-1.8230702) q[3];
sx q[3];
rz(1.9510423) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6337223) q[0];
sx q[0];
rz(-0.81955925) q[0];
sx q[0];
rz(0.67935294) q[0];
rz(2.9044115) q[1];
sx q[1];
rz(-0.91865426) q[1];
sx q[1];
rz(-2.0749345) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0268008) q[0];
sx q[0];
rz(-1.566883) q[0];
sx q[0];
rz(0.51030178) q[0];
rz(-pi) q[1];
rz(2.1751498) q[2];
sx q[2];
rz(-1.6487412) q[2];
sx q[2];
rz(-1.9029531) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.9257366) q[1];
sx q[1];
rz(-1.0514823) q[1];
sx q[1];
rz(0.14332736) q[1];
rz(-0.92654631) q[3];
sx q[3];
rz(-2.6253581) q[3];
sx q[3];
rz(-2.981319) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.8827768) q[2];
sx q[2];
rz(-1.2315742) q[2];
sx q[2];
rz(0.41137496) q[2];
rz(1.2808293) q[3];
sx q[3];
rz(-2.4954093) q[3];
sx q[3];
rz(2.083174) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48788747) q[0];
sx q[0];
rz(-1.7937086) q[0];
sx q[0];
rz(0.93061289) q[0];
rz(-0.46512428) q[1];
sx q[1];
rz(-0.5828751) q[1];
sx q[1];
rz(-1.4177657) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0865308) q[0];
sx q[0];
rz(-1.8095922) q[0];
sx q[0];
rz(-2.5361674) q[0];
rz(-pi) q[1];
rz(-0.64196824) q[2];
sx q[2];
rz(-2.1479508) q[2];
sx q[2];
rz(-2.746144) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.237933) q[1];
sx q[1];
rz(-1.9114019) q[1];
sx q[1];
rz(-1.8800859) q[1];
rz(-pi) q[2];
x q[2];
rz(0.11563372) q[3];
sx q[3];
rz(-2.1468785) q[3];
sx q[3];
rz(0.86365684) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.174939) q[2];
sx q[2];
rz(-2.311309) q[2];
sx q[2];
rz(-1.9261599) q[2];
rz(-0.71839607) q[3];
sx q[3];
rz(-1.6738439) q[3];
sx q[3];
rz(0.64232701) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65828085) q[0];
sx q[0];
rz(-2.719306) q[0];
sx q[0];
rz(-1.8102113) q[0];
rz(2.4347958) q[1];
sx q[1];
rz(-1.7168609) q[1];
sx q[1];
rz(-1.4250863) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0646149) q[0];
sx q[0];
rz(-1.1277024) q[0];
sx q[0];
rz(0.71378543) q[0];
rz(-pi) q[1];
rz(-0.91412506) q[2];
sx q[2];
rz(-1.6222092) q[2];
sx q[2];
rz(2.35319) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.4689961) q[1];
sx q[1];
rz(-1.1332268) q[1];
sx q[1];
rz(-2.5621179) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0277068) q[3];
sx q[3];
rz(-0.82172365) q[3];
sx q[3];
rz(2.3253019) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.66296545) q[2];
sx q[2];
rz(-2.4586283) q[2];
sx q[2];
rz(-1.9662201) q[2];
rz(3.1183682) q[3];
sx q[3];
rz(-2.569779) q[3];
sx q[3];
rz(0.2636675) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5671253) q[0];
sx q[0];
rz(-2.1463285) q[0];
sx q[0];
rz(1.2405201) q[0];
rz(2.4868838) q[1];
sx q[1];
rz(-1.2798825) q[1];
sx q[1];
rz(1.9725694) q[1];
rz(-3.0289247) q[2];
sx q[2];
rz(-1.7074399) q[2];
sx q[2];
rz(1.8497333) q[2];
rz(0.67778604) q[3];
sx q[3];
rz(-0.98193632) q[3];
sx q[3];
rz(2.2241398) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
