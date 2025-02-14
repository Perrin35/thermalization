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
rz(0.95247954) q[1];
sx q[1];
rz(-1.6713961) q[1];
sx q[1];
rz(-2.267946) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7352493) q[0];
sx q[0];
rz(-2.2910718) q[0];
sx q[0];
rz(-0.27168897) q[0];
rz(-pi) q[1];
x q[1];
rz(0.35189907) q[2];
sx q[2];
rz(-1.9939559) q[2];
sx q[2];
rz(2.9911135) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.80471604) q[1];
sx q[1];
rz(-1.0226486) q[1];
sx q[1];
rz(-3.0744661) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1659032) q[3];
sx q[3];
rz(-1.8244074) q[3];
sx q[3];
rz(0.028279956) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.9690507) q[2];
sx q[2];
rz(-1.3323063) q[2];
sx q[2];
rz(1.2782028) q[2];
rz(-1.5365907) q[3];
sx q[3];
rz(-1.9032109) q[3];
sx q[3];
rz(2.8240375) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4441967) q[0];
sx q[0];
rz(-0.2187271) q[0];
sx q[0];
rz(0.87795192) q[0];
rz(0.76389337) q[1];
sx q[1];
rz(-2.0593675) q[1];
sx q[1];
rz(-1.1306995) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3618523) q[0];
sx q[0];
rz(-1.5105783) q[0];
sx q[0];
rz(-1.4152933) q[0];
rz(1.133059) q[2];
sx q[2];
rz(-1.7113051) q[2];
sx q[2];
rz(2.2214691) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.4476026) q[1];
sx q[1];
rz(-0.75335767) q[1];
sx q[1];
rz(-1.5521469) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.14826397) q[3];
sx q[3];
rz(-2.5103323) q[3];
sx q[3];
rz(0.12693044) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.6905602) q[2];
sx q[2];
rz(-1.5105543) q[2];
sx q[2];
rz(0.65563273) q[2];
rz(1.9111274) q[3];
sx q[3];
rz(-0.96223193) q[3];
sx q[3];
rz(0.83278304) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7448298) q[0];
sx q[0];
rz(-1.9704882) q[0];
sx q[0];
rz(0.50286621) q[0];
rz(-0.14547959) q[1];
sx q[1];
rz(-1.611404) q[1];
sx q[1];
rz(-1.9810289) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7877951) q[0];
sx q[0];
rz(-2.702454) q[0];
sx q[0];
rz(2.3697858) q[0];
rz(2.7918589) q[2];
sx q[2];
rz(-1.4473414) q[2];
sx q[2];
rz(-2.6598005) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.38830259) q[1];
sx q[1];
rz(-0.75725746) q[1];
sx q[1];
rz(0.99885916) q[1];
rz(-pi) q[2];
x q[2];
rz(1.840402) q[3];
sx q[3];
rz(-1.2343711) q[3];
sx q[3];
rz(-2.2158282) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.9460556) q[2];
sx q[2];
rz(-1.035752) q[2];
sx q[2];
rz(-0.76988402) q[2];
rz(2.5464673) q[3];
sx q[3];
rz(-0.49134058) q[3];
sx q[3];
rz(-2.6814521) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.50764099) q[0];
sx q[0];
rz(-2.0946298) q[0];
sx q[0];
rz(-1.3557583) q[0];
rz(0.62375623) q[1];
sx q[1];
rz(-1.535894) q[1];
sx q[1];
rz(0.49130586) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89739931) q[0];
sx q[0];
rz(-1.7241052) q[0];
sx q[0];
rz(0.53863948) q[0];
rz(-pi) q[1];
x q[1];
rz(0.63734881) q[2];
sx q[2];
rz(-2.7985811) q[2];
sx q[2];
rz(0.98877871) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.92190336) q[1];
sx q[1];
rz(-1.2619201) q[1];
sx q[1];
rz(2.3050453) q[1];
x q[2];
rz(1.0503429) q[3];
sx q[3];
rz(-0.97176516) q[3];
sx q[3];
rz(1.9584306) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.9451311) q[2];
sx q[2];
rz(-1.6170231) q[2];
sx q[2];
rz(0.54836908) q[2];
rz(1.3301001) q[3];
sx q[3];
rz(-0.29812223) q[3];
sx q[3];
rz(0.33838457) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
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
rz(2.2746628) q[0];
sx q[0];
rz(-2.6882956) q[0];
sx q[0];
rz(1.0484265) q[0];
rz(-0.10920814) q[1];
sx q[1];
rz(-1.5501153) q[1];
sx q[1];
rz(2.0352667) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8745426) q[0];
sx q[0];
rz(-1.4629034) q[0];
sx q[0];
rz(-1.3496198) q[0];
rz(-pi) q[1];
x q[1];
rz(0.335131) q[2];
sx q[2];
rz(-0.8742399) q[2];
sx q[2];
rz(-2.1256688) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.2611003) q[1];
sx q[1];
rz(-1.2478117) q[1];
sx q[1];
rz(-0.1552945) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.901793) q[3];
sx q[3];
rz(-1.8604401) q[3];
sx q[3];
rz(-2.0568796) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.76935524) q[2];
sx q[2];
rz(-1.2662338) q[2];
sx q[2];
rz(0.20277578) q[2];
rz(-0.77962223) q[3];
sx q[3];
rz(-2.0519665) q[3];
sx q[3];
rz(-0.37105086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62696537) q[0];
sx q[0];
rz(-0.045128673) q[0];
sx q[0];
rz(-0.87920642) q[0];
rz(-0.59576398) q[1];
sx q[1];
rz(-2.2674982) q[1];
sx q[1];
rz(-1.2557868) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1308252) q[0];
sx q[0];
rz(-1.4786451) q[0];
sx q[0];
rz(1.6238201) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5686092) q[2];
sx q[2];
rz(-1.8610916) q[2];
sx q[2];
rz(-0.03183768) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.25099157) q[1];
sx q[1];
rz(-1.0724147) q[1];
sx q[1];
rz(-1.8283286) q[1];
rz(-pi) q[2];
x q[2];
rz(1.321927) q[3];
sx q[3];
rz(-0.96967317) q[3];
sx q[3];
rz(2.7271276) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.51861989) q[2];
sx q[2];
rz(-0.33766782) q[2];
sx q[2];
rz(-1.6597623) q[2];
rz(1.4431813) q[3];
sx q[3];
rz(-1.7549763) q[3];
sx q[3];
rz(-1.1186918) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54565322) q[0];
sx q[0];
rz(-0.82939363) q[0];
sx q[0];
rz(0.19350061) q[0];
rz(-0.045348383) q[1];
sx q[1];
rz(-1.1646611) q[1];
sx q[1];
rz(0.68101105) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4880002) q[0];
sx q[0];
rz(-0.18463273) q[0];
sx q[0];
rz(-2.6499278) q[0];
rz(-pi) q[1];
rz(1.0424642) q[2];
sx q[2];
rz(-1.3501423) q[2];
sx q[2];
rz(-2.200732) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.4894708) q[1];
sx q[1];
rz(-0.014111405) q[1];
sx q[1];
rz(-1.997214) q[1];
rz(-pi) q[2];
rz(-2.7549823) q[3];
sx q[3];
rz(-1.6865239) q[3];
sx q[3];
rz(0.9936617) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.27651522) q[2];
sx q[2];
rz(-2.4962208) q[2];
sx q[2];
rz(1.3481677) q[2];
rz(1.7776141) q[3];
sx q[3];
rz(-1.8230702) q[3];
sx q[3];
rz(-1.9510423) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6337223) q[0];
sx q[0];
rz(-0.81955925) q[0];
sx q[0];
rz(-0.67935294) q[0];
rz(-2.9044115) q[1];
sx q[1];
rz(-0.91865426) q[1];
sx q[1];
rz(2.0749345) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54618597) q[0];
sx q[0];
rz(-2.0810938) q[0];
sx q[0];
rz(1.575281) q[0];
rz(-pi) q[1];
x q[1];
rz(0.96644281) q[2];
sx q[2];
rz(-1.4928515) q[2];
sx q[2];
rz(-1.9029531) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.2158561) q[1];
sx q[1];
rz(-2.0901103) q[1];
sx q[1];
rz(-0.14332736) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8130624) q[3];
sx q[3];
rz(-1.1650929) q[3];
sx q[3];
rz(0.87268396) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.25881585) q[2];
sx q[2];
rz(-1.2315742) q[2];
sx q[2];
rz(-0.41137496) q[2];
rz(1.2808293) q[3];
sx q[3];
rz(-0.6461834) q[3];
sx q[3];
rz(-2.083174) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
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
rz(2.6537052) q[0];
sx q[0];
rz(-1.7937086) q[0];
sx q[0];
rz(-2.2109798) q[0];
rz(2.6764684) q[1];
sx q[1];
rz(-0.5828751) q[1];
sx q[1];
rz(-1.4177657) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3534622) q[0];
sx q[0];
rz(-0.98488082) q[0];
sx q[0];
rz(-1.2829554) q[0];
rz(2.3143461) q[2];
sx q[2];
rz(-2.3066024) q[2];
sx q[2];
rz(-0.5448676) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.5807996) q[1];
sx q[1];
rz(-1.2798112) q[1];
sx q[1];
rz(-0.35620226) q[1];
rz(-0.11563372) q[3];
sx q[3];
rz(-0.99471417) q[3];
sx q[3];
rz(-2.2779358) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.96665367) q[2];
sx q[2];
rz(-0.83028364) q[2];
sx q[2];
rz(1.2154328) q[2];
rz(-0.71839607) q[3];
sx q[3];
rz(-1.4677488) q[3];
sx q[3];
rz(2.4992656) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4833118) q[0];
sx q[0];
rz(-0.42228666) q[0];
sx q[0];
rz(1.8102113) q[0];
rz(-2.4347958) q[1];
sx q[1];
rz(-1.7168609) q[1];
sx q[1];
rz(-1.7165064) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0646149) q[0];
sx q[0];
rz(-1.1277024) q[0];
sx q[0];
rz(-2.4278072) q[0];
x q[1];
rz(-1.6548884) q[2];
sx q[2];
rz(-0.65838366) q[2];
sx q[2];
rz(0.7158196) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.37230587) q[1];
sx q[1];
rz(-2.0897749) q[1];
sx q[1];
rz(-1.0610046) q[1];
rz(0.44300191) q[3];
sx q[3];
rz(-0.85369977) q[3];
sx q[3];
rz(-0.19099654) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.4786272) q[2];
sx q[2];
rz(-0.68296432) q[2];
sx q[2];
rz(1.1753725) q[2];
rz(-0.023224467) q[3];
sx q[3];
rz(-2.569779) q[3];
sx q[3];
rz(0.2636675) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5671253) q[0];
sx q[0];
rz(-2.1463285) q[0];
sx q[0];
rz(1.2405201) q[0];
rz(-0.65470882) q[1];
sx q[1];
rz(-1.2798825) q[1];
sx q[1];
rz(1.9725694) q[1];
rz(0.11266795) q[2];
sx q[2];
rz(-1.7074399) q[2];
sx q[2];
rz(1.8497333) q[2];
rz(-0.8169218) q[3];
sx q[3];
rz(-2.2754442) q[3];
sx q[3];
rz(1.2572921) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
