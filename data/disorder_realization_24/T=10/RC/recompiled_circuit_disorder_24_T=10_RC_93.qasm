OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.55968094) q[0];
sx q[0];
rz(2.0547325) q[0];
sx q[0];
rz(7.6261043) q[0];
rz(1.5496594) q[1];
sx q[1];
rz(-0.078443371) q[1];
sx q[1];
rz(-0.6426386) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0770744) q[0];
sx q[0];
rz(-1.0318349) q[0];
sx q[0];
rz(1.6032739) q[0];
rz(-0.50813714) q[2];
sx q[2];
rz(-1.4027486) q[2];
sx q[2];
rz(2.7009168) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.77036422) q[1];
sx q[1];
rz(-2.9712354) q[1];
sx q[1];
rz(0.78070663) q[1];
rz(-pi) q[2];
rz(-2.9075165) q[3];
sx q[3];
rz(-1.5353068) q[3];
sx q[3];
rz(-1.8518098) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.6671483) q[2];
sx q[2];
rz(-2.4953304) q[2];
sx q[2];
rz(-1.2791963) q[2];
rz(2.4228418) q[3];
sx q[3];
rz(-1.5712534) q[3];
sx q[3];
rz(-0.32354245) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6779125) q[0];
sx q[0];
rz(-1.7872515) q[0];
sx q[0];
rz(-2.1610778) q[0];
rz(-2.9837043) q[1];
sx q[1];
rz(-1.6442559) q[1];
sx q[1];
rz(0.78871361) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0448488) q[0];
sx q[0];
rz(-1.6941438) q[0];
sx q[0];
rz(0.11476536) q[0];
x q[1];
rz(-0.12077232) q[2];
sx q[2];
rz(-2.883652) q[2];
sx q[2];
rz(-2.068012) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.84856725) q[1];
sx q[1];
rz(-1.2663406) q[1];
sx q[1];
rz(0.57363631) q[1];
rz(-pi) q[2];
rz(-2.3716281) q[3];
sx q[3];
rz(-1.0125481) q[3];
sx q[3];
rz(-1.528873) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.7754037) q[2];
sx q[2];
rz(-2.0663694) q[2];
sx q[2];
rz(-0.186084) q[2];
rz(-2.4880593) q[3];
sx q[3];
rz(-1.3862405) q[3];
sx q[3];
rz(-0.11793605) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9300951) q[0];
sx q[0];
rz(-2.1989172) q[0];
sx q[0];
rz(-0.63013664) q[0];
rz(-3.0139626) q[1];
sx q[1];
rz(-2.4928513) q[1];
sx q[1];
rz(0.72174597) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5436514) q[0];
sx q[0];
rz(-1.3184034) q[0];
sx q[0];
rz(2.4734205) q[0];
rz(-pi) q[1];
rz(-0.53775215) q[2];
sx q[2];
rz(-0.76965145) q[2];
sx q[2];
rz(2.4485181) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.1303781) q[1];
sx q[1];
rz(-1.5323258) q[1];
sx q[1];
rz(-0.4442261) q[1];
rz(2.8053022) q[3];
sx q[3];
rz(-0.97067562) q[3];
sx q[3];
rz(2.1093413) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.3815986) q[2];
sx q[2];
rz(-0.95278946) q[2];
sx q[2];
rz(-0.90908137) q[2];
rz(0.51820731) q[3];
sx q[3];
rz(-1.0639233) q[3];
sx q[3];
rz(-2.853493) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5575314) q[0];
sx q[0];
rz(-2.4052305) q[0];
sx q[0];
rz(-3.0990565) q[0];
rz(-0.78035367) q[1];
sx q[1];
rz(-2.6413554) q[1];
sx q[1];
rz(0.76400486) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1799058) q[0];
sx q[0];
rz(-0.67877239) q[0];
sx q[0];
rz(-1.5508482) q[0];
x q[1];
rz(1.6971223) q[2];
sx q[2];
rz(-2.2115123) q[2];
sx q[2];
rz(-1.390552) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.82356794) q[1];
sx q[1];
rz(-1.9293702) q[1];
sx q[1];
rz(1.093868) q[1];
rz(-pi) q[2];
rz(2.0352827) q[3];
sx q[3];
rz(-0.96040695) q[3];
sx q[3];
rz(2.8697517) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.8445231) q[2];
sx q[2];
rz(-1.2244747) q[2];
sx q[2];
rz(0.53304535) q[2];
rz(2.879203) q[3];
sx q[3];
rz(-2.5049987) q[3];
sx q[3];
rz(-0.46245241) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9168636) q[0];
sx q[0];
rz(-2.841195) q[0];
sx q[0];
rz(-1.2438783) q[0];
rz(-2.9149756) q[1];
sx q[1];
rz(-0.77426568) q[1];
sx q[1];
rz(-0.40333834) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6809083) q[0];
sx q[0];
rz(-1.4049238) q[0];
sx q[0];
rz(-0.18780639) q[0];
x q[1];
rz(-2.4558805) q[2];
sx q[2];
rz(-1.2843686) q[2];
sx q[2];
rz(-0.098766947) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.7008608) q[1];
sx q[1];
rz(-1.7248389) q[1];
sx q[1];
rz(2.6101019) q[1];
x q[2];
rz(0.67540695) q[3];
sx q[3];
rz(-1.1408148) q[3];
sx q[3];
rz(0.59795415) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.8205745) q[2];
sx q[2];
rz(-1.0514739) q[2];
sx q[2];
rz(0.78249758) q[2];
rz(-1.1123505) q[3];
sx q[3];
rz(-0.4370884) q[3];
sx q[3];
rz(-2.1330244) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7111506) q[0];
sx q[0];
rz(-0.40200457) q[0];
sx q[0];
rz(-0.45561403) q[0];
rz(3.0420711) q[1];
sx q[1];
rz(-2.0030622) q[1];
sx q[1];
rz(-3.1351556) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.422056) q[0];
sx q[0];
rz(-2.3483843) q[0];
sx q[0];
rz(-0.35736812) q[0];
rz(3.0354584) q[2];
sx q[2];
rz(-1.1960293) q[2];
sx q[2];
rz(2.8449164) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.2575063) q[1];
sx q[1];
rz(-1.3370561) q[1];
sx q[1];
rz(1.9952378) q[1];
rz(2.3901859) q[3];
sx q[3];
rz(-0.91114984) q[3];
sx q[3];
rz(0.71099647) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.7727938) q[2];
sx q[2];
rz(-1.6364748) q[2];
sx q[2];
rz(-2.2951365) q[2];
rz(0.99772292) q[3];
sx q[3];
rz(-0.80934757) q[3];
sx q[3];
rz(-1.458582) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.098175123) q[0];
sx q[0];
rz(-2.257453) q[0];
sx q[0];
rz(-0.055710677) q[0];
rz(-2.3588691) q[1];
sx q[1];
rz(-1.1306154) q[1];
sx q[1];
rz(1.1605211) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6877277) q[0];
sx q[0];
rz(-1.223432) q[0];
sx q[0];
rz(-1.0930644) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2968282) q[2];
sx q[2];
rz(-0.65462199) q[2];
sx q[2];
rz(2.1824333) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.5934505) q[1];
sx q[1];
rz(-1.1957809) q[1];
sx q[1];
rz(-0.96199357) q[1];
rz(1.9556932) q[3];
sx q[3];
rz(-1.2166096) q[3];
sx q[3];
rz(-0.76364309) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.0760076) q[2];
sx q[2];
rz(-2.2202754) q[2];
sx q[2];
rz(0.78545061) q[2];
rz(0.75585946) q[3];
sx q[3];
rz(-1.8463585) q[3];
sx q[3];
rz(0.41539645) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3128368) q[0];
sx q[0];
rz(-0.090933196) q[0];
sx q[0];
rz(-1.1428517) q[0];
rz(1.8354592) q[1];
sx q[1];
rz(-1.9530692) q[1];
sx q[1];
rz(-2.725504) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4131665) q[0];
sx q[0];
rz(-0.51983716) q[0];
sx q[0];
rz(-1.9782938) q[0];
rz(-pi) q[1];
rz(-0.14256723) q[2];
sx q[2];
rz(-2.0349742) q[2];
sx q[2];
rz(0.67205059) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.9302952) q[1];
sx q[1];
rz(-0.43356178) q[1];
sx q[1];
rz(0.3890721) q[1];
rz(-pi) q[2];
rz(1.3564524) q[3];
sx q[3];
rz(-2.6594901) q[3];
sx q[3];
rz(-1.4516423) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.0351506) q[2];
sx q[2];
rz(-1.687655) q[2];
sx q[2];
rz(-2.8184334) q[2];
rz(0.20251003) q[3];
sx q[3];
rz(-1.860362) q[3];
sx q[3];
rz(0.57730738) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
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
rz(-0.37339661) q[0];
sx q[0];
rz(-0.62921262) q[0];
sx q[0];
rz(1.1449822) q[0];
rz(-1.1960944) q[1];
sx q[1];
rz(-2.9856666) q[1];
sx q[1];
rz(-2.6224565) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8828744) q[0];
sx q[0];
rz(-1.1943733) q[0];
sx q[0];
rz(-2.9988078) q[0];
x q[1];
rz(1.8066508) q[2];
sx q[2];
rz(-1.7867076) q[2];
sx q[2];
rz(-2.4049135) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(3.0874333) q[1];
sx q[1];
rz(-2.3490153) q[1];
sx q[1];
rz(1.6967609) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5095652) q[3];
sx q[3];
rz(-2.8497189) q[3];
sx q[3];
rz(-2.8458418) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.3141994) q[2];
sx q[2];
rz(-1.9284356) q[2];
sx q[2];
rz(0.040977565) q[2];
rz(-2.273902) q[3];
sx q[3];
rz(-2.4882081) q[3];
sx q[3];
rz(2.6303671) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.39559078) q[0];
sx q[0];
rz(-0.87930167) q[0];
sx q[0];
rz(1.6145153) q[0];
rz(-1.4279667) q[1];
sx q[1];
rz(-2.7797996) q[1];
sx q[1];
rz(-1.6428927) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.41560995) q[0];
sx q[0];
rz(-1.590953) q[0];
sx q[0];
rz(-1.495342) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8685568) q[2];
sx q[2];
rz(-0.16212633) q[2];
sx q[2];
rz(0.33725421) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.5012706) q[1];
sx q[1];
rz(-0.81201279) q[1];
sx q[1];
rz(-2.0444319) q[1];
rz(-pi) q[2];
rz(-1.4478217) q[3];
sx q[3];
rz(-1.7719367) q[3];
sx q[3];
rz(2.5589383) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.3830118) q[2];
sx q[2];
rz(-2.0643533) q[2];
sx q[2];
rz(2.995058) q[2];
rz(-0.81418973) q[3];
sx q[3];
rz(-0.79629961) q[3];
sx q[3];
rz(2.7248785) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2232589) q[0];
sx q[0];
rz(-1.2467361) q[0];
sx q[0];
rz(0.99714739) q[0];
rz(1.8756443) q[1];
sx q[1];
rz(-1.0026149) q[1];
sx q[1];
rz(1.2276585) q[1];
rz(2.0064034) q[2];
sx q[2];
rz(-2.8786567) q[2];
sx q[2];
rz(-1.565956) q[2];
rz(0.90445789) q[3];
sx q[3];
rz(-2.1790128) q[3];
sx q[3];
rz(2.0603767) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
