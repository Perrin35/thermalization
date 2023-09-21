OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.38874415) q[0];
sx q[0];
rz(3.677877) q[0];
sx q[0];
rz(10.372547) q[0];
rz(1.8127958) q[1];
sx q[1];
rz(-1.2674018) q[1];
sx q[1];
rz(-1.0277494) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9337024) q[0];
sx q[0];
rz(-0.37651248) q[0];
sx q[0];
rz(-3.0794789) q[0];
rz(-pi) q[1];
x q[1];
rz(2.920354) q[2];
sx q[2];
rz(-1.0916296) q[2];
sx q[2];
rz(-2.0146807) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.4878405) q[1];
sx q[1];
rz(-0.90812212) q[1];
sx q[1];
rz(-2.0648048) q[1];
rz(-pi) q[2];
x q[2];
rz(0.58524744) q[3];
sx q[3];
rz(-0.71551502) q[3];
sx q[3];
rz(1.1701208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.0119005) q[2];
sx q[2];
rz(-1.7069858) q[2];
sx q[2];
rz(-2.091308) q[2];
rz(-1.1132647) q[3];
sx q[3];
rz(-2.2498825) q[3];
sx q[3];
rz(-3.0734857) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
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
rz(3.0691836) q[0];
sx q[0];
rz(-1.8962815) q[0];
sx q[0];
rz(0.29775277) q[0];
rz(0.61966664) q[1];
sx q[1];
rz(-1.0071808) q[1];
sx q[1];
rz(1.108095) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4257658) q[0];
sx q[0];
rz(-0.96164385) q[0];
sx q[0];
rz(-2.8858375) q[0];
x q[1];
rz(-2.7116398) q[2];
sx q[2];
rz(-1.0064831) q[2];
sx q[2];
rz(-2.058409) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.260028) q[1];
sx q[1];
rz(-1.9238872) q[1];
sx q[1];
rz(-0.41207037) q[1];
x q[2];
rz(2.8262029) q[3];
sx q[3];
rz(-0.72761256) q[3];
sx q[3];
rz(1.7244347) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.6790598) q[2];
sx q[2];
rz(-0.97683895) q[2];
sx q[2];
rz(-2.1662946) q[2];
rz(-2.2235928) q[3];
sx q[3];
rz(-1.2851597) q[3];
sx q[3];
rz(0.29618922) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.179203) q[0];
sx q[0];
rz(-2.2900892) q[0];
sx q[0];
rz(0.54291022) q[0];
rz(-0.88223282) q[1];
sx q[1];
rz(-2.0062607) q[1];
sx q[1];
rz(0.96484819) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0438) q[0];
sx q[0];
rz(-1.291073) q[0];
sx q[0];
rz(0.33491896) q[0];
rz(1.6134878) q[2];
sx q[2];
rz(-1.8674769) q[2];
sx q[2];
rz(-1.3670849) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.0286897) q[1];
sx q[1];
rz(-1.5064872) q[1];
sx q[1];
rz(1.4494004) q[1];
x q[2];
rz(0.19823719) q[3];
sx q[3];
rz(-1.1612411) q[3];
sx q[3];
rz(2.4026681) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.09459153) q[2];
sx q[2];
rz(-2.5307405) q[2];
sx q[2];
rz(1.1331406) q[2];
rz(-0.23162332) q[3];
sx q[3];
rz(-1.8685721) q[3];
sx q[3];
rz(2.384322) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.27947458) q[0];
sx q[0];
rz(-3.1311488) q[0];
sx q[0];
rz(-1.7650771) q[0];
rz(0.51849413) q[1];
sx q[1];
rz(-1.8771749) q[1];
sx q[1];
rz(0.24212295) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3926485) q[0];
sx q[0];
rz(-2.2376275) q[0];
sx q[0];
rz(-1.3632266) q[0];
rz(-pi) q[1];
rz(1.5577199) q[2];
sx q[2];
rz(-2.0003951) q[2];
sx q[2];
rz(1.0193046) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.8072847) q[1];
sx q[1];
rz(-1.1242928) q[1];
sx q[1];
rz(-1.7590894) q[1];
rz(-0.41800349) q[3];
sx q[3];
rz(-1.6998708) q[3];
sx q[3];
rz(-1.2152745) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.018192856) q[2];
sx q[2];
rz(-0.96863666) q[2];
sx q[2];
rz(-0.094853178) q[2];
rz(-1.3421966) q[3];
sx q[3];
rz(-1.3972524) q[3];
sx q[3];
rz(0.43911394) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3577394) q[0];
sx q[0];
rz(-1.8308324) q[0];
sx q[0];
rz(0.36079303) q[0];
rz(1.3882673) q[1];
sx q[1];
rz(-1.8107982) q[1];
sx q[1];
rz(1.1345908) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15123385) q[0];
sx q[0];
rz(-0.84999527) q[0];
sx q[0];
rz(2.4609341) q[0];
rz(-pi) q[1];
x q[1];
rz(0.88419948) q[2];
sx q[2];
rz(-0.53495416) q[2];
sx q[2];
rz(-1.3605489) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.1281631) q[1];
sx q[1];
rz(-1.1268106) q[1];
sx q[1];
rz(-1.3982989) q[1];
x q[2];
rz(2.9633425) q[3];
sx q[3];
rz(-2.6986487) q[3];
sx q[3];
rz(0.49803842) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.0052884) q[2];
sx q[2];
rz(-1.7207928) q[2];
sx q[2];
rz(0.57265442) q[2];
rz(0.92875656) q[3];
sx q[3];
rz(-0.52162617) q[3];
sx q[3];
rz(-1.9740392) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.3271493) q[0];
sx q[0];
rz(-2.0454018) q[0];
sx q[0];
rz(-1.2493398) q[0];
rz(1.918474) q[1];
sx q[1];
rz(-1.616281) q[1];
sx q[1];
rz(-1.9893601) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6845247) q[0];
sx q[0];
rz(-0.10604924) q[0];
sx q[0];
rz(1.690879) q[0];
rz(-pi) q[1];
rz(0.67201519) q[2];
sx q[2];
rz(-0.61215559) q[2];
sx q[2];
rz(-1.9415346) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.8730948) q[1];
sx q[1];
rz(-0.22646204) q[1];
sx q[1];
rz(2.6576659) q[1];
rz(-pi) q[2];
rz(0.96529393) q[3];
sx q[3];
rz(-0.98635841) q[3];
sx q[3];
rz(0.30129978) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.46889177) q[2];
sx q[2];
rz(-1.9273309) q[2];
sx q[2];
rz(0.99096283) q[2];
rz(-2.4937566) q[3];
sx q[3];
rz(-0.9655374) q[3];
sx q[3];
rz(-2.794054) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25794849) q[0];
sx q[0];
rz(-2.914496) q[0];
sx q[0];
rz(-3.0793072) q[0];
rz(0.1858055) q[1];
sx q[1];
rz(-1.4567016) q[1];
sx q[1];
rz(0.3947765) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.018935238) q[0];
sx q[0];
rz(-1.1766953) q[0];
sx q[0];
rz(1.8271441) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4418418) q[2];
sx q[2];
rz(-0.47669461) q[2];
sx q[2];
rz(-1.2303908) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.3697606) q[1];
sx q[1];
rz(-0.50060111) q[1];
sx q[1];
rz(2.6066149) q[1];
rz(-1.6077605) q[3];
sx q[3];
rz(-1.7241524) q[3];
sx q[3];
rz(-1.5330029) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.64849598) q[2];
sx q[2];
rz(-0.33005565) q[2];
sx q[2];
rz(-0.27302343) q[2];
rz(1.3027044) q[3];
sx q[3];
rz(-1.3132934) q[3];
sx q[3];
rz(0.31204143) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7438695) q[0];
sx q[0];
rz(-2.0070772) q[0];
sx q[0];
rz(-1.8564818) q[0];
rz(-1.6400736) q[1];
sx q[1];
rz(-1.39095) q[1];
sx q[1];
rz(1.3407019) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96269755) q[0];
sx q[0];
rz(-2.6216051) q[0];
sx q[0];
rz(-0.87332256) q[0];
rz(-pi) q[1];
rz(2.4457744) q[2];
sx q[2];
rz(-1.6098445) q[2];
sx q[2];
rz(-0.85862904) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.50780523) q[1];
sx q[1];
rz(-1.9314737) q[1];
sx q[1];
rz(0.37640576) q[1];
rz(-pi) q[2];
rz(-2.0128485) q[3];
sx q[3];
rz(-1.8574517) q[3];
sx q[3];
rz(-0.096404508) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.1061219) q[2];
sx q[2];
rz(-0.80646986) q[2];
sx q[2];
rz(1.0236053) q[2];
rz(2.9566531) q[3];
sx q[3];
rz(-2.7513294) q[3];
sx q[3];
rz(-0.24188724) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0446562) q[0];
sx q[0];
rz(-2.1450295) q[0];
sx q[0];
rz(1.5203083) q[0];
rz(-2.8114491) q[1];
sx q[1];
rz(-1.9338927) q[1];
sx q[1];
rz(0.83713371) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5569816) q[0];
sx q[0];
rz(-0.71174445) q[0];
sx q[0];
rz(2.8296489) q[0];
rz(-pi) q[1];
rz(-2.7021072) q[2];
sx q[2];
rz(-0.62289933) q[2];
sx q[2];
rz(-0.57657951) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.7652119) q[1];
sx q[1];
rz(-1.4188671) q[1];
sx q[1];
rz(-3.0357749) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8505441) q[3];
sx q[3];
rz(-2.0488727) q[3];
sx q[3];
rz(-0.92791286) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.5403486) q[2];
sx q[2];
rz(-2.0566172) q[2];
sx q[2];
rz(-0.17962757) q[2];
rz(2.1458697) q[3];
sx q[3];
rz(-1.8928173) q[3];
sx q[3];
rz(1.8306336) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3778465) q[0];
sx q[0];
rz(-2.7959931) q[0];
sx q[0];
rz(-1.0572222) q[0];
rz(0.10748848) q[1];
sx q[1];
rz(-1.2534393) q[1];
sx q[1];
rz(-2.1616139) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7769023) q[0];
sx q[0];
rz(-1.6216462) q[0];
sx q[0];
rz(-1.6965673) q[0];
rz(-pi) q[1];
rz(2.6238742) q[2];
sx q[2];
rz(-2.0132408) q[2];
sx q[2];
rz(2.8874318) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.500463) q[1];
sx q[1];
rz(-2.6424721) q[1];
sx q[1];
rz(-1.0541381) q[1];
rz(-pi) q[2];
rz(1.2782709) q[3];
sx q[3];
rz(-0.33698002) q[3];
sx q[3];
rz(1.302759) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.8367299) q[2];
sx q[2];
rz(-1.2048081) q[2];
sx q[2];
rz(-1.0277964) q[2];
rz(1.3868388) q[3];
sx q[3];
rz(-1.3091062) q[3];
sx q[3];
rz(-0.28361472) q[3];
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
x q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86823157) q[0];
sx q[0];
rz(-1.0366476) q[0];
sx q[0];
rz(-1.114053) q[0];
rz(-2.3095619) q[1];
sx q[1];
rz(-0.46453005) q[1];
sx q[1];
rz(0.66418905) q[1];
rz(-1.2912512) q[2];
sx q[2];
rz(-0.57689473) q[2];
sx q[2];
rz(2.9070791) q[2];
rz(-1.9527312) q[3];
sx q[3];
rz(-1.100913) q[3];
sx q[3];
rz(2.037896) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
