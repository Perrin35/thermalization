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
rz(-0.66145575) q[0];
sx q[0];
rz(-0.28132004) q[0];
sx q[0];
rz(-1.9982279) q[0];
rz(-2.4885664) q[1];
sx q[1];
rz(-1.8206568) q[1];
sx q[1];
rz(3.1261669) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2020039) q[0];
sx q[0];
rz(-2.1103854) q[0];
sx q[0];
rz(2.5651776) q[0];
rz(-pi) q[1];
rz(-2.562343) q[2];
sx q[2];
rz(-2.0655895) q[2];
sx q[2];
rz(0.65828568) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.8566638) q[1];
sx q[1];
rz(-2.2010815) q[1];
sx q[1];
rz(-2.6394415) q[1];
rz(-pi) q[2];
rz(-0.53807844) q[3];
sx q[3];
rz(-1.4002677) q[3];
sx q[3];
rz(1.9146843) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.780484) q[2];
sx q[2];
rz(-3.1321654) q[2];
sx q[2];
rz(2.0375605) q[2];
rz(0.55936724) q[3];
sx q[3];
rz(-2.1909824) q[3];
sx q[3];
rz(-2.2573788) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1169432) q[0];
sx q[0];
rz(-2.556612) q[0];
sx q[0];
rz(-2.2351433) q[0];
rz(2.8001884) q[1];
sx q[1];
rz(-0.87541348) q[1];
sx q[1];
rz(-0.34814775) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9885555) q[0];
sx q[0];
rz(-1.8796024) q[0];
sx q[0];
rz(-1.3111809) q[0];
rz(-1.6115336) q[2];
sx q[2];
rz(-1.5968644) q[2];
sx q[2];
rz(-0.099310087) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.22869884) q[1];
sx q[1];
rz(-1.4958819) q[1];
sx q[1];
rz(2.1829437) q[1];
rz(-pi) q[2];
rz(2.313638) q[3];
sx q[3];
rz(-1.3997692) q[3];
sx q[3];
rz(0.010874484) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.97588313) q[2];
sx q[2];
rz(-2.7996863) q[2];
sx q[2];
rz(0.085414097) q[2];
rz(-0.092770569) q[3];
sx q[3];
rz(-0.86920357) q[3];
sx q[3];
rz(-0.87576491) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.77466011) q[0];
sx q[0];
rz(-2.7800738) q[0];
sx q[0];
rz(2.7969978) q[0];
rz(-1.5498281) q[1];
sx q[1];
rz(-1.7235618) q[1];
sx q[1];
rz(-1.6835015) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.944779) q[0];
sx q[0];
rz(-0.85374248) q[0];
sx q[0];
rz(-1.6193006) q[0];
rz(-pi) q[1];
rz(0.05942101) q[2];
sx q[2];
rz(-2.0869617) q[2];
sx q[2];
rz(0.29250654) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.5815253) q[1];
sx q[1];
rz(-1.4017644) q[1];
sx q[1];
rz(1.324111) q[1];
rz(-0.52234274) q[3];
sx q[3];
rz(-1.038658) q[3];
sx q[3];
rz(1.6381611) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.6452667) q[2];
sx q[2];
rz(-0.93924773) q[2];
sx q[2];
rz(2.0110896) q[2];
rz(-0.94275236) q[3];
sx q[3];
rz(-1.0470942) q[3];
sx q[3];
rz(1.2070791) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
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
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7154253) q[0];
sx q[0];
rz(-1.796145) q[0];
sx q[0];
rz(-1.6834393) q[0];
rz(-1.0484877) q[1];
sx q[1];
rz(-1.6494992) q[1];
sx q[1];
rz(-0.47007158) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1753569) q[0];
sx q[0];
rz(-2.5585173) q[0];
sx q[0];
rz(-0.52578853) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0993217) q[2];
sx q[2];
rz(-1.8026581) q[2];
sx q[2];
rz(3.0523155) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.9841135) q[1];
sx q[1];
rz(-1.762885) q[1];
sx q[1];
rz(0.68208154) q[1];
x q[2];
rz(0.80959971) q[3];
sx q[3];
rz(-0.86000809) q[3];
sx q[3];
rz(0.87464911) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.5423535) q[2];
sx q[2];
rz(-0.75672823) q[2];
sx q[2];
rz(-2.2139464) q[2];
rz(1.8574235) q[3];
sx q[3];
rz(-1.410306) q[3];
sx q[3];
rz(0.94902432) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1272142) q[0];
sx q[0];
rz(-0.40507409) q[0];
sx q[0];
rz(-0.48144427) q[0];
rz(-0.97775793) q[1];
sx q[1];
rz(-0.6441741) q[1];
sx q[1];
rz(-2.0668623) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46512544) q[0];
sx q[0];
rz(-2.6338086) q[0];
sx q[0];
rz(-0.98498084) q[0];
rz(-0.75047173) q[2];
sx q[2];
rz(-0.40130645) q[2];
sx q[2];
rz(1.4913781) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.8930298) q[1];
sx q[1];
rz(-2.6932635) q[1];
sx q[1];
rz(2.1488229) q[1];
rz(-pi) q[2];
rz(1.8874218) q[3];
sx q[3];
rz(-1.9434557) q[3];
sx q[3];
rz(-2.0697442) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.0398728) q[2];
sx q[2];
rz(-2.0247255) q[2];
sx q[2];
rz(-0.57363415) q[2];
rz(1.6576069) q[3];
sx q[3];
rz(-1.0480169) q[3];
sx q[3];
rz(0.60605961) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1178591) q[0];
sx q[0];
rz(-0.25256279) q[0];
sx q[0];
rz(-0.31164393) q[0];
rz(0.32870865) q[1];
sx q[1];
rz(-1.8054104) q[1];
sx q[1];
rz(0.74105826) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8097654) q[0];
sx q[0];
rz(-0.96967319) q[0];
sx q[0];
rz(2.580216) q[0];
x q[1];
rz(-1.7996721) q[2];
sx q[2];
rz(-1.6308349) q[2];
sx q[2];
rz(2.6031074) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.7007826) q[1];
sx q[1];
rz(-2.4572671) q[1];
sx q[1];
rz(1.6709063) q[1];
x q[2];
rz(1.5719919) q[3];
sx q[3];
rz(-2.5107267) q[3];
sx q[3];
rz(-1.4470755) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.3662423) q[2];
sx q[2];
rz(-0.58997184) q[2];
sx q[2];
rz(0.74014202) q[2];
rz(2.7548693) q[3];
sx q[3];
rz(-2.4155278) q[3];
sx q[3];
rz(2.6630785) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
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
rz(1.9750403) q[0];
sx q[0];
rz(-2.1599025) q[0];
sx q[0];
rz(-0.24328406) q[0];
rz(-2.2443306) q[1];
sx q[1];
rz(-1.5586531) q[1];
sx q[1];
rz(0.74403393) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1420741) q[0];
sx q[0];
rz(-1.653695) q[0];
sx q[0];
rz(-2.7510452) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2076188) q[2];
sx q[2];
rz(-1.2609856) q[2];
sx q[2];
rz(1.241809) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.1330382) q[1];
sx q[1];
rz(-0.65804505) q[1];
sx q[1];
rz(-2.1997019) q[1];
rz(-pi) q[2];
rz(-1.2598557) q[3];
sx q[3];
rz(-2.5937383) q[3];
sx q[3];
rz(2.1825298) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.1870785) q[2];
sx q[2];
rz(-2.3176471) q[2];
sx q[2];
rz(-0.0079060923) q[2];
rz(-1.4963957) q[3];
sx q[3];
rz(-3.0683066) q[3];
sx q[3];
rz(0.50905281) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1169443) q[0];
sx q[0];
rz(-0.032289676) q[0];
sx q[0];
rz(-0.13667983) q[0];
rz(-1.7928803) q[1];
sx q[1];
rz(-1.8486479) q[1];
sx q[1];
rz(2.6332556) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1404554) q[0];
sx q[0];
rz(-1.9656202) q[0];
sx q[0];
rz(3.0009901) q[0];
rz(0.29859297) q[2];
sx q[2];
rz(-2.2246839) q[2];
sx q[2];
rz(-0.25748006) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.66022422) q[1];
sx q[1];
rz(-1.9268039) q[1];
sx q[1];
rz(-0.86467177) q[1];
rz(-1.4517205) q[3];
sx q[3];
rz(-1.649646) q[3];
sx q[3];
rz(2.598473) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.4623744) q[2];
sx q[2];
rz(-2.4269673) q[2];
sx q[2];
rz(-3.083631) q[2];
rz(-2.1884749) q[3];
sx q[3];
rz(-2.0942196) q[3];
sx q[3];
rz(2.5061149) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6111095) q[0];
sx q[0];
rz(-1.4169175) q[0];
sx q[0];
rz(2.2755151) q[0];
rz(-1.7762314) q[1];
sx q[1];
rz(-1.583464) q[1];
sx q[1];
rz(-0.80397111) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.98483368) q[0];
sx q[0];
rz(-1.8560243) q[0];
sx q[0];
rz(-2.3063176) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3525904) q[2];
sx q[2];
rz(-1.8937292) q[2];
sx q[2];
rz(0.6168405) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.70529363) q[1];
sx q[1];
rz(-1.9914522) q[1];
sx q[1];
rz(2.3068025) q[1];
rz(1.9744121) q[3];
sx q[3];
rz(-1.5909373) q[3];
sx q[3];
rz(-2.8668364) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.08635252) q[2];
sx q[2];
rz(-0.1487727) q[2];
sx q[2];
rz(-1.0521592) q[2];
rz(0.85938984) q[3];
sx q[3];
rz(-2.342577) q[3];
sx q[3];
rz(-2.1990282) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4683485) q[0];
sx q[0];
rz(-2.4788661) q[0];
sx q[0];
rz(-2.7224139) q[0];
rz(3.1255417) q[1];
sx q[1];
rz(-1.586986) q[1];
sx q[1];
rz(-3.0174461) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0483748) q[0];
sx q[0];
rz(-1.4182874) q[0];
sx q[0];
rz(2.1270063) q[0];
x q[1];
rz(-2.2474278) q[2];
sx q[2];
rz(-2.8196206) q[2];
sx q[2];
rz(-1.8338667) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(3.1087703) q[1];
sx q[1];
rz(-1.747393) q[1];
sx q[1];
rz(-0.25417166) q[1];
rz(-pi) q[2];
rz(2.9591363) q[3];
sx q[3];
rz(-0.92853755) q[3];
sx q[3];
rz(-0.69891847) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.9496256) q[2];
sx q[2];
rz(-0.73835915) q[2];
sx q[2];
rz(0.16853608) q[2];
rz(-1.656823) q[3];
sx q[3];
rz(-2.6717581) q[3];
sx q[3];
rz(-0.43009871) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
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
rz(-2.1930502) q[0];
sx q[0];
rz(-0.47946231) q[0];
sx q[0];
rz(-0.23647501) q[0];
rz(-0.82462689) q[1];
sx q[1];
rz(-1.5390479) q[1];
sx q[1];
rz(1.8631757) q[1];
rz(2.5599418) q[2];
sx q[2];
rz(-2.2764475) q[2];
sx q[2];
rz(1.225133) q[2];
rz(3.0807224) q[3];
sx q[3];
rz(-1.0247598) q[3];
sx q[3];
rz(1.8714874) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
