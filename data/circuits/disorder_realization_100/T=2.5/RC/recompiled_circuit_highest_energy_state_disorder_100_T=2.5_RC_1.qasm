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
rz(2.9313791) q[0];
sx q[0];
rz(-0.35603324) q[0];
sx q[0];
rz(-1.5176679) q[0];
rz(-1.3104982) q[1];
sx q[1];
rz(-0.85135353) q[1];
sx q[1];
rz(0.89827615) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1318645) q[0];
sx q[0];
rz(-0.72355958) q[0];
sx q[0];
rz(1.0148763) q[0];
rz(-pi) q[1];
rz(2.064205) q[2];
sx q[2];
rz(-1.8283683) q[2];
sx q[2];
rz(0.34573629) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.5643257) q[1];
sx q[1];
rz(-2.942406) q[1];
sx q[1];
rz(2.729862) q[1];
x q[2];
rz(-3.1329915) q[3];
sx q[3];
rz(-1.7349173) q[3];
sx q[3];
rz(-0.65879956) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.5400759) q[2];
sx q[2];
rz(-2.8105152) q[2];
sx q[2];
rz(-2.4430821) q[2];
rz(1.6554333) q[3];
sx q[3];
rz(-2.555002) q[3];
sx q[3];
rz(1.1521888) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0965213) q[0];
sx q[0];
rz(-2.4424545) q[0];
sx q[0];
rz(-0.29749468) q[0];
rz(-0.43111626) q[1];
sx q[1];
rz(-0.86687207) q[1];
sx q[1];
rz(-1.8284304) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0348059) q[0];
sx q[0];
rz(-1.9138096) q[0];
sx q[0];
rz(-2.4136132) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.976491) q[2];
sx q[2];
rz(-1.5303474) q[2];
sx q[2];
rz(-1.8581529) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.6152108) q[1];
sx q[1];
rz(-1.5632249) q[1];
sx q[1];
rz(-2.3442299) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8937469) q[3];
sx q[3];
rz(-1.447926) q[3];
sx q[3];
rz(-2.7932699) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.3678652) q[2];
sx q[2];
rz(-0.44168681) q[2];
sx q[2];
rz(-2.3750677) q[2];
rz(-0.68486989) q[3];
sx q[3];
rz(-1.6133512) q[3];
sx q[3];
rz(0.49055704) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3016475) q[0];
sx q[0];
rz(-1.3293581) q[0];
sx q[0];
rz(2.8369821) q[0];
rz(-2.1412762) q[1];
sx q[1];
rz(-1.1023003) q[1];
sx q[1];
rz(2.2432378) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4340916) q[0];
sx q[0];
rz(-1.487949) q[0];
sx q[0];
rz(1.0368686) q[0];
rz(-pi) q[1];
rz(2.8091891) q[2];
sx q[2];
rz(-0.89611182) q[2];
sx q[2];
rz(1.589366) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.5542986) q[1];
sx q[1];
rz(-1.0815433) q[1];
sx q[1];
rz(0.6502519) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7285181) q[3];
sx q[3];
rz(-2.6135859) q[3];
sx q[3];
rz(-2.4977893) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.0072713) q[2];
sx q[2];
rz(-1.4192702) q[2];
sx q[2];
rz(-2.814494) q[2];
rz(-1.6359811) q[3];
sx q[3];
rz(-1.6749195) q[3];
sx q[3];
rz(-1.3350284) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
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
rz(-0.026767749) q[0];
sx q[0];
rz(-2.6752495) q[0];
sx q[0];
rz(-0.53792167) q[0];
rz(-2.3665358) q[1];
sx q[1];
rz(-1.8102976) q[1];
sx q[1];
rz(2.7937826) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8224768) q[0];
sx q[0];
rz(-0.70981423) q[0];
sx q[0];
rz(0.7546203) q[0];
rz(-1.6762907) q[2];
sx q[2];
rz(-1.6755793) q[2];
sx q[2];
rz(-1.3380761) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.8533142) q[1];
sx q[1];
rz(-2.2396259) q[1];
sx q[1];
rz(-2.1200256) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1600445) q[3];
sx q[3];
rz(-0.42489811) q[3];
sx q[3];
rz(-1.2228257) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.20370087) q[2];
sx q[2];
rz(-2.4163279) q[2];
sx q[2];
rz(-2.4791278) q[2];
rz(1.0558111) q[3];
sx q[3];
rz(-0.31848389) q[3];
sx q[3];
rz(1.2663579) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.98589677) q[0];
sx q[0];
rz(-0.53845423) q[0];
sx q[0];
rz(2.1180617) q[0];
rz(2.0899978) q[1];
sx q[1];
rz(-1.4902427) q[1];
sx q[1];
rz(3.1255299) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15323205) q[0];
sx q[0];
rz(-1.705959) q[0];
sx q[0];
rz(-1.022781) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.93463411) q[2];
sx q[2];
rz(-0.98630465) q[2];
sx q[2];
rz(2.5713657) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.93960436) q[1];
sx q[1];
rz(-1.7919722) q[1];
sx q[1];
rz(2.6084234) q[1];
rz(1.2993126) q[3];
sx q[3];
rz(-1.6405655) q[3];
sx q[3];
rz(1.2742364) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.2992799) q[2];
sx q[2];
rz(-0.88025847) q[2];
sx q[2];
rz(0.24097815) q[2];
rz(-2.0321417) q[3];
sx q[3];
rz(-1.2883319) q[3];
sx q[3];
rz(-0.39458767) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1269162) q[0];
sx q[0];
rz(-0.83634818) q[0];
sx q[0];
rz(-2.8447004) q[0];
rz(-1.9706005) q[1];
sx q[1];
rz(-1.820727) q[1];
sx q[1];
rz(0.5353294) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7479582) q[0];
sx q[0];
rz(-1.7553333) q[0];
sx q[0];
rz(-2.9519597) q[0];
rz(-2.8487327) q[2];
sx q[2];
rz(-0.83838481) q[2];
sx q[2];
rz(2.5590961) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.1254309) q[1];
sx q[1];
rz(-2.6088073) q[1];
sx q[1];
rz(1.1440117) q[1];
x q[2];
rz(-0.84821852) q[3];
sx q[3];
rz(-2.9322538) q[3];
sx q[3];
rz(2.6959553) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.31763306) q[2];
sx q[2];
rz(-2.3667658) q[2];
sx q[2];
rz(2.5114457) q[2];
rz(1.0050425) q[3];
sx q[3];
rz(-2.9229087) q[3];
sx q[3];
rz(0.69766831) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.180535) q[0];
sx q[0];
rz(-0.068004161) q[0];
sx q[0];
rz(-1.8649752) q[0];
rz(-1.0446154) q[1];
sx q[1];
rz(-1.7216543) q[1];
sx q[1];
rz(-0.23342625) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0208527) q[0];
sx q[0];
rz(-2.4197612) q[0];
sx q[0];
rz(1.5068677) q[0];
x q[1];
rz(1.8769924) q[2];
sx q[2];
rz(-0.95149916) q[2];
sx q[2];
rz(-0.77250749) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.9809493) q[1];
sx q[1];
rz(-0.82349524) q[1];
sx q[1];
rz(-1.6232729) q[1];
rz(-0.10226055) q[3];
sx q[3];
rz(-2.5610211) q[3];
sx q[3];
rz(-1.5796661) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.476568) q[2];
sx q[2];
rz(-1.2580405) q[2];
sx q[2];
rz(0.27102077) q[2];
rz(-2.792231) q[3];
sx q[3];
rz(-2.2237399) q[3];
sx q[3];
rz(2.5640986) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2006328) q[0];
sx q[0];
rz(-0.9372434) q[0];
sx q[0];
rz(1.8120793) q[0];
rz(0.26893523) q[1];
sx q[1];
rz(-2.4766998) q[1];
sx q[1];
rz(1.7609133) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.63686812) q[0];
sx q[0];
rz(-2.8177525) q[0];
sx q[0];
rz(-2.25405) q[0];
x q[1];
rz(-2.5455238) q[2];
sx q[2];
rz(-1.5567994) q[2];
sx q[2];
rz(-0.67971855) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.5561292) q[1];
sx q[1];
rz(-1.7335725) q[1];
sx q[1];
rz(2.1750431) q[1];
x q[2];
rz(0.21694048) q[3];
sx q[3];
rz(-0.57121459) q[3];
sx q[3];
rz(-1.2071213) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.5844476) q[2];
sx q[2];
rz(-2.2717805) q[2];
sx q[2];
rz(-2.4073041) q[2];
rz(1.0555438) q[3];
sx q[3];
rz(-1.5922092) q[3];
sx q[3];
rz(-0.9744823) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.664481) q[0];
sx q[0];
rz(-0.25150126) q[0];
sx q[0];
rz(-0.86790458) q[0];
rz(-1.3033298) q[1];
sx q[1];
rz(-1.5918599) q[1];
sx q[1];
rz(-2.910639) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0030913) q[0];
sx q[0];
rz(-1.0338817) q[0];
sx q[0];
rz(-2.0559539) q[0];
rz(-pi) q[1];
rz(-2.7074293) q[2];
sx q[2];
rz(-1.4462398) q[2];
sx q[2];
rz(-1.9366154) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.26511525) q[1];
sx q[1];
rz(-0.5101706) q[1];
sx q[1];
rz(2.6835346) q[1];
rz(1.5695464) q[3];
sx q[3];
rz(-2.2428572) q[3];
sx q[3];
rz(-1.6015893) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.2086198) q[2];
sx q[2];
rz(-1.4459556) q[2];
sx q[2];
rz(-2.7033973) q[2];
rz(-0.70133251) q[3];
sx q[3];
rz(-2.0739136) q[3];
sx q[3];
rz(-0.35593885) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1040795) q[0];
sx q[0];
rz(-2.6689745) q[0];
sx q[0];
rz(-2.0613101) q[0];
rz(-1.0516306) q[1];
sx q[1];
rz(-1.1311421) q[1];
sx q[1];
rz(-1.0952605) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8488377) q[0];
sx q[0];
rz(-1.0612604) q[0];
sx q[0];
rz(1.2408907) q[0];
x q[1];
rz(-3.0313086) q[2];
sx q[2];
rz(-0.37026893) q[2];
sx q[2];
rz(0.25843378) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.238823) q[1];
sx q[1];
rz(-0.84759694) q[1];
sx q[1];
rz(2.9937503) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5445537) q[3];
sx q[3];
rz(-1.1445012) q[3];
sx q[3];
rz(-0.36620127) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-3.0206535) q[2];
sx q[2];
rz(-1.2660657) q[2];
sx q[2];
rz(-0.90547639) q[2];
rz(-0.71546537) q[3];
sx q[3];
rz(-2.4162636) q[3];
sx q[3];
rz(-2.4874617) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
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
rz(-2.7249971) q[0];
sx q[0];
rz(-1.409197) q[0];
sx q[0];
rz(-2.5247164) q[0];
rz(2.0116518) q[1];
sx q[1];
rz(-1.8604953) q[1];
sx q[1];
rz(-1.4166191) q[1];
rz(0.8413419) q[2];
sx q[2];
rz(-2.0615933) q[2];
sx q[2];
rz(1.9398873) q[2];
rz(-0.94983421) q[3];
sx q[3];
rz(-2.8845308) q[3];
sx q[3];
rz(-1.844248) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
