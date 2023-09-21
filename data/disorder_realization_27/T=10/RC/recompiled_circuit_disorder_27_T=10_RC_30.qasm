OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(3.1383837) q[0];
sx q[0];
rz(-2.9870343) q[0];
sx q[0];
rz(-0.69252339) q[0];
rz(-1.2094296) q[1];
sx q[1];
rz(-1.8930607) q[1];
sx q[1];
rz(-1.7564397) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9557578) q[0];
sx q[0];
rz(-1.8896566) q[0];
sx q[0];
rz(2.3715109) q[0];
rz(-pi) q[1];
rz(-1.295747) q[2];
sx q[2];
rz(-2.1056471) q[2];
sx q[2];
rz(-1.5616852) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.79917158) q[1];
sx q[1];
rz(-2.85596) q[1];
sx q[1];
rz(-2.530453) q[1];
x q[2];
rz(2.5296506) q[3];
sx q[3];
rz(-0.7512593) q[3];
sx q[3];
rz(0.9179759) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.8866855) q[2];
sx q[2];
rz(-0.79780769) q[2];
sx q[2];
rz(-2.936426) q[2];
rz(-2.3702879) q[3];
sx q[3];
rz(-2.3588534) q[3];
sx q[3];
rz(-2.0390959) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7339864) q[0];
sx q[0];
rz(-0.74626958) q[0];
sx q[0];
rz(-0.45390391) q[0];
rz(2.1167963) q[1];
sx q[1];
rz(-2.7339934) q[1];
sx q[1];
rz(-1.227238) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62118976) q[0];
sx q[0];
rz(-1.495201) q[0];
sx q[0];
rz(1.7002506) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4279459) q[2];
sx q[2];
rz(-0.64672856) q[2];
sx q[2];
rz(-1.5681859) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.7770575) q[1];
sx q[1];
rz(-0.39452663) q[1];
sx q[1];
rz(-0.72655861) q[1];
rz(1.4144054) q[3];
sx q[3];
rz(-1.3723433) q[3];
sx q[3];
rz(-0.95201492) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(3.1002905) q[2];
sx q[2];
rz(-1.1854478) q[2];
sx q[2];
rz(-2.5741637) q[2];
rz(-0.36519095) q[3];
sx q[3];
rz(-1.7285715) q[3];
sx q[3];
rz(2.173483) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48297468) q[0];
sx q[0];
rz(-0.56476074) q[0];
sx q[0];
rz(-0.89865249) q[0];
rz(-2.1458416) q[1];
sx q[1];
rz(-1.5581222) q[1];
sx q[1];
rz(0.333289) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.22526564) q[0];
sx q[0];
rz(-2.3258381) q[0];
sx q[0];
rz(-2.4832721) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3938445) q[2];
sx q[2];
rz(-1.6672009) q[2];
sx q[2];
rz(-0.63444885) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.052735141) q[1];
sx q[1];
rz(-0.6328859) q[1];
sx q[1];
rz(-1.3120552) q[1];
x q[2];
rz(2.7379052) q[3];
sx q[3];
rz(-1.0567769) q[3];
sx q[3];
rz(-2.6727563) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.4553392) q[2];
sx q[2];
rz(-1.7909966) q[2];
sx q[2];
rz(-1.9906445) q[2];
rz(0.84093705) q[3];
sx q[3];
rz(-2.0261814) q[3];
sx q[3];
rz(-1.9870728) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6999917) q[0];
sx q[0];
rz(-1.5804407) q[0];
sx q[0];
rz(2.4568795) q[0];
rz(2.1060064) q[1];
sx q[1];
rz(-0.5077478) q[1];
sx q[1];
rz(-1.205014) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.372615) q[0];
sx q[0];
rz(-2.2616771) q[0];
sx q[0];
rz(-0.0432424) q[0];
rz(-pi) q[1];
rz(-1.8393458) q[2];
sx q[2];
rz(-0.66187243) q[2];
sx q[2];
rz(-0.72999398) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.4971784) q[1];
sx q[1];
rz(-0.89343151) q[1];
sx q[1];
rz(0.35269423) q[1];
rz(-pi) q[2];
rz(-0.071675008) q[3];
sx q[3];
rz(-0.86719162) q[3];
sx q[3];
rz(3.049831) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.24923199) q[2];
sx q[2];
rz(-1.426733) q[2];
sx q[2];
rz(0.37115804) q[2];
rz(-1.7403729) q[3];
sx q[3];
rz(-0.6597844) q[3];
sx q[3];
rz(-1.1192809) q[3];
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
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.054759653) q[0];
sx q[0];
rz(-2.355447) q[0];
sx q[0];
rz(3.0084685) q[0];
rz(-2.1482824) q[1];
sx q[1];
rz(-1.7555833) q[1];
sx q[1];
rz(-2.5865119) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.426794) q[0];
sx q[0];
rz(-1.5580651) q[0];
sx q[0];
rz(1.8861594) q[0];
x q[1];
rz(1.416989) q[2];
sx q[2];
rz(-0.4193192) q[2];
sx q[2];
rz(-2.9749982) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.75025573) q[1];
sx q[1];
rz(-1.5899982) q[1];
sx q[1];
rz(2.0160497) q[1];
rz(-pi) q[2];
rz(1.5186148) q[3];
sx q[3];
rz(-1.6515886) q[3];
sx q[3];
rz(-2.3063456) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.83539) q[2];
sx q[2];
rz(-1.0058879) q[2];
sx q[2];
rz(3.0026657) q[2];
rz(-0.94240087) q[3];
sx q[3];
rz(-1.5709632) q[3];
sx q[3];
rz(2.5901103) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.54365629) q[0];
sx q[0];
rz(-0.56977001) q[0];
sx q[0];
rz(-0.57975769) q[0];
rz(-3.014091) q[1];
sx q[1];
rz(-1.189905) q[1];
sx q[1];
rz(-1.6019843) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.416076) q[0];
sx q[0];
rz(-2.2726739) q[0];
sx q[0];
rz(-2.8743125) q[0];
rz(-pi) q[1];
rz(0.33467218) q[2];
sx q[2];
rz(-3.0034608) q[2];
sx q[2];
rz(1.3484671) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.08338883) q[1];
sx q[1];
rz(-1.997588) q[1];
sx q[1];
rz(-0.31174387) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.29069889) q[3];
sx q[3];
rz(-2.2346367) q[3];
sx q[3];
rz(1.3568527) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.55398983) q[2];
sx q[2];
rz(-0.25045276) q[2];
sx q[2];
rz(-0.26947752) q[2];
rz(0.23412165) q[3];
sx q[3];
rz(-2.6168489) q[3];
sx q[3];
rz(0.060119303) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44678974) q[0];
sx q[0];
rz(-0.61820522) q[0];
sx q[0];
rz(-0.68429464) q[0];
rz(0.11958312) q[1];
sx q[1];
rz(-1.2922829) q[1];
sx q[1];
rz(-2.6228242) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2001901) q[0];
sx q[0];
rz(-1.4369643) q[0];
sx q[0];
rz(-0.83394136) q[0];
rz(2.9154645) q[2];
sx q[2];
rz(-0.78352189) q[2];
sx q[2];
rz(-2.4353611) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.5468532) q[1];
sx q[1];
rz(-0.93187983) q[1];
sx q[1];
rz(-0.98313318) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.1154237) q[3];
sx q[3];
rz(-1.9713638) q[3];
sx q[3];
rz(2.8699584) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.2542904) q[2];
sx q[2];
rz(-1.7742949) q[2];
sx q[2];
rz(2.5781412) q[2];
rz(3.0900132) q[3];
sx q[3];
rz(-2.0023465) q[3];
sx q[3];
rz(0.95190597) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96034399) q[0];
sx q[0];
rz(-0.5287756) q[0];
sx q[0];
rz(1.7425591) q[0];
rz(-0.78701204) q[1];
sx q[1];
rz(-1.1287289) q[1];
sx q[1];
rz(-0.74434892) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3314914) q[0];
sx q[0];
rz(-1.6252675) q[0];
sx q[0];
rz(-2.9936552) q[0];
rz(-pi) q[1];
rz(2.0853945) q[2];
sx q[2];
rz(-1.4375763) q[2];
sx q[2];
rz(1.0366057) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.3202312) q[1];
sx q[1];
rz(-1.3470955) q[1];
sx q[1];
rz(0.51364586) q[1];
rz(1.5472502) q[3];
sx q[3];
rz(-2.0252725) q[3];
sx q[3];
rz(-2.8794895) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.4259592) q[2];
sx q[2];
rz(-1.7461494) q[2];
sx q[2];
rz(-2.5320833) q[2];
rz(-2.4842747) q[3];
sx q[3];
rz(-0.64353639) q[3];
sx q[3];
rz(-0.26143423) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1881926) q[0];
sx q[0];
rz(-3.0472026) q[0];
sx q[0];
rz(1.6375861) q[0];
rz(1.9001182) q[1];
sx q[1];
rz(-1.1499317) q[1];
sx q[1];
rz(-0.77493587) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3459754) q[0];
sx q[0];
rz(-0.91751639) q[0];
sx q[0];
rz(0.4972636) q[0];
x q[1];
rz(-2.0140892) q[2];
sx q[2];
rz(-1.4294251) q[2];
sx q[2];
rz(-3.1183426) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.6248524) q[1];
sx q[1];
rz(-1.7588741) q[1];
sx q[1];
rz(0.46051689) q[1];
rz(-pi) q[2];
rz(0.14258607) q[3];
sx q[3];
rz(-1.7877868) q[3];
sx q[3];
rz(-2.2246974) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.9541786) q[2];
sx q[2];
rz(-2.9118907) q[2];
sx q[2];
rz(2.9837218) q[2];
rz(1.212451) q[3];
sx q[3];
rz(-2.082943) q[3];
sx q[3];
rz(1.8635748) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0697486) q[0];
sx q[0];
rz(-0.97244111) q[0];
sx q[0];
rz(-2.9272595) q[0];
rz(2.4841323) q[1];
sx q[1];
rz(-2.9174556) q[1];
sx q[1];
rz(-1.0459895) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25046529) q[0];
sx q[0];
rz(-1.9292826) q[0];
sx q[0];
rz(-1.6842708) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0984512) q[2];
sx q[2];
rz(-0.59141814) q[2];
sx q[2];
rz(0.45515781) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.9606049) q[1];
sx q[1];
rz(-1.029403) q[1];
sx q[1];
rz(2.7156746) q[1];
x q[2];
rz(-1.5973741) q[3];
sx q[3];
rz(-1.9760895) q[3];
sx q[3];
rz(1.6534896) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.41632286) q[2];
sx q[2];
rz(-3.0815093) q[2];
sx q[2];
rz(1.0894758) q[2];
rz(-1.5754835) q[3];
sx q[3];
rz(-1.8982866) q[3];
sx q[3];
rz(2.4889448) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2789223) q[0];
sx q[0];
rz(-2.537732) q[0];
sx q[0];
rz(-2.296007) q[0];
rz(-1.6090341) q[1];
sx q[1];
rz(-1.6747723) q[1];
sx q[1];
rz(2.0369045) q[1];
rz(2.7325148) q[2];
sx q[2];
rz(-1.8553875) q[2];
sx q[2];
rz(2.6864048) q[2];
rz(3.0388721) q[3];
sx q[3];
rz(-2.5892047) q[3];
sx q[3];
rz(1.8451167) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
