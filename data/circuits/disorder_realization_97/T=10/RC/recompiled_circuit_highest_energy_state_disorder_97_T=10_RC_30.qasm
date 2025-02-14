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
rz(2.661854) q[0];
sx q[0];
rz(5.148611) q[0];
sx q[0];
rz(10.892286) q[0];
rz(1.6788586) q[1];
sx q[1];
rz(0.94574133) q[1];
sx q[1];
rz(8.9206817) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5228793) q[0];
sx q[0];
rz(-2.0815432) q[0];
sx q[0];
rz(-0.63453959) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7376704) q[2];
sx q[2];
rz(-1.9144399) q[2];
sx q[2];
rz(-2.0640896) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.9047743) q[1];
sx q[1];
rz(-2.1874551) q[1];
sx q[1];
rz(-0.84994933) q[1];
rz(0.13926718) q[3];
sx q[3];
rz(-0.33815171) q[3];
sx q[3];
rz(2.5168745) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.41245875) q[2];
sx q[2];
rz(-1.2428186) q[2];
sx q[2];
rz(-0.18623713) q[2];
rz(-1.3650182) q[3];
sx q[3];
rz(-0.61187196) q[3];
sx q[3];
rz(1.2046825) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.81011009) q[0];
sx q[0];
rz(-2.6832566) q[0];
sx q[0];
rz(0.052074281) q[0];
rz(-2.1049818) q[1];
sx q[1];
rz(-0.60984817) q[1];
sx q[1];
rz(-0.61526543) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2589073) q[0];
sx q[0];
rz(-1.8947766) q[0];
sx q[0];
rz(0.86470226) q[0];
x q[1];
rz(-1.7024666) q[2];
sx q[2];
rz(-1.3352009) q[2];
sx q[2];
rz(2.5680096) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.4897291) q[1];
sx q[1];
rz(-2.4861838) q[1];
sx q[1];
rz(1.4784527) q[1];
rz(-pi) q[2];
x q[2];
rz(0.21188696) q[3];
sx q[3];
rz(-0.78915262) q[3];
sx q[3];
rz(2.0921538) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.692824) q[2];
sx q[2];
rz(-1.245446) q[2];
sx q[2];
rz(-1.2358865) q[2];
rz(-2.0125194) q[3];
sx q[3];
rz(-0.90714199) q[3];
sx q[3];
rz(-2.3045519) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6579984) q[0];
sx q[0];
rz(-2.3880385) q[0];
sx q[0];
rz(-1.3760706) q[0];
rz(2.6013069) q[1];
sx q[1];
rz(-1.0397725) q[1];
sx q[1];
rz(-1.6859863) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7175563) q[0];
sx q[0];
rz(-0.40923318) q[0];
sx q[0];
rz(0.17828973) q[0];
x q[1];
rz(1.4810782) q[2];
sx q[2];
rz(-1.3680581) q[2];
sx q[2];
rz(-0.81344542) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.27764717) q[1];
sx q[1];
rz(-1.2182115) q[1];
sx q[1];
rz(-2.2318798) q[1];
rz(-pi) q[2];
rz(-0.40850477) q[3];
sx q[3];
rz(-0.91445476) q[3];
sx q[3];
rz(1.0711627) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.49440631) q[2];
sx q[2];
rz(-2.8072085) q[2];
sx q[2];
rz(1.6953267) q[2];
rz(-1.9485731) q[3];
sx q[3];
rz(-2.1665067) q[3];
sx q[3];
rz(1.6993935) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0541662) q[0];
sx q[0];
rz(-2.8717201) q[0];
sx q[0];
rz(-0.42359459) q[0];
rz(-2.1845747) q[1];
sx q[1];
rz(-1.7901763) q[1];
sx q[1];
rz(-3.0439923) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1621057) q[0];
sx q[0];
rz(-0.66363664) q[0];
sx q[0];
rz(-0.62941334) q[0];
rz(-pi) q[1];
rz(-0.1431485) q[2];
sx q[2];
rz(-2.2586056) q[2];
sx q[2];
rz(2.0542415) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.7047119) q[1];
sx q[1];
rz(-1.2578674) q[1];
sx q[1];
rz(1.9347316) q[1];
rz(-pi) q[2];
x q[2];
rz(2.76582) q[3];
sx q[3];
rz(-0.76533459) q[3];
sx q[3];
rz(-0.96804141) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.5264954) q[2];
sx q[2];
rz(-2.6690833) q[2];
sx q[2];
rz(-0.16461593) q[2];
rz(-0.047867157) q[3];
sx q[3];
rz(-1.7367626) q[3];
sx q[3];
rz(-0.78711787) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1614302) q[0];
sx q[0];
rz(-0.4466559) q[0];
sx q[0];
rz(1.584378) q[0];
rz(-0.36852512) q[1];
sx q[1];
rz(-1.8199814) q[1];
sx q[1];
rz(1.6065067) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6045348) q[0];
sx q[0];
rz(-0.7365444) q[0];
sx q[0];
rz(-2.3307072) q[0];
rz(-0.82442762) q[2];
sx q[2];
rz(-1.8608837) q[2];
sx q[2];
rz(2.1739391) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.2531703) q[1];
sx q[1];
rz(-1.0895604) q[1];
sx q[1];
rz(-0.28216823) q[1];
x q[2];
rz(-1.2822578) q[3];
sx q[3];
rz(-1.6464982) q[3];
sx q[3];
rz(2.5385365) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.2313472) q[2];
sx q[2];
rz(-2.7635837) q[2];
sx q[2];
rz(2.0965915) q[2];
rz(0.16799489) q[3];
sx q[3];
rz(-1.955227) q[3];
sx q[3];
rz(0.35429889) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1217693) q[0];
sx q[0];
rz(-2.1124117) q[0];
sx q[0];
rz(-1.9371012) q[0];
rz(0.80351859) q[1];
sx q[1];
rz(-0.81356994) q[1];
sx q[1];
rz(-2.4258851) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2432118) q[0];
sx q[0];
rz(-1.0168494) q[0];
sx q[0];
rz(0.57147632) q[0];
rz(1.7693232) q[2];
sx q[2];
rz(-1.1305222) q[2];
sx q[2];
rz(-2.0581051) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.0625851) q[1];
sx q[1];
rz(-0.91168303) q[1];
sx q[1];
rz(1.7138033) q[1];
x q[2];
rz(-0.29362595) q[3];
sx q[3];
rz(-2.01247) q[3];
sx q[3];
rz(-2.3229118) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.7314926) q[2];
sx q[2];
rz(-0.64971739) q[2];
sx q[2];
rz(1.043172) q[2];
rz(-3.0336174) q[3];
sx q[3];
rz(-2.2004746) q[3];
sx q[3];
rz(-2.030355) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9589979) q[0];
sx q[0];
rz(-1.3866871) q[0];
sx q[0];
rz(-1.2249655) q[0];
rz(0.23172465) q[1];
sx q[1];
rz(-0.8871612) q[1];
sx q[1];
rz(-1.2506332) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.38076048) q[0];
sx q[0];
rz(-0.75801132) q[0];
sx q[0];
rz(-0.17802989) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.022103775) q[2];
sx q[2];
rz(-1.3726418) q[2];
sx q[2];
rz(-1.0304067) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.9385443) q[1];
sx q[1];
rz(-0.37751679) q[1];
sx q[1];
rz(2.50752) q[1];
rz(-pi) q[2];
rz(-0.020718109) q[3];
sx q[3];
rz(-2.585034) q[3];
sx q[3];
rz(-2.8902365) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(3.0687678) q[2];
sx q[2];
rz(-0.67049694) q[2];
sx q[2];
rz(3.0094299) q[2];
rz(-2.4459548) q[3];
sx q[3];
rz(-1.1824824) q[3];
sx q[3];
rz(2.3343991) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
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
rz(0.19145963) q[0];
sx q[0];
rz(-1.5810672) q[0];
sx q[0];
rz(-0.70924846) q[0];
rz(1.6356155) q[1];
sx q[1];
rz(-1.9874856) q[1];
sx q[1];
rz(-2.0567315) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1902311) q[0];
sx q[0];
rz(-1.6282646) q[0];
sx q[0];
rz(-0.60110737) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3368767) q[2];
sx q[2];
rz(-1.6931117) q[2];
sx q[2];
rz(1.027077) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.9429837) q[1];
sx q[1];
rz(-1.3169022) q[1];
sx q[1];
rz(0.47603807) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3832757) q[3];
sx q[3];
rz(-0.061358364) q[3];
sx q[3];
rz(0.96422577) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.6214577) q[2];
sx q[2];
rz(-1.0026714) q[2];
sx q[2];
rz(-1.3168859) q[2];
rz(-0.78553158) q[3];
sx q[3];
rz(-1.9436049) q[3];
sx q[3];
rz(2.7151916) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.921628) q[0];
sx q[0];
rz(-1.6554609) q[0];
sx q[0];
rz(-2.7332136) q[0];
rz(-2.1829103) q[1];
sx q[1];
rz(-0.32902333) q[1];
sx q[1];
rz(1.6201409) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7470541) q[0];
sx q[0];
rz(-1.9346721) q[0];
sx q[0];
rz(-0.25157148) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0651363) q[2];
sx q[2];
rz(-2.5707013) q[2];
sx q[2];
rz(2.5236208) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.1679528) q[1];
sx q[1];
rz(-0.79881891) q[1];
sx q[1];
rz(2.1583907) q[1];
x q[2];
rz(1.0318087) q[3];
sx q[3];
rz(-2.781311) q[3];
sx q[3];
rz(-0.96794766) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.3756322) q[2];
sx q[2];
rz(-1.1124632) q[2];
sx q[2];
rz(-2.5755303) q[2];
rz(1.3426956) q[3];
sx q[3];
rz(-0.415396) q[3];
sx q[3];
rz(-0.28877637) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.600243) q[0];
sx q[0];
rz(-0.039027795) q[0];
sx q[0];
rz(1.4254697) q[0];
rz(-1.0936945) q[1];
sx q[1];
rz(-2.2192571) q[1];
sx q[1];
rz(2.7519382) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3411639) q[0];
sx q[0];
rz(-0.50398705) q[0];
sx q[0];
rz(-2.1227073) q[0];
rz(-pi) q[1];
rz(2.1946743) q[2];
sx q[2];
rz(-1.7830666) q[2];
sx q[2];
rz(-0.52516261) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-3.006778) q[1];
sx q[1];
rz(-1.14535) q[1];
sx q[1];
rz(-3.1158034) q[1];
rz(-pi) q[2];
rz(-1.2663307) q[3];
sx q[3];
rz(-1.0926477) q[3];
sx q[3];
rz(-1.696451) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.23715544) q[2];
sx q[2];
rz(-0.5390141) q[2];
sx q[2];
rz(-1.1971486) q[2];
rz(-0.14687471) q[3];
sx q[3];
rz(-2.4683888) q[3];
sx q[3];
rz(-2.1541514) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7214397) q[0];
sx q[0];
rz(-1.0316105) q[0];
sx q[0];
rz(1.6461865) q[0];
rz(2.738476) q[1];
sx q[1];
rz(-1.3251726) q[1];
sx q[1];
rz(-1.6641738) q[1];
rz(1.4473128) q[2];
sx q[2];
rz(-2.1421943) q[2];
sx q[2];
rz(-1.2669834) q[2];
rz(-1.758214) q[3];
sx q[3];
rz(-1.1915177) q[3];
sx q[3];
rz(2.4080924) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
