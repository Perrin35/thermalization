OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.91650668) q[0];
sx q[0];
rz(-2.9986311) q[0];
sx q[0];
rz(1.4885055) q[0];
rz(-2.0090964) q[1];
sx q[1];
rz(5.8512591) q[1];
sx q[1];
rz(6.9195256) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4637488) q[0];
sx q[0];
rz(-1.1882458) q[0];
sx q[0];
rz(-0.93884672) q[0];
rz(1.9540399) q[2];
sx q[2];
rz(-1.9343209) q[2];
sx q[2];
rz(0.55628796) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.591037) q[1];
sx q[1];
rz(-1.4466969) q[1];
sx q[1];
rz(-1.9684588) q[1];
x q[2];
rz(-2.431178) q[3];
sx q[3];
rz(-2.1542366) q[3];
sx q[3];
rz(-1.5855011) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.1119614) q[2];
sx q[2];
rz(-2.0165069) q[2];
sx q[2];
rz(-0.65043989) q[2];
rz(1.8247617) q[3];
sx q[3];
rz(-1.7951169) q[3];
sx q[3];
rz(3.0397547) q[3];
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
x q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0874852) q[0];
sx q[0];
rz(-2.5547145) q[0];
sx q[0];
rz(0.45463872) q[0];
rz(0.023160402) q[1];
sx q[1];
rz(-1.4440447) q[1];
sx q[1];
rz(-0.32618162) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5989203) q[0];
sx q[0];
rz(-0.70264757) q[0];
sx q[0];
rz(-2.513078) q[0];
rz(-2.1670879) q[2];
sx q[2];
rz(-2.6428004) q[2];
sx q[2];
rz(1.6287664) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.61039576) q[1];
sx q[1];
rz(-0.67761868) q[1];
sx q[1];
rz(-2.4536831) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0241859) q[3];
sx q[3];
rz(-2.2281584) q[3];
sx q[3];
rz(-0.18169345) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(3.0252016) q[2];
sx q[2];
rz(-1.9390743) q[2];
sx q[2];
rz(-2.1873059) q[2];
rz(1.8918234) q[3];
sx q[3];
rz(-0.3229177) q[3];
sx q[3];
rz(-3.1402816) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8702451) q[0];
sx q[0];
rz(-1.1529237) q[0];
sx q[0];
rz(-1.1767607) q[0];
rz(-2.7765043) q[1];
sx q[1];
rz(-1.8389713) q[1];
sx q[1];
rz(-1.5286068) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5534262) q[0];
sx q[0];
rz(-1.5575101) q[0];
sx q[0];
rz(1.5499328) q[0];
rz(-pi) q[1];
rz(-1.1272548) q[2];
sx q[2];
rz(-1.4993641) q[2];
sx q[2];
rz(-2.254666) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.4499719) q[1];
sx q[1];
rz(-2.9868538) q[1];
sx q[1];
rz(1.7697033) q[1];
rz(0.69832506) q[3];
sx q[3];
rz(-3.0520682) q[3];
sx q[3];
rz(-1.9036628) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.8927346) q[2];
sx q[2];
rz(-0.55344075) q[2];
sx q[2];
rz(-1.152732) q[2];
rz(1.24498) q[3];
sx q[3];
rz(-2.1608519) q[3];
sx q[3];
rz(-2.8583756) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.98274851) q[0];
sx q[0];
rz(-11*pi/12) q[0];
sx q[0];
rz(-0.68956476) q[0];
rz(-0.025029643) q[1];
sx q[1];
rz(-1.6233147) q[1];
sx q[1];
rz(1.7207346) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11224225) q[0];
sx q[0];
rz(-0.94932244) q[0];
sx q[0];
rz(1.6981324) q[0];
x q[1];
rz(1.1406607) q[2];
sx q[2];
rz(-1.4073512) q[2];
sx q[2];
rz(0.30980936) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.6384699) q[1];
sx q[1];
rz(-1.5984189) q[1];
sx q[1];
rz(-2.6383328) q[1];
x q[2];
rz(-0.33735621) q[3];
sx q[3];
rz(-1.3135664) q[3];
sx q[3];
rz(-1.6599865) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.932852) q[2];
sx q[2];
rz(-2.9661861) q[2];
sx q[2];
rz(-1.1669195) q[2];
rz(-0.44421998) q[3];
sx q[3];
rz(-1.2362213) q[3];
sx q[3];
rz(-2.8319632) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2883478) q[0];
sx q[0];
rz(-1.3695559) q[0];
sx q[0];
rz(2.7254768) q[0];
rz(0.5101282) q[1];
sx q[1];
rz(-0.77113873) q[1];
sx q[1];
rz(0.77521926) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.090714723) q[0];
sx q[0];
rz(-0.64656252) q[0];
sx q[0];
rz(-0.334686) q[0];
rz(-pi) q[1];
rz(-2.7363766) q[2];
sx q[2];
rz(-1.7381258) q[2];
sx q[2];
rz(2.2115603) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.0803804) q[1];
sx q[1];
rz(-2.9728372) q[1];
sx q[1];
rz(-2.3024998) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9479423) q[3];
sx q[3];
rz(-2.1324649) q[3];
sx q[3];
rz(0.18581729) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.9267209) q[2];
sx q[2];
rz(-2.1472609) q[2];
sx q[2];
rz(-2.1057687) q[2];
rz(0.18946798) q[3];
sx q[3];
rz(-2.6697956) q[3];
sx q[3];
rz(2.6389879) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8871317) q[0];
sx q[0];
rz(-0.04763617) q[0];
sx q[0];
rz(-2.7857842) q[0];
rz(1.4954781) q[1];
sx q[1];
rz(-0.9551841) q[1];
sx q[1];
rz(-1.1411508) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4116652) q[0];
sx q[0];
rz(-2.228274) q[0];
sx q[0];
rz(-0.2667) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5579434) q[2];
sx q[2];
rz(-2.5299978) q[2];
sx q[2];
rz(2.0919959) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.43737632) q[1];
sx q[1];
rz(-0.17568888) q[1];
sx q[1];
rz(-2.499627) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5969971) q[3];
sx q[3];
rz(-0.43403927) q[3];
sx q[3];
rz(-1.3846014) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.63153875) q[2];
sx q[2];
rz(-0.71439356) q[2];
sx q[2];
rz(-1.7047403) q[2];
rz(-2.0884183) q[3];
sx q[3];
rz(-1.5318233) q[3];
sx q[3];
rz(-2.6385782) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.18008867) q[0];
sx q[0];
rz(-2.8453974) q[0];
sx q[0];
rz(0.12251138) q[0];
rz(-0.65613121) q[1];
sx q[1];
rz(-1.2686814) q[1];
sx q[1];
rz(-1.8001385) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.91514523) q[0];
sx q[0];
rz(-1.6065803) q[0];
sx q[0];
rz(-3.0962178) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7846863) q[2];
sx q[2];
rz(-0.30115899) q[2];
sx q[2];
rz(-0.44277175) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.16640284) q[1];
sx q[1];
rz(-2.4619048) q[1];
sx q[1];
rz(1.1389334) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4571601) q[3];
sx q[3];
rz(-1.1294147) q[3];
sx q[3];
rz(-1.4809004) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.5682257) q[2];
sx q[2];
rz(-1.2218916) q[2];
sx q[2];
rz(-1.4637671) q[2];
rz(-0.73417869) q[3];
sx q[3];
rz(-0.28407431) q[3];
sx q[3];
rz(-1.8370139) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7952591) q[0];
sx q[0];
rz(-1.3672375) q[0];
sx q[0];
rz(2.6829868) q[0];
rz(2.1268225) q[1];
sx q[1];
rz(-1.9472803) q[1];
sx q[1];
rz(-2.0789304) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5073481) q[0];
sx q[0];
rz(-1.5131823) q[0];
sx q[0];
rz(-0.16648023) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3331679) q[2];
sx q[2];
rz(-1.9701517) q[2];
sx q[2];
rz(-2.6724919) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.1717689) q[1];
sx q[1];
rz(-0.14482982) q[1];
sx q[1];
rz(-1.1317731) q[1];
rz(-pi) q[2];
rz(1.0513145) q[3];
sx q[3];
rz(-1.8534582) q[3];
sx q[3];
rz(1.1804898) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.25166544) q[2];
sx q[2];
rz(-1.7532316) q[2];
sx q[2];
rz(1.7903222) q[2];
rz(-1.9225559) q[3];
sx q[3];
rz(-2.7128897) q[3];
sx q[3];
rz(-0.8333227) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.079000533) q[0];
sx q[0];
rz(-1.3840249) q[0];
sx q[0];
rz(2.2341527) q[0];
rz(-0.22706789) q[1];
sx q[1];
rz(-1.0457958) q[1];
sx q[1];
rz(0.89920941) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9470422) q[0];
sx q[0];
rz(-2.3029746) q[0];
sx q[0];
rz(-2.28449) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0631457) q[2];
sx q[2];
rz(-1.4885934) q[2];
sx q[2];
rz(1.7771099) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.9252852) q[1];
sx q[1];
rz(-1.7561963) q[1];
sx q[1];
rz(2.0865296) q[1];
rz(1.6257004) q[3];
sx q[3];
rz(-0.74885741) q[3];
sx q[3];
rz(-2.1532173) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.1850618) q[2];
sx q[2];
rz(-1.1659634) q[2];
sx q[2];
rz(0.60066191) q[2];
rz(2.0447842) q[3];
sx q[3];
rz(-2.5513702) q[3];
sx q[3];
rz(2.2685331) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.7507062) q[0];
sx q[0];
rz(-1.0989256) q[0];
sx q[0];
rz(-0.98439687) q[0];
rz(0.4920494) q[1];
sx q[1];
rz(-1.4130519) q[1];
sx q[1];
rz(1.9459928) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9162684) q[0];
sx q[0];
rz(-1.9231503) q[0];
sx q[0];
rz(2.5114162) q[0];
rz(-pi) q[1];
rz(0.49585988) q[2];
sx q[2];
rz(-0.95520077) q[2];
sx q[2];
rz(0.95814182) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(3.1164536) q[1];
sx q[1];
rz(-0.39158106) q[1];
sx q[1];
rz(1.7844329) q[1];
x q[2];
rz(2.0298084) q[3];
sx q[3];
rz(-0.85235607) q[3];
sx q[3];
rz(1.8000923) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.055858) q[2];
sx q[2];
rz(-1.5208289) q[2];
sx q[2];
rz(-2.1585507) q[2];
rz(2.8899657) q[3];
sx q[3];
rz(-2.7139137) q[3];
sx q[3];
rz(2.1380641) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
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
rz(-0.49878237) q[0];
sx q[0];
rz(-0.97220535) q[0];
sx q[0];
rz(-2.4844949) q[0];
rz(1.2596754) q[1];
sx q[1];
rz(-1.4114264) q[1];
sx q[1];
rz(-2.2687601) q[1];
rz(-0.85763422) q[2];
sx q[2];
rz(-1.7339755) q[2];
sx q[2];
rz(-0.82790464) q[2];
rz(2.8079001) q[3];
sx q[3];
rz(-1.8751388) q[3];
sx q[3];
rz(1.4061385) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
