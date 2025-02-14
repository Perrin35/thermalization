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
rz(0.089601547) q[0];
sx q[0];
rz(-0.90553415) q[0];
sx q[0];
rz(2.9517458) q[0];
rz(-0.39659652) q[1];
sx q[1];
rz(-2.8009669) q[1];
sx q[1];
rz(0.61000282) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6515577) q[0];
sx q[0];
rz(-1.8076979) q[0];
sx q[0];
rz(1.8238571) q[0];
rz(-pi) q[1];
rz(0.18596639) q[2];
sx q[2];
rz(-1.1347316) q[2];
sx q[2];
rz(0.036967335) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.09874092) q[1];
sx q[1];
rz(-1.0611649) q[1];
sx q[1];
rz(-2.3364728) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0656432) q[3];
sx q[3];
rz(-0.54022721) q[3];
sx q[3];
rz(-1.6151419) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.6923328) q[2];
sx q[2];
rz(-1.924694) q[2];
sx q[2];
rz(2.8999691) q[2];
rz(3.0786476) q[3];
sx q[3];
rz(-1.2017622) q[3];
sx q[3];
rz(-0.81301779) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6908506) q[0];
sx q[0];
rz(-1.2797322) q[0];
sx q[0];
rz(-2.8513841) q[0];
rz(-0.33503512) q[1];
sx q[1];
rz(-2.2033043) q[1];
sx q[1];
rz(2.8163574) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7656416) q[0];
sx q[0];
rz(-2.0360569) q[0];
sx q[0];
rz(1.6164533) q[0];
rz(-pi) q[1];
x q[1];
rz(0.26518719) q[2];
sx q[2];
rz(-2.2329805) q[2];
sx q[2];
rz(-1.7005077) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.27694459) q[1];
sx q[1];
rz(-1.7200302) q[1];
sx q[1];
rz(1.7407038) q[1];
x q[2];
rz(0.70230647) q[3];
sx q[3];
rz(-0.78506535) q[3];
sx q[3];
rz(0.64083523) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.30514303) q[2];
sx q[2];
rz(-2.2025351) q[2];
sx q[2];
rz(-1.1391501) q[2];
rz(-2.3818453) q[3];
sx q[3];
rz(-0.097948827) q[3];
sx q[3];
rz(0.37794149) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.29001319) q[0];
sx q[0];
rz(-0.72999287) q[0];
sx q[0];
rz(2.3082025) q[0];
rz(-3.1381798) q[1];
sx q[1];
rz(-0.5178057) q[1];
sx q[1];
rz(-1.980967) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.072642751) q[0];
sx q[0];
rz(-2.8153689) q[0];
sx q[0];
rz(1.6618927) q[0];
rz(-1.7951444) q[2];
sx q[2];
rz(-2.2065711) q[2];
sx q[2];
rz(1.1368542) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.034027122) q[1];
sx q[1];
rz(-2.9068397) q[1];
sx q[1];
rz(1.5723133) q[1];
rz(-pi) q[2];
rz(2.0187601) q[3];
sx q[3];
rz(-2.7395757) q[3];
sx q[3];
rz(2.1340919) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.87355906) q[2];
sx q[2];
rz(-1.4159091) q[2];
sx q[2];
rz(-0.1669008) q[2];
rz(-1.9514826) q[3];
sx q[3];
rz(-2.8968865) q[3];
sx q[3];
rz(-2.0580097) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13389182) q[0];
sx q[0];
rz(-1.5950483) q[0];
sx q[0];
rz(-2.290945) q[0];
rz(-0.8029241) q[1];
sx q[1];
rz(-1.9839169) q[1];
sx q[1];
rz(-2.0096774) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.23861215) q[0];
sx q[0];
rz(-1.0762666) q[0];
sx q[0];
rz(0.20574768) q[0];
x q[1];
rz(0.19271199) q[2];
sx q[2];
rz(-1.2817749) q[2];
sx q[2];
rz(-0.98622903) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.67805659) q[1];
sx q[1];
rz(-2.4812323) q[1];
sx q[1];
rz(0.44594958) q[1];
x q[2];
rz(2.6904039) q[3];
sx q[3];
rz(-0.78253096) q[3];
sx q[3];
rz(0.53451371) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.5256646) q[2];
sx q[2];
rz(-1.6388845) q[2];
sx q[2];
rz(-2.8611503) q[2];
rz(2.7403455) q[3];
sx q[3];
rz(-0.29956996) q[3];
sx q[3];
rz(-1.7846599) q[3];
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
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2111557) q[0];
sx q[0];
rz(-1.7481952) q[0];
sx q[0];
rz(2.9537971) q[0];
rz(-0.64741778) q[1];
sx q[1];
rz(-0.96962601) q[1];
sx q[1];
rz(2.5513249) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85344106) q[0];
sx q[0];
rz(-0.36124215) q[0];
sx q[0];
rz(-1.4807184) q[0];
x q[1];
rz(0.78863229) q[2];
sx q[2];
rz(-1.8344387) q[2];
sx q[2];
rz(-0.90197998) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.62377159) q[1];
sx q[1];
rz(-1.5989027) q[1];
sx q[1];
rz(-3.0610316) q[1];
rz(-pi) q[2];
x q[2];
rz(0.66450232) q[3];
sx q[3];
rz(-1.4034617) q[3];
sx q[3];
rz(-1.4897084) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.65842998) q[2];
sx q[2];
rz(-1.1700609) q[2];
sx q[2];
rz(-2.9675193) q[2];
rz(-0.46711323) q[3];
sx q[3];
rz(-0.73254782) q[3];
sx q[3];
rz(-2.9113801) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1751404) q[0];
sx q[0];
rz(-1.3850965) q[0];
sx q[0];
rz(1.0873644) q[0];
rz(2.9804136) q[1];
sx q[1];
rz(-1.168707) q[1];
sx q[1];
rz(0.93343121) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3573489) q[0];
sx q[0];
rz(-2.0406561) q[0];
sx q[0];
rz(-2.78746) q[0];
rz(0.12045355) q[2];
sx q[2];
rz(-1.9423123) q[2];
sx q[2];
rz(0.62830433) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.63042917) q[1];
sx q[1];
rz(-1.9952718) q[1];
sx q[1];
rz(1.9349348) q[1];
x q[2];
rz(2.5640413) q[3];
sx q[3];
rz(-2.1114919) q[3];
sx q[3];
rz(0.78839984) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.33438385) q[2];
sx q[2];
rz(-1.7024567) q[2];
sx q[2];
rz(1.8857694) q[2];
rz(0.0028336023) q[3];
sx q[3];
rz(-2.5815559) q[3];
sx q[3];
rz(-0.66633666) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.46347076) q[0];
sx q[0];
rz(-3.1246298) q[0];
sx q[0];
rz(2.8609138) q[0];
rz(-0.30411389) q[1];
sx q[1];
rz(-2.3801443) q[1];
sx q[1];
rz(-1.4842518) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4481135) q[0];
sx q[0];
rz(-1.5511314) q[0];
sx q[0];
rz(-3.124191) q[0];
rz(-pi) q[1];
x q[1];
rz(0.78042729) q[2];
sx q[2];
rz(-1.6882511) q[2];
sx q[2];
rz(-0.800363) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.74976512) q[1];
sx q[1];
rz(-0.76376643) q[1];
sx q[1];
rz(0.78322021) q[1];
rz(-pi) q[2];
rz(1.901504) q[3];
sx q[3];
rz(-2.7484012) q[3];
sx q[3];
rz(1.1101983) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.41875276) q[2];
sx q[2];
rz(-0.80416983) q[2];
sx q[2];
rz(-1.0214825) q[2];
rz(0.57182765) q[3];
sx q[3];
rz(-0.56838667) q[3];
sx q[3];
rz(-0.9182601) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.21917139) q[0];
sx q[0];
rz(-2.2497441) q[0];
sx q[0];
rz(0.047155596) q[0];
rz(-1.4157408) q[1];
sx q[1];
rz(-1.9784617) q[1];
sx q[1];
rz(-0.39173752) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4120868) q[0];
sx q[0];
rz(-1.1177309) q[0];
sx q[0];
rz(1.5231251) q[0];
rz(1.91003) q[2];
sx q[2];
rz(-1.3922786) q[2];
sx q[2];
rz(-2.0925131) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.22931598) q[1];
sx q[1];
rz(-0.01688334) q[1];
sx q[1];
rz(-0.90516488) q[1];
x q[2];
rz(0.39288951) q[3];
sx q[3];
rz(-1.3706012) q[3];
sx q[3];
rz(1.6416711) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.2262816) q[2];
sx q[2];
rz(-0.87258029) q[2];
sx q[2];
rz(-0.20381168) q[2];
rz(-0.63621825) q[3];
sx q[3];
rz(-1.7812984) q[3];
sx q[3];
rz(-0.067559592) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5930475) q[0];
sx q[0];
rz(-0.021634463) q[0];
sx q[0];
rz(2.0245323) q[0];
rz(1.2720269) q[1];
sx q[1];
rz(-2.2950324) q[1];
sx q[1];
rz(-0.55714947) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2761189) q[0];
sx q[0];
rz(-2.0003002) q[0];
sx q[0];
rz(-0.19709023) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.093634886) q[2];
sx q[2];
rz(-1.4057004) q[2];
sx q[2];
rz(1.5993725) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.8073358) q[1];
sx q[1];
rz(-2.7149221) q[1];
sx q[1];
rz(-0.92324275) q[1];
rz(-pi) q[2];
rz(-2.689834) q[3];
sx q[3];
rz(-2.5521899) q[3];
sx q[3];
rz(1.9077993) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7316526) q[2];
sx q[2];
rz(-2.9059548) q[2];
sx q[2];
rz(1.635599) q[2];
rz(2.951494) q[3];
sx q[3];
rz(-0.59571576) q[3];
sx q[3];
rz(3.0715517) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6079123) q[0];
sx q[0];
rz(-0.30617014) q[0];
sx q[0];
rz(2.876907) q[0];
rz(0.39868042) q[1];
sx q[1];
rz(-2.61187) q[1];
sx q[1];
rz(2.9369489) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5335352) q[0];
sx q[0];
rz(-1.7711955) q[0];
sx q[0];
rz(2.1222369) q[0];
rz(-pi) q[1];
rz(-2.2024353) q[2];
sx q[2];
rz(-1.3654786) q[2];
sx q[2];
rz(0.0042631342) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.8890587) q[1];
sx q[1];
rz(-1.6315347) q[1];
sx q[1];
rz(0.82734682) q[1];
x q[2];
rz(-0.27575043) q[3];
sx q[3];
rz(-2.6480452) q[3];
sx q[3];
rz(-0.38024732) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.22145049) q[2];
sx q[2];
rz(-1.7998989) q[2];
sx q[2];
rz(-1.7500925) q[2];
rz(-0.57735389) q[3];
sx q[3];
rz(-2.5632863) q[3];
sx q[3];
rz(2.3277843) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.53265424) q[0];
sx q[0];
rz(-1.895099) q[0];
sx q[0];
rz(1.1005703) q[0];
rz(-2.7990923) q[1];
sx q[1];
rz(-0.96440146) q[1];
sx q[1];
rz(-1.0920116) q[1];
rz(-1.2630812) q[2];
sx q[2];
rz(-1.7379496) q[2];
sx q[2];
rz(-1.3676277) q[2];
rz(-0.58811515) q[3];
sx q[3];
rz(-0.93265231) q[3];
sx q[3];
rz(1.5476641) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
