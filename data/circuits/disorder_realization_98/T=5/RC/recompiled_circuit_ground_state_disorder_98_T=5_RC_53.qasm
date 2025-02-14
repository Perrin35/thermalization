OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.40795657) q[0];
sx q[0];
rz(-2.9562558) q[0];
sx q[0];
rz(1.6428525) q[0];
rz(-0.19417956) q[1];
sx q[1];
rz(-0.3436389) q[1];
sx q[1];
rz(-0.37771168) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2687338) q[0];
sx q[0];
rz(-1.3493363) q[0];
sx q[0];
rz(-0.24747699) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.333856) q[2];
sx q[2];
rz(-2.5662072) q[2];
sx q[2];
rz(0.45562109) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.9895674) q[1];
sx q[1];
rz(-1.5366239) q[1];
sx q[1];
rz(-0.65712838) q[1];
rz(-pi) q[2];
x q[2];
rz(1.463428) q[3];
sx q[3];
rz(-2.8753548) q[3];
sx q[3];
rz(-0.050198089) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.3898042) q[2];
sx q[2];
rz(-1.3499539) q[2];
sx q[2];
rz(-2.8765615) q[2];
rz(0.42936471) q[3];
sx q[3];
rz(-1.8931484) q[3];
sx q[3];
rz(-0.00092367729) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65381831) q[0];
sx q[0];
rz(-2.9980897) q[0];
sx q[0];
rz(-1.8408467) q[0];
rz(3.0744413) q[1];
sx q[1];
rz(-2.2323699) q[1];
sx q[1];
rz(0.86004177) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81102449) q[0];
sx q[0];
rz(-2.1386792) q[0];
sx q[0];
rz(-2.3080243) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3798736) q[2];
sx q[2];
rz(-0.8702308) q[2];
sx q[2];
rz(-0.127244) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.94446429) q[1];
sx q[1];
rz(-1.010506) q[1];
sx q[1];
rz(-0.48546882) q[1];
rz(0.18386358) q[3];
sx q[3];
rz(-1.2753092) q[3];
sx q[3];
rz(-1.9744524) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.63090515) q[2];
sx q[2];
rz(-0.61669934) q[2];
sx q[2];
rz(2.5751233) q[2];
rz(0.72930068) q[3];
sx q[3];
rz(-2.2972378) q[3];
sx q[3];
rz(-2.9384889) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1931964) q[0];
sx q[0];
rz(-1.0862792) q[0];
sx q[0];
rz(-2.7753944) q[0];
rz(-1.5858448) q[1];
sx q[1];
rz(-1.0428753) q[1];
sx q[1];
rz(3.0911456) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9589466) q[0];
sx q[0];
rz(-3.052127) q[0];
sx q[0];
rz(-2.0561809) q[0];
rz(0.68784662) q[2];
sx q[2];
rz(-1.1769933) q[2];
sx q[2];
rz(1.1325768) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.0512426) q[1];
sx q[1];
rz(-1.7900677) q[1];
sx q[1];
rz(-1.7419113) q[1];
rz(-pi) q[2];
rz(1.8837711) q[3];
sx q[3];
rz(-2.69776) q[3];
sx q[3];
rz(2.7601506) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.067165) q[2];
sx q[2];
rz(-1.0600435) q[2];
sx q[2];
rz(1.3107497) q[2];
rz(-1.358076) q[3];
sx q[3];
rz(-2.570593) q[3];
sx q[3];
rz(1.3091492) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8686789) q[0];
sx q[0];
rz(-0.22628117) q[0];
sx q[0];
rz(-1.7568461) q[0];
rz(1.9028496) q[1];
sx q[1];
rz(-1.9107995) q[1];
sx q[1];
rz(-1.967427) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.56856241) q[0];
sx q[0];
rz(-1.996604) q[0];
sx q[0];
rz(2.9338475) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.27710813) q[2];
sx q[2];
rz(-1.8113675) q[2];
sx q[2];
rz(-1.6024557) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.2432855) q[1];
sx q[1];
rz(-1.8970058) q[1];
sx q[1];
rz(-1.0350569) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0405868) q[3];
sx q[3];
rz(-2.3755666) q[3];
sx q[3];
rz(-1.211973) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.7063286) q[2];
sx q[2];
rz(-1.3202983) q[2];
sx q[2];
rz(3.1026133) q[2];
rz(0.6428166) q[3];
sx q[3];
rz(-1.0182074) q[3];
sx q[3];
rz(0.62079287) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3103127) q[0];
sx q[0];
rz(-1.0142925) q[0];
sx q[0];
rz(-3.1211299) q[0];
rz(2.9488355) q[1];
sx q[1];
rz(-0.61734504) q[1];
sx q[1];
rz(-1.9761168) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.504731) q[0];
sx q[0];
rz(-0.76380542) q[0];
sx q[0];
rz(2.5434407) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4679409) q[2];
sx q[2];
rz(-1.21017) q[2];
sx q[2];
rz(-1.1471105) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.9679035) q[1];
sx q[1];
rz(-0.5040796) q[1];
sx q[1];
rz(-1.5289115) q[1];
rz(-pi) q[2];
rz(1.9475157) q[3];
sx q[3];
rz(-2.93749) q[3];
sx q[3];
rz(-2.3538176) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.3141025) q[2];
sx q[2];
rz(-2.1496488) q[2];
sx q[2];
rz(-2.6143383) q[2];
rz(0.54316795) q[3];
sx q[3];
rz(-2.4542464) q[3];
sx q[3];
rz(-0.34272042) q[3];
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
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0602033) q[0];
sx q[0];
rz(-2.9855766) q[0];
sx q[0];
rz(-2.7348837) q[0];
rz(-1.1609062) q[1];
sx q[1];
rz(-0.91861594) q[1];
sx q[1];
rz(2.1312174) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5766616) q[0];
sx q[0];
rz(-2.1157678) q[0];
sx q[0];
rz(0.14133639) q[0];
rz(-pi) q[1];
x q[1];
rz(0.40596227) q[2];
sx q[2];
rz(-2.7926707) q[2];
sx q[2];
rz(-2.5418848) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.6076489) q[1];
sx q[1];
rz(-1.9393117) q[1];
sx q[1];
rz(1.7979969) q[1];
rz(0.80548894) q[3];
sx q[3];
rz(-1.7542766) q[3];
sx q[3];
rz(3.0940521) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.0747718) q[2];
sx q[2];
rz(-1.8596884) q[2];
sx q[2];
rz(-2.1066966) q[2];
rz(1.3567989) q[3];
sx q[3];
rz(-2.5467338) q[3];
sx q[3];
rz(2.518173) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9676301) q[0];
sx q[0];
rz(-1.6169463) q[0];
sx q[0];
rz(-0.3279283) q[0];
rz(-0.11218849) q[1];
sx q[1];
rz(-1.9393549) q[1];
sx q[1];
rz(-2.1690538) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2043861) q[0];
sx q[0];
rz(-0.68336672) q[0];
sx q[0];
rz(2.1478189) q[0];
rz(-0.99268408) q[2];
sx q[2];
rz(-2.6277866) q[2];
sx q[2];
rz(0.93573278) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.33148872) q[1];
sx q[1];
rz(-1.5520025) q[1];
sx q[1];
rz(2.5522425) q[1];
rz(1.2054382) q[3];
sx q[3];
rz(-0.97627538) q[3];
sx q[3];
rz(-0.48156092) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.605725) q[2];
sx q[2];
rz(-0.87656993) q[2];
sx q[2];
rz(-2.4884339) q[2];
rz(0.43290916) q[3];
sx q[3];
rz(-2.631729) q[3];
sx q[3];
rz(-2.9292817) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9369649) q[0];
sx q[0];
rz(-2.300394) q[0];
sx q[0];
rz(-2.9168108) q[0];
rz(2.3460491) q[1];
sx q[1];
rz(-1.3139775) q[1];
sx q[1];
rz(-2.0733817) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.58086181) q[0];
sx q[0];
rz(-2.1819072) q[0];
sx q[0];
rz(-0.9471425) q[0];
rz(-pi) q[1];
rz(-0.96282962) q[2];
sx q[2];
rz(-0.89917937) q[2];
sx q[2];
rz(1.9266018) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.873899) q[1];
sx q[1];
rz(-1.5963449) q[1];
sx q[1];
rz(-2.4539095) q[1];
rz(-pi) q[2];
rz(-1.3146888) q[3];
sx q[3];
rz(-0.22897069) q[3];
sx q[3];
rz(3.1310905) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.75605679) q[2];
sx q[2];
rz(-1.7588561) q[2];
sx q[2];
rz(0.63560152) q[2];
rz(0.33291891) q[3];
sx q[3];
rz(-2.1300485) q[3];
sx q[3];
rz(-0.46271589) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3391089) q[0];
sx q[0];
rz(-0.026263069) q[0];
sx q[0];
rz(0.12301692) q[0];
rz(1.3975551) q[1];
sx q[1];
rz(-1.3879644) q[1];
sx q[1];
rz(0.77973286) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5194979) q[0];
sx q[0];
rz(-0.72185282) q[0];
sx q[0];
rz(1.9869542) q[0];
rz(0.17895584) q[2];
sx q[2];
rz(-0.6935941) q[2];
sx q[2];
rz(-3.0998942) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.8416948) q[1];
sx q[1];
rz(-1.3720241) q[1];
sx q[1];
rz(1.9511306) q[1];
x q[2];
rz(-0.10752435) q[3];
sx q[3];
rz(-1.4443099) q[3];
sx q[3];
rz(-1.7582939) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.4745549) q[2];
sx q[2];
rz(-1.8755308) q[2];
sx q[2];
rz(-3.0511268) q[2];
rz(-0.41680923) q[3];
sx q[3];
rz(-0.79206812) q[3];
sx q[3];
rz(2.1478103) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4823293) q[0];
sx q[0];
rz(-1.3220795) q[0];
sx q[0];
rz(0.089476712) q[0];
rz(-0.80884519) q[1];
sx q[1];
rz(-2.6409812) q[1];
sx q[1];
rz(2.083875) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.63454483) q[0];
sx q[0];
rz(-1.9814361) q[0];
sx q[0];
rz(-1.9837186) q[0];
rz(3.1359768) q[2];
sx q[2];
rz(-2.1227909) q[2];
sx q[2];
rz(-2.698026) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.6597418) q[1];
sx q[1];
rz(-0.54164499) q[1];
sx q[1];
rz(-0.84195824) q[1];
rz(3.1154409) q[3];
sx q[3];
rz(-0.75687486) q[3];
sx q[3];
rz(2.7084086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.7384501) q[2];
sx q[2];
rz(-0.53692997) q[2];
sx q[2];
rz(0.37330791) q[2];
rz(2.8849854) q[3];
sx q[3];
rz(-1.4213057) q[3];
sx q[3];
rz(0.074450113) q[3];
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
rz(-pi/2) q[3];
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
rz(-1.2557209) q[0];
sx q[0];
rz(-1.5626361) q[0];
sx q[0];
rz(1.7241021) q[0];
rz(-1.4295084) q[1];
sx q[1];
rz(-2.7919339) q[1];
sx q[1];
rz(-2.3333593) q[1];
rz(-0.019349426) q[2];
sx q[2];
rz(-2.6826998) q[2];
sx q[2];
rz(-1.5738375) q[2];
rz(3.1028845) q[3];
sx q[3];
rz(-0.9699655) q[3];
sx q[3];
rz(-1.7326596) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
