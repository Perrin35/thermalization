OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.70513201) q[0];
sx q[0];
rz(6.8350514) q[0];
sx q[0];
rz(9.4466136) q[0];
rz(2.7472189) q[1];
sx q[1];
rz(-1.4596649) q[1];
sx q[1];
rz(2.9266761) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37786814) q[0];
sx q[0];
rz(-1.2963429) q[0];
sx q[0];
rz(-0.26835693) q[0];
x q[1];
rz(0.46918842) q[2];
sx q[2];
rz(-1.708963) q[2];
sx q[2];
rz(-2.5772622) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.6282181) q[1];
sx q[1];
rz(-2.3623423) q[1];
sx q[1];
rz(2.2378504) q[1];
rz(-pi) q[2];
rz(2.5932556) q[3];
sx q[3];
rz(-0.79474802) q[3];
sx q[3];
rz(0.37219513) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.73137838) q[2];
sx q[2];
rz(-1.6822858) q[2];
sx q[2];
rz(0.56420502) q[2];
rz(1.7764067) q[3];
sx q[3];
rz(-2.6919638) q[3];
sx q[3];
rz(-1.8723429) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0974225) q[0];
sx q[0];
rz(-1.928227) q[0];
sx q[0];
rz(0.92798293) q[0];
rz(1.1652975) q[1];
sx q[1];
rz(-1.5382643) q[1];
sx q[1];
rz(-2.2448418) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39820406) q[0];
sx q[0];
rz(-2.2687015) q[0];
sx q[0];
rz(-2.8845805) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4977658) q[2];
sx q[2];
rz(-1.739193) q[2];
sx q[2];
rz(-1.8254691) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.0838544) q[1];
sx q[1];
rz(-1.6946304) q[1];
sx q[1];
rz(-0.79353516) q[1];
x q[2];
rz(-2.3223022) q[3];
sx q[3];
rz(-1.2580401) q[3];
sx q[3];
rz(2.6889338) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.26560489) q[2];
sx q[2];
rz(-2.6066055) q[2];
sx q[2];
rz(1.0401475) q[2];
rz(1.6863719) q[3];
sx q[3];
rz(-1.8486332) q[3];
sx q[3];
rz(-0.35475981) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4002157) q[0];
sx q[0];
rz(-0.57755661) q[0];
sx q[0];
rz(1.0282015) q[0];
rz(-2.0630515) q[1];
sx q[1];
rz(-0.56285793) q[1];
sx q[1];
rz(2.7064586) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2745278) q[0];
sx q[0];
rz(-1.9221677) q[0];
sx q[0];
rz(1.4562796) q[0];
x q[1];
rz(2.8912796) q[2];
sx q[2];
rz(-1.4280983) q[2];
sx q[2];
rz(0.38562361) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.5907026) q[1];
sx q[1];
rz(-2.0520376) q[1];
sx q[1];
rz(-0.18443702) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9153254) q[3];
sx q[3];
rz(-0.83762729) q[3];
sx q[3];
rz(-1.6813577) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.6083287) q[2];
sx q[2];
rz(-1.8412795) q[2];
sx q[2];
rz(2.8386774) q[2];
rz(-1.8164002) q[3];
sx q[3];
rz(-1.9830827) q[3];
sx q[3];
rz(-3.0505676) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9451697) q[0];
sx q[0];
rz(-1.4099932) q[0];
sx q[0];
rz(-2.2241425) q[0];
rz(2.4687185) q[1];
sx q[1];
rz(-2.0560975) q[1];
sx q[1];
rz(-0.26487574) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.53752758) q[0];
sx q[0];
rz(-0.91188216) q[0];
sx q[0];
rz(-1.5327246) q[0];
x q[1];
rz(-2.3189544) q[2];
sx q[2];
rz(-2.7441141) q[2];
sx q[2];
rz(-0.17695225) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.5502364) q[1];
sx q[1];
rz(-1.9758421) q[1];
sx q[1];
rz(0.77781271) q[1];
rz(-pi) q[2];
rz(2.3817252) q[3];
sx q[3];
rz(-0.93899512) q[3];
sx q[3];
rz(-0.85216537) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.778487) q[2];
sx q[2];
rz(-2.6532756) q[2];
sx q[2];
rz(-1.5765566) q[2];
rz(-2.1145084) q[3];
sx q[3];
rz(-2.4001207) q[3];
sx q[3];
rz(1.1013793) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89001369) q[0];
sx q[0];
rz(-0.13680923) q[0];
sx q[0];
rz(-2.662861) q[0];
rz(2.1084673) q[1];
sx q[1];
rz(-0.9712351) q[1];
sx q[1];
rz(-2.1889401) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6081776) q[0];
sx q[0];
rz(-2.1153643) q[0];
sx q[0];
rz(-1.4334701) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7258337) q[2];
sx q[2];
rz(-1.2649049) q[2];
sx q[2];
rz(-1.8005467) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.488727) q[1];
sx q[1];
rz(-1.1796724) q[1];
sx q[1];
rz(0.25686849) q[1];
rz(-pi) q[2];
rz(1.5289375) q[3];
sx q[3];
rz(-1.0930982) q[3];
sx q[3];
rz(-0.65934138) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.4218563) q[2];
sx q[2];
rz(-0.36281261) q[2];
sx q[2];
rz(-2.6110113) q[2];
rz(1.7355708) q[3];
sx q[3];
rz(-2.0134182) q[3];
sx q[3];
rz(0.83166844) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.56753165) q[0];
sx q[0];
rz(-1.644779) q[0];
sx q[0];
rz(1.6249599) q[0];
rz(-1.8364871) q[1];
sx q[1];
rz(-1.790698) q[1];
sx q[1];
rz(-2.9690202) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1369551) q[0];
sx q[0];
rz(-1.1543659) q[0];
sx q[0];
rz(-0.01491551) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0270945) q[2];
sx q[2];
rz(-1.8048394) q[2];
sx q[2];
rz(-2.8161088) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-3.1291618) q[1];
sx q[1];
rz(-0.42814246) q[1];
sx q[1];
rz(1.3151602) q[1];
rz(0.31130143) q[3];
sx q[3];
rz(-0.65400306) q[3];
sx q[3];
rz(-1.9037387) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.83795786) q[2];
sx q[2];
rz(-1.4540318) q[2];
sx q[2];
rz(-2.0149569) q[2];
rz(2.3593694) q[3];
sx q[3];
rz(-1.2354847) q[3];
sx q[3];
rz(1.8036028) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.85309) q[0];
sx q[0];
rz(-2.8350916) q[0];
sx q[0];
rz(-0.66147584) q[0];
rz(-0.96039564) q[1];
sx q[1];
rz(-1.7405225) q[1];
sx q[1];
rz(-0.75659928) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.21515439) q[0];
sx q[0];
rz(-1.4190136) q[0];
sx q[0];
rz(0.093869165) q[0];
x q[1];
rz(2.7034764) q[2];
sx q[2];
rz(-0.44898673) q[2];
sx q[2];
rz(-0.24030906) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.8826897) q[1];
sx q[1];
rz(-1.3414624) q[1];
sx q[1];
rz(-2.980466) q[1];
rz(-pi) q[2];
rz(-1.6301304) q[3];
sx q[3];
rz(-1.3795128) q[3];
sx q[3];
rz(0.64185601) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.0044272) q[2];
sx q[2];
rz(-0.1903154) q[2];
sx q[2];
rz(-2.9471617) q[2];
rz(2.2284609) q[3];
sx q[3];
rz(-1.7539932) q[3];
sx q[3];
rz(2.156179) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-0.59654355) q[0];
sx q[0];
rz(-0.61674917) q[0];
sx q[0];
rz(0.066666691) q[0];
rz(2.8170259) q[1];
sx q[1];
rz(-1.5044183) q[1];
sx q[1];
rz(0.98888046) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.37799997) q[0];
sx q[0];
rz(-2.6110296) q[0];
sx q[0];
rz(-2.2647122) q[0];
x q[1];
rz(3.0481911) q[2];
sx q[2];
rz(-1.5548692) q[2];
sx q[2];
rz(-1.8708558) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.40438548) q[1];
sx q[1];
rz(-0.42006856) q[1];
sx q[1];
rz(-1.406548) q[1];
x q[2];
rz(0.75108053) q[3];
sx q[3];
rz(-2.4820648) q[3];
sx q[3];
rz(-2.2758323) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.4618335) q[2];
sx q[2];
rz(-2.2448848) q[2];
sx q[2];
rz(0.40763339) q[2];
rz(-0.76861012) q[3];
sx q[3];
rz(-1.3137484) q[3];
sx q[3];
rz(2.2176946) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.75893629) q[0];
sx q[0];
rz(-1.8396682) q[0];
sx q[0];
rz(-2.5323903) q[0];
rz(0.095104782) q[1];
sx q[1];
rz(-1.2520049) q[1];
sx q[1];
rz(-0.87337714) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.72520032) q[0];
sx q[0];
rz(-1.7500449) q[0];
sx q[0];
rz(0.068530131) q[0];
x q[1];
rz(2.3960605) q[2];
sx q[2];
rz(-0.19956707) q[2];
sx q[2];
rz(1.0996639) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.59233353) q[1];
sx q[1];
rz(-0.98456406) q[1];
sx q[1];
rz(1.9418282) q[1];
x q[2];
rz(-2.7018413) q[3];
sx q[3];
rz(-1.0381178) q[3];
sx q[3];
rz(0.38910481) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.0218899) q[2];
sx q[2];
rz(-2.3535574) q[2];
sx q[2];
rz(-2.4592887) q[2];
rz(-0.37426379) q[3];
sx q[3];
rz(-1.5687317) q[3];
sx q[3];
rz(-2.9746829) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69797126) q[0];
sx q[0];
rz(-1.781783) q[0];
sx q[0];
rz(-0.95296729) q[0];
rz(-0.8264181) q[1];
sx q[1];
rz(-0.73917878) q[1];
sx q[1];
rz(-1.3964765) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.29165927) q[0];
sx q[0];
rz(-0.33432654) q[0];
sx q[0];
rz(0.26419421) q[0];
x q[1];
rz(-2.4834677) q[2];
sx q[2];
rz(-1.3326367) q[2];
sx q[2];
rz(0.70891526) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.64618387) q[1];
sx q[1];
rz(-2.7135661) q[1];
sx q[1];
rz(2.3133548) q[1];
rz(-pi) q[2];
rz(0.20874899) q[3];
sx q[3];
rz(-0.26780805) q[3];
sx q[3];
rz(2.1684614) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.46618) q[2];
sx q[2];
rz(-2.7853577) q[2];
sx q[2];
rz(-2.9818025) q[2];
rz(2.8397078) q[3];
sx q[3];
rz(-0.92697898) q[3];
sx q[3];
rz(-2.8543499) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.044534279) q[0];
sx q[0];
rz(-2.4659768) q[0];
sx q[0];
rz(1.5855047) q[0];
rz(0.13327577) q[1];
sx q[1];
rz(-1.6242846) q[1];
sx q[1];
rz(-0.12856738) q[1];
rz(-2.9524654) q[2];
sx q[2];
rz(-0.94869877) q[2];
sx q[2];
rz(-1.000065) q[2];
rz(-1.521048) q[3];
sx q[3];
rz(-1.0377025) q[3];
sx q[3];
rz(2.514537) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
