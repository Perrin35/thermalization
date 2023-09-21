OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.70541731) q[0];
sx q[0];
rz(-2.5751312) q[0];
sx q[0];
rz(-0.17106549) q[0];
rz(1.45362) q[1];
sx q[1];
rz(-0.34314081) q[1];
sx q[1];
rz(1.8106102) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1216461) q[0];
sx q[0];
rz(-2.2197147) q[0];
sx q[0];
rz(-1.5050423) q[0];
rz(-pi) q[1];
rz(1.5491484) q[2];
sx q[2];
rz(-2.0865555) q[2];
sx q[2];
rz(-1.6025008) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.0806182) q[1];
sx q[1];
rz(-1.3197101) q[1];
sx q[1];
rz(-0.017647839) q[1];
rz(-0.64294502) q[3];
sx q[3];
rz(-1.391198) q[3];
sx q[3];
rz(2.5760866) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.3657637) q[2];
sx q[2];
rz(-0.77314955) q[2];
sx q[2];
rz(1.3295056) q[2];
rz(-1.4154411) q[3];
sx q[3];
rz(-1.5664145) q[3];
sx q[3];
rz(-2.9392488) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1502894) q[0];
sx q[0];
rz(-2.5973899) q[0];
sx q[0];
rz(1.6888899) q[0];
rz(2.9406722) q[1];
sx q[1];
rz(-2.0566437) q[1];
sx q[1];
rz(-0.30028775) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.18314221) q[0];
sx q[0];
rz(-0.97226769) q[0];
sx q[0];
rz(1.4248225) q[0];
rz(-2.5710201) q[2];
sx q[2];
rz(-1.8405481) q[2];
sx q[2];
rz(-2.3938993) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.092204658) q[1];
sx q[1];
rz(-2.7687679) q[1];
sx q[1];
rz(0.91918175) q[1];
rz(-pi) q[2];
rz(0.16365666) q[3];
sx q[3];
rz(-1.7573342) q[3];
sx q[3];
rz(1.2638307) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.369027) q[2];
sx q[2];
rz(-2.7894661) q[2];
sx q[2];
rz(-2.1035813) q[2];
rz(2.299262) q[3];
sx q[3];
rz(-1.4146283) q[3];
sx q[3];
rz(-2.7712908) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5791941) q[0];
sx q[0];
rz(-2.6694522) q[0];
sx q[0];
rz(2.753479) q[0];
rz(-0.072470486) q[1];
sx q[1];
rz(-1.7150755) q[1];
sx q[1];
rz(-2.8252576) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4916723) q[0];
sx q[0];
rz(-0.24501093) q[0];
sx q[0];
rz(-1.578376) q[0];
rz(0.57245589) q[2];
sx q[2];
rz(-1.2081895) q[2];
sx q[2];
rz(2.2628757) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.0132732) q[1];
sx q[1];
rz(-1.1228021) q[1];
sx q[1];
rz(-2.0112579) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9100788) q[3];
sx q[3];
rz(-1.2959891) q[3];
sx q[3];
rz(0.85698444) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.0324273) q[2];
sx q[2];
rz(-0.18523231) q[2];
sx q[2];
rz(-2.9807828) q[2];
rz(-0.12604776) q[3];
sx q[3];
rz(-1.7536609) q[3];
sx q[3];
rz(2.1179312) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-1.9008824) q[0];
sx q[0];
rz(-0.70403376) q[0];
sx q[0];
rz(2.3098992) q[0];
rz(-1.7968934) q[1];
sx q[1];
rz(-2.1422377) q[1];
sx q[1];
rz(1.5136738) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9261525) q[0];
sx q[0];
rz(-2.2427796) q[0];
sx q[0];
rz(-1.6295022) q[0];
rz(-pi) q[1];
rz(-3.0604826) q[2];
sx q[2];
rz(-1.0081648) q[2];
sx q[2];
rz(-0.62603355) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.1308243) q[1];
sx q[1];
rz(-1.1225892) q[1];
sx q[1];
rz(3.1229449) q[1];
x q[2];
rz(-2.4098445) q[3];
sx q[3];
rz(-1.2536612) q[3];
sx q[3];
rz(-1.5665311) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.7164798) q[2];
sx q[2];
rz(-1.0093062) q[2];
sx q[2];
rz(-1.224219) q[2];
rz(-0.41695693) q[3];
sx q[3];
rz(-1.0427534) q[3];
sx q[3];
rz(-0.52880374) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.058218) q[0];
sx q[0];
rz(-0.26979065) q[0];
sx q[0];
rz(-1.8378687) q[0];
rz(-2.6858792) q[1];
sx q[1];
rz(-0.2625176) q[1];
sx q[1];
rz(-3.0900893) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7464375) q[0];
sx q[0];
rz(-1.6704541) q[0];
sx q[0];
rz(1.6874203) q[0];
x q[1];
rz(0.15437834) q[2];
sx q[2];
rz(-1.8533684) q[2];
sx q[2];
rz(2.4992361) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.67700779) q[1];
sx q[1];
rz(-2.0160463) q[1];
sx q[1];
rz(0.17226179) q[1];
rz(2.4078835) q[3];
sx q[3];
rz(-1.6624311) q[3];
sx q[3];
rz(1.1250145) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.45433989) q[2];
sx q[2];
rz(-1.7752825) q[2];
sx q[2];
rz(-2.5435737) q[2];
rz(0.55001843) q[3];
sx q[3];
rz(-0.23967448) q[3];
sx q[3];
rz(-2.0050744) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-2.0710058) q[0];
sx q[0];
rz(-0.07659176) q[0];
sx q[0];
rz(2.9396074) q[0];
rz(-2.1760991) q[1];
sx q[1];
rz(-1.0243203) q[1];
sx q[1];
rz(0.083267033) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.14295386) q[0];
sx q[0];
rz(-2.8790701) q[0];
sx q[0];
rz(-2.4130164) q[0];
rz(-2.8667738) q[2];
sx q[2];
rz(-1.8289369) q[2];
sx q[2];
rz(-1.2061662) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.1042085) q[1];
sx q[1];
rz(-1.6777778) q[1];
sx q[1];
rz(1.2030829) q[1];
rz(1.2721328) q[3];
sx q[3];
rz(-2.1080058) q[3];
sx q[3];
rz(-0.99029535) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-3.0319556) q[2];
sx q[2];
rz(-1.6836616) q[2];
sx q[2];
rz(1.7588245) q[2];
rz(-2.8921195) q[3];
sx q[3];
rz(-1.8981257) q[3];
sx q[3];
rz(1.680254) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4270585) q[0];
sx q[0];
rz(-2.336851) q[0];
sx q[0];
rz(-2.2447341) q[0];
rz(-0.25009051) q[1];
sx q[1];
rz(-2.9607594) q[1];
sx q[1];
rz(2.0239963) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5898949) q[0];
sx q[0];
rz(-2.4015744) q[0];
sx q[0];
rz(-0.19703534) q[0];
rz(-pi) q[1];
rz(0.859185) q[2];
sx q[2];
rz(-1.4118328) q[2];
sx q[2];
rz(-2.8189916) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.1114137) q[1];
sx q[1];
rz(-1.3707268) q[1];
sx q[1];
rz(2.5740037) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8459122) q[3];
sx q[3];
rz(-1.0529622) q[3];
sx q[3];
rz(-2.8003611) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.08427944) q[2];
sx q[2];
rz(-1.7417615) q[2];
sx q[2];
rz(2.9220667) q[2];
rz(0.38671842) q[3];
sx q[3];
rz(-0.76824776) q[3];
sx q[3];
rz(2.1462671) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.052208386) q[0];
sx q[0];
rz(-1.5445671) q[0];
sx q[0];
rz(0.21729939) q[0];
rz(-0.030933881) q[1];
sx q[1];
rz(-2.5035281) q[1];
sx q[1];
rz(-2.4826179) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0156408) q[0];
sx q[0];
rz(-0.77227393) q[0];
sx q[0];
rz(3.1053567) q[0];
rz(-pi) q[1];
rz(-2.4202004) q[2];
sx q[2];
rz(-1.0288236) q[2];
sx q[2];
rz(1.1162356) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.79170376) q[1];
sx q[1];
rz(-0.83194299) q[1];
sx q[1];
rz(-1.8658584) q[1];
rz(-pi) q[2];
x q[2];
rz(0.92774763) q[3];
sx q[3];
rz(-1.5156931) q[3];
sx q[3];
rz(-0.4004713) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.4387681) q[2];
sx q[2];
rz(-1.3805026) q[2];
sx q[2];
rz(-2.8213815) q[2];
rz(-1.5252339) q[3];
sx q[3];
rz(-1.9459008) q[3];
sx q[3];
rz(-1.2493791) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3808688) q[0];
sx q[0];
rz(-0.18146935) q[0];
sx q[0];
rz(-2.4587801) q[0];
rz(1.0702417) q[1];
sx q[1];
rz(-1.8990592) q[1];
sx q[1];
rz(-1.9314996) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.72737981) q[0];
sx q[0];
rz(-2.7276037) q[0];
sx q[0];
rz(2.1075641) q[0];
rz(-pi) q[1];
rz(-1.2145043) q[2];
sx q[2];
rz(-2.5062343) q[2];
sx q[2];
rz(2.2941342) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.853211) q[1];
sx q[1];
rz(-2.2182811) q[1];
sx q[1];
rz(-2.7700469) q[1];
rz(-pi) q[2];
rz(2.6415778) q[3];
sx q[3];
rz(-0.68272907) q[3];
sx q[3];
rz(0.4993712) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.5658297) q[2];
sx q[2];
rz(-2.1255707) q[2];
sx q[2];
rz(2.5749717) q[2];
rz(-2.2120655) q[3];
sx q[3];
rz(-2.1598787) q[3];
sx q[3];
rz(-1.7738336) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.728445) q[0];
sx q[0];
rz(-2.8840273) q[0];
sx q[0];
rz(1.7096827) q[0];
rz(2.6152949) q[1];
sx q[1];
rz(-0.52733517) q[1];
sx q[1];
rz(0.79968232) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.66946) q[0];
sx q[0];
rz(-1.2233943) q[0];
sx q[0];
rz(-0.090263788) q[0];
x q[1];
rz(2.9808447) q[2];
sx q[2];
rz(-1.6511917) q[2];
sx q[2];
rz(-0.71842566) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.666115) q[1];
sx q[1];
rz(-1.0628504) q[1];
sx q[1];
rz(0.39132262) q[1];
rz(-pi) q[2];
rz(1.4950072) q[3];
sx q[3];
rz(-1.7160551) q[3];
sx q[3];
rz(2.4491572) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.2852823) q[2];
sx q[2];
rz(-0.44390634) q[2];
sx q[2];
rz(-0.69520673) q[2];
rz(-0.55784145) q[3];
sx q[3];
rz(-1.7772243) q[3];
sx q[3];
rz(-0.69303304) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0859062) q[0];
sx q[0];
rz(-1.9491371) q[0];
sx q[0];
rz(0.51167713) q[0];
rz(1.862539) q[1];
sx q[1];
rz(-2.3976354) q[1];
sx q[1];
rz(2.4639113) q[1];
rz(2.8056801) q[2];
sx q[2];
rz(-1.5716142) q[2];
sx q[2];
rz(-1.4128528) q[2];
rz(1.4064895) q[3];
sx q[3];
rz(-0.77948979) q[3];
sx q[3];
rz(-0.43850552) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];