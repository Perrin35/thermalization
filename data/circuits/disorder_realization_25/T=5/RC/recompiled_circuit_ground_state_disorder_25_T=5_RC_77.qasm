OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.8018262) q[0];
sx q[0];
rz(-2.5751994) q[0];
sx q[0];
rz(1.9422148) q[0];
rz(2.3102923) q[1];
sx q[1];
rz(-1.6142694) q[1];
sx q[1];
rz(-2.6418614) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.45699994) q[0];
sx q[0];
rz(-1.3641806) q[0];
sx q[0];
rz(2.3159136) q[0];
rz(1.209916) q[2];
sx q[2];
rz(-1.3832051) q[2];
sx q[2];
rz(2.9071992) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.7669406) q[1];
sx q[1];
rz(-2.3071676) q[1];
sx q[1];
rz(1.9689421) q[1];
rz(-pi) q[2];
rz(-1.5152895) q[3];
sx q[3];
rz(-0.71904564) q[3];
sx q[3];
rz(-2.2628502) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.1685593) q[2];
sx q[2];
rz(-1.7567822) q[2];
sx q[2];
rz(2.2200269) q[2];
rz(1.5714931) q[3];
sx q[3];
rz(-0.67155963) q[3];
sx q[3];
rz(-0.69935548) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8916931) q[0];
sx q[0];
rz(-2.1136916) q[0];
sx q[0];
rz(1.6722884) q[0];
rz(-1.6647313) q[1];
sx q[1];
rz(-1.3427443) q[1];
sx q[1];
rz(-0.5161759) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0196386) q[0];
sx q[0];
rz(-0.40651822) q[0];
sx q[0];
rz(-1.3234786) q[0];
rz(-1.9024114) q[2];
sx q[2];
rz(-1.8027163) q[2];
sx q[2];
rz(0.9882016) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.33715926) q[1];
sx q[1];
rz(-2.8808631) q[1];
sx q[1];
rz(0.27267898) q[1];
rz(-2.1352876) q[3];
sx q[3];
rz(-0.7745452) q[3];
sx q[3];
rz(-1.8902799) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.76253647) q[2];
sx q[2];
rz(-1.9409337) q[2];
sx q[2];
rz(0.29736796) q[2];
rz(-0.51658806) q[3];
sx q[3];
rz(-1.1597495) q[3];
sx q[3];
rz(2.5166701) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49279889) q[0];
sx q[0];
rz(-2.6060217) q[0];
sx q[0];
rz(0.18185644) q[0];
rz(-2.6975373) q[1];
sx q[1];
rz(-2.2477138) q[1];
sx q[1];
rz(-0.081238834) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59653283) q[0];
sx q[0];
rz(-1.9659916) q[0];
sx q[0];
rz(-1.6020653) q[0];
x q[1];
rz(2.8625254) q[2];
sx q[2];
rz(-2.2418602) q[2];
sx q[2];
rz(-2.8299675) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.013147203) q[1];
sx q[1];
rz(-1.1025635) q[1];
sx q[1];
rz(0.20856671) q[1];
rz(-1.0100097) q[3];
sx q[3];
rz(-1.6688235) q[3];
sx q[3];
rz(-1.3658226) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.2150725) q[2];
sx q[2];
rz(-0.36538616) q[2];
sx q[2];
rz(-0.60058769) q[2];
rz(1.8961204) q[3];
sx q[3];
rz(-0.54491091) q[3];
sx q[3];
rz(3.1020402) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
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
rz(-1.1323701) q[0];
sx q[0];
rz(-2.4339269) q[0];
sx q[0];
rz(-0.02455499) q[0];
rz(0.37697667) q[1];
sx q[1];
rz(-1.8507277) q[1];
sx q[1];
rz(-0.36184186) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.50411669) q[0];
sx q[0];
rz(-1.8808805) q[0];
sx q[0];
rz(-1.4092567) q[0];
rz(-pi) q[1];
rz(1.9458188) q[2];
sx q[2];
rz(-1.6437569) q[2];
sx q[2];
rz(-0.53158497) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.3256073) q[1];
sx q[1];
rz(-0.66130844) q[1];
sx q[1];
rz(2.8490752) q[1];
rz(0.41976069) q[3];
sx q[3];
rz(-1.1294522) q[3];
sx q[3];
rz(0.10548681) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.61818608) q[2];
sx q[2];
rz(-1.2647102) q[2];
sx q[2];
rz(0.82320172) q[2];
rz(0.19436714) q[3];
sx q[3];
rz(-2.7394962) q[3];
sx q[3];
rz(-0.91485867) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0219367) q[0];
sx q[0];
rz(-1.6031665) q[0];
sx q[0];
rz(0.62171474) q[0];
rz(2.3960528) q[1];
sx q[1];
rz(-2.3108683) q[1];
sx q[1];
rz(0.088836975) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6357036) q[0];
sx q[0];
rz(-0.14378366) q[0];
sx q[0];
rz(0.48322387) q[0];
rz(1.9000824) q[2];
sx q[2];
rz(-1.3111787) q[2];
sx q[2];
rz(3.064472) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.34804071) q[1];
sx q[1];
rz(-1.8489535) q[1];
sx q[1];
rz(-2.3955961) q[1];
x q[2];
rz(1.2337429) q[3];
sx q[3];
rz(-1.7255069) q[3];
sx q[3];
rz(-1.0708933) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.9125774) q[2];
sx q[2];
rz(-1.6385767) q[2];
sx q[2];
rz(-1.6055321) q[2];
rz(2.33365) q[3];
sx q[3];
rz(-0.84351051) q[3];
sx q[3];
rz(1.9690751) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.87869969) q[0];
sx q[0];
rz(-1.4233339) q[0];
sx q[0];
rz(0.95243564) q[0];
rz(1.7772504) q[1];
sx q[1];
rz(-1.9352501) q[1];
sx q[1];
rz(-0.84699026) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1035396) q[0];
sx q[0];
rz(-2.969739) q[0];
sx q[0];
rz(-0.58809963) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.32520099) q[2];
sx q[2];
rz(-2.3210827) q[2];
sx q[2];
rz(-2.6186117) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.0274732) q[1];
sx q[1];
rz(-2.4714258) q[1];
sx q[1];
rz(-0.7188188) q[1];
rz(-pi) q[2];
rz(3.0383864) q[3];
sx q[3];
rz(-1.5563193) q[3];
sx q[3];
rz(-0.38179427) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.18818894) q[2];
sx q[2];
rz(-1.0651411) q[2];
sx q[2];
rz(-2.7597001) q[2];
rz(2.0415107) q[3];
sx q[3];
rz(-0.81172687) q[3];
sx q[3];
rz(2.7303117) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5683658) q[0];
sx q[0];
rz(-1.8170284) q[0];
sx q[0];
rz(0.49760094) q[0];
rz(0.35500232) q[1];
sx q[1];
rz(-1.4811131) q[1];
sx q[1];
rz(-2.3752046) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3610182) q[0];
sx q[0];
rz(-2.0832591) q[0];
sx q[0];
rz(1.6105525) q[0];
x q[1];
rz(-2.1649579) q[2];
sx q[2];
rz(-2.4554002) q[2];
sx q[2];
rz(2.2399069) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.97624) q[1];
sx q[1];
rz(-1.2818977) q[1];
sx q[1];
rz(1.5084748) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3192407) q[3];
sx q[3];
rz(-2.2330604) q[3];
sx q[3];
rz(1.7098302) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.7327205) q[2];
sx q[2];
rz(-1.6072105) q[2];
sx q[2];
rz(1.6498148) q[2];
rz(-2.0626227) q[3];
sx q[3];
rz(-1.3078728) q[3];
sx q[3];
rz(0.61416793) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
sx q[3];
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
rz(0.33477467) q[0];
sx q[0];
rz(-0.65009999) q[0];
sx q[0];
rz(-2.7244869) q[0];
rz(3.0879703) q[1];
sx q[1];
rz(-1.6157179) q[1];
sx q[1];
rz(-2.0448304) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3552719) q[0];
sx q[0];
rz(-1.6910166) q[0];
sx q[0];
rz(-2.1646196) q[0];
x q[1];
rz(0.057889537) q[2];
sx q[2];
rz(-1.1107003) q[2];
sx q[2];
rz(1.636508) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.5973401) q[1];
sx q[1];
rz(-1.5637923) q[1];
sx q[1];
rz(0.16503741) q[1];
rz(-pi) q[2];
x q[2];
rz(0.18713672) q[3];
sx q[3];
rz(-0.68681648) q[3];
sx q[3];
rz(2.2568964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.9216448) q[2];
sx q[2];
rz(-1.767445) q[2];
sx q[2];
rz(2.5448223) q[2];
rz(0.88045949) q[3];
sx q[3];
rz(-1.7170649) q[3];
sx q[3];
rz(-2.5318291) q[3];
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
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.45744687) q[0];
sx q[0];
rz(-0.54231751) q[0];
sx q[0];
rz(-0.31998262) q[0];
rz(-0.34842247) q[1];
sx q[1];
rz(-0.44487822) q[1];
sx q[1];
rz(-2.6382823) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1419174) q[0];
sx q[0];
rz(-1.7192695) q[0];
sx q[0];
rz(-3.0485118) q[0];
rz(-pi) q[1];
rz(0.82292436) q[2];
sx q[2];
rz(-1.1852929) q[2];
sx q[2];
rz(-2.4679135) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.7855808) q[1];
sx q[1];
rz(-1.4921489) q[1];
sx q[1];
rz(-1.1376343) q[1];
rz(-pi) q[2];
x q[2];
rz(0.21280913) q[3];
sx q[3];
rz(-2.6104038) q[3];
sx q[3];
rz(2.2399898) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.2549501) q[2];
sx q[2];
rz(-0.4526259) q[2];
sx q[2];
rz(2.5060999) q[2];
rz(-1.6057711) q[3];
sx q[3];
rz(-1.7255892) q[3];
sx q[3];
rz(-3.0543069) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.6814293) q[0];
sx q[0];
rz(-0.59578139) q[0];
sx q[0];
rz(1.7728565) q[0];
rz(0.5303371) q[1];
sx q[1];
rz(-1.4884721) q[1];
sx q[1];
rz(-0.77290767) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3221098) q[0];
sx q[0];
rz(-2.3272996) q[0];
sx q[0];
rz(2.5952475) q[0];
rz(-pi) q[1];
rz(0.79276104) q[2];
sx q[2];
rz(-1.498641) q[2];
sx q[2];
rz(2.4506086) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.7557395) q[1];
sx q[1];
rz(-0.47702152) q[1];
sx q[1];
rz(-1.9128837) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.021455) q[3];
sx q[3];
rz(-1.8348331) q[3];
sx q[3];
rz(-2.9024189) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.7698001) q[2];
sx q[2];
rz(-2.1830406) q[2];
sx q[2];
rz(0.27251631) q[2];
rz(2.6767327) q[3];
sx q[3];
rz(-1.2010937) q[3];
sx q[3];
rz(-0.71729898) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40806017) q[0];
sx q[0];
rz(-1.536674) q[0];
sx q[0];
rz(1.6996171) q[0];
rz(1.5245262) q[1];
sx q[1];
rz(-0.99823141) q[1];
sx q[1];
rz(0.29481606) q[1];
rz(1.8275558) q[2];
sx q[2];
rz(-1.8770367) q[2];
sx q[2];
rz(-0.0035280037) q[2];
rz(0.22188998) q[3];
sx q[3];
rz(-2.0035012) q[3];
sx q[3];
rz(1.1825081) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
