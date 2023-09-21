OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.36800185) q[0];
sx q[0];
rz(-0.79080963) q[0];
sx q[0];
rz(-2.8074582) q[0];
rz(-0.45733991) q[1];
sx q[1];
rz(5.338905) q[1];
sx q[1];
rz(10.64325) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.89361184) q[0];
sx q[0];
rz(-2.3900744) q[0];
sx q[0];
rz(0.052608842) q[0];
x q[1];
rz(-2.676744) q[2];
sx q[2];
rz(-1.5343622) q[2];
sx q[2];
rz(-2.6236617) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.4301181) q[1];
sx q[1];
rz(-2.4415486) q[1];
sx q[1];
rz(-2.3858566) q[1];
rz(2.9524607) q[3];
sx q[3];
rz(-1.3526219) q[3];
sx q[3];
rz(1.8060341) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.5228287) q[2];
sx q[2];
rz(-0.4814119) q[2];
sx q[2];
rz(-2.5640326) q[2];
rz(-1.9918359) q[3];
sx q[3];
rz(-1.3883608) q[3];
sx q[3];
rz(-2.4770588) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1537271) q[0];
sx q[0];
rz(-0.58652121) q[0];
sx q[0];
rz(0.38744774) q[0];
rz(-0.93915835) q[1];
sx q[1];
rz(-0.99717957) q[1];
sx q[1];
rz(-1.4025677) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99968602) q[0];
sx q[0];
rz(-1.5667331) q[0];
sx q[0];
rz(3.1121029) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8005373) q[2];
sx q[2];
rz(-2.555072) q[2];
sx q[2];
rz(1.7797433) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.2991997) q[1];
sx q[1];
rz(-2.0626915) q[1];
sx q[1];
rz(2.030034) q[1];
rz(-pi) q[2];
rz(3.0121147) q[3];
sx q[3];
rz(-2.5054512) q[3];
sx q[3];
rz(0.57487088) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.42276057) q[2];
sx q[2];
rz(-1.3964802) q[2];
sx q[2];
rz(2.823901) q[2];
rz(-0.20673949) q[3];
sx q[3];
rz(-0.59967774) q[3];
sx q[3];
rz(2.3247705) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7725672) q[0];
sx q[0];
rz(-1.439753) q[0];
sx q[0];
rz(1.4136219) q[0];
rz(-2.6638022) q[1];
sx q[1];
rz(-1.7910035) q[1];
sx q[1];
rz(2.7405222) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0553592) q[0];
sx q[0];
rz(-1.8952574) q[0];
sx q[0];
rz(0.1089036) q[0];
rz(2.4972649) q[2];
sx q[2];
rz(-1.7228848) q[2];
sx q[2];
rz(-0.49567859) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.229278) q[1];
sx q[1];
rz(-0.73017987) q[1];
sx q[1];
rz(-3.1320523) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9455455) q[3];
sx q[3];
rz(-1.9571597) q[3];
sx q[3];
rz(0.26323174) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.59427375) q[2];
sx q[2];
rz(-1.5218364) q[2];
sx q[2];
rz(-2.5857914) q[2];
rz(-2.1650971) q[3];
sx q[3];
rz(-2.591811) q[3];
sx q[3];
rz(0.78021375) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.34898409) q[0];
sx q[0];
rz(-1.5190834) q[0];
sx q[0];
rz(-1.4439616) q[0];
rz(1.5199039) q[1];
sx q[1];
rz(-2.4855721) q[1];
sx q[1];
rz(2.8881883) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69617535) q[0];
sx q[0];
rz(-1.4679969) q[0];
sx q[0];
rz(1.0284543) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5379982) q[2];
sx q[2];
rz(-2.6555853) q[2];
sx q[2];
rz(2.3151468) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.036382347) q[1];
sx q[1];
rz(-1.4745592) q[1];
sx q[1];
rz(0.63016816) q[1];
x q[2];
rz(-1.4686618) q[3];
sx q[3];
rz(-0.92389744) q[3];
sx q[3];
rz(-1.1598066) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.110934) q[2];
sx q[2];
rz(-1.3867644) q[2];
sx q[2];
rz(-0.78732642) q[2];
rz(0.91286719) q[3];
sx q[3];
rz(-2.392231) q[3];
sx q[3];
rz(2.1319938) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.426429) q[0];
sx q[0];
rz(-0.63868317) q[0];
sx q[0];
rz(3.0786247) q[0];
rz(3.0175623) q[1];
sx q[1];
rz(-2.3359559) q[1];
sx q[1];
rz(-0.45809349) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9034018) q[0];
sx q[0];
rz(-1.5684677) q[0];
sx q[0];
rz(-2.6989614) q[0];
x q[1];
rz(2.6890254) q[2];
sx q[2];
rz(-1.0107702) q[2];
sx q[2];
rz(0.82816154) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.3921567) q[1];
sx q[1];
rz(-1.6877618) q[1];
sx q[1];
rz(-0.82526308) q[1];
rz(1.0701418) q[3];
sx q[3];
rz(-1.4874914) q[3];
sx q[3];
rz(1.3798151) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.9690341) q[2];
sx q[2];
rz(-2.2183552) q[2];
sx q[2];
rz(0.22949533) q[2];
rz(0.0028006639) q[3];
sx q[3];
rz(-2.2715748) q[3];
sx q[3];
rz(1.3389448) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5620419) q[0];
sx q[0];
rz(-0.27844772) q[0];
sx q[0];
rz(-0.91947412) q[0];
rz(-0.062285034) q[1];
sx q[1];
rz(-2.1376164) q[1];
sx q[1];
rz(1.8744291) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.29612449) q[0];
sx q[0];
rz(-0.86588973) q[0];
sx q[0];
rz(-2.092093) q[0];
rz(-pi) q[1];
rz(2.3976372) q[2];
sx q[2];
rz(-1.4641054) q[2];
sx q[2];
rz(0.26847408) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.25176469) q[1];
sx q[1];
rz(-1.7675195) q[1];
sx q[1];
rz(1.9038492) q[1];
x q[2];
rz(2.8980428) q[3];
sx q[3];
rz(-1.1034414) q[3];
sx q[3];
rz(-1.5598284) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.9617812) q[2];
sx q[2];
rz(-2.2634025) q[2];
sx q[2];
rz(-0.58376694) q[2];
rz(-0.70872712) q[3];
sx q[3];
rz(-1.3112336) q[3];
sx q[3];
rz(3.1183929) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.07847438) q[0];
sx q[0];
rz(-2.6874976) q[0];
sx q[0];
rz(1.0725347) q[0];
rz(2.5947) q[1];
sx q[1];
rz(-1.2439589) q[1];
sx q[1];
rz(-1.1118719) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3355013) q[0];
sx q[0];
rz(-1.4772692) q[0];
sx q[0];
rz(-2.0557687) q[0];
rz(-pi) q[1];
rz(-2.578031) q[2];
sx q[2];
rz(-2.5346018) q[2];
sx q[2];
rz(2.3129472) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.3355616) q[1];
sx q[1];
rz(-0.30004382) q[1];
sx q[1];
rz(-2.8801444) q[1];
rz(-pi) q[2];
rz(0.38970077) q[3];
sx q[3];
rz(-2.0412363) q[3];
sx q[3];
rz(0.54159347) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.303858) q[2];
sx q[2];
rz(-1.4759109) q[2];
sx q[2];
rz(-1.973935) q[2];
rz(1.5363103) q[3];
sx q[3];
rz(-1.4669908) q[3];
sx q[3];
rz(1.7355828) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4090356) q[0];
sx q[0];
rz(-1.5719825) q[0];
sx q[0];
rz(2.4107966) q[0];
rz(-0.90019512) q[1];
sx q[1];
rz(-0.80454818) q[1];
sx q[1];
rz(2.3866167) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2728111) q[0];
sx q[0];
rz(-0.73505721) q[0];
sx q[0];
rz(1.3520157) q[0];
rz(2.5078012) q[2];
sx q[2];
rz(-2.3049424) q[2];
sx q[2];
rz(-1.6939236) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.8440486) q[1];
sx q[1];
rz(-0.58931671) q[1];
sx q[1];
rz(0.20357666) q[1];
x q[2];
rz(2.2488392) q[3];
sx q[3];
rz(-2.4340981) q[3];
sx q[3];
rz(2.1219818) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.3770611) q[2];
sx q[2];
rz(-1.7691282) q[2];
sx q[2];
rz(0.63684741) q[2];
rz(0.26646715) q[3];
sx q[3];
rz(-1.0576495) q[3];
sx q[3];
rz(-1.5554957) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-2.4144142) q[0];
sx q[0];
rz(-1.1295015) q[0];
sx q[0];
rz(-1.138858) q[0];
rz(-0.75421929) q[1];
sx q[1];
rz(-0.33640877) q[1];
sx q[1];
rz(-0.019502217) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0855904) q[0];
sx q[0];
rz(-1.3645571) q[0];
sx q[0];
rz(-0.081736728) q[0];
rz(-pi) q[1];
x q[1];
rz(0.15140622) q[2];
sx q[2];
rz(-1.3726241) q[2];
sx q[2];
rz(-2.1422447) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.86917415) q[1];
sx q[1];
rz(-2.4915016) q[1];
sx q[1];
rz(-1.8723349) q[1];
x q[2];
rz(0.1653413) q[3];
sx q[3];
rz(-1.8792361) q[3];
sx q[3];
rz(0.43499085) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-3.1372244) q[2];
sx q[2];
rz(-1.4164111) q[2];
sx q[2];
rz(2.2926889) q[2];
rz(-2.7539339) q[3];
sx q[3];
rz(-1.1281697) q[3];
sx q[3];
rz(-1.6000115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
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
rz(1.7983109) q[0];
sx q[0];
rz(-0.16769519) q[0];
sx q[0];
rz(-2.6570901) q[0];
rz(-1.3867406) q[1];
sx q[1];
rz(-1.4258899) q[1];
sx q[1];
rz(1.1482931) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.590811) q[0];
sx q[0];
rz(-1.592144) q[0];
sx q[0];
rz(0.14069964) q[0];
rz(-pi) q[1];
rz(-1.3482773) q[2];
sx q[2];
rz(-2.1902124) q[2];
sx q[2];
rz(-2.7811108) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.53459586) q[1];
sx q[1];
rz(-1.7417223) q[1];
sx q[1];
rz(2.2439438) q[1];
x q[2];
rz(-2.9989472) q[3];
sx q[3];
rz(-2.7890165) q[3];
sx q[3];
rz(-2.8838317) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.91223532) q[2];
sx q[2];
rz(-1.2939913) q[2];
sx q[2];
rz(-2.7764376) q[2];
rz(0.12864104) q[3];
sx q[3];
rz(-1.235685) q[3];
sx q[3];
rz(0.45583367) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0317595) q[0];
sx q[0];
rz(-0.84072996) q[0];
sx q[0];
rz(1.6054556) q[0];
rz(0.96314349) q[1];
sx q[1];
rz(-1.8704725) q[1];
sx q[1];
rz(2.0830547) q[1];
rz(1.1744432) q[2];
sx q[2];
rz(-2.195993) q[2];
sx q[2];
rz(3.1396951) q[2];
rz(3.1046042) q[3];
sx q[3];
rz(-2.130641) q[3];
sx q[3];
rz(-0.11161042) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
