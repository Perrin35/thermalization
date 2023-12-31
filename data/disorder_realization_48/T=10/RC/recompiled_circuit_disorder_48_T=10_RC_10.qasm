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
rz(-0.9442803) q[1];
sx q[1];
rz(1.2184719) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.96555644) q[0];
sx q[0];
rz(-0.82057014) q[0];
sx q[0];
rz(1.6198938) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.081131447) q[2];
sx q[2];
rz(-0.46617026) q[2];
sx q[2];
rz(-2.1612338) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.54108809) q[1];
sx q[1];
rz(-2.0588015) q[1];
sx q[1];
rz(-1.0469251) q[1];
rz(-2.9524607) q[3];
sx q[3];
rz(-1.7889708) q[3];
sx q[3];
rz(1.8060341) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.618764) q[2];
sx q[2];
rz(-0.4814119) q[2];
sx q[2];
rz(2.5640326) q[2];
rz(1.1497568) q[3];
sx q[3];
rz(-1.3883608) q[3];
sx q[3];
rz(0.66453385) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.98786551) q[0];
sx q[0];
rz(-2.5550714) q[0];
sx q[0];
rz(-0.38744774) q[0];
rz(-2.2024343) q[1];
sx q[1];
rz(-0.99717957) q[1];
sx q[1];
rz(1.4025677) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.70799202) q[0];
sx q[0];
rz(-0.029768243) q[0];
sx q[0];
rz(3.004651) q[0];
rz(-0.34105532) q[2];
sx q[2];
rz(-2.555072) q[2];
sx q[2];
rz(-1.7797433) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-3.1075588) q[1];
sx q[1];
rz(-2.4817913) q[1];
sx q[1];
rz(0.69114139) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.12947793) q[3];
sx q[3];
rz(-0.63614142) q[3];
sx q[3];
rz(-0.57487088) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.42276057) q[2];
sx q[2];
rz(-1.3964802) q[2];
sx q[2];
rz(0.31769162) q[2];
rz(-2.9348532) q[3];
sx q[3];
rz(-2.5419149) q[3];
sx q[3];
rz(-0.81682214) q[3];
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
x q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7725672) q[0];
sx q[0];
rz(-1.439753) q[0];
sx q[0];
rz(1.4136219) q[0];
rz(-0.47779045) q[1];
sx q[1];
rz(-1.3505892) q[1];
sx q[1];
rz(-0.40107045) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6221878) q[0];
sx q[0];
rz(-1.6739968) q[0];
sx q[0];
rz(1.2445356) q[0];
rz(-pi) q[1];
x q[1];
rz(0.6443278) q[2];
sx q[2];
rz(-1.4187078) q[2];
sx q[2];
rz(-0.49567859) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.8071825) q[1];
sx q[1];
rz(-1.5644329) q[1];
sx q[1];
rz(-0.73015726) q[1];
x q[2];
rz(1.1776351) q[3];
sx q[3];
rz(-1.3893681) q[3];
sx q[3];
rz(1.7593311) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.59427375) q[2];
sx q[2];
rz(-1.6197562) q[2];
sx q[2];
rz(2.5857914) q[2];
rz(-2.1650971) q[3];
sx q[3];
rz(-0.54978168) q[3];
sx q[3];
rz(2.3613789) q[3];
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
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.34898409) q[0];
sx q[0];
rz(-1.5190834) q[0];
sx q[0];
rz(1.4439616) q[0];
rz(-1.6216888) q[1];
sx q[1];
rz(-2.4855721) q[1];
sx q[1];
rz(2.8881883) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.69617535) q[0];
sx q[0];
rz(-1.4679969) q[0];
sx q[0];
rz(2.1131383) q[0];
x q[1];
rz(3.1242712) q[2];
sx q[2];
rz(-2.0565196) q[2];
sx q[2];
rz(2.3522365) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.464444) q[1];
sx q[1];
rz(-2.1975937) q[1];
sx q[1];
rz(1.6897175) q[1];
rz(-0.6494135) q[3];
sx q[3];
rz(-1.6522437) q[3];
sx q[3];
rz(-0.34929517) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.0306586) q[2];
sx q[2];
rz(-1.3867644) q[2];
sx q[2];
rz(-0.78732642) q[2];
rz(2.2287255) q[3];
sx q[3];
rz(-2.392231) q[3];
sx q[3];
rz(-2.1319938) q[3];
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
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71516365) q[0];
sx q[0];
rz(-2.5029095) q[0];
sx q[0];
rz(3.0786247) q[0];
rz(0.12403034) q[1];
sx q[1];
rz(-2.3359559) q[1];
sx q[1];
rz(0.45809349) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3315017) q[0];
sx q[0];
rz(-1.1281663) q[0];
sx q[0];
rz(1.5733733) q[0];
rz(-pi) q[1];
x q[1];
rz(0.45256726) q[2];
sx q[2];
rz(-1.0107702) q[2];
sx q[2];
rz(-0.82816154) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.8372247) q[1];
sx q[1];
rz(-2.3886884) q[1];
sx q[1];
rz(1.7423082) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7430274) q[3];
sx q[3];
rz(-0.50695626) q[3];
sx q[3];
rz(2.7996922) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.9690341) q[2];
sx q[2];
rz(-0.92323747) q[2];
sx q[2];
rz(-0.22949533) q[2];
rz(0.0028006639) q[3];
sx q[3];
rz(-2.2715748) q[3];
sx q[3];
rz(1.3389448) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5620419) q[0];
sx q[0];
rz(-0.27844772) q[0];
sx q[0];
rz(-0.91947412) q[0];
rz(-3.0793076) q[1];
sx q[1];
rz(-2.1376164) q[1];
sx q[1];
rz(1.2671635) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4290659) q[0];
sx q[0];
rz(-0.8493087) q[0];
sx q[0];
rz(2.6119786) q[0];
rz(-0.74395545) q[2];
sx q[2];
rz(-1.6774872) q[2];
sx q[2];
rz(2.8731186) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.2515182) q[1];
sx q[1];
rz(-1.8971844) q[1];
sx q[1];
rz(-0.20784394) q[1];
rz(1.1250161) q[3];
sx q[3];
rz(-0.52281724) q[3];
sx q[3];
rz(1.0558053) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.9617812) q[2];
sx q[2];
rz(-0.87819019) q[2];
sx q[2];
rz(2.5578257) q[2];
rz(-2.4328655) q[3];
sx q[3];
rz(-1.3112336) q[3];
sx q[3];
rz(-3.1183929) q[3];
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
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.07847438) q[0];
sx q[0];
rz(-0.45409504) q[0];
sx q[0];
rz(-2.069058) q[0];
rz(2.5947) q[1];
sx q[1];
rz(-1.8976338) q[1];
sx q[1];
rz(-2.0297208) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.81387732) q[0];
sx q[0];
rz(-1.0881256) q[0];
sx q[0];
rz(3.0359603) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6107437) q[2];
sx q[2];
rz(-1.2611654) q[2];
sx q[2];
rz(-1.9206778) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.062473) q[1];
sx q[1];
rz(-1.2812496) q[1];
sx q[1];
rz(1.4909966) q[1];
x q[2];
rz(-0.92915793) q[3];
sx q[3];
rz(-0.60141364) q[3];
sx q[3];
rz(1.2777559) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.303858) q[2];
sx q[2];
rz(-1.4759109) q[2];
sx q[2];
rz(-1.1676577) q[2];
rz(1.5363103) q[3];
sx q[3];
rz(-1.6746018) q[3];
sx q[3];
rz(1.4060098) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4090356) q[0];
sx q[0];
rz(-1.5719825) q[0];
sx q[0];
rz(-0.73079601) q[0];
rz(0.90019512) q[1];
sx q[1];
rz(-2.3370445) q[1];
sx q[1];
rz(-0.75497595) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2728111) q[0];
sx q[0];
rz(-2.4065354) q[0];
sx q[0];
rz(-1.3520157) q[0];
x q[1];
rz(-2.5078012) q[2];
sx q[2];
rz(-2.3049424) q[2];
sx q[2];
rz(1.6939236) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.6983812) q[1];
sx q[1];
rz(-1.6834007) q[1];
sx q[1];
rz(-0.57971445) q[1];
x q[2];
rz(-0.98324361) q[3];
sx q[3];
rz(-1.1508815) q[3];
sx q[3];
rz(-1.1004694) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.3770611) q[2];
sx q[2];
rz(-1.7691282) q[2];
sx q[2];
rz(0.63684741) q[2];
rz(2.8751255) q[3];
sx q[3];
rz(-2.0839432) q[3];
sx q[3];
rz(1.586097) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4144142) q[0];
sx q[0];
rz(-2.0120912) q[0];
sx q[0];
rz(-2.0027347) q[0];
rz(-2.3873734) q[1];
sx q[1];
rz(-2.8051839) q[1];
sx q[1];
rz(-0.019502217) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0560023) q[0];
sx q[0];
rz(-1.7770355) q[0];
sx q[0];
rz(3.0598559) q[0];
rz(0.92658652) q[2];
sx q[2];
rz(-0.24878657) q[2];
sx q[2];
rz(-2.8015346) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.86917415) q[1];
sx q[1];
rz(-2.4915016) q[1];
sx q[1];
rz(1.8723349) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8832302) q[3];
sx q[3];
rz(-1.7282681) q[3];
sx q[3];
rz(-1.0851932) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.0043682178) q[2];
sx q[2];
rz(-1.4164111) q[2];
sx q[2];
rz(-0.84890378) q[2];
rz(0.38765872) q[3];
sx q[3];
rz(-2.013423) q[3];
sx q[3];
rz(1.6000115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
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
rz(-1.9932995) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.590811) q[0];
sx q[0];
rz(-1.592144) q[0];
sx q[0];
rz(-3.000893) q[0];
rz(-pi) q[1];
rz(-0.30015595) q[2];
sx q[2];
rz(-0.65320063) q[2];
sx q[2];
rz(2.4094827) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.3155568) q[1];
sx q[1];
rz(-2.450374) q[1];
sx q[1];
rz(-1.8408937) q[1];
x q[2];
rz(-0.14264543) q[3];
sx q[3];
rz(-2.7890165) q[3];
sx q[3];
rz(2.8838317) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.91223532) q[2];
sx q[2];
rz(-1.2939913) q[2];
sx q[2];
rz(0.36515507) q[2];
rz(3.0129516) q[3];
sx q[3];
rz(-1.9059076) q[3];
sx q[3];
rz(-2.685759) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
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
rz(-0.1098332) q[0];
sx q[0];
rz(-2.3008627) q[0];
sx q[0];
rz(-1.536137) q[0];
rz(0.96314349) q[1];
sx q[1];
rz(-1.8704725) q[1];
sx q[1];
rz(2.0830547) q[1];
rz(-0.49114901) q[2];
sx q[2];
rz(-0.72577234) q[2];
sx q[2];
rz(2.5189248) q[2];
rz(2.1309489) q[3];
sx q[3];
rz(-1.602136) q[3];
sx q[3];
rz(-1.7020561) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
