OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.5947333) q[0];
sx q[0];
rz(-1.5164627) q[0];
sx q[0];
rz(-2.8773142) q[0];
rz(-0.9737941) q[1];
sx q[1];
rz(5.073054) q[1];
sx q[1];
rz(10.160025) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3448479) q[0];
sx q[0];
rz(-1.42294) q[0];
sx q[0];
rz(-2.4792838) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0779607) q[2];
sx q[2];
rz(-2.1858366) q[2];
sx q[2];
rz(-2.087649) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.4525675) q[1];
sx q[1];
rz(-1.2847632) q[1];
sx q[1];
rz(2.9806115) q[1];
x q[2];
rz(-2.5611476) q[3];
sx q[3];
rz(-1.547303) q[3];
sx q[3];
rz(1.1593727) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.4720817) q[2];
sx q[2];
rz(-1.8410204) q[2];
sx q[2];
rz(-1.1038587) q[2];
rz(1.2708698) q[3];
sx q[3];
rz(-1.9138252) q[3];
sx q[3];
rz(0.27403533) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1141777) q[0];
sx q[0];
rz(-2.6129621) q[0];
sx q[0];
rz(-0.43637481) q[0];
rz(-0.46288681) q[1];
sx q[1];
rz(-2.1040237) q[1];
sx q[1];
rz(2.8754821) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9929745) q[0];
sx q[0];
rz(-0.87478144) q[0];
sx q[0];
rz(-0.50177411) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1462609) q[2];
sx q[2];
rz(-1.7082214) q[2];
sx q[2];
rz(2.1910138) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.38763603) q[1];
sx q[1];
rz(-2.1093183) q[1];
sx q[1];
rz(-2.3073879) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.9640785) q[3];
sx q[3];
rz(-1.7214516) q[3];
sx q[3];
rz(1.3552624) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.46488547) q[2];
sx q[2];
rz(-1.2767982) q[2];
sx q[2];
rz(2.6300988) q[2];
rz(-2.3320847) q[3];
sx q[3];
rz(-1.5313238) q[3];
sx q[3];
rz(2.8539343) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.70616102) q[0];
sx q[0];
rz(-1.5478739) q[0];
sx q[0];
rz(-2.2128552) q[0];
rz(1.4061032) q[1];
sx q[1];
rz(-2.4415253) q[1];
sx q[1];
rz(1.3471289) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.72235332) q[0];
sx q[0];
rz(-2.0273211) q[0];
sx q[0];
rz(1.2786352) q[0];
rz(-pi) q[1];
rz(-2.9834619) q[2];
sx q[2];
rz(-1.8526544) q[2];
sx q[2];
rz(0.49027298) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.0128847) q[1];
sx q[1];
rz(-1.2816267) q[1];
sx q[1];
rz(-2.1916926) q[1];
rz(-pi) q[2];
rz(-2.2778355) q[3];
sx q[3];
rz(-1.6861526) q[3];
sx q[3];
rz(2.399721) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.1469664) q[2];
sx q[2];
rz(-1.3866084) q[2];
sx q[2];
rz(-1.5220801) q[2];
rz(0.26432031) q[3];
sx q[3];
rz(-2.1051354) q[3];
sx q[3];
rz(2.6456397) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3494444) q[0];
sx q[0];
rz(-1.9932207) q[0];
sx q[0];
rz(-2.175892) q[0];
rz(-2.4194338) q[1];
sx q[1];
rz(-1.5042217) q[1];
sx q[1];
rz(2.5818363) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9352919) q[0];
sx q[0];
rz(-1.1456523) q[0];
sx q[0];
rz(1.6598808) q[0];
rz(-2.3085924) q[2];
sx q[2];
rz(-0.90655316) q[2];
sx q[2];
rz(-0.26495648) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.30448118) q[1];
sx q[1];
rz(-1.5251535) q[1];
sx q[1];
rz(-0.95506217) q[1];
rz(0.78919952) q[3];
sx q[3];
rz(-1.4025098) q[3];
sx q[3];
rz(2.5326953) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.42797783) q[2];
sx q[2];
rz(-1.6327991) q[2];
sx q[2];
rz(1.4245865) q[2];
rz(0.26040855) q[3];
sx q[3];
rz(-1.3924761) q[3];
sx q[3];
rz(-2.7105455) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.150862) q[0];
sx q[0];
rz(-1.6435511) q[0];
sx q[0];
rz(-0.36636233) q[0];
rz(1.5953966) q[1];
sx q[1];
rz(-2.5876744) q[1];
sx q[1];
rz(2.7979134) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.30003798) q[0];
sx q[0];
rz(-2.2404788) q[0];
sx q[0];
rz(-2.4276053) q[0];
rz(2.5630066) q[2];
sx q[2];
rz(-0.95539504) q[2];
sx q[2];
rz(-1.7631284) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.9586473) q[1];
sx q[1];
rz(-1.4970386) q[1];
sx q[1];
rz(2.069509) q[1];
rz(-pi) q[2];
rz(-1.5991391) q[3];
sx q[3];
rz(-2.4975371) q[3];
sx q[3];
rz(-1.6213662) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.3917824) q[2];
sx q[2];
rz(-0.85931531) q[2];
sx q[2];
rz(3.0878477) q[2];
rz(1.7371477) q[3];
sx q[3];
rz(-0.497538) q[3];
sx q[3];
rz(2.8500309) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6546201) q[0];
sx q[0];
rz(-2.4916861) q[0];
sx q[0];
rz(-0.75575954) q[0];
rz(0.02515633) q[1];
sx q[1];
rz(-2.2143366) q[1];
sx q[1];
rz(-0.25973928) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.040859) q[0];
sx q[0];
rz(-1.20964) q[0];
sx q[0];
rz(0.33613899) q[0];
x q[1];
rz(-3.1218006) q[2];
sx q[2];
rz(-1.9070101) q[2];
sx q[2];
rz(2.3078231) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.0670843) q[1];
sx q[1];
rz(-1.2423406) q[1];
sx q[1];
rz(-2.902608) q[1];
rz(-pi) q[2];
x q[2];
rz(1.608058) q[3];
sx q[3];
rz(-2.691949) q[3];
sx q[3];
rz(2.8763308) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.6529237) q[2];
sx q[2];
rz(-0.80604625) q[2];
sx q[2];
rz(0.10061131) q[2];
rz(-0.18209022) q[3];
sx q[3];
rz(-0.87825769) q[3];
sx q[3];
rz(1.3476868) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9942193) q[0];
sx q[0];
rz(-1.7213151) q[0];
sx q[0];
rz(0.31016645) q[0];
rz(0.50225964) q[1];
sx q[1];
rz(-2.6175833) q[1];
sx q[1];
rz(-2.535634) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.80752559) q[0];
sx q[0];
rz(-1.345732) q[0];
sx q[0];
rz(2.5146757) q[0];
x q[1];
rz(-0.59553869) q[2];
sx q[2];
rz(-2.7936802) q[2];
sx q[2];
rz(2.4728342) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.76304945) q[1];
sx q[1];
rz(-2.0705283) q[1];
sx q[1];
rz(-0.39079697) q[1];
rz(-pi) q[2];
rz(2.3378387) q[3];
sx q[3];
rz(-0.70138068) q[3];
sx q[3];
rz(-2.3826117) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.9014152) q[2];
sx q[2];
rz(-1.9577953) q[2];
sx q[2];
rz(-0.85285464) q[2];
rz(-1.3700221) q[3];
sx q[3];
rz(-1.6829237) q[3];
sx q[3];
rz(2.9366233) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9119499) q[0];
sx q[0];
rz(-0.59597534) q[0];
sx q[0];
rz(-1.461331) q[0];
rz(-1.7386859) q[1];
sx q[1];
rz(-0.97424126) q[1];
sx q[1];
rz(-0.064037474) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1916699) q[0];
sx q[0];
rz(-0.64282286) q[0];
sx q[0];
rz(1.886801) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0520347) q[2];
sx q[2];
rz(-2.5737737) q[2];
sx q[2];
rz(-2.0702814) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.1053972) q[1];
sx q[1];
rz(-2.7846249) q[1];
sx q[1];
rz(3.0371975) q[1];
x q[2];
rz(1.6927034) q[3];
sx q[3];
rz(-0.31414437) q[3];
sx q[3];
rz(-1.4172518) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.091207592) q[2];
sx q[2];
rz(-2.51077) q[2];
sx q[2];
rz(-1.5554265) q[2];
rz(-0.88820109) q[3];
sx q[3];
rz(-1.9017838) q[3];
sx q[3];
rz(-2.2122038) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9973307) q[0];
sx q[0];
rz(-2.0781131) q[0];
sx q[0];
rz(-1.7096747) q[0];
rz(-2.572708) q[1];
sx q[1];
rz(-2.6061997) q[1];
sx q[1];
rz(-2.0137537) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8059212) q[0];
sx q[0];
rz(-1.7195042) q[0];
sx q[0];
rz(-1.554603) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.94294195) q[2];
sx q[2];
rz(-2.0161511) q[2];
sx q[2];
rz(-3.359059e-05) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.20967083) q[1];
sx q[1];
rz(-1.5591963) q[1];
sx q[1];
rz(-1.329122) q[1];
x q[2];
rz(-1.0117709) q[3];
sx q[3];
rz(-1.8852919) q[3];
sx q[3];
rz(0.77034706) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.52788064) q[2];
sx q[2];
rz(-0.85313672) q[2];
sx q[2];
rz(2.9679427) q[2];
rz(0.33637834) q[3];
sx q[3];
rz(-1.9189546) q[3];
sx q[3];
rz(-2.7500847) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7062475) q[0];
sx q[0];
rz(-0.58723891) q[0];
sx q[0];
rz(1.6760814) q[0];
rz(-2.3174875) q[1];
sx q[1];
rz(-1.5287639) q[1];
sx q[1];
rz(-0.5724268) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7596282) q[0];
sx q[0];
rz(-1.0257162) q[0];
sx q[0];
rz(-0.99960021) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0296713) q[2];
sx q[2];
rz(-1.8177114) q[2];
sx q[2];
rz(-1.0690451) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(3.010667) q[1];
sx q[1];
rz(-1.7867242) q[1];
sx q[1];
rz(0.24433498) q[1];
rz(-pi) q[2];
rz(2.8045373) q[3];
sx q[3];
rz(-0.73738499) q[3];
sx q[3];
rz(2.9205703) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3616025) q[2];
sx q[2];
rz(-2.6670691) q[2];
sx q[2];
rz(-2.6043716) q[2];
rz(1.0572664) q[3];
sx q[3];
rz(-0.89151645) q[3];
sx q[3];
rz(2.4479772) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2873516) q[0];
sx q[0];
rz(-1.1653405) q[0];
sx q[0];
rz(-1.5821138) q[0];
rz(-0.46335012) q[1];
sx q[1];
rz(-0.87711038) q[1];
sx q[1];
rz(-1.6323485) q[1];
rz(2.2255185) q[2];
sx q[2];
rz(-1.4929885) q[2];
sx q[2];
rz(2.3103466) q[2];
rz(-3.0795931) q[3];
sx q[3];
rz(-2.0580895) q[3];
sx q[3];
rz(2.5930391) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
