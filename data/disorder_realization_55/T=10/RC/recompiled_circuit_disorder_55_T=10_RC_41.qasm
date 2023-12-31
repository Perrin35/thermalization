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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.41266325) q[0];
sx q[0];
rz(-0.67617765) q[0];
sx q[0];
rz(0.23763188) q[0];
x q[1];
rz(-2.5392883) q[2];
sx q[2];
rz(-2.3659083) q[2];
sx q[2];
rz(1.3210981) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.6890251) q[1];
sx q[1];
rz(-1.8568294) q[1];
sx q[1];
rz(2.9806115) q[1];
rz(3.0987708) q[3];
sx q[3];
rz(-2.560727) q[3];
sx q[3];
rz(-0.37561852) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.66951093) q[2];
sx q[2];
rz(-1.8410204) q[2];
sx q[2];
rz(2.0377339) q[2];
rz(-1.8707229) q[3];
sx q[3];
rz(-1.9138252) q[3];
sx q[3];
rz(-2.8675573) q[3];
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
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0274149) q[0];
sx q[0];
rz(-2.6129621) q[0];
sx q[0];
rz(0.43637481) q[0];
rz(2.6787058) q[1];
sx q[1];
rz(-1.0375689) q[1];
sx q[1];
rz(-2.8754821) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55914315) q[0];
sx q[0];
rz(-2.3087915) q[0];
sx q[0];
rz(1.0484496) q[0];
rz(-2.1462609) q[2];
sx q[2];
rz(-1.7082214) q[2];
sx q[2];
rz(2.1910138) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.3937711) q[1];
sx q[1];
rz(-2.1854679) q[1];
sx q[1];
rz(0.6786896) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.9640785) q[3];
sx q[3];
rz(-1.420141) q[3];
sx q[3];
rz(-1.3552624) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.6767072) q[2];
sx q[2];
rz(-1.2767982) q[2];
sx q[2];
rz(-2.6300988) q[2];
rz(-0.809508) q[3];
sx q[3];
rz(-1.5313238) q[3];
sx q[3];
rz(-2.8539343) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.70616102) q[0];
sx q[0];
rz(-1.5478739) q[0];
sx q[0];
rz(0.92873746) q[0];
rz(1.4061032) q[1];
sx q[1];
rz(-0.7000674) q[1];
sx q[1];
rz(1.7944638) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1613306) q[0];
sx q[0];
rz(-1.3093003) q[0];
sx q[0];
rz(0.47388347) q[0];
rz(-pi) q[1];
rz(0.15813078) q[2];
sx q[2];
rz(-1.2889382) q[2];
sx q[2];
rz(2.6513197) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.9007064) q[1];
sx q[1];
rz(-2.1623003) q[1];
sx q[1];
rz(2.790931) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9903528) q[3];
sx q[3];
rz(-2.2721707) q[3];
sx q[3];
rz(-0.73089862) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.1469664) q[2];
sx q[2];
rz(-1.3866084) q[2];
sx q[2];
rz(1.5220801) q[2];
rz(2.8772723) q[3];
sx q[3];
rz(-2.1051354) q[3];
sx q[3];
rz(-2.6456397) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.79214823) q[0];
sx q[0];
rz(-1.9932207) q[0];
sx q[0];
rz(2.175892) q[0];
rz(0.72215885) q[1];
sx q[1];
rz(-1.5042217) q[1];
sx q[1];
rz(2.5818363) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1346261) q[0];
sx q[0];
rz(-2.7077733) q[0];
sx q[0];
rz(-0.19402786) q[0];
x q[1];
rz(-2.431805) q[2];
sx q[2];
rz(-2.1927532) q[2];
sx q[2];
rz(2.4316535) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.2985845) q[1];
sx q[1];
rz(-2.1857939) q[1];
sx q[1];
rz(-3.0857012) q[1];
rz(-0.23493725) q[3];
sx q[3];
rz(-0.80312356) q[3];
sx q[3];
rz(-1.1266176) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.42797783) q[2];
sx q[2];
rz(-1.5087936) q[2];
sx q[2];
rz(1.7170061) q[2];
rz(2.8811841) q[3];
sx q[3];
rz(-1.7491165) q[3];
sx q[3];
rz(-2.7105455) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.150862) q[0];
sx q[0];
rz(-1.6435511) q[0];
sx q[0];
rz(2.7752303) q[0];
rz(1.5953966) q[1];
sx q[1];
rz(-2.5876744) q[1];
sx q[1];
rz(-0.34367925) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8415547) q[0];
sx q[0];
rz(-2.2404788) q[0];
sx q[0];
rz(0.71398736) q[0];
rz(-pi) q[1];
rz(2.5630066) q[2];
sx q[2];
rz(-0.95539504) q[2];
sx q[2];
rz(1.3784642) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.1829454) q[1];
sx q[1];
rz(-1.4970386) q[1];
sx q[1];
rz(-1.0720836) q[1];
rz(-pi) q[2];
rz(-1.5991391) q[3];
sx q[3];
rz(-0.64405555) q[3];
sx q[3];
rz(1.6213662) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.3917824) q[2];
sx q[2];
rz(-2.2822773) q[2];
sx q[2];
rz(-0.053744944) q[2];
rz(-1.7371477) q[3];
sx q[3];
rz(-2.6440547) q[3];
sx q[3];
rz(-0.29156175) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6546201) q[0];
sx q[0];
rz(-0.64990652) q[0];
sx q[0];
rz(-0.75575954) q[0];
rz(-0.02515633) q[1];
sx q[1];
rz(-2.2143366) q[1];
sx q[1];
rz(-2.8818534) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.040859) q[0];
sx q[0];
rz(-1.20964) q[0];
sx q[0];
rz(-0.33613899) q[0];
rz(-pi) q[1];
rz(-0.01979205) q[2];
sx q[2];
rz(-1.9070101) q[2];
sx q[2];
rz(0.83376955) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.427634) q[1];
sx q[1];
rz(-2.7379588) q[1];
sx q[1];
rz(0.9637109) q[1];
rz(-pi) q[2];
rz(1.5335347) q[3];
sx q[3];
rz(-2.691949) q[3];
sx q[3];
rz(0.26526181) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.48866895) q[2];
sx q[2];
rz(-0.80604625) q[2];
sx q[2];
rz(-0.10061131) q[2];
rz(-2.9595024) q[3];
sx q[3];
rz(-0.87825769) q[3];
sx q[3];
rz(1.7939059) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1473734) q[0];
sx q[0];
rz(-1.4202776) q[0];
sx q[0];
rz(2.8314262) q[0];
rz(-2.639333) q[1];
sx q[1];
rz(-0.52400932) q[1];
sx q[1];
rz(2.535634) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.92354846) q[0];
sx q[0];
rz(-2.179562) q[0];
sx q[0];
rz(1.2952842) q[0];
rz(-pi) q[1];
rz(1.3700968) q[2];
sx q[2];
rz(-1.2846652) q[2];
sx q[2];
rz(1.8482894) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.0026605) q[1];
sx q[1];
rz(-1.9117038) q[1];
sx q[1];
rz(2.1041811) q[1];
x q[2];
rz(-0.53020729) q[3];
sx q[3];
rz(-1.0876417) q[3];
sx q[3];
rz(2.9999441) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.24017748) q[2];
sx q[2];
rz(-1.1837974) q[2];
sx q[2];
rz(2.288738) q[2];
rz(-1.3700221) q[3];
sx q[3];
rz(-1.6829237) q[3];
sx q[3];
rz(2.9366233) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
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
rz(1.9119499) q[0];
sx q[0];
rz(-2.5456173) q[0];
sx q[0];
rz(-1.6802616) q[0];
rz(1.4029067) q[1];
sx q[1];
rz(-0.97424126) q[1];
sx q[1];
rz(-0.064037474) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.56211573) q[0];
sx q[0];
rz(-0.96456438) q[0];
sx q[0];
rz(-2.9129145) q[0];
rz(1.0649101) q[2];
sx q[2];
rz(-1.8407028) q[2];
sx q[2];
rz(0.050886521) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.21676) q[1];
sx q[1];
rz(-1.2158582) q[1];
sx q[1];
rz(-1.6096398) q[1];
rz(-pi) q[2];
x q[2];
rz(0.039489432) q[3];
sx q[3];
rz(-1.2590623) q[3];
sx q[3];
rz(-1.5962275) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-3.0503851) q[2];
sx q[2];
rz(-0.63082266) q[2];
sx q[2];
rz(-1.5861661) q[2];
rz(0.88820109) q[3];
sx q[3];
rz(-1.9017838) q[3];
sx q[3];
rz(-0.92938882) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9973307) q[0];
sx q[0];
rz(-2.0781131) q[0];
sx q[0];
rz(-1.7096747) q[0];
rz(-0.56888467) q[1];
sx q[1];
rz(-0.535393) q[1];
sx q[1];
rz(1.127839) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.226798) q[0];
sx q[0];
rz(-2.9920122) q[0];
sx q[0];
rz(3.0339255) q[0];
rz(-2.1986507) q[2];
sx q[2];
rz(-1.1254416) q[2];
sx q[2];
rz(-3.359059e-05) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.3582663) q[1];
sx q[1];
rz(-1.8124541) q[1];
sx q[1];
rz(-0.011947167) q[1];
rz(-pi) q[2];
rz(1.0206251) q[3];
sx q[3];
rz(-0.63311011) q[3];
sx q[3];
rz(1.8819295) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.52788064) q[2];
sx q[2];
rz(-0.85313672) q[2];
sx q[2];
rz(-2.9679427) q[2];
rz(0.33637834) q[3];
sx q[3];
rz(-1.222638) q[3];
sx q[3];
rz(2.7500847) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7062475) q[0];
sx q[0];
rz(-2.5543537) q[0];
sx q[0];
rz(-1.6760814) q[0];
rz(-0.82410518) q[1];
sx q[1];
rz(-1.5287639) q[1];
sx q[1];
rz(-2.5691659) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3819645) q[0];
sx q[0];
rz(-2.1158764) q[0];
sx q[0];
rz(0.99960021) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.053439) q[2];
sx q[2];
rz(-0.51689076) q[2];
sx q[2];
rz(3.0992532) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(3.010667) q[1];
sx q[1];
rz(-1.7867242) q[1];
sx q[1];
rz(0.24433498) q[1];
x q[2];
rz(0.33705538) q[3];
sx q[3];
rz(-2.4042077) q[3];
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
rz(-0.47452351) q[2];
sx q[2];
rz(2.6043716) q[2];
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
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2873516) q[0];
sx q[0];
rz(-1.1653405) q[0];
sx q[0];
rz(-1.5821138) q[0];
rz(2.6782425) q[1];
sx q[1];
rz(-0.87711038) q[1];
sx q[1];
rz(-1.6323485) q[1];
rz(1.6981381) q[2];
sx q[2];
rz(-0.6586532) q[2];
sx q[2];
rz(-2.3011343) q[2];
rz(1.4544009) q[3];
sx q[3];
rz(-2.6506861) q[3];
sx q[3];
rz(-0.41674137) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
