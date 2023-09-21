OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.54685932) q[0];
sx q[0];
rz(-1.62513) q[0];
sx q[0];
rz(-0.2642785) q[0];
rz(-0.9737941) q[1];
sx q[1];
rz(-1.2101313) q[1];
sx q[1];
rz(0.73524737) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7967448) q[0];
sx q[0];
rz(-1.7186527) q[0];
sx q[0];
rz(0.66230886) q[0];
rz(-pi) q[1];
rz(-1.0636319) q[2];
sx q[2];
rz(-0.9557561) q[2];
sx q[2];
rz(2.087649) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.4525675) q[1];
sx q[1];
rz(-1.2847632) q[1];
sx q[1];
rz(2.9806115) q[1];
x q[2];
rz(-3.0987708) q[3];
sx q[3];
rz(-2.560727) q[3];
sx q[3];
rz(-2.7659741) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.4720817) q[2];
sx q[2];
rz(-1.8410204) q[2];
sx q[2];
rz(2.0377339) q[2];
rz(1.8707229) q[3];
sx q[3];
rz(-1.9138252) q[3];
sx q[3];
rz(2.8675573) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1141777) q[0];
sx q[0];
rz(-2.6129621) q[0];
sx q[0];
rz(-0.43637481) q[0];
rz(-2.6787058) q[1];
sx q[1];
rz(-2.1040237) q[1];
sx q[1];
rz(0.26611051) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5824495) q[0];
sx q[0];
rz(-0.83280116) q[0];
sx q[0];
rz(1.0484496) q[0];
rz(-0.16337784) q[2];
sx q[2];
rz(-2.1401569) q[2];
sx q[2];
rz(-0.53158224) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.3937711) q[1];
sx q[1];
rz(-2.1854679) q[1];
sx q[1];
rz(-0.6786896) q[1];
x q[2];
rz(2.9588685) q[3];
sx q[3];
rz(-0.97191873) q[3];
sx q[3];
rz(-3.0298508) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.46488547) q[2];
sx q[2];
rz(-1.8647944) q[2];
sx q[2];
rz(0.51149386) q[2];
rz(-0.809508) q[3];
sx q[3];
rz(-1.6102689) q[3];
sx q[3];
rz(2.8539343) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4354316) q[0];
sx q[0];
rz(-1.5937188) q[0];
sx q[0];
rz(0.92873746) q[0];
rz(1.7354895) q[1];
sx q[1];
rz(-2.4415253) q[1];
sx q[1];
rz(-1.3471289) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0179694) q[0];
sx q[0];
rz(-0.53640134) q[0];
sx q[0];
rz(0.53039741) q[0];
rz(-1.8560266) q[2];
sx q[2];
rz(-1.4189548) q[2];
sx q[2];
rz(1.0361995) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.82145065) q[1];
sx q[1];
rz(-2.4648033) q[1];
sx q[1];
rz(2.0435145) q[1];
rz(0.15123983) q[3];
sx q[3];
rz(-0.86942196) q[3];
sx q[3];
rz(2.410694) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.1469664) q[2];
sx q[2];
rz(-1.3866084) q[2];
sx q[2];
rz(-1.5220801) q[2];
rz(-0.26432031) q[3];
sx q[3];
rz(-2.1051354) q[3];
sx q[3];
rz(0.49595293) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(-0.79214823) q[0];
sx q[0];
rz(-1.148372) q[0];
sx q[0];
rz(2.175892) q[0];
rz(-0.72215885) q[1];
sx q[1];
rz(-1.5042217) q[1];
sx q[1];
rz(-2.5818363) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.0069665466) q[0];
sx q[0];
rz(-0.43381938) q[0];
sx q[0];
rz(-0.19402786) q[0];
x q[1];
rz(-0.70978769) q[2];
sx q[2];
rz(-0.94883942) q[2];
sx q[2];
rz(2.4316535) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.30448118) q[1];
sx q[1];
rz(-1.5251535) q[1];
sx q[1];
rz(0.95506217) q[1];
rz(-pi) q[2];
rz(-2.3523931) q[3];
sx q[3];
rz(-1.4025098) q[3];
sx q[3];
rz(-0.60889739) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.42797783) q[2];
sx q[2];
rz(-1.6327991) q[2];
sx q[2];
rz(1.7170061) q[2];
rz(0.26040855) q[3];
sx q[3];
rz(-1.3924761) q[3];
sx q[3];
rz(-2.7105455) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.99073064) q[0];
sx q[0];
rz(-1.4980415) q[0];
sx q[0];
rz(2.7752303) q[0];
rz(1.5461961) q[1];
sx q[1];
rz(-0.55391824) q[1];
sx q[1];
rz(2.7979134) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2491964) q[0];
sx q[0];
rz(-0.93660347) q[0];
sx q[0];
rz(-0.87974324) q[0];
x q[1];
rz(-0.86954388) q[2];
sx q[2];
rz(-1.1079259) q[2];
sx q[2];
rz(-2.9733544) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.1829454) q[1];
sx q[1];
rz(-1.6445541) q[1];
sx q[1];
rz(2.069509) q[1];
rz(-pi) q[2];
rz(1.5991391) q[3];
sx q[3];
rz(-0.64405555) q[3];
sx q[3];
rz(-1.6213662) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.7498103) q[2];
sx q[2];
rz(-2.2822773) q[2];
sx q[2];
rz(-3.0878477) q[2];
rz(-1.404445) q[3];
sx q[3];
rz(-0.497538) q[3];
sx q[3];
rz(2.8500309) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6546201) q[0];
sx q[0];
rz(-0.64990652) q[0];
sx q[0];
rz(2.3858331) q[0];
rz(0.02515633) q[1];
sx q[1];
rz(-2.2143366) q[1];
sx q[1];
rz(-0.25973928) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.320967) q[0];
sx q[0];
rz(-0.48829406) q[0];
sx q[0];
rz(0.85296209) q[0];
rz(1.2345215) q[2];
sx q[2];
rz(-1.5521126) q[2];
sx q[2];
rz(-0.74355723) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.0670843) q[1];
sx q[1];
rz(-1.2423406) q[1];
sx q[1];
rz(0.23898464) q[1];
rz(-3.1236157) q[3];
sx q[3];
rz(-1.1214876) q[3];
sx q[3];
rz(2.9177005) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.48866895) q[2];
sx q[2];
rz(-2.3355464) q[2];
sx q[2];
rz(3.0409813) q[2];
rz(2.9595024) q[3];
sx q[3];
rz(-2.263335) q[3];
sx q[3];
rz(-1.3476868) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9942193) q[0];
sx q[0];
rz(-1.7213151) q[0];
sx q[0];
rz(-0.31016645) q[0];
rz(-2.639333) q[1];
sx q[1];
rz(-2.6175833) q[1];
sx q[1];
rz(-2.535634) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46442407) q[0];
sx q[0];
rz(-2.4806528) q[0];
sx q[0];
rz(2.7695157) q[0];
x q[1];
rz(-1.7714959) q[2];
sx q[2];
rz(-1.2846652) q[2];
sx q[2];
rz(1.8482894) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.0026605) q[1];
sx q[1];
rz(-1.2298889) q[1];
sx q[1];
rz(-2.1041811) q[1];
rz(-pi) q[2];
rz(-0.53020729) q[3];
sx q[3];
rz(-2.0539509) q[3];
sx q[3];
rz(-2.9999441) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.9014152) q[2];
sx q[2];
rz(-1.9577953) q[2];
sx q[2];
rz(-2.288738) q[2];
rz(-1.7715706) q[3];
sx q[3];
rz(-1.6829237) q[3];
sx q[3];
rz(-2.9366233) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9119499) q[0];
sx q[0];
rz(-0.59597534) q[0];
sx q[0];
rz(-1.461331) q[0];
rz(-1.4029067) q[1];
sx q[1];
rz(-2.1673514) q[1];
sx q[1];
rz(3.0775552) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.94992276) q[0];
sx q[0];
rz(-2.4987698) q[0];
sx q[0];
rz(-1.886801) q[0];
rz(-pi) q[1];
rz(1.0520347) q[2];
sx q[2];
rz(-0.56781893) q[2];
sx q[2];
rz(-2.0702814) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.36754164) q[1];
sx q[1];
rz(-1.6072175) q[1];
sx q[1];
rz(0.35518412) q[1];
x q[2];
rz(-1.8827581) q[3];
sx q[3];
rz(-1.5332111) q[3];
sx q[3];
rz(3.104044) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(3.0503851) q[2];
sx q[2];
rz(-0.63082266) q[2];
sx q[2];
rz(-1.5554265) q[2];
rz(0.88820109) q[3];
sx q[3];
rz(-1.2398088) q[3];
sx q[3];
rz(0.92938882) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.14426194) q[0];
sx q[0];
rz(-1.0634796) q[0];
sx q[0];
rz(-1.7096747) q[0];
rz(-0.56888467) q[1];
sx q[1];
rz(-2.6061997) q[1];
sx q[1];
rz(-1.127839) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9040684) q[0];
sx q[0];
rz(-1.5868109) q[0];
sx q[0];
rz(-2.9928656) q[0];
rz(-pi) q[1];
x q[1];
rz(0.88840975) q[2];
sx q[2];
rz(-2.3896304) q[2];
sx q[2];
rz(2.106452) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.7334426) q[1];
sx q[1];
rz(-0.24194716) q[1];
sx q[1];
rz(-1.5223632) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0206251) q[3];
sx q[3];
rz(-0.63311011) q[3];
sx q[3];
rz(1.8819295) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.52788064) q[2];
sx q[2];
rz(-0.85313672) q[2];
sx q[2];
rz(0.17364994) q[2];
rz(0.33637834) q[3];
sx q[3];
rz(-1.222638) q[3];
sx q[3];
rz(-0.39150795) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4353452) q[0];
sx q[0];
rz(-0.58723891) q[0];
sx q[0];
rz(-1.4655112) q[0];
rz(2.3174875) q[1];
sx q[1];
rz(-1.5287639) q[1];
sx q[1];
rz(0.5724268) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3819645) q[0];
sx q[0];
rz(-2.1158764) q[0];
sx q[0];
rz(2.1419924) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0881537) q[2];
sx q[2];
rz(-2.6247019) q[2];
sx q[2];
rz(0.042339485) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.6483557) q[1];
sx q[1];
rz(-1.3322468) q[1];
sx q[1];
rz(-1.7931213) q[1];
x q[2];
rz(-2.8045373) q[3];
sx q[3];
rz(-0.73738499) q[3];
sx q[3];
rz(0.22102236) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.77999014) q[2];
sx q[2];
rz(-0.47452351) q[2];
sx q[2];
rz(-0.53722107) q[2];
rz(-1.0572664) q[3];
sx q[3];
rz(-2.2500762) q[3];
sx q[3];
rz(-0.69361544) q[3];
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
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
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
rz(-0.097975227) q[2];
sx q[2];
rz(-0.91839472) q[2];
sx q[2];
rz(0.67994946) q[2];
rz(-1.0827071) q[3];
sx q[3];
rz(-1.516021) q[3];
sx q[3];
rz(1.0513023) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
