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
rz(2.1677986) q[1];
sx q[1];
rz(-1.9314613) q[1];
sx q[1];
rz(2.4063453) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7289294) q[0];
sx q[0];
rz(-0.67617765) q[0];
sx q[0];
rz(0.23763188) q[0];
x q[1];
rz(2.4618857) q[2];
sx q[2];
rz(-1.1628816) q[2];
sx q[2];
rz(-2.9349875) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.3056065) q[1];
sx q[1];
rz(-1.7251833) q[1];
sx q[1];
rz(1.2812213) q[1];
rz(-3.0987708) q[3];
sx q[3];
rz(-0.58086568) q[3];
sx q[3];
rz(-0.37561852) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.4720817) q[2];
sx q[2];
rz(-1.3005723) q[2];
sx q[2];
rz(1.1038587) q[2];
rz(-1.8707229) q[3];
sx q[3];
rz(-1.2277675) q[3];
sx q[3];
rz(2.8675573) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0274149) q[0];
sx q[0];
rz(-0.52863055) q[0];
sx q[0];
rz(-0.43637481) q[0];
rz(-2.6787058) q[1];
sx q[1];
rz(-2.1040237) q[1];
sx q[1];
rz(0.26611051) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9929745) q[0];
sx q[0];
rz(-0.87478144) q[0];
sx q[0];
rz(-2.6398185) q[0];
x q[1];
rz(2.9782148) q[2];
sx q[2];
rz(-1.0014357) q[2];
sx q[2];
rz(0.53158224) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.38763603) q[1];
sx q[1];
rz(-1.0322744) q[1];
sx q[1];
rz(0.83420475) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3105884) q[3];
sx q[3];
rz(-0.62285138) q[3];
sx q[3];
rz(0.42850307) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.46488547) q[2];
sx q[2];
rz(-1.2767982) q[2];
sx q[2];
rz(-2.6300988) q[2];
rz(-0.809508) q[3];
sx q[3];
rz(-1.6102689) q[3];
sx q[3];
rz(-0.28765837) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4354316) q[0];
sx q[0];
rz(-1.5937188) q[0];
sx q[0];
rz(-0.92873746) q[0];
rz(-1.7354895) q[1];
sx q[1];
rz(-0.7000674) q[1];
sx q[1];
rz(1.7944638) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4192393) q[0];
sx q[0];
rz(-2.0273211) q[0];
sx q[0];
rz(1.2786352) q[0];
rz(-pi) q[1];
rz(-2.9834619) q[2];
sx q[2];
rz(-1.2889382) q[2];
sx q[2];
rz(2.6513197) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.128708) q[1];
sx q[1];
rz(-1.8599659) q[1];
sx q[1];
rz(2.1916926) q[1];
rz(-pi) q[2];
rz(0.8637572) q[3];
sx q[3];
rz(-1.45544) q[3];
sx q[3];
rz(-2.399721) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.1469664) q[2];
sx q[2];
rz(-1.7549843) q[2];
sx q[2];
rz(1.5220801) q[2];
rz(-0.26432031) q[3];
sx q[3];
rz(-2.1051354) q[3];
sx q[3];
rz(0.49595293) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3494444) q[0];
sx q[0];
rz(-1.9932207) q[0];
sx q[0];
rz(2.175892) q[0];
rz(-0.72215885) q[1];
sx q[1];
rz(-1.5042217) q[1];
sx q[1];
rz(-2.5818363) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.0069665466) q[0];
sx q[0];
rz(-2.7077733) q[0];
sx q[0];
rz(0.19402786) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3085924) q[2];
sx q[2];
rz(-0.90655316) q[2];
sx q[2];
rz(-0.26495648) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.8371115) q[1];
sx q[1];
rz(-1.5251535) q[1];
sx q[1];
rz(2.1865305) q[1];
x q[2];
rz(0.78919952) q[3];
sx q[3];
rz(-1.4025098) q[3];
sx q[3];
rz(-0.60889739) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.7136148) q[2];
sx q[2];
rz(-1.5087936) q[2];
sx q[2];
rz(-1.7170061) q[2];
rz(0.26040855) q[3];
sx q[3];
rz(-1.3924761) q[3];
sx q[3];
rz(-2.7105455) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8923963) q[0];
sx q[0];
rz(-2.2049892) q[0];
sx q[0];
rz(-2.2618494) q[0];
rz(0.57858606) q[2];
sx q[2];
rz(-2.1861976) q[2];
sx q[2];
rz(1.3784642) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.5695565) q[1];
sx q[1];
rz(-2.0680288) q[1];
sx q[1];
rz(-3.0576502) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5991391) q[3];
sx q[3];
rz(-2.4975371) q[3];
sx q[3];
rz(-1.6213662) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.7498103) q[2];
sx q[2];
rz(-2.2822773) q[2];
sx q[2];
rz(3.0878477) q[2];
rz(-1.404445) q[3];
sx q[3];
rz(-0.497538) q[3];
sx q[3];
rz(2.8500309) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
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
rz(-0.4869726) q[0];
sx q[0];
rz(-2.4916861) q[0];
sx q[0];
rz(-0.75575954) q[0];
rz(-3.1164363) q[1];
sx q[1];
rz(-0.92725602) q[1];
sx q[1];
rz(-2.8818534) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4070968) q[0];
sx q[0];
rz(-1.2571113) q[0];
sx q[0];
rz(-1.19019) q[0];
rz(-pi) q[1];
rz(1.6273646) q[2];
sx q[2];
rz(-2.8048189) q[2];
sx q[2];
rz(0.77384225) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.4252792) q[1];
sx q[1];
rz(-1.3448161) q[1];
sx q[1];
rz(1.9081566) q[1];
x q[2];
rz(-0.017976956) q[3];
sx q[3];
rz(-2.020105) q[3];
sx q[3];
rz(-0.22389212) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.48866895) q[2];
sx q[2];
rz(-0.80604625) q[2];
sx q[2];
rz(-0.10061131) q[2];
rz(2.9595024) q[3];
sx q[3];
rz(-2.263335) q[3];
sx q[3];
rz(1.7939059) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9942193) q[0];
sx q[0];
rz(-1.4202776) q[0];
sx q[0];
rz(2.8314262) q[0];
rz(-0.50225964) q[1];
sx q[1];
rz(-2.6175833) q[1];
sx q[1];
rz(-0.60595864) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.46442407) q[0];
sx q[0];
rz(-0.66093984) q[0];
sx q[0];
rz(2.7695157) q[0];
rz(-pi) q[1];
rz(-2.546054) q[2];
sx q[2];
rz(-2.7936802) q[2];
sx q[2];
rz(-2.4728342) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.1389321) q[1];
sx q[1];
rz(-1.9117038) q[1];
sx q[1];
rz(1.0374116) q[1];
rz(2.1171655) q[3];
sx q[3];
rz(-1.1064648) q[3];
sx q[3];
rz(-1.446561) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.9014152) q[2];
sx q[2];
rz(-1.1837974) q[2];
sx q[2];
rz(2.288738) q[2];
rz(1.7715706) q[3];
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
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9119499) q[0];
sx q[0];
rz(-0.59597534) q[0];
sx q[0];
rz(1.461331) q[0];
rz(1.7386859) q[1];
sx q[1];
rz(-2.1673514) q[1];
sx q[1];
rz(-0.064037474) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.56211573) q[0];
sx q[0];
rz(-0.96456438) q[0];
sx q[0];
rz(-2.9129145) q[0];
rz(0.30631752) q[2];
sx q[2];
rz(-1.0848572) q[2];
sx q[2];
rz(-1.475032) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.0361955) q[1];
sx q[1];
rz(-0.35696778) q[1];
sx q[1];
rz(-0.10439519) q[1];
rz(1.6927034) q[3];
sx q[3];
rz(-0.31414437) q[3];
sx q[3];
rz(1.7243408) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(3.0503851) q[2];
sx q[2];
rz(-0.63082266) q[2];
sx q[2];
rz(1.5861661) q[2];
rz(2.2533916) q[3];
sx q[3];
rz(-1.2398088) q[3];
sx q[3];
rz(2.2122038) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14426194) q[0];
sx q[0];
rz(-2.0781131) q[0];
sx q[0];
rz(-1.7096747) q[0];
rz(2.572708) q[1];
sx q[1];
rz(-2.6061997) q[1];
sx q[1];
rz(-1.127839) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3356714) q[0];
sx q[0];
rz(-1.7195042) q[0];
sx q[0];
rz(1.554603) q[0];
rz(-pi) q[1];
rz(2.608689) q[2];
sx q[2];
rz(-2.1295296) q[2];
sx q[2];
rz(1.8738054) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.9319218) q[1];
sx q[1];
rz(-1.5823963) q[1];
sx q[1];
rz(1.329122) q[1];
x q[2];
rz(1.0206251) q[3];
sx q[3];
rz(-0.63311011) q[3];
sx q[3];
rz(-1.2596631) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.613712) q[2];
sx q[2];
rz(-2.2884559) q[2];
sx q[2];
rz(-2.9679427) q[2];
rz(-2.8052143) q[3];
sx q[3];
rz(-1.9189546) q[3];
sx q[3];
rz(-2.7500847) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[2];
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
rz(0.82410518) q[1];
sx q[1];
rz(-1.5287639) q[1];
sx q[1];
rz(2.5691659) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.51047072) q[0];
sx q[0];
rz(-2.0513751) q[0];
sx q[0];
rz(-2.5170588) q[0];
rz(-pi) q[1];
x q[1];
rz(0.27406759) q[2];
sx q[2];
rz(-2.0147418) q[2];
sx q[2];
rz(2.519671) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.4116284) q[1];
sx q[1];
rz(-0.3246383) q[1];
sx q[1];
rz(-2.4050729) q[1];
rz(-pi) q[2];
rz(0.33705538) q[3];
sx q[3];
rz(-2.4042077) q[3];
sx q[3];
rz(2.9205703) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.3616025) q[2];
sx q[2];
rz(-2.6670691) q[2];
sx q[2];
rz(2.6043716) q[2];
rz(1.0572664) q[3];
sx q[3];
rz(-0.89151645) q[3];
sx q[3];
rz(-0.69361544) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(-1.854241) q[0];
sx q[0];
rz(-1.1653405) q[0];
sx q[0];
rz(-1.5821138) q[0];
rz(-0.46335012) q[1];
sx q[1];
rz(-0.87711038) q[1];
sx q[1];
rz(-1.6323485) q[1];
rz(-1.6981381) q[2];
sx q[2];
rz(-2.4829395) q[2];
sx q[2];
rz(0.84045835) q[2];
rz(0.061999576) q[3];
sx q[3];
rz(-2.0580895) q[3];
sx q[3];
rz(2.5930391) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];