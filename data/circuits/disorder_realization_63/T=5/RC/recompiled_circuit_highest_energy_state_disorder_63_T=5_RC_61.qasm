OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.0492078) q[0];
sx q[0];
rz(7.0424289) q[0];
sx q[0];
rz(9.0483604) q[0];
rz(-1.5837826) q[1];
sx q[1];
rz(4.2111068) q[1];
sx q[1];
rz(11.61011) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.21951708) q[0];
sx q[0];
rz(-1.886036) q[0];
sx q[0];
rz(-2.2275777) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.42266385) q[2];
sx q[2];
rz(-1.7476805) q[2];
sx q[2];
rz(1.6852578) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.53102436) q[1];
sx q[1];
rz(-0.7501561) q[1];
sx q[1];
rz(-1.3917188) q[1];
rz(-pi) q[2];
rz(1.3322863) q[3];
sx q[3];
rz(-1.8347036) q[3];
sx q[3];
rz(-1.1562386) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.3236397) q[2];
sx q[2];
rz(-2.2223667) q[2];
sx q[2];
rz(1.2032571) q[2];
rz(1.0407) q[3];
sx q[3];
rz(-2.2668656) q[3];
sx q[3];
rz(-0.11999764) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0284718) q[0];
sx q[0];
rz(-2.5682243) q[0];
sx q[0];
rz(-2.0290802) q[0];
rz(1.6587229) q[1];
sx q[1];
rz(-0.590938) q[1];
sx q[1];
rz(0.15731752) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.33799115) q[0];
sx q[0];
rz(-0.71283075) q[0];
sx q[0];
rz(2.8719851) q[0];
x q[1];
rz(0.21944616) q[2];
sx q[2];
rz(-2.7601961) q[2];
sx q[2];
rz(1.4666605) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.78029481) q[1];
sx q[1];
rz(-2.2289702) q[1];
sx q[1];
rz(2.8315413) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8202845) q[3];
sx q[3];
rz(-2.9071593) q[3];
sx q[3];
rz(-0.050416273) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.6436254) q[2];
sx q[2];
rz(-2.6622541) q[2];
sx q[2];
rz(0.41110006) q[2];
rz(0.57605612) q[3];
sx q[3];
rz(-2.1274302) q[3];
sx q[3];
rz(-0.99803734) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3882947) q[0];
sx q[0];
rz(-1.0030712) q[0];
sx q[0];
rz(0.74111795) q[0];
rz(1.4499715) q[1];
sx q[1];
rz(-1.8284109) q[1];
sx q[1];
rz(2.5695739) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.30878371) q[0];
sx q[0];
rz(-2.7833412) q[0];
sx q[0];
rz(-2.1705296) q[0];
rz(-pi) q[1];
rz(-2.4837844) q[2];
sx q[2];
rz(-0.33051046) q[2];
sx q[2];
rz(-3.0140439) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.3254889) q[1];
sx q[1];
rz(-1.1469335) q[1];
sx q[1];
rz(0.44386379) q[1];
x q[2];
rz(-0.97870632) q[3];
sx q[3];
rz(-0.39798576) q[3];
sx q[3];
rz(0.0090816895) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.584562) q[2];
sx q[2];
rz(-1.9494373) q[2];
sx q[2];
rz(-0.7977879) q[2];
rz(-1.3095464) q[3];
sx q[3];
rz(-1.8833501) q[3];
sx q[3];
rz(-1.7358739) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0107467) q[0];
sx q[0];
rz(-2.2444785) q[0];
sx q[0];
rz(2.77453) q[0];
rz(1.9233854) q[1];
sx q[1];
rz(-1.4116849) q[1];
sx q[1];
rz(-2.9201115) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9629134) q[0];
sx q[0];
rz(-1.7811462) q[0];
sx q[0];
rz(1.5014929) q[0];
rz(-pi) q[1];
x q[1];
rz(0.367737) q[2];
sx q[2];
rz(-1.6984816) q[2];
sx q[2];
rz(0.62386306) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.023497009) q[1];
sx q[1];
rz(-2.6129236) q[1];
sx q[1];
rz(0.65331991) q[1];
x q[2];
rz(-0.52677299) q[3];
sx q[3];
rz(-1.793141) q[3];
sx q[3];
rz(2.6759345) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.14747846) q[2];
sx q[2];
rz(-1.1722112) q[2];
sx q[2];
rz(1.3383024) q[2];
rz(-2.639468) q[3];
sx q[3];
rz(-0.10052557) q[3];
sx q[3];
rz(0.24371915) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62622708) q[0];
sx q[0];
rz(-0.60972649) q[0];
sx q[0];
rz(-1.5355661) q[0];
rz(-0.045479927) q[1];
sx q[1];
rz(-1.3810424) q[1];
sx q[1];
rz(2.6754726) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.202318) q[0];
sx q[0];
rz(-2.8554248) q[0];
sx q[0];
rz(1.0658468) q[0];
rz(-0.32987288) q[2];
sx q[2];
rz(-0.95392841) q[2];
sx q[2];
rz(0.62451476) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.5369617) q[1];
sx q[1];
rz(-1.5535181) q[1];
sx q[1];
rz(-3.0503154) q[1];
rz(-pi) q[2];
rz(-0.97753559) q[3];
sx q[3];
rz(-2.1338507) q[3];
sx q[3];
rz(-2.7138591) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.8972299) q[2];
sx q[2];
rz(-1.676061) q[2];
sx q[2];
rz(2.3587295) q[2];
rz(0.016294567) q[3];
sx q[3];
rz(-1.8330845) q[3];
sx q[3];
rz(-1.0619987) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8472854) q[0];
sx q[0];
rz(-0.48536244) q[0];
sx q[0];
rz(2.7622727) q[0];
rz(-2.8324221) q[1];
sx q[1];
rz(-1.199017) q[1];
sx q[1];
rz(0.22383037) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4113385) q[0];
sx q[0];
rz(-0.77654949) q[0];
sx q[0];
rz(2.0221439) q[0];
rz(1.0439453) q[2];
sx q[2];
rz(-2.9211524) q[2];
sx q[2];
rz(-1.8048087) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.47403112) q[1];
sx q[1];
rz(-1.4646475) q[1];
sx q[1];
rz(0.1072704) q[1];
rz(-0.70586127) q[3];
sx q[3];
rz(-2.0076723) q[3];
sx q[3];
rz(1.0918416) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.0763187) q[2];
sx q[2];
rz(-1.167401) q[2];
sx q[2];
rz(-0.25263146) q[2];
rz(-0.97918716) q[3];
sx q[3];
rz(-0.88117176) q[3];
sx q[3];
rz(0.34669909) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.52794367) q[0];
sx q[0];
rz(-0.76512965) q[0];
sx q[0];
rz(-1.2921523) q[0];
rz(2.6076803) q[1];
sx q[1];
rz(-1.4521867) q[1];
sx q[1];
rz(2.7814878) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.502141) q[0];
sx q[0];
rz(-1.2107673) q[0];
sx q[0];
rz(-3.0675824) q[0];
x q[1];
rz(-2.656237) q[2];
sx q[2];
rz(-1.2872641) q[2];
sx q[2];
rz(-2.527371) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.26038995) q[1];
sx q[1];
rz(-0.61169988) q[1];
sx q[1];
rz(-2.7124497) q[1];
rz(-pi) q[2];
rz(1.4680181) q[3];
sx q[3];
rz(-1.53781) q[3];
sx q[3];
rz(-1.7793293) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.1503633) q[2];
sx q[2];
rz(-1.5183307) q[2];
sx q[2];
rz(-0.72597996) q[2];
rz(2.7109072) q[3];
sx q[3];
rz(-2.3613598) q[3];
sx q[3];
rz(-0.45346692) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.78366572) q[0];
sx q[0];
rz(-1.6538606) q[0];
sx q[0];
rz(-2.3554392) q[0];
rz(0.17732009) q[1];
sx q[1];
rz(-2.0568078) q[1];
sx q[1];
rz(-2.1254983) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.02064146) q[0];
sx q[0];
rz(-0.30123152) q[0];
sx q[0];
rz(1.2094638) q[0];
rz(-2.9295278) q[2];
sx q[2];
rz(-1.4432505) q[2];
sx q[2];
rz(0.86033347) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.0018277) q[1];
sx q[1];
rz(-1.5878126) q[1];
sx q[1];
rz(-0.66405343) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.75288624) q[3];
sx q[3];
rz(-1.8532003) q[3];
sx q[3];
rz(1.2644067) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.5440386) q[2];
sx q[2];
rz(-0.95557135) q[2];
sx q[2];
rz(-2.3717144) q[2];
rz(0.17635135) q[3];
sx q[3];
rz(-1.6917112) q[3];
sx q[3];
rz(0.31340733) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65799323) q[0];
sx q[0];
rz(-0.083881065) q[0];
sx q[0];
rz(-1.0461079) q[0];
rz(-1.5484035) q[1];
sx q[1];
rz(-1.7240588) q[1];
sx q[1];
rz(-1.9200602) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0964805) q[0];
sx q[0];
rz(-1.6754221) q[0];
sx q[0];
rz(-1.2136995) q[0];
rz(2.9905926) q[2];
sx q[2];
rz(-1.8859474) q[2];
sx q[2];
rz(-1.1282819) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.56961774) q[1];
sx q[1];
rz(-1.3571661) q[1];
sx q[1];
rz(-2.9353082) q[1];
x q[2];
rz(-3.0296589) q[3];
sx q[3];
rz(-3.0977664) q[3];
sx q[3];
rz(-2.6153713) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.017642411) q[2];
sx q[2];
rz(-1.7240179) q[2];
sx q[2];
rz(1.9564015) q[2];
rz(-2.8017398) q[3];
sx q[3];
rz(-0.9797107) q[3];
sx q[3];
rz(3.0237696) q[3];
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
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3621984) q[0];
sx q[0];
rz(-1.3267936) q[0];
sx q[0];
rz(-2.4438044) q[0];
rz(2.4367874) q[1];
sx q[1];
rz(-2.0185399) q[1];
sx q[1];
rz(0.25308457) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11624434) q[0];
sx q[0];
rz(-1.3379813) q[0];
sx q[0];
rz(-0.076007621) q[0];
x q[1];
rz(-2.4867851) q[2];
sx q[2];
rz(-1.9960446) q[2];
sx q[2];
rz(-3.0523093) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.733787) q[1];
sx q[1];
rz(-1.1152667) q[1];
sx q[1];
rz(-2.2940919) q[1];
rz(-pi) q[2];
rz(-0.28530085) q[3];
sx q[3];
rz(-2.3904475) q[3];
sx q[3];
rz(1.4073262) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.42085984) q[2];
sx q[2];
rz(-1.3497817) q[2];
sx q[2];
rz(-0.68812686) q[2];
rz(0.94528919) q[3];
sx q[3];
rz(-1.1752081) q[3];
sx q[3];
rz(-0.92760408) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-2.5505493) q[0];
sx q[0];
rz(-1.3591546) q[0];
sx q[0];
rz(-1.2426283) q[0];
rz(0.91724829) q[1];
sx q[1];
rz(-1.4731673) q[1];
sx q[1];
rz(0.57631667) q[1];
rz(1.2082214) q[2];
sx q[2];
rz(-2.1167663) q[2];
sx q[2];
rz(2.703985) q[2];
rz(-1.1444202) q[3];
sx q[3];
rz(-1.4833428) q[3];
sx q[3];
rz(-1.1946027) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
