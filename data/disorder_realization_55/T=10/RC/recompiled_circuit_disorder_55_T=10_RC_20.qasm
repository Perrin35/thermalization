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
rz(-0.73524737) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.11156946) q[0];
sx q[0];
rz(-0.91696793) q[0];
sx q[0];
rz(-1.7574739) q[0];
x q[1];
rz(0.67970694) q[2];
sx q[2];
rz(-1.1628816) q[2];
sx q[2];
rz(2.9349875) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.3056065) q[1];
sx q[1];
rz(-1.7251833) q[1];
sx q[1];
rz(1.8603714) q[1];
rz(-pi) q[2];
rz(0.04282184) q[3];
sx q[3];
rz(-0.58086568) q[3];
sx q[3];
rz(-0.37561852) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.66951093) q[2];
sx q[2];
rz(-1.8410204) q[2];
sx q[2];
rz(-2.0377339) q[2];
rz(-1.8707229) q[3];
sx q[3];
rz(-1.9138252) q[3];
sx q[3];
rz(0.27403533) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1141777) q[0];
sx q[0];
rz(-0.52863055) q[0];
sx q[0];
rz(2.7052178) q[0];
rz(0.46288681) q[1];
sx q[1];
rz(-2.1040237) q[1];
sx q[1];
rz(0.26611051) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.14861815) q[0];
sx q[0];
rz(-0.87478144) q[0];
sx q[0];
rz(0.50177411) q[0];
x q[1];
rz(1.3219464) q[2];
sx q[2];
rz(-0.58983931) q[2];
sx q[2];
rz(2.313254) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.7539566) q[1];
sx q[1];
rz(-2.1093183) q[1];
sx q[1];
rz(-2.3073879) q[1];
rz(-pi) q[2];
rz(-1.3105884) q[3];
sx q[3];
rz(-0.62285138) q[3];
sx q[3];
rz(-2.7130896) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.6767072) q[2];
sx q[2];
rz(-1.8647944) q[2];
sx q[2];
rz(-0.51149386) q[2];
rz(-0.809508) q[3];
sx q[3];
rz(-1.5313238) q[3];
sx q[3];
rz(-2.8539343) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.70616102) q[0];
sx q[0];
rz(-1.5478739) q[0];
sx q[0];
rz(-2.2128552) q[0];
rz(-1.4061032) q[1];
sx q[1];
rz(-2.4415253) q[1];
sx q[1];
rz(-1.3471289) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.12362326) q[0];
sx q[0];
rz(-2.6051913) q[0];
sx q[0];
rz(-0.53039741) q[0];
rz(2.0688829) q[2];
sx q[2];
rz(-0.32215298) q[2];
sx q[2];
rz(-2.1307532) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.320142) q[1];
sx q[1];
rz(-2.4648033) q[1];
sx q[1];
rz(1.0980781) q[1];
x q[2];
rz(2.2778355) q[3];
sx q[3];
rz(-1.45544) q[3];
sx q[3];
rz(2.399721) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.9946263) q[2];
sx q[2];
rz(-1.3866084) q[2];
sx q[2];
rz(1.5220801) q[2];
rz(-2.8772723) q[3];
sx q[3];
rz(-2.1051354) q[3];
sx q[3];
rz(-0.49595293) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3494444) q[0];
sx q[0];
rz(-1.9932207) q[0];
sx q[0];
rz(2.175892) q[0];
rz(2.4194338) q[1];
sx q[1];
rz(-1.637371) q[1];
sx q[1];
rz(2.5818363) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4013195) q[0];
sx q[0];
rz(-1.6519321) q[0];
sx q[0];
rz(2.7149537) q[0];
x q[1];
rz(-0.83300029) q[2];
sx q[2];
rz(-2.2350395) q[2];
sx q[2];
rz(2.8766362) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.30448118) q[1];
sx q[1];
rz(-1.6164391) q[1];
sx q[1];
rz(-0.95506217) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.78919952) q[3];
sx q[3];
rz(-1.4025098) q[3];
sx q[3];
rz(0.60889739) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.42797783) q[2];
sx q[2];
rz(-1.5087936) q[2];
sx q[2];
rz(1.7170061) q[2];
rz(-0.26040855) q[3];
sx q[3];
rz(-1.3924761) q[3];
sx q[3];
rz(2.7105455) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.150862) q[0];
sx q[0];
rz(-1.4980415) q[0];
sx q[0];
rz(-0.36636233) q[0];
rz(1.5461961) q[1];
sx q[1];
rz(-0.55391824) q[1];
sx q[1];
rz(2.7979134) q[1];
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
x q[1];
rz(-2.5630066) q[2];
sx q[2];
rz(-2.1861976) q[2];
sx q[2];
rz(1.3784642) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.3949563) q[1];
sx q[1];
rz(-2.6379105) q[1];
sx q[1];
rz(-1.4175182) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.021275612) q[3];
sx q[3];
rz(-2.2145503) q[3];
sx q[3];
rz(1.6568041) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.7498103) q[2];
sx q[2];
rz(-0.85931531) q[2];
sx q[2];
rz(3.0878477) q[2];
rz(-1.7371477) q[3];
sx q[3];
rz(-0.497538) q[3];
sx q[3];
rz(0.29156175) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4869726) q[0];
sx q[0];
rz(-0.64990652) q[0];
sx q[0];
rz(-2.3858331) q[0];
rz(-0.02515633) q[1];
sx q[1];
rz(-0.92725602) q[1];
sx q[1];
rz(2.8818534) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1007337) q[0];
sx q[0];
rz(-1.9319527) q[0];
sx q[0];
rz(0.33613899) q[0];
x q[1];
rz(0.01979205) q[2];
sx q[2];
rz(-1.9070101) q[2];
sx q[2];
rz(2.3078231) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.7139587) q[1];
sx q[1];
rz(-0.40363388) q[1];
sx q[1];
rz(0.9637109) q[1];
rz(-pi) q[2];
rz(1.608058) q[3];
sx q[3];
rz(-2.691949) q[3];
sx q[3];
rz(2.8763308) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.48866895) q[2];
sx q[2];
rz(-2.3355464) q[2];
sx q[2];
rz(-3.0409813) q[2];
rz(2.9595024) q[3];
sx q[3];
rz(-0.87825769) q[3];
sx q[3];
rz(-1.7939059) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9942193) q[0];
sx q[0];
rz(-1.7213151) q[0];
sx q[0];
rz(0.31016645) q[0];
rz(-0.50225964) q[1];
sx q[1];
rz(-2.6175833) q[1];
sx q[1];
rz(-0.60595864) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.92354846) q[0];
sx q[0];
rz(-2.179562) q[0];
sx q[0];
rz(1.8463085) q[0];
rz(-pi) q[1];
rz(2.546054) q[2];
sx q[2];
rz(-0.34791246) q[2];
sx q[2];
rz(-2.4728342) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(3.088674) q[1];
sx q[1];
rz(-2.5175736) q[1];
sx q[1];
rz(0.96159972) q[1];
x q[2];
rz(-1.0244272) q[3];
sx q[3];
rz(-1.1064648) q[3];
sx q[3];
rz(1.6950316) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.24017748) q[2];
sx q[2];
rz(-1.1837974) q[2];
sx q[2];
rz(-2.288738) q[2];
rz(1.7715706) q[3];
sx q[3];
rz(-1.4586689) q[3];
sx q[3];
rz(-2.9366233) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
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
rz(3.0775552) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2647576) q[0];
sx q[0];
rz(-1.3834073) q[0];
sx q[0];
rz(0.95215709) q[0];
rz(-1.0520347) q[2];
sx q[2];
rz(-0.56781893) q[2];
sx q[2];
rz(-1.0713112) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.36754164) q[1];
sx q[1];
rz(-1.6072175) q[1];
sx q[1];
rz(-0.35518412) q[1];
rz(1.4488892) q[3];
sx q[3];
rz(-2.8274483) q[3];
sx q[3];
rz(-1.4172518) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.0503851) q[2];
sx q[2];
rz(-0.63082266) q[2];
sx q[2];
rz(-1.5861661) q[2];
rz(2.2533916) q[3];
sx q[3];
rz(-1.9017838) q[3];
sx q[3];
rz(-2.2122038) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.14426194) q[0];
sx q[0];
rz(-2.0781131) q[0];
sx q[0];
rz(-1.4319179) q[0];
rz(0.56888467) q[1];
sx q[1];
rz(-0.535393) q[1];
sx q[1];
rz(-1.127839) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2375243) q[0];
sx q[0];
rz(-1.5547817) q[0];
sx q[0];
rz(-2.9928656) q[0];
rz(-0.88840975) q[2];
sx q[2];
rz(-0.75196224) q[2];
sx q[2];
rz(-1.0351406) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.4081501) q[1];
sx q[1];
rz(-2.8996455) q[1];
sx q[1];
rz(-1.5223632) q[1];
rz(-pi) q[2];
rz(-0.3663775) q[3];
sx q[3];
rz(-1.0421703) q[3];
sx q[3];
rz(-2.5322994) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.52788064) q[2];
sx q[2];
rz(-2.2884559) q[2];
sx q[2];
rz(-0.17364994) q[2];
rz(-0.33637834) q[3];
sx q[3];
rz(-1.222638) q[3];
sx q[3];
rz(-2.7500847) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4353452) q[0];
sx q[0];
rz(-0.58723891) q[0];
sx q[0];
rz(1.4655112) q[0];
rz(0.82410518) q[1];
sx q[1];
rz(-1.6128287) q[1];
sx q[1];
rz(0.5724268) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6515598) q[0];
sx q[0];
rz(-2.3736554) q[0];
sx q[0];
rz(0.72816531) q[0];
x q[1];
rz(-1.1119214) q[2];
sx q[2];
rz(-1.8177114) q[2];
sx q[2];
rz(2.0725476) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.72996424) q[1];
sx q[1];
rz(-2.8169544) q[1];
sx q[1];
rz(2.4050729) q[1];
rz(-pi) q[2];
rz(1.8626067) q[3];
sx q[3];
rz(-0.88340532) q[3];
sx q[3];
rz(0.22112267) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.77999014) q[2];
sx q[2];
rz(-0.47452351) q[2];
sx q[2];
rz(-0.53722107) q[2];
rz(-2.0843263) q[3];
sx q[3];
rz(-2.2500762) q[3];
sx q[3];
rz(0.69361544) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.854241) q[0];
sx q[0];
rz(-1.1653405) q[0];
sx q[0];
rz(-1.5821138) q[0];
rz(2.6782425) q[1];
sx q[1];
rz(-0.87711038) q[1];
sx q[1];
rz(-1.6323485) q[1];
rz(-0.91607416) q[2];
sx q[2];
rz(-1.4929885) q[2];
sx q[2];
rz(2.3103466) q[2];
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
