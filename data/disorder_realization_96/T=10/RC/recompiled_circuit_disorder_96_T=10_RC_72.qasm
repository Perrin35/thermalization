OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-3.0918026) q[0];
sx q[0];
rz(-3.0135305) q[0];
sx q[0];
rz(-0.81737104) q[0];
rz(0.983239) q[1];
sx q[1];
rz(2.6020738) q[1];
sx q[1];
rz(10.625216) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6450206) q[0];
sx q[0];
rz(-1.8121769) q[0];
sx q[0];
rz(2.7215331) q[0];
rz(0.13237662) q[2];
sx q[2];
rz(-1.0356324) q[2];
sx q[2];
rz(-1.3995427) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.5719205) q[1];
sx q[1];
rz(-1.8351646) q[1];
sx q[1];
rz(2.4155248) q[1];
x q[2];
rz(1.8957542) q[3];
sx q[3];
rz(-2.339139) q[3];
sx q[3];
rz(-1.5047531) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.1203221) q[2];
sx q[2];
rz(-0.56818429) q[2];
sx q[2];
rz(1.5585287) q[2];
rz(-2.1448686) q[3];
sx q[3];
rz(-0.45209) q[3];
sx q[3];
rz(0.42580095) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5581756) q[0];
sx q[0];
rz(-2.3893864) q[0];
sx q[0];
rz(0.054071991) q[0];
rz(1.1955098) q[1];
sx q[1];
rz(-1.0369438) q[1];
sx q[1];
rz(0.53584677) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8230096) q[0];
sx q[0];
rz(-0.59360015) q[0];
sx q[0];
rz(1.7943322) q[0];
rz(-pi) q[1];
x q[1];
rz(0.29157721) q[2];
sx q[2];
rz(-0.73086408) q[2];
sx q[2];
rz(-0.72327327) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.6260813) q[1];
sx q[1];
rz(-1.1647845) q[1];
sx q[1];
rz(-1.2938234) q[1];
rz(-pi) q[2];
rz(1.7934947) q[3];
sx q[3];
rz(-2.9950812) q[3];
sx q[3];
rz(0.033586249) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.1198931) q[2];
sx q[2];
rz(-1.0571486) q[2];
sx q[2];
rz(-2.1832441) q[2];
rz(3.0751394) q[3];
sx q[3];
rz(-1.5840014) q[3];
sx q[3];
rz(2.691793) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7217343) q[0];
sx q[0];
rz(-1.9323213) q[0];
sx q[0];
rz(-0.15047519) q[0];
rz(-0.45723215) q[1];
sx q[1];
rz(-0.91232863) q[1];
sx q[1];
rz(3.1157852) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.74224) q[0];
sx q[0];
rz(-0.25063801) q[0];
sx q[0];
rz(1.4901194) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.089695887) q[2];
sx q[2];
rz(-1.2604453) q[2];
sx q[2];
rz(-2.2657564) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.4986213) q[1];
sx q[1];
rz(-1.9830623) q[1];
sx q[1];
rz(2.4436414) q[1];
x q[2];
rz(1.2529536) q[3];
sx q[3];
rz(-0.40099537) q[3];
sx q[3];
rz(2.3272115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0187443) q[2];
sx q[2];
rz(-0.3652502) q[2];
sx q[2];
rz(-0.5853816) q[2];
rz(2.9600926) q[3];
sx q[3];
rz(-1.8094962) q[3];
sx q[3];
rz(1.5766778) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(2.240775) q[0];
sx q[0];
rz(-0.62269354) q[0];
sx q[0];
rz(0.17661072) q[0];
rz(-2.2606842) q[1];
sx q[1];
rz(-2.0753588) q[1];
sx q[1];
rz(-0.53612971) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49318424) q[0];
sx q[0];
rz(-0.83183653) q[0];
sx q[0];
rz(-1.0113082) q[0];
rz(-pi) q[1];
rz(-2.9085607) q[2];
sx q[2];
rz(-2.8985902) q[2];
sx q[2];
rz(-2.255893) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.94169468) q[1];
sx q[1];
rz(-0.75921339) q[1];
sx q[1];
rz(0.56337507) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2767775) q[3];
sx q[3];
rz(-2.2259568) q[3];
sx q[3];
rz(-2.7340207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.6716016) q[2];
sx q[2];
rz(-1.4208379) q[2];
sx q[2];
rz(-2.0969351) q[2];
rz(-2.4345543) q[3];
sx q[3];
rz(-2.1576594) q[3];
sx q[3];
rz(-2.7846591) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9064643) q[0];
sx q[0];
rz(-1.3236073) q[0];
sx q[0];
rz(1.0513603) q[0];
rz(-1.4936739) q[1];
sx q[1];
rz(-2.5525679) q[1];
sx q[1];
rz(-0.043118127) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9263822) q[0];
sx q[0];
rz(-2.1676817) q[0];
sx q[0];
rz(-0.2290639) q[0];
rz(1.0518603) q[2];
sx q[2];
rz(-1.8199925) q[2];
sx q[2];
rz(-0.11617004) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.4700714) q[1];
sx q[1];
rz(-2.7579691) q[1];
sx q[1];
rz(2.3684711) q[1];
rz(-pi) q[2];
x q[2];
rz(2.99302) q[3];
sx q[3];
rz(-0.79509495) q[3];
sx q[3];
rz(-3.1032004) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.12895) q[2];
sx q[2];
rz(-0.49762112) q[2];
sx q[2];
rz(2.7894003) q[2];
rz(-0.59018618) q[3];
sx q[3];
rz(-0.47368172) q[3];
sx q[3];
rz(-0.56110704) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-2.6234289) q[0];
sx q[0];
rz(-1.2151027) q[0];
sx q[0];
rz(1.1556926) q[0];
rz(2.39134) q[1];
sx q[1];
rz(-2.2022088) q[1];
sx q[1];
rz(1.0587143) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6569865) q[0];
sx q[0];
rz(-1.779042) q[0];
sx q[0];
rz(0.44021846) q[0];
rz(-pi) q[1];
x q[1];
rz(0.58745678) q[2];
sx q[2];
rz(-2.537478) q[2];
sx q[2];
rz(-2.6641252) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.25710479) q[1];
sx q[1];
rz(-0.37933644) q[1];
sx q[1];
rz(2.2124955) q[1];
x q[2];
rz(-0.70763564) q[3];
sx q[3];
rz(-1.71873) q[3];
sx q[3];
rz(1.9572452) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.47026149) q[2];
sx q[2];
rz(-1.7241717) q[2];
sx q[2];
rz(1.2188101) q[2];
rz(-1.9865215) q[3];
sx q[3];
rz(-2.9030436) q[3];
sx q[3];
rz(1.7003805) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0999775) q[0];
sx q[0];
rz(-1.8247373) q[0];
sx q[0];
rz(-1.6301427) q[0];
rz(-1.3776243) q[1];
sx q[1];
rz(-2.8306077) q[1];
sx q[1];
rz(0.84164936) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3232988) q[0];
sx q[0];
rz(-0.57384402) q[0];
sx q[0];
rz(-2.7159575) q[0];
rz(-pi) q[1];
rz(0.26288962) q[2];
sx q[2];
rz(-0.6436231) q[2];
sx q[2];
rz(-0.4292683) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.51534286) q[1];
sx q[1];
rz(-1.5222349) q[1];
sx q[1];
rz(-1.4944782) q[1];
x q[2];
rz(0.21861403) q[3];
sx q[3];
rz(-0.84158763) q[3];
sx q[3];
rz(1.5900172) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.001361751) q[2];
sx q[2];
rz(-1.7732239) q[2];
sx q[2];
rz(3.0916396) q[2];
rz(2.4800381) q[3];
sx q[3];
rz(-2.619132) q[3];
sx q[3];
rz(2.9522827) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2426303) q[0];
sx q[0];
rz(-1.0268509) q[0];
sx q[0];
rz(-1.7393973) q[0];
rz(-3.0461123) q[1];
sx q[1];
rz(-1.9752558) q[1];
sx q[1];
rz(-2.7239674) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2554889) q[0];
sx q[0];
rz(-2.4574453) q[0];
sx q[0];
rz(-0.72649254) q[0];
rz(-pi) q[1];
rz(-1.2049963) q[2];
sx q[2];
rz(-1.6779643) q[2];
sx q[2];
rz(1.1366213) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.3216746) q[1];
sx q[1];
rz(-2.3707317) q[1];
sx q[1];
rz(-2.2714771) q[1];
x q[2];
rz(0.91248625) q[3];
sx q[3];
rz(-2.6819326) q[3];
sx q[3];
rz(1.0115136) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.79545704) q[2];
sx q[2];
rz(-1.4435578) q[2];
sx q[2];
rz(-2.0987089) q[2];
rz(2.4677094) q[3];
sx q[3];
rz(-1.4833114) q[3];
sx q[3];
rz(-2.2414482) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3089356) q[0];
sx q[0];
rz(-1.4082264) q[0];
sx q[0];
rz(2.5119264) q[0];
rz(-0.57485238) q[1];
sx q[1];
rz(-1.31153) q[1];
sx q[1];
rz(-0.94690698) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4862343) q[0];
sx q[0];
rz(-1.2709193) q[0];
sx q[0];
rz(-0.34098682) q[0];
rz(-pi) q[1];
rz(0.097054585) q[2];
sx q[2];
rz(-1.5311054) q[2];
sx q[2];
rz(1.5511712) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.8967594) q[1];
sx q[1];
rz(-1.4866801) q[1];
sx q[1];
rz(2.8068078) q[1];
x q[2];
rz(2.9506748) q[3];
sx q[3];
rz(-2.0058161) q[3];
sx q[3];
rz(-0.85489475) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.56069121) q[2];
sx q[2];
rz(-1.1297444) q[2];
sx q[2];
rz(1.8927195) q[2];
rz(2.4272264) q[3];
sx q[3];
rz(-1.2759821) q[3];
sx q[3];
rz(-0.26556382) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(-2.0666075) q[0];
sx q[0];
rz(-2.5191436) q[0];
sx q[0];
rz(2.4865436) q[0];
rz(0.89637268) q[1];
sx q[1];
rz(-0.90677774) q[1];
sx q[1];
rz(2.4972829) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5376741) q[0];
sx q[0];
rz(-0.18225056) q[0];
sx q[0];
rz(0.3084348) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3078493) q[2];
sx q[2];
rz(-2.408228) q[2];
sx q[2];
rz(2.3761689) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.58562467) q[1];
sx q[1];
rz(-1.854419) q[1];
sx q[1];
rz(-0.3507627) q[1];
x q[2];
rz(1.0993016) q[3];
sx q[3];
rz(-1.1700556) q[3];
sx q[3];
rz(0.61939643) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.7908988) q[2];
sx q[2];
rz(-1.6692946) q[2];
sx q[2];
rz(2.541686) q[2];
rz(-0.89896262) q[3];
sx q[3];
rz(-0.18342429) q[3];
sx q[3];
rz(1.3658587) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.29466378) q[0];
sx q[0];
rz(-1.8142721) q[0];
sx q[0];
rz(-0.66692764) q[0];
rz(0.22944336) q[1];
sx q[1];
rz(-0.89090092) q[1];
sx q[1];
rz(0.13577239) q[1];
rz(-1.3953801) q[2];
sx q[2];
rz(-2.5098364) q[2];
sx q[2];
rz(-2.9927158) q[2];
rz(-2.2453528) q[3];
sx q[3];
rz(-1.2590209) q[3];
sx q[3];
rz(-0.38929064) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
