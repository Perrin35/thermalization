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
rz(1.7841568) q[0];
sx q[0];
rz(-0.62499243) q[0];
sx q[0];
rz(-1.619119) q[0];
rz(-1.9982665) q[1];
sx q[1];
rz(-0.63653094) q[1];
sx q[1];
rz(1.358939) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9516884) q[0];
sx q[0];
rz(-1.0190599) q[0];
sx q[0];
rz(2.8008528) q[0];
rz(2.8463616) q[2];
sx q[2];
rz(-1.9257987) q[2];
sx q[2];
rz(0.4685185) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.2083019) q[1];
sx q[1];
rz(-1.8027026) q[1];
sx q[1];
rz(1.8342706) q[1];
rz(-pi) q[2];
rz(-1.3728849) q[3];
sx q[3];
rz(-2.6555914) q[3];
sx q[3];
rz(-0.81102351) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.4437359) q[2];
sx q[2];
rz(-0.82252684) q[2];
sx q[2];
rz(2.6653384) q[2];
rz(-0.13992986) q[3];
sx q[3];
rz(-1.6401451) q[3];
sx q[3];
rz(-3.068315) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(-0.78854617) q[0];
sx q[0];
rz(-1.6206425) q[0];
sx q[0];
rz(-0.36343685) q[0];
rz(0.67047554) q[1];
sx q[1];
rz(-1.8164219) q[1];
sx q[1];
rz(-2.061981) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.99936679) q[0];
sx q[0];
rz(-1.8086664) q[0];
sx q[0];
rz(-1.0269985) q[0];
rz(-pi) q[1];
rz(-2.1035467) q[2];
sx q[2];
rz(-1.6831301) q[2];
sx q[2];
rz(2.7137604) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.9448589) q[1];
sx q[1];
rz(-2.8962039) q[1];
sx q[1];
rz(1.3887482) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7538746) q[3];
sx q[3];
rz(-0.61655515) q[3];
sx q[3];
rz(1.316837) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.67538992) q[2];
sx q[2];
rz(-0.73489302) q[2];
sx q[2];
rz(2.7030763) q[2];
rz(-1.0112666) q[3];
sx q[3];
rz(-2.4562685) q[3];
sx q[3];
rz(1.4810168) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(0.48280516) q[0];
sx q[0];
rz(-2.6130982) q[0];
sx q[0];
rz(-0.82861376) q[0];
rz(1.8939182) q[1];
sx q[1];
rz(-1.9745461) q[1];
sx q[1];
rz(0.25965986) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.22406507) q[0];
sx q[0];
rz(-2.4026786) q[0];
sx q[0];
rz(1.1631484) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9405878) q[2];
sx q[2];
rz(-1.2134873) q[2];
sx q[2];
rz(1.0543038) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.5071755) q[1];
sx q[1];
rz(-2.2440662) q[1];
sx q[1];
rz(-2.5310764) q[1];
rz(-1.1961706) q[3];
sx q[3];
rz(-0.68794477) q[3];
sx q[3];
rz(1.7619531) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.9344249) q[2];
sx q[2];
rz(-1.5639037) q[2];
sx q[2];
rz(2.3210607) q[2];
rz(0.82646838) q[3];
sx q[3];
rz(-0.99244899) q[3];
sx q[3];
rz(-1.8291399) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8823223) q[0];
sx q[0];
rz(-0.81040183) q[0];
sx q[0];
rz(-2.209254) q[0];
rz(-1.3177634) q[1];
sx q[1];
rz(-1.980314) q[1];
sx q[1];
rz(3.0188149) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2185635) q[0];
sx q[0];
rz(-1.6549131) q[0];
sx q[0];
rz(-0.53972678) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8343646) q[2];
sx q[2];
rz(-0.84103742) q[2];
sx q[2];
rz(2.407069) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.1689366) q[1];
sx q[1];
rz(-2.3984809) q[1];
sx q[1];
rz(-1.572322) q[1];
rz(-pi) q[2];
rz(-0.56468954) q[3];
sx q[3];
rz(-1.3211234) q[3];
sx q[3];
rz(-2.3752729) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.5968898) q[2];
sx q[2];
rz(-2.2201316) q[2];
sx q[2];
rz(-0.60453647) q[2];
rz(0.70945159) q[3];
sx q[3];
rz(-0.76990288) q[3];
sx q[3];
rz(-2.3179222) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9619047) q[0];
sx q[0];
rz(-0.12374319) q[0];
sx q[0];
rz(1.7748348) q[0];
rz(0.99864117) q[1];
sx q[1];
rz(-0.81592453) q[1];
sx q[1];
rz(0.5407812) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3872469) q[0];
sx q[0];
rz(-1.7779113) q[0];
sx q[0];
rz(-1.4282754) q[0];
x q[1];
rz(0.86414637) q[2];
sx q[2];
rz(-0.96071488) q[2];
sx q[2];
rz(-2.6928201) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.04356727) q[1];
sx q[1];
rz(-0.35570733) q[1];
sx q[1];
rz(-0.58009899) q[1];
rz(1.389796) q[3];
sx q[3];
rz(-2.1297526) q[3];
sx q[3];
rz(1.4829092) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.6775386) q[2];
sx q[2];
rz(-2.3276734) q[2];
sx q[2];
rz(0.82826725) q[2];
rz(-1.6895435) q[3];
sx q[3];
rz(-1.4938846) q[3];
sx q[3];
rz(-0.90788466) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.117711) q[0];
sx q[0];
rz(-1.1957059) q[0];
sx q[0];
rz(-1.1718132) q[0];
rz(2.4614629) q[1];
sx q[1];
rz(-1.2260194) q[1];
sx q[1];
rz(-2.5853058) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2883462) q[0];
sx q[0];
rz(-2.6814406) q[0];
sx q[0];
rz(2.1723738) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7173217) q[2];
sx q[2];
rz(-1.5735448) q[2];
sx q[2];
rz(-0.78080746) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.7032675) q[1];
sx q[1];
rz(-2.5653953) q[1];
sx q[1];
rz(-3.0238815) q[1];
rz(-pi) q[2];
rz(1.0917967) q[3];
sx q[3];
rz(-2.1508039) q[3];
sx q[3];
rz(-2.403808) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.80643001) q[2];
sx q[2];
rz(-1.282629) q[2];
sx q[2];
rz(2.9365017) q[2];
rz(-2.3388376) q[3];
sx q[3];
rz(-1.1057248) q[3];
sx q[3];
rz(0.55454379) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.71451181) q[0];
sx q[0];
rz(-1.0733805) q[0];
sx q[0];
rz(2.9141973) q[0];
rz(2.3511476) q[1];
sx q[1];
rz(-0.77722725) q[1];
sx q[1];
rz(-2.3048185) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6444708) q[0];
sx q[0];
rz(-0.56915149) q[0];
sx q[0];
rz(-2.7272124) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6322962) q[2];
sx q[2];
rz(-2.4894538) q[2];
sx q[2];
rz(-0.76667128) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.5414801) q[1];
sx q[1];
rz(-1.2573693) q[1];
sx q[1];
rz(2.7164943) q[1];
rz(2.7526189) q[3];
sx q[3];
rz(-0.79998652) q[3];
sx q[3];
rz(0.98158132) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.5493912) q[2];
sx q[2];
rz(-1.6624007) q[2];
sx q[2];
rz(2.5679892) q[2];
rz(0.32522374) q[3];
sx q[3];
rz(-0.89541382) q[3];
sx q[3];
rz(-0.4044683) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10930546) q[0];
sx q[0];
rz(-2.9635297) q[0];
sx q[0];
rz(-2.4626515) q[0];
rz(0.13488723) q[1];
sx q[1];
rz(-2.0811681) q[1];
sx q[1];
rz(1.4237684) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7622691) q[0];
sx q[0];
rz(-2.0998619) q[0];
sx q[0];
rz(0.53405116) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.55439722) q[2];
sx q[2];
rz(-0.94292484) q[2];
sx q[2];
rz(0.45931057) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.0404741) q[1];
sx q[1];
rz(-0.76381874) q[1];
sx q[1];
rz(-1.5274672) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1072949) q[3];
sx q[3];
rz(-2.6141659) q[3];
sx q[3];
rz(0.81125439) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.5474825) q[2];
sx q[2];
rz(-2.2308733) q[2];
sx q[2];
rz(2.9098848) q[2];
rz(1.0101275) q[3];
sx q[3];
rz(-1.908327) q[3];
sx q[3];
rz(2.2505984) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6758839) q[0];
sx q[0];
rz(-0.75513419) q[0];
sx q[0];
rz(2.6278507) q[0];
rz(-1.7087917) q[1];
sx q[1];
rz(-1.2756196) q[1];
sx q[1];
rz(-1.5076216) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5861478) q[0];
sx q[0];
rz(-2.9990104) q[0];
sx q[0];
rz(-0.34839387) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5812923) q[2];
sx q[2];
rz(-2.8447084) q[2];
sx q[2];
rz(2.1765944) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.3097174) q[1];
sx q[1];
rz(-2.9621705) q[1];
sx q[1];
rz(-0.76283263) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8225372) q[3];
sx q[3];
rz(-1.5466403) q[3];
sx q[3];
rz(-2.5509953) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.0214009) q[2];
sx q[2];
rz(-0.84453619) q[2];
sx q[2];
rz(2.8821442) q[2];
rz(-1.9868959) q[3];
sx q[3];
rz(-2.3754933) q[3];
sx q[3];
rz(-1.8416789) q[3];
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
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4907289) q[0];
sx q[0];
rz(-1.6614953) q[0];
sx q[0];
rz(-1.4890626) q[0];
rz(-2.5015855) q[1];
sx q[1];
rz(-2.044544) q[1];
sx q[1];
rz(-2.1515813) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31371597) q[0];
sx q[0];
rz(-1.4826688) q[0];
sx q[0];
rz(-2.9966596) q[0];
x q[1];
rz(2.3732164) q[2];
sx q[2];
rz(-1.2207165) q[2];
sx q[2];
rz(-0.96736747) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.4943018) q[1];
sx q[1];
rz(-1.9875981) q[1];
sx q[1];
rz(-0.17636645) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3481989) q[3];
sx q[3];
rz(-1.5951459) q[3];
sx q[3];
rz(-0.81741945) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.2085569) q[2];
sx q[2];
rz(-1.2988337) q[2];
sx q[2];
rz(-2.9162858) q[2];
rz(-2.2202668) q[3];
sx q[3];
rz(-2.7511629) q[3];
sx q[3];
rz(3.0909753) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(3.0113572) q[0];
sx q[0];
rz(-0.36778944) q[0];
sx q[0];
rz(0.854048) q[0];
rz(1.3182974) q[1];
sx q[1];
rz(-1.5250991) q[1];
sx q[1];
rz(-0.93223882) q[1];
rz(-1.83056) q[2];
sx q[2];
rz(-1.6648035) q[2];
sx q[2];
rz(1.8307508) q[2];
rz(2.1193567) q[3];
sx q[3];
rz(-1.4132186) q[3];
sx q[3];
rz(1.3506387) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
