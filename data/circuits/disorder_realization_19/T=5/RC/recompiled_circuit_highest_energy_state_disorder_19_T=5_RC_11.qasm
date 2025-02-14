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
rz(0.49993604) q[0];
sx q[0];
rz(1.9498107) q[0];
sx q[0];
rz(10.796588) q[0];
rz(-1.491188) q[1];
sx q[1];
rz(-2.2743382) q[1];
sx q[1];
rz(2.7452918) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.16184743) q[0];
sx q[0];
rz(-1.7801741) q[0];
sx q[0];
rz(-1.6518948) q[0];
rz(-pi) q[1];
rz(1.9464067) q[2];
sx q[2];
rz(-1.4109138) q[2];
sx q[2];
rz(2.218046) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.22452488) q[1];
sx q[1];
rz(-0.5366508) q[1];
sx q[1];
rz(1.7557675) q[1];
x q[2];
rz(-0.39537786) q[3];
sx q[3];
rz(-1.8413183) q[3];
sx q[3];
rz(2.6879355) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.0013915) q[2];
sx q[2];
rz(-2.3346257) q[2];
sx q[2];
rz(-1.3414471) q[2];
rz(2.0724824) q[3];
sx q[3];
rz(-0.1461229) q[3];
sx q[3];
rz(1.7599546) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.57217252) q[0];
sx q[0];
rz(-2.2241346) q[0];
sx q[0];
rz(0.65895748) q[0];
rz(3.068889) q[1];
sx q[1];
rz(-2.0854988) q[1];
sx q[1];
rz(-1.2603849) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7911071) q[0];
sx q[0];
rz(-1.1189776) q[0];
sx q[0];
rz(0.090403948) q[0];
rz(-1.057823) q[2];
sx q[2];
rz(-1.1758846) q[2];
sx q[2];
rz(2.9213443) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.94174332) q[1];
sx q[1];
rz(-1.3056271) q[1];
sx q[1];
rz(1.5035692) q[1];
rz(1.7409513) q[3];
sx q[3];
rz(-2.5239787) q[3];
sx q[3];
rz(-3.1132085) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.8947123) q[2];
sx q[2];
rz(-0.74328819) q[2];
sx q[2];
rz(-1.7903719) q[2];
rz(2.5249935) q[3];
sx q[3];
rz(-1.0298046) q[3];
sx q[3];
rz(-0.29115796) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
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
rz(-2.5093444) q[0];
sx q[0];
rz(-0.46724874) q[0];
sx q[0];
rz(-2.9314991) q[0];
rz(-3.0585152) q[1];
sx q[1];
rz(-1.2385085) q[1];
sx q[1];
rz(3.074379) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.38131648) q[0];
sx q[0];
rz(-0.018122321) q[0];
sx q[0];
rz(-2.4118547) q[0];
rz(-pi) q[1];
rz(1.1804586) q[2];
sx q[2];
rz(-1.3403041) q[2];
sx q[2];
rz(-2.0978438) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.9720502) q[1];
sx q[1];
rz(-1.2360555) q[1];
sx q[1];
rz(-2.2856345) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4806858) q[3];
sx q[3];
rz(-2.2743974) q[3];
sx q[3];
rz(1.7357602) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.82789603) q[2];
sx q[2];
rz(-1.6500429) q[2];
sx q[2];
rz(-1.776604) q[2];
rz(2.7374173) q[3];
sx q[3];
rz(-1.8242691) q[3];
sx q[3];
rz(-1.9425758) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86554027) q[0];
sx q[0];
rz(-1.4006389) q[0];
sx q[0];
rz(-2.6633967) q[0];
rz(-0.95282355) q[1];
sx q[1];
rz(-0.15004221) q[1];
sx q[1];
rz(-2.731954) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3296369) q[0];
sx q[0];
rz(-2.4853737) q[0];
sx q[0];
rz(-3.0983221) q[0];
rz(-pi) q[1];
x q[1];
rz(0.19721376) q[2];
sx q[2];
rz(-2.6773415) q[2];
sx q[2];
rz(1.0997538) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.8781902) q[1];
sx q[1];
rz(-0.95118427) q[1];
sx q[1];
rz(-0.81677572) q[1];
x q[2];
rz(1.9799819) q[3];
sx q[3];
rz(-2.7936862) q[3];
sx q[3];
rz(-2.6285841) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.5459368) q[2];
sx q[2];
rz(-1.7031534) q[2];
sx q[2];
rz(-2.0036009) q[2];
rz(2.1465837) q[3];
sx q[3];
rz(-2.5840839) q[3];
sx q[3];
rz(3.0549684) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.397641) q[0];
sx q[0];
rz(-1.3084363) q[0];
sx q[0];
rz(-2.8635136) q[0];
rz(1.2644348) q[1];
sx q[1];
rz(-1.8587298) q[1];
sx q[1];
rz(-1.5778731) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.57347689) q[0];
sx q[0];
rz(-1.9327123) q[0];
sx q[0];
rz(3.0058577) q[0];
rz(-2.7817338) q[2];
sx q[2];
rz(-2.3198876) q[2];
sx q[2];
rz(0.51279678) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.3914096) q[1];
sx q[1];
rz(-0.58779189) q[1];
sx q[1];
rz(-0.21601917) q[1];
x q[2];
rz(1.7298572) q[3];
sx q[3];
rz(-1.4791227) q[3];
sx q[3];
rz(1.5605469) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.2084674) q[2];
sx q[2];
rz(-1.5310023) q[2];
sx q[2];
rz(2.7739286) q[2];
rz(-2.0465046) q[3];
sx q[3];
rz(-1.0235267) q[3];
sx q[3];
rz(-0.17525214) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4472189) q[0];
sx q[0];
rz(-0.29964888) q[0];
sx q[0];
rz(-1.3391986) q[0];
rz(-0.12204349) q[1];
sx q[1];
rz(-1.3587147) q[1];
sx q[1];
rz(0.36062127) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9741518) q[0];
sx q[0];
rz(-1.1888778) q[0];
sx q[0];
rz(0.72446211) q[0];
rz(-pi) q[1];
rz(-1.4819988) q[2];
sx q[2];
rz(-2.7200903) q[2];
sx q[2];
rz(2.5777566) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.0039978) q[1];
sx q[1];
rz(-1.2892525) q[1];
sx q[1];
rz(2.0291734) q[1];
rz(-0.020645647) q[3];
sx q[3];
rz(-2.8886578) q[3];
sx q[3];
rz(-1.7684574) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.3660761) q[2];
sx q[2];
rz(-1.4630432) q[2];
sx q[2];
rz(2.0029081) q[2];
rz(0.25389296) q[3];
sx q[3];
rz(-0.64520276) q[3];
sx q[3];
rz(2.3937288) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2426386) q[0];
sx q[0];
rz(-2.2518318) q[0];
sx q[0];
rz(-2.1658072) q[0];
rz(1.1294533) q[1];
sx q[1];
rz(-0.44810805) q[1];
sx q[1];
rz(-1.3536369) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24714805) q[0];
sx q[0];
rz(-1.2784875) q[0];
sx q[0];
rz(-2.4185989) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1466402) q[2];
sx q[2];
rz(-1.643702) q[2];
sx q[2];
rz(0.50114252) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.2372386) q[1];
sx q[1];
rz(-1.905958) q[1];
sx q[1];
rz(1.7101014) q[1];
rz(1.9951172) q[3];
sx q[3];
rz(-2.0034308) q[3];
sx q[3];
rz(-0.53782767) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.56753105) q[2];
sx q[2];
rz(-0.83063829) q[2];
sx q[2];
rz(0.91721025) q[2];
rz(-1.5926825) q[3];
sx q[3];
rz(-2.8062688) q[3];
sx q[3];
rz(-2.3622021) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9354644) q[0];
sx q[0];
rz(-1.520282) q[0];
sx q[0];
rz(-0.024854831) q[0];
rz(1.8553597) q[1];
sx q[1];
rz(-1.356946) q[1];
sx q[1];
rz(1.6110427) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2516305) q[0];
sx q[0];
rz(-1.1306835) q[0];
sx q[0];
rz(-1.1391231) q[0];
rz(2.8299061) q[2];
sx q[2];
rz(-0.71060813) q[2];
sx q[2];
rz(-0.59618261) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.74421404) q[1];
sx q[1];
rz(-1.2894003) q[1];
sx q[1];
rz(-1.6069657) q[1];
rz(1.7814662) q[3];
sx q[3];
rz(-1.4349185) q[3];
sx q[3];
rz(-1.5158412) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.020869104) q[2];
sx q[2];
rz(-1.0039696) q[2];
sx q[2];
rz(-1.2845767) q[2];
rz(2.8203216) q[3];
sx q[3];
rz(-2.4925241) q[3];
sx q[3];
rz(0.2215213) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1247509) q[0];
sx q[0];
rz(-1.1469954) q[0];
sx q[0];
rz(2.2220213) q[0];
rz(2.549767) q[1];
sx q[1];
rz(-1.0706736) q[1];
sx q[1];
rz(0.33669534) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.29558173) q[0];
sx q[0];
rz(-0.3008371) q[0];
sx q[0];
rz(1.1016125) q[0];
rz(-pi) q[1];
rz(2.7783103) q[2];
sx q[2];
rz(-1.6840625) q[2];
sx q[2];
rz(-0.47081468) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.7408161) q[1];
sx q[1];
rz(-2.0753521) q[1];
sx q[1];
rz(-1.7187121) q[1];
rz(0.15922861) q[3];
sx q[3];
rz(-1.319229) q[3];
sx q[3];
rz(0.23345527) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.1153974) q[2];
sx q[2];
rz(-2.7112609) q[2];
sx q[2];
rz(-0.78286147) q[2];
rz(3.03249) q[3];
sx q[3];
rz(-1.0166054) q[3];
sx q[3];
rz(1.2469863) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.03610177) q[0];
sx q[0];
rz(-1.6721268) q[0];
sx q[0];
rz(2.2176657) q[0];
rz(-0.52664122) q[1];
sx q[1];
rz(-1.8662607) q[1];
sx q[1];
rz(2.142876) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17820464) q[0];
sx q[0];
rz(-2.9792157) q[0];
sx q[0];
rz(2.8078662) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2800824) q[2];
sx q[2];
rz(-1.6502585) q[2];
sx q[2];
rz(3.0018011) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.0622039) q[1];
sx q[1];
rz(-2.1583589) q[1];
sx q[1];
rz(-0.04157898) q[1];
rz(-pi) q[2];
rz(2.6941006) q[3];
sx q[3];
rz(-2.4119141) q[3];
sx q[3];
rz(0.11451572) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.0052789) q[2];
sx q[2];
rz(-2.2515991) q[2];
sx q[2];
rz(-2.8078553) q[2];
rz(-0.8606832) q[3];
sx q[3];
rz(-1.4849097) q[3];
sx q[3];
rz(3.1349643) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6776047) q[0];
sx q[0];
rz(-1.1900359) q[0];
sx q[0];
rz(0.66587454) q[0];
rz(2.5869276) q[1];
sx q[1];
rz(-1.278109) q[1];
sx q[1];
rz(0.14229933) q[1];
rz(1.8205504) q[2];
sx q[2];
rz(-1.142923) q[2];
sx q[2];
rz(-0.10071071) q[2];
rz(0.71698112) q[3];
sx q[3];
rz(-2.3752799) q[3];
sx q[3];
rz(0.3175288) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
