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
rz(0.61252874) q[0];
sx q[0];
rz(-2.5084485) q[0];
sx q[0];
rz(0.36344114) q[0];
rz(2.5545622) q[1];
sx q[1];
rz(-0.5293923) q[1];
sx q[1];
rz(1.6338978) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3359409) q[0];
sx q[0];
rz(-1.664289) q[0];
sx q[0];
rz(0.050430817) q[0];
rz(0.58552758) q[2];
sx q[2];
rz(-2.3351933) q[2];
sx q[2];
rz(-1.6438649) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.9161497) q[1];
sx q[1];
rz(-0.60497153) q[1];
sx q[1];
rz(-1.9523282) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5938183) q[3];
sx q[3];
rz(-2.2507077) q[3];
sx q[3];
rz(0.31991968) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.2734566) q[2];
sx q[2];
rz(-0.09082219) q[2];
sx q[2];
rz(-2.691213) q[2];
rz(3.0716589) q[3];
sx q[3];
rz(-1.115333) q[3];
sx q[3];
rz(2.6993338) q[3];
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
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0792585) q[0];
sx q[0];
rz(-2.313995) q[0];
sx q[0];
rz(-2.8424679) q[0];
rz(-0.66080719) q[1];
sx q[1];
rz(-1.1102763) q[1];
sx q[1];
rz(0.92346907) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.35645928) q[0];
sx q[0];
rz(-1.74455) q[0];
sx q[0];
rz(1.2453402) q[0];
x q[1];
rz(-2.9060676) q[2];
sx q[2];
rz(-1.7851637) q[2];
sx q[2];
rz(1.9247687) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.077775309) q[1];
sx q[1];
rz(-0.37250054) q[1];
sx q[1];
rz(0.17613303) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9808867) q[3];
sx q[3];
rz(-2.514127) q[3];
sx q[3];
rz(-0.80998176) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.6580761) q[2];
sx q[2];
rz(-0.4008171) q[2];
sx q[2];
rz(0.083902396) q[2];
rz(-0.99533254) q[3];
sx q[3];
rz(-1.1570802) q[3];
sx q[3];
rz(1.5108494) q[3];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.29626572) q[0];
sx q[0];
rz(-1.9060598) q[0];
sx q[0];
rz(-0.96257019) q[0];
rz(-2.8341017) q[1];
sx q[1];
rz(-2.2566785) q[1];
sx q[1];
rz(-2.9226551) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6097267) q[0];
sx q[0];
rz(-1.6403188) q[0];
sx q[0];
rz(-2.7951434) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.89246558) q[2];
sx q[2];
rz(-1.7545059) q[2];
sx q[2];
rz(-0.018863686) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.57812095) q[1];
sx q[1];
rz(-1.5679661) q[1];
sx q[1];
rz(1.1712537) q[1];
rz(2.2454065) q[3];
sx q[3];
rz(-2.0749669) q[3];
sx q[3];
rz(-1.653513) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.683005) q[2];
sx q[2];
rz(-1.747749) q[2];
sx q[2];
rz(1.5031507) q[2];
rz(0.25034869) q[3];
sx q[3];
rz(-0.72385794) q[3];
sx q[3];
rz(0.19846465) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4918764) q[0];
sx q[0];
rz(-2.8174077) q[0];
sx q[0];
rz(-2.8247996) q[0];
rz(-3.1253539) q[1];
sx q[1];
rz(-0.79252487) q[1];
sx q[1];
rz(-2.8932755) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.028848948) q[0];
sx q[0];
rz(-0.95579493) q[0];
sx q[0];
rz(-0.93407913) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4961719) q[2];
sx q[2];
rz(-2.0368825) q[2];
sx q[2];
rz(1.8209195) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.8115719) q[1];
sx q[1];
rz(-1.4452126) q[1];
sx q[1];
rz(-2.8126393) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.85217441) q[3];
sx q[3];
rz(-0.2874473) q[3];
sx q[3];
rz(-1.804627) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.877964) q[2];
sx q[2];
rz(-3.0402277) q[2];
sx q[2];
rz(-1.9975837) q[2];
rz(2.3794203) q[3];
sx q[3];
rz(-0.52505571) q[3];
sx q[3];
rz(0.99153668) q[3];
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
sx q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2852729) q[0];
sx q[0];
rz(-0.090949051) q[0];
sx q[0];
rz(2.0891304) q[0];
rz(0.4902803) q[1];
sx q[1];
rz(-1.0888381) q[1];
sx q[1];
rz(3.0163684) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1220684) q[0];
sx q[0];
rz(-0.53036371) q[0];
sx q[0];
rz(0.85774647) q[0];
rz(-pi) q[1];
rz(-0.54791656) q[2];
sx q[2];
rz(-1.8109908) q[2];
sx q[2];
rz(1.554516) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.0908806) q[1];
sx q[1];
rz(-1.6428324) q[1];
sx q[1];
rz(1.3204095) q[1];
x q[2];
rz(-1.8242314) q[3];
sx q[3];
rz(-2.4134688) q[3];
sx q[3];
rz(0.65092939) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.48309717) q[2];
sx q[2];
rz(-0.53762895) q[2];
sx q[2];
rz(-1.5527234) q[2];
rz(3.0850278) q[3];
sx q[3];
rz(-1.6274692) q[3];
sx q[3];
rz(-0.27730832) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5796984) q[0];
sx q[0];
rz(-2.6772006) q[0];
sx q[0];
rz(0.42771801) q[0];
rz(-0.5411717) q[1];
sx q[1];
rz(-0.41053694) q[1];
sx q[1];
rz(-1.3637095) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.83374524) q[0];
sx q[0];
rz(-2.5302093) q[0];
sx q[0];
rz(-1.7627384) q[0];
x q[1];
rz(-3.0156216) q[2];
sx q[2];
rz(-0.4301542) q[2];
sx q[2];
rz(0.90829721) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.1178255) q[1];
sx q[1];
rz(-2.3741873) q[1];
sx q[1];
rz(-1.1249705) q[1];
rz(-pi) q[2];
rz(1.545752) q[3];
sx q[3];
rz(-0.4129172) q[3];
sx q[3];
rz(3.1144547) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-3.0987739) q[2];
sx q[2];
rz(-2.6218178) q[2];
sx q[2];
rz(0.40936142) q[2];
rz(-0.21971075) q[3];
sx q[3];
rz(-1.5922056) q[3];
sx q[3];
rz(-1.6519914) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4396502) q[0];
sx q[0];
rz(-0.42124978) q[0];
sx q[0];
rz(0.14081328) q[0];
rz(-2.6048062) q[1];
sx q[1];
rz(-1.355492) q[1];
sx q[1];
rz(2.7649194) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0842782) q[0];
sx q[0];
rz(-1.8207826) q[0];
sx q[0];
rz(-2.9765073) q[0];
rz(-pi) q[1];
x q[1];
rz(2.656761) q[2];
sx q[2];
rz(-0.73532447) q[2];
sx q[2];
rz(1.5905274) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.72195327) q[1];
sx q[1];
rz(-0.7925539) q[1];
sx q[1];
rz(-2.9591317) q[1];
rz(1.4145109) q[3];
sx q[3];
rz(-1.9816074) q[3];
sx q[3];
rz(2.1235091) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.0385527) q[2];
sx q[2];
rz(-1.1662177) q[2];
sx q[2];
rz(-0.80123365) q[2];
rz(-1.0771105) q[3];
sx q[3];
rz(-2.884765) q[3];
sx q[3];
rz(-0.1424772) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
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
rz(2.899026) q[0];
sx q[0];
rz(-2.9809256) q[0];
sx q[0];
rz(-0.94108474) q[0];
rz(-0.89648992) q[1];
sx q[1];
rz(-1.7231548) q[1];
sx q[1];
rz(-0.81472188) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.488637) q[0];
sx q[0];
rz(-1.205648) q[0];
sx q[0];
rz(-2.7333951) q[0];
rz(-pi) q[1];
rz(-2.9707387) q[2];
sx q[2];
rz(-1.5207982) q[2];
sx q[2];
rz(1.4871666) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.43489698) q[1];
sx q[1];
rz(-2.416947) q[1];
sx q[1];
rz(-2.0963293) q[1];
rz(-pi) q[2];
rz(1.0227459) q[3];
sx q[3];
rz(-2.4999385) q[3];
sx q[3];
rz(-1.9985355) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.5251289) q[2];
sx q[2];
rz(-0.63462555) q[2];
sx q[2];
rz(-2.5647707) q[2];
rz(0.47354627) q[3];
sx q[3];
rz(-1.1487995) q[3];
sx q[3];
rz(0.043341652) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.43805495) q[0];
sx q[0];
rz(-0.75551581) q[0];
sx q[0];
rz(-3.0536861) q[0];
rz(-1.5880623) q[1];
sx q[1];
rz(-0.5843662) q[1];
sx q[1];
rz(2.8839674) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15749685) q[0];
sx q[0];
rz(-1.983674) q[0];
sx q[0];
rz(0.12992125) q[0];
rz(-pi) q[1];
rz(1.6682569) q[2];
sx q[2];
rz(-0.62964338) q[2];
sx q[2];
rz(2.7435945) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.1779319) q[1];
sx q[1];
rz(-2.6448634) q[1];
sx q[1];
rz(-0.50951783) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.99434234) q[3];
sx q[3];
rz(-1.0691158) q[3];
sx q[3];
rz(-0.27015206) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.075500175) q[2];
sx q[2];
rz(-1.4614212) q[2];
sx q[2];
rz(3.063805) q[2];
rz(-1.0545688) q[3];
sx q[3];
rz(-0.74020296) q[3];
sx q[3];
rz(-2.4906702) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8973549) q[0];
sx q[0];
rz(-3.1081508) q[0];
sx q[0];
rz(-0.43575409) q[0];
rz(-3.1089973) q[1];
sx q[1];
rz(-2.4041912) q[1];
sx q[1];
rz(1.3479412) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7028026) q[0];
sx q[0];
rz(-2.3910671) q[0];
sx q[0];
rz(0.010396718) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7392883) q[2];
sx q[2];
rz(-0.56395203) q[2];
sx q[2];
rz(2.5072797) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.0950198) q[1];
sx q[1];
rz(-0.038674861) q[1];
sx q[1];
rz(1.2016618) q[1];
rz(-pi) q[2];
rz(-0.3147821) q[3];
sx q[3];
rz(-1.3482932) q[3];
sx q[3];
rz(-0.088929847) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.2160448) q[2];
sx q[2];
rz(-1.0341158) q[2];
sx q[2];
rz(2.3852868) q[2];
rz(-0.082152724) q[3];
sx q[3];
rz(-0.3923471) q[3];
sx q[3];
rz(-2.4503585) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.77770752) q[0];
sx q[0];
rz(-0.73962279) q[0];
sx q[0];
rz(2.2864322) q[0];
rz(-1.1126407) q[1];
sx q[1];
rz(-2.009195) q[1];
sx q[1];
rz(2.7809273) q[1];
rz(2.1786244) q[2];
sx q[2];
rz(-0.30207015) q[2];
sx q[2];
rz(1.3312726) q[2];
rz(2.2459946) q[3];
sx q[3];
rz(-0.67586312) q[3];
sx q[3];
rz(-2.6296658) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
