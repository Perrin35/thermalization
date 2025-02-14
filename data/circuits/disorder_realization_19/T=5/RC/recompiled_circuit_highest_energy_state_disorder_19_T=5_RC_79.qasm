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
rz(-2.6416566) q[0];
sx q[0];
rz(-1.9498107) q[0];
sx q[0];
rz(-1.7697822) q[0];
rz(1.6504047) q[1];
sx q[1];
rz(-0.86725441) q[1];
sx q[1];
rz(0.39630085) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4258408) q[0];
sx q[0];
rz(-1.6501199) q[0];
sx q[0];
rz(-0.21004814) q[0];
rz(-pi) q[1];
rz(1.9849586) q[2];
sx q[2];
rz(-2.7348619) q[2];
sx q[2];
rz(-1.0309645) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.635879) q[1];
sx q[1];
rz(-1.6649655) q[1];
sx q[1];
rz(-2.0999184) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2788775) q[3];
sx q[3];
rz(-1.1905498) q[3];
sx q[3];
rz(1.0060665) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.0013915) q[2];
sx q[2];
rz(-0.80696693) q[2];
sx q[2];
rz(1.8001455) q[2];
rz(-1.0691102) q[3];
sx q[3];
rz(-0.1461229) q[3];
sx q[3];
rz(1.7599546) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
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
rz(0.57217252) q[0];
sx q[0];
rz(-2.2241346) q[0];
sx q[0];
rz(0.65895748) q[0];
rz(3.068889) q[1];
sx q[1];
rz(-2.0854988) q[1];
sx q[1];
rz(1.8812077) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.35048553) q[0];
sx q[0];
rz(-2.022615) q[0];
sx q[0];
rz(3.0511887) q[0];
rz(-2.6953893) q[2];
sx q[2];
rz(-1.1006736) q[2];
sx q[2];
rz(2.0044495) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.9483794) q[1];
sx q[1];
rz(-0.27336379) q[1];
sx q[1];
rz(-0.24250077) q[1];
rz(-pi) q[2];
rz(-0.96003344) q[3];
sx q[3];
rz(-1.472578) q[3];
sx q[3];
rz(-1.7383505) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.2468804) q[2];
sx q[2];
rz(-2.3983045) q[2];
sx q[2];
rz(1.7903719) q[2];
rz(-2.5249935) q[3];
sx q[3];
rz(-1.0298046) q[3];
sx q[3];
rz(-2.8504347) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-0.63224822) q[0];
sx q[0];
rz(-0.46724874) q[0];
sx q[0];
rz(-0.21009357) q[0];
rz(-0.083077438) q[1];
sx q[1];
rz(-1.9030842) q[1];
sx q[1];
rz(3.074379) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7930896) q[0];
sx q[0];
rz(-1.5843035) q[0];
sx q[0];
rz(-1.5828787) q[0];
rz(-pi) q[1];
rz(-2.1234197) q[2];
sx q[2];
rz(-2.6913096) q[2];
sx q[2];
rz(1.0339111) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.9720502) q[1];
sx q[1];
rz(-1.2360555) q[1];
sx q[1];
rz(2.2856345) q[1];
rz(0.10566575) q[3];
sx q[3];
rz(-0.70836954) q[3];
sx q[3];
rz(1.8745223) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.3136966) q[2];
sx q[2];
rz(-1.6500429) q[2];
sx q[2];
rz(-1.776604) q[2];
rz(2.7374173) q[3];
sx q[3];
rz(-1.3173236) q[3];
sx q[3];
rz(-1.1990168) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86554027) q[0];
sx q[0];
rz(-1.7409538) q[0];
sx q[0];
rz(-0.47819594) q[0];
rz(-0.95282355) q[1];
sx q[1];
rz(-0.15004221) q[1];
sx q[1];
rz(-2.731954) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9347246) q[0];
sx q[0];
rz(-1.5971916) q[0];
sx q[0];
rz(2.4858263) q[0];
x q[1];
rz(-2.6851321) q[2];
sx q[2];
rz(-1.6586411) q[2];
sx q[2];
rz(2.8473471) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.8781902) q[1];
sx q[1];
rz(-0.95118427) q[1];
sx q[1];
rz(2.3248169) q[1];
rz(-pi) q[2];
rz(1.9799819) q[3];
sx q[3];
rz(-2.7936862) q[3];
sx q[3];
rz(-2.6285841) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.59565583) q[2];
sx q[2];
rz(-1.4384392) q[2];
sx q[2];
rz(-2.0036009) q[2];
rz(-0.99500895) q[3];
sx q[3];
rz(-2.5840839) q[3];
sx q[3];
rz(-0.086624302) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.397641) q[0];
sx q[0];
rz(-1.3084363) q[0];
sx q[0];
rz(0.27807903) q[0];
rz(-1.8771578) q[1];
sx q[1];
rz(-1.2828628) q[1];
sx q[1];
rz(1.5778731) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5681158) q[0];
sx q[0];
rz(-1.2088803) q[0];
sx q[0];
rz(-0.135735) q[0];
rz(-0.78863849) q[2];
sx q[2];
rz(-1.83162) q[2];
sx q[2];
rz(-2.3344085) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.75018308) q[1];
sx q[1];
rz(-0.58779189) q[1];
sx q[1];
rz(2.9255735) q[1];
rz(-1.4117354) q[3];
sx q[3];
rz(-1.4791227) q[3];
sx q[3];
rz(1.5605469) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.9331253) q[2];
sx q[2];
rz(-1.6105904) q[2];
sx q[2];
rz(0.36766407) q[2];
rz(1.095088) q[3];
sx q[3];
rz(-2.118066) q[3];
sx q[3];
rz(-2.9663405) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6943738) q[0];
sx q[0];
rz(-0.29964888) q[0];
sx q[0];
rz(1.8023941) q[0];
rz(0.12204349) q[1];
sx q[1];
rz(-1.782878) q[1];
sx q[1];
rz(-2.7809714) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.339645) q[0];
sx q[0];
rz(-2.3390798) q[0];
sx q[0];
rz(2.5967477) q[0];
rz(-pi) q[1];
rz(0.039741411) q[2];
sx q[2];
rz(-1.1510599) q[2];
sx q[2];
rz(2.6750203) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.569446) q[1];
sx q[1];
rz(-1.1317557) q[1];
sx q[1];
rz(-0.31198685) q[1];
rz(-pi) q[2];
rz(0.020645647) q[3];
sx q[3];
rz(-0.25293487) q[3];
sx q[3];
rz(1.3731352) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.7755166) q[2];
sx q[2];
rz(-1.6785494) q[2];
sx q[2];
rz(-2.0029081) q[2];
rz(-2.8876997) q[3];
sx q[3];
rz(-2.4963899) q[3];
sx q[3];
rz(0.74786389) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.89895407) q[0];
sx q[0];
rz(-0.88976088) q[0];
sx q[0];
rz(-0.97578543) q[0];
rz(-1.1294533) q[1];
sx q[1];
rz(-0.44810805) q[1];
sx q[1];
rz(1.3536369) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0669374) q[0];
sx q[0];
rz(-2.2569509) q[0];
sx q[0];
rz(-1.1891436) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9949525) q[2];
sx q[2];
rz(-1.4978906) q[2];
sx q[2];
rz(2.6404501) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.2372386) q[1];
sx q[1];
rz(-1.905958) q[1];
sx q[1];
rz(-1.4314913) q[1];
x q[2];
rz(0.72809345) q[3];
sx q[3];
rz(-2.545176) q[3];
sx q[3];
rz(2.8567258) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.5740616) q[2];
sx q[2];
rz(-2.3109544) q[2];
sx q[2];
rz(0.91721025) q[2];
rz(-1.5489102) q[3];
sx q[3];
rz(-2.8062688) q[3];
sx q[3];
rz(-0.77939051) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9354644) q[0];
sx q[0];
rz(-1.6213106) q[0];
sx q[0];
rz(3.1167378) q[0];
rz(1.8553597) q[1];
sx q[1];
rz(-1.356946) q[1];
sx q[1];
rz(-1.53055) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0760374) q[0];
sx q[0];
rz(-2.5352056) q[0];
sx q[0];
rz(-2.4151912) q[0];
rz(-pi) q[1];
rz(-0.68629219) q[2];
sx q[2];
rz(-1.3694021) q[2];
sx q[2];
rz(-1.9274595) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.3973786) q[1];
sx q[1];
rz(-1.8521924) q[1];
sx q[1];
rz(-1.534627) q[1];
x q[2];
rz(-1.7814662) q[3];
sx q[3];
rz(-1.4349185) q[3];
sx q[3];
rz(-1.6257515) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.020869104) q[2];
sx q[2];
rz(-2.1376231) q[2];
sx q[2];
rz(1.2845767) q[2];
rz(0.32127109) q[3];
sx q[3];
rz(-0.64906859) q[3];
sx q[3];
rz(-2.9200714) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1247509) q[0];
sx q[0];
rz(-1.9945972) q[0];
sx q[0];
rz(0.91957134) q[0];
rz(-2.549767) q[1];
sx q[1];
rz(-1.0706736) q[1];
sx q[1];
rz(2.8048973) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.82430912) q[0];
sx q[0];
rz(-1.4364087) q[0];
sx q[0];
rz(-1.8407673) q[0];
rz(1.4496964) q[2];
sx q[2];
rz(-1.2099488) q[2];
sx q[2];
rz(1.0570414) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.40077651) q[1];
sx q[1];
rz(-2.0753521) q[1];
sx q[1];
rz(1.7187121) q[1];
rz(2.982364) q[3];
sx q[3];
rz(-1.319229) q[3];
sx q[3];
rz(2.9081374) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.0261953) q[2];
sx q[2];
rz(-2.7112609) q[2];
sx q[2];
rz(0.78286147) q[2];
rz(-3.03249) q[3];
sx q[3];
rz(-2.1249873) q[3];
sx q[3];
rz(-1.8946064) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1054909) q[0];
sx q[0];
rz(-1.4694659) q[0];
sx q[0];
rz(-0.92392695) q[0];
rz(0.52664122) q[1];
sx q[1];
rz(-1.8662607) q[1];
sx q[1];
rz(-2.142876) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.625542) q[0];
sx q[0];
rz(-1.417451) q[0];
sx q[0];
rz(1.6244066) q[0];
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
rz(-2.0793887) q[1];
sx q[1];
rz(-2.1583589) q[1];
sx q[1];
rz(-3.1000137) q[1];
rz(-1.9400334) q[3];
sx q[3];
rz(-2.2155361) q[3];
sx q[3];
rz(0.45756868) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.1363137) q[2];
sx q[2];
rz(-0.88999358) q[2];
sx q[2];
rz(0.33373731) q[2];
rz(-2.2809095) q[3];
sx q[3];
rz(-1.656683) q[3];
sx q[3];
rz(3.1349643) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6776047) q[0];
sx q[0];
rz(-1.9515568) q[0];
sx q[0];
rz(-2.4757181) q[0];
rz(-0.55466501) q[1];
sx q[1];
rz(-1.278109) q[1];
sx q[1];
rz(0.14229933) q[1];
rz(-1.8205504) q[2];
sx q[2];
rz(-1.9986696) q[2];
sx q[2];
rz(3.0408819) q[2];
rz(-2.1347682) q[3];
sx q[3];
rz(-2.1208616) q[3];
sx q[3];
rz(-1.943945) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
