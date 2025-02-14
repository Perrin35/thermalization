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
rz(-0.91312042) q[0];
sx q[0];
rz(-1.0364391) q[0];
sx q[0];
rz(1.6057462) q[0];
rz(-2.4802471) q[1];
sx q[1];
rz(-1.9547434) q[1];
sx q[1];
rz(-0.43651906) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.8494328) q[0];
sx q[0];
rz(-1.8884553) q[0];
sx q[0];
rz(0.10411291) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4869957) q[2];
sx q[2];
rz(-0.67383781) q[2];
sx q[2];
rz(1.0985653) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.3665665) q[1];
sx q[1];
rz(-2.4373551) q[1];
sx q[1];
rz(1.5077059) q[1];
x q[2];
rz(-2.2546886) q[3];
sx q[3];
rz(-1.4330079) q[3];
sx q[3];
rz(-2.6215214) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.46301699) q[2];
sx q[2];
rz(-2.7028694) q[2];
sx q[2];
rz(-0.92301816) q[2];
rz(-3.0758514) q[3];
sx q[3];
rz(-1.3275423) q[3];
sx q[3];
rz(1.8008697) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1317257) q[0];
sx q[0];
rz(-2.0781524) q[0];
sx q[0];
rz(-0.03751066) q[0];
rz(0.58049479) q[1];
sx q[1];
rz(-0.66941222) q[1];
sx q[1];
rz(-0.69211778) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.53216171) q[0];
sx q[0];
rz(-1.2464644) q[0];
sx q[0];
rz(2.956809) q[0];
rz(-2.0522444) q[2];
sx q[2];
rz(-0.55608597) q[2];
sx q[2];
rz(2.5646445) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.0040652635) q[1];
sx q[1];
rz(-2.1430227) q[1];
sx q[1];
rz(-0.96478421) q[1];
x q[2];
rz(0.75196079) q[3];
sx q[3];
rz(-1.4436566) q[3];
sx q[3];
rz(0.76691915) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.6633501) q[2];
sx q[2];
rz(-1.9448091) q[2];
sx q[2];
rz(-1.7459858) q[2];
rz(-2.9546402) q[3];
sx q[3];
rz(-1.9713277) q[3];
sx q[3];
rz(-0.54005867) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66889399) q[0];
sx q[0];
rz(-2.5085594) q[0];
sx q[0];
rz(-2.582666) q[0];
rz(0.82866296) q[1];
sx q[1];
rz(-1.5907954) q[1];
sx q[1];
rz(2.8899946) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54788816) q[0];
sx q[0];
rz(-1.3798837) q[0];
sx q[0];
rz(2.6949791) q[0];
x q[1];
rz(0.34464406) q[2];
sx q[2];
rz(-2.0953669) q[2];
sx q[2];
rz(1.3087496) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.14519599) q[1];
sx q[1];
rz(-0.72219488) q[1];
sx q[1];
rz(1.1678371) q[1];
rz(-pi) q[2];
rz(0.91355998) q[3];
sx q[3];
rz(-1.0042448) q[3];
sx q[3];
rz(-0.65627894) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.805213) q[2];
sx q[2];
rz(-0.42625913) q[2];
sx q[2];
rz(1.4317929) q[2];
rz(-0.21670565) q[3];
sx q[3];
rz(-1.6101937) q[3];
sx q[3];
rz(-1.7823035) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25877229) q[0];
sx q[0];
rz(-2.6348305) q[0];
sx q[0];
rz(-1.2910507) q[0];
rz(-2.3123815) q[1];
sx q[1];
rz(-2.0334838) q[1];
sx q[1];
rz(2.1270027) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.050796789) q[0];
sx q[0];
rz(-1.8136588) q[0];
sx q[0];
rz(-0.54625752) q[0];
rz(-pi) q[1];
rz(-0.5710602) q[2];
sx q[2];
rz(-0.64415613) q[2];
sx q[2];
rz(1.8824492) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.3922799) q[1];
sx q[1];
rz(-1.5661513) q[1];
sx q[1];
rz(2.1455812) q[1];
rz(-pi) q[2];
rz(1.4271554) q[3];
sx q[3];
rz(-2.3206986) q[3];
sx q[3];
rz(0.088975541) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.1318876) q[2];
sx q[2];
rz(-1.9405126) q[2];
sx q[2];
rz(2.3809872) q[2];
rz(0.46918121) q[3];
sx q[3];
rz(-1.4696591) q[3];
sx q[3];
rz(0.73452264) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0483265) q[0];
sx q[0];
rz(-0.51437298) q[0];
sx q[0];
rz(2.560428) q[0];
rz(1.5900853) q[1];
sx q[1];
rz(-2.4947512) q[1];
sx q[1];
rz(-2.4466628) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7030554) q[0];
sx q[0];
rz(-1.3378654) q[0];
sx q[0];
rz(1.8700048) q[0];
rz(0.30384906) q[2];
sx q[2];
rz(-0.30551791) q[2];
sx q[2];
rz(1.6696827) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.4798812) q[1];
sx q[1];
rz(-1.6146361) q[1];
sx q[1];
rz(2.5978894) q[1];
rz(0.23628037) q[3];
sx q[3];
rz(-1.7485011) q[3];
sx q[3];
rz(1.1760933) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.86357) q[2];
sx q[2];
rz(-1.5333971) q[2];
sx q[2];
rz(-2.8387873) q[2];
rz(1.7152202) q[3];
sx q[3];
rz(-2.7808166) q[3];
sx q[3];
rz(-2.5572131) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.67702174) q[0];
sx q[0];
rz(-0.40957054) q[0];
sx q[0];
rz(-0.69497481) q[0];
rz(2.8547844) q[1];
sx q[1];
rz(-0.60568714) q[1];
sx q[1];
rz(-1.8668176) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5329842) q[0];
sx q[0];
rz(-1.8814981) q[0];
sx q[0];
rz(-1.0644738) q[0];
rz(-pi) q[1];
rz(-2.5996501) q[2];
sx q[2];
rz(-1.1620318) q[2];
sx q[2];
rz(0.33318502) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.50837738) q[1];
sx q[1];
rz(-1.214809) q[1];
sx q[1];
rz(-0.87320019) q[1];
rz(-pi) q[2];
rz(-2.5596095) q[3];
sx q[3];
rz(-1.3701539) q[3];
sx q[3];
rz(1.1021922) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.0605165) q[2];
sx q[2];
rz(-0.21848564) q[2];
sx q[2];
rz(-2.1742353) q[2];
rz(0.71990144) q[3];
sx q[3];
rz(-0.94341174) q[3];
sx q[3];
rz(1.8657743) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28441456) q[0];
sx q[0];
rz(-3.1107749) q[0];
sx q[0];
rz(1.4469294) q[0];
rz(-0.46724304) q[1];
sx q[1];
rz(-1.2930361) q[1];
sx q[1];
rz(1.1955059) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0366096) q[0];
sx q[0];
rz(-0.65147841) q[0];
sx q[0];
rz(2.3250513) q[0];
x q[1];
rz(-1.906997) q[2];
sx q[2];
rz(-0.78374353) q[2];
sx q[2];
rz(-0.23434336) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.0250281) q[1];
sx q[1];
rz(-1.0724853) q[1];
sx q[1];
rz(3.0710941) q[1];
x q[2];
rz(-1.6067664) q[3];
sx q[3];
rz(-1.5580873) q[3];
sx q[3];
rz(0.1606687) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.82470992) q[2];
sx q[2];
rz(-0.83493817) q[2];
sx q[2];
rz(-0.52428025) q[2];
rz(-0.25238642) q[3];
sx q[3];
rz(-0.10656825) q[3];
sx q[3];
rz(3.0805123) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.97279945) q[0];
sx q[0];
rz(-1.1058829) q[0];
sx q[0];
rz(-1.9148781) q[0];
rz(0.75617689) q[1];
sx q[1];
rz(-2.1935479) q[1];
sx q[1];
rz(0.39047584) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.062138012) q[0];
sx q[0];
rz(-1.3083757) q[0];
sx q[0];
rz(-0.45824893) q[0];
x q[1];
rz(-3.0608589) q[2];
sx q[2];
rz(-0.44534031) q[2];
sx q[2];
rz(2.3824218) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.23747102) q[1];
sx q[1];
rz(-1.3573705) q[1];
sx q[1];
rz(-1.2280812) q[1];
rz(-pi) q[2];
x q[2];
rz(0.49440439) q[3];
sx q[3];
rz(-2.1526511) q[3];
sx q[3];
rz(-1.0487919) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.32435027) q[2];
sx q[2];
rz(-0.70049006) q[2];
sx q[2];
rz(1.3520757) q[2];
rz(-2.9663441) q[3];
sx q[3];
rz(-1.8027571) q[3];
sx q[3];
rz(1.9372743) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7290203) q[0];
sx q[0];
rz(-0.46171284) q[0];
sx q[0];
rz(2.4793258) q[0];
rz(2.5114255) q[1];
sx q[1];
rz(-1.2799193) q[1];
sx q[1];
rz(0.23552775) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.61155897) q[0];
sx q[0];
rz(-1.7445824) q[0];
sx q[0];
rz(-0.21428971) q[0];
rz(-pi) q[1];
rz(0.56925242) q[2];
sx q[2];
rz(-1.6496837) q[2];
sx q[2];
rz(1.1752664) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.77068608) q[1];
sx q[1];
rz(-1.5471493) q[1];
sx q[1];
rz(-1.6321983) q[1];
rz(-2.2673549) q[3];
sx q[3];
rz(-1.6739419) q[3];
sx q[3];
rz(-2.646605) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.0887289) q[2];
sx q[2];
rz(-2.1910618) q[2];
sx q[2];
rz(-2.5843184) q[2];
rz(2.3142464) q[3];
sx q[3];
rz(-1.6297623) q[3];
sx q[3];
rz(1.5000337) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2130704) q[0];
sx q[0];
rz(-2.3445573) q[0];
sx q[0];
rz(0.56761566) q[0];
rz(2.7396743) q[1];
sx q[1];
rz(-2.4979976) q[1];
sx q[1];
rz(1.3622805) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3338254) q[0];
sx q[0];
rz(-2.5596715) q[0];
sx q[0];
rz(-2.3533217) q[0];
x q[1];
rz(-2.882233) q[2];
sx q[2];
rz(-2.1742714) q[2];
sx q[2];
rz(-0.063613907) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.6475923) q[1];
sx q[1];
rz(-2.5060182) q[1];
sx q[1];
rz(-1.7032436) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.79093334) q[3];
sx q[3];
rz(-2.2245527) q[3];
sx q[3];
rz(2.4733651) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.42084971) q[2];
sx q[2];
rz(-1.9385312) q[2];
sx q[2];
rz(-2.8813349) q[2];
rz(-2.4506954) q[3];
sx q[3];
rz(-2.2862209) q[3];
sx q[3];
rz(1.6183841) q[3];
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
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8061168) q[0];
sx q[0];
rz(-1.5609043) q[0];
sx q[0];
rz(1.5326395) q[0];
rz(2.7083022) q[1];
sx q[1];
rz(-1.8186124) q[1];
sx q[1];
rz(-1.1846452) q[1];
rz(-2.4240906) q[2];
sx q[2];
rz(-2.2926471) q[2];
sx q[2];
rz(-2.9903464) q[2];
rz(-1.4337486) q[3];
sx q[3];
rz(-2.3564586) q[3];
sx q[3];
rz(-2.3740507) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
