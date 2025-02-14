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
rz(-2.8895145) q[0];
sx q[0];
rz(-0.57555389) q[0];
sx q[0];
rz(0.12230305) q[0];
rz(-2.0343434) q[1];
sx q[1];
rz(-2.1802433) q[1];
sx q[1];
rz(0.038318757) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3275546) q[0];
sx q[0];
rz(-0.71056847) q[0];
sx q[0];
rz(2.1585805) q[0];
x q[1];
rz(-1.6898481) q[2];
sx q[2];
rz(-1.6887293) q[2];
sx q[2];
rz(-2.8182507) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.1120834) q[1];
sx q[1];
rz(-1.2473628) q[1];
sx q[1];
rz(-0.98677633) q[1];
rz(-pi) q[2];
rz(-1.5502717) q[3];
sx q[3];
rz(-1.8181268) q[3];
sx q[3];
rz(2.1065245) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.64500874) q[2];
sx q[2];
rz(-1.714548) q[2];
sx q[2];
rz(1.4487779) q[2];
rz(0.19041666) q[3];
sx q[3];
rz(-1.0551635) q[3];
sx q[3];
rz(0.14934389) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9377624) q[0];
sx q[0];
rz(-1.8730524) q[0];
sx q[0];
rz(-2.8835836) q[0];
rz(1.5500655) q[1];
sx q[1];
rz(-1.9015046) q[1];
sx q[1];
rz(0.16664997) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.193522) q[0];
sx q[0];
rz(-1.6662473) q[0];
sx q[0];
rz(0.58346099) q[0];
rz(-pi) q[1];
x q[1];
rz(0.62549297) q[2];
sx q[2];
rz(-2.6776367) q[2];
sx q[2];
rz(-2.6044012) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.0936792) q[1];
sx q[1];
rz(-2.4166345) q[1];
sx q[1];
rz(-2.9525501) q[1];
x q[2];
rz(2.3349668) q[3];
sx q[3];
rz(-2.4203033) q[3];
sx q[3];
rz(1.2794577) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.0698645) q[2];
sx q[2];
rz(-1.122415) q[2];
sx q[2];
rz(1.2321164) q[2];
rz(-2.2169436) q[3];
sx q[3];
rz(-2.2053714) q[3];
sx q[3];
rz(2.0360951) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.6919747) q[0];
sx q[0];
rz(-1.6716577) q[0];
sx q[0];
rz(-3.1396507) q[0];
rz(3.107403) q[1];
sx q[1];
rz(-1.2119774) q[1];
sx q[1];
rz(1.5984104) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40638518) q[0];
sx q[0];
rz(-1.8282561) q[0];
sx q[0];
rz(-0.46528146) q[0];
rz(1.1748741) q[2];
sx q[2];
rz(-0.87475077) q[2];
sx q[2];
rz(-0.093122236) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.48422563) q[1];
sx q[1];
rz(-1.5316841) q[1];
sx q[1];
rz(0.35376057) q[1];
rz(-pi) q[2];
rz(2.8974124) q[3];
sx q[3];
rz(-0.93726087) q[3];
sx q[3];
rz(-1.7297431) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.2446642) q[2];
sx q[2];
rz(-2.7074773) q[2];
sx q[2];
rz(-2.1042692) q[2];
rz(2.0364929) q[3];
sx q[3];
rz(-1.5367855) q[3];
sx q[3];
rz(1.1387811) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26745519) q[0];
sx q[0];
rz(-2.656811) q[0];
sx q[0];
rz(1.4111891) q[0];
rz(2.5022068) q[1];
sx q[1];
rz(-1.6849898) q[1];
sx q[1];
rz(-0.04714084) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.61091238) q[0];
sx q[0];
rz(-2.147104) q[0];
sx q[0];
rz(1.9082101) q[0];
x q[1];
rz(-2.669692) q[2];
sx q[2];
rz(-1.7444897) q[2];
sx q[2];
rz(-1.7431517) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.47430965) q[1];
sx q[1];
rz(-1.9243) q[1];
sx q[1];
rz(-2.1992285) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2940801) q[3];
sx q[3];
rz(-1.3467595) q[3];
sx q[3];
rz(-0.23875313) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.3186657) q[2];
sx q[2];
rz(-1.1456127) q[2];
sx q[2];
rz(1.9226496) q[2];
rz(-2.8259891) q[3];
sx q[3];
rz(-3.0867519) q[3];
sx q[3];
rz(0.1489197) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8197935) q[0];
sx q[0];
rz(-0.55235523) q[0];
sx q[0];
rz(2.7929982) q[0];
rz(1.2527342) q[1];
sx q[1];
rz(-1.9786973) q[1];
sx q[1];
rz(1.8399651) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0785261) q[0];
sx q[0];
rz(-1.491703) q[0];
sx q[0];
rz(2.3344759) q[0];
x q[1];
rz(1.1548066) q[2];
sx q[2];
rz(-2.092053) q[2];
sx q[2];
rz(-0.74300569) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.0685454) q[1];
sx q[1];
rz(-1.3497735) q[1];
sx q[1];
rz(-0.18494341) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0318448) q[3];
sx q[3];
rz(-2.5205118) q[3];
sx q[3];
rz(-0.10824848) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.9584413) q[2];
sx q[2];
rz(-2.8330467) q[2];
sx q[2];
rz(1.6661673) q[2];
rz(-2.4557377) q[3];
sx q[3];
rz(-2.1619022) q[3];
sx q[3];
rz(-0.53387749) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1118065) q[0];
sx q[0];
rz(-0.15203467) q[0];
sx q[0];
rz(1.6492122) q[0];
rz(-2.1624508) q[1];
sx q[1];
rz(-1.5359595) q[1];
sx q[1];
rz(-1.917256) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2164152) q[0];
sx q[0];
rz(-1.5402216) q[0];
sx q[0];
rz(0.19801099) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5923205) q[2];
sx q[2];
rz(-1.7584929) q[2];
sx q[2];
rz(2.930738) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.27367941) q[1];
sx q[1];
rz(-0.95703322) q[1];
sx q[1];
rz(3.1263208) q[1];
x q[2];
rz(-0.94131366) q[3];
sx q[3];
rz(-1.286003) q[3];
sx q[3];
rz(-2.3465057) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.41902038) q[2];
sx q[2];
rz(-1.4973065) q[2];
sx q[2];
rz(0.9160308) q[2];
rz(-1.7570868) q[3];
sx q[3];
rz(-1.8839096) q[3];
sx q[3];
rz(-2.4345583) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
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
rz(-0.0093805669) q[0];
sx q[0];
rz(-3.0896602) q[0];
sx q[0];
rz(2.1873271) q[0];
rz(2.7047899) q[1];
sx q[1];
rz(-1.5635468) q[1];
sx q[1];
rz(0.11016914) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.018579114) q[0];
sx q[0];
rz(-1.8948104) q[0];
sx q[0];
rz(1.3186245) q[0];
x q[1];
rz(-1.0344347) q[2];
sx q[2];
rz(-2.1948994) q[2];
sx q[2];
rz(-2.1339061) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.4272144) q[1];
sx q[1];
rz(-2.7958779) q[1];
sx q[1];
rz(-1.5835254) q[1];
rz(-2.8875868) q[3];
sx q[3];
rz(-2.0585103) q[3];
sx q[3];
rz(2.8539477) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.43180141) q[2];
sx q[2];
rz(-3.1352037) q[2];
sx q[2];
rz(2.1978281) q[2];
rz(-0.56378311) q[3];
sx q[3];
rz(-2.0457334) q[3];
sx q[3];
rz(-2.0311267) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3540045) q[0];
sx q[0];
rz(-0.28689757) q[0];
sx q[0];
rz(-0.34550825) q[0];
rz(-1.0954789) q[1];
sx q[1];
rz(-1.4778719) q[1];
sx q[1];
rz(1.5798205) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6248572) q[0];
sx q[0];
rz(-0.78291946) q[0];
sx q[0];
rz(1.6153512) q[0];
rz(-pi) q[1];
rz(-2.2279635) q[2];
sx q[2];
rz(-1.6198748) q[2];
sx q[2];
rz(2.9947481) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.53198481) q[1];
sx q[1];
rz(-2.4188305) q[1];
sx q[1];
rz(-1.0015798) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0828219) q[3];
sx q[3];
rz(-0.49859014) q[3];
sx q[3];
rz(0.22496794) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.76274189) q[2];
sx q[2];
rz(-0.26143917) q[2];
sx q[2];
rz(-1.4307107) q[2];
rz(0.53449574) q[3];
sx q[3];
rz(-1.675324) q[3];
sx q[3];
rz(-2.1889595) q[3];
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
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5296103) q[0];
sx q[0];
rz(-3.072325) q[0];
sx q[0];
rz(-1.8310504) q[0];
rz(2.2546841) q[1];
sx q[1];
rz(-2.3014018) q[1];
sx q[1];
rz(-1.2695674) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7879155) q[0];
sx q[0];
rz(-1.4424818) q[0];
sx q[0];
rz(-1.2169891) q[0];
rz(-2.9863556) q[2];
sx q[2];
rz(-1.198215) q[2];
sx q[2];
rz(2.9205204) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.7901526) q[1];
sx q[1];
rz(-1.8428486) q[1];
sx q[1];
rz(0.2403918) q[1];
rz(-0.19143243) q[3];
sx q[3];
rz(-0.20729724) q[3];
sx q[3];
rz(1.8199004) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.612192) q[2];
sx q[2];
rz(-1.8058913) q[2];
sx q[2];
rz(-2.9086435) q[2];
rz(-1.285078) q[3];
sx q[3];
rz(-2.464747) q[3];
sx q[3];
rz(-2.7555833) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9832298) q[0];
sx q[0];
rz(-0.60929275) q[0];
sx q[0];
rz(-0.097271517) q[0];
rz(-2.3376047) q[1];
sx q[1];
rz(-0.709788) q[1];
sx q[1];
rz(0.21673094) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0618503) q[0];
sx q[0];
rz(-1.6260132) q[0];
sx q[0];
rz(1.7415857) q[0];
rz(-pi) q[1];
x q[1];
rz(0.31483908) q[2];
sx q[2];
rz(-0.14893571) q[2];
sx q[2];
rz(0.51477942) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.3236802) q[1];
sx q[1];
rz(-2.2127832) q[1];
sx q[1];
rz(-2.2873408) q[1];
x q[2];
rz(2.4794934) q[3];
sx q[3];
rz(-1.2939014) q[3];
sx q[3];
rz(-0.84018842) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.29297605) q[2];
sx q[2];
rz(-2.8605707) q[2];
sx q[2];
rz(0.49523735) q[2];
rz(-1.5589335) q[3];
sx q[3];
rz(-2.0571183) q[3];
sx q[3];
rz(-2.9787279) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
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
rz(-2.0804629) q[0];
sx q[0];
rz(-1.6292138) q[0];
sx q[0];
rz(-0.52393352) q[0];
rz(-1.0864661) q[1];
sx q[1];
rz(-2.6230984) q[1];
sx q[1];
rz(-2.7883504) q[1];
rz(-1.3475781) q[2];
sx q[2];
rz(-0.032608727) q[2];
sx q[2];
rz(-0.41667117) q[2];
rz(-2.2829655) q[3];
sx q[3];
rz(-1.6391907) q[3];
sx q[3];
rz(-1.0427604) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
