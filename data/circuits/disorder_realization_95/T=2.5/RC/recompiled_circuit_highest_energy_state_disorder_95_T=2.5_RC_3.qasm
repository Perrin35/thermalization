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
rz(-2.4565826) q[0];
sx q[0];
rz(-0.69184875) q[0];
sx q[0];
rz(-0.43842167) q[0];
rz(3.0216079) q[1];
sx q[1];
rz(3.383145) q[1];
sx q[1];
rz(8.4278843) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0783421) q[0];
sx q[0];
rz(-1.7228427) q[0];
sx q[0];
rz(0.087910533) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6108151) q[2];
sx q[2];
rz(-1.3005031) q[2];
sx q[2];
rz(1.1129462) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.01693657) q[1];
sx q[1];
rz(-1.1775832) q[1];
sx q[1];
rz(-0.26166088) q[1];
rz(-1.3658872) q[3];
sx q[3];
rz(-1.3384322) q[3];
sx q[3];
rz(3.0099677) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.74242705) q[2];
sx q[2];
rz(-2.4319067) q[2];
sx q[2];
rz(2.4281375) q[2];
rz(-0.82792884) q[3];
sx q[3];
rz(-1.4917587) q[3];
sx q[3];
rz(2.2152065) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4543318) q[0];
sx q[0];
rz(-0.61606544) q[0];
sx q[0];
rz(1.2607505) q[0];
rz(-2.8997391) q[1];
sx q[1];
rz(-1.8459903) q[1];
sx q[1];
rz(-0.58580011) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0457092) q[0];
sx q[0];
rz(-1.2564141) q[0];
sx q[0];
rz(0.16736253) q[0];
rz(-1.5795779) q[2];
sx q[2];
rz(-1.9987371) q[2];
sx q[2];
rz(-1.2195171) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-3.0351981) q[1];
sx q[1];
rz(-1.2751586) q[1];
sx q[1];
rz(0.8080478) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5590267) q[3];
sx q[3];
rz(-1.5705207) q[3];
sx q[3];
rz(-2.2451412) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.62023097) q[2];
sx q[2];
rz(-0.62180454) q[2];
sx q[2];
rz(-0.28548959) q[2];
rz(-0.96389687) q[3];
sx q[3];
rz(-2.4745092) q[3];
sx q[3];
rz(2.9662761) q[3];
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
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6730839) q[0];
sx q[0];
rz(-2.5991169) q[0];
sx q[0];
rz(-2.8520404) q[0];
rz(-0.30329224) q[1];
sx q[1];
rz(-1.2723609) q[1];
sx q[1];
rz(1.2264651) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68115409) q[0];
sx q[0];
rz(-1.9931721) q[0];
sx q[0];
rz(2.284221) q[0];
x q[1];
rz(-0.3704088) q[2];
sx q[2];
rz(-0.61282149) q[2];
sx q[2];
rz(-2.1343663) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.94680601) q[1];
sx q[1];
rz(-2.1420519) q[1];
sx q[1];
rz(-2.8033546) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.11762186) q[3];
sx q[3];
rz(-0.78946992) q[3];
sx q[3];
rz(2.5817613) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.6931285) q[2];
sx q[2];
rz(-0.99702865) q[2];
sx q[2];
rz(-2.2779951) q[2];
rz(3.038285) q[3];
sx q[3];
rz(-1.7938675) q[3];
sx q[3];
rz(-0.13532713) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1358262) q[0];
sx q[0];
rz(-2.5209881) q[0];
sx q[0];
rz(0.045106877) q[0];
rz(-0.82334423) q[1];
sx q[1];
rz(-0.38213676) q[1];
sx q[1];
rz(-0.31450048) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.088011) q[0];
sx q[0];
rz(-1.6550261) q[0];
sx q[0];
rz(0.65902577) q[0];
rz(2.7240924) q[2];
sx q[2];
rz(-1.0733777) q[2];
sx q[2];
rz(-0.50499798) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-3.1005974) q[1];
sx q[1];
rz(-1.081859) q[1];
sx q[1];
rz(-1.3616427) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.32415402) q[3];
sx q[3];
rz(-1.3676893) q[3];
sx q[3];
rz(-2.5517011) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.0663674) q[2];
sx q[2];
rz(-1.4868569) q[2];
sx q[2];
rz(0.60869795) q[2];
rz(-1.4501976) q[3];
sx q[3];
rz(-0.83289731) q[3];
sx q[3];
rz(-2.6628185) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4920014) q[0];
sx q[0];
rz(-2.2719125) q[0];
sx q[0];
rz(-0.087652303) q[0];
rz(1.8805257) q[1];
sx q[1];
rz(-1.7365716) q[1];
sx q[1];
rz(0.87439775) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0387303) q[0];
sx q[0];
rz(-2.4695354) q[0];
sx q[0];
rz(1.5182537) q[0];
rz(-pi) q[1];
rz(0.034406729) q[2];
sx q[2];
rz(-2.161247) q[2];
sx q[2];
rz(1.4077983) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.72681422) q[1];
sx q[1];
rz(-2.2764479) q[1];
sx q[1];
rz(0.7945286) q[1];
rz(-2.9196872) q[3];
sx q[3];
rz(-2.1067224) q[3];
sx q[3];
rz(1.3886896) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.7586691) q[2];
sx q[2];
rz(-2.1047968) q[2];
sx q[2];
rz(-2.5385638) q[2];
rz(2.1488819) q[3];
sx q[3];
rz(-2.755136) q[3];
sx q[3];
rz(2.4934798) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8822766) q[0];
sx q[0];
rz(-2.3969605) q[0];
sx q[0];
rz(-0.69389206) q[0];
rz(-2.4248185) q[1];
sx q[1];
rz(-2.0718772) q[1];
sx q[1];
rz(-0.96376354) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1070654) q[0];
sx q[0];
rz(-1.1363159) q[0];
sx q[0];
rz(-1.700842) q[0];
rz(-pi) q[1];
rz(-2.9593857) q[2];
sx q[2];
rz(-1.4002396) q[2];
sx q[2];
rz(2.272416) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.9020667) q[1];
sx q[1];
rz(-2.0888302) q[1];
sx q[1];
rz(-2.7076028) q[1];
rz(1.5429872) q[3];
sx q[3];
rz(-2.1429792) q[3];
sx q[3];
rz(-1.5554626) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.11513772) q[2];
sx q[2];
rz(-0.28885767) q[2];
sx q[2];
rz(-3.0333983) q[2];
rz(2.1859956) q[3];
sx q[3];
rz(-0.038766131) q[3];
sx q[3];
rz(2.5819216) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.74681246) q[0];
sx q[0];
rz(-2.96947) q[0];
sx q[0];
rz(-3.0346003) q[0];
rz(-2.6394898) q[1];
sx q[1];
rz(-1.6306449) q[1];
sx q[1];
rz(-2.9923901) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0923) q[0];
sx q[0];
rz(-0.23173103) q[0];
sx q[0];
rz(2.2845303) q[0];
rz(-2.4154941) q[2];
sx q[2];
rz(-0.70961414) q[2];
sx q[2];
rz(-2.6626056) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.97905695) q[1];
sx q[1];
rz(-1.874089) q[1];
sx q[1];
rz(-0.83502646) q[1];
rz(-pi) q[2];
rz(-0.42909166) q[3];
sx q[3];
rz(-1.7011931) q[3];
sx q[3];
rz(-0.71500378) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.2627829) q[2];
sx q[2];
rz(-2.171319) q[2];
sx q[2];
rz(-2.523017) q[2];
rz(0.68459073) q[3];
sx q[3];
rz(-0.22817831) q[3];
sx q[3];
rz(2.1555307) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8743643) q[0];
sx q[0];
rz(-1.769913) q[0];
sx q[0];
rz(-0.57269639) q[0];
rz(1.0193846) q[1];
sx q[1];
rz(-0.23713325) q[1];
sx q[1];
rz(3.0651029) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9534308) q[0];
sx q[0];
rz(-2.5929586) q[0];
sx q[0];
rz(-1.3004957) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8440153) q[2];
sx q[2];
rz(-1.3598765) q[2];
sx q[2];
rz(0.95266137) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.9323378) q[1];
sx q[1];
rz(-0.9088974) q[1];
sx q[1];
rz(-0.13929892) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0763909) q[3];
sx q[3];
rz(-2.3870584) q[3];
sx q[3];
rz(0.86598321) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.0387592) q[2];
sx q[2];
rz(-2.2566654) q[2];
sx q[2];
rz(-0.49212512) q[2];
rz(2.7627908) q[3];
sx q[3];
rz(-0.4327966) q[3];
sx q[3];
rz(-2.2545599) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.73961306) q[0];
sx q[0];
rz(-0.48960632) q[0];
sx q[0];
rz(0.42994764) q[0];
rz(2.6509189) q[1];
sx q[1];
rz(-2.6538167) q[1];
sx q[1];
rz(-2.1388163) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.96129629) q[0];
sx q[0];
rz(-1.8334098) q[0];
sx q[0];
rz(1.0680593) q[0];
x q[1];
rz(1.6677708) q[2];
sx q[2];
rz(-0.67296689) q[2];
sx q[2];
rz(0.99373465) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.44081894) q[1];
sx q[1];
rz(-2.6731468) q[1];
sx q[1];
rz(-0.10611315) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4409695) q[3];
sx q[3];
rz(-1.2081606) q[3];
sx q[3];
rz(0.36421916) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.7542725) q[2];
sx q[2];
rz(-2.5405799) q[2];
sx q[2];
rz(2.4499272) q[2];
rz(0.12868853) q[3];
sx q[3];
rz(-1.5494989) q[3];
sx q[3];
rz(0.58840978) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.16828951) q[0];
sx q[0];
rz(-3.0423218) q[0];
sx q[0];
rz(-2.5183992) q[0];
rz(-0.19459952) q[1];
sx q[1];
rz(-2.01229) q[1];
sx q[1];
rz(0.87619877) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0191747) q[0];
sx q[0];
rz(-0.12309531) q[0];
sx q[0];
rz(-1.9341296) q[0];
x q[1];
rz(-0.74157115) q[2];
sx q[2];
rz(-1.4585765) q[2];
sx q[2];
rz(0.73597344) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.9651803) q[1];
sx q[1];
rz(-1.1438055) q[1];
sx q[1];
rz(-1.7199442) q[1];
rz(-1.9153173) q[3];
sx q[3];
rz(-2.1117941) q[3];
sx q[3];
rz(1.8581543) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.1759922) q[2];
sx q[2];
rz(-2.8922562) q[2];
sx q[2];
rz(2.5052137) q[2];
rz(-1.5198358) q[3];
sx q[3];
rz(-0.88080019) q[3];
sx q[3];
rz(2.9145068) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(2.8272921) q[0];
sx q[0];
rz(-1.5416332) q[0];
sx q[0];
rz(2.736349) q[0];
rz(-1.4018519) q[1];
sx q[1];
rz(-1.6442465) q[1];
sx q[1];
rz(0.13784611) q[1];
rz(-2.5194216) q[2];
sx q[2];
rz(-1.6710186) q[2];
sx q[2];
rz(0.98115151) q[2];
rz(-0.53608175) q[3];
sx q[3];
rz(-1.5403219) q[3];
sx q[3];
rz(2.2684682) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
