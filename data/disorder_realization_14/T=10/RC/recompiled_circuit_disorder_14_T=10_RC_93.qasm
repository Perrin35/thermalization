OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.68552652) q[0];
sx q[0];
rz(-2.752562) q[0];
sx q[0];
rz(-2.2580137) q[0];
rz(3.1318624) q[1];
sx q[1];
rz(-1.6844123) q[1];
sx q[1];
rz(-1.943346) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1238026) q[0];
sx q[0];
rz(-1.4746116) q[0];
sx q[0];
rz(-0.17734806) q[0];
rz(-2.8134402) q[2];
sx q[2];
rz(-0.80291623) q[2];
sx q[2];
rz(2.3053055) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.27904305) q[1];
sx q[1];
rz(-2.0837796) q[1];
sx q[1];
rz(-1.1898477) q[1];
x q[2];
rz(1.8331068) q[3];
sx q[3];
rz(-0.17671083) q[3];
sx q[3];
rz(-1.3486598) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.9464232) q[2];
sx q[2];
rz(-0.98313466) q[2];
sx q[2];
rz(2.9602489) q[2];
rz(-2.8803853) q[3];
sx q[3];
rz(-1.2657335) q[3];
sx q[3];
rz(-2.3852824) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3023892) q[0];
sx q[0];
rz(-0.2897245) q[0];
sx q[0];
rz(2.7547577) q[0];
rz(0.50239262) q[1];
sx q[1];
rz(-2.1680809) q[1];
sx q[1];
rz(-1.5997255) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5889266) q[0];
sx q[0];
rz(-0.98836556) q[0];
sx q[0];
rz(-1.0731339) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2453169) q[2];
sx q[2];
rz(-0.4782593) q[2];
sx q[2];
rz(1.9667728) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.6501573) q[1];
sx q[1];
rz(-3.1045034) q[1];
sx q[1];
rz(2.8703719) q[1];
rz(-2.6582791) q[3];
sx q[3];
rz(-1.9938853) q[3];
sx q[3];
rz(1.2276633) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.60275045) q[2];
sx q[2];
rz(-2.2637612) q[2];
sx q[2];
rz(-1.7269469) q[2];
rz(0.85033068) q[3];
sx q[3];
rz(-0.43262216) q[3];
sx q[3];
rz(-1.458228) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26281115) q[0];
sx q[0];
rz(-1.6101863) q[0];
sx q[0];
rz(-2.5313654) q[0];
rz(-1.8521076) q[1];
sx q[1];
rz(-0.97924966) q[1];
sx q[1];
rz(-0.99197018) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.85342825) q[0];
sx q[0];
rz(-2.5531904) q[0];
sx q[0];
rz(1.9073652) q[0];
rz(-pi) q[1];
rz(-2.8875071) q[2];
sx q[2];
rz(-2.457329) q[2];
sx q[2];
rz(-2.4917045) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.2838973) q[1];
sx q[1];
rz(-1.1313492) q[1];
sx q[1];
rz(1.6698014) q[1];
rz(-pi) q[2];
rz(-2.6716258) q[3];
sx q[3];
rz(-0.50438687) q[3];
sx q[3];
rz(-2.4117163) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.38301864) q[2];
sx q[2];
rz(-0.16123161) q[2];
sx q[2];
rz(2.8733011) q[2];
rz(2.7475083) q[3];
sx q[3];
rz(-1.2309309) q[3];
sx q[3];
rz(0.18850732) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2878993) q[0];
sx q[0];
rz(-0.51369602) q[0];
sx q[0];
rz(-2.7365141) q[0];
rz(-0.69008094) q[1];
sx q[1];
rz(-1.9837374) q[1];
sx q[1];
rz(2.4437723) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.47397428) q[0];
sx q[0];
rz(-1.4809429) q[0];
sx q[0];
rz(3.1180179) q[0];
rz(1.9525098) q[2];
sx q[2];
rz(-1.0587947) q[2];
sx q[2];
rz(2.2005759) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.2600425) q[1];
sx q[1];
rz(-1.3870194) q[1];
sx q[1];
rz(1.0764513) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1504193) q[3];
sx q[3];
rz(-1.8432518) q[3];
sx q[3];
rz(-3.0416995) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.71965376) q[2];
sx q[2];
rz(-1.3976588) q[2];
sx q[2];
rz(1.5092124) q[2];
rz(0.40431067) q[3];
sx q[3];
rz(-2.4590838) q[3];
sx q[3];
rz(1.6633165) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25800911) q[0];
sx q[0];
rz(-1.3556577) q[0];
sx q[0];
rz(2.1160545) q[0];
rz(2.569596) q[1];
sx q[1];
rz(-1.0943741) q[1];
sx q[1];
rz(2.5122723) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48110163) q[0];
sx q[0];
rz(-1.8697303) q[0];
sx q[0];
rz(-2.7116508) q[0];
rz(-1.8848558) q[2];
sx q[2];
rz(-0.70054189) q[2];
sx q[2];
rz(-0.26740057) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.8992936) q[1];
sx q[1];
rz(-1.0596024) q[1];
sx q[1];
rz(-1.4392698) q[1];
rz(2.1719645) q[3];
sx q[3];
rz(-1.789635) q[3];
sx q[3];
rz(2.9588483) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.5489674) q[2];
sx q[2];
rz(-2.7719345) q[2];
sx q[2];
rz(2.7992451) q[2];
rz(1.7957548) q[3];
sx q[3];
rz(-1.6941518) q[3];
sx q[3];
rz(-2.9455744) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65258566) q[0];
sx q[0];
rz(-1.2359897) q[0];
sx q[0];
rz(0.18519369) q[0];
rz(-1.7350896) q[1];
sx q[1];
rz(-1.0909189) q[1];
sx q[1];
rz(-1.3669744) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0104116) q[0];
sx q[0];
rz(-0.62958065) q[0];
sx q[0];
rz(-0.38139831) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8390859) q[2];
sx q[2];
rz(-2.2777646) q[2];
sx q[2];
rz(2.6940341) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.5609315) q[1];
sx q[1];
rz(-1.7226189) q[1];
sx q[1];
rz(0.092756943) q[1];
rz(-pi) q[2];
rz(1.9983415) q[3];
sx q[3];
rz(-1.0435836) q[3];
sx q[3];
rz(-2.9851819) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.8967445) q[2];
sx q[2];
rz(-0.5534133) q[2];
sx q[2];
rz(-0.883376) q[2];
rz(-1.4128489) q[3];
sx q[3];
rz(-2.4491375) q[3];
sx q[3];
rz(0.81378716) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25154034) q[0];
sx q[0];
rz(-0.10245704) q[0];
sx q[0];
rz(1.2782156) q[0];
rz(-3.1037519) q[1];
sx q[1];
rz(-2.3262639) q[1];
sx q[1];
rz(1.7657123) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28988518) q[0];
sx q[0];
rz(-2.2635248) q[0];
sx q[0];
rz(-0.38498199) q[0];
rz(-pi) q[1];
rz(1.6872348) q[2];
sx q[2];
rz(-0.47633007) q[2];
sx q[2];
rz(-1.4637228) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.5232877) q[1];
sx q[1];
rz(-1.4555281) q[1];
sx q[1];
rz(-1.6769888) q[1];
rz(-pi) q[2];
rz(-1.2405346) q[3];
sx q[3];
rz(-1.6263279) q[3];
sx q[3];
rz(-2.8020669) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.6841131) q[2];
sx q[2];
rz(-0.71762466) q[2];
sx q[2];
rz(0.73105556) q[2];
rz(0.11080065) q[3];
sx q[3];
rz(-1.585107) q[3];
sx q[3];
rz(-2.4462162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5450127) q[0];
sx q[0];
rz(-0.87485635) q[0];
sx q[0];
rz(0.73356432) q[0];
rz(0.60797524) q[1];
sx q[1];
rz(-1.9476451) q[1];
sx q[1];
rz(0.2342934) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9600784) q[0];
sx q[0];
rz(-2.5900822) q[0];
sx q[0];
rz(-1.0612723) q[0];
rz(-pi) q[1];
rz(-2.2393164) q[2];
sx q[2];
rz(-2.1736645) q[2];
sx q[2];
rz(-0.55842802) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.5768847) q[1];
sx q[1];
rz(-1.7600093) q[1];
sx q[1];
rz(-2.3343759) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0571363) q[3];
sx q[3];
rz(-2.497695) q[3];
sx q[3];
rz(2.7542496) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(3.0155448) q[2];
sx q[2];
rz(-1.8323106) q[2];
sx q[2];
rz(-1.7162494) q[2];
rz(1.4632633) q[3];
sx q[3];
rz(-2.3571456) q[3];
sx q[3];
rz(-1.4956168) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.54206806) q[0];
sx q[0];
rz(-0.33518377) q[0];
sx q[0];
rz(1.9482127) q[0];
rz(1.2606196) q[1];
sx q[1];
rz(-1.7766989) q[1];
sx q[1];
rz(0.9448005) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1244922) q[0];
sx q[0];
rz(-1.1860577) q[0];
sx q[0];
rz(-0.7229294) q[0];
x q[1];
rz(-2.8464727) q[2];
sx q[2];
rz(-1.1913538) q[2];
sx q[2];
rz(1.6377359) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.8560801) q[1];
sx q[1];
rz(-1.7417007) q[1];
sx q[1];
rz(2.8588365) q[1];
rz(-pi) q[2];
rz(-1.2148083) q[3];
sx q[3];
rz(-1.755852) q[3];
sx q[3];
rz(-1.9389648) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.9251359) q[2];
sx q[2];
rz(-1.0711121) q[2];
sx q[2];
rz(2.8533868) q[2];
rz(-2.6618585) q[3];
sx q[3];
rz(-1.0498485) q[3];
sx q[3];
rz(1.5079927) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0297246) q[0];
sx q[0];
rz(-2.8653963) q[0];
sx q[0];
rz(0.9129886) q[0];
rz(-0.37462014) q[1];
sx q[1];
rz(-1.7381784) q[1];
sx q[1];
rz(0.8909117) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4961632) q[0];
sx q[0];
rz(-1.2362288) q[0];
sx q[0];
rz(1.8339001) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.183379) q[2];
sx q[2];
rz(-1.3365067) q[2];
sx q[2];
rz(-0.9226375) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.3112463) q[1];
sx q[1];
rz(-1.3724569) q[1];
sx q[1];
rz(0.34459181) q[1];
x q[2];
rz(0.18817801) q[3];
sx q[3];
rz(-1.4965701) q[3];
sx q[3];
rz(1.0159514) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.9051819) q[2];
sx q[2];
rz(-0.85835251) q[2];
sx q[2];
rz(2.6399844) q[2];
rz(-1.2891399) q[3];
sx q[3];
rz(-1.4533549) q[3];
sx q[3];
rz(-0.44617173) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.63507737) q[0];
sx q[0];
rz(-1.7000533) q[0];
sx q[0];
rz(0.62361367) q[0];
rz(2.0093909) q[1];
sx q[1];
rz(-2.3846346) q[1];
sx q[1];
rz(0.089288575) q[1];
rz(-2.6272527) q[2];
sx q[2];
rz(-2.8375576) q[2];
sx q[2];
rz(0.67630771) q[2];
rz(3.0922079) q[3];
sx q[3];
rz(-2.1168843) q[3];
sx q[3];
rz(-2.4612853) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
