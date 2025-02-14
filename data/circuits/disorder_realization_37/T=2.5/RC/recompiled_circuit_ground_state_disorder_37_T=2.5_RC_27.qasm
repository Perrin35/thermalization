OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.3545544) q[0];
sx q[0];
rz(0.7631425) q[0];
sx q[0];
rz(7.2416303) q[0];
rz(-2.7148442) q[1];
sx q[1];
rz(-2.7172105) q[1];
sx q[1];
rz(0.97270614) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5066507) q[0];
sx q[0];
rz(-1.4774862) q[0];
sx q[0];
rz(-0.20185457) q[0];
rz(-pi) q[1];
x q[1];
rz(0.17225687) q[2];
sx q[2];
rz(-1.63802) q[2];
sx q[2];
rz(-2.7329815) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.7643775) q[1];
sx q[1];
rz(-0.7684349) q[1];
sx q[1];
rz(0.37698104) q[1];
rz(-pi) q[2];
rz(2.8593721) q[3];
sx q[3];
rz(-1.422554) q[3];
sx q[3];
rz(-2.6765598) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-3.0613609) q[2];
sx q[2];
rz(-1.1876567) q[2];
sx q[2];
rz(2.3415671) q[2];
rz(-3.0112265) q[3];
sx q[3];
rz(-2.3733807) q[3];
sx q[3];
rz(0.31726328) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.80113634) q[0];
sx q[0];
rz(-1.9213333) q[0];
sx q[0];
rz(-0.98943797) q[0];
rz(-0.89775741) q[1];
sx q[1];
rz(-1.0549301) q[1];
sx q[1];
rz(-2.3602233) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5637276) q[0];
sx q[0];
rz(-1.0672309) q[0];
sx q[0];
rz(2.9211852) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4825253) q[2];
sx q[2];
rz(-1.3376184) q[2];
sx q[2];
rz(-0.424343) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.20618379) q[1];
sx q[1];
rz(-2.0061135) q[1];
sx q[1];
rz(1.0401506) q[1];
rz(-pi) q[2];
rz(-1.2739232) q[3];
sx q[3];
rz(-1.9386995) q[3];
sx q[3];
rz(-1.5987336) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(3.117131) q[2];
sx q[2];
rz(-1.7639561) q[2];
sx q[2];
rz(-2.3632939) q[2];
rz(-0.69027573) q[3];
sx q[3];
rz(-2.2199151) q[3];
sx q[3];
rz(1.5604304) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
rz(0.86163259) q[0];
sx q[0];
rz(-2.3065688) q[0];
sx q[0];
rz(0.41197187) q[0];
rz(0.95097923) q[1];
sx q[1];
rz(-0.6528267) q[1];
sx q[1];
rz(-1.9297809) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.73746956) q[0];
sx q[0];
rz(-1.5323191) q[0];
sx q[0];
rz(3.1079556) q[0];
x q[1];
rz(1.1583382) q[2];
sx q[2];
rz(-0.24246267) q[2];
sx q[2];
rz(-1.2810156) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.131457) q[1];
sx q[1];
rz(-1.5431511) q[1];
sx q[1];
rz(2.4316039) q[1];
rz(-1.4097628) q[3];
sx q[3];
rz(-1.9317683) q[3];
sx q[3];
rz(1.4396813) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.1872306) q[2];
sx q[2];
rz(-0.53202859) q[2];
sx q[2];
rz(1.6383891) q[2];
rz(-0.040180834) q[3];
sx q[3];
rz(-2.2153722) q[3];
sx q[3];
rz(-0.23475501) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9784341) q[0];
sx q[0];
rz(-1.6272767) q[0];
sx q[0];
rz(-0.12983313) q[0];
rz(1.8561329) q[1];
sx q[1];
rz(-1.9655656) q[1];
sx q[1];
rz(1.0349549) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1945788) q[0];
sx q[0];
rz(-1.1348551) q[0];
sx q[0];
rz(2.2313124) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7096388) q[2];
sx q[2];
rz(-2.4731257) q[2];
sx q[2];
rz(-0.61651361) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.6432727) q[1];
sx q[1];
rz(-2.8788212) q[1];
sx q[1];
rz(-2.8064125) q[1];
rz(-1.0134936) q[3];
sx q[3];
rz(-2.3689007) q[3];
sx q[3];
rz(2.0183636) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.3639565) q[2];
sx q[2];
rz(-1.3866321) q[2];
sx q[2];
rz(-0.094495471) q[2];
rz(-2.7206521) q[3];
sx q[3];
rz(-1.7178444) q[3];
sx q[3];
rz(-1.4360992) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7218551) q[0];
sx q[0];
rz(-0.57191816) q[0];
sx q[0];
rz(3.1255334) q[0];
rz(0.84469604) q[1];
sx q[1];
rz(-2.7222996) q[1];
sx q[1];
rz(0.90763456) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1354271) q[0];
sx q[0];
rz(-2.3434397) q[0];
sx q[0];
rz(-2.8291875) q[0];
rz(-pi) q[1];
rz(-2.0727022) q[2];
sx q[2];
rz(-0.47119432) q[2];
sx q[2];
rz(-2.1387177) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.076489623) q[1];
sx q[1];
rz(-1.8999083) q[1];
sx q[1];
rz(2.0338716) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5441555) q[3];
sx q[3];
rz(-1.2268042) q[3];
sx q[3];
rz(2.2709803) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.57546651) q[2];
sx q[2];
rz(-2.0822272) q[2];
sx q[2];
rz(0.063684138) q[2];
rz(0.56695402) q[3];
sx q[3];
rz(-2.3072115) q[3];
sx q[3];
rz(-0.59757346) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31396922) q[0];
sx q[0];
rz(-0.14179985) q[0];
sx q[0];
rz(1.1578479) q[0];
rz(-1.4843548) q[1];
sx q[1];
rz(-1.4446222) q[1];
sx q[1];
rz(2.8725913) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4542352) q[0];
sx q[0];
rz(-1.4299862) q[0];
sx q[0];
rz(-1.6160377) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8431191) q[2];
sx q[2];
rz(-1.6629873) q[2];
sx q[2];
rz(2.628679) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.288876) q[1];
sx q[1];
rz(-1.3914795) q[1];
sx q[1];
rz(-0.10097558) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7888467) q[3];
sx q[3];
rz(-1.9761128) q[3];
sx q[3];
rz(0.97783711) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.3645662) q[2];
sx q[2];
rz(-2.624161) q[2];
sx q[2];
rz(-2.2606134) q[2];
rz(3.0173054) q[3];
sx q[3];
rz(-1.4275987) q[3];
sx q[3];
rz(-0.65829268) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.86150259) q[0];
sx q[0];
rz(-2.8741591) q[0];
sx q[0];
rz(-2.4524443) q[0];
rz(0.46734494) q[1];
sx q[1];
rz(-2.5655949) q[1];
sx q[1];
rz(-1.0980094) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.91749597) q[0];
sx q[0];
rz(-2.9584329) q[0];
sx q[0];
rz(0.44395019) q[0];
x q[1];
rz(-1.271209) q[2];
sx q[2];
rz(-2.0986084) q[2];
sx q[2];
rz(2.6115548) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.5913691) q[1];
sx q[1];
rz(-2.6003631) q[1];
sx q[1];
rz(3.079097) q[1];
rz(-0.50366537) q[3];
sx q[3];
rz(-2.1379469) q[3];
sx q[3];
rz(-0.43225542) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.565862) q[2];
sx q[2];
rz(-0.87811676) q[2];
sx q[2];
rz(-0.6305387) q[2];
rz(-2.5343043) q[3];
sx q[3];
rz(-0.90689617) q[3];
sx q[3];
rz(-1.0926532) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31834114) q[0];
sx q[0];
rz(-0.14597758) q[0];
sx q[0];
rz(1.7952221) q[0];
rz(1.775555) q[1];
sx q[1];
rz(-0.64907688) q[1];
sx q[1];
rz(2.5128561) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7766984) q[0];
sx q[0];
rz(-1.3140251) q[0];
sx q[0];
rz(2.5904694) q[0];
rz(-pi) q[1];
rz(-0.95405719) q[2];
sx q[2];
rz(-0.77191691) q[2];
sx q[2];
rz(1.3961015) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.4688411) q[1];
sx q[1];
rz(-0.76460014) q[1];
sx q[1];
rz(0.73672626) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7698742) q[3];
sx q[3];
rz(-1.645483) q[3];
sx q[3];
rz(1.5499566) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.2748572) q[2];
sx q[2];
rz(-1.3254415) q[2];
sx q[2];
rz(-2.2692915) q[2];
rz(-2.9288779) q[3];
sx q[3];
rz(-1.3958967) q[3];
sx q[3];
rz(0.58342903) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.23671959) q[0];
sx q[0];
rz(-1.6107591) q[0];
sx q[0];
rz(1.2637631) q[0];
rz(1.0172552) q[1];
sx q[1];
rz(-2.2433498) q[1];
sx q[1];
rz(0.57473007) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0145469) q[0];
sx q[0];
rz(-1.0406063) q[0];
sx q[0];
rz(1.4958025) q[0];
rz(-pi) q[1];
rz(-2.9024486) q[2];
sx q[2];
rz(-0.84296339) q[2];
sx q[2];
rz(0.69918888) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.3077449) q[1];
sx q[1];
rz(-2.7258432) q[1];
sx q[1];
rz(-0.47885696) q[1];
rz(-pi) q[2];
rz(1.0887126) q[3];
sx q[3];
rz(-1.3441372) q[3];
sx q[3];
rz(2.4149293) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.7175674) q[2];
sx q[2];
rz(-1.8342476) q[2];
sx q[2];
rz(3.1094816) q[2];
rz(-1.9854246) q[3];
sx q[3];
rz(-0.38438946) q[3];
sx q[3];
rz(1.2110075) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0403989) q[0];
sx q[0];
rz(-0.54553425) q[0];
sx q[0];
rz(1.1241166) q[0];
rz(-0.54565564) q[1];
sx q[1];
rz(-2.7898495) q[1];
sx q[1];
rz(2.5901332) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7246147) q[0];
sx q[0];
rz(-1.7584086) q[0];
sx q[0];
rz(-3.0155988) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7573171) q[2];
sx q[2];
rz(-0.2787481) q[2];
sx q[2];
rz(-2.2997625) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.3797326) q[1];
sx q[1];
rz(-1.1071536) q[1];
sx q[1];
rz(-2.597888) q[1];
rz(-pi) q[2];
rz(1.7343246) q[3];
sx q[3];
rz(-1.8732605) q[3];
sx q[3];
rz(2.1238984) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.55190279) q[2];
sx q[2];
rz(-0.51108131) q[2];
sx q[2];
rz(1.7299293) q[2];
rz(-2.3949413) q[3];
sx q[3];
rz(-1.6801497) q[3];
sx q[3];
rz(-0.97486973) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.830037) q[0];
sx q[0];
rz(-1.5811601) q[0];
sx q[0];
rz(1.5720221) q[0];
rz(-0.39881067) q[1];
sx q[1];
rz(-1.8067982) q[1];
sx q[1];
rz(-1.6511818) q[1];
rz(-1.082303) q[2];
sx q[2];
rz(-1.563896) q[2];
sx q[2];
rz(0.36120216) q[2];
rz(-1.87181) q[3];
sx q[3];
rz(-2.3475937) q[3];
sx q[3];
rz(0.54318843) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
