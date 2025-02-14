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
rz(1.1627816) q[0];
sx q[0];
rz(5.3428234) q[0];
sx q[0];
rz(12.760395) q[0];
rz(2.5972875) q[1];
sx q[1];
rz(-1.5736009) q[1];
sx q[1];
rz(0.1967217) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.17795632) q[0];
sx q[0];
rz(-0.58695176) q[0];
sx q[0];
rz(-1.3474083) q[0];
rz(-pi) q[1];
x q[1];
rz(0.098539515) q[2];
sx q[2];
rz(-0.218501) q[2];
sx q[2];
rz(-0.69536415) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.3595092) q[1];
sx q[1];
rz(-1.8386786) q[1];
sx q[1];
rz(-2.8102576) q[1];
rz(-pi) q[2];
rz(1.807735) q[3];
sx q[3];
rz(-0.58962017) q[3];
sx q[3];
rz(0.20649621) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.0804312) q[2];
sx q[2];
rz(-1.72074) q[2];
sx q[2];
rz(0.013193456) q[2];
rz(1.0641998) q[3];
sx q[3];
rz(-1.0366169) q[3];
sx q[3];
rz(2.9707151) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
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
rz(-1.6212293) q[0];
sx q[0];
rz(-0.55084387) q[0];
sx q[0];
rz(-0.4826104) q[0];
rz(-2.6699325) q[1];
sx q[1];
rz(-0.99716798) q[1];
sx q[1];
rz(-2.4536123) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2575098) q[0];
sx q[0];
rz(-2.7796352) q[0];
sx q[0];
rz(-2.5779526) q[0];
rz(1.8652164) q[2];
sx q[2];
rz(-2.3465119) q[2];
sx q[2];
rz(0.15995041) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(3.1222306) q[1];
sx q[1];
rz(-1.5270342) q[1];
sx q[1];
rz(-2.6701609) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7741995) q[3];
sx q[3];
rz(-1.6546094) q[3];
sx q[3];
rz(2.0544685) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.017435) q[2];
sx q[2];
rz(-1.0729125) q[2];
sx q[2];
rz(-2.9166481) q[2];
rz(-1.0791091) q[3];
sx q[3];
rz(-2.2033117) q[3];
sx q[3];
rz(-2.6911531) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9709622) q[0];
sx q[0];
rz(-0.58929515) q[0];
sx q[0];
rz(-1.851409) q[0];
rz(1.8633441) q[1];
sx q[1];
rz(-1.1561013) q[1];
sx q[1];
rz(1.2826756) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0541924) q[0];
sx q[0];
rz(-0.71119672) q[0];
sx q[0];
rz(-2.1367226) q[0];
rz(-pi) q[1];
rz(0.12357269) q[2];
sx q[2];
rz(-1.4915492) q[2];
sx q[2];
rz(2.045897) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.6911192) q[1];
sx q[1];
rz(-1.5509924) q[1];
sx q[1];
rz(3.1234689) q[1];
rz(-0.76126666) q[3];
sx q[3];
rz(-2.1652615) q[3];
sx q[3];
rz(-2.6192384) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.67905417) q[2];
sx q[2];
rz(-1.5306229) q[2];
sx q[2];
rz(-2.6666717) q[2];
rz(-1.9654407) q[3];
sx q[3];
rz(-2.2852496) q[3];
sx q[3];
rz(2.8600051) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7684286) q[0];
sx q[0];
rz(-1.7730862) q[0];
sx q[0];
rz(-0.1678634) q[0];
rz(-0.31790512) q[1];
sx q[1];
rz(-1.8524421) q[1];
sx q[1];
rz(1.0344523) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8694591) q[0];
sx q[0];
rz(-1.8328761) q[0];
sx q[0];
rz(3.1054405) q[0];
rz(-2.690763) q[2];
sx q[2];
rz(-0.73720142) q[2];
sx q[2];
rz(-0.99861713) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.75986162) q[1];
sx q[1];
rz(-1.1118366) q[1];
sx q[1];
rz(2.3242287) q[1];
rz(-pi) q[2];
rz(1.4698704) q[3];
sx q[3];
rz(-0.96138326) q[3];
sx q[3];
rz(0.17796365) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.1075403) q[2];
sx q[2];
rz(-0.62959051) q[2];
sx q[2];
rz(2.8752987) q[2];
rz(-2.9722424) q[3];
sx q[3];
rz(-1.063187) q[3];
sx q[3];
rz(-2.9625986) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.87194815) q[0];
sx q[0];
rz(-2.9565411) q[0];
sx q[0];
rz(1.8484775) q[0];
rz(3.0710908) q[1];
sx q[1];
rz(-0.87424707) q[1];
sx q[1];
rz(0.806113) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1925114) q[0];
sx q[0];
rz(-1.3638138) q[0];
sx q[0];
rz(0.1271609) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7197554) q[2];
sx q[2];
rz(-1.6276013) q[2];
sx q[2];
rz(-0.42213747) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.55743417) q[1];
sx q[1];
rz(-2.5100252) q[1];
sx q[1];
rz(-0.33051349) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6025196) q[3];
sx q[3];
rz(-1.381449) q[3];
sx q[3];
rz(0.93261485) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.9301794) q[2];
sx q[2];
rz(-1.4750007) q[2];
sx q[2];
rz(-2.831947) q[2];
rz(2.6304701) q[3];
sx q[3];
rz(-2.5680254) q[3];
sx q[3];
rz(2.3206594) q[3];
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
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9822134) q[0];
sx q[0];
rz(-2.5143304) q[0];
sx q[0];
rz(-0.43701592) q[0];
rz(0.45477319) q[1];
sx q[1];
rz(-2.3450856) q[1];
sx q[1];
rz(-2.2589267) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2734786) q[0];
sx q[0];
rz(-1.4820288) q[0];
sx q[0];
rz(-1.300439) q[0];
x q[1];
rz(-1.4910554) q[2];
sx q[2];
rz(-0.3796595) q[2];
sx q[2];
rz(0.14282957) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.5875579) q[1];
sx q[1];
rz(-1.9022501) q[1];
sx q[1];
rz(-2.0786839) q[1];
rz(-2.168247) q[3];
sx q[3];
rz(-2.4472467) q[3];
sx q[3];
rz(2.2996817) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.63558811) q[2];
sx q[2];
rz(-1.4281861) q[2];
sx q[2];
rz(2.7006855) q[2];
rz(-2.0697557) q[3];
sx q[3];
rz(-2.8863638) q[3];
sx q[3];
rz(-2.5409839) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.014024409) q[0];
sx q[0];
rz(-1.2034282) q[0];
sx q[0];
rz(-1.8593651) q[0];
rz(-2.0558689) q[1];
sx q[1];
rz(-1.6174569) q[1];
sx q[1];
rz(-1.5132743) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31103872) q[0];
sx q[0];
rz(-2.3966602) q[0];
sx q[0];
rz(5*pi/9) q[0];
rz(-pi) q[1];
rz(1.2499181) q[2];
sx q[2];
rz(-0.89628637) q[2];
sx q[2];
rz(-2.6303076) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.46678695) q[1];
sx q[1];
rz(-1.1937831) q[1];
sx q[1];
rz(0.88175718) q[1];
rz(-pi) q[2];
rz(1.9508771) q[3];
sx q[3];
rz(-2.2909559) q[3];
sx q[3];
rz(-2.5270568) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.47161237) q[2];
sx q[2];
rz(-0.70576224) q[2];
sx q[2];
rz(-0.82028779) q[2];
rz(-2.7465076) q[3];
sx q[3];
rz(-0.79212752) q[3];
sx q[3];
rz(-2.5274966) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.87913269) q[0];
sx q[0];
rz(-1.6430055) q[0];
sx q[0];
rz(2.5347128) q[0];
rz(-0.34582368) q[1];
sx q[1];
rz(-1.1864097) q[1];
sx q[1];
rz(0.80723673) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5257634) q[0];
sx q[0];
rz(-2.7836972) q[0];
sx q[0];
rz(-1.8496022) q[0];
x q[1];
rz(-2.8736595) q[2];
sx q[2];
rz(-1.180928) q[2];
sx q[2];
rz(-0.34696503) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.2507627) q[1];
sx q[1];
rz(-1.3034262) q[1];
sx q[1];
rz(1.7946589) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1764507) q[3];
sx q[3];
rz(-2.0803841) q[3];
sx q[3];
rz(-1.6965564) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.81749934) q[2];
sx q[2];
rz(-1.1439415) q[2];
sx q[2];
rz(0.49638805) q[2];
rz(3.056622) q[3];
sx q[3];
rz(-0.43046633) q[3];
sx q[3];
rz(-2.5607204) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1519415) q[0];
sx q[0];
rz(-1.4284416) q[0];
sx q[0];
rz(-2.7570643) q[0];
rz(-0.25228581) q[1];
sx q[1];
rz(-1.3832046) q[1];
sx q[1];
rz(-1.5581473) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4760732) q[0];
sx q[0];
rz(-1.5613149) q[0];
sx q[0];
rz(0.46561636) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7436888) q[2];
sx q[2];
rz(-2.845484) q[2];
sx q[2];
rz(0.95740151) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.059824316) q[1];
sx q[1];
rz(-1.5579954) q[1];
sx q[1];
rz(-2.1600385) q[1];
rz(-pi) q[2];
rz(2.2966398) q[3];
sx q[3];
rz(-2.264177) q[3];
sx q[3];
rz(-1.4385731) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.3462476) q[2];
sx q[2];
rz(-1.8833505) q[2];
sx q[2];
rz(-2.0015049) q[2];
rz(-1.358486) q[3];
sx q[3];
rz(-1.2525109) q[3];
sx q[3];
rz(-2.8235249) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.43050218) q[0];
sx q[0];
rz(-1.6999812) q[0];
sx q[0];
rz(0.39518133) q[0];
rz(1.9197865) q[1];
sx q[1];
rz(-0.8539353) q[1];
sx q[1];
rz(2.2538259) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2632653) q[0];
sx q[0];
rz(-2.3774231) q[0];
sx q[0];
rz(-2.7394638) q[0];
x q[1];
rz(-0.0035462499) q[2];
sx q[2];
rz(-1.0958091) q[2];
sx q[2];
rz(1.3925683) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.99988538) q[1];
sx q[1];
rz(-1.5962684) q[1];
sx q[1];
rz(1.9266121) q[1];
rz(-pi) q[2];
rz(-1.2492845) q[3];
sx q[3];
rz(-0.32673937) q[3];
sx q[3];
rz(0.48130166) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.6937574) q[2];
sx q[2];
rz(-1.8303904) q[2];
sx q[2];
rz(-2.5698404) q[2];
rz(-1.5104431) q[3];
sx q[3];
rz(-1.4331199) q[3];
sx q[3];
rz(1.3408287) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0157264) q[0];
sx q[0];
rz(-1.3073574) q[0];
sx q[0];
rz(0.14412185) q[0];
rz(1.1852599) q[1];
sx q[1];
rz(-2.2846501) q[1];
sx q[1];
rz(1.2774998) q[1];
rz(0.38707564) q[2];
sx q[2];
rz(-1.7995926) q[2];
sx q[2];
rz(2.8298701) q[2];
rz(-1.2492467) q[3];
sx q[3];
rz(-0.95148181) q[3];
sx q[3];
rz(-0.8368809) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
