OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.78703824) q[0];
sx q[0];
rz(-0.7631425) q[0];
sx q[0];
rz(0.95844498) q[0];
rz(0.42674843) q[1];
sx q[1];
rz(-0.42438212) q[1];
sx q[1];
rz(-0.97270614) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6329968) q[0];
sx q[0];
rz(-2.9194814) q[0];
sx q[0];
rz(2.7048777) q[0];
rz(-pi) q[1];
rz(-0.3742674) q[2];
sx q[2];
rz(-2.9568045) q[2];
sx q[2];
rz(-1.6109465) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.083746759) q[1];
sx q[1];
rz(-1.3120756) q[1];
sx q[1];
rz(2.4094635) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6493727) q[3];
sx q[3];
rz(-2.8237298) q[3];
sx q[3];
rz(1.5647056) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-3.0613609) q[2];
sx q[2];
rz(-1.1876567) q[2];
sx q[2];
rz(-0.80002552) q[2];
rz(-3.0112265) q[3];
sx q[3];
rz(-0.76821199) q[3];
sx q[3];
rz(-0.31726328) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3404563) q[0];
sx q[0];
rz(-1.2202593) q[0];
sx q[0];
rz(-0.98943797) q[0];
rz(2.2438352) q[1];
sx q[1];
rz(-2.0866626) q[1];
sx q[1];
rz(2.3602233) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1290478) q[0];
sx q[0];
rz(-2.5957286) q[0];
sx q[0];
rz(1.9485628) q[0];
rz(-pi) q[1];
rz(1.6590674) q[2];
sx q[2];
rz(-1.8039743) q[2];
sx q[2];
rz(2.7172497) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.9877517) q[1];
sx q[1];
rz(-0.67286009) q[1];
sx q[1];
rz(-0.82760906) q[1];
rz(-pi) q[2];
rz(-2.7584293) q[3];
sx q[3];
rz(-1.8472611) q[3];
sx q[3];
rz(-0.081646669) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.024461688) q[2];
sx q[2];
rz(-1.3776366) q[2];
sx q[2];
rz(0.77829877) q[2];
rz(0.69027573) q[3];
sx q[3];
rz(-0.92167753) q[3];
sx q[3];
rz(1.5604304) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.86163259) q[0];
sx q[0];
rz(-0.83502382) q[0];
sx q[0];
rz(2.7296208) q[0];
rz(-2.1906134) q[1];
sx q[1];
rz(-0.6528267) q[1];
sx q[1];
rz(1.2118118) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4561583) q[0];
sx q[0];
rz(-3.0904909) q[0];
sx q[0];
rz(-0.85275485) q[0];
x q[1];
rz(1.7936208) q[2];
sx q[2];
rz(-1.6671902) q[2];
sx q[2];
rz(-2.4501462) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.6784994) q[1];
sx q[1];
rz(-2.2804567) q[1];
sx q[1];
rz(1.5343496) q[1];
rz(-2.7762967) q[3];
sx q[3];
rz(-1.72137) q[3];
sx q[3];
rz(0.07380658) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.1872306) q[2];
sx q[2];
rz(-2.6095641) q[2];
sx q[2];
rz(1.6383891) q[2];
rz(0.040180834) q[3];
sx q[3];
rz(-0.92622042) q[3];
sx q[3];
rz(-0.23475501) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9784341) q[0];
sx q[0];
rz(-1.6272767) q[0];
sx q[0];
rz(3.0117595) q[0];
rz(1.8561329) q[1];
sx q[1];
rz(-1.1760271) q[1];
sx q[1];
rz(2.1066378) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.059211123) q[0];
sx q[0];
rz(-2.1605345) q[0];
sx q[0];
rz(-0.53296169) q[0];
x q[1];
rz(0.1088684) q[2];
sx q[2];
rz(-0.90990674) q[2];
sx q[2];
rz(-2.3488597) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.49832) q[1];
sx q[1];
rz(-2.8788212) q[1];
sx q[1];
rz(0.33518016) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.47607039) q[3];
sx q[3];
rz(-0.93671533) q[3];
sx q[3];
rz(0.40704029) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.7776362) q[2];
sx q[2];
rz(-1.3866321) q[2];
sx q[2];
rz(-3.0470972) q[2];
rz(2.7206521) q[3];
sx q[3];
rz(-1.4237483) q[3];
sx q[3];
rz(-1.4360992) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7218551) q[0];
sx q[0];
rz(-0.57191816) q[0];
sx q[0];
rz(0.016059248) q[0];
rz(0.84469604) q[1];
sx q[1];
rz(-2.7222996) q[1];
sx q[1];
rz(0.90763456) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3428872) q[0];
sx q[0];
rz(-1.3488975) q[0];
sx q[0];
rz(-2.368244) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0688905) q[2];
sx q[2];
rz(-0.47119432) q[2];
sx q[2];
rz(1.0028749) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.4872949) q[1];
sx q[1];
rz(-2.0072486) q[1];
sx q[1];
rz(0.36466332) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0673712) q[3];
sx q[3];
rz(-0.34498131) q[3];
sx q[3];
rz(-0.94946194) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.57546651) q[2];
sx q[2];
rz(-1.0593654) q[2];
sx q[2];
rz(-0.063684138) q[2];
rz(-2.5746386) q[3];
sx q[3];
rz(-2.3072115) q[3];
sx q[3];
rz(2.5440192) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(2.8276234) q[0];
sx q[0];
rz(-0.14179985) q[0];
sx q[0];
rz(1.9837448) q[0];
rz(-1.4843548) q[1];
sx q[1];
rz(-1.6969705) q[1];
sx q[1];
rz(-2.8725913) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6873574) q[0];
sx q[0];
rz(-1.7116065) q[0];
sx q[0];
rz(-1.5255549) q[0];
rz(-pi) q[1];
rz(-1.2984736) q[2];
sx q[2];
rz(-1.6629873) q[2];
sx q[2];
rz(-2.628679) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.8416031) q[1];
sx q[1];
rz(-1.4714452) q[1];
sx q[1];
rz(-1.7510115) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.352746) q[3];
sx q[3];
rz(-1.1654799) q[3];
sx q[3];
rz(-2.1637555) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.77702648) q[2];
sx q[2];
rz(-0.51743162) q[2];
sx q[2];
rz(0.88097921) q[2];
rz(3.0173054) q[3];
sx q[3];
rz(-1.7139939) q[3];
sx q[3];
rz(0.65829268) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2800901) q[0];
sx q[0];
rz(-0.26743356) q[0];
sx q[0];
rz(-0.68914831) q[0];
rz(2.6742477) q[1];
sx q[1];
rz(-0.57599774) q[1];
sx q[1];
rz(2.0435832) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9257346) q[0];
sx q[0];
rz(-1.6491062) q[0];
sx q[0];
rz(-0.16574482) q[0];
rz(-2.5937366) q[2];
sx q[2];
rz(-1.3129873) q[2];
sx q[2];
rz(-1.9465035) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-3.0674379) q[1];
sx q[1];
rz(-1.6029781) q[1];
sx q[1];
rz(0.54036714) q[1];
rz(-pi) q[2];
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
rz(-2.2634759) q[2];
sx q[2];
rz(-2.511054) q[2];
rz(0.60728836) q[3];
sx q[3];
rz(-0.90689617) q[3];
sx q[3];
rz(-1.0926532) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(0.31834114) q[0];
sx q[0];
rz(-2.9956151) q[0];
sx q[0];
rz(-1.3463705) q[0];
rz(-1.3660376) q[1];
sx q[1];
rz(-0.64907688) q[1];
sx q[1];
rz(-0.62873658) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3607488) q[0];
sx q[0];
rz(-2.1018902) q[0];
sx q[0];
rz(1.2718334) q[0];
rz(-pi) q[1];
rz(-2.2419078) q[2];
sx q[2];
rz(-1.9860528) q[2];
sx q[2];
rz(-0.64476162) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.4688411) q[1];
sx q[1];
rz(-2.3769925) q[1];
sx q[1];
rz(-0.73672626) q[1];
rz(0.076185779) q[3];
sx q[3];
rz(-1.7693118) q[3];
sx q[3];
rz(0.0057868231) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.8667355) q[2];
sx q[2];
rz(-1.3254415) q[2];
sx q[2];
rz(-2.2692915) q[2];
rz(2.9288779) q[3];
sx q[3];
rz(-1.3958967) q[3];
sx q[3];
rz(-0.58342903) q[3];
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
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9048731) q[0];
sx q[0];
rz(-1.6107591) q[0];
sx q[0];
rz(-1.8778296) q[0];
rz(2.1243375) q[1];
sx q[1];
rz(-2.2433498) q[1];
sx q[1];
rz(-0.57473007) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0145469) q[0];
sx q[0];
rz(-1.0406063) q[0];
sx q[0];
rz(-1.4958025) q[0];
rz(1.310964) q[2];
sx q[2];
rz(-0.75922478) q[2];
sx q[2];
rz(0.34789839) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.8239261) q[1];
sx q[1];
rz(-1.2041908) q[1];
sx q[1];
rz(1.7714785) q[1];
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
rz(-0.42402521) q[2];
sx q[2];
rz(-1.307345) q[2];
sx q[2];
rz(-3.1094816) q[2];
rz(1.1561681) q[3];
sx q[3];
rz(-2.7572032) q[3];
sx q[3];
rz(1.9305852) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0403989) q[0];
sx q[0];
rz(-0.54553425) q[0];
sx q[0];
rz(2.0174761) q[0];
rz(2.595937) q[1];
sx q[1];
rz(-0.35174313) q[1];
sx q[1];
rz(0.55145946) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.41697793) q[0];
sx q[0];
rz(-1.383184) q[0];
sx q[0];
rz(-0.12599385) q[0];
rz(-0.053023382) q[2];
sx q[2];
rz(-1.2970088) q[2];
sx q[2];
rz(-1.0356569) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.76186) q[1];
sx q[1];
rz(-1.1071536) q[1];
sx q[1];
rz(2.597888) q[1];
x q[2];
rz(-2.6607108) q[3];
sx q[3];
rz(-2.7989498) q[3];
sx q[3];
rz(0.51183701) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.5896899) q[2];
sx q[2];
rz(-0.51108131) q[2];
sx q[2];
rz(-1.7299293) q[2];
rz(2.3949413) q[3];
sx q[3];
rz(-1.6801497) q[3];
sx q[3];
rz(-2.1667229) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.830037) q[0];
sx q[0];
rz(-1.5604326) q[0];
sx q[0];
rz(-1.5695705) q[0];
rz(-0.39881067) q[1];
sx q[1];
rz(-1.8067982) q[1];
sx q[1];
rz(-1.6511818) q[1];
rz(-2.0592897) q[2];
sx q[2];
rz(-1.5776967) q[2];
sx q[2];
rz(-2.7803905) q[2];
rz(1.87181) q[3];
sx q[3];
rz(-0.79399899) q[3];
sx q[3];
rz(-2.5984042) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
