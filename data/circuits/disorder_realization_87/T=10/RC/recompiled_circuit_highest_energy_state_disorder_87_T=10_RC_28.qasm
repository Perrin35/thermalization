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
rz(1.4662161) q[0];
sx q[0];
rz(-0.94453064) q[0];
sx q[0];
rz(0.20139995) q[0];
rz(-2.0693076) q[1];
sx q[1];
rz(-2.3787002) q[1];
sx q[1];
rz(1.720517) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11036377) q[0];
sx q[0];
rz(-2.5284934) q[0];
sx q[0];
rz(1.4754063) q[0];
rz(-pi) q[1];
rz(0.95008738) q[2];
sx q[2];
rz(-2.3664775) q[2];
sx q[2];
rz(0.24573869) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.0396467) q[1];
sx q[1];
rz(-1.1803487) q[1];
sx q[1];
rz(-0.70266188) q[1];
x q[2];
rz(2.8086923) q[3];
sx q[3];
rz(-1.1984636) q[3];
sx q[3];
rz(-1.2941192) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.034915514) q[2];
sx q[2];
rz(-2.3372529) q[2];
sx q[2];
rz(-2.8625281) q[2];
rz(-1.8883102) q[3];
sx q[3];
rz(-2.8918355) q[3];
sx q[3];
rz(2.9899924) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6107553) q[0];
sx q[0];
rz(-2.294367) q[0];
sx q[0];
rz(-2.8952059) q[0];
rz(1.1950182) q[1];
sx q[1];
rz(-2.2476826) q[1];
sx q[1];
rz(1.5066719) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4621861) q[0];
sx q[0];
rz(-0.52216086) q[0];
sx q[0];
rz(-1.0617816) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4343403) q[2];
sx q[2];
rz(-2.8359957) q[2];
sx q[2];
rz(-1.9836677) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.23562787) q[1];
sx q[1];
rz(-1.2081138) q[1];
sx q[1];
rz(-3.079804) q[1];
x q[2];
rz(0.43466062) q[3];
sx q[3];
rz(-2.0512329) q[3];
sx q[3];
rz(-2.2312589) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.7404777) q[2];
sx q[2];
rz(-0.97492188) q[2];
sx q[2];
rz(1.9096036) q[2];
rz(2.039382) q[3];
sx q[3];
rz(-0.004318459) q[3];
sx q[3];
rz(1.3365041) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44322893) q[0];
sx q[0];
rz(-2.4446428) q[0];
sx q[0];
rz(2.2337636) q[0];
rz(-0.56697956) q[1];
sx q[1];
rz(-2.1588529) q[1];
sx q[1];
rz(0.10279113) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.028057) q[0];
sx q[0];
rz(-0.60657078) q[0];
sx q[0];
rz(-2.7675178) q[0];
x q[1];
rz(3.0605761) q[2];
sx q[2];
rz(-1.6747432) q[2];
sx q[2];
rz(2.8250717) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.56270617) q[1];
sx q[1];
rz(-1.1670615) q[1];
sx q[1];
rz(0.63376928) q[1];
rz(2.0538883) q[3];
sx q[3];
rz(-2.5317041) q[3];
sx q[3];
rz(2.1263378) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.28321442) q[2];
sx q[2];
rz(-0.80867043) q[2];
sx q[2];
rz(-1.4374479) q[2];
rz(-0.39250675) q[3];
sx q[3];
rz(-1.1350574) q[3];
sx q[3];
rz(0.2984305) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0933541) q[0];
sx q[0];
rz(-1.3627351) q[0];
sx q[0];
rz(-1.1887953) q[0];
rz(-0.28611723) q[1];
sx q[1];
rz(-0.61758271) q[1];
sx q[1];
rz(-1.6292705) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8522569) q[0];
sx q[0];
rz(-1.8484383) q[0];
sx q[0];
rz(-1.0278948) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3740923) q[2];
sx q[2];
rz(-1.8165117) q[2];
sx q[2];
rz(-1.9227288) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.0551853) q[1];
sx q[1];
rz(-2.1595104) q[1];
sx q[1];
rz(2.9179395) q[1];
x q[2];
rz(1.063594) q[3];
sx q[3];
rz(-1.385687) q[3];
sx q[3];
rz(2.0349658) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.7901223) q[2];
sx q[2];
rz(-1.0834379) q[2];
sx q[2];
rz(-0.67698395) q[2];
rz(1.2466768) q[3];
sx q[3];
rz(-1.8691984) q[3];
sx q[3];
rz(-1.6251132) q[3];
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
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.93382728) q[0];
sx q[0];
rz(-2.8849869) q[0];
sx q[0];
rz(0.024913464) q[0];
rz(1.0554396) q[1];
sx q[1];
rz(-2.1565304) q[1];
sx q[1];
rz(-2.8581462) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0830529) q[0];
sx q[0];
rz(-1.3136567) q[0];
sx q[0];
rz(-1.8399946) q[0];
rz(-2.1589958) q[2];
sx q[2];
rz(-0.23410205) q[2];
sx q[2];
rz(-0.38123075) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.39549669) q[1];
sx q[1];
rz(-1.0701792) q[1];
sx q[1];
rz(0.034828111) q[1];
x q[2];
rz(-2.2437566) q[3];
sx q[3];
rz(-1.3011609) q[3];
sx q[3];
rz(-1.8561038) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.31263605) q[2];
sx q[2];
rz(-0.41746155) q[2];
sx q[2];
rz(0.32583315) q[2];
rz(-1.2827986) q[3];
sx q[3];
rz(-1.157434) q[3];
sx q[3];
rz(1.5939943) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7376937) q[0];
sx q[0];
rz(-1.6588545) q[0];
sx q[0];
rz(0.37495908) q[0];
rz(0.47863475) q[1];
sx q[1];
rz(-2.4691983) q[1];
sx q[1];
rz(-2.7714444) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2787196) q[0];
sx q[0];
rz(-1.5939043) q[0];
sx q[0];
rz(1.5824759) q[0];
x q[1];
rz(-0.13932087) q[2];
sx q[2];
rz(-1.3874386) q[2];
sx q[2];
rz(-2.574711) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.3016175) q[1];
sx q[1];
rz(-0.8704005) q[1];
sx q[1];
rz(2.5216334) q[1];
rz(-pi) q[2];
rz(0.42780723) q[3];
sx q[3];
rz(-1.5662346) q[3];
sx q[3];
rz(-3.004194) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.7438573) q[2];
sx q[2];
rz(-2.9501259) q[2];
sx q[2];
rz(-1.7605304) q[2];
rz(-1.0487652) q[3];
sx q[3];
rz(-1.7966813) q[3];
sx q[3];
rz(-0.63961187) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25310707) q[0];
sx q[0];
rz(-0.83919224) q[0];
sx q[0];
rz(-2.2174368) q[0];
rz(-2.2249075) q[1];
sx q[1];
rz(-1.0662096) q[1];
sx q[1];
rz(-1.4013269) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9144672) q[0];
sx q[0];
rz(-2.0642529) q[0];
sx q[0];
rz(1.5052133) q[0];
rz(-pi) q[1];
rz(0.22221128) q[2];
sx q[2];
rz(-2.4743818) q[2];
sx q[2];
rz(1.3363584) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.37937134) q[1];
sx q[1];
rz(-1.4109542) q[1];
sx q[1];
rz(-0.74121468) q[1];
rz(-pi) q[2];
rz(-1.9643149) q[3];
sx q[3];
rz(-0.52214185) q[3];
sx q[3];
rz(0.71557426) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.2283198) q[2];
sx q[2];
rz(-1.6797804) q[2];
sx q[2];
rz(1.224996) q[2];
rz(0.0082958881) q[3];
sx q[3];
rz(-0.010992916) q[3];
sx q[3];
rz(-2.2744961) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55815721) q[0];
sx q[0];
rz(-2.3112516) q[0];
sx q[0];
rz(-0.48026568) q[0];
rz(-2.9971314) q[1];
sx q[1];
rz(-2.2339349) q[1];
sx q[1];
rz(2.2772148) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7378993) q[0];
sx q[0];
rz(-1.4072352) q[0];
sx q[0];
rz(-3.0944194) q[0];
rz(-pi) q[1];
x q[1];
rz(0.094674663) q[2];
sx q[2];
rz(-1.6045158) q[2];
sx q[2];
rz(-1.124749) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.0675812) q[1];
sx q[1];
rz(-1.0339104) q[1];
sx q[1];
rz(-2.8606961) q[1];
x q[2];
rz(0.79473991) q[3];
sx q[3];
rz(-0.38215548) q[3];
sx q[3];
rz(-2.4433745) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.9353443) q[2];
sx q[2];
rz(-2.1312921) q[2];
sx q[2];
rz(-0.63507357) q[2];
rz(2.5687929) q[3];
sx q[3];
rz(-1.7856995) q[3];
sx q[3];
rz(-2.2662381) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11005814) q[0];
sx q[0];
rz(-0.77339554) q[0];
sx q[0];
rz(-3.0173259) q[0];
rz(0.36525137) q[1];
sx q[1];
rz(-1.4271913) q[1];
sx q[1];
rz(-1.431538) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3076403) q[0];
sx q[0];
rz(-0.26276428) q[0];
sx q[0];
rz(-2.3961179) q[0];
rz(-pi) q[1];
rz(0.65406873) q[2];
sx q[2];
rz(-1.1032915) q[2];
sx q[2];
rz(1.9273749) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.46043432) q[1];
sx q[1];
rz(-1.2660813) q[1];
sx q[1];
rz(-2.4439993) q[1];
rz(0.67691524) q[3];
sx q[3];
rz(-0.64161086) q[3];
sx q[3];
rz(0.42719242) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.26620904) q[2];
sx q[2];
rz(-2.2201846) q[2];
sx q[2];
rz(1.9688152) q[2];
rz(1.4069936) q[3];
sx q[3];
rz(-1.809779) q[3];
sx q[3];
rz(1.9269358) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5491972) q[0];
sx q[0];
rz(-1.1347436) q[0];
sx q[0];
rz(0.59463516) q[0];
rz(1.0058962) q[1];
sx q[1];
rz(-1.9391831) q[1];
sx q[1];
rz(-0.88919052) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62458663) q[0];
sx q[0];
rz(-1.3505624) q[0];
sx q[0];
rz(3.0896679) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0709549) q[2];
sx q[2];
rz(-1.1114235) q[2];
sx q[2];
rz(0.86462155) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.31510776) q[1];
sx q[1];
rz(-1.4764305) q[1];
sx q[1];
rz(0.06006518) q[1];
rz(-pi) q[2];
rz(2.6252006) q[3];
sx q[3];
rz(-1.603873) q[3];
sx q[3];
rz(-0.75504485) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.1903926) q[2];
sx q[2];
rz(-1.1756281) q[2];
sx q[2];
rz(0.54538837) q[2];
rz(-0.96327463) q[3];
sx q[3];
rz(-1.1759718) q[3];
sx q[3];
rz(-1.6776599) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49190285) q[0];
sx q[0];
rz(-2.2408673) q[0];
sx q[0];
rz(1.3229205) q[0];
rz(-0.84359618) q[1];
sx q[1];
rz(-1.0782764) q[1];
sx q[1];
rz(-1.2596399) q[1];
rz(2.9424473) q[2];
sx q[2];
rz(-2.0024588) q[2];
sx q[2];
rz(-1.3695516) q[2];
rz(2.4424408) q[3];
sx q[3];
rz(-0.71157645) q[3];
sx q[3];
rz(-0.012307766) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
