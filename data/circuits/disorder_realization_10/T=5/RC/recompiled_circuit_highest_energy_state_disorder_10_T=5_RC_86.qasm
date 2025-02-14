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
rz(-1.7614814) q[0];
sx q[0];
rz(7.8742134) q[0];
sx q[0];
rz(5.9013517) q[0];
rz(2.7497357) q[1];
sx q[1];
rz(-1.5500103) q[1];
sx q[1];
rz(2.2171807) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59413737) q[0];
sx q[0];
rz(-2.1004803) q[0];
sx q[0];
rz(-0.87181566) q[0];
rz(-pi) q[1];
rz(0.37021356) q[2];
sx q[2];
rz(-1.5523124) q[2];
sx q[2];
rz(-0.12272515) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.83340911) q[1];
sx q[1];
rz(-1.4410958) q[1];
sx q[1];
rz(-2.7144578) q[1];
rz(-pi) q[2];
rz(-1.0409058) q[3];
sx q[3];
rz(-1.232551) q[3];
sx q[3];
rz(1.4640946) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.9367289) q[2];
sx q[2];
rz(-1.2932581) q[2];
sx q[2];
rz(-1.5667685) q[2];
rz(1.1664248) q[3];
sx q[3];
rz(-1.0650485) q[3];
sx q[3];
rz(-2.4488357) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.67398706) q[0];
sx q[0];
rz(-0.80802149) q[0];
sx q[0];
rz(-1.1760733) q[0];
rz(-0.067497079) q[1];
sx q[1];
rz(-2.564023) q[1];
sx q[1];
rz(0.31879058) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2683811) q[0];
sx q[0];
rz(-0.86641524) q[0];
sx q[0];
rz(-1.9834788) q[0];
rz(-1.3347049) q[2];
sx q[2];
rz(-1.1523048) q[2];
sx q[2];
rz(0.22007569) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.31736703) q[1];
sx q[1];
rz(-1.1928802) q[1];
sx q[1];
rz(-0.96940689) q[1];
x q[2];
rz(-1.8570921) q[3];
sx q[3];
rz(-1.4501374) q[3];
sx q[3];
rz(0.21721043) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.0789644) q[2];
sx q[2];
rz(-2.2603409) q[2];
sx q[2];
rz(2.6641565) q[2];
rz(1.7834974) q[3];
sx q[3];
rz(-0.6529468) q[3];
sx q[3];
rz(3.131598) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39182144) q[0];
sx q[0];
rz(-1.3359767) q[0];
sx q[0];
rz(1.1392449) q[0];
rz(0.40621743) q[1];
sx q[1];
rz(-1.8832877) q[1];
sx q[1];
rz(2.4635945) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8157522) q[0];
sx q[0];
rz(-1.6597972) q[0];
sx q[0];
rz(2.1584216) q[0];
rz(2.1966372) q[2];
sx q[2];
rz(-1.8650818) q[2];
sx q[2];
rz(-1.4465894) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.9443612) q[1];
sx q[1];
rz(-2.0294831) q[1];
sx q[1];
rz(-0.58802559) q[1];
rz(-2.0800679) q[3];
sx q[3];
rz(-1.8647927) q[3];
sx q[3];
rz(-1.7727838) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.71538007) q[2];
sx q[2];
rz(-1.6343626) q[2];
sx q[2];
rz(2.0646084) q[2];
rz(1.2601323) q[3];
sx q[3];
rz(-1.0271122) q[3];
sx q[3];
rz(2.8857359) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
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
rz(-0.28427163) q[0];
sx q[0];
rz(-1.9673328) q[0];
sx q[0];
rz(2.381109) q[0];
rz(-2.7898232) q[1];
sx q[1];
rz(-2.4302509) q[1];
sx q[1];
rz(0.15028353) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.12211299) q[0];
sx q[0];
rz(-1.5688591) q[0];
sx q[0];
rz(-1.5741072) q[0];
x q[1];
rz(1.9801122) q[2];
sx q[2];
rz(-2.7285353) q[2];
sx q[2];
rz(-1.5395791) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.6071668) q[1];
sx q[1];
rz(-0.98984209) q[1];
sx q[1];
rz(2.7423285) q[1];
rz(1.3129381) q[3];
sx q[3];
rz(-1.1262731) q[3];
sx q[3];
rz(0.049098102) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.030423) q[2];
sx q[2];
rz(-2.9141278) q[2];
sx q[2];
rz(2.5566768) q[2];
rz(-0.065718204) q[3];
sx q[3];
rz(-1.1394371) q[3];
sx q[3];
rz(2.0871302) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86972648) q[0];
sx q[0];
rz(-1.2249185) q[0];
sx q[0];
rz(-0.72750339) q[0];
rz(-1.8456521) q[1];
sx q[1];
rz(-1.0546874) q[1];
sx q[1];
rz(-2.4268699) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9172404) q[0];
sx q[0];
rz(-2.9091798) q[0];
sx q[0];
rz(-0.31347855) q[0];
rz(-pi) q[1];
rz(-0.3008414) q[2];
sx q[2];
rz(-1.46035) q[2];
sx q[2];
rz(0.12953239) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.5446075) q[1];
sx q[1];
rz(-1.3860354) q[1];
sx q[1];
rz(-1.9562592) q[1];
x q[2];
rz(-2.7477164) q[3];
sx q[3];
rz(-2.1112313) q[3];
sx q[3];
rz(-0.5462164) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.1456445) q[2];
sx q[2];
rz(-1.6137292) q[2];
sx q[2];
rz(1.145251) q[2];
rz(-0.20259419) q[3];
sx q[3];
rz(-0.069772094) q[3];
sx q[3];
rz(-2.8128459) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8458493) q[0];
sx q[0];
rz(-0.29243094) q[0];
sx q[0];
rz(-2.9737245) q[0];
rz(-0.56495086) q[1];
sx q[1];
rz(-1.047537) q[1];
sx q[1];
rz(1.3289183) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.57083) q[0];
sx q[0];
rz(-2.0432297) q[0];
sx q[0];
rz(-0.60586141) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8064427) q[2];
sx q[2];
rz(-0.14832917) q[2];
sx q[2];
rz(-0.89491612) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.5766474) q[1];
sx q[1];
rz(-2.7912346) q[1];
sx q[1];
rz(2.3957361) q[1];
x q[2];
rz(0.82002016) q[3];
sx q[3];
rz(-2.8201172) q[3];
sx q[3];
rz(0.12382759) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.38587511) q[2];
sx q[2];
rz(-1.2473325) q[2];
sx q[2];
rz(0.34783777) q[2];
rz(2.8268585) q[3];
sx q[3];
rz(-1.0957402) q[3];
sx q[3];
rz(1.1381963) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1494074) q[0];
sx q[0];
rz(-1.5062165) q[0];
sx q[0];
rz(2.1904679) q[0];
rz(0.29257193) q[1];
sx q[1];
rz(-2.2508299) q[1];
sx q[1];
rz(-2.1975885) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4571465) q[0];
sx q[0];
rz(-1.5862042) q[0];
sx q[0];
rz(1.2627361) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3500309) q[2];
sx q[2];
rz(-0.63523173) q[2];
sx q[2];
rz(-0.61295618) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.85691014) q[1];
sx q[1];
rz(-1.3976685) q[1];
sx q[1];
rz(3.1186652) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1066426) q[3];
sx q[3];
rz(-2.6831045) q[3];
sx q[3];
rz(1.5036086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.49147478) q[2];
sx q[2];
rz(-2.398141) q[2];
sx q[2];
rz(-1.4415119) q[2];
rz(-0.024638351) q[3];
sx q[3];
rz(-1.325343) q[3];
sx q[3];
rz(0.98954454) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.80855075) q[0];
sx q[0];
rz(-1.7311743) q[0];
sx q[0];
rz(3.0035875) q[0];
rz(-2.0637312) q[1];
sx q[1];
rz(-1.6777918) q[1];
sx q[1];
rz(0.93625751) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9247583) q[0];
sx q[0];
rz(-0.48473293) q[0];
sx q[0];
rz(0.68883552) q[0];
rz(-pi) q[1];
rz(0.9844043) q[2];
sx q[2];
rz(-0.83293646) q[2];
sx q[2];
rz(-0.10228233) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.030979217) q[1];
sx q[1];
rz(-2.958221) q[1];
sx q[1];
rz(-1.6513837) q[1];
rz(-2.8209723) q[3];
sx q[3];
rz(-0.71133876) q[3];
sx q[3];
rz(2.5742755) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.998698) q[2];
sx q[2];
rz(-1.2724718) q[2];
sx q[2];
rz(1.8227089) q[2];
rz(-3.0190492) q[3];
sx q[3];
rz(-0.69869852) q[3];
sx q[3];
rz(-2.7114649) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1634624) q[0];
sx q[0];
rz(-0.4438816) q[0];
sx q[0];
rz(-0.36460707) q[0];
rz(1.3847903) q[1];
sx q[1];
rz(-2.2099647) q[1];
sx q[1];
rz(0.91420954) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.6624266) q[0];
sx q[0];
rz(-1.437133) q[0];
sx q[0];
rz(-1.1396618) q[0];
rz(-pi) q[1];
rz(2.4070074) q[2];
sx q[2];
rz(-1.5365639) q[2];
sx q[2];
rz(-1.5058668) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.57371512) q[1];
sx q[1];
rz(-2.2063631) q[1];
sx q[1];
rz(-1.741841) q[1];
rz(-pi) q[2];
rz(-0.54085853) q[3];
sx q[3];
rz(-1.4424994) q[3];
sx q[3];
rz(3.1311664) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.8387575) q[2];
sx q[2];
rz(-2.5573533) q[2];
sx q[2];
rz(-1.5128822) q[2];
rz(-0.58879876) q[3];
sx q[3];
rz(-1.9353119) q[3];
sx q[3];
rz(2.485386) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47422472) q[0];
sx q[0];
rz(-0.56741699) q[0];
sx q[0];
rz(-2.8571416) q[0];
rz(-2.4868763) q[1];
sx q[1];
rz(-1.7660716) q[1];
sx q[1];
rz(0.42025748) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66061879) q[0];
sx q[0];
rz(-1.7848564) q[0];
sx q[0];
rz(-1.1461036) q[0];
rz(-pi) q[1];
rz(2.0345694) q[2];
sx q[2];
rz(-2.3056185) q[2];
sx q[2];
rz(-0.42586621) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-3.0033231) q[1];
sx q[1];
rz(-1.4821293) q[1];
sx q[1];
rz(1.6900464) q[1];
x q[2];
rz(0.52892567) q[3];
sx q[3];
rz(-0.63415895) q[3];
sx q[3];
rz(2.8932539) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.1493211) q[2];
sx q[2];
rz(-2.6565266) q[2];
sx q[2];
rz(-2.1743656) q[2];
rz(-2.1238756) q[3];
sx q[3];
rz(-1.8571721) q[3];
sx q[3];
rz(0.73219901) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.16348542) q[0];
sx q[0];
rz(-2.0757984) q[0];
sx q[0];
rz(2.1708873) q[0];
rz(0.12374395) q[1];
sx q[1];
rz(-2.1950304) q[1];
sx q[1];
rz(-1.310941) q[1];
rz(-0.32044784) q[2];
sx q[2];
rz(-1.6978227) q[2];
sx q[2];
rz(2.6200079) q[2];
rz(2.8295201) q[3];
sx q[3];
rz(-1.5371116) q[3];
sx q[3];
rz(-0.1826382) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
