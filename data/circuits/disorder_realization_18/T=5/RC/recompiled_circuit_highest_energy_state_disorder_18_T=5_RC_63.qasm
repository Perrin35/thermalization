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
rz(0.87822479) q[0];
sx q[0];
rz(3.8642519) q[0];
sx q[0];
rz(10.470358) q[0];
rz(0.0027520952) q[1];
sx q[1];
rz(1.6834919) q[1];
sx q[1];
rz(10.145664) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.79851914) q[0];
sx q[0];
rz(-1.4785754) q[0];
sx q[0];
rz(-0.8999805) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3119705) q[2];
sx q[2];
rz(-0.19779655) q[2];
sx q[2];
rz(-0.91413864) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.1926177) q[1];
sx q[1];
rz(-1.1921391) q[1];
sx q[1];
rz(1.8693792) q[1];
rz(-2.3985596) q[3];
sx q[3];
rz(-1.6662906) q[3];
sx q[3];
rz(-2.9955289) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.8568933) q[2];
sx q[2];
rz(-1.928494) q[2];
sx q[2];
rz(2.6846679) q[2];
rz(2.4871155) q[3];
sx q[3];
rz(-0.5623397) q[3];
sx q[3];
rz(0.42403179) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65861312) q[0];
sx q[0];
rz(-0.35372666) q[0];
sx q[0];
rz(-0.60802996) q[0];
rz(1.1990625) q[1];
sx q[1];
rz(-0.46747318) q[1];
sx q[1];
rz(-2.3466568) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2450144) q[0];
sx q[0];
rz(-0.45077401) q[0];
sx q[0];
rz(0.71433492) q[0];
x q[1];
rz(-2.8744187) q[2];
sx q[2];
rz(-1.5985711) q[2];
sx q[2];
rz(-2.6849314) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.39917377) q[1];
sx q[1];
rz(-1.292242) q[1];
sx q[1];
rz(2.2225077) q[1];
rz(2.7057418) q[3];
sx q[3];
rz(-2.4727949) q[3];
sx q[3];
rz(-2.8613402) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.6752211) q[2];
sx q[2];
rz(-2.0872865) q[2];
sx q[2];
rz(1.4804473) q[2];
rz(2.772707) q[3];
sx q[3];
rz(-1.2494913) q[3];
sx q[3];
rz(2.4108385) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.22496255) q[0];
sx q[0];
rz(-2.1376305) q[0];
sx q[0];
rz(0.68600255) q[0];
rz(-1.2432159) q[1];
sx q[1];
rz(-2.9153283) q[1];
sx q[1];
rz(-0.6023947) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7083654) q[0];
sx q[0];
rz(-2.4658563) q[0];
sx q[0];
rz(2.2883718) q[0];
rz(-pi) q[1];
x q[1];
rz(3.1092293) q[2];
sx q[2];
rz(-1.5302858) q[2];
sx q[2];
rz(1.4537741) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.2418306) q[1];
sx q[1];
rz(-1.6679523) q[1];
sx q[1];
rz(2.1400356) q[1];
rz(1.1021578) q[3];
sx q[3];
rz(-1.2775565) q[3];
sx q[3];
rz(1.6654429) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.7556222) q[2];
sx q[2];
rz(-1.3034416) q[2];
sx q[2];
rz(-1.2516359) q[2];
rz(-2.5899467) q[3];
sx q[3];
rz(-0.14277221) q[3];
sx q[3];
rz(-0.84757203) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.30152339) q[0];
sx q[0];
rz(-1.2395549) q[0];
sx q[0];
rz(3.0384592) q[0];
rz(-0.81941191) q[1];
sx q[1];
rz(-2.2401919) q[1];
sx q[1];
rz(-0.55319667) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5846636) q[0];
sx q[0];
rz(-1.9346666) q[0];
sx q[0];
rz(-2.8332024) q[0];
rz(-pi) q[1];
rz(-2.791094) q[2];
sx q[2];
rz(-1.0523049) q[2];
sx q[2];
rz(1.1120591) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.4955289) q[1];
sx q[1];
rz(-0.49139443) q[1];
sx q[1];
rz(1.8629406) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0023213) q[3];
sx q[3];
rz(-1.6628886) q[3];
sx q[3];
rz(1.0669277) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.72054046) q[2];
sx q[2];
rz(-0.74137509) q[2];
sx q[2];
rz(0.53483024) q[2];
rz(-2.363291) q[3];
sx q[3];
rz(-0.72746593) q[3];
sx q[3];
rz(0.21540575) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68742037) q[0];
sx q[0];
rz(-1.9947616) q[0];
sx q[0];
rz(-2.4046894) q[0];
rz(-3.0829644) q[1];
sx q[1];
rz(-0.77719378) q[1];
sx q[1];
rz(-2.5545205) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.042087) q[0];
sx q[0];
rz(-1.6027959) q[0];
sx q[0];
rz(2.2048897) q[0];
x q[1];
rz(-0.063063085) q[2];
sx q[2];
rz(-1.3418256) q[2];
sx q[2];
rz(-1.582043) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.88564155) q[1];
sx q[1];
rz(-1.6033065) q[1];
sx q[1];
rz(-1.2684275) q[1];
rz(-pi) q[2];
x q[2];
rz(0.30549099) q[3];
sx q[3];
rz(-1.8754212) q[3];
sx q[3];
rz(-2.6412499) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.8495142) q[2];
sx q[2];
rz(-2.2167315) q[2];
sx q[2];
rz(-0.64797956) q[2];
rz(0.61602151) q[3];
sx q[3];
rz(-1.5024622) q[3];
sx q[3];
rz(-3.093921) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7306526) q[0];
sx q[0];
rz(-2.3766282) q[0];
sx q[0];
rz(-2.1867645) q[0];
rz(-3.094589) q[1];
sx q[1];
rz(-0.67263293) q[1];
sx q[1];
rz(-2.1254553) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.34441062) q[0];
sx q[0];
rz(-1.1007995) q[0];
sx q[0];
rz(-0.60352602) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.418265) q[2];
sx q[2];
rz(-1.1473341) q[2];
sx q[2];
rz(-2.9996769) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.816427) q[1];
sx q[1];
rz(-1.0652958) q[1];
sx q[1];
rz(-1.5542404) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.627907) q[3];
sx q[3];
rz(-0.94739238) q[3];
sx q[3];
rz(0.36950612) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.59812349) q[2];
sx q[2];
rz(-1.8757966) q[2];
sx q[2];
rz(-0.050696105) q[2];
rz(-3.0905881) q[3];
sx q[3];
rz(-1.8404605) q[3];
sx q[3];
rz(0.72371975) q[3];
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
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8766668) q[0];
sx q[0];
rz(-2.6267138) q[0];
sx q[0];
rz(1.6616954) q[0];
rz(-2.0199846) q[1];
sx q[1];
rz(-1.2622967) q[1];
sx q[1];
rz(-2.6720537) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9037989) q[0];
sx q[0];
rz(-1.5241166) q[0];
sx q[0];
rz(1.4633798) q[0];
x q[1];
rz(1.1781349) q[2];
sx q[2];
rz(-1.715065) q[2];
sx q[2];
rz(-0.96111125) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.28227678) q[1];
sx q[1];
rz(-2.6398602) q[1];
sx q[1];
rz(-1.938107) q[1];
x q[2];
rz(-2.363405) q[3];
sx q[3];
rz(-2.0594849) q[3];
sx q[3];
rz(-2.5358939) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.40662128) q[2];
sx q[2];
rz(-1.5903641) q[2];
sx q[2];
rz(-2.0242019) q[2];
rz(-1.0920852) q[3];
sx q[3];
rz(-1.4039682) q[3];
sx q[3];
rz(2.7058069) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0206873) q[0];
sx q[0];
rz(-1.1130604) q[0];
sx q[0];
rz(-3.0764965) q[0];
rz(-2.8096325) q[1];
sx q[1];
rz(-1.8354225) q[1];
sx q[1];
rz(-1.8665705) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41439636) q[0];
sx q[0];
rz(-2.9119618) q[0];
sx q[0];
rz(-1.5288607) q[0];
rz(-1.128007) q[2];
sx q[2];
rz(-1.4472359) q[2];
sx q[2];
rz(2.3905001) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.4091585) q[1];
sx q[1];
rz(-2.2526764) q[1];
sx q[1];
rz(-1.5049008) q[1];
rz(-1.3887819) q[3];
sx q[3];
rz(-0.776178) q[3];
sx q[3];
rz(2.9266199) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.28978213) q[2];
sx q[2];
rz(-1.449457) q[2];
sx q[2];
rz(2.1688478) q[2];
rz(-0.70074493) q[3];
sx q[3];
rz(-2.2795129) q[3];
sx q[3];
rz(-2.3166482) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7026611) q[0];
sx q[0];
rz(-2.7362566) q[0];
sx q[0];
rz(2.7324556) q[0];
rz(1.1459972) q[1];
sx q[1];
rz(-1.0640249) q[1];
sx q[1];
rz(-1.8786059) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.43009297) q[0];
sx q[0];
rz(-1.5706675) q[0];
sx q[0];
rz(1.5789728) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8645889) q[2];
sx q[2];
rz(-1.190349) q[2];
sx q[2];
rz(1.6858619) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.2330875) q[1];
sx q[1];
rz(-0.89492866) q[1];
sx q[1];
rz(0.29260537) q[1];
rz(-pi) q[2];
rz(1.658012) q[3];
sx q[3];
rz(-2.1391807) q[3];
sx q[3];
rz(1.6001112) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.24141773) q[2];
sx q[2];
rz(-1.9200385) q[2];
sx q[2];
rz(-2.9202666) q[2];
rz(-2.6545702) q[3];
sx q[3];
rz(-2.2723787) q[3];
sx q[3];
rz(-1.2041436) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0601592) q[0];
sx q[0];
rz(-0.24516036) q[0];
sx q[0];
rz(1.2231476) q[0];
rz(3.0606015) q[1];
sx q[1];
rz(-1.6209737) q[1];
sx q[1];
rz(1.5222668) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8232097) q[0];
sx q[0];
rz(-1.902474) q[0];
sx q[0];
rz(1.7451502) q[0];
rz(-0.98098251) q[2];
sx q[2];
rz(-2.207617) q[2];
sx q[2];
rz(2.3825519) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-3.1159525) q[1];
sx q[1];
rz(-1.2668477) q[1];
sx q[1];
rz(2.9362684) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.31701458) q[3];
sx q[3];
rz(-0.51643378) q[3];
sx q[3];
rz(1.3610193) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(3.0167375) q[2];
sx q[2];
rz(-1.7502461) q[2];
sx q[2];
rz(2.4511231) q[2];
rz(0.27103439) q[3];
sx q[3];
rz(-0.46509585) q[3];
sx q[3];
rz(-1.5439699) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6799714) q[0];
sx q[0];
rz(-0.93628708) q[0];
sx q[0];
rz(-0.17967889) q[0];
rz(-0.24196729) q[1];
sx q[1];
rz(-0.89121834) q[1];
sx q[1];
rz(-1.2527087) q[1];
rz(-0.42226462) q[2];
sx q[2];
rz(-1.3801856) q[2];
sx q[2];
rz(-1.2587067) q[2];
rz(-0.88884066) q[3];
sx q[3];
rz(-1.325229) q[3];
sx q[3];
rz(-1.6390507) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
