OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.2607245) q[0];
sx q[0];
rz(-0.27624929) q[0];
sx q[0];
rz(0.84375381) q[0];
rz(1.002797) q[1];
sx q[1];
rz(3.9763713) q[1];
sx q[1];
rz(3.6741771) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40520378) q[0];
sx q[0];
rz(-1.7199992) q[0];
sx q[0];
rz(0.91513855) q[0];
rz(-pi) q[1];
rz(-0.5101846) q[2];
sx q[2];
rz(-1.4483588) q[2];
sx q[2];
rz(-0.40529006) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.3644173) q[1];
sx q[1];
rz(-2.2580986) q[1];
sx q[1];
rz(0.19608325) q[1];
x q[2];
rz(1.9483564) q[3];
sx q[3];
rz(-1.9185563) q[3];
sx q[3];
rz(2.7703066) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.3089932) q[2];
sx q[2];
rz(-1.4048046) q[2];
sx q[2];
rz(-0.12283202) q[2];
rz(0.041736688) q[3];
sx q[3];
rz(-1.5580956) q[3];
sx q[3];
rz(2.6532555) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
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
rz(-1.1500583) q[0];
sx q[0];
rz(-2.8128862) q[0];
sx q[0];
rz(2.026189) q[0];
rz(-0.54488048) q[1];
sx q[1];
rz(-1.7439525) q[1];
sx q[1];
rz(2.5741408) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.79115907) q[0];
sx q[0];
rz(-1.6164403) q[0];
sx q[0];
rz(1.6045536) q[0];
x q[1];
rz(0.77623244) q[2];
sx q[2];
rz(-1.5934332) q[2];
sx q[2];
rz(-0.55720383) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.9055674) q[1];
sx q[1];
rz(-1.4705338) q[1];
sx q[1];
rz(-2.5699993) q[1];
rz(-0.48226815) q[3];
sx q[3];
rz(-1.3025138) q[3];
sx q[3];
rz(2.0577459) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.6313717) q[2];
sx q[2];
rz(-0.12237445) q[2];
sx q[2];
rz(0.1046293) q[2];
rz(2.6369324) q[3];
sx q[3];
rz(-1.3480836) q[3];
sx q[3];
rz(-2.4095355) q[3];
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
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1402682) q[0];
sx q[0];
rz(-0.19135419) q[0];
sx q[0];
rz(-1.485317) q[0];
rz(2.5775919) q[1];
sx q[1];
rz(-1.0537078) q[1];
sx q[1];
rz(-2.3174813) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4553412) q[0];
sx q[0];
rz(-1.7353962) q[0];
sx q[0];
rz(2.1061939) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5512087) q[2];
sx q[2];
rz(-2.2548642) q[2];
sx q[2];
rz(1.7192529) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.1396795) q[1];
sx q[1];
rz(-1.1458108) q[1];
sx q[1];
rz(-1.3496467) q[1];
x q[2];
rz(-3.0715354) q[3];
sx q[3];
rz(-1.7825216) q[3];
sx q[3];
rz(-2.5222561) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.092502681) q[2];
sx q[2];
rz(-0.29128942) q[2];
sx q[2];
rz(1.606344) q[2];
rz(2.2384079) q[3];
sx q[3];
rz(-1.3245557) q[3];
sx q[3];
rz(-2.7170392) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20151888) q[0];
sx q[0];
rz(-2.1551082) q[0];
sx q[0];
rz(-0.41341138) q[0];
rz(-0.19551936) q[1];
sx q[1];
rz(-2.1045411) q[1];
sx q[1];
rz(1.4541385) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66178807) q[0];
sx q[0];
rz(-2.256752) q[0];
sx q[0];
rz(-2.9770026) q[0];
x q[1];
rz(0.35572534) q[2];
sx q[2];
rz(-1.7497831) q[2];
sx q[2];
rz(2.9236874) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.98764963) q[1];
sx q[1];
rz(-2.5763303) q[1];
sx q[1];
rz(1.0260737) q[1];
rz(-pi) q[2];
rz(0.69444434) q[3];
sx q[3];
rz(-1.4141091) q[3];
sx q[3];
rz(0.82686916) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.2426408) q[2];
sx q[2];
rz(-0.61598888) q[2];
sx q[2];
rz(-1.8757437) q[2];
rz(-0.85593456) q[3];
sx q[3];
rz(-1.3866235) q[3];
sx q[3];
rz(1.2148733) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4796546) q[0];
sx q[0];
rz(-2.9582773) q[0];
sx q[0];
rz(-1.7228175) q[0];
rz(1.7105557) q[1];
sx q[1];
rz(-1.6173671) q[1];
sx q[1];
rz(-2.4180744) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2884525) q[0];
sx q[0];
rz(-2.5306764) q[0];
sx q[0];
rz(-0.31129654) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3695643) q[2];
sx q[2];
rz(-0.92185045) q[2];
sx q[2];
rz(-1.2452425) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.4505087) q[1];
sx q[1];
rz(-1.5732068) q[1];
sx q[1];
rz(-0.86402135) q[1];
x q[2];
rz(-2.789322) q[3];
sx q[3];
rz(-1.5639502) q[3];
sx q[3];
rz(-0.37876836) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.7001223) q[2];
sx q[2];
rz(-2.1258326) q[2];
sx q[2];
rz(-2.7777242) q[2];
rz(1.3077334) q[3];
sx q[3];
rz(-1.6677083) q[3];
sx q[3];
rz(-1.7893121) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7135007) q[0];
sx q[0];
rz(-2.2245753) q[0];
sx q[0];
rz(-0.8999024) q[0];
rz(1.046754) q[1];
sx q[1];
rz(-0.37938198) q[1];
sx q[1];
rz(-2.5894763) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9166342) q[0];
sx q[0];
rz(-0.33946589) q[0];
sx q[0];
rz(3.0566264) q[0];
rz(-pi) q[1];
x q[1];
rz(0.61814536) q[2];
sx q[2];
rz(-0.80354881) q[2];
sx q[2];
rz(1.8185563) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.49657288) q[1];
sx q[1];
rz(-1.7420118) q[1];
sx q[1];
rz(3.0497996) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5512929) q[3];
sx q[3];
rz(-1.5029782) q[3];
sx q[3];
rz(-2.6026311) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.8206574) q[2];
sx q[2];
rz(-2.4652017) q[2];
sx q[2];
rz(-0.19860849) q[2];
rz(-1.1292388) q[3];
sx q[3];
rz(-1.4720474) q[3];
sx q[3];
rz(0.78013295) q[3];
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
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8010913) q[0];
sx q[0];
rz(-1.3022364) q[0];
sx q[0];
rz(0.71402016) q[0];
rz(-1.2227614) q[1];
sx q[1];
rz(-0.87632767) q[1];
sx q[1];
rz(0.27539918) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3396153) q[0];
sx q[0];
rz(-2.6875949) q[0];
sx q[0];
rz(2.0699487) q[0];
x q[1];
rz(2.0049663) q[2];
sx q[2];
rz(-2.482802) q[2];
sx q[2];
rz(0.24497797) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(3.1186953) q[1];
sx q[1];
rz(-2.0400908) q[1];
sx q[1];
rz(1.0245299) q[1];
rz(0.72526284) q[3];
sx q[3];
rz(-1.5607087) q[3];
sx q[3];
rz(0.46848224) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.3068646) q[2];
sx q[2];
rz(-1.6708259) q[2];
sx q[2];
rz(1.4677706) q[2];
rz(-1.62014) q[3];
sx q[3];
rz(-2.455267) q[3];
sx q[3];
rz(-2.2036688) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(-2.9891605) q[0];
sx q[0];
rz(-0.91355649) q[0];
sx q[0];
rz(-2.1754225) q[0];
rz(2.3518111) q[1];
sx q[1];
rz(-2.8826394) q[1];
sx q[1];
rz(0.97377473) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1570996) q[0];
sx q[0];
rz(-0.2737845) q[0];
sx q[0];
rz(-2.0272001) q[0];
rz(-0.42736407) q[2];
sx q[2];
rz(-1.9081429) q[2];
sx q[2];
rz(0.34572476) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.3605347) q[1];
sx q[1];
rz(-1.2670749) q[1];
sx q[1];
rz(-2.2045434) q[1];
rz(-pi) q[2];
rz(3.0816541) q[3];
sx q[3];
rz(-2.0459243) q[3];
sx q[3];
rz(-2.5157473) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.3070273) q[2];
sx q[2];
rz(-1.4296738) q[2];
sx q[2];
rz(2.5359421) q[2];
rz(-1.0904795) q[3];
sx q[3];
rz(-2.8993789) q[3];
sx q[3];
rz(-0.098701326) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.403648) q[0];
sx q[0];
rz(-0.048796766) q[0];
sx q[0];
rz(2.5700289) q[0];
rz(2.7538708) q[1];
sx q[1];
rz(-2.3523836) q[1];
sx q[1];
rz(-2.2472084) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8799202) q[0];
sx q[0];
rz(-2.399647) q[0];
sx q[0];
rz(0.80961734) q[0];
rz(3.0249075) q[2];
sx q[2];
rz(-1.0067847) q[2];
sx q[2];
rz(2.2224768) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.9891178) q[1];
sx q[1];
rz(-0.47429171) q[1];
sx q[1];
rz(-1.0204587) q[1];
rz(-pi) q[2];
rz(-2.6329805) q[3];
sx q[3];
rz(-1.6392574) q[3];
sx q[3];
rz(-0.48229846) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.1423219) q[2];
sx q[2];
rz(-0.88185507) q[2];
sx q[2];
rz(0.87232653) q[2];
rz(-3.0669323) q[3];
sx q[3];
rz(-1.229769) q[3];
sx q[3];
rz(-2.5653896) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1416624) q[0];
sx q[0];
rz(-0.44121989) q[0];
sx q[0];
rz(0.73750752) q[0];
rz(-0.69955889) q[1];
sx q[1];
rz(-2.12205) q[1];
sx q[1];
rz(2.2950744) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0695037) q[0];
sx q[0];
rz(-0.78323019) q[0];
sx q[0];
rz(0.58132078) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9969127) q[2];
sx q[2];
rz(-1.5240545) q[2];
sx q[2];
rz(0.5668723) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(3.032585) q[1];
sx q[1];
rz(-1.7081339) q[1];
sx q[1];
rz(2.0020383) q[1];
rz(-pi) q[2];
x q[2];
rz(1.75714) q[3];
sx q[3];
rz(-1.1769503) q[3];
sx q[3];
rz(0.20993103) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.49125853) q[2];
sx q[2];
rz(-1.9776055) q[2];
sx q[2];
rz(2.7872046) q[2];
rz(-2.3645511) q[3];
sx q[3];
rz(-2.5208426) q[3];
sx q[3];
rz(2.0508155) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1914094) q[0];
sx q[0];
rz(-1.8671028) q[0];
sx q[0];
rz(-2.6381459) q[0];
rz(-1.3632111) q[1];
sx q[1];
rz(-0.53047219) q[1];
sx q[1];
rz(-0.06906876) q[1];
rz(-2.2835636) q[2];
sx q[2];
rz(-1.4893441) q[2];
sx q[2];
rz(-2.2841397) q[2];
rz(2.0263528) q[3];
sx q[3];
rz(-1.2610049) q[3];
sx q[3];
rz(-0.42381248) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
