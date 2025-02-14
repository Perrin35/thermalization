OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.8808682) q[0];
sx q[0];
rz(-2.8653434) q[0];
sx q[0];
rz(-0.84375381) q[0];
rz(-2.1387956) q[1];
sx q[1];
rz(-0.83477867) q[1];
sx q[1];
rz(-0.53258449) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7363889) q[0];
sx q[0];
rz(-1.4215934) q[0];
sx q[0];
rz(-0.91513855) q[0];
rz(-pi) q[1];
rz(2.8947484) q[2];
sx q[2];
rz(-0.52340639) q[2];
sx q[2];
rz(-2.190965) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.0805894) q[1];
sx q[1];
rz(-2.4312651) q[1];
sx q[1];
rz(-1.3377473) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1932362) q[3];
sx q[3];
rz(-1.2230363) q[3];
sx q[3];
rz(-0.37128604) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.8325995) q[2];
sx q[2];
rz(-1.7367881) q[2];
sx q[2];
rz(-3.0187606) q[2];
rz(-3.099856) q[3];
sx q[3];
rz(-1.583497) q[3];
sx q[3];
rz(0.48833716) q[3];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9915344) q[0];
sx q[0];
rz(-0.3287065) q[0];
sx q[0];
rz(-2.026189) q[0];
rz(-0.54488048) q[1];
sx q[1];
rz(-1.7439525) q[1];
sx q[1];
rz(-0.56745183) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3504336) q[0];
sx q[0];
rz(-1.6164403) q[0];
sx q[0];
rz(1.537039) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.602515) q[2];
sx q[2];
rz(-0.79481541) q[2];
sx q[2];
rz(-1.0358126) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.9055674) q[1];
sx q[1];
rz(-1.6710588) q[1];
sx q[1];
rz(0.57159337) q[1];
x q[2];
rz(-2.6065234) q[3];
sx q[3];
rz(-2.5949083) q[3];
sx q[3];
rz(0.018230326) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.6313717) q[2];
sx q[2];
rz(-0.12237445) q[2];
sx q[2];
rz(0.1046293) q[2];
rz(-0.50466022) q[3];
sx q[3];
rz(-1.3480836) q[3];
sx q[3];
rz(-2.4095355) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1402682) q[0];
sx q[0];
rz(-0.19135419) q[0];
sx q[0];
rz(-1.485317) q[0];
rz(0.56400076) q[1];
sx q[1];
rz(-2.0878849) q[1];
sx q[1];
rz(-2.3174813) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1230319) q[0];
sx q[0];
rz(-2.0981952) q[0];
sx q[0];
rz(2.9508136) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3468909) q[2];
sx q[2];
rz(-2.0168755) q[2];
sx q[2];
rz(-0.54903713) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.0019131) q[1];
sx q[1];
rz(-1.9957819) q[1];
sx q[1];
rz(-1.3496467) q[1];
rz(-1.885627) q[3];
sx q[3];
rz(-2.9187429) q[3];
sx q[3];
rz(0.94161445) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.092502681) q[2];
sx q[2];
rz(-2.8503032) q[2];
sx q[2];
rz(1.5352486) q[2];
rz(2.2384079) q[3];
sx q[3];
rz(-1.817037) q[3];
sx q[3];
rz(-0.42455348) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9400738) q[0];
sx q[0];
rz(-2.1551082) q[0];
sx q[0];
rz(2.7281813) q[0];
rz(2.9460733) q[1];
sx q[1];
rz(-1.0370516) q[1];
sx q[1];
rz(1.6874541) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2233596) q[0];
sx q[0];
rz(-2.4392895) q[0];
sx q[0];
rz(1.7683517) q[0];
x q[1];
rz(-1.7614582) q[2];
sx q[2];
rz(-1.2209998) q[2];
sx q[2];
rz(1.7226532) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.5315872) q[1];
sx q[1];
rz(-2.0466699) q[1];
sx q[1];
rz(0.31756084) q[1];
rz(1.3680254) q[3];
sx q[3];
rz(-0.88651171) q[3];
sx q[3];
rz(2.2684284) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.89895189) q[2];
sx q[2];
rz(-2.5256038) q[2];
sx q[2];
rz(-1.8757437) q[2];
rz(2.2856581) q[3];
sx q[3];
rz(-1.7549691) q[3];
sx q[3];
rz(-1.2148733) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.6619381) q[0];
sx q[0];
rz(-0.18331535) q[0];
sx q[0];
rz(1.7228175) q[0];
rz(-1.4310369) q[1];
sx q[1];
rz(-1.6173671) q[1];
sx q[1];
rz(0.72351825) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.60155) q[0];
sx q[0];
rz(-1.7474084) q[0];
sx q[0];
rz(-2.5536182) q[0];
rz(2.3143598) q[2];
sx q[2];
rz(-0.9632573) q[2];
sx q[2];
rz(2.9116258) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.11746487) q[1];
sx q[1];
rz(-2.4348143) q[1];
sx q[1];
rz(-1.5745082) q[1];
x q[2];
rz(0.35227065) q[3];
sx q[3];
rz(-1.5639502) q[3];
sx q[3];
rz(2.7628243) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.4414703) q[2];
sx q[2];
rz(-1.0157601) q[2];
sx q[2];
rz(2.7777242) q[2];
rz(-1.8338592) q[3];
sx q[3];
rz(-1.6677083) q[3];
sx q[3];
rz(-1.7893121) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7135007) q[0];
sx q[0];
rz(-2.2245753) q[0];
sx q[0];
rz(2.2416903) q[0];
rz(-2.0948386) q[1];
sx q[1];
rz(-2.7622107) q[1];
sx q[1];
rz(2.5894763) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.22495843) q[0];
sx q[0];
rz(-2.8021268) q[0];
sx q[0];
rz(3.0566264) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4399567) q[2];
sx q[2];
rz(-2.0011099) q[2];
sx q[2];
rz(-0.21077354) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.0585409) q[1];
sx q[1];
rz(-1.4803491) q[1];
sx q[1];
rz(-1.3988711) q[1];
rz(-1.6523727) q[3];
sx q[3];
rz(-2.1595567) q[3];
sx q[3];
rz(-2.0643864) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.8206574) q[2];
sx q[2];
rz(-0.67639095) q[2];
sx q[2];
rz(-2.9429842) q[2];
rz(-2.0123539) q[3];
sx q[3];
rz(-1.6695453) q[3];
sx q[3];
rz(0.78013295) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3405014) q[0];
sx q[0];
rz(-1.3022364) q[0];
sx q[0];
rz(0.71402016) q[0];
rz(1.9188312) q[1];
sx q[1];
rz(-0.87632767) q[1];
sx q[1];
rz(0.27539918) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8283683) q[0];
sx q[0];
rz(-1.3592915) q[0];
sx q[0];
rz(-1.9755892) q[0];
rz(-pi) q[1];
rz(-1.1366264) q[2];
sx q[2];
rz(-2.482802) q[2];
sx q[2];
rz(-2.8966147) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.95412817) q[1];
sx q[1];
rz(-2.4373332) q[1];
sx q[1];
rz(-0.7974979) q[1];
rz(-0.015206788) q[3];
sx q[3];
rz(-2.4162724) q[3];
sx q[3];
rz(2.050658) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.3068646) q[2];
sx q[2];
rz(-1.6708259) q[2];
sx q[2];
rz(-1.6738221) q[2];
rz(-1.62014) q[3];
sx q[3];
rz(-0.68632564) q[3];
sx q[3];
rz(2.2036688) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9891605) q[0];
sx q[0];
rz(-2.2280362) q[0];
sx q[0];
rz(0.96617019) q[0];
rz(-2.3518111) q[1];
sx q[1];
rz(-2.8826394) q[1];
sx q[1];
rz(2.1678179) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0278661) q[0];
sx q[0];
rz(-1.4513512) q[0];
sx q[0];
rz(-1.3238504) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.70232161) q[2];
sx q[2];
rz(-2.6036545) q[2];
sx q[2];
rz(2.5449716) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.3605347) q[1];
sx q[1];
rz(-1.8745178) q[1];
sx q[1];
rz(-0.93704929) q[1];
rz(-pi) q[2];
rz(2.0466558) q[3];
sx q[3];
rz(-1.5175036) q[3];
sx q[3];
rz(-2.1691967) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.3070273) q[2];
sx q[2];
rz(-1.4296738) q[2];
sx q[2];
rz(-2.5359421) q[2];
rz(-2.0511131) q[3];
sx q[3];
rz(-0.24221371) q[3];
sx q[3];
rz(-0.098701326) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7379446) q[0];
sx q[0];
rz(-0.048796766) q[0];
sx q[0];
rz(-0.57156372) q[0];
rz(0.38772186) q[1];
sx q[1];
rz(-2.3523836) q[1];
sx q[1];
rz(-0.8943843) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26167241) q[0];
sx q[0];
rz(-0.74194569) q[0];
sx q[0];
rz(-2.3319753) q[0];
rz(-pi) q[1];
rz(1.7528083) q[2];
sx q[2];
rz(-0.57467206) q[2];
sx q[2];
rz(2.4383309) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.2236261) q[1];
sx q[1];
rz(-1.8119748) q[1];
sx q[1];
rz(-1.1583223) q[1];
rz(-pi) q[2];
rz(0.13988755) q[3];
sx q[3];
rz(-2.6287946) q[3];
sx q[3];
rz(0.96641216) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.1423219) q[2];
sx q[2];
rz(-2.2597376) q[2];
sx q[2];
rz(-0.87232653) q[2];
rz(0.074660389) q[3];
sx q[3];
rz(-1.9118237) q[3];
sx q[3];
rz(2.5653896) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1416624) q[0];
sx q[0];
rz(-0.44121989) q[0];
sx q[0];
rz(-0.73750752) q[0];
rz(0.69955889) q[1];
sx q[1];
rz(-2.12205) q[1];
sx q[1];
rz(0.84651822) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32444137) q[0];
sx q[0];
rz(-2.2014509) q[0];
sx q[0];
rz(2.0711428) q[0];
rz(-1.4581095) q[2];
sx q[2];
rz(-0.42851617) q[2];
sx q[2];
rz(-0.90135114) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(3.032585) q[1];
sx q[1];
rz(-1.7081339) q[1];
sx q[1];
rz(-1.1395543) q[1];
rz(0.40007253) q[3];
sx q[3];
rz(-1.7427252) q[3];
sx q[3];
rz(-1.4330869) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.49125853) q[2];
sx q[2];
rz(-1.1639872) q[2];
sx q[2];
rz(2.7872046) q[2];
rz(0.77704159) q[3];
sx q[3];
rz(-2.5208426) q[3];
sx q[3];
rz(2.0508155) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.95018321) q[0];
sx q[0];
rz(-1.2744899) q[0];
sx q[0];
rz(0.50344678) q[0];
rz(-1.7783816) q[1];
sx q[1];
rz(-2.6111205) q[1];
sx q[1];
rz(3.0725239) q[1];
rz(0.10748482) q[2];
sx q[2];
rz(-0.86089118) q[2];
sx q[2];
rz(-0.64313342) q[2];
rz(-0.94181413) q[3];
sx q[3];
rz(-0.54472803) q[3];
sx q[3];
rz(1.7036078) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
