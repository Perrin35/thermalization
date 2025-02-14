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
rz(1.002797) q[1];
sx q[1];
rz(-2.306814) q[1];
sx q[1];
rz(0.53258449) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40520378) q[0];
sx q[0];
rz(-1.4215934) q[0];
sx q[0];
rz(-2.2264541) q[0];
rz(-2.8947484) q[2];
sx q[2];
rz(-2.6181863) q[2];
sx q[2];
rz(0.95062765) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.0610032) q[1];
sx q[1];
rz(-0.71032754) q[1];
sx q[1];
rz(1.8038453) q[1];
x q[2];
rz(-0.37181969) q[3];
sx q[3];
rz(-1.216868) q[3];
sx q[3];
rz(-1.33385) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.3089932) q[2];
sx q[2];
rz(-1.7367881) q[2];
sx q[2];
rz(-0.12283202) q[2];
rz(3.099856) q[3];
sx q[3];
rz(-1.583497) q[3];
sx q[3];
rz(-0.48833716) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
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
rz(-1.1500583) q[0];
sx q[0];
rz(-2.8128862) q[0];
sx q[0];
rz(-1.1154037) q[0];
rz(-0.54488048) q[1];
sx q[1];
rz(-1.3976401) q[1];
sx q[1];
rz(-2.5741408) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3504336) q[0];
sx q[0];
rz(-1.5251524) q[0];
sx q[0];
rz(1.6045536) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3653602) q[2];
sx q[2];
rz(-1.5934332) q[2];
sx q[2];
rz(2.5843888) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.8711155) q[1];
sx q[1];
rz(-1.0024299) q[1];
sx q[1];
rz(1.6898443) q[1];
rz(-2.6065234) q[3];
sx q[3];
rz(-2.5949083) q[3];
sx q[3];
rz(0.018230326) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.6313717) q[2];
sx q[2];
rz(-3.0192182) q[2];
sx q[2];
rz(3.0369634) q[2];
rz(2.6369324) q[3];
sx q[3];
rz(-1.3480836) q[3];
sx q[3];
rz(0.73205718) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0013244) q[0];
sx q[0];
rz(-0.19135419) q[0];
sx q[0];
rz(-1.6562756) q[0];
rz(0.56400076) q[1];
sx q[1];
rz(-2.0878849) q[1];
sx q[1];
rz(-2.3174813) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1230319) q[0];
sx q[0];
rz(-1.0433974) q[0];
sx q[0];
rz(2.9508136) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1698157) q[2];
sx q[2];
rz(-0.87131558) q[2];
sx q[2];
rz(0.60764193) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.6388986) q[1];
sx q[1];
rz(-2.6656409) q[1];
sx q[1];
rz(2.6902728) q[1];
rz(-pi) q[2];
x q[2];
rz(0.070057292) q[3];
sx q[3];
rz(-1.3590711) q[3];
sx q[3];
rz(2.5222561) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(3.04909) q[2];
sx q[2];
rz(-2.8503032) q[2];
sx q[2];
rz(1.5352486) q[2];
rz(0.90318471) q[3];
sx q[3];
rz(-1.3245557) q[3];
sx q[3];
rz(2.7170392) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.20151888) q[0];
sx q[0];
rz(-2.1551082) q[0];
sx q[0];
rz(0.41341138) q[0];
rz(2.9460733) q[1];
sx q[1];
rz(-1.0370516) q[1];
sx q[1];
rz(-1.4541385) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.91823309) q[0];
sx q[0];
rz(-2.4392895) q[0];
sx q[0];
rz(-1.7683517) q[0];
x q[1];
rz(1.3801345) q[2];
sx q[2];
rz(-1.9205928) q[2];
sx q[2];
rz(1.4189394) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.98764963) q[1];
sx q[1];
rz(-2.5763303) q[1];
sx q[1];
rz(2.115519) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3680254) q[3];
sx q[3];
rz(-2.2550809) q[3];
sx q[3];
rz(0.87316421) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.2426408) q[2];
sx q[2];
rz(-0.61598888) q[2];
sx q[2];
rz(1.2658489) q[2];
rz(-2.2856581) q[3];
sx q[3];
rz(-1.3866235) q[3];
sx q[3];
rz(1.9267193) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.6619381) q[0];
sx q[0];
rz(-2.9582773) q[0];
sx q[0];
rz(1.7228175) q[0];
rz(-1.4310369) q[1];
sx q[1];
rz(-1.6173671) q[1];
sx q[1];
rz(-2.4180744) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8531401) q[0];
sx q[0];
rz(-0.61091629) q[0];
sx q[0];
rz(-0.31129654) q[0];
rz(0.75690143) q[2];
sx q[2];
rz(-0.98148966) q[2];
sx q[2];
rz(2.2841931) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-3.0192467) q[1];
sx q[1];
rz(-0.86402383) q[1];
sx q[1];
rz(-3.138423) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.789322) q[3];
sx q[3];
rz(-1.5776424) q[3];
sx q[3];
rz(0.37876836) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.4414703) q[2];
sx q[2];
rz(-2.1258326) q[2];
sx q[2];
rz(0.36386841) q[2];
rz(-1.3077334) q[3];
sx q[3];
rz(-1.6677083) q[3];
sx q[3];
rz(-1.3522805) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7135007) q[0];
sx q[0];
rz(-0.91701737) q[0];
sx q[0];
rz(-2.2416903) q[0];
rz(-2.0948386) q[1];
sx q[1];
rz(-0.37938198) q[1];
sx q[1];
rz(0.55211639) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.715616) q[0];
sx q[0];
rz(-1.5425342) q[0];
sx q[0];
rz(-2.8032599) q[0];
rz(-2.5234473) q[2];
sx q[2];
rz(-2.3380438) q[2];
sx q[2];
rz(-1.8185563) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.99191715) q[1];
sx q[1];
rz(-2.9475355) q[1];
sx q[1];
rz(2.0582576) q[1];
rz(-1.6523727) q[3];
sx q[3];
rz(-0.98203595) q[3];
sx q[3];
rz(-1.0772062) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.3209352) q[2];
sx q[2];
rz(-0.67639095) q[2];
sx q[2];
rz(-2.9429842) q[2];
rz(2.0123539) q[3];
sx q[3];
rz(-1.4720474) q[3];
sx q[3];
rz(0.78013295) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(2.8010913) q[0];
sx q[0];
rz(-1.8393562) q[0];
sx q[0];
rz(0.71402016) q[0];
rz(-1.2227614) q[1];
sx q[1];
rz(-0.87632767) q[1];
sx q[1];
rz(-2.8661935) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3396153) q[0];
sx q[0];
rz(-0.45399779) q[0];
sx q[0];
rz(1.071644) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.95850079) q[2];
sx q[2];
rz(-1.3103518) q[2];
sx q[2];
rz(-2.1672003) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.279567) q[1];
sx q[1];
rz(-2.0525524) q[1];
sx q[1];
rz(2.606009) q[1];
rz(-pi) q[2];
rz(-3.1263859) q[3];
sx q[3];
rz(-2.4162724) q[3];
sx q[3];
rz(1.0909347) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.8347281) q[2];
sx q[2];
rz(-1.4707668) q[2];
sx q[2];
rz(1.6738221) q[2];
rz(1.62014) q[3];
sx q[3];
rz(-2.455267) q[3];
sx q[3];
rz(-0.93792382) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
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
rz(-0.78978157) q[1];
sx q[1];
rz(-2.8826394) q[1];
sx q[1];
rz(0.97377473) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6286929) q[0];
sx q[0];
rz(-1.8159465) q[0];
sx q[0];
rz(3.0184477) q[0];
rz(2.7142286) q[2];
sx q[2];
rz(-1.2334497) q[2];
sx q[2];
rz(-0.34572476) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.5681947) q[1];
sx q[1];
rz(-0.97028448) q[1];
sx q[1];
rz(-0.37094122) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0949368) q[3];
sx q[3];
rz(-1.5175036) q[3];
sx q[3];
rz(2.1691967) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.3070273) q[2];
sx q[2];
rz(-1.4296738) q[2];
sx q[2];
rz(-0.60565051) q[2];
rz(2.0511131) q[3];
sx q[3];
rz(-2.8993789) q[3];
sx q[3];
rz(3.0428913) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7379446) q[0];
sx q[0];
rz(-3.0927959) q[0];
sx q[0];
rz(-2.5700289) q[0];
rz(0.38772186) q[1];
sx q[1];
rz(-2.3523836) q[1];
sx q[1];
rz(-0.8943843) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8799202) q[0];
sx q[0];
rz(-2.399647) q[0];
sx q[0];
rz(0.80961734) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0036976) q[2];
sx q[2];
rz(-1.4722479) q[2];
sx q[2];
rz(-0.71426094) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.91796658) q[1];
sx q[1];
rz(-1.8119748) q[1];
sx q[1];
rz(-1.1583223) q[1];
x q[2];
rz(1.4924516) q[3];
sx q[3];
rz(-1.0634897) q[3];
sx q[3];
rz(-1.1266192) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.9992708) q[2];
sx q[2];
rz(-0.88185507) q[2];
sx q[2];
rz(2.2692661) q[2];
rz(-0.074660389) q[3];
sx q[3];
rz(-1.229769) q[3];
sx q[3];
rz(-0.57620302) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9999303) q[0];
sx q[0];
rz(-2.7003728) q[0];
sx q[0];
rz(-0.73750752) q[0];
rz(2.4420338) q[1];
sx q[1];
rz(-1.0195426) q[1];
sx q[1];
rz(0.84651822) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8171513) q[0];
sx q[0];
rz(-2.2014509) q[0];
sx q[0];
rz(-1.0704499) q[0];
rz(-pi) q[1];
rz(1.6834832) q[2];
sx q[2];
rz(-0.42851617) q[2];
sx q[2];
rz(2.2402415) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.6168931) q[1];
sx q[1];
rz(-1.1438826) q[1];
sx q[1];
rz(0.15097832) q[1];
rz(-pi) q[2];
rz(0.40007253) q[3];
sx q[3];
rz(-1.3988675) q[3];
sx q[3];
rz(-1.7085058) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.49125853) q[2];
sx q[2];
rz(-1.1639872) q[2];
sx q[2];
rz(0.35438806) q[2];
rz(2.3645511) q[3];
sx q[3];
rz(-2.5208426) q[3];
sx q[3];
rz(-2.0508155) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.95018321) q[0];
sx q[0];
rz(-1.8671028) q[0];
sx q[0];
rz(-2.6381459) q[0];
rz(1.7783816) q[1];
sx q[1];
rz(-0.53047219) q[1];
sx q[1];
rz(-0.06906876) q[1];
rz(1.6949881) q[2];
sx q[2];
rz(-2.4249981) q[2];
sx q[2];
rz(-0.80719765) q[2];
rz(1.1152399) q[3];
sx q[3];
rz(-1.8805877) q[3];
sx q[3];
rz(2.7177802) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
