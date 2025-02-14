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
rz(2.702873) q[0];
sx q[0];
rz(-0.60957849) q[0];
sx q[0];
rz(-0.69019812) q[0];
rz(1.26735) q[1];
sx q[1];
rz(-0.48935088) q[1];
sx q[1];
rz(2.7971921) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1699088) q[0];
sx q[0];
rz(-0.2580041) q[0];
sx q[0];
rz(-2.1841316) q[0];
x q[1];
rz(-1.8163258) q[2];
sx q[2];
rz(-1.9787162) q[2];
sx q[2];
rz(-1.7657042) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.5459402) q[1];
sx q[1];
rz(-1.3134688) q[1];
sx q[1];
rz(0.20471065) q[1];
rz(-1.2574168) q[3];
sx q[3];
rz(-2.628577) q[3];
sx q[3];
rz(0.080400217) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.0979746) q[2];
sx q[2];
rz(-2.2100885) q[2];
sx q[2];
rz(-2.0841058) q[2];
rz(2.1438694) q[3];
sx q[3];
rz(-1.3910553) q[3];
sx q[3];
rz(-1.1725496) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69720307) q[0];
sx q[0];
rz(-0.81962219) q[0];
sx q[0];
rz(2.3890553) q[0];
rz(1.2731816) q[1];
sx q[1];
rz(-1.9877142) q[1];
sx q[1];
rz(-2.2862327) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68356791) q[0];
sx q[0];
rz(-0.23580256) q[0];
sx q[0];
rz(-2.0033512) q[0];
rz(-0.94903058) q[2];
sx q[2];
rz(-1.3675642) q[2];
sx q[2];
rz(-1.0946314) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.7400186) q[1];
sx q[1];
rz(-0.42045004) q[1];
sx q[1];
rz(0.40989532) q[1];
rz(-pi) q[2];
rz(1.5454888) q[3];
sx q[3];
rz(-1.3719333) q[3];
sx q[3];
rz(-1.0024827) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.37318834) q[2];
sx q[2];
rz(-1.5311925) q[2];
sx q[2];
rz(-1.9739523) q[2];
rz(0.097955616) q[3];
sx q[3];
rz(-0.28147134) q[3];
sx q[3];
rz(1.9907985) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
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
rz(2.7975467) q[0];
sx q[0];
rz(-2.5757289) q[0];
sx q[0];
rz(1.6746445) q[0];
rz(-3.103718) q[1];
sx q[1];
rz(-1.9297618) q[1];
sx q[1];
rz(1.1143335) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.21143154) q[0];
sx q[0];
rz(-2.7142148) q[0];
sx q[0];
rz(1.2965802) q[0];
rz(-pi) q[1];
rz(-0.55338316) q[2];
sx q[2];
rz(-1.1107365) q[2];
sx q[2];
rz(1.9529238) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.73015651) q[1];
sx q[1];
rz(-1.1056285) q[1];
sx q[1];
rz(-0.48480715) q[1];
rz(-0.30348482) q[3];
sx q[3];
rz(-1.0110098) q[3];
sx q[3];
rz(1.8523077) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.7797839) q[2];
sx q[2];
rz(-2.5461758) q[2];
sx q[2];
rz(-2.6105866) q[2];
rz(0.94046721) q[3];
sx q[3];
rz(-2.2898424) q[3];
sx q[3];
rz(-1.6865591) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4130212) q[0];
sx q[0];
rz(-2.72609) q[0];
sx q[0];
rz(-2.3492133) q[0];
rz(-0.13262311) q[1];
sx q[1];
rz(-2.4515929) q[1];
sx q[1];
rz(-2.2161868) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1559832) q[0];
sx q[0];
rz(-1.4155843) q[0];
sx q[0];
rz(1.7585836) q[0];
rz(-pi) q[1];
rz(1.119198) q[2];
sx q[2];
rz(-1.5182505) q[2];
sx q[2];
rz(0.53942142) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.99605152) q[1];
sx q[1];
rz(-0.62623238) q[1];
sx q[1];
rz(1.9171417) q[1];
x q[2];
rz(2.6230249) q[3];
sx q[3];
rz(-1.858497) q[3];
sx q[3];
rz(0.30687919) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.7856019) q[2];
sx q[2];
rz(-1.2987368) q[2];
sx q[2];
rz(1.7388657) q[2];
rz(-1.8897024) q[3];
sx q[3];
rz(-1.3024878) q[3];
sx q[3];
rz(0.29786626) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4414325) q[0];
sx q[0];
rz(-2.1692363) q[0];
sx q[0];
rz(-2.1527619) q[0];
rz(-1.6255469) q[1];
sx q[1];
rz(-1.6700309) q[1];
sx q[1];
rz(0.54738799) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1158041) q[0];
sx q[0];
rz(-2.0594547) q[0];
sx q[0];
rz(1.5343094) q[0];
rz(-2.9314032) q[2];
sx q[2];
rz(-2.7700305) q[2];
sx q[2];
rz(-2.4331869) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.061364) q[1];
sx q[1];
rz(-1.1669817) q[1];
sx q[1];
rz(-2.3263127) q[1];
x q[2];
rz(-0.27606729) q[3];
sx q[3];
rz(-1.4763586) q[3];
sx q[3];
rz(-0.51968473) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.0872515) q[2];
sx q[2];
rz(-2.6332899) q[2];
sx q[2];
rz(0.15092078) q[2];
rz(-1.6357251) q[3];
sx q[3];
rz(-1.3003636) q[3];
sx q[3];
rz(-1.4668363) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
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
rz(-2.6557789) q[0];
sx q[0];
rz(-1.1968311) q[0];
sx q[0];
rz(-2.6659513) q[0];
rz(-2.7173243) q[1];
sx q[1];
rz(-1.2356267) q[1];
sx q[1];
rz(1.7686527) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7141372) q[0];
sx q[0];
rz(-1.3081453) q[0];
sx q[0];
rz(-1.3930687) q[0];
x q[1];
rz(0.67694725) q[2];
sx q[2];
rz(-0.71093762) q[2];
sx q[2];
rz(2.9589911) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.9109505) q[1];
sx q[1];
rz(-1.9021209) q[1];
sx q[1];
rz(1.001292) q[1];
rz(3.025029) q[3];
sx q[3];
rz(-2.0102083) q[3];
sx q[3];
rz(-1.9842465) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.0515392) q[2];
sx q[2];
rz(-0.97344437) q[2];
sx q[2];
rz(0.092183979) q[2];
rz(-1.9696382) q[3];
sx q[3];
rz(-2.6086174) q[3];
sx q[3];
rz(-2.2251825) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6072657) q[0];
sx q[0];
rz(-0.83616513) q[0];
sx q[0];
rz(2.109206) q[0];
rz(3.0119925) q[1];
sx q[1];
rz(-1.7584636) q[1];
sx q[1];
rz(1.3547156) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1537405) q[0];
sx q[0];
rz(-2.3153785) q[0];
sx q[0];
rz(0.67954845) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9117891) q[2];
sx q[2];
rz(-1.908506) q[2];
sx q[2];
rz(-0.72125283) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.9216283) q[1];
sx q[1];
rz(-2.3458743) q[1];
sx q[1];
rz(-2.398215) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0951429) q[3];
sx q[3];
rz(-1.5009592) q[3];
sx q[3];
rz(-1.2418408) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.9356392) q[2];
sx q[2];
rz(-0.26198584) q[2];
sx q[2];
rz(1.4865173) q[2];
rz(-0.67241159) q[3];
sx q[3];
rz(-2.0629864) q[3];
sx q[3];
rz(0.068517223) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2023778) q[0];
sx q[0];
rz(-0.12513932) q[0];
sx q[0];
rz(2.8676721) q[0];
rz(2.9778453) q[1];
sx q[1];
rz(-1.3122908) q[1];
sx q[1];
rz(0.40072498) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76213259) q[0];
sx q[0];
rz(-0.5665938) q[0];
sx q[0];
rz(-2.6521126) q[0];
x q[1];
rz(1.2873307) q[2];
sx q[2];
rz(-0.41850433) q[2];
sx q[2];
rz(-2.037627) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.5916047) q[1];
sx q[1];
rz(-0.89877909) q[1];
sx q[1];
rz(3.1179291) q[1];
x q[2];
rz(-2.2473281) q[3];
sx q[3];
rz(-2.6163963) q[3];
sx q[3];
rz(2.0390455) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.2694232) q[2];
sx q[2];
rz(-0.17435208) q[2];
sx q[2];
rz(-1.6677469) q[2];
rz(-0.10146865) q[3];
sx q[3];
rz(-1.2227819) q[3];
sx q[3];
rz(-1.7469223) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2673016) q[0];
sx q[0];
rz(-0.75031459) q[0];
sx q[0];
rz(0.0016203298) q[0];
rz(2.5770309) q[1];
sx q[1];
rz(-1.3357342) q[1];
sx q[1];
rz(-0.50216215) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2005794) q[0];
sx q[0];
rz(-2.1920125) q[0];
sx q[0];
rz(2.4401779) q[0];
rz(-1.9829168) q[2];
sx q[2];
rz(-2.0860163) q[2];
sx q[2];
rz(0.56319046) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.50599033) q[1];
sx q[1];
rz(-2.5924304) q[1];
sx q[1];
rz(-0.378774) q[1];
rz(2.929964) q[3];
sx q[3];
rz(-1.0628848) q[3];
sx q[3];
rz(-0.25903364) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(3.0994215) q[2];
sx q[2];
rz(-1.9440938) q[2];
sx q[2];
rz(1.2449123) q[2];
rz(-0.20243195) q[3];
sx q[3];
rz(-1.7499685) q[3];
sx q[3];
rz(2.1421053) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5779483) q[0];
sx q[0];
rz(-0.36643323) q[0];
sx q[0];
rz(-1.7250489) q[0];
rz(-1.1997403) q[1];
sx q[1];
rz(-1.1734633) q[1];
sx q[1];
rz(0.27552584) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4352132) q[0];
sx q[0];
rz(-1.9322898) q[0];
sx q[0];
rz(-0.3630017) q[0];
rz(-1.4813384) q[2];
sx q[2];
rz(-1.2209792) q[2];
sx q[2];
rz(1.1773156) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.72807377) q[1];
sx q[1];
rz(-2.6797543) q[1];
sx q[1];
rz(-1.6175265) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.18003612) q[3];
sx q[3];
rz(-2.535365) q[3];
sx q[3];
rz(-1.1763371) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.766091) q[2];
sx q[2];
rz(-1.6771064) q[2];
sx q[2];
rz(2.7601833) q[2];
rz(2.8219847) q[3];
sx q[3];
rz(-0.8173129) q[3];
sx q[3];
rz(0.6680502) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-1.9850591) q[0];
sx q[0];
rz(-1.1548797) q[0];
sx q[0];
rz(-0.67697939) q[0];
rz(1.001724) q[1];
sx q[1];
rz(-1.4717419) q[1];
sx q[1];
rz(-0.91632661) q[1];
rz(0.15663319) q[2];
sx q[2];
rz(-2.4703783) q[2];
sx q[2];
rz(-3.0117161) q[2];
rz(-1.6826717) q[3];
sx q[3];
rz(-2.2266012) q[3];
sx q[3];
rz(0.83960017) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
