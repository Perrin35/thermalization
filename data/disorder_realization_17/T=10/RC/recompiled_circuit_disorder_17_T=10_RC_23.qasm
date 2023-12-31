OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.9085812) q[0];
sx q[0];
rz(-1.9549978) q[0];
sx q[0];
rz(-1.1458122) q[0];
rz(0.97283483) q[1];
sx q[1];
rz(-1.4714779) q[1];
sx q[1];
rz(0.29247984) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.16380331) q[0];
sx q[0];
rz(-1.6262494) q[0];
sx q[0];
rz(-1.9652912) q[0];
x q[1];
rz(0.77941676) q[2];
sx q[2];
rz(-1.5343108) q[2];
sx q[2];
rz(1.2365637) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.7293207) q[1];
sx q[1];
rz(-2.1643157) q[1];
sx q[1];
rz(-3.0862942) q[1];
rz(-pi) q[2];
rz(1.0969767) q[3];
sx q[3];
rz(-1.572527) q[3];
sx q[3];
rz(-0.8918106) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.73074377) q[2];
sx q[2];
rz(-1.0935254) q[2];
sx q[2];
rz(0.4494108) q[2];
rz(-2.4959026) q[3];
sx q[3];
rz(-0.68104762) q[3];
sx q[3];
rz(1.8507563) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(2.9280424) q[0];
sx q[0];
rz(-0.88590652) q[0];
sx q[0];
rz(2.399562) q[0];
rz(1.6702601) q[1];
sx q[1];
rz(-2.4447618) q[1];
sx q[1];
rz(-1.3630294) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.33107685) q[0];
sx q[0];
rz(-1.9609945) q[0];
sx q[0];
rz(-0.66731989) q[0];
rz(0.40740168) q[2];
sx q[2];
rz(-0.50440895) q[2];
sx q[2];
rz(-2.2638869) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.14582536) q[1];
sx q[1];
rz(-1.9103423) q[1];
sx q[1];
rz(1.550436) q[1];
rz(-pi) q[2];
x q[2];
rz(0.19248776) q[3];
sx q[3];
rz(-1.3093595) q[3];
sx q[3];
rz(0.9769494) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.30119511) q[2];
sx q[2];
rz(-2.6448554) q[2];
sx q[2];
rz(1.2724686) q[2];
rz(2.8043591) q[3];
sx q[3];
rz(-1.1703346) q[3];
sx q[3];
rz(-0.01005323) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0740046) q[0];
sx q[0];
rz(-2.8019866) q[0];
sx q[0];
rz(2.8787676) q[0];
rz(-0.35573959) q[1];
sx q[1];
rz(-1.4590615) q[1];
sx q[1];
rz(1.9062818) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4609769) q[0];
sx q[0];
rz(-1.5605698) q[0];
sx q[0];
rz(3.1264722) q[0];
x q[1];
rz(1.6791061) q[2];
sx q[2];
rz(-0.24178594) q[2];
sx q[2];
rz(2.0568648) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.0318109) q[1];
sx q[1];
rz(-2.8019252) q[1];
sx q[1];
rz(0.32082816) q[1];
rz(-1.9800817) q[3];
sx q[3];
rz(-0.73771362) q[3];
sx q[3];
rz(-1.1324594) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.1245023) q[2];
sx q[2];
rz(-1.3829145) q[2];
sx q[2];
rz(2.9612605) q[2];
rz(2.5056433) q[3];
sx q[3];
rz(-1.162581) q[3];
sx q[3];
rz(2.7699871) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3746049) q[0];
sx q[0];
rz(-2.7174482) q[0];
sx q[0];
rz(-2.4225127) q[0];
rz(0.51302296) q[1];
sx q[1];
rz(-1.9388371) q[1];
sx q[1];
rz(2.1069353) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4661515) q[0];
sx q[0];
rz(-0.062128566) q[0];
sx q[0];
rz(-2.8176869) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.047495202) q[2];
sx q[2];
rz(-1.7058027) q[2];
sx q[2];
rz(-0.53802711) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.6369578) q[1];
sx q[1];
rz(-0.89244288) q[1];
sx q[1];
rz(1.9370609) q[1];
rz(-pi) q[2];
rz(0.33200522) q[3];
sx q[3];
rz(-2.6453291) q[3];
sx q[3];
rz(2.9018324) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.2085312) q[2];
sx q[2];
rz(-2.8716817) q[2];
sx q[2];
rz(-0.2362403) q[2];
rz(-1.9832206) q[3];
sx q[3];
rz(-2.1393496) q[3];
sx q[3];
rz(-0.85197824) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1119969) q[0];
sx q[0];
rz(-1.010226) q[0];
sx q[0];
rz(0.90233666) q[0];
rz(0.92102712) q[1];
sx q[1];
rz(-0.61704707) q[1];
sx q[1];
rz(2.5193118) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3023923) q[0];
sx q[0];
rz(-1.7532187) q[0];
sx q[0];
rz(1.2829885) q[0];
rz(-pi) q[1];
rz(0.64825443) q[2];
sx q[2];
rz(-1.6443242) q[2];
sx q[2];
rz(-1.8258121) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.8473377) q[1];
sx q[1];
rz(-1.6604742) q[1];
sx q[1];
rz(-1.7404106) q[1];
rz(-pi) q[2];
rz(-1.1774109) q[3];
sx q[3];
rz(-2.0188361) q[3];
sx q[3];
rz(-0.084669948) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.92419147) q[2];
sx q[2];
rz(-1.8175586) q[2];
sx q[2];
rz(3.1203111) q[2];
rz(1.4422669) q[3];
sx q[3];
rz(-2.4271624) q[3];
sx q[3];
rz(-2.2309979) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6570046) q[0];
sx q[0];
rz(-0.5963043) q[0];
sx q[0];
rz(3.0969627) q[0];
rz(1.0021098) q[1];
sx q[1];
rz(-2.1410746) q[1];
sx q[1];
rz(-3.1071641) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9653942) q[0];
sx q[0];
rz(-1.584443) q[0];
sx q[0];
rz(-1.3854881) q[0];
rz(-0.12182932) q[2];
sx q[2];
rz(-1.0708772) q[2];
sx q[2];
rz(-0.027539754) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.5163102) q[1];
sx q[1];
rz(-1.0370266) q[1];
sx q[1];
rz(2.7620402) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4061635) q[3];
sx q[3];
rz(-2.2693686) q[3];
sx q[3];
rz(-2.5005831) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.78952152) q[2];
sx q[2];
rz(-1.8692724) q[2];
sx q[2];
rz(0.96674031) q[2];
rz(0.012185193) q[3];
sx q[3];
rz(-2.2120078) q[3];
sx q[3];
rz(-3.0325082) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1352585) q[0];
sx q[0];
rz(-2.874458) q[0];
sx q[0];
rz(-0.99408856) q[0];
rz(-0.055123568) q[1];
sx q[1];
rz(-1.8693285) q[1];
sx q[1];
rz(2.4023043) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2544884) q[0];
sx q[0];
rz(-1.4758037) q[0];
sx q[0];
rz(-1.4332921) q[0];
x q[1];
rz(-2.0544858) q[2];
sx q[2];
rz(-2.0927883) q[2];
sx q[2];
rz(2.1453478) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.90034396) q[1];
sx q[1];
rz(-2.5785891) q[1];
sx q[1];
rz(-2.7808933) q[1];
rz(1.8592632) q[3];
sx q[3];
rz(-2.1707284) q[3];
sx q[3];
rz(-0.44655061) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.58549515) q[2];
sx q[2];
rz(-1.7332417) q[2];
sx q[2];
rz(0.56376702) q[2];
rz(2.901315) q[3];
sx q[3];
rz(-2.6085745) q[3];
sx q[3];
rz(1.3747922) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.69650841) q[0];
sx q[0];
rz(-0.17906469) q[0];
sx q[0];
rz(2.5575496) q[0];
rz(0.42075992) q[1];
sx q[1];
rz(-2.1003484) q[1];
sx q[1];
rz(-2.8631794) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.70008343) q[0];
sx q[0];
rz(-2.1372876) q[0];
sx q[0];
rz(0.14871116) q[0];
rz(-2.7114696) q[2];
sx q[2];
rz(-0.42818907) q[2];
sx q[2];
rz(2.6870514) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.486057) q[1];
sx q[1];
rz(-1.2206519) q[1];
sx q[1];
rz(2.9085367) q[1];
rz(-1.1823468) q[3];
sx q[3];
rz(-2.5591345) q[3];
sx q[3];
rz(2.9782481) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.2723096) q[2];
sx q[2];
rz(-0.81988207) q[2];
sx q[2];
rz(0.38796866) q[2];
rz(-1.3502454) q[3];
sx q[3];
rz(-2.6781121) q[3];
sx q[3];
rz(-0.2750245) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.85912722) q[0];
sx q[0];
rz(-2.9282741) q[0];
sx q[0];
rz(2.4097089) q[0];
rz(0.16648079) q[1];
sx q[1];
rz(-2.6961168) q[1];
sx q[1];
rz(-0.073721185) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8548944) q[0];
sx q[0];
rz(-2.694283) q[0];
sx q[0];
rz(1.3351424) q[0];
rz(-1.6905626) q[2];
sx q[2];
rz(-1.4904067) q[2];
sx q[2];
rz(0.62270852) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.8970866) q[1];
sx q[1];
rz(-0.58468854) q[1];
sx q[1];
rz(3.0737682) q[1];
x q[2];
rz(-0.57314408) q[3];
sx q[3];
rz(-1.6640275) q[3];
sx q[3];
rz(-1.4179109) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.59447294) q[2];
sx q[2];
rz(-2.9056845) q[2];
sx q[2];
rz(0.42868844) q[2];
rz(-1.3171014) q[3];
sx q[3];
rz(-1.3573815) q[3];
sx q[3];
rz(1.930442) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
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
rz(-0.73228943) q[0];
sx q[0];
rz(-1.6560873) q[0];
sx q[0];
rz(2.0987341) q[0];
rz(1.6304784) q[1];
sx q[1];
rz(-1.7621367) q[1];
sx q[1];
rz(-1.7932549) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3955584) q[0];
sx q[0];
rz(-1.4575882) q[0];
sx q[0];
rz(-0.025870196) q[0];
rz(-1.6595608) q[2];
sx q[2];
rz(-0.73854337) q[2];
sx q[2];
rz(-2.2752787) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.1765602) q[1];
sx q[1];
rz(-1.6415879) q[1];
sx q[1];
rz(-1.895158) q[1];
rz(-0.24095778) q[3];
sx q[3];
rz(-1.7562508) q[3];
sx q[3];
rz(-2.0897739) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.80031359) q[2];
sx q[2];
rz(-2.6276402) q[2];
sx q[2];
rz(-0.094816118) q[2];
rz(1.173165) q[3];
sx q[3];
rz(-1.7404107) q[3];
sx q[3];
rz(0.15869424) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.512758) q[0];
sx q[0];
rz(-2.6678968) q[0];
sx q[0];
rz(-2.0976023) q[0];
rz(-1.6013153) q[1];
sx q[1];
rz(-1.5581144) q[1];
sx q[1];
rz(1.5099572) q[1];
rz(-1.2039456) q[2];
sx q[2];
rz(-0.22782142) q[2];
sx q[2];
rz(0.08624764) q[2];
rz(0.52900984) q[3];
sx q[3];
rz(-0.88376868) q[3];
sx q[3];
rz(0.39467011) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
