OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.2330115) q[0];
sx q[0];
rz(-1.1865948) q[0];
sx q[0];
rz(-1.9957805) q[0];
rz(-2.1687578) q[1];
sx q[1];
rz(-1.6701148) q[1];
sx q[1];
rz(2.8491128) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8669517) q[0];
sx q[0];
rz(-0.39817087) q[0];
sx q[0];
rz(1.7142332) q[0];
x q[1];
rz(2.3621759) q[2];
sx q[2];
rz(-1.6072818) q[2];
sx q[2];
rz(1.2365637) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.8279743) q[1];
sx q[1];
rz(-2.5458114) q[1];
sx q[1];
rz(-1.4890563) q[1];
rz(0.0019449751) q[3];
sx q[3];
rz(-1.0969775) q[3];
sx q[3];
rz(2.4634944) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.73074377) q[2];
sx q[2];
rz(-2.0480672) q[2];
sx q[2];
rz(0.4494108) q[2];
rz(2.4959026) q[3];
sx q[3];
rz(-2.460545) q[3];
sx q[3];
rz(-1.2908363) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.21355024) q[0];
sx q[0];
rz(-0.88590652) q[0];
sx q[0];
rz(-2.399562) q[0];
rz(1.6702601) q[1];
sx q[1];
rz(-0.69683087) q[1];
sx q[1];
rz(1.3630294) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5308967) q[0];
sx q[0];
rz(-0.96141059) q[0];
sx q[0];
rz(1.0884398) q[0];
rz(-pi) q[1];
rz(0.40740168) q[2];
sx q[2];
rz(-2.6371837) q[2];
sx q[2];
rz(-0.87770578) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.14582536) q[1];
sx q[1];
rz(-1.9103423) q[1];
sx q[1];
rz(1.550436) q[1];
x q[2];
rz(1.8369254) q[3];
sx q[3];
rz(-1.3849272) q[3];
sx q[3];
rz(0.64418018) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.8403975) q[2];
sx q[2];
rz(-0.49673721) q[2];
sx q[2];
rz(-1.2724686) q[2];
rz(0.33723351) q[3];
sx q[3];
rz(-1.971258) q[3];
sx q[3];
rz(3.1315394) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.067588016) q[0];
sx q[0];
rz(-0.33960605) q[0];
sx q[0];
rz(2.8787676) q[0];
rz(-0.35573959) q[1];
sx q[1];
rz(-1.4590615) q[1];
sx q[1];
rz(1.9062818) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4848029) q[0];
sx q[0];
rz(-0.018253837) q[0];
sx q[0];
rz(-0.59469964) q[0];
rz(-pi) q[1];
rz(-0.026651816) q[2];
sx q[2];
rz(-1.3304552) q[2];
sx q[2];
rz(-1.196256) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.23535145) q[1];
sx q[1];
rz(-1.6760577) q[1];
sx q[1];
rz(2.8180442) q[1];
rz(1.1615109) q[3];
sx q[3];
rz(-0.73771362) q[3];
sx q[3];
rz(2.0091332) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.0170903) q[2];
sx q[2];
rz(-1.7586781) q[2];
sx q[2];
rz(-2.9612605) q[2];
rz(-0.63594937) q[3];
sx q[3];
rz(-1.162581) q[3];
sx q[3];
rz(-0.37160555) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76698774) q[0];
sx q[0];
rz(-0.42414442) q[0];
sx q[0];
rz(2.4225127) q[0];
rz(2.6285697) q[1];
sx q[1];
rz(-1.2027556) q[1];
sx q[1];
rz(2.1069353) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4661515) q[0];
sx q[0];
rz(-3.0794641) q[0];
sx q[0];
rz(-2.8176869) q[0];
rz(-pi) q[1];
x q[1];
rz(0.047495202) q[2];
sx q[2];
rz(-1.7058027) q[2];
sx q[2];
rz(0.53802711) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.0532916) q[1];
sx q[1];
rz(-0.75690714) q[1];
sx q[1];
rz(-2.7234368) q[1];
x q[2];
rz(1.3961117) q[3];
sx q[3];
rz(-2.0377199) q[3];
sx q[3];
rz(0.13388453) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.93306142) q[2];
sx q[2];
rz(-2.8716817) q[2];
sx q[2];
rz(-0.2362403) q[2];
rz(1.158372) q[3];
sx q[3];
rz(-2.1393496) q[3];
sx q[3];
rz(-0.85197824) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1119969) q[0];
sx q[0];
rz(-2.1313666) q[0];
sx q[0];
rz(2.239256) q[0];
rz(0.92102712) q[1];
sx q[1];
rz(-2.5245456) q[1];
sx q[1];
rz(-2.5193118) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3233933) q[0];
sx q[0];
rz(-2.8022031) q[0];
sx q[0];
rz(-2.1470977) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4933382) q[2];
sx q[2];
rz(-1.4972685) q[2];
sx q[2];
rz(-1.3157805) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.29187782) q[1];
sx q[1];
rz(-1.4018702) q[1];
sx q[1];
rz(-3.0506163) q[1];
rz(2.6617674) q[3];
sx q[3];
rz(-1.2180426) q[3];
sx q[3];
rz(-1.6640116) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.92419147) q[2];
sx q[2];
rz(-1.3240341) q[2];
sx q[2];
rz(-3.1203111) q[2];
rz(1.6993258) q[3];
sx q[3];
rz(-0.71443021) q[3];
sx q[3];
rz(0.91059476) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4845881) q[0];
sx q[0];
rz(-2.5452884) q[0];
sx q[0];
rz(-3.0969627) q[0];
rz(-2.1394829) q[1];
sx q[1];
rz(-1.0005181) q[1];
sx q[1];
rz(-0.034428509) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9653942) q[0];
sx q[0];
rz(-1.584443) q[0];
sx q[0];
rz(1.7561046) q[0];
rz(1.0677412) q[2];
sx q[2];
rz(-1.6776553) q[2];
sx q[2];
rz(-1.5397132) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.5163102) q[1];
sx q[1];
rz(-1.0370266) q[1];
sx q[1];
rz(-2.7620402) q[1];
x q[2];
rz(-2.4363082) q[3];
sx q[3];
rz(-1.6966288) q[3];
sx q[3];
rz(0.82334405) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.78952152) q[2];
sx q[2];
rz(-1.2723203) q[2];
sx q[2];
rz(-2.1748523) q[2];
rz(-0.012185193) q[3];
sx q[3];
rz(-2.2120078) q[3];
sx q[3];
rz(-0.10908443) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0063342) q[0];
sx q[0];
rz(-2.874458) q[0];
sx q[0];
rz(0.99408856) q[0];
rz(-3.0864691) q[1];
sx q[1];
rz(-1.2722641) q[1];
sx q[1];
rz(-0.73928839) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.71541876) q[0];
sx q[0];
rz(-2.9746375) q[0];
sx q[0];
rz(-0.963361) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0544858) q[2];
sx q[2];
rz(-1.0488044) q[2];
sx q[2];
rz(-0.99624485) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.1623605) q[1];
sx q[1];
rz(-1.3812961) q[1];
sx q[1];
rz(-0.53342553) q[1];
x q[2];
rz(2.7474653) q[3];
sx q[3];
rz(-2.4836802) q[3];
sx q[3];
rz(0.037362785) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.5560975) q[2];
sx q[2];
rz(-1.7332417) q[2];
sx q[2];
rz(2.5778256) q[2];
rz(-0.24027763) q[3];
sx q[3];
rz(-2.6085745) q[3];
sx q[3];
rz(-1.7668004) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69650841) q[0];
sx q[0];
rz(-2.962528) q[0];
sx q[0];
rz(2.5575496) q[0];
rz(2.7208327) q[1];
sx q[1];
rz(-2.1003484) q[1];
sx q[1];
rz(2.8631794) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1906492) q[0];
sx q[0];
rz(-1.4454495) q[0];
sx q[0];
rz(0.9992674) q[0];
rz(-0.43012302) q[2];
sx q[2];
rz(-0.42818907) q[2];
sx q[2];
rz(0.45454121) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.091374) q[1];
sx q[1];
rz(-2.723657) q[1];
sx q[1];
rz(1.0068847) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0233381) q[3];
sx q[3];
rz(-1.3609144) q[3];
sx q[3];
rz(1.4048214) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.86928308) q[2];
sx q[2];
rz(-2.3217106) q[2];
sx q[2];
rz(0.38796866) q[2];
rz(-1.7913473) q[3];
sx q[3];
rz(-0.46348059) q[3];
sx q[3];
rz(2.8665682) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2824654) q[0];
sx q[0];
rz(-0.21331856) q[0];
sx q[0];
rz(2.4097089) q[0];
rz(2.9751119) q[1];
sx q[1];
rz(-2.6961168) q[1];
sx q[1];
rz(0.073721185) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28669824) q[0];
sx q[0];
rz(-0.44730967) q[0];
sx q[0];
rz(-1.8064503) q[0];
x q[1];
rz(0.080967112) q[2];
sx q[2];
rz(-1.4514187) q[2];
sx q[2];
rz(-0.93842426) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.32578711) q[1];
sx q[1];
rz(-0.98762883) q[1];
sx q[1];
rz(-1.5259685) q[1];
rz(2.5684486) q[3];
sx q[3];
rz(-1.4775652) q[3];
sx q[3];
rz(-1.7236818) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.5471197) q[2];
sx q[2];
rz(-2.9056845) q[2];
sx q[2];
rz(-2.7129042) q[2];
rz(-1.8244913) q[3];
sx q[3];
rz(-1.3573815) q[3];
sx q[3];
rz(-1.930442) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.73228943) q[0];
sx q[0];
rz(-1.4855054) q[0];
sx q[0];
rz(-1.0428585) q[0];
rz(-1.5111142) q[1];
sx q[1];
rz(-1.379456) q[1];
sx q[1];
rz(1.7932549) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.82768517) q[0];
sx q[0];
rz(-1.5450918) q[0];
sx q[0];
rz(1.4575507) q[0];
rz(-pi) q[1];
rz(2.3073763) q[2];
sx q[2];
rz(-1.5110821) q[2];
sx q[2];
rz(-2.5028253) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.18689449) q[1];
sx q[1];
rz(-2.8098626) q[1];
sx q[1];
rz(-1.7897254) q[1];
rz(1.3799558) q[3];
sx q[3];
rz(-1.8075426) q[3];
sx q[3];
rz(-2.5773347) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.3412791) q[2];
sx q[2];
rz(-2.6276402) q[2];
sx q[2];
rz(3.0467765) q[2];
rz(-1.173165) q[3];
sx q[3];
rz(-1.4011819) q[3];
sx q[3];
rz(-2.9828984) q[3];
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
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
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
rz(-1.7839292) q[2];
sx q[2];
rz(-1.6518946) q[2];
sx q[2];
rz(1.2988731) q[2];
rz(0.81090609) q[3];
sx q[3];
rz(-1.1699642) q[3];
sx q[3];
rz(-1.5311833) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];