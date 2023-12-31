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
rz(-2.1687578) q[1];
sx q[1];
rz(-1.6701148) q[1];
sx q[1];
rz(2.8491128) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7115294) q[0];
sx q[0];
rz(-1.9646514) q[0];
sx q[0];
rz(-0.060056134) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5195261) q[2];
sx q[2];
rz(-2.3495557) q[2];
sx q[2];
rz(2.771332) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.18947345) q[1];
sx q[1];
rz(-1.6166302) q[1];
sx q[1];
rz(-0.97656753) q[1];
rz(-3.1396477) q[3];
sx q[3];
rz(-2.0446152) q[3];
sx q[3];
rz(-2.4634944) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.4108489) q[2];
sx q[2];
rz(-1.0935254) q[2];
sx q[2];
rz(-2.6921819) q[2];
rz(0.64569008) q[3];
sx q[3];
rz(-0.68104762) q[3];
sx q[3];
rz(1.8507563) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
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
rz(-1.6702601) q[1];
sx q[1];
rz(-2.4447618) q[1];
sx q[1];
rz(1.3630294) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8105158) q[0];
sx q[0];
rz(-1.9609945) q[0];
sx q[0];
rz(-0.66731989) q[0];
rz(-1.7861373) q[2];
sx q[2];
rz(-2.0305579) q[2];
sx q[2];
rz(2.721867) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.431753) q[1];
sx q[1];
rz(-1.589994) q[1];
sx q[1];
rz(0.33961105) q[1];
x q[2];
rz(-1.3046673) q[3];
sx q[3];
rz(-1.3849272) q[3];
sx q[3];
rz(0.64418018) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.30119511) q[2];
sx q[2];
rz(-2.6448554) q[2];
sx q[2];
rz(-1.8691241) q[2];
rz(0.33723351) q[3];
sx q[3];
rz(-1.971258) q[3];
sx q[3];
rz(3.1315394) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0740046) q[0];
sx q[0];
rz(-0.33960605) q[0];
sx q[0];
rz(-0.2628251) q[0];
rz(0.35573959) q[1];
sx q[1];
rz(-1.6825312) q[1];
sx q[1];
rz(1.9062818) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10997406) q[0];
sx q[0];
rz(-1.5556766) q[0];
sx q[0];
rz(-1.5605687) q[0];
x q[1];
rz(1.6791061) q[2];
sx q[2];
rz(-0.24178594) q[2];
sx q[2];
rz(-1.0847278) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.0318109) q[1];
sx q[1];
rz(-0.33966741) q[1];
sx q[1];
rz(2.8207645) q[1];
x q[2];
rz(-0.87576207) q[3];
sx q[3];
rz(-1.8417629) q[3];
sx q[3];
rz(0.7489487) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.0170903) q[2];
sx q[2];
rz(-1.3829145) q[2];
sx q[2];
rz(0.18033218) q[2];
rz(2.5056433) q[3];
sx q[3];
rz(-1.9790117) q[3];
sx q[3];
rz(0.37160555) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3746049) q[0];
sx q[0];
rz(-2.7174482) q[0];
sx q[0];
rz(-2.4225127) q[0];
rz(-0.51302296) q[1];
sx q[1];
rz(-1.2027556) q[1];
sx q[1];
rz(-1.0346574) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2186787) q[0];
sx q[0];
rz(-1.5905587) q[0];
sx q[0];
rz(3.0826871) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7059533) q[2];
sx q[2];
rz(-1.5237336) q[2];
sx q[2];
rz(-2.1152209) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.3116341) q[1];
sx q[1];
rz(-1.8533851) q[1];
sx q[1];
rz(-2.4294873) q[1];
rz(-pi) q[2];
rz(-1.3961117) q[3];
sx q[3];
rz(-2.0377199) q[3];
sx q[3];
rz(-0.13388453) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.2085312) q[2];
sx q[2];
rz(-2.8716817) q[2];
sx q[2];
rz(-2.9053524) q[2];
rz(1.9832206) q[3];
sx q[3];
rz(-1.002243) q[3];
sx q[3];
rz(-0.85197824) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.029595705) q[0];
sx q[0];
rz(-1.010226) q[0];
sx q[0];
rz(2.239256) q[0];
rz(2.2205655) q[1];
sx q[1];
rz(-2.5245456) q[1];
sx q[1];
rz(-0.62228084) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.21475269) q[0];
sx q[0];
rz(-1.8536957) q[0];
sx q[0];
rz(-0.19006417) q[0];
rz(-0.64825443) q[2];
sx q[2];
rz(-1.4972685) q[2];
sx q[2];
rz(1.3157805) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.29187782) q[1];
sx q[1];
rz(-1.7397225) q[1];
sx q[1];
rz(-3.0506163) q[1];
x q[2];
rz(-1.1774109) q[3];
sx q[3];
rz(-2.0188361) q[3];
sx q[3];
rz(3.0569227) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.92419147) q[2];
sx q[2];
rz(-1.3240341) q[2];
sx q[2];
rz(3.1203111) q[2];
rz(-1.6993258) q[3];
sx q[3];
rz(-2.4271624) q[3];
sx q[3];
rz(0.91059476) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6570046) q[0];
sx q[0];
rz(-2.5452884) q[0];
sx q[0];
rz(-0.044629991) q[0];
rz(-2.1394829) q[1];
sx q[1];
rz(-1.0005181) q[1];
sx q[1];
rz(-0.034428509) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1761985) q[0];
sx q[0];
rz(-1.584443) q[0];
sx q[0];
rz(-1.7561046) q[0];
rz(-pi) q[1];
rz(0.12182932) q[2];
sx q[2];
rz(-2.0707154) q[2];
sx q[2];
rz(-0.027539754) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.9958482) q[1];
sx q[1];
rz(-1.8954344) q[1];
sx q[1];
rz(2.1374628) q[1];
rz(-pi) q[2];
rz(2.9488726) q[3];
sx q[3];
rz(-0.7145213) q[3];
sx q[3];
rz(2.2477828) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.78952152) q[2];
sx q[2];
rz(-1.2723203) q[2];
sx q[2];
rz(2.1748523) q[2];
rz(3.1294075) q[3];
sx q[3];
rz(-2.2120078) q[3];
sx q[3];
rz(-0.10908443) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1352585) q[0];
sx q[0];
rz(-2.874458) q[0];
sx q[0];
rz(0.99408856) q[0];
rz(0.055123568) q[1];
sx q[1];
rz(-1.8693285) q[1];
sx q[1];
rz(0.73928839) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2544884) q[0];
sx q[0];
rz(-1.6657889) q[0];
sx q[0];
rz(-1.7083005) q[0];
rz(-pi) q[1];
rz(0.67989345) q[2];
sx q[2];
rz(-0.69603622) q[2];
sx q[2];
rz(1.807715) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.1623605) q[1];
sx q[1];
rz(-1.7602966) q[1];
sx q[1];
rz(-0.53342553) q[1];
rz(0.39412734) q[3];
sx q[3];
rz(-2.4836802) q[3];
sx q[3];
rz(3.1042299) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.5560975) q[2];
sx q[2];
rz(-1.408351) q[2];
sx q[2];
rz(-0.56376702) q[2];
rz(0.24027763) q[3];
sx q[3];
rz(-0.53301817) q[3];
sx q[3];
rz(1.3747922) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4450842) q[0];
sx q[0];
rz(-2.962528) q[0];
sx q[0];
rz(-0.58404303) q[0];
rz(-0.42075992) q[1];
sx q[1];
rz(-2.1003484) q[1];
sx q[1];
rz(2.8631794) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42785545) q[0];
sx q[0];
rz(-2.5579778) q[0];
sx q[0];
rz(-1.3419271) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7483447) q[2];
sx q[2];
rz(-1.7448145) q[2];
sx q[2];
rz(-2.4207123) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.091374) q[1];
sx q[1];
rz(-2.723657) q[1];
sx q[1];
rz(-2.1347079) q[1];
rz(-pi) q[2];
x q[2];
rz(2.897103) q[3];
sx q[3];
rz(-2.1049307) q[3];
sx q[3];
rz(0.29230803) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.86928308) q[2];
sx q[2];
rz(-2.3217106) q[2];
sx q[2];
rz(0.38796866) q[2];
rz(-1.3502454) q[3];
sx q[3];
rz(-0.46348059) q[3];
sx q[3];
rz(0.2750245) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.85912722) q[0];
sx q[0];
rz(-0.21331856) q[0];
sx q[0];
rz(0.7318837) q[0];
rz(-2.9751119) q[1];
sx q[1];
rz(-2.6961168) q[1];
sx q[1];
rz(-0.073721185) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28669824) q[0];
sx q[0];
rz(-0.44730967) q[0];
sx q[0];
rz(1.3351424) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.45103) q[2];
sx q[2];
rz(-1.6511859) q[2];
sx q[2];
rz(0.62270852) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.8970866) q[1];
sx q[1];
rz(-2.5569041) q[1];
sx q[1];
rz(3.0737682) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.681626) q[3];
sx q[3];
rz(-2.14114) q[3];
sx q[3];
rz(-0.092872083) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.5471197) q[2];
sx q[2];
rz(-2.9056845) q[2];
sx q[2];
rz(2.7129042) q[2];
rz(1.3171014) q[3];
sx q[3];
rz(-1.7842112) q[3];
sx q[3];
rz(-1.2111506) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4093032) q[0];
sx q[0];
rz(-1.6560873) q[0];
sx q[0];
rz(2.0987341) q[0];
rz(1.5111142) q[1];
sx q[1];
rz(-1.379456) q[1];
sx q[1];
rz(-1.7932549) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6207328) q[0];
sx q[0];
rz(-3.0254786) q[0];
sx q[0];
rz(-1.7945047) q[0];
rz(-pi) q[1];
rz(3.0610602) q[2];
sx q[2];
rz(-0.83583145) q[2];
sx q[2];
rz(-0.9860952) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.9546982) q[1];
sx q[1];
rz(-0.3317301) q[1];
sx q[1];
rz(1.3518672) q[1];
rz(-pi) q[2];
rz(-2.9006349) q[3];
sx q[3];
rz(-1.7562508) q[3];
sx q[3];
rz(2.0897739) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.3412791) q[2];
sx q[2];
rz(-0.51395243) q[2];
sx q[2];
rz(-3.0467765) q[2];
rz(-1.173165) q[3];
sx q[3];
rz(-1.7404107) q[3];
sx q[3];
rz(2.9828984) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6288347) q[0];
sx q[0];
rz(-0.47369581) q[0];
sx q[0];
rz(1.0439903) q[0];
rz(1.5402773) q[1];
sx q[1];
rz(-1.5581144) q[1];
sx q[1];
rz(1.5099572) q[1];
rz(0.082967233) q[2];
sx q[2];
rz(-1.7832179) q[2];
sx q[2];
rz(-0.28945343) q[2];
rz(2.6125828) q[3];
sx q[3];
rz(-2.257824) q[3];
sx q[3];
rz(-2.7469225) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
