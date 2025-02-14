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
rz(-2.4189334) q[0];
sx q[0];
rz(1.0455796) q[0];
rz(0.0027520952) q[1];
sx q[1];
rz(-1.4581008) q[1];
sx q[1];
rz(2.4207065) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3430735) q[0];
sx q[0];
rz(-1.6630173) q[0];
sx q[0];
rz(-0.8999805) q[0];
rz(1.7621668) q[2];
sx q[2];
rz(-1.6211133) q[2];
sx q[2];
rz(-2.2309395) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.406479) q[1];
sx q[1];
rz(-1.8476474) q[1];
sx q[1];
rz(-2.7471354) q[1];
rz(-pi) q[2];
rz(-1.4414532) q[3];
sx q[3];
rz(-2.3096519) q[3];
sx q[3];
rz(-1.8042313) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.8568933) q[2];
sx q[2];
rz(-1.2130986) q[2];
sx q[2];
rz(2.6846679) q[2];
rz(0.65447718) q[3];
sx q[3];
rz(-0.5623397) q[3];
sx q[3];
rz(-0.42403179) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65861312) q[0];
sx q[0];
rz(-2.787866) q[0];
sx q[0];
rz(2.5335627) q[0];
rz(-1.9425302) q[1];
sx q[1];
rz(-2.6741195) q[1];
sx q[1];
rz(-0.79493585) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6632961) q[0];
sx q[0];
rz(-1.9062055) q[0];
sx q[0];
rz(1.2637423) q[0];
rz(-1.5995922) q[2];
sx q[2];
rz(-1.8378647) q[2];
sx q[2];
rz(-1.1065337) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.1767456) q[1];
sx q[1];
rz(-2.1934192) q[1];
sx q[1];
rz(-2.7962901) q[1];
rz(-pi) q[2];
rz(-0.43585082) q[3];
sx q[3];
rz(-0.66879771) q[3];
sx q[3];
rz(2.8613402) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.46637154) q[2];
sx q[2];
rz(-2.0872865) q[2];
sx q[2];
rz(1.6611453) q[2];
rz(-2.772707) q[3];
sx q[3];
rz(-1.8921013) q[3];
sx q[3];
rz(-0.7307542) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9166301) q[0];
sx q[0];
rz(-2.1376305) q[0];
sx q[0];
rz(0.68600255) q[0];
rz(-1.8983768) q[1];
sx q[1];
rz(-0.22626433) q[1];
sx q[1];
rz(-0.6023947) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7083654) q[0];
sx q[0];
rz(-2.4658563) q[0];
sx q[0];
rz(2.2883718) q[0];
x q[1];
rz(-1.5302646) q[2];
sx q[2];
rz(-1.5384595) q[2];
sx q[2];
rz(3.0258816) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.73302669) q[1];
sx q[1];
rz(-1.0045718) q[1];
sx q[1];
rz(0.11519656) q[1];
rz(-pi) q[2];
rz(-2.0394348) q[3];
sx q[3];
rz(-1.8640362) q[3];
sx q[3];
rz(1.4761497) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.7556222) q[2];
sx q[2];
rz(-1.3034416) q[2];
sx q[2];
rz(1.2516359) q[2];
rz(-0.55164591) q[3];
sx q[3];
rz(-0.14277221) q[3];
sx q[3];
rz(0.84757203) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8400693) q[0];
sx q[0];
rz(-1.9020377) q[0];
sx q[0];
rz(3.0384592) q[0];
rz(-2.3221807) q[1];
sx q[1];
rz(-0.90140072) q[1];
sx q[1];
rz(-0.55319667) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5846636) q[0];
sx q[0];
rz(-1.2069261) q[0];
sx q[0];
rz(-0.30839021) q[0];
rz(-pi) q[1];
rz(2.791094) q[2];
sx q[2];
rz(-1.0523049) q[2];
sx q[2];
rz(2.0295335) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.4955289) q[1];
sx q[1];
rz(-2.6501982) q[1];
sx q[1];
rz(1.8629406) q[1];
rz(-pi) q[2];
rz(-1.1392713) q[3];
sx q[3];
rz(-1.6628886) q[3];
sx q[3];
rz(-1.0669277) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.72054046) q[2];
sx q[2];
rz(-2.4002176) q[2];
sx q[2];
rz(-2.6067624) q[2];
rz(-0.77830166) q[3];
sx q[3];
rz(-2.4141267) q[3];
sx q[3];
rz(0.21540575) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4541723) q[0];
sx q[0];
rz(-1.9947616) q[0];
sx q[0];
rz(0.73690328) q[0];
rz(-0.058628254) q[1];
sx q[1];
rz(-2.3643989) q[1];
sx q[1];
rz(-2.5545205) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5694052) q[0];
sx q[0];
rz(-2.5068034) q[0];
sx q[0];
rz(1.6247747) q[0];
rz(0.063063085) q[2];
sx q[2];
rz(-1.7997671) q[2];
sx q[2];
rz(-1.582043) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.3526094) q[1];
sx q[1];
rz(-2.8375344) q[1];
sx q[1];
rz(1.4620146) q[1];
rz(-pi) q[2];
rz(1.2523417) q[3];
sx q[3];
rz(-1.8618004) q[3];
sx q[3];
rz(2.165447) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.8495142) q[2];
sx q[2];
rz(-2.2167315) q[2];
sx q[2];
rz(-2.4936131) q[2];
rz(-2.5255711) q[3];
sx q[3];
rz(-1.6391305) q[3];
sx q[3];
rz(-0.047671635) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7306526) q[0];
sx q[0];
rz(-2.3766282) q[0];
sx q[0];
rz(2.1867645) q[0];
rz(-3.094589) q[1];
sx q[1];
rz(-2.4689597) q[1];
sx q[1];
rz(2.1254553) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.797182) q[0];
sx q[0];
rz(-1.1007995) q[0];
sx q[0];
rz(2.5380666) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1121641) q[2];
sx q[2];
rz(-0.9230136) q[2];
sx q[2];
rz(-1.7769396) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.8879442) q[1];
sx q[1];
rz(-1.5563111) q[1];
sx q[1];
rz(-2.6360341) q[1];
rz(-pi) q[2];
rz(-3.0623779) q[3];
sx q[3];
rz(-0.62566775) q[3];
sx q[3];
rz(0.46712671) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.59812349) q[2];
sx q[2];
rz(-1.2657961) q[2];
sx q[2];
rz(3.0908965) q[2];
rz(-3.0905881) q[3];
sx q[3];
rz(-1.3011322) q[3];
sx q[3];
rz(-0.72371975) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8766668) q[0];
sx q[0];
rz(-0.51487881) q[0];
sx q[0];
rz(1.6616954) q[0];
rz(-2.0199846) q[1];
sx q[1];
rz(-1.879296) q[1];
sx q[1];
rz(-0.46953896) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.74140841) q[0];
sx q[0];
rz(-3.0245076) q[0];
sx q[0];
rz(-1.9817121) q[0];
x q[1];
rz(1.93365) q[2];
sx q[2];
rz(-2.7245625) q[2];
sx q[2];
rz(-0.94379683) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.6958298) q[1];
sx q[1];
rz(-1.1053021) q[1];
sx q[1];
rz(2.9470934) q[1];
rz(2.4933563) q[3];
sx q[3];
rz(-2.2508374) q[3];
sx q[3];
rz(1.4095962) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.7349714) q[2];
sx q[2];
rz(-1.5512286) q[2];
sx q[2];
rz(1.1173908) q[2];
rz(-2.0495074) q[3];
sx q[3];
rz(-1.7376244) q[3];
sx q[3];
rz(-0.4357858) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
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
rz(-1.0206873) q[0];
sx q[0];
rz(-2.0285323) q[0];
sx q[0];
rz(0.06509617) q[0];
rz(0.33196017) q[1];
sx q[1];
rz(-1.8354225) q[1];
sx q[1];
rz(-1.8665705) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41439636) q[0];
sx q[0];
rz(-2.9119618) q[0];
sx q[0];
rz(-1.5288607) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0135856) q[2];
sx q[2];
rz(-1.6943568) q[2];
sx q[2];
rz(-0.75109252) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.4091585) q[1];
sx q[1];
rz(-2.2526764) q[1];
sx q[1];
rz(1.5049008) q[1];
rz(-pi) q[2];
rz(-1.7528107) q[3];
sx q[3];
rz(-2.3654147) q[3];
sx q[3];
rz(2.9266199) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.8518105) q[2];
sx q[2];
rz(-1.6921356) q[2];
sx q[2];
rz(0.97274485) q[2];
rz(2.4408477) q[3];
sx q[3];
rz(-2.2795129) q[3];
sx q[3];
rz(-2.3166482) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7026611) q[0];
sx q[0];
rz(-2.7362566) q[0];
sx q[0];
rz(-2.7324556) q[0];
rz(-1.9955955) q[1];
sx q[1];
rz(-2.0775677) q[1];
sx q[1];
rz(-1.2629868) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1564558) q[0];
sx q[0];
rz(-3.1334152) q[0];
sx q[0];
rz(-1.5865492) q[0];
rz(2.1705698) q[2];
sx q[2];
rz(-0.46659887) q[2];
sx q[2];
rz(-0.80228892) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.9085051) q[1];
sx q[1];
rz(-0.89492866) q[1];
sx q[1];
rz(0.29260537) q[1];
x q[2];
rz(3.0060482) q[3];
sx q[3];
rz(-2.5672847) q[3];
sx q[3];
rz(1.7611461) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.9001749) q[2];
sx q[2];
rz(-1.9200385) q[2];
sx q[2];
rz(2.9202666) q[2];
rz(0.48702249) q[3];
sx q[3];
rz(-2.2723787) q[3];
sx q[3];
rz(1.937449) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0601592) q[0];
sx q[0];
rz(-2.8964323) q[0];
sx q[0];
rz(-1.2231476) q[0];
rz(3.0606015) q[1];
sx q[1];
rz(-1.5206189) q[1];
sx q[1];
rz(-1.5222668) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3190614) q[0];
sx q[0];
rz(-0.37320787) q[0];
sx q[0];
rz(-2.6749956) q[0];
rz(-pi) q[1];
rz(2.4967983) q[2];
sx q[2];
rz(-0.83907467) q[2];
sx q[2];
rz(-0.085345833) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.63357089) q[1];
sx q[1];
rz(-2.7765818) q[1];
sx q[1];
rz(0.99442039) q[1];
rz(-pi) q[2];
rz(-0.31701458) q[3];
sx q[3];
rz(-2.6251589) q[3];
sx q[3];
rz(1.7805733) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.12485518) q[2];
sx q[2];
rz(-1.3913466) q[2];
sx q[2];
rz(-0.69046956) q[2];
rz(-0.27103439) q[3];
sx q[3];
rz(-2.6764968) q[3];
sx q[3];
rz(-1.5439699) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4616213) q[0];
sx q[0];
rz(-2.2053056) q[0];
sx q[0];
rz(2.9619138) q[0];
rz(2.8996254) q[1];
sx q[1];
rz(-0.89121834) q[1];
sx q[1];
rz(-1.2527087) q[1];
rz(2.719328) q[2];
sx q[2];
rz(-1.3801856) q[2];
sx q[2];
rz(-1.2587067) q[2];
rz(-1.1923424) q[3];
sx q[3];
rz(-0.71810953) q[3];
sx q[3];
rz(0.22267846) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
