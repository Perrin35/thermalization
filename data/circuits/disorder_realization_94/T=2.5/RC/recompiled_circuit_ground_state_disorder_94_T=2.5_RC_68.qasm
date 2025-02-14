OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.2296978) q[0];
sx q[0];
rz(-1.8954281) q[0];
sx q[0];
rz(1.5204313) q[0];
rz(-2.7819832) q[1];
sx q[1];
rz(-0.26878992) q[1];
sx q[1];
rz(0.51528817) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.775151) q[0];
sx q[0];
rz(-0.40344813) q[0];
sx q[0];
rz(2.6029384) q[0];
rz(1.8727539) q[2];
sx q[2];
rz(-2.0965323) q[2];
sx q[2];
rz(0.56509226) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.8003502) q[1];
sx q[1];
rz(-1.9669232) q[1];
sx q[1];
rz(1.800283) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4024967) q[3];
sx q[3];
rz(-1.0217654) q[3];
sx q[3];
rz(-2.9573553) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.9109853) q[2];
sx q[2];
rz(-1.5585941) q[2];
sx q[2];
rz(0.67414635) q[2];
rz(-2.9437183) q[3];
sx q[3];
rz(-1.9015692) q[3];
sx q[3];
rz(1.6363293) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3182217) q[0];
sx q[0];
rz(-2.8405393) q[0];
sx q[0];
rz(0.2163042) q[0];
rz(3.0220616) q[1];
sx q[1];
rz(-2.0152338) q[1];
sx q[1];
rz(1.4215887) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7455121) q[0];
sx q[0];
rz(-2.0876309) q[0];
sx q[0];
rz(-0.091334657) q[0];
rz(2.6720409) q[2];
sx q[2];
rz(-0.44347635) q[2];
sx q[2];
rz(1.8592905) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.6338809) q[1];
sx q[1];
rz(-1.8844023) q[1];
sx q[1];
rz(-2.0416946) q[1];
rz(-pi) q[2];
rz(-0.64038527) q[3];
sx q[3];
rz(-2.3043568) q[3];
sx q[3];
rz(-0.094410019) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.8372832) q[2];
sx q[2];
rz(-1.9459566) q[2];
sx q[2];
rz(-0.99011123) q[2];
rz(0.75634161) q[3];
sx q[3];
rz(-0.6476616) q[3];
sx q[3];
rz(0.17624779) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15683098) q[0];
sx q[0];
rz(-2.3909843) q[0];
sx q[0];
rz(2.005715) q[0];
rz(-2.1319977) q[1];
sx q[1];
rz(-1.986809) q[1];
sx q[1];
rz(2.6957767) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8404482) q[0];
sx q[0];
rz(-3.0160286) q[0];
sx q[0];
rz(2.156394) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.03568825) q[2];
sx q[2];
rz(-0.70377398) q[2];
sx q[2];
rz(1.698871) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.7103084) q[1];
sx q[1];
rz(-0.7529707) q[1];
sx q[1];
rz(-2.1953039) q[1];
rz(-pi) q[2];
x q[2];
rz(0.11196158) q[3];
sx q[3];
rz(-0.44977934) q[3];
sx q[3];
rz(1.8729608) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.3710215) q[2];
sx q[2];
rz(-2.1393445) q[2];
sx q[2];
rz(-2.501343) q[2];
rz(0.69784969) q[3];
sx q[3];
rz(-1.6777439) q[3];
sx q[3];
rz(2.8016134) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3014389) q[0];
sx q[0];
rz(-0.88742632) q[0];
sx q[0];
rz(-1.3288757) q[0];
rz(1.9245194) q[1];
sx q[1];
rz(-0.71958676) q[1];
sx q[1];
rz(0.77635366) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8912131) q[0];
sx q[0];
rz(-0.73085143) q[0];
sx q[0];
rz(-2.5688085) q[0];
rz(-0.66993454) q[2];
sx q[2];
rz(-2.4342854) q[2];
sx q[2];
rz(-2.7279127) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.2393155) q[1];
sx q[1];
rz(-1.8493686) q[1];
sx q[1];
rz(0.011812731) q[1];
x q[2];
rz(2.4569974) q[3];
sx q[3];
rz(-1.5719613) q[3];
sx q[3];
rz(-0.62546414) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.2886469) q[2];
sx q[2];
rz(-0.32175803) q[2];
sx q[2];
rz(1.9256437) q[2];
rz(1.1207885) q[3];
sx q[3];
rz(-1.2633163) q[3];
sx q[3];
rz(-2.2585675) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(2.9057366) q[0];
sx q[0];
rz(-0.97705066) q[0];
sx q[0];
rz(2.6781154) q[0];
rz(1.0889168) q[1];
sx q[1];
rz(-1.8810279) q[1];
sx q[1];
rz(-1.3190528) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.70196296) q[0];
sx q[0];
rz(-1.3101398) q[0];
sx q[0];
rz(-1.8219276) q[0];
x q[1];
rz(-2.3438101) q[2];
sx q[2];
rz(-1.6897276) q[2];
sx q[2];
rz(-0.2695131) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.0214349) q[1];
sx q[1];
rz(-1.6579227) q[1];
sx q[1];
rz(2.1539262) q[1];
rz(1.1518258) q[3];
sx q[3];
rz(-1.7363461) q[3];
sx q[3];
rz(2.3922684) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.9702381) q[2];
sx q[2];
rz(-0.78533185) q[2];
sx q[2];
rz(2.5023517) q[2];
rz(-0.81131896) q[3];
sx q[3];
rz(-2.6526484) q[3];
sx q[3];
rz(1.5343687) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2163579) q[0];
sx q[0];
rz(-2.7282867) q[0];
sx q[0];
rz(0.21892029) q[0];
rz(1.9256598) q[1];
sx q[1];
rz(-0.464012) q[1];
sx q[1];
rz(-1.4930412) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7812778) q[0];
sx q[0];
rz(-2.6584315) q[0];
sx q[0];
rz(0.29668087) q[0];
x q[1];
rz(-0.41755192) q[2];
sx q[2];
rz(-1.5680015) q[2];
sx q[2];
rz(-2.8101943) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.5801489) q[1];
sx q[1];
rz(-2.3030048) q[1];
sx q[1];
rz(-1.5484018) q[1];
rz(-pi) q[2];
rz(0.099148765) q[3];
sx q[3];
rz(-1.5739023) q[3];
sx q[3];
rz(-0.24156027) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.1047989) q[2];
sx q[2];
rz(-1.3497738) q[2];
sx q[2];
rz(1.3118504) q[2];
rz(1.5051684) q[3];
sx q[3];
rz(-1.2994095) q[3];
sx q[3];
rz(2.6817491) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
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
rz(-0.63603193) q[0];
sx q[0];
rz(-2.6513031) q[0];
sx q[0];
rz(1.7946515) q[0];
rz(1.8147644) q[1];
sx q[1];
rz(-2.3074) q[1];
sx q[1];
rz(-0.38536513) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6198719) q[0];
sx q[0];
rz(-1.6024029) q[0];
sx q[0];
rz(-0.35913976) q[0];
rz(-pi) q[1];
rz(2.1009675) q[2];
sx q[2];
rz(-1.1405924) q[2];
sx q[2];
rz(-0.29597798) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.24033156) q[1];
sx q[1];
rz(-1.5442532) q[1];
sx q[1];
rz(0.44601299) q[1];
x q[2];
rz(-0.47355424) q[3];
sx q[3];
rz(-1.9723168) q[3];
sx q[3];
rz(-0.36916379) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.68606004) q[2];
sx q[2];
rz(-2.3251688) q[2];
sx q[2];
rz(-1.653999) q[2];
rz(-0.16573302) q[3];
sx q[3];
rz(-2.3159852) q[3];
sx q[3];
rz(-2.9674271) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1626749) q[0];
sx q[0];
rz(-0.62541494) q[0];
sx q[0];
rz(0.22853525) q[0];
rz(0.31271115) q[1];
sx q[1];
rz(-2.2622908) q[1];
sx q[1];
rz(1.3794587) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40423597) q[0];
sx q[0];
rz(-1.9725058) q[0];
sx q[0];
rz(1.6549003) q[0];
x q[1];
rz(-1.7422471) q[2];
sx q[2];
rz(-1.2510692) q[2];
sx q[2];
rz(2.7227744) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.5818565) q[1];
sx q[1];
rz(-0.35071555) q[1];
sx q[1];
rz(-1.6718256) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5614235) q[3];
sx q[3];
rz(-2.3456367) q[3];
sx q[3];
rz(-0.43886504) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.10869965) q[2];
sx q[2];
rz(-2.097082) q[2];
sx q[2];
rz(2.1006987) q[2];
rz(2.6563472) q[3];
sx q[3];
rz(-1.3121366) q[3];
sx q[3];
rz(2.7237039) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8682206) q[0];
sx q[0];
rz(-0.98216787) q[0];
sx q[0];
rz(-2.3305273) q[0];
rz(-1.9048994) q[1];
sx q[1];
rz(-2.2752454) q[1];
sx q[1];
rz(2.9467357) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9962353) q[0];
sx q[0];
rz(-0.35439098) q[0];
sx q[0];
rz(-2.3217391) q[0];
rz(-2.5783875) q[2];
sx q[2];
rz(-0.99894542) q[2];
sx q[2];
rz(0.88746136) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.5681618) q[1];
sx q[1];
rz(-0.57933148) q[1];
sx q[1];
rz(0.90888494) q[1];
rz(-2.5833292) q[3];
sx q[3];
rz(-2.1462064) q[3];
sx q[3];
rz(0.093274506) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.98086944) q[2];
sx q[2];
rz(-1.8941433) q[2];
sx q[2];
rz(-0.54171872) q[2];
rz(1.8187652) q[3];
sx q[3];
rz(-1.3397237) q[3];
sx q[3];
rz(0.26255783) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59170359) q[0];
sx q[0];
rz(-1.6240969) q[0];
sx q[0];
rz(0.24895915) q[0];
rz(-1.9242363) q[1];
sx q[1];
rz(-1.267642) q[1];
sx q[1];
rz(0.41044661) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3623724) q[0];
sx q[0];
rz(-1.1817314) q[0];
sx q[0];
rz(2.7011407) q[0];
rz(-1.6617695) q[2];
sx q[2];
rz(-1.5814648) q[2];
sx q[2];
rz(1.7144698) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.96406781) q[1];
sx q[1];
rz(-0.52492889) q[1];
sx q[1];
rz(2.0711511) q[1];
rz(-0.66907042) q[3];
sx q[3];
rz(-0.72815547) q[3];
sx q[3];
rz(-1.5012342) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.88811389) q[2];
sx q[2];
rz(-1.1682744) q[2];
sx q[2];
rz(2.0237563) q[2];
rz(-3.0228293) q[3];
sx q[3];
rz(-2.4517086) q[3];
sx q[3];
rz(-1.9411055) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(0.73040199) q[0];
sx q[0];
rz(-1.7350736) q[0];
sx q[0];
rz(-0.34179678) q[0];
rz(0.047601184) q[1];
sx q[1];
rz(-1.2536512) q[1];
sx q[1];
rz(1.8258078) q[1];
rz(0.81472266) q[2];
sx q[2];
rz(-1.9152894) q[2];
sx q[2];
rz(0.49283129) q[2];
rz(0.66524617) q[3];
sx q[3];
rz(-0.86626296) q[3];
sx q[3];
rz(-2.5588425) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
