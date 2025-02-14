OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.0889283) q[0];
sx q[0];
rz(5.1533617) q[0];
sx q[0];
rz(12.552153) q[0];
rz(0.05834236) q[1];
sx q[1];
rz(-2.3925233) q[1];
sx q[1];
rz(-2.5383389) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39122836) q[0];
sx q[0];
rz(-2.3846013) q[0];
sx q[0];
rz(-0.52387497) q[0];
rz(-pi) q[1];
rz(-2.2255664) q[2];
sx q[2];
rz(-1.5226622) q[2];
sx q[2];
rz(2.2544238) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.1849564) q[1];
sx q[1];
rz(-2.2549964) q[1];
sx q[1];
rz(-1.7723899) q[1];
rz(-pi) q[2];
rz(0.18666191) q[3];
sx q[3];
rz(-2.2687654) q[3];
sx q[3];
rz(-1.9912793) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.3732036) q[2];
sx q[2];
rz(-1.6037805) q[2];
sx q[2];
rz(-1.0962037) q[2];
rz(-1.0783892) q[3];
sx q[3];
rz(-2.6069141) q[3];
sx q[3];
rz(2.6064742) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.32644367) q[0];
sx q[0];
rz(-1.0468227) q[0];
sx q[0];
rz(0.4775508) q[0];
rz(2.1304456) q[1];
sx q[1];
rz(-2.5944581) q[1];
sx q[1];
rz(-1.0844885) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3352184) q[0];
sx q[0];
rz(-1.9382297) q[0];
sx q[0];
rz(-0.12560496) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6376375) q[2];
sx q[2];
rz(-2.7324841) q[2];
sx q[2];
rz(0.44873777) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.4762759) q[1];
sx q[1];
rz(-1.3937104) q[1];
sx q[1];
rz(0.55426532) q[1];
rz(-2.4020477) q[3];
sx q[3];
rz(-1.7266183) q[3];
sx q[3];
rz(-2.9062944) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.5227205) q[2];
sx q[2];
rz(-2.7174157) q[2];
sx q[2];
rz(-2.0514533) q[2];
rz(2.5777396) q[3];
sx q[3];
rz(-2.4520051) q[3];
sx q[3];
rz(3.1088945) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76145935) q[0];
sx q[0];
rz(-0.021012336) q[0];
sx q[0];
rz(1.5367966) q[0];
rz(-1.8236632) q[1];
sx q[1];
rz(-1.0937966) q[1];
sx q[1];
rz(-0.39753786) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1434162) q[0];
sx q[0];
rz(-2.2097124) q[0];
sx q[0];
rz(2.8692416) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0084559) q[2];
sx q[2];
rz(-2.5175885) q[2];
sx q[2];
rz(-0.74988264) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.2469625) q[1];
sx q[1];
rz(-1.8254465) q[1];
sx q[1];
rz(-0.92856426) q[1];
rz(2.378355) q[3];
sx q[3];
rz(-0.67585603) q[3];
sx q[3];
rz(1.8873029) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.2354108) q[2];
sx q[2];
rz(-1.9637039) q[2];
sx q[2];
rz(0.57514352) q[2];
rz(-0.67342526) q[3];
sx q[3];
rz(-1.6005102) q[3];
sx q[3];
rz(1.5967691) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6509318) q[0];
sx q[0];
rz(-2.0199825) q[0];
sx q[0];
rz(-2.2326873) q[0];
rz(3.124369) q[1];
sx q[1];
rz(-2.6135018) q[1];
sx q[1];
rz(-0.012103279) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72649137) q[0];
sx q[0];
rz(-2.4612777) q[0];
sx q[0];
rz(-0.52796396) q[0];
x q[1];
rz(-1.034017) q[2];
sx q[2];
rz(-0.99269789) q[2];
sx q[2];
rz(0.54852099) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.7722655) q[1];
sx q[1];
rz(-2.4924534) q[1];
sx q[1];
rz(1.9499798) q[1];
rz(-0.34194591) q[3];
sx q[3];
rz(-3.0071335) q[3];
sx q[3];
rz(-0.97245261) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.072731344) q[2];
sx q[2];
rz(-0.27421811) q[2];
sx q[2];
rz(-2.7552628) q[2];
rz(-0.38935152) q[3];
sx q[3];
rz(-1.6081622) q[3];
sx q[3];
rz(-2.1336335) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7858081) q[0];
sx q[0];
rz(-2.4171827) q[0];
sx q[0];
rz(-0.32387787) q[0];
rz(2.7549699) q[1];
sx q[1];
rz(-1.7522248) q[1];
sx q[1];
rz(1.5868384) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1607035) q[0];
sx q[0];
rz(-1.5063154) q[0];
sx q[0];
rz(-2.4095201) q[0];
x q[1];
rz(2.1065169) q[2];
sx q[2];
rz(-1.2457827) q[2];
sx q[2];
rz(-0.49301389) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.28858057) q[1];
sx q[1];
rz(-1.1094571) q[1];
sx q[1];
rz(0.71801676) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3426476) q[3];
sx q[3];
rz(-1.9513592) q[3];
sx q[3];
rz(2.9137672) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.7524449) q[2];
sx q[2];
rz(-0.84263313) q[2];
sx q[2];
rz(-3.0774934) q[2];
rz(2.5373503) q[3];
sx q[3];
rz(-1.7490381) q[3];
sx q[3];
rz(-1.2324246) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9391249) q[0];
sx q[0];
rz(-2.5305643) q[0];
sx q[0];
rz(0.87919277) q[0];
rz(2.7912256) q[1];
sx q[1];
rz(-1.3060952) q[1];
sx q[1];
rz(2.9894357) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8504753) q[0];
sx q[0];
rz(-0.88248173) q[0];
sx q[0];
rz(1.9242213) q[0];
rz(-2.9486676) q[2];
sx q[2];
rz(-2.0461651) q[2];
sx q[2];
rz(1.0032636) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.8250936) q[1];
sx q[1];
rz(-0.48227019) q[1];
sx q[1];
rz(2.6939166) q[1];
rz(1.1760487) q[3];
sx q[3];
rz(-1.1541894) q[3];
sx q[3];
rz(-0.28232251) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.0396314) q[2];
sx q[2];
rz(-2.8122718) q[2];
sx q[2];
rz(-0.24001089) q[2];
rz(-2.4890684) q[3];
sx q[3];
rz(-1.1381166) q[3];
sx q[3];
rz(-1.8664546) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0535102) q[0];
sx q[0];
rz(-1.2367915) q[0];
sx q[0];
rz(-0.85813338) q[0];
rz(-1.7543322) q[1];
sx q[1];
rz(-0.45448449) q[1];
sx q[1];
rz(2.8087356) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.34073601) q[0];
sx q[0];
rz(-1.8854257) q[0];
sx q[0];
rz(0.97140177) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5120686) q[2];
sx q[2];
rz(-1.8203805) q[2];
sx q[2];
rz(1.6317612) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.8695527) q[1];
sx q[1];
rz(-1.7764047) q[1];
sx q[1];
rz(-2.7133792) q[1];
rz(-pi) q[2];
x q[2];
rz(0.44705963) q[3];
sx q[3];
rz(-0.70526988) q[3];
sx q[3];
rz(-1.3652319) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.7994999) q[2];
sx q[2];
rz(-0.19793333) q[2];
sx q[2];
rz(1.0510772) q[2];
rz(1.6445232) q[3];
sx q[3];
rz(-1.1727419) q[3];
sx q[3];
rz(-0.79819775) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
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
rz(0.77338162) q[0];
sx q[0];
rz(-2.1147275) q[0];
sx q[0];
rz(-0.89212242) q[0];
rz(0.46961531) q[1];
sx q[1];
rz(-0.62925595) q[1];
sx q[1];
rz(-2.3687252) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5443916) q[0];
sx q[0];
rz(-2.3344451) q[0];
sx q[0];
rz(-2.3788484) q[0];
rz(-pi) q[1];
rz(0.30420423) q[2];
sx q[2];
rz(-0.72663621) q[2];
sx q[2];
rz(1.8223423) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.6438021) q[1];
sx q[1];
rz(-2.0333159) q[1];
sx q[1];
rz(-1.2116572) q[1];
rz(0.77247932) q[3];
sx q[3];
rz(-2.8101343) q[3];
sx q[3];
rz(0.14642388) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.5272687) q[2];
sx q[2];
rz(-1.9169151) q[2];
sx q[2];
rz(-0.69250715) q[2];
rz(2.6912189) q[3];
sx q[3];
rz(-1.9362484) q[3];
sx q[3];
rz(-1.6374121) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5699128) q[0];
sx q[0];
rz(-0.79140651) q[0];
sx q[0];
rz(-2.1929542) q[0];
rz(1.0665077) q[1];
sx q[1];
rz(-2.152161) q[1];
sx q[1];
rz(1.271064) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1656837) q[0];
sx q[0];
rz(-1.4780439) q[0];
sx q[0];
rz(-1.776987) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.022106013) q[2];
sx q[2];
rz(-1.5367295) q[2];
sx q[2];
rz(0.23512041) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.5258057) q[1];
sx q[1];
rz(-2.5005625) q[1];
sx q[1];
rz(1.7980144) q[1];
rz(2.8715114) q[3];
sx q[3];
rz(-2.684786) q[3];
sx q[3];
rz(-1.0246667) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.54242414) q[2];
sx q[2];
rz(-1.8343265) q[2];
sx q[2];
rz(-3.001281) q[2];
rz(-3.1091651) q[3];
sx q[3];
rz(-0.92493886) q[3];
sx q[3];
rz(1.2062581) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5133544) q[0];
sx q[0];
rz(-1.4435377) q[0];
sx q[0];
rz(-2.3640609) q[0];
rz(-2.2041722) q[1];
sx q[1];
rz(-1.9098858) q[1];
sx q[1];
rz(-0.44313988) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69618689) q[0];
sx q[0];
rz(-1.6517795) q[0];
sx q[0];
rz(3.0565673) q[0];
rz(-pi) q[1];
rz(-1.1045157) q[2];
sx q[2];
rz(-0.65090042) q[2];
sx q[2];
rz(-2.4492531) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.7305371) q[1];
sx q[1];
rz(-1.2403508) q[1];
sx q[1];
rz(2.4871268) q[1];
rz(-2.0432161) q[3];
sx q[3];
rz(-2.7948423) q[3];
sx q[3];
rz(-0.43671331) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0810658) q[2];
sx q[2];
rz(-2.492283) q[2];
sx q[2];
rz(-3.013986) q[2];
rz(-2.8429032) q[3];
sx q[3];
rz(-1.9689711) q[3];
sx q[3];
rz(-2.7711788) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.61113197) q[0];
sx q[0];
rz(-1.8831384) q[0];
sx q[0];
rz(0.25181121) q[0];
rz(2.5313189) q[1];
sx q[1];
rz(-0.88519575) q[1];
sx q[1];
rz(-2.0559678) q[1];
rz(1.4971785) q[2];
sx q[2];
rz(-2.4215019) q[2];
sx q[2];
rz(-0.99046594) q[2];
rz(2.0213303) q[3];
sx q[3];
rz(-1.7623175) q[3];
sx q[3];
rz(0.13468066) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
