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
rz(0.74066585) q[0];
sx q[0];
rz(1.8494777) q[0];
sx q[0];
rz(10.194095) q[0];
rz(1.5953335) q[1];
sx q[1];
rz(-2.9268664) q[1];
sx q[1];
rz(-0.62155849) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4165933) q[0];
sx q[0];
rz(-2.160797) q[0];
sx q[0];
rz(-0.9585462) q[0];
x q[1];
rz(-1.1530959) q[2];
sx q[2];
rz(-1.5547032) q[2];
sx q[2];
rz(-2.180445) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.5947508) q[1];
sx q[1];
rz(-1.8545281) q[1];
sx q[1];
rz(-0.8304412) q[1];
x q[2];
rz(1.6509455) q[3];
sx q[3];
rz(-1.4372131) q[3];
sx q[3];
rz(2.7770673) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.2810716) q[2];
sx q[2];
rz(-0.29025429) q[2];
sx q[2];
rz(2.2104635) q[2];
rz(-2.9562601) q[3];
sx q[3];
rz(-1.0328181) q[3];
sx q[3];
rz(-1.4144271) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0463878) q[0];
sx q[0];
rz(-0.33400184) q[0];
sx q[0];
rz(2.1686676) q[0];
rz(-2.3880549) q[1];
sx q[1];
rz(-2.2851508) q[1];
sx q[1];
rz(0.10261745) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0282183) q[0];
sx q[0];
rz(-1.4732142) q[0];
sx q[0];
rz(-0.085208864) q[0];
rz(2.5061812) q[2];
sx q[2];
rz(-0.16588587) q[2];
sx q[2];
rz(0.36402853) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.30604369) q[1];
sx q[1];
rz(-1.4144327) q[1];
sx q[1];
rz(-0.083680731) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1903193) q[3];
sx q[3];
rz(-0.35934908) q[3];
sx q[3];
rz(-0.41216601) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.7332581) q[2];
sx q[2];
rz(-1.4380941) q[2];
sx q[2];
rz(0.56780887) q[2];
rz(-2.7693977) q[3];
sx q[3];
rz(-1.3293543) q[3];
sx q[3];
rz(1.6224434) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5837625) q[0];
sx q[0];
rz(-0.74302858) q[0];
sx q[0];
rz(1.2416191) q[0];
rz(-2.325233) q[1];
sx q[1];
rz(-2.7528449) q[1];
sx q[1];
rz(-2.1519318) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96134391) q[0];
sx q[0];
rz(-1.3575866) q[0];
sx q[0];
rz(-0.38589392) q[0];
rz(-0.092421345) q[2];
sx q[2];
rz(-1.4105317) q[2];
sx q[2];
rz(-1.607464) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.951713) q[1];
sx q[1];
rz(-0.44766176) q[1];
sx q[1];
rz(-2.7192804) q[1];
rz(-pi) q[2];
rz(1.9450467) q[3];
sx q[3];
rz(-1.47577) q[3];
sx q[3];
rz(-2.414741) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.376754) q[2];
sx q[2];
rz(-1.9003442) q[2];
sx q[2];
rz(-1.6468916) q[2];
rz(2.3946848) q[3];
sx q[3];
rz(-2.5266095) q[3];
sx q[3];
rz(0.99228215) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(0.84585369) q[0];
sx q[0];
rz(-2.5339412) q[0];
sx q[0];
rz(-0.94957748) q[0];
rz(-1.1664561) q[1];
sx q[1];
rz(-1.8828853) q[1];
sx q[1];
rz(2.9393401) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.14213053) q[0];
sx q[0];
rz(-0.058555257) q[0];
sx q[0];
rz(0.55051319) q[0];
rz(-pi) q[1];
rz(-2.9390088) q[2];
sx q[2];
rz(-0.69949818) q[2];
sx q[2];
rz(-2.0650748) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.7864322) q[1];
sx q[1];
rz(-2.0789685) q[1];
sx q[1];
rz(-0.90759214) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1465591) q[3];
sx q[3];
rz(-2.0821794) q[3];
sx q[3];
rz(0.56727876) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.64968455) q[2];
sx q[2];
rz(-1.4338355) q[2];
sx q[2];
rz(1.3362308) q[2];
rz(2.6797471) q[3];
sx q[3];
rz(-1.0577842) q[3];
sx q[3];
rz(1.0285146) q[3];
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
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.860054) q[0];
sx q[0];
rz(-3.0335732) q[0];
sx q[0];
rz(-2.7463013) q[0];
rz(-2.1832502) q[1];
sx q[1];
rz(-1.3805362) q[1];
sx q[1];
rz(-2.7470632) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1720264) q[0];
sx q[0];
rz(-2.8752669) q[0];
sx q[0];
rz(-1.6380493) q[0];
x q[1];
rz(-0.5455234) q[2];
sx q[2];
rz(-0.50456556) q[2];
sx q[2];
rz(-1.2009753) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.3580618) q[1];
sx q[1];
rz(-0.79800843) q[1];
sx q[1];
rz(-1.111387) q[1];
rz(-pi) q[2];
x q[2];
rz(1.35704) q[3];
sx q[3];
rz(-0.96352623) q[3];
sx q[3];
rz(2.7979545) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.40256527) q[2];
sx q[2];
rz(-1.6141011) q[2];
sx q[2];
rz(-1.019545) q[2];
rz(-1.3022276) q[3];
sx q[3];
rz(-3.0890833) q[3];
sx q[3];
rz(-2.4359865) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(2.5604945) q[0];
sx q[0];
rz(-0.56466931) q[0];
sx q[0];
rz(-1.0127006) q[0];
rz(-0.46214354) q[1];
sx q[1];
rz(-1.4116762) q[1];
sx q[1];
rz(0.046646811) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3585535) q[0];
sx q[0];
rz(-1.5601242) q[0];
sx q[0];
rz(1.4477167) q[0];
rz(2.279206) q[2];
sx q[2];
rz(-0.38225952) q[2];
sx q[2];
rz(1.6406989) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.44429) q[1];
sx q[1];
rz(-1.9026766) q[1];
sx q[1];
rz(2.1730971) q[1];
rz(-pi) q[2];
rz(-1.1936452) q[3];
sx q[3];
rz(-1.5932249) q[3];
sx q[3];
rz(0.42322657) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.3095653) q[2];
sx q[2];
rz(-0.088531606) q[2];
sx q[2];
rz(0.95477611) q[2];
rz(0.66983062) q[3];
sx q[3];
rz(-1.1845183) q[3];
sx q[3];
rz(-0.36439103) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9949812) q[0];
sx q[0];
rz(-2.8611188) q[0];
sx q[0];
rz(-1.9914419) q[0];
rz(-3.0448044) q[1];
sx q[1];
rz(-2.0014747) q[1];
sx q[1];
rz(1.6614301) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.87931765) q[0];
sx q[0];
rz(-1.1918544) q[0];
sx q[0];
rz(2.3728601) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3315013) q[2];
sx q[2];
rz(-1.8836438) q[2];
sx q[2];
rz(-0.32101908) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.11058775) q[1];
sx q[1];
rz(-1.8613986) q[1];
sx q[1];
rz(-0.1433934) q[1];
x q[2];
rz(2.0050762) q[3];
sx q[3];
rz(-1.7107114) q[3];
sx q[3];
rz(2.851448) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4074771) q[2];
sx q[2];
rz(-2.263767) q[2];
sx q[2];
rz(2.5131098) q[2];
rz(2.8001522) q[3];
sx q[3];
rz(-2.1554558) q[3];
sx q[3];
rz(-0.83610523) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2899365) q[0];
sx q[0];
rz(-0.66806942) q[0];
sx q[0];
rz(-3.0776899) q[0];
rz(0.54197657) q[1];
sx q[1];
rz(-2.0109476) q[1];
sx q[1];
rz(-2.7625387) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0444894) q[0];
sx q[0];
rz(-1.7428723) q[0];
sx q[0];
rz(1.8015566) q[0];
rz(-1.7070243) q[2];
sx q[2];
rz(-0.76644015) q[2];
sx q[2];
rz(2.5640783) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.29330373) q[1];
sx q[1];
rz(-2.1142936) q[1];
sx q[1];
rz(-1.5889542) q[1];
x q[2];
rz(0.97207119) q[3];
sx q[3];
rz(-2.4643371) q[3];
sx q[3];
rz(-0.87358356) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.1698251) q[2];
sx q[2];
rz(-2.0242033) q[2];
sx q[2];
rz(0.28905147) q[2];
rz(-1.0158094) q[3];
sx q[3];
rz(-0.93520516) q[3];
sx q[3];
rz(-2.9142006) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10580258) q[0];
sx q[0];
rz(-2.7547024) q[0];
sx q[0];
rz(-2.8111358) q[0];
rz(-0.48078787) q[1];
sx q[1];
rz(-1.5590706) q[1];
sx q[1];
rz(-0.50723433) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.98927724) q[0];
sx q[0];
rz(-1.0535272) q[0];
sx q[0];
rz(0.74801561) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.66659303) q[2];
sx q[2];
rz(-2.5799516) q[2];
sx q[2];
rz(-0.14065841) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.54785905) q[1];
sx q[1];
rz(-2.6630109) q[1];
sx q[1];
rz(-1.0793769) q[1];
x q[2];
rz(1.6303924) q[3];
sx q[3];
rz(-0.97138849) q[3];
sx q[3];
rz(2.5580278) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.5759739) q[2];
sx q[2];
rz(-2.2810292) q[2];
sx q[2];
rz(-0.38798517) q[2];
rz(-0.062813736) q[3];
sx q[3];
rz(-2.4020436) q[3];
sx q[3];
rz(-1.1886103) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
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
rz(-1.0509309) q[0];
sx q[0];
rz(-3.1212786) q[0];
sx q[0];
rz(-3.0860197) q[0];
rz(-2.9203501) q[1];
sx q[1];
rz(-2.1317) q[1];
sx q[1];
rz(-1.6411068) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9433946) q[0];
sx q[0];
rz(-1.0695116) q[0];
sx q[0];
rz(-2.8265165) q[0];
x q[1];
rz(3.0048641) q[2];
sx q[2];
rz(-1.6199582) q[2];
sx q[2];
rz(1.8135742) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.21162578) q[1];
sx q[1];
rz(-1.9504419) q[1];
sx q[1];
rz(-0.34348653) q[1];
rz(-pi) q[2];
x q[2];
rz(2.290912) q[3];
sx q[3];
rz(-1.097957) q[3];
sx q[3];
rz(-2.8306876) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.9556433) q[2];
sx q[2];
rz(-0.187749) q[2];
sx q[2];
rz(-2.032568) q[2];
rz(3.1254613) q[3];
sx q[3];
rz(-1.4343836) q[3];
sx q[3];
rz(2.6764638) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
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
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3817417) q[0];
sx q[0];
rz(-2.1325337) q[0];
sx q[0];
rz(2.7251563) q[0];
rz(0.79553678) q[1];
sx q[1];
rz(-1.3858613) q[1];
sx q[1];
rz(-0.081079986) q[1];
rz(-3.0007985) q[2];
sx q[2];
rz(-1.2296321) q[2];
sx q[2];
rz(2.932569) q[2];
rz(-2.6017009) q[3];
sx q[3];
rz(-0.85689982) q[3];
sx q[3];
rz(0.11107437) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
