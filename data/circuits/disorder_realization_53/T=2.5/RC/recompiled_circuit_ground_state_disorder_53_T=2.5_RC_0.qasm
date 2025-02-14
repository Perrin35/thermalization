OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.1528435) q[0];
sx q[0];
rz(3.710521) q[0];
sx q[0];
rz(13.518128) q[0];
rz(3.9495502) q[1];
sx q[1];
rz(6.1679975) q[1];
sx q[1];
rz(8.4313784) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0885926) q[0];
sx q[0];
rz(-1.3198084) q[0];
sx q[0];
rz(-0.95337501) q[0];
rz(-pi) q[1];
rz(-0.98007085) q[2];
sx q[2];
rz(-0.23071846) q[2];
sx q[2];
rz(1.427111) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.6019878) q[1];
sx q[1];
rz(-1.0430416) q[1];
sx q[1];
rz(-0.62202203) q[1];
rz(-1.7073329) q[3];
sx q[3];
rz(-2.2280558) q[3];
sx q[3];
rz(3.0252489) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.50977612) q[2];
sx q[2];
rz(-0.6002554) q[2];
sx q[2];
rz(2.996345) q[2];
rz(2.1692569) q[3];
sx q[3];
rz(-1.626868) q[3];
sx q[3];
rz(-1.2783031) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-1.7250605) q[0];
sx q[0];
rz(-2.5711377) q[0];
sx q[0];
rz(-2.3544627) q[0];
rz(-2.9540673) q[1];
sx q[1];
rz(-1.7555321) q[1];
sx q[1];
rz(-3.131391) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.37885168) q[0];
sx q[0];
rz(-1.6798899) q[0];
sx q[0];
rz(-0.0082882546) q[0];
rz(1.446063) q[2];
sx q[2];
rz(-0.48791781) q[2];
sx q[2];
rz(-3.0782521) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.9157313) q[1];
sx q[1];
rz(-0.26090947) q[1];
sx q[1];
rz(0.4930851) q[1];
rz(-0.83008978) q[3];
sx q[3];
rz(-1.4674868) q[3];
sx q[3];
rz(1.938397) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.06112222) q[2];
sx q[2];
rz(-0.11961131) q[2];
sx q[2];
rz(-1.97869) q[2];
rz(1.2969147) q[3];
sx q[3];
rz(-1.4328522) q[3];
sx q[3];
rz(0.0059303693) q[3];
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
rz(1.6892683) q[0];
sx q[0];
rz(-2.5732714) q[0];
sx q[0];
rz(-2.7962621) q[0];
rz(-3.0863702) q[1];
sx q[1];
rz(-0.60631141) q[1];
sx q[1];
rz(-0.20923722) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.064406618) q[0];
sx q[0];
rz(-1.4923054) q[0];
sx q[0];
rz(-0.85759832) q[0];
x q[1];
rz(-1.8183858) q[2];
sx q[2];
rz(-0.074562975) q[2];
sx q[2];
rz(-2.8207795) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.0303231) q[1];
sx q[1];
rz(-1.7220338) q[1];
sx q[1];
rz(-2.4200685) q[1];
rz(-pi) q[2];
rz(2.8992527) q[3];
sx q[3];
rz(-2.3018357) q[3];
sx q[3];
rz(3.1104607) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.0914586) q[2];
sx q[2];
rz(-1.4860934) q[2];
sx q[2];
rz(0.44198188) q[2];
rz(-0.4778536) q[3];
sx q[3];
rz(-1.7171532) q[3];
sx q[3];
rz(1.8817687) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0100937) q[0];
sx q[0];
rz(-2.5690014) q[0];
sx q[0];
rz(1.2226489) q[0];
rz(-1.6652416) q[1];
sx q[1];
rz(-2.1189225) q[1];
sx q[1];
rz(0.30278912) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0593392) q[0];
sx q[0];
rz(-0.48294623) q[0];
sx q[0];
rz(2.5046136) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1533621) q[2];
sx q[2];
rz(-0.086518651) q[2];
sx q[2];
rz(0.38245538) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.5975133) q[1];
sx q[1];
rz(-1.5248393) q[1];
sx q[1];
rz(1.2567029) q[1];
rz(0.92118951) q[3];
sx q[3];
rz(-2.0944893) q[3];
sx q[3];
rz(2.7648448) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.0629695) q[2];
sx q[2];
rz(-1.4457694) q[2];
sx q[2];
rz(-1.3485738) q[2];
rz(3.1283227) q[3];
sx q[3];
rz(-0.66418663) q[3];
sx q[3];
rz(-1.9224904) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.211798) q[0];
sx q[0];
rz(-3.0199265) q[0];
sx q[0];
rz(1.1965055) q[0];
rz(-2.3579146) q[1];
sx q[1];
rz(-1.0107026) q[1];
sx q[1];
rz(0.7799305) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9570628) q[0];
sx q[0];
rz(-1.3100782) q[0];
sx q[0];
rz(-2.3924559) q[0];
rz(-3.0230396) q[2];
sx q[2];
rz(-2.701512) q[2];
sx q[2];
rz(-1.5660182) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(3.0454516) q[1];
sx q[1];
rz(-2.2287031) q[1];
sx q[1];
rz(1.0784574) q[1];
rz(-pi) q[2];
rz(-2.6013589) q[3];
sx q[3];
rz(-1.295394) q[3];
sx q[3];
rz(1.6364975) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.1448867) q[2];
sx q[2];
rz(-0.11652623) q[2];
sx q[2];
rz(3.0987926) q[2];
rz(-1.5329817) q[3];
sx q[3];
rz(-1.7973085) q[3];
sx q[3];
rz(-0.90095055) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
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
rz(-0.75339371) q[0];
sx q[0];
rz(-0.94069427) q[0];
sx q[0];
rz(0.70145506) q[0];
rz(-0.88297168) q[1];
sx q[1];
rz(-1.3621623) q[1];
sx q[1];
rz(-2.0515474) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0018679) q[0];
sx q[0];
rz(-1.9650978) q[0];
sx q[0];
rz(0.35625881) q[0];
rz(-pi) q[1];
rz(-0.50925635) q[2];
sx q[2];
rz(-2.351077) q[2];
sx q[2];
rz(-2.5495286) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.54451128) q[1];
sx q[1];
rz(-1.4284572) q[1];
sx q[1];
rz(-0.049171731) q[1];
x q[2];
rz(-1.8102856) q[3];
sx q[3];
rz(-1.4109932) q[3];
sx q[3];
rz(-1.0996475) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.8604454) q[2];
sx q[2];
rz(-0.88938418) q[2];
sx q[2];
rz(-1.3775728) q[2];
rz(0.095976202) q[3];
sx q[3];
rz(-2.4317661) q[3];
sx q[3];
rz(2.6095552) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
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
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1808566) q[0];
sx q[0];
rz(-2.2611698) q[0];
sx q[0];
rz(-1.9218504) q[0];
rz(2.5324054) q[1];
sx q[1];
rz(-1.6222619) q[1];
sx q[1];
rz(2.6763197) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1022961) q[0];
sx q[0];
rz(-1.0900153) q[0];
sx q[0];
rz(-2.395438) q[0];
rz(2.005079) q[2];
sx q[2];
rz(-0.65616504) q[2];
sx q[2];
rz(-0.16193552) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.45256685) q[1];
sx q[1];
rz(-1.5264319) q[1];
sx q[1];
rz(1.1242799) q[1];
x q[2];
rz(0.19962387) q[3];
sx q[3];
rz(-2.3559442) q[3];
sx q[3];
rz(0.35405418) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.6294127) q[2];
sx q[2];
rz(-2.8359154) q[2];
sx q[2];
rz(0.42210397) q[2];
rz(3.1230538) q[3];
sx q[3];
rz(-1.5471349) q[3];
sx q[3];
rz(2.5277933) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2043532) q[0];
sx q[0];
rz(-2.6528907) q[0];
sx q[0];
rz(-2.3204284) q[0];
rz(-0.69860727) q[1];
sx q[1];
rz(-2.3804074) q[1];
sx q[1];
rz(3.1331114) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9152731) q[0];
sx q[0];
rz(-2.2157359) q[0];
sx q[0];
rz(-1.9553095) q[0];
x q[1];
rz(0.39706612) q[2];
sx q[2];
rz(-1.0806568) q[2];
sx q[2];
rz(-1.9352143) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.20336452) q[1];
sx q[1];
rz(-1.3044323) q[1];
sx q[1];
rz(-2.3530234) q[1];
x q[2];
rz(0.92681411) q[3];
sx q[3];
rz(-1.5379526) q[3];
sx q[3];
rz(0.72073267) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.2681793) q[2];
sx q[2];
rz(-3.029533) q[2];
sx q[2];
rz(-2.6123135) q[2];
rz(-1.4389634) q[3];
sx q[3];
rz(-1.8301423) q[3];
sx q[3];
rz(-2.8992991) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
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
rz(2.5731803) q[0];
sx q[0];
rz(-0.72887623) q[0];
sx q[0];
rz(-2.6045784) q[0];
rz(0.88045398) q[1];
sx q[1];
rz(-1.8388137) q[1];
sx q[1];
rz(2.3939078) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0466008) q[0];
sx q[0];
rz(-1.23833) q[0];
sx q[0];
rz(2.4103122) q[0];
x q[1];
rz(2.7319261) q[2];
sx q[2];
rz(-0.94575277) q[2];
sx q[2];
rz(-1.1898439) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.6612271) q[1];
sx q[1];
rz(-1.2496619) q[1];
sx q[1];
rz(-2.8868746) q[1];
rz(1.8219833) q[3];
sx q[3];
rz(-0.97339857) q[3];
sx q[3];
rz(0.26439127) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.487454) q[2];
sx q[2];
rz(-2.4345001) q[2];
sx q[2];
rz(0.96997619) q[2];
rz(-1.6449432) q[3];
sx q[3];
rz(-2.1275438) q[3];
sx q[3];
rz(-1.0223201) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3269761) q[0];
sx q[0];
rz(-1.6642445) q[0];
sx q[0];
rz(0.60272637) q[0];
rz(0.0043491443) q[1];
sx q[1];
rz(-1.3163722) q[1];
sx q[1];
rz(2.0306921) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0051188) q[0];
sx q[0];
rz(-0.92347758) q[0];
sx q[0];
rz(0.74764436) q[0];
x q[1];
rz(-0.87238042) q[2];
sx q[2];
rz(-1.9920298) q[2];
sx q[2];
rz(0.55544686) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.7676684) q[1];
sx q[1];
rz(-1.6925188) q[1];
sx q[1];
rz(2.5472104) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6708665) q[3];
sx q[3];
rz(-1.8672872) q[3];
sx q[3];
rz(-2.61946) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.5000308) q[2];
sx q[2];
rz(-1.2189453) q[2];
sx q[2];
rz(-2.5489589) q[2];
rz(2.5617013) q[3];
sx q[3];
rz(-0.65120828) q[3];
sx q[3];
rz(-1.7013223) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.92689571) q[0];
sx q[0];
rz(-1.0746645) q[0];
sx q[0];
rz(-0.97434531) q[0];
rz(-0.018085619) q[1];
sx q[1];
rz(-0.34307243) q[1];
sx q[1];
rz(-3.051563) q[1];
rz(-1.8263578) q[2];
sx q[2];
rz(-1.2719874) q[2];
sx q[2];
rz(-0.86612305) q[2];
rz(-0.52094372) q[3];
sx q[3];
rz(-0.81650618) q[3];
sx q[3];
rz(3.0405844) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
