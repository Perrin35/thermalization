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
rz(2.2284722) q[0];
sx q[0];
rz(-2.1051536) q[0];
sx q[0];
rz(1.5358465) q[0];
rz(-2.4802471) q[1];
sx q[1];
rz(-1.9547434) q[1];
sx q[1];
rz(-0.43651906) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.75398895) q[0];
sx q[0];
rz(-1.6696829) q[0];
sx q[0];
rz(-1.2515227) q[0];
rz(1.1182734) q[2];
sx q[2];
rz(-2.088639) q[2];
sx q[2];
rz(2.81942) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.4492717) q[1];
sx q[1];
rz(-0.8682478) q[1];
sx q[1];
rz(3.0880804) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7868452) q[3];
sx q[3];
rz(-2.4461544) q[3];
sx q[3];
rz(-1.2176633) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.46301699) q[2];
sx q[2];
rz(-2.7028694) q[2];
sx q[2];
rz(0.92301816) q[2];
rz(0.065741278) q[3];
sx q[3];
rz(-1.3275423) q[3];
sx q[3];
rz(-1.340723) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.0098669212) q[0];
sx q[0];
rz(-2.0781524) q[0];
sx q[0];
rz(-3.104082) q[0];
rz(2.5610979) q[1];
sx q[1];
rz(-2.4721804) q[1];
sx q[1];
rz(2.4494749) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1624533) q[0];
sx q[0];
rz(-1.7458436) q[0];
sx q[0];
rz(-1.2412423) q[0];
x q[1];
rz(0.28022556) q[2];
sx q[2];
rz(-1.0839274) q[2];
sx q[2];
rz(3.1161199) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-3.1375274) q[1];
sx q[1];
rz(-2.1430227) q[1];
sx q[1];
rz(-2.1768084) q[1];
rz(-2.9565937) q[3];
sx q[3];
rz(-0.76055148) q[3];
sx q[3];
rz(-0.669125) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.47824255) q[2];
sx q[2];
rz(-1.1967836) q[2];
sx q[2];
rz(1.3956068) q[2];
rz(0.1869525) q[3];
sx q[3];
rz(-1.170265) q[3];
sx q[3];
rz(0.54005867) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66889399) q[0];
sx q[0];
rz(-2.5085594) q[0];
sx q[0];
rz(-2.582666) q[0];
rz(-0.82866296) q[1];
sx q[1];
rz(-1.5507973) q[1];
sx q[1];
rz(-0.25159803) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.64564451) q[0];
sx q[0];
rz(-2.6584266) q[0];
sx q[0];
rz(0.42074411) q[0];
rz(-pi) q[1];
rz(1.0195658) q[2];
sx q[2];
rz(-1.8675659) q[2];
sx q[2];
rz(3.0574329) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.025574) q[1];
sx q[1];
rz(-1.3085828) q[1];
sx q[1];
rz(2.2518454) q[1];
rz(0.91355998) q[3];
sx q[3];
rz(-2.1373478) q[3];
sx q[3];
rz(-2.4853137) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.33637968) q[2];
sx q[2];
rz(-0.42625913) q[2];
sx q[2];
rz(1.7097998) q[2];
rz(0.21670565) q[3];
sx q[3];
rz(-1.6101937) q[3];
sx q[3];
rz(1.7823035) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25877229) q[0];
sx q[0];
rz(-2.6348305) q[0];
sx q[0];
rz(-1.2910507) q[0];
rz(-2.3123815) q[1];
sx q[1];
rz(-1.1081089) q[1];
sx q[1];
rz(-2.1270027) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0907959) q[0];
sx q[0];
rz(-1.8136588) q[0];
sx q[0];
rz(-2.5953351) q[0];
x q[1];
rz(-1.9564187) q[2];
sx q[2];
rz(-2.1004538) q[2];
sx q[2];
rz(-0.5822863) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.7493128) q[1];
sx q[1];
rz(-1.5661513) q[1];
sx q[1];
rz(-0.99601142) q[1];
rz(-1.7144373) q[3];
sx q[3];
rz(-0.82089409) q[3];
sx q[3];
rz(-0.088975541) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.1318876) q[2];
sx q[2];
rz(-1.9405126) q[2];
sx q[2];
rz(-2.3809872) q[2];
rz(-0.46918121) q[3];
sx q[3];
rz(-1.4696591) q[3];
sx q[3];
rz(-0.73452264) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0932662) q[0];
sx q[0];
rz(-0.51437298) q[0];
sx q[0];
rz(0.58116466) q[0];
rz(-1.5515074) q[1];
sx q[1];
rz(-0.64684144) q[1];
sx q[1];
rz(2.4466628) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0804176) q[0];
sx q[0];
rz(-1.8616849) q[0];
sx q[0];
rz(-2.8982452) q[0];
rz(-pi) q[1];
rz(-2.8492691) q[2];
sx q[2];
rz(-1.6609123) q[2];
sx q[2];
rz(0.19168774) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(3.123053) q[1];
sx q[1];
rz(-0.54529069) q[1];
sx q[1];
rz(0.084597708) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.23628037) q[3];
sx q[3];
rz(-1.3930916) q[3];
sx q[3];
rz(-1.9654993) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.27802262) q[2];
sx q[2];
rz(-1.6081955) q[2];
sx q[2];
rz(-0.30280534) q[2];
rz(-1.4263724) q[3];
sx q[3];
rz(-0.36077603) q[3];
sx q[3];
rz(-0.58437955) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.67702174) q[0];
sx q[0];
rz(-0.40957054) q[0];
sx q[0];
rz(-2.4466178) q[0];
rz(-0.28680828) q[1];
sx q[1];
rz(-2.5359055) q[1];
sx q[1];
rz(1.8668176) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5329842) q[0];
sx q[0];
rz(-1.8814981) q[0];
sx q[0];
rz(-1.0644738) q[0];
rz(-pi) q[1];
x q[1];
rz(0.54194258) q[2];
sx q[2];
rz(-1.9795609) q[2];
sx q[2];
rz(2.8084076) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.4732237) q[1];
sx q[1];
rz(-2.372207) q[1];
sx q[1];
rz(2.0955057) q[1];
rz(2.7872179) q[3];
sx q[3];
rz(-0.61180002) q[3];
sx q[3];
rz(-2.967088) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.0810762) q[2];
sx q[2];
rz(-2.923107) q[2];
sx q[2];
rz(-0.9673574) q[2];
rz(-0.71990144) q[3];
sx q[3];
rz(-0.94341174) q[3];
sx q[3];
rz(1.2758183) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(0.28441456) q[0];
sx q[0];
rz(-3.1107749) q[0];
sx q[0];
rz(1.6946633) q[0];
rz(0.46724304) q[1];
sx q[1];
rz(-1.2930361) q[1];
sx q[1];
rz(1.9460868) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0341114) q[0];
sx q[0];
rz(-1.1426289) q[0];
sx q[0];
rz(2.0780245) q[0];
x q[1];
rz(1.906997) q[2];
sx q[2];
rz(-0.78374353) q[2];
sx q[2];
rz(-2.9072493) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.171716) q[1];
sx q[1];
rz(-2.6387353) q[1];
sx q[1];
rz(1.4420532) q[1];
rz(-pi) q[2];
rz(1.9105114) q[3];
sx q[3];
rz(-0.038148316) q[3];
sx q[3];
rz(-1.3919786) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.82470992) q[2];
sx q[2];
rz(-0.83493817) q[2];
sx q[2];
rz(0.52428025) q[2];
rz(-2.8892062) q[3];
sx q[3];
rz(-3.0350244) q[3];
sx q[3];
rz(3.0805123) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1687932) q[0];
sx q[0];
rz(-1.1058829) q[0];
sx q[0];
rz(1.9148781) q[0];
rz(-2.3854158) q[1];
sx q[1];
rz(-2.1935479) q[1];
sx q[1];
rz(-2.7511168) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7602063) q[0];
sx q[0];
rz(-1.1293656) q[0];
sx q[0];
rz(1.2797829) q[0];
rz(1.5323213) q[2];
sx q[2];
rz(-2.0145825) q[2];
sx q[2];
rz(-2.2930068) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.272534) q[1];
sx q[1];
rz(-2.7400937) q[1];
sx q[1];
rz(0.99799207) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.94589767) q[3];
sx q[3];
rz(-2.3970553) q[3];
sx q[3];
rz(0.27305279) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.8172424) q[2];
sx q[2];
rz(-0.70049006) q[2];
sx q[2];
rz(1.3520757) q[2];
rz(2.9663441) q[3];
sx q[3];
rz(-1.8027571) q[3];
sx q[3];
rz(1.2043183) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[3];
rz(-pi/2) q[3];
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
rz(-1.4125724) q[0];
sx q[0];
rz(-2.6798798) q[0];
sx q[0];
rz(-2.4793258) q[0];
rz(-2.5114255) q[1];
sx q[1];
rz(-1.2799193) q[1];
sx q[1];
rz(-0.23552775) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6306289) q[0];
sx q[0];
rz(-2.8665344) q[0];
sx q[0];
rz(2.4514626) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6643737) q[2];
sx q[2];
rz(-2.1380599) q[2];
sx q[2];
rz(2.7964489) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.3709066) q[1];
sx q[1];
rz(-1.5471493) q[1];
sx q[1];
rz(1.5093944) q[1];
rz(-pi) q[2];
rz(-0.13413818) q[3];
sx q[3];
rz(-0.87867498) q[3];
sx q[3];
rz(-1.1617171) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.0887289) q[2];
sx q[2];
rz(-2.1910618) q[2];
sx q[2];
rz(0.55727422) q[2];
rz(-0.82734621) q[3];
sx q[3];
rz(-1.5118303) q[3];
sx q[3];
rz(-1.5000337) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2130704) q[0];
sx q[0];
rz(-0.7970354) q[0];
sx q[0];
rz(2.573977) q[0];
rz(2.7396743) q[1];
sx q[1];
rz(-2.4979976) q[1];
sx q[1];
rz(1.3622805) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6797736) q[0];
sx q[0];
rz(-1.9711694) q[0];
sx q[0];
rz(-0.43433614) q[0];
rz(-pi) q[1];
rz(2.1902607) q[2];
sx q[2];
rz(-1.7835604) q[2];
sx q[2];
rz(1.484953) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.6475923) q[1];
sx q[1];
rz(-2.5060182) q[1];
sx q[1];
rz(-1.7032436) q[1];
rz(-pi) q[2];
rz(-0.79093334) q[3];
sx q[3];
rz(-0.91703992) q[3];
sx q[3];
rz(-2.4733651) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.42084971) q[2];
sx q[2];
rz(-1.2030615) q[2];
sx q[2];
rz(-2.8813349) q[2];
rz(2.4506954) q[3];
sx q[3];
rz(-0.8553718) q[3];
sx q[3];
rz(-1.5232085) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8061168) q[0];
sx q[0];
rz(-1.5609043) q[0];
sx q[0];
rz(1.5326395) q[0];
rz(0.43329049) q[1];
sx q[1];
rz(-1.3229803) q[1];
sx q[1];
rz(1.9569474) q[1];
rz(-0.717502) q[2];
sx q[2];
rz(-0.84894553) q[2];
sx q[2];
rz(0.15124627) q[2];
rz(0.79037249) q[3];
sx q[3];
rz(-1.4740667) q[3];
sx q[3];
rz(2.4355751) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
