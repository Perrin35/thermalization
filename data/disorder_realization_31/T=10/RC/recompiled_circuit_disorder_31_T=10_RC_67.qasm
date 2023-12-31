OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.6150317) q[0];
sx q[0];
rz(-0.57305133) q[0];
sx q[0];
rz(0.84258643) q[0];
rz(2.1057582) q[1];
sx q[1];
rz(8.3254568) q[1];
sx q[1];
rz(7.96666) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.18544491) q[0];
sx q[0];
rz(-1.3584134) q[0];
sx q[0];
rz(-0.74622112) q[0];
rz(-pi) q[1];
x q[1];
rz(3.1137142) q[2];
sx q[2];
rz(-2.5275702) q[2];
sx q[2];
rz(0.30464722) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.018651389) q[1];
sx q[1];
rz(-0.39995799) q[1];
sx q[1];
rz(0.33756983) q[1];
rz(-pi) q[2];
rz(1.4487212) q[3];
sx q[3];
rz(-0.97994084) q[3];
sx q[3];
rz(-1.6216244) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.52790102) q[2];
sx q[2];
rz(-1.0062904) q[2];
sx q[2];
rz(-2.9620985) q[2];
rz(1.2256631) q[3];
sx q[3];
rz(-1.3464728) q[3];
sx q[3];
rz(-0.82204449) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3935788) q[0];
sx q[0];
rz(-0.8809692) q[0];
sx q[0];
rz(-0.32546145) q[0];
rz(-1.356396) q[1];
sx q[1];
rz(-1.0486832) q[1];
sx q[1];
rz(-1.9869841) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1240631) q[0];
sx q[0];
rz(-1.5895827) q[0];
sx q[0];
rz(1.5536874) q[0];
x q[1];
rz(-2.7484659) q[2];
sx q[2];
rz(-0.98193411) q[2];
sx q[2];
rz(-1.7413505) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.11101152) q[1];
sx q[1];
rz(-1.6554553) q[1];
sx q[1];
rz(-2.3762523) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9566544) q[3];
sx q[3];
rz(-1.8293081) q[3];
sx q[3];
rz(-1.0227433) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.4521728) q[2];
sx q[2];
rz(-1.8916811) q[2];
sx q[2];
rz(0.88341218) q[2];
rz(-2.6702821) q[3];
sx q[3];
rz(-1.4383957) q[3];
sx q[3];
rz(-0.78770351) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8283591) q[0];
sx q[0];
rz(-1.6468843) q[0];
sx q[0];
rz(1.6261684) q[0];
rz(2.5405163) q[1];
sx q[1];
rz(-2.5939012) q[1];
sx q[1];
rz(-1.0916969) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1356782) q[0];
sx q[0];
rz(-2.1164829) q[0];
sx q[0];
rz(2.2360327) q[0];
rz(-2.2268779) q[2];
sx q[2];
rz(-1.9064184) q[2];
sx q[2];
rz(-0.4682954) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.25635438) q[1];
sx q[1];
rz(-1.9372205) q[1];
sx q[1];
rz(2.3430235) q[1];
x q[2];
rz(3.0136209) q[3];
sx q[3];
rz(-1.6347173) q[3];
sx q[3];
rz(-1.2325866) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.8213356) q[2];
sx q[2];
rz(-0.50575033) q[2];
sx q[2];
rz(0.88095218) q[2];
rz(-1.7679924) q[3];
sx q[3];
rz(-1.6146086) q[3];
sx q[3];
rz(-1.0176456) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.83051935) q[0];
sx q[0];
rz(-1.3922465) q[0];
sx q[0];
rz(-0.4367035) q[0];
rz(2.9084335) q[1];
sx q[1];
rz(-1.2522839) q[1];
sx q[1];
rz(-2.8312347) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.719602) q[0];
sx q[0];
rz(-1.4388226) q[0];
sx q[0];
rz(-2.7601526) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.15375806) q[2];
sx q[2];
rz(-2.4506844) q[2];
sx q[2];
rz(-1.549987) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.0536641) q[1];
sx q[1];
rz(-1.9287319) q[1];
sx q[1];
rz(-0.019836516) q[1];
x q[2];
rz(-1.849732) q[3];
sx q[3];
rz(-0.12860563) q[3];
sx q[3];
rz(-2.0898553) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-3.0115396) q[2];
sx q[2];
rz(-2.4235642) q[2];
sx q[2];
rz(1.0774353) q[2];
rz(3.0854026) q[3];
sx q[3];
rz(-0.63779938) q[3];
sx q[3];
rz(1.5475387) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.189165) q[0];
sx q[0];
rz(-1.0452894) q[0];
sx q[0];
rz(-0.24965832) q[0];
rz(1.5769618) q[1];
sx q[1];
rz(-2.3639634) q[1];
sx q[1];
rz(2.2713984) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91709671) q[0];
sx q[0];
rz(-1.3670237) q[0];
sx q[0];
rz(2.8208371) q[0];
rz(-pi) q[1];
rz(2.9179847) q[2];
sx q[2];
rz(-0.70906559) q[2];
sx q[2];
rz(2.2774334) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.40904564) q[1];
sx q[1];
rz(-0.88102075) q[1];
sx q[1];
rz(1.9802666) q[1];
rz(-0.16964511) q[3];
sx q[3];
rz(-1.0153474) q[3];
sx q[3];
rz(0.58562216) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.27328086) q[2];
sx q[2];
rz(-1.8153278) q[2];
sx q[2];
rz(-0.67374054) q[2];
rz(-0.30361787) q[3];
sx q[3];
rz(-1.9165336) q[3];
sx q[3];
rz(-1.822086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7917787) q[0];
sx q[0];
rz(-2.204201) q[0];
sx q[0];
rz(-0.2579903) q[0];
rz(-2.7164283) q[1];
sx q[1];
rz(-2.185967) q[1];
sx q[1];
rz(-1.4917096) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.30116044) q[0];
sx q[0];
rz(-1.4727117) q[0];
sx q[0];
rz(-0.43099404) q[0];
rz(-1.1509402) q[2];
sx q[2];
rz(-2.206415) q[2];
sx q[2];
rz(2.2395526) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.41960934) q[1];
sx q[1];
rz(-2.1188542) q[1];
sx q[1];
rz(-0.67668681) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.10156472) q[3];
sx q[3];
rz(-1.2079117) q[3];
sx q[3];
rz(0.18728072) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.1288746) q[2];
sx q[2];
rz(-2.1779163) q[2];
sx q[2];
rz(0.091726124) q[2];
rz(-2.2979459) q[3];
sx q[3];
rz(-0.97976145) q[3];
sx q[3];
rz(0.89404026) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3180852) q[0];
sx q[0];
rz(-1.2473236) q[0];
sx q[0];
rz(-2.7303625) q[0];
rz(-2.2757018) q[1];
sx q[1];
rz(-2.829268) q[1];
sx q[1];
rz(0.033989865) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86173979) q[0];
sx q[0];
rz(-1.732677) q[0];
sx q[0];
rz(0.78761657) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2074239) q[2];
sx q[2];
rz(-1.120943) q[2];
sx q[2];
rz(2.5069782) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.6778292) q[1];
sx q[1];
rz(-0.67968183) q[1];
sx q[1];
rz(1.6512647) q[1];
x q[2];
rz(2.9339318) q[3];
sx q[3];
rz(-0.62519473) q[3];
sx q[3];
rz(-3.0561662) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.5380481) q[2];
sx q[2];
rz(-0.59843439) q[2];
sx q[2];
rz(-0.87654385) q[2];
rz(-0.34902469) q[3];
sx q[3];
rz(-1.2004431) q[3];
sx q[3];
rz(2.9984737) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2440764) q[0];
sx q[0];
rz(-1.709047) q[0];
sx q[0];
rz(-2.7476655) q[0];
rz(-2.774033) q[1];
sx q[1];
rz(-1.7575248) q[1];
sx q[1];
rz(1.6961018) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1720393) q[0];
sx q[0];
rz(-1.5759828) q[0];
sx q[0];
rz(1.1374377) q[0];
rz(-pi) q[1];
x q[1];
rz(0.58480279) q[2];
sx q[2];
rz(-2.1091166) q[2];
sx q[2];
rz(2.1540097) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.98390025) q[1];
sx q[1];
rz(-2.1196113) q[1];
sx q[1];
rz(-1.7556612) q[1];
rz(-pi) q[2];
rz(-0.74650713) q[3];
sx q[3];
rz(-2.3254447) q[3];
sx q[3];
rz(-1.5619123) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.2945071) q[2];
sx q[2];
rz(-0.89035788) q[2];
sx q[2];
rz(-2.7344446) q[2];
rz(-1.6242705) q[3];
sx q[3];
rz(-1.1573236) q[3];
sx q[3];
rz(-2.8919162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3354934) q[0];
sx q[0];
rz(-2.6265916) q[0];
sx q[0];
rz(1.2517713) q[0];
rz(2.4720526) q[1];
sx q[1];
rz(-1.1839097) q[1];
sx q[1];
rz(-0.30977419) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99684925) q[0];
sx q[0];
rz(-2.092917) q[0];
sx q[0];
rz(-1.7168619) q[0];
rz(1.5179713) q[2];
sx q[2];
rz(-2.9644358) q[2];
sx q[2];
rz(-1.0747386) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.0032776912) q[1];
sx q[1];
rz(-1.365108) q[1];
sx q[1];
rz(0.059271952) q[1];
x q[2];
rz(-0.36848948) q[3];
sx q[3];
rz(-2.512305) q[3];
sx q[3];
rz(-0.027241782) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.70242515) q[2];
sx q[2];
rz(-0.71321407) q[2];
sx q[2];
rz(-1.2072198) q[2];
rz(1.0369982) q[3];
sx q[3];
rz(-1.8959277) q[3];
sx q[3];
rz(2.4859378) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0913775) q[0];
sx q[0];
rz(-1.3239048) q[0];
sx q[0];
rz(1.2058831) q[0];
rz(-0.58569113) q[1];
sx q[1];
rz(-1.0810477) q[1];
sx q[1];
rz(1.4996128) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8382032) q[0];
sx q[0];
rz(-1.8614385) q[0];
sx q[0];
rz(0.10586664) q[0];
x q[1];
rz(0.76403107) q[2];
sx q[2];
rz(-1.521763) q[2];
sx q[2];
rz(1.7977835) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.75913945) q[1];
sx q[1];
rz(-2.0331953) q[1];
sx q[1];
rz(-0.23666246) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5851192) q[3];
sx q[3];
rz(-0.17281547) q[3];
sx q[3];
rz(-2.1567791) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.2575834) q[2];
sx q[2];
rz(-1.7929701) q[2];
sx q[2];
rz(-0.94669) q[2];
rz(2.7729014) q[3];
sx q[3];
rz(-1.5654516) q[3];
sx q[3];
rz(2.6856016) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7832227) q[0];
sx q[0];
rz(-1.9932278) q[0];
sx q[0];
rz(2.7182462) q[0];
rz(-0.070925698) q[1];
sx q[1];
rz(-1.4535041) q[1];
sx q[1];
rz(2.8765875) q[1];
rz(-1.0735738) q[2];
sx q[2];
rz(-1.9855269) q[2];
sx q[2];
rz(-1.2863458) q[2];
rz(-0.72017097) q[3];
sx q[3];
rz(-1.1838893) q[3];
sx q[3];
rz(0.67223687) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
