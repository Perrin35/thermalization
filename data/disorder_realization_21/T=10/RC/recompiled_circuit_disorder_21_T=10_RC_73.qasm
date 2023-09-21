OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.055846) q[0];
sx q[0];
rz(-3.0598109) q[0];
sx q[0];
rz(2.6401289) q[0];
rz(1.4986829) q[1];
sx q[1];
rz(-2.745435) q[1];
sx q[1];
rz(-0.3224386) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8034536) q[0];
sx q[0];
rz(-2.8087466) q[0];
sx q[0];
rz(1.3071609) q[0];
rz(2.5311004) q[2];
sx q[2];
rz(-2.1571026) q[2];
sx q[2];
rz(-2.5551978) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.8391708) q[1];
sx q[1];
rz(-1.9170554) q[1];
sx q[1];
rz(2.195921) q[1];
rz(2.1665855) q[3];
sx q[3];
rz(-2.4132204) q[3];
sx q[3];
rz(-1.7155852) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.50513187) q[2];
sx q[2];
rz(-2.5487066) q[2];
sx q[2];
rz(0.55603975) q[2];
rz(0.83267823) q[3];
sx q[3];
rz(-1.6502389) q[3];
sx q[3];
rz(-2.1957943) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6933724) q[0];
sx q[0];
rz(-1.4602666) q[0];
sx q[0];
rz(-0.15727501) q[0];
rz(-0.26113025) q[1];
sx q[1];
rz(-1.3477247) q[1];
sx q[1];
rz(-0.10903407) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5058274) q[0];
sx q[0];
rz(-2.8553537) q[0];
sx q[0];
rz(0.92606996) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.63983812) q[2];
sx q[2];
rz(-0.32391641) q[2];
sx q[2];
rz(1.6537635) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.038139548) q[1];
sx q[1];
rz(-0.98587576) q[1];
sx q[1];
rz(2.1409047) q[1];
x q[2];
rz(-1.7440967) q[3];
sx q[3];
rz(-1.0565851) q[3];
sx q[3];
rz(1.3483931) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.9033501) q[2];
sx q[2];
rz(-1.1652596) q[2];
sx q[2];
rz(1.2634574) q[2];
rz(-2.8144828) q[3];
sx q[3];
rz(-1.5644904) q[3];
sx q[3];
rz(-1.2143149) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0771714) q[0];
sx q[0];
rz(-3.0922958) q[0];
sx q[0];
rz(-1.7984614) q[0];
rz(2.893977) q[1];
sx q[1];
rz(-0.74664465) q[1];
sx q[1];
rz(0.48167357) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0945064) q[0];
sx q[0];
rz(-1.1990093) q[0];
sx q[0];
rz(2.0949754) q[0];
x q[1];
rz(-1.1142251) q[2];
sx q[2];
rz(-0.66514981) q[2];
sx q[2];
rz(-0.84012023) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.18322769) q[1];
sx q[1];
rz(-2.1124358) q[1];
sx q[1];
rz(2.4135713) q[1];
x q[2];
rz(0.84677245) q[3];
sx q[3];
rz(-2.0623042) q[3];
sx q[3];
rz(1.5566952) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.3383011) q[2];
sx q[2];
rz(-2.3241966) q[2];
sx q[2];
rz(-2.6417007) q[2];
rz(0.56097427) q[3];
sx q[3];
rz(-1.2596954) q[3];
sx q[3];
rz(-1.6803754) q[3];
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
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.19701476) q[0];
sx q[0];
rz(-2.975583) q[0];
sx q[0];
rz(0.552185) q[0];
rz(-1.588297) q[1];
sx q[1];
rz(-2.242656) q[1];
sx q[1];
rz(1.8968556) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9024076) q[0];
sx q[0];
rz(-0.33948487) q[0];
sx q[0];
rz(-0.0048588077) q[0];
x q[1];
rz(-0.99889465) q[2];
sx q[2];
rz(-0.23795393) q[2];
sx q[2];
rz(-1.9902802) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.2326395) q[1];
sx q[1];
rz(-2.5563572) q[1];
sx q[1];
rz(-1.1802243) q[1];
rz(-0.57085412) q[3];
sx q[3];
rz(-1.7613162) q[3];
sx q[3];
rz(1.1754456) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.84919471) q[2];
sx q[2];
rz(-1.2595824) q[2];
sx q[2];
rz(1.1506895) q[2];
rz(-1.4771279) q[3];
sx q[3];
rz(-1.632558) q[3];
sx q[3];
rz(0.48294827) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(3.0634336) q[0];
sx q[0];
rz(-0.76197356) q[0];
sx q[0];
rz(3.0601236) q[0];
rz(-3.07913) q[1];
sx q[1];
rz(-2.000258) q[1];
sx q[1];
rz(-1.6385471) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7446639) q[0];
sx q[0];
rz(-1.443112) q[0];
sx q[0];
rz(-0.98608195) q[0];
x q[1];
rz(1.6790752) q[2];
sx q[2];
rz(-1.6533972) q[2];
sx q[2];
rz(0.29512197) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(3.0283708) q[1];
sx q[1];
rz(-1.2591259) q[1];
sx q[1];
rz(0.10722864) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8893858) q[3];
sx q[3];
rz(-1.0831523) q[3];
sx q[3];
rz(-0.82061758) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.5082671) q[2];
sx q[2];
rz(-2.1990364) q[2];
sx q[2];
rz(1.2379237) q[2];
rz(-2.0189019) q[3];
sx q[3];
rz(-2.4653547) q[3];
sx q[3];
rz(2.6200263) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6102819) q[0];
sx q[0];
rz(-2.1370482) q[0];
sx q[0];
rz(-2.8748728) q[0];
rz(-2.5807014) q[1];
sx q[1];
rz(-1.2979049) q[1];
sx q[1];
rz(-0.7985324) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0621588) q[0];
sx q[0];
rz(-1.2957797) q[0];
sx q[0];
rz(-2.9955203) q[0];
rz(-2.550225) q[2];
sx q[2];
rz(-2.5632576) q[2];
sx q[2];
rz(-0.27507281) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.861607) q[1];
sx q[1];
rz(-2.1415188) q[1];
sx q[1];
rz(3.0793889) q[1];
rz(-0.9540326) q[3];
sx q[3];
rz(-0.87152374) q[3];
sx q[3];
rz(-1.9642533) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.9810527) q[2];
sx q[2];
rz(-1.8926228) q[2];
sx q[2];
rz(0.36671656) q[2];
rz(1.2612873) q[3];
sx q[3];
rz(-1.6882608) q[3];
sx q[3];
rz(0.096207531) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40903184) q[0];
sx q[0];
rz(-2.2213187) q[0];
sx q[0];
rz(0.60638705) q[0];
rz(-2.9442893) q[1];
sx q[1];
rz(-1.1261255) q[1];
sx q[1];
rz(2.6775449) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7921917) q[0];
sx q[0];
rz(-2.1878625) q[0];
sx q[0];
rz(-0.52857907) q[0];
rz(-pi) q[1];
x q[1];
rz(3.1214141) q[2];
sx q[2];
rz(-1.8115461) q[2];
sx q[2];
rz(1.9764331) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(3.0276427) q[1];
sx q[1];
rz(-1.0877891) q[1];
sx q[1];
rz(-0.34160683) q[1];
rz(-pi) q[2];
x q[2];
rz(0.23630996) q[3];
sx q[3];
rz(-1.650562) q[3];
sx q[3];
rz(-0.41111708) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.93418926) q[2];
sx q[2];
rz(-1.0031909) q[2];
sx q[2];
rz(-2.8835473) q[2];
rz(1.9559654) q[3];
sx q[3];
rz(-1.5153171) q[3];
sx q[3];
rz(0.0035088249) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.19514062) q[0];
sx q[0];
rz(-1.8608681) q[0];
sx q[0];
rz(-0.38129693) q[0];
rz(3.0463468) q[1];
sx q[1];
rz(-0.97243273) q[1];
sx q[1];
rz(1.7000748) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7521116) q[0];
sx q[0];
rz(-1.6342667) q[0];
sx q[0];
rz(-1.5556637) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.61261119) q[2];
sx q[2];
rz(-1.6022041) q[2];
sx q[2];
rz(-2.999246) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(3.0563593) q[1];
sx q[1];
rz(-1.7006405) q[1];
sx q[1];
rz(1.4443881) q[1];
rz(-pi) q[2];
rz(-3.0184047) q[3];
sx q[3];
rz(-2.8052969) q[3];
sx q[3];
rz(-1.5598701) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.1429446) q[2];
sx q[2];
rz(-2.7286077) q[2];
sx q[2];
rz(2.9150035) q[2];
rz(2.6930124) q[3];
sx q[3];
rz(-1.5357163) q[3];
sx q[3];
rz(0.8297689) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3806234) q[0];
sx q[0];
rz(-0.7614823) q[0];
sx q[0];
rz(-1.3990078) q[0];
rz(-2.8245068) q[1];
sx q[1];
rz(-1.4750907) q[1];
sx q[1];
rz(-0.98659602) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.164455) q[0];
sx q[0];
rz(-1.4129606) q[0];
sx q[0];
rz(2.6437003) q[0];
rz(-pi) q[1];
rz(-0.53194745) q[2];
sx q[2];
rz(-2.1956458) q[2];
sx q[2];
rz(0.22235409) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.0341067) q[1];
sx q[1];
rz(-1.7120546) q[1];
sx q[1];
rz(-0.60955255) q[1];
rz(-pi) q[2];
x q[2];
rz(0.17296965) q[3];
sx q[3];
rz(-1.9860455) q[3];
sx q[3];
rz(0.017661974) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.2150779) q[2];
sx q[2];
rz(-2.4145917) q[2];
sx q[2];
rz(-0.40965664) q[2];
rz(2.8783197) q[3];
sx q[3];
rz(-1.3116838) q[3];
sx q[3];
rz(0.51945654) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39919329) q[0];
sx q[0];
rz(-3.0629459) q[0];
sx q[0];
rz(1.4051399) q[0];
rz(0.82110226) q[1];
sx q[1];
rz(-2.2228873) q[1];
sx q[1];
rz(-1.7260889) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.50354276) q[0];
sx q[0];
rz(-1.3347795) q[0];
sx q[0];
rz(0.34061265) q[0];
x q[1];
rz(1.2070933) q[2];
sx q[2];
rz(-1.3377681) q[2];
sx q[2];
rz(1.2365014) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.532383) q[1];
sx q[1];
rz(-1.7222952) q[1];
sx q[1];
rz(1.7864368) q[1];
rz(-1.2763001) q[3];
sx q[3];
rz(-2.6101972) q[3];
sx q[3];
rz(-2.6607799) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.4225509) q[2];
sx q[2];
rz(-0.30095235) q[2];
sx q[2];
rz(0.12410513) q[2];
rz(2.1758046) q[3];
sx q[3];
rz(-1.6671168) q[3];
sx q[3];
rz(2.4805099) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
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
rz(2.538095) q[0];
sx q[0];
rz(-0.24833831) q[0];
sx q[0];
rz(-0.86059358) q[0];
rz(-2.8339236) q[1];
sx q[1];
rz(-1.888702) q[1];
sx q[1];
rz(-1.9370334) q[1];
rz(-2.5406191) q[2];
sx q[2];
rz(-0.31243159) q[2];
sx q[2];
rz(2.5675788) q[2];
rz(0.23871213) q[3];
sx q[3];
rz(-2.5720027) q[3];
sx q[3];
rz(-3.0726074) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];