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
rz(3.5377503) q[1];
sx q[1];
rz(9.1023394) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3381391) q[0];
sx q[0];
rz(-0.33284602) q[0];
sx q[0];
rz(-1.3071609) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2828322) q[2];
sx q[2];
rz(-2.3220064) q[2];
sx q[2];
rz(-0.31529266) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.8391708) q[1];
sx q[1];
rz(-1.2245373) q[1];
sx q[1];
rz(2.195921) q[1];
rz(-2.206771) q[3];
sx q[3];
rz(-1.1879731) q[3];
sx q[3];
rz(0.32360199) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.50513187) q[2];
sx q[2];
rz(-2.5487066) q[2];
sx q[2];
rz(0.55603975) q[2];
rz(-2.3089144) q[3];
sx q[3];
rz(-1.6502389) q[3];
sx q[3];
rz(0.94579831) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6933724) q[0];
sx q[0];
rz(-1.4602666) q[0];
sx q[0];
rz(0.15727501) q[0];
rz(2.8804624) q[1];
sx q[1];
rz(-1.3477247) q[1];
sx q[1];
rz(3.0325586) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.581668) q[0];
sx q[0];
rz(-1.7413057) q[0];
sx q[0];
rz(-1.3397564) q[0];
x q[1];
rz(-0.63983812) q[2];
sx q[2];
rz(-0.32391641) q[2];
sx q[2];
rz(1.6537635) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.038139548) q[1];
sx q[1];
rz(-2.1557169) q[1];
sx q[1];
rz(1.000688) q[1];
x q[2];
rz(0.29626366) q[3];
sx q[3];
rz(-2.6014572) q[3];
sx q[3];
rz(-1.6903282) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.2382425) q[2];
sx q[2];
rz(-1.976333) q[2];
sx q[2];
rz(-1.8781352) q[2];
rz(-0.3271099) q[3];
sx q[3];
rz(-1.5644904) q[3];
sx q[3];
rz(1.2143149) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.064421244) q[0];
sx q[0];
rz(-0.049296878) q[0];
sx q[0];
rz(-1.7984614) q[0];
rz(2.893977) q[1];
sx q[1];
rz(-0.74664465) q[1];
sx q[1];
rz(-2.6599191) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1789078) q[0];
sx q[0];
rz(-0.6324397) q[0];
sx q[0];
rz(-0.90895598) q[0];
rz(0.33294296) q[2];
sx q[2];
rz(-0.98368401) q[2];
sx q[2];
rz(-1.743403) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.18322769) q[1];
sx q[1];
rz(-2.1124358) q[1];
sx q[1];
rz(-0.72802131) q[1];
x q[2];
rz(-0.62044575) q[3];
sx q[3];
rz(-0.9471604) q[3];
sx q[3];
rz(-0.40943957) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.3383011) q[2];
sx q[2];
rz(-2.3241966) q[2];
sx q[2];
rz(2.6417007) q[2];
rz(-2.5806184) q[3];
sx q[3];
rz(-1.8818972) q[3];
sx q[3];
rz(1.6803754) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.19701476) q[0];
sx q[0];
rz(-0.16600969) q[0];
sx q[0];
rz(-0.552185) q[0];
rz(1.588297) q[1];
sx q[1];
rz(-0.89893666) q[1];
sx q[1];
rz(1.8968556) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8053999) q[0];
sx q[0];
rz(-1.5724143) q[0];
sx q[0];
rz(-0.33948116) q[0];
x q[1];
rz(-3.0110637) q[2];
sx q[2];
rz(-1.7703238) q[2];
sx q[2];
rz(-2.5752657) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(3.1340449) q[1];
sx q[1];
rz(-1.3589077) q[1];
sx q[1];
rz(1.0210387) q[1];
rz(-pi) q[2];
rz(1.3454868) q[3];
sx q[3];
rz(-1.011519) q[3];
sx q[3];
rz(2.6252281) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.2923979) q[2];
sx q[2];
rz(-1.8820102) q[2];
sx q[2];
rz(-1.1506895) q[2];
rz(-1.4771279) q[3];
sx q[3];
rz(-1.632558) q[3];
sx q[3];
rz(0.48294827) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0634336) q[0];
sx q[0];
rz(-0.76197356) q[0];
sx q[0];
rz(0.081469014) q[0];
rz(-0.062462656) q[1];
sx q[1];
rz(-1.1413347) q[1];
sx q[1];
rz(-1.6385471) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3969288) q[0];
sx q[0];
rz(-1.443112) q[0];
sx q[0];
rz(0.98608195) q[0];
rz(2.224515) q[2];
sx q[2];
rz(-0.13609016) q[2];
sx q[2];
rz(1.9249141) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.11322184) q[1];
sx q[1];
rz(-1.2591259) q[1];
sx q[1];
rz(0.10722864) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2522069) q[3];
sx q[3];
rz(-2.0584403) q[3];
sx q[3];
rz(-2.3209751) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.6333255) q[2];
sx q[2];
rz(-2.1990364) q[2];
sx q[2];
rz(-1.903669) q[2];
rz(2.0189019) q[3];
sx q[3];
rz(-0.676238) q[3];
sx q[3];
rz(-0.52156633) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5313107) q[0];
sx q[0];
rz(-2.1370482) q[0];
sx q[0];
rz(0.26671985) q[0];
rz(0.56089127) q[1];
sx q[1];
rz(-1.8436878) q[1];
sx q[1];
rz(-2.3430603) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7245367) q[0];
sx q[0];
rz(-2.8310611) q[0];
sx q[0];
rz(1.0945555) q[0];
rz(-1.919826) q[2];
sx q[2];
rz(-1.0997699) q[2];
sx q[2];
rz(2.7407321) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.27998566) q[1];
sx q[1];
rz(-1.0000739) q[1];
sx q[1];
rz(0.062203783) q[1];
rz(-pi) q[2];
rz(0.60243209) q[3];
sx q[3];
rz(-0.89655399) q[3];
sx q[3];
rz(-0.34365052) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.9810527) q[2];
sx q[2];
rz(-1.2489698) q[2];
sx q[2];
rz(2.7748761) q[2];
rz(-1.2612873) q[3];
sx q[3];
rz(-1.4533318) q[3];
sx q[3];
rz(0.096207531) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40903184) q[0];
sx q[0];
rz(-2.2213187) q[0];
sx q[0];
rz(-2.5352056) q[0];
rz(-2.9442893) q[1];
sx q[1];
rz(-1.1261255) q[1];
sx q[1];
rz(2.6775449) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.34940091) q[0];
sx q[0];
rz(-2.1878625) q[0];
sx q[0];
rz(2.6130136) q[0];
rz(-pi) q[1];
rz(-1.4887965) q[2];
sx q[2];
rz(-0.24157761) q[2];
sx q[2];
rz(1.2496003) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(3.0276427) q[1];
sx q[1];
rz(-2.0538035) q[1];
sx q[1];
rz(0.34160683) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4887605) q[3];
sx q[3];
rz(-1.335252) q[3];
sx q[3];
rz(-1.9627278) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.2074034) q[2];
sx q[2];
rz(-2.1384017) q[2];
sx q[2];
rz(0.25804538) q[2];
rz(1.9559654) q[3];
sx q[3];
rz(-1.6262755) q[3];
sx q[3];
rz(-0.0035088249) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.19514062) q[0];
sx q[0];
rz(-1.2807245) q[0];
sx q[0];
rz(-0.38129693) q[0];
rz(0.095245846) q[1];
sx q[1];
rz(-2.1691599) q[1];
sx q[1];
rz(1.7000748) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9863319) q[0];
sx q[0];
rz(-3.0763456) q[0];
sx q[0];
rz(2.9078527) q[0];
x q[1];
rz(-0.054585329) q[2];
sx q[2];
rz(-0.61331257) q[2];
sx q[2];
rz(1.3837981) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.86130202) q[1];
sx q[1];
rz(-2.960627) q[1];
sx q[1];
rz(2.3737532) q[1];
rz(-pi) q[2];
rz(2.8076595) q[3];
sx q[3];
rz(-1.5302368) q[3];
sx q[3];
rz(-3.0143152) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.1429446) q[2];
sx q[2];
rz(-2.7286077) q[2];
sx q[2];
rz(-2.9150035) q[2];
rz(2.6930124) q[3];
sx q[3];
rz(-1.6058763) q[3];
sx q[3];
rz(-0.8297689) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3806234) q[0];
sx q[0];
rz(-0.7614823) q[0];
sx q[0];
rz(1.7425849) q[0];
rz(2.8245068) q[1];
sx q[1];
rz(-1.6665019) q[1];
sx q[1];
rz(2.1549966) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.9771377) q[0];
sx q[0];
rz(-1.4129606) q[0];
sx q[0];
rz(2.6437003) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.87395845) q[2];
sx q[2];
rz(-1.9947589) q[2];
sx q[2];
rz(1.4615814) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.36531891) q[1];
sx q[1];
rz(-2.1734108) q[1];
sx q[1];
rz(1.3990632) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1499693) q[3];
sx q[3];
rz(-1.7289366) q[3];
sx q[3];
rz(-1.4827673) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.9265147) q[2];
sx q[2];
rz(-0.72700095) q[2];
sx q[2];
rz(0.40965664) q[2];
rz(-0.26327291) q[3];
sx q[3];
rz(-1.3116838) q[3];
sx q[3];
rz(0.51945654) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7423994) q[0];
sx q[0];
rz(-0.078646794) q[0];
sx q[0];
rz(1.4051399) q[0];
rz(-2.3204904) q[1];
sx q[1];
rz(-2.2228873) q[1];
sx q[1];
rz(-1.7260889) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.50354276) q[0];
sx q[0];
rz(-1.8068131) q[0];
sx q[0];
rz(2.80098) q[0];
rz(1.9344994) q[2];
sx q[2];
rz(-1.3377681) q[2];
sx q[2];
rz(1.9050913) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.0714598) q[1];
sx q[1];
rz(-1.3576641) q[1];
sx q[1];
rz(-2.9865584) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2763001) q[3];
sx q[3];
rz(-2.6101972) q[3];
sx q[3];
rz(2.6607799) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.4225509) q[2];
sx q[2];
rz(-2.8406403) q[2];
sx q[2];
rz(-3.0174875) q[2];
rz(0.96578807) q[3];
sx q[3];
rz(-1.4744759) q[3];
sx q[3];
rz(2.4805099) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.60349764) q[0];
sx q[0];
rz(-0.24833831) q[0];
sx q[0];
rz(-0.86059358) q[0];
rz(0.30766906) q[1];
sx q[1];
rz(-1.888702) q[1];
sx q[1];
rz(-1.9370334) q[1];
rz(-1.3901426) q[2];
sx q[2];
rz(-1.8271108) q[2];
sx q[2];
rz(1.9432632) q[2];
rz(-2.5849871) q[3];
sx q[3];
rz(-1.6986596) q[3];
sx q[3];
rz(-1.2996775) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
