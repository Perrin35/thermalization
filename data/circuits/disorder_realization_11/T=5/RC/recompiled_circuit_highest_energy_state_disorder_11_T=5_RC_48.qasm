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
rz(1.4229245) q[0];
sx q[0];
rz(-2.0473502) q[0];
sx q[0];
rz(-0.25804582) q[0];
rz(2.1482422) q[1];
sx q[1];
rz(-1.300783) q[1];
sx q[1];
rz(0.25564495) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6171744) q[0];
sx q[0];
rz(-2.1416695) q[0];
sx q[0];
rz(-0.85321315) q[0];
rz(0.3886224) q[2];
sx q[2];
rz(-1.5124694) q[2];
sx q[2];
rz(-1.7494534) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.5140443) q[1];
sx q[1];
rz(-1.5020292) q[1];
sx q[1];
rz(-1.7538944) q[1];
rz(-0.49778865) q[3];
sx q[3];
rz(-1.3827795) q[3];
sx q[3];
rz(1.1527485) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.51097441) q[2];
sx q[2];
rz(-1.3773842) q[2];
sx q[2];
rz(1.1154491) q[2];
rz(-0.014178064) q[3];
sx q[3];
rz(-1.7970128) q[3];
sx q[3];
rz(0.43629638) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2071335) q[0];
sx q[0];
rz(-2.5196228) q[0];
sx q[0];
rz(-1.2472664) q[0];
rz(2.6994052) q[1];
sx q[1];
rz(-1.7555883) q[1];
sx q[1];
rz(2.3242548) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2984262) q[0];
sx q[0];
rz(-1.048537) q[0];
sx q[0];
rz(2.9199615) q[0];
rz(-pi) q[1];
rz(-2.1139718) q[2];
sx q[2];
rz(-0.50293844) q[2];
sx q[2];
rz(-0.78928141) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.0603588) q[1];
sx q[1];
rz(-2.4801755) q[1];
sx q[1];
rz(-0.27137406) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4659791) q[3];
sx q[3];
rz(-0.88425501) q[3];
sx q[3];
rz(3.0385449) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.3680129) q[2];
sx q[2];
rz(-2.7648338) q[2];
sx q[2];
rz(0.22029857) q[2];
rz(-1.1822654) q[3];
sx q[3];
rz(-1.9672829) q[3];
sx q[3];
rz(-0.063145414) q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0910864) q[0];
sx q[0];
rz(-1.3171221) q[0];
sx q[0];
rz(2.0029946) q[0];
rz(0.41172045) q[1];
sx q[1];
rz(-0.76985923) q[1];
sx q[1];
rz(-1.0134816) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0385376) q[0];
sx q[0];
rz(-2.8881209) q[0];
sx q[0];
rz(-0.059132476) q[0];
rz(1.2960245) q[2];
sx q[2];
rz(-2.6856075) q[2];
sx q[2];
rz(2.6528461) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.1057844) q[1];
sx q[1];
rz(-2.2757451) q[1];
sx q[1];
rz(2.1501599) q[1];
rz(0.8441505) q[3];
sx q[3];
rz(-1.6682079) q[3];
sx q[3];
rz(1.0507492) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.98075214) q[2];
sx q[2];
rz(-1.1557121) q[2];
sx q[2];
rz(1.2551003) q[2];
rz(0.27215019) q[3];
sx q[3];
rz(-2.3641219) q[3];
sx q[3];
rz(-0.14061558) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.960152) q[0];
sx q[0];
rz(-0.93638268) q[0];
sx q[0];
rz(0.61035672) q[0];
rz(-2.4329674) q[1];
sx q[1];
rz(-1.0162063) q[1];
sx q[1];
rz(-1.5708539) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8886531) q[0];
sx q[0];
rz(-1.4331237) q[0];
sx q[0];
rz(0.15451365) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6897292) q[2];
sx q[2];
rz(-1.8557252) q[2];
sx q[2];
rz(-2.1708084) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.849087) q[1];
sx q[1];
rz(-1.5226814) q[1];
sx q[1];
rz(2.8699257) q[1];
rz(-2.3409178) q[3];
sx q[3];
rz(-1.281637) q[3];
sx q[3];
rz(-1.9137933) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.9307956) q[2];
sx q[2];
rz(-1.0127298) q[2];
sx q[2];
rz(-2.5861758) q[2];
rz(-0.72426116) q[3];
sx q[3];
rz(-1.9925502) q[3];
sx q[3];
rz(-0.65653062) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.637218) q[0];
sx q[0];
rz(-1.9509622) q[0];
sx q[0];
rz(2.7727238) q[0];
rz(-0.24636191) q[1];
sx q[1];
rz(-1.8135704) q[1];
sx q[1];
rz(-1.7074283) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96077641) q[0];
sx q[0];
rz(-0.81396507) q[0];
sx q[0];
rz(-1.5778944) q[0];
rz(2.6961961) q[2];
sx q[2];
rz(-1.4120308) q[2];
sx q[2];
rz(-2.7743055) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.37232698) q[1];
sx q[1];
rz(-2.7397836) q[1];
sx q[1];
rz(0.010784464) q[1];
rz(-3.0415972) q[3];
sx q[3];
rz(-1.7971562) q[3];
sx q[3];
rz(0.55734998) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.4904334) q[2];
sx q[2];
rz(-1.8305402) q[2];
sx q[2];
rz(2.0595713) q[2];
rz(-0.79484445) q[3];
sx q[3];
rz(-1.4719897) q[3];
sx q[3];
rz(1.2580416) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
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
rz(0.27570462) q[0];
sx q[0];
rz(-3.0401433) q[0];
sx q[0];
rz(-0.24359447) q[0];
rz(0.98681915) q[1];
sx q[1];
rz(-0.74816626) q[1];
sx q[1];
rz(-2.2056244) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2631131) q[0];
sx q[0];
rz(-1.2498444) q[0];
sx q[0];
rz(1.1522989) q[0];
rz(-pi) q[1];
x q[1];
rz(0.82359969) q[2];
sx q[2];
rz(-1.3098426) q[2];
sx q[2];
rz(-1.5023155) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.3159891) q[1];
sx q[1];
rz(-2.2554419) q[1];
sx q[1];
rz(-2.0988093) q[1];
rz(-pi) q[2];
rz(-0.97096393) q[3];
sx q[3];
rz(-1.55313) q[3];
sx q[3];
rz(-2.3223557) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.662107) q[2];
sx q[2];
rz(-2.0027436) q[2];
sx q[2];
rz(-0.34995079) q[2];
rz(2.5838666) q[3];
sx q[3];
rz(-1.2099268) q[3];
sx q[3];
rz(0.85132712) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6996985) q[0];
sx q[0];
rz(-2.5040369) q[0];
sx q[0];
rz(0.30174524) q[0];
rz(-2.8217577) q[1];
sx q[1];
rz(-1.6022857) q[1];
sx q[1];
rz(2.005827) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5556896) q[0];
sx q[0];
rz(-1.5169889) q[0];
sx q[0];
rz(-2.4967628) q[0];
x q[1];
rz(-2.6276845) q[2];
sx q[2];
rz(-1.0036047) q[2];
sx q[2];
rz(-2.6197768) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.9965976) q[1];
sx q[1];
rz(-0.97255822) q[1];
sx q[1];
rz(-0.096417565) q[1];
rz(-0.79085042) q[3];
sx q[3];
rz(-1.6501657) q[3];
sx q[3];
rz(-0.96109238) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.4096421) q[2];
sx q[2];
rz(-2.4739517) q[2];
sx q[2];
rz(-0.14275924) q[2];
rz(1.7536633) q[3];
sx q[3];
rz(-1.9526491) q[3];
sx q[3];
rz(2.4587542) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1801572) q[0];
sx q[0];
rz(-0.016594369) q[0];
sx q[0];
rz(-0.55602443) q[0];
rz(2.9122638) q[1];
sx q[1];
rz(-1.2622958) q[1];
sx q[1];
rz(1.75846) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1616104) q[0];
sx q[0];
rz(-1.5357657) q[0];
sx q[0];
rz(-2.4039414) q[0];
rz(0.11923125) q[2];
sx q[2];
rz(-2.6802353) q[2];
sx q[2];
rz(2.2329494) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.8236716) q[1];
sx q[1];
rz(-1.2814199) q[1];
sx q[1];
rz(1.1422864) q[1];
rz(-pi) q[2];
rz(1.9541627) q[3];
sx q[3];
rz(-1.995242) q[3];
sx q[3];
rz(-0.92680537) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.51656276) q[2];
sx q[2];
rz(-0.27250686) q[2];
sx q[2];
rz(-1.9005091) q[2];
rz(-3.0789913) q[3];
sx q[3];
rz(-1.4227941) q[3];
sx q[3];
rz(0.82714287) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
sx q[3];
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
rz(-1.0144219) q[0];
sx q[0];
rz(-1.7272471) q[0];
sx q[0];
rz(2.993809) q[0];
rz(-2.6328909) q[1];
sx q[1];
rz(-1.3721507) q[1];
sx q[1];
rz(2.7630189) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7258647) q[0];
sx q[0];
rz(-1.8047389) q[0];
sx q[0];
rz(-0.42629231) q[0];
x q[1];
rz(0.97855391) q[2];
sx q[2];
rz(-1.2662953) q[2];
sx q[2];
rz(2.1155807) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.3562153) q[1];
sx q[1];
rz(-1.2357724) q[1];
sx q[1];
rz(2.8122693) q[1];
x q[2];
rz(-1.6949953) q[3];
sx q[3];
rz(-2.327965) q[3];
sx q[3];
rz(-0.99909335) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.1754237) q[2];
sx q[2];
rz(-2.1899026) q[2];
sx q[2];
rz(-1.8625205) q[2];
rz(1.7133948) q[3];
sx q[3];
rz(-2.1396075) q[3];
sx q[3];
rz(2.9960347) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
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
rz(-1.3928423) q[0];
sx q[0];
rz(-0.57806438) q[0];
sx q[0];
rz(-3.0774935) q[0];
rz(-2.6129258) q[1];
sx q[1];
rz(-1.8730947) q[1];
sx q[1];
rz(2.166523) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9720358) q[0];
sx q[0];
rz(-1.0825048) q[0];
sx q[0];
rz(-1.8711062) q[0];
rz(-pi) q[1];
rz(-2.5717402) q[2];
sx q[2];
rz(-1.0625216) q[2];
sx q[2];
rz(0.13392553) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.9916315) q[1];
sx q[1];
rz(-1.35222) q[1];
sx q[1];
rz(0.44709713) q[1];
rz(-pi) q[2];
rz(-2.3894839) q[3];
sx q[3];
rz(-0.32882133) q[3];
sx q[3];
rz(-0.11550918) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.9958682) q[2];
sx q[2];
rz(-2.4805562) q[2];
sx q[2];
rz(-0.60689849) q[2];
rz(2.8912344) q[3];
sx q[3];
rz(-1.4162049) q[3];
sx q[3];
rz(2.6355766) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.09457) q[0];
sx q[0];
rz(-1.2170412) q[0];
sx q[0];
rz(1.8893597) q[0];
rz(0.49867123) q[1];
sx q[1];
rz(-1.55232) q[1];
sx q[1];
rz(1.3943863) q[1];
rz(-2.2398938) q[2];
sx q[2];
rz(-0.82868262) q[2];
sx q[2];
rz(2.5862624) q[2];
rz(-1.3648894) q[3];
sx q[3];
rz(-1.6905224) q[3];
sx q[3];
rz(-0.57991309) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
