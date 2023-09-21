OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.9392202) q[0];
sx q[0];
rz(-0.4063172) q[0];
sx q[0];
rz(0.82011861) q[0];
rz(-0.36110538) q[1];
sx q[1];
rz(0.63280025) q[1];
sx q[1];
rz(11.735698) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.75738534) q[0];
sx q[0];
rz(-1.2641347) q[0];
sx q[0];
rz(-0.39461179) q[0];
rz(-pi) q[1];
rz(1.0933502) q[2];
sx q[2];
rz(-2.2508143) q[2];
sx q[2];
rz(0.29104656) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.0349866) q[1];
sx q[1];
rz(-2.258746) q[1];
sx q[1];
rz(2.7290542) q[1];
x q[2];
rz(-1.9740231) q[3];
sx q[3];
rz(-1.9864169) q[3];
sx q[3];
rz(-0.88413873) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.43710199) q[2];
sx q[2];
rz(-0.4318684) q[2];
sx q[2];
rz(-0.1201771) q[2];
rz(-1.1581356) q[3];
sx q[3];
rz(-1.7430256) q[3];
sx q[3];
rz(2.2226298) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.08081089) q[0];
sx q[0];
rz(-1.3588384) q[0];
sx q[0];
rz(2.2297915) q[0];
rz(2.3520825) q[1];
sx q[1];
rz(-2.1531838) q[1];
sx q[1];
rz(-0.3266913) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.806843) q[0];
sx q[0];
rz(-2.202889) q[0];
sx q[0];
rz(0.64081162) q[0];
rz(-pi) q[1];
rz(-0.15906449) q[2];
sx q[2];
rz(-1.1317562) q[2];
sx q[2];
rz(2.1829) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.5888728) q[1];
sx q[1];
rz(-2.1761804) q[1];
sx q[1];
rz(-2.6436716) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4665524) q[3];
sx q[3];
rz(-0.6904656) q[3];
sx q[3];
rz(-3.1363827) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.6313173) q[2];
sx q[2];
rz(-0.77284208) q[2];
sx q[2];
rz(-1.2871683) q[2];
rz(-3.0316947) q[3];
sx q[3];
rz(-1.4108312) q[3];
sx q[3];
rz(-1.7597642) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.54863769) q[0];
sx q[0];
rz(-0.73919636) q[0];
sx q[0];
rz(-0.32989311) q[0];
rz(-0.27711162) q[1];
sx q[1];
rz(-1.8246633) q[1];
sx q[1];
rz(-2.0842016) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49845093) q[0];
sx q[0];
rz(-2.8006449) q[0];
sx q[0];
rz(1.6324415) q[0];
rz(-pi) q[1];
rz(0.51809394) q[2];
sx q[2];
rz(-0.38803852) q[2];
sx q[2];
rz(2.5990017) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.33144618) q[1];
sx q[1];
rz(-0.76008893) q[1];
sx q[1];
rz(2.0600832) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7626761) q[3];
sx q[3];
rz(-1.9865611) q[3];
sx q[3];
rz(-0.92431812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.3588336) q[2];
sx q[2];
rz(-1.1153355) q[2];
sx q[2];
rz(1.3624297) q[2];
rz(-0.6247012) q[3];
sx q[3];
rz(-2.0420859) q[3];
sx q[3];
rz(0.81956285) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3574922) q[0];
sx q[0];
rz(-0.52629137) q[0];
sx q[0];
rz(-0.5350565) q[0];
rz(-1.1401945) q[1];
sx q[1];
rz(-1.813872) q[1];
sx q[1];
rz(-2.9761956) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9388158) q[0];
sx q[0];
rz(-1.446615) q[0];
sx q[0];
rz(1.3846272) q[0];
rz(-0.8382767) q[2];
sx q[2];
rz(-0.24176134) q[2];
sx q[2];
rz(-1.4571112) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.4277028) q[1];
sx q[1];
rz(-2.1117003) q[1];
sx q[1];
rz(1.9546024) q[1];
rz(-1.0159675) q[3];
sx q[3];
rz(-1.3362243) q[3];
sx q[3];
rz(2.6382584) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.7455204) q[2];
sx q[2];
rz(-1.3362276) q[2];
sx q[2];
rz(-0.20544927) q[2];
rz(1.127634) q[3];
sx q[3];
rz(-1.9879568) q[3];
sx q[3];
rz(1.2566465) q[3];
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
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4797392) q[0];
sx q[0];
rz(-2.2275708) q[0];
sx q[0];
rz(2.7868295) q[0];
rz(1.9873437) q[1];
sx q[1];
rz(-0.92461363) q[1];
sx q[1];
rz(2.9096471) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35271586) q[0];
sx q[0];
rz(-2.0577601) q[0];
sx q[0];
rz(2.6608174) q[0];
rz(-1.4300214) q[2];
sx q[2];
rz(-2.2614334) q[2];
sx q[2];
rz(1.8292793) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.9310589) q[1];
sx q[1];
rz(-1.3506883) q[1];
sx q[1];
rz(-1.348043) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7730745) q[3];
sx q[3];
rz(-1.5889865) q[3];
sx q[3];
rz(-0.25572488) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.0014687) q[2];
sx q[2];
rz(-1.3342369) q[2];
sx q[2];
rz(-1.2403963) q[2];
rz(-0.59605789) q[3];
sx q[3];
rz(-1.3052992) q[3];
sx q[3];
rz(1.3034472) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0971138) q[0];
sx q[0];
rz(-0.070274027) q[0];
sx q[0];
rz(0.20275673) q[0];
rz(0.98908201) q[1];
sx q[1];
rz(-1.443807) q[1];
sx q[1];
rz(2.1441377) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7339242) q[0];
sx q[0];
rz(-1.7051538) q[0];
sx q[0];
rz(2.2785447) q[0];
rz(-pi) q[1];
rz(2.7315797) q[2];
sx q[2];
rz(-2.3668681) q[2];
sx q[2];
rz(1.1020401) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.5207386) q[1];
sx q[1];
rz(-1.3315017) q[1];
sx q[1];
rz(1.1362856) q[1];
x q[2];
rz(-1.0783844) q[3];
sx q[3];
rz(-0.78740722) q[3];
sx q[3];
rz(1.2315962) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.9138907) q[2];
sx q[2];
rz(-1.1662377) q[2];
sx q[2];
rz(2.690199) q[2];
rz(2.732892) q[3];
sx q[3];
rz(-1.6025851) q[3];
sx q[3];
rz(1.2020948) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.30615) q[0];
sx q[0];
rz(-0.4040443) q[0];
sx q[0];
rz(2.5174482) q[0];
rz(1.5165326) q[1];
sx q[1];
rz(-2.8846278) q[1];
sx q[1];
rz(-2.5278032) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3261624) q[0];
sx q[0];
rz(-0.76061941) q[0];
sx q[0];
rz(-1.8421696) q[0];
rz(-0.21653793) q[2];
sx q[2];
rz(-0.29007402) q[2];
sx q[2];
rz(-1.9921583) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.9871414) q[1];
sx q[1];
rz(-0.11949355) q[1];
sx q[1];
rz(2.6277072) q[1];
rz(0.11717637) q[3];
sx q[3];
rz(-2.0332608) q[3];
sx q[3];
rz(1.5069435) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.43341407) q[2];
sx q[2];
rz(-2.2257979) q[2];
sx q[2];
rz(1.1748574) q[2];
rz(2.5332149) q[3];
sx q[3];
rz(-1.7374246) q[3];
sx q[3];
rz(1.7181989) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90010086) q[0];
sx q[0];
rz(-1.8429723) q[0];
sx q[0];
rz(-0.4883782) q[0];
rz(1.6237367) q[1];
sx q[1];
rz(-1.7428215) q[1];
sx q[1];
rz(0.98446313) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1293837) q[0];
sx q[0];
rz(-2.1166271) q[0];
sx q[0];
rz(-1.6145541) q[0];
x q[1];
rz(1.9756873) q[2];
sx q[2];
rz(-2.3003909) q[2];
sx q[2];
rz(2.0617495) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.3417495) q[1];
sx q[1];
rz(-1.5768331) q[1];
sx q[1];
rz(0.36532613) q[1];
rz(0.47655388) q[3];
sx q[3];
rz(-1.2911951) q[3];
sx q[3];
rz(0.24147803) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.4349334) q[2];
sx q[2];
rz(-1.443053) q[2];
sx q[2];
rz(-2.6064176) q[2];
rz(1.0501856) q[3];
sx q[3];
rz(-1.340056) q[3];
sx q[3];
rz(2.1323269) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.500279) q[0];
sx q[0];
rz(-2.2135493) q[0];
sx q[0];
rz(-2.4556659) q[0];
rz(-2.7507239) q[1];
sx q[1];
rz(-0.94376826) q[1];
sx q[1];
rz(0.92591441) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6305011) q[0];
sx q[0];
rz(-1.4303659) q[0];
sx q[0];
rz(-1.3569843) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.8530059) q[2];
sx q[2];
rz(-0.91663137) q[2];
sx q[2];
rz(-2.5329563) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.029286413) q[1];
sx q[1];
rz(-2.6767113) q[1];
sx q[1];
rz(-1.2390562) q[1];
x q[2];
rz(-1.6894475) q[3];
sx q[3];
rz(-1.7735529) q[3];
sx q[3];
rz(-0.53393902) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.19568504) q[2];
sx q[2];
rz(-2.6699799) q[2];
sx q[2];
rz(2.6055028) q[2];
rz(2.7219971) q[3];
sx q[3];
rz(-2.5511238) q[3];
sx q[3];
rz(-1.6528116) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0560028) q[0];
sx q[0];
rz(-2.7828126) q[0];
sx q[0];
rz(0.39500239) q[0];
rz(-1.6292054) q[1];
sx q[1];
rz(-1.2724266) q[1];
sx q[1];
rz(1.013247) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.17908827) q[0];
sx q[0];
rz(-2.8694186) q[0];
sx q[0];
rz(-2.1129235) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7878694) q[2];
sx q[2];
rz(-1.5295267) q[2];
sx q[2];
rz(-1.313414) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.0777178) q[1];
sx q[1];
rz(-2.7091654) q[1];
sx q[1];
rz(0.23250154) q[1];
x q[2];
rz(0.89684422) q[3];
sx q[3];
rz(-0.92771155) q[3];
sx q[3];
rz(-2.602592) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.6932678) q[2];
sx q[2];
rz(-1.6222745) q[2];
sx q[2];
rz(-2.3821793) q[2];
rz(1.7761207) q[3];
sx q[3];
rz(-1.1080192) q[3];
sx q[3];
rz(2.571648) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7286745) q[0];
sx q[0];
rz(-1.8771794) q[0];
sx q[0];
rz(-1.9369453) q[0];
rz(-0.57327523) q[1];
sx q[1];
rz(-1.234006) q[1];
sx q[1];
rz(-1.3201859) q[1];
rz(-1.7309932) q[2];
sx q[2];
rz(-2.6916531) q[2];
sx q[2];
rz(2.0315363) q[2];
rz(2.2686601) q[3];
sx q[3];
rz(-1.6193661) q[3];
sx q[3];
rz(-1.247874) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
