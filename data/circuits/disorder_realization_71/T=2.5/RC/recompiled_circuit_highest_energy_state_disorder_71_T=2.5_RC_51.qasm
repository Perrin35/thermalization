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
rz(-2.5833997) q[0];
sx q[0];
rz(-0.71891958) q[0];
sx q[0];
rz(-0.67607003) q[0];
rz(-2.8683635) q[1];
sx q[1];
rz(-0.9571119) q[1];
sx q[1];
rz(2.2303384) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8000923) q[0];
sx q[0];
rz(-1.5908634) q[0];
sx q[0];
rz(-1.5485199) q[0];
x q[1];
rz(0.62556556) q[2];
sx q[2];
rz(-1.5766597) q[2];
sx q[2];
rz(-1.0926343) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.75685749) q[1];
sx q[1];
rz(-2.5674106) q[1];
sx q[1];
rz(-0.059348182) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3852547) q[3];
sx q[3];
rz(-1.66258) q[3];
sx q[3];
rz(0.97912195) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.91997826) q[2];
sx q[2];
rz(-2.6875434) q[2];
sx q[2];
rz(2.3294219) q[2];
rz(-3.0657366) q[3];
sx q[3];
rz(-0.64801884) q[3];
sx q[3];
rz(2.5465452) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49038637) q[0];
sx q[0];
rz(-0.47845978) q[0];
sx q[0];
rz(-1.867021) q[0];
rz(-1.5050585) q[1];
sx q[1];
rz(-0.36205629) q[1];
sx q[1];
rz(1.9777745) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4951492) q[0];
sx q[0];
rz(-2.4730485) q[0];
sx q[0];
rz(0.77495452) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6357119) q[2];
sx q[2];
rz(-1.9569279) q[2];
sx q[2];
rz(-2.1979111) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.5608985) q[1];
sx q[1];
rz(-0.95798641) q[1];
sx q[1];
rz(1.8263019) q[1];
rz(1.1436361) q[3];
sx q[3];
rz(-1.8619976) q[3];
sx q[3];
rz(-1.8607651) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.290648) q[2];
sx q[2];
rz(-3.1182351) q[2];
sx q[2];
rz(1.5604431) q[2];
rz(-0.23049878) q[3];
sx q[3];
rz(-1.0483619) q[3];
sx q[3];
rz(0.44052625) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5480963) q[0];
sx q[0];
rz(-0.31262147) q[0];
sx q[0];
rz(3.0189959) q[0];
rz(0.7971881) q[1];
sx q[1];
rz(-1.0120579) q[1];
sx q[1];
rz(-0.0880934) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0038550475) q[0];
sx q[0];
rz(-2.3311576) q[0];
sx q[0];
rz(1.9268981) q[0];
rz(-pi) q[1];
rz(-1.2253907) q[2];
sx q[2];
rz(-0.24125762) q[2];
sx q[2];
rz(3.0221107) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.28411394) q[1];
sx q[1];
rz(-1.4339674) q[1];
sx q[1];
rz(-1.72422) q[1];
rz(0.54666211) q[3];
sx q[3];
rz(-0.61184363) q[3];
sx q[3];
rz(0.33239588) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.88283551) q[2];
sx q[2];
rz(-1.5828524) q[2];
sx q[2];
rz(-1.3239512) q[2];
rz(-2.6435408) q[3];
sx q[3];
rz(-0.96314722) q[3];
sx q[3];
rz(2.3948885) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28904706) q[0];
sx q[0];
rz(-2.9883224) q[0];
sx q[0];
rz(-2.9943941) q[0];
rz(-0.64951253) q[1];
sx q[1];
rz(-1.5107061) q[1];
sx q[1];
rz(1.46896) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.29028673) q[0];
sx q[0];
rz(-1.3477579) q[0];
sx q[0];
rz(-0.0073435535) q[0];
rz(-pi) q[1];
rz(2.9584453) q[2];
sx q[2];
rz(-1.1467271) q[2];
sx q[2];
rz(-1.9522425) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.6522347) q[1];
sx q[1];
rz(-1.918004) q[1];
sx q[1];
rz(2.2330797) q[1];
x q[2];
rz(0.62375237) q[3];
sx q[3];
rz(-2.459708) q[3];
sx q[3];
rz(0.15450134) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.2088251) q[2];
sx q[2];
rz(-2.6471477) q[2];
sx q[2];
rz(0.24173582) q[2];
rz(-0.73511165) q[3];
sx q[3];
rz(-1.7416411) q[3];
sx q[3];
rz(-0.66489768) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6786137) q[0];
sx q[0];
rz(-2.0942056) q[0];
sx q[0];
rz(0.12272923) q[0];
rz(-2.2562476) q[1];
sx q[1];
rz(-1.0095936) q[1];
sx q[1];
rz(-0.56040323) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3353676) q[0];
sx q[0];
rz(-1.0606375) q[0];
sx q[0];
rz(-0.2038184) q[0];
x q[1];
rz(-0.29717314) q[2];
sx q[2];
rz(-1.3577596) q[2];
sx q[2];
rz(-1.7293255) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.082277) q[1];
sx q[1];
rz(-2.6059285) q[1];
sx q[1];
rz(-2.0945063) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.4462279) q[3];
sx q[3];
rz(-2.0110907) q[3];
sx q[3];
rz(1.4780731) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.6178599) q[2];
sx q[2];
rz(-0.79404074) q[2];
sx q[2];
rz(0.20475556) q[2];
rz(2.2023885) q[3];
sx q[3];
rz(-1.771628) q[3];
sx q[3];
rz(-0.088976629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.27696779) q[0];
sx q[0];
rz(-0.11700103) q[0];
sx q[0];
rz(-0.65698874) q[0];
rz(-0.14343801) q[1];
sx q[1];
rz(-1.5564432) q[1];
sx q[1];
rz(-2.4030446) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6929106) q[0];
sx q[0];
rz(-1.4553207) q[0];
sx q[0];
rz(2.4666952) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0792909) q[2];
sx q[2];
rz(-1.3710944) q[2];
sx q[2];
rz(-2.8815567) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.025102928) q[1];
sx q[1];
rz(-1.2560351) q[1];
sx q[1];
rz(0.88309755) q[1];
rz(-1.012771) q[3];
sx q[3];
rz(-0.9415365) q[3];
sx q[3];
rz(-2.2037639) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.8831545) q[2];
sx q[2];
rz(-0.65885764) q[2];
sx q[2];
rz(0.0030041791) q[2];
rz(-0.7891807) q[3];
sx q[3];
rz(-0.3843669) q[3];
sx q[3];
rz(1.9743617) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.060519144) q[0];
sx q[0];
rz(-2.7923212) q[0];
sx q[0];
rz(0.14807598) q[0];
rz(0.5332467) q[1];
sx q[1];
rz(-2.8956469) q[1];
sx q[1];
rz(1.5124793) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0740967) q[0];
sx q[0];
rz(-1.3036911) q[0];
sx q[0];
rz(-0.47988923) q[0];
rz(-0.92502906) q[2];
sx q[2];
rz(-0.86849493) q[2];
sx q[2];
rz(1.4214732) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.80485) q[1];
sx q[1];
rz(-2.7798493) q[1];
sx q[1];
rz(1.9182349) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3324206) q[3];
sx q[3];
rz(-2.4052103) q[3];
sx q[3];
rz(0.045698085) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.6988354) q[2];
sx q[2];
rz(-1.7008702) q[2];
sx q[2];
rz(-0.093078144) q[2];
rz(-0.062151521) q[3];
sx q[3];
rz(-2.744894) q[3];
sx q[3];
rz(1.2232346) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1208948) q[0];
sx q[0];
rz(-2.5466205) q[0];
sx q[0];
rz(0.66463941) q[0];
rz(1.7918034) q[1];
sx q[1];
rz(-2.2638958) q[1];
sx q[1];
rz(0.68914366) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1375547) q[0];
sx q[0];
rz(-1.327164) q[0];
sx q[0];
rz(1.451127) q[0];
rz(0.87920636) q[2];
sx q[2];
rz(-1.9636646) q[2];
sx q[2];
rz(2.526432) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-3.0943567) q[1];
sx q[1];
rz(-2.7942889) q[1];
sx q[1];
rz(-0.44012614) q[1];
x q[2];
rz(2.7921121) q[3];
sx q[3];
rz(-0.17519028) q[3];
sx q[3];
rz(-2.3153806) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.050921116) q[2];
sx q[2];
rz(-0.65616578) q[2];
sx q[2];
rz(1.3760759) q[2];
rz(1.225166) q[3];
sx q[3];
rz(-0.71198946) q[3];
sx q[3];
rz(-0.041697748) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0390778) q[0];
sx q[0];
rz(-0.044487655) q[0];
sx q[0];
rz(-2.5503889) q[0];
rz(2.8996331) q[1];
sx q[1];
rz(-1.0958902) q[1];
sx q[1];
rz(0.42344365) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5862783) q[0];
sx q[0];
rz(-1.8722157) q[0];
sx q[0];
rz(-2.8239646) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9482163) q[2];
sx q[2];
rz(-0.49623734) q[2];
sx q[2];
rz(1.9691182) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.6953227) q[1];
sx q[1];
rz(-2.5188418) q[1];
sx q[1];
rz(2.0957309) q[1];
x q[2];
rz(-0.54355346) q[3];
sx q[3];
rz(-2.4494236) q[3];
sx q[3];
rz(-2.8080733) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.55687904) q[2];
sx q[2];
rz(-0.60606474) q[2];
sx q[2];
rz(1.1663743) q[2];
rz(-1.1276468) q[3];
sx q[3];
rz(-2.1731264) q[3];
sx q[3];
rz(-0.52496547) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.128085) q[0];
sx q[0];
rz(-2.7078244) q[0];
sx q[0];
rz(-0.67310131) q[0];
rz(-0.85064864) q[1];
sx q[1];
rz(-1.7885845) q[1];
sx q[1];
rz(-2.8180715) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.464768) q[0];
sx q[0];
rz(-0.45545211) q[0];
sx q[0];
rz(-2.1663329) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.27574972) q[2];
sx q[2];
rz(-1.1080896) q[2];
sx q[2];
rz(-1.5327061) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.4648351) q[1];
sx q[1];
rz(-0.82993648) q[1];
sx q[1];
rz(2.6840032) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5786641) q[3];
sx q[3];
rz(-2.1546225) q[3];
sx q[3];
rz(1.5848094) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.9234151) q[2];
sx q[2];
rz(-2.9604993) q[2];
sx q[2];
rz(0.072048135) q[2];
rz(-2.1273023) q[3];
sx q[3];
rz(-0.91599661) q[3];
sx q[3];
rz(-0.54916507) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2608248) q[0];
sx q[0];
rz(-1.4243955) q[0];
sx q[0];
rz(-2.0212174) q[0];
rz(-2.8063759) q[1];
sx q[1];
rz(-0.64246476) q[1];
sx q[1];
rz(-1.1186218) q[1];
rz(0.14090385) q[2];
sx q[2];
rz(-2.0456516) q[2];
sx q[2];
rz(2.4608154) q[2];
rz(2.9670197) q[3];
sx q[3];
rz(-2.9963507) q[3];
sx q[3];
rz(0.33490845) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
