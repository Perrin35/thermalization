OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.6634231) q[0];
sx q[0];
rz(-2.2278251) q[0];
sx q[0];
rz(-1.0650286) q[0];
rz(1.6051259) q[1];
sx q[1];
rz(-1.1054339) q[1];
sx q[1];
rz(0.096535834) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8138344) q[0];
sx q[0];
rz(-1.3655435) q[0];
sx q[0];
rz(1.1428409) q[0];
rz(-pi) q[1];
rz(-0.41806721) q[2];
sx q[2];
rz(-2.5804248) q[2];
sx q[2];
rz(0.18218606) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.9241544) q[1];
sx q[1];
rz(-0.70637843) q[1];
sx q[1];
rz(-0.46831727) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4289342) q[3];
sx q[3];
rz(-2.5569698) q[3];
sx q[3];
rz(-1.4277924) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.0194633) q[2];
sx q[2];
rz(-0.18522842) q[2];
sx q[2];
rz(-1.8786001) q[2];
rz(2.7428135) q[3];
sx q[3];
rz(-2.0006657) q[3];
sx q[3];
rz(-1.0293452) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1102092) q[0];
sx q[0];
rz(-0.6321913) q[0];
sx q[0];
rz(1.6803918) q[0];
rz(-0.070143135) q[1];
sx q[1];
rz(-1.5927916) q[1];
sx q[1];
rz(-2.7102914) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3628497) q[0];
sx q[0];
rz(-1.6322789) q[0];
sx q[0];
rz(-0.009470609) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8441418) q[2];
sx q[2];
rz(-1.97142) q[2];
sx q[2];
rz(1.9971776) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.13557735) q[1];
sx q[1];
rz(-0.90432465) q[1];
sx q[1];
rz(2.2462387) q[1];
rz(-pi) q[2];
rz(2.7106695) q[3];
sx q[3];
rz(-0.049073372) q[3];
sx q[3];
rz(-1.8762833) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.18616072) q[2];
sx q[2];
rz(-1.7308851) q[2];
sx q[2];
rz(0.53039941) q[2];
rz(2.0576599) q[3];
sx q[3];
rz(-2.0518186) q[3];
sx q[3];
rz(-1.8918022) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0386061) q[0];
sx q[0];
rz(-2.5654721) q[0];
sx q[0];
rz(1.9392133) q[0];
rz(-2.1309958) q[1];
sx q[1];
rz(-1.3253515) q[1];
sx q[1];
rz(0.030166322) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.3836244) q[0];
sx q[0];
rz(-0.080648184) q[0];
sx q[0];
rz(-0.6225125) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.518345) q[2];
sx q[2];
rz(-1.5546335) q[2];
sx q[2];
rz(-2.9132622) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.846993) q[1];
sx q[1];
rz(-1.9140118) q[1];
sx q[1];
rz(-1.4354475) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4143589) q[3];
sx q[3];
rz(-0.57830252) q[3];
sx q[3];
rz(2.3012637) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.5304337) q[2];
sx q[2];
rz(-2.9035089) q[2];
sx q[2];
rz(-2.7670624) q[2];
rz(1.1206867) q[3];
sx q[3];
rz(-1.6569542) q[3];
sx q[3];
rz(0.69168958) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.94367868) q[0];
sx q[0];
rz(-1.8489842) q[0];
sx q[0];
rz(-1.9669272) q[0];
rz(1.3150175) q[1];
sx q[1];
rz(-2.0349793) q[1];
sx q[1];
rz(0.16474251) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.574802) q[0];
sx q[0];
rz(-1.0079449) q[0];
sx q[0];
rz(-0.99250162) q[0];
rz(-pi) q[1];
rz(-2.6283662) q[2];
sx q[2];
rz(-1.7145707) q[2];
sx q[2];
rz(-1.9699186) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.10197266) q[1];
sx q[1];
rz(-1.0821187) q[1];
sx q[1];
rz(-0.38445977) q[1];
rz(0.03892558) q[3];
sx q[3];
rz(-1.4994326) q[3];
sx q[3];
rz(-1.7046649) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.4508007) q[2];
sx q[2];
rz(-1.5361293) q[2];
sx q[2];
rz(-0.50626051) q[2];
rz(2.6311503) q[3];
sx q[3];
rz(-2.3531395) q[3];
sx q[3];
rz(2.3300664) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6798169) q[0];
sx q[0];
rz(-0.33648574) q[0];
sx q[0];
rz(0.80999723) q[0];
rz(-3.0432155) q[1];
sx q[1];
rz(-0.38946238) q[1];
sx q[1];
rz(2.817165) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2330889) q[0];
sx q[0];
rz(-1.5727186) q[0];
sx q[0];
rz(-1.6932827) q[0];
x q[1];
rz(2.5971707) q[2];
sx q[2];
rz(-2.8834497) q[2];
sx q[2];
rz(-0.38636696) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.1871029) q[1];
sx q[1];
rz(-1.3924161) q[1];
sx q[1];
rz(2.5599077) q[1];
rz(0.47629997) q[3];
sx q[3];
rz(-1.7576964) q[3];
sx q[3];
rz(-0.66038654) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.7847932) q[2];
sx q[2];
rz(-1.5609799) q[2];
sx q[2];
rz(2.5830833) q[2];
rz(-1.2950581) q[3];
sx q[3];
rz(-1.7241155) q[3];
sx q[3];
rz(-2.7891125) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2297939) q[0];
sx q[0];
rz(-0.42775446) q[0];
sx q[0];
rz(1.3707004) q[0];
rz(1.7478583) q[1];
sx q[1];
rz(-2.1262028) q[1];
sx q[1];
rz(-0.14399354) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32329473) q[0];
sx q[0];
rz(-0.12579623) q[0];
sx q[0];
rz(-0.97387858) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1156111) q[2];
sx q[2];
rz(-1.7894288) q[2];
sx q[2];
rz(-1.3022547) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.2768133) q[1];
sx q[1];
rz(-2.1739156) q[1];
sx q[1];
rz(0.13637194) q[1];
rz(-pi) q[2];
rz(2.6089588) q[3];
sx q[3];
rz(-0.71511474) q[3];
sx q[3];
rz(-1.2330513) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.41649848) q[2];
sx q[2];
rz(-0.35606733) q[2];
sx q[2];
rz(2.1605675) q[2];
rz(2.6128926) q[3];
sx q[3];
rz(-1.7873535) q[3];
sx q[3];
rz(-2.5813812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.33787802) q[0];
sx q[0];
rz(-2.3382472) q[0];
sx q[0];
rz(-1.2038318) q[0];
rz(-3.1375258) q[1];
sx q[1];
rz(-1.4264359) q[1];
sx q[1];
rz(2.8003661) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1004286) q[0];
sx q[0];
rz(-2.251314) q[0];
sx q[0];
rz(-1.0170223) q[0];
x q[1];
rz(1.5725102) q[2];
sx q[2];
rz(-0.69134635) q[2];
sx q[2];
rz(-1.5087663) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.8123521) q[1];
sx q[1];
rz(-1.7197548) q[1];
sx q[1];
rz(-1.5801058) q[1];
rz(-pi) q[2];
rz(-2.7507095) q[3];
sx q[3];
rz(-2.7733186) q[3];
sx q[3];
rz(-2.0892909) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.61871201) q[2];
sx q[2];
rz(-1.2423542) q[2];
sx q[2];
rz(0.27113554) q[2];
rz(1.6297657) q[3];
sx q[3];
rz(-2.1745493) q[3];
sx q[3];
rz(-2.7395774) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4790344) q[0];
sx q[0];
rz(-2.5738578) q[0];
sx q[0];
rz(-2.973279) q[0];
rz(1.0100826) q[1];
sx q[1];
rz(-1.964566) q[1];
sx q[1];
rz(-3.041306) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1948457) q[0];
sx q[0];
rz(-1.2887717) q[0];
sx q[0];
rz(2.8213324) q[0];
rz(-pi) q[1];
rz(0.73302539) q[2];
sx q[2];
rz(-1.4372281) q[2];
sx q[2];
rz(2.095993) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.6319233) q[1];
sx q[1];
rz(-0.20317101) q[1];
sx q[1];
rz(-0.12170273) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1264621) q[3];
sx q[3];
rz(-2.0499886) q[3];
sx q[3];
rz(1.5492619) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.4526796) q[2];
sx q[2];
rz(-2.6426688) q[2];
sx q[2];
rz(-1.720724) q[2];
rz(1.0639327) q[3];
sx q[3];
rz(-1.7199793) q[3];
sx q[3];
rz(-0.091778278) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38561472) q[0];
sx q[0];
rz(-1.6843963) q[0];
sx q[0];
rz(-3.1047367) q[0];
rz(-1.256975) q[1];
sx q[1];
rz(-1.5024065) q[1];
sx q[1];
rz(1.771155) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.12808558) q[0];
sx q[0];
rz(-2.8261746) q[0];
sx q[0];
rz(2.0924139) q[0];
rz(-pi) q[1];
rz(1.7212935) q[2];
sx q[2];
rz(-1.2062819) q[2];
sx q[2];
rz(0.30752674) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.47342604) q[1];
sx q[1];
rz(-1.2915732) q[1];
sx q[1];
rz(0.88069005) q[1];
rz(-pi) q[2];
rz(-1.4716205) q[3];
sx q[3];
rz(-0.41891019) q[3];
sx q[3];
rz(0.3813972) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.86159414) q[2];
sx q[2];
rz(-1.7656606) q[2];
sx q[2];
rz(0.17246788) q[2];
rz(2.5837768) q[3];
sx q[3];
rz(-0.54111257) q[3];
sx q[3];
rz(-1.4403758) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7114792) q[0];
sx q[0];
rz(-2.5028296) q[0];
sx q[0];
rz(0.3399671) q[0];
rz(-2.4760447) q[1];
sx q[1];
rz(-1.7219209) q[1];
sx q[1];
rz(1.3131712) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42178828) q[0];
sx q[0];
rz(-1.3911085) q[0];
sx q[0];
rz(-1.6186886) q[0];
rz(-0.27965172) q[2];
sx q[2];
rz(-2.8707829) q[2];
sx q[2];
rz(-0.93284399) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.0413805) q[1];
sx q[1];
rz(-1.3330701) q[1];
sx q[1];
rz(-1.8820534) q[1];
rz(0.35452133) q[3];
sx q[3];
rz(-0.42663048) q[3];
sx q[3];
rz(1.420295) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.36249915) q[2];
sx q[2];
rz(-0.43467793) q[2];
sx q[2];
rz(-2.4707826) q[2];
rz(0.32495156) q[3];
sx q[3];
rz(-2.3281125) q[3];
sx q[3];
rz(-2.1285098) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.308607) q[0];
sx q[0];
rz(-1.8600464) q[0];
sx q[0];
rz(-1.6639584) q[0];
rz(-1.0302522) q[1];
sx q[1];
rz(-1.4719084) q[1];
sx q[1];
rz(1.7869064) q[1];
rz(-0.5253039) q[2];
sx q[2];
rz(-1.1266194) q[2];
sx q[2];
rz(2.5818679) q[2];
rz(0.58086953) q[3];
sx q[3];
rz(-1.8884251) q[3];
sx q[3];
rz(-2.6407218) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
