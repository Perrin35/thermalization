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
rz(0.53606755) q[0];
sx q[0];
rz(-1.9789088) q[0];
sx q[0];
rz(-0.43098488) q[0];
rz(-5.9302063) q[1];
sx q[1];
rz(7.5063345) q[1];
sx q[1];
rz(5.8533668) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0551222) q[0];
sx q[0];
rz(-1.4807379) q[0];
sx q[0];
rz(-0.49205972) q[0];
rz(-2.0031646) q[2];
sx q[2];
rz(-1.5562487) q[2];
sx q[2];
rz(2.5231139) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.7215772) q[1];
sx q[1];
rz(-0.41999451) q[1];
sx q[1];
rz(2.6073522) q[1];
rz(-pi) q[2];
rz(0.63685959) q[3];
sx q[3];
rz(-0.59266337) q[3];
sx q[3];
rz(-0.67837472) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.063804403) q[2];
sx q[2];
rz(-0.52738515) q[2];
sx q[2];
rz(1.0112313) q[2];
rz(-1.0412591) q[3];
sx q[3];
rz(-2.6285089) q[3];
sx q[3];
rz(-0.43963715) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1083199) q[0];
sx q[0];
rz(-2.1450873) q[0];
sx q[0];
rz(-2.394115) q[0];
rz(-0.29391089) q[1];
sx q[1];
rz(-1.1312609) q[1];
sx q[1];
rz(-2.8864313) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.19758148) q[0];
sx q[0];
rz(-1.4048049) q[0];
sx q[0];
rz(-2.152488) q[0];
rz(-pi) q[1];
rz(1.3841278) q[2];
sx q[2];
rz(-1.871284) q[2];
sx q[2];
rz(1.0050424) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.0415062) q[1];
sx q[1];
rz(-1.1250682) q[1];
sx q[1];
rz(-2.2959832) q[1];
rz(-pi) q[2];
rz(2.1952755) q[3];
sx q[3];
rz(-1.7229652) q[3];
sx q[3];
rz(2.5077903) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.1110288) q[2];
sx q[2];
rz(-3.0869637) q[2];
sx q[2];
rz(2.2500706) q[2];
rz(0.24957481) q[3];
sx q[3];
rz(-2.2587743) q[3];
sx q[3];
rz(-0.3961302) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4075277) q[0];
sx q[0];
rz(-2.0052818) q[0];
sx q[0];
rz(3.0974645) q[0];
rz(-1.2491501) q[1];
sx q[1];
rz(-0.42611486) q[1];
sx q[1];
rz(2.6557907) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7120948) q[0];
sx q[0];
rz(-2.9154239) q[0];
sx q[0];
rz(0.087271376) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2052231) q[2];
sx q[2];
rz(-0.23072019) q[2];
sx q[2];
rz(1.7796734) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.9456182) q[1];
sx q[1];
rz(-0.83688078) q[1];
sx q[1];
rz(2.3952146) q[1];
rz(-2.3084675) q[3];
sx q[3];
rz(-0.55748788) q[3];
sx q[3];
rz(-1.1380029) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.8140063) q[2];
sx q[2];
rz(-2.0531211) q[2];
sx q[2];
rz(-2.4468454) q[2];
rz(-2.5658549) q[3];
sx q[3];
rz(-1.8703987) q[3];
sx q[3];
rz(-1.7339138) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(-2.1908252) q[0];
sx q[0];
rz(-0.29470834) q[0];
sx q[0];
rz(0.29275352) q[0];
rz(2.5726908) q[1];
sx q[1];
rz(-2.3141373) q[1];
sx q[1];
rz(-0.84691602) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6349061) q[0];
sx q[0];
rz(-1.4932079) q[0];
sx q[0];
rz(-1.9798748) q[0];
rz(-3.0693502) q[2];
sx q[2];
rz(-0.90257593) q[2];
sx q[2];
rz(-0.39230686) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.9255815) q[1];
sx q[1];
rz(-0.30864172) q[1];
sx q[1];
rz(0.68283783) q[1];
x q[2];
rz(-1.9027223) q[3];
sx q[3];
rz(-0.20461018) q[3];
sx q[3];
rz(2.0505035) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.7920821) q[2];
sx q[2];
rz(-1.8831207) q[2];
sx q[2];
rz(0.2087896) q[2];
rz(2.4236692) q[3];
sx q[3];
rz(-0.28306285) q[3];
sx q[3];
rz(2.2010402) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1013041) q[0];
sx q[0];
rz(-0.52241075) q[0];
sx q[0];
rz(-2.8029602) q[0];
rz(0.70308095) q[1];
sx q[1];
rz(-0.73892006) q[1];
sx q[1];
rz(-0.83831659) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.188952) q[0];
sx q[0];
rz(-1.8243655) q[0];
sx q[0];
rz(2.6978561) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9432326) q[2];
sx q[2];
rz(-2.3075929) q[2];
sx q[2];
rz(-1.5892346) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.053239659) q[1];
sx q[1];
rz(-0.87118282) q[1];
sx q[1];
rz(3.0979473) q[1];
x q[2];
rz(-2.5637066) q[3];
sx q[3];
rz(-0.69212038) q[3];
sx q[3];
rz(-2.1731165) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.4579953) q[2];
sx q[2];
rz(-1.0883051) q[2];
sx q[2];
rz(-1.0028769) q[2];
rz(-0.14497997) q[3];
sx q[3];
rz(-3.0478015) q[3];
sx q[3];
rz(-3.0785479) q[3];
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
rz(1.018464) q[0];
sx q[0];
rz(-0.58582425) q[0];
sx q[0];
rz(2.1783094) q[0];
rz(2.9604498) q[1];
sx q[1];
rz(-0.89814848) q[1];
sx q[1];
rz(-1.9021665) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2116263) q[0];
sx q[0];
rz(-0.74290448) q[0];
sx q[0];
rz(0.56908281) q[0];
x q[1];
rz(-0.507429) q[2];
sx q[2];
rz(-2.6118738) q[2];
sx q[2];
rz(1.5629753) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.9037348) q[1];
sx q[1];
rz(-1.1958201) q[1];
sx q[1];
rz(1.9720248) q[1];
rz(-pi) q[2];
rz(-1.7224947) q[3];
sx q[3];
rz(-1.6983508) q[3];
sx q[3];
rz(-0.90326819) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.2565101) q[2];
sx q[2];
rz(-2.2233267) q[2];
sx q[2];
rz(2.2086842) q[2];
rz(1.4126623) q[3];
sx q[3];
rz(-1.4224854) q[3];
sx q[3];
rz(-0.2093813) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4466208) q[0];
sx q[0];
rz(-0.3681204) q[0];
sx q[0];
rz(2.3387261) q[0];
rz(2.39957) q[1];
sx q[1];
rz(-1.4890393) q[1];
sx q[1];
rz(-0.41025695) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8703394) q[0];
sx q[0];
rz(-1.2294719) q[0];
sx q[0];
rz(2.5786311) q[0];
rz(-pi) q[1];
rz(2.3773277) q[2];
sx q[2];
rz(-1.6558803) q[2];
sx q[2];
rz(-2.6972023) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.8852119) q[1];
sx q[1];
rz(-1.7325914) q[1];
sx q[1];
rz(0.075116861) q[1];
x q[2];
rz(-1.0954082) q[3];
sx q[3];
rz(-0.96964004) q[3];
sx q[3];
rz(0.22829311) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.86847574) q[2];
sx q[2];
rz(-1.5584402) q[2];
sx q[2];
rz(-2.9203019) q[2];
rz(2.7224329) q[3];
sx q[3];
rz(-0.48210382) q[3];
sx q[3];
rz(-2.6656849) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.23685037) q[0];
sx q[0];
rz(-2.1433647) q[0];
sx q[0];
rz(1.6118443) q[0];
rz(-1.8086241) q[1];
sx q[1];
rz(-2.1905441) q[1];
sx q[1];
rz(-1.453368) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26448956) q[0];
sx q[0];
rz(-2.1801688) q[0];
sx q[0];
rz(2.2731056) q[0];
x q[1];
rz(-2.9506893) q[2];
sx q[2];
rz(-1.4530572) q[2];
sx q[2];
rz(0.28798739) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.61078) q[1];
sx q[1];
rz(-1.4207834) q[1];
sx q[1];
rz(0.71041815) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2923041) q[3];
sx q[3];
rz(-1.7207816) q[3];
sx q[3];
rz(-1.5453135) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.86614418) q[2];
sx q[2];
rz(-0.25442213) q[2];
sx q[2];
rz(-0.010738372) q[2];
rz(2.8209316) q[3];
sx q[3];
rz(-1.843822) q[3];
sx q[3];
rz(2.5519154) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7084932) q[0];
sx q[0];
rz(-0.95516094) q[0];
sx q[0];
rz(-2.6668715) q[0];
rz(-0.2933329) q[1];
sx q[1];
rz(-1.1055929) q[1];
sx q[1];
rz(-1.6620103) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.06874456) q[0];
sx q[0];
rz(-1.5340051) q[0];
sx q[0];
rz(-1.6108914) q[0];
rz(-pi) q[1];
rz(-2.4655645) q[2];
sx q[2];
rz(-1.3814298) q[2];
sx q[2];
rz(0.87847906) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.41614) q[1];
sx q[1];
rz(-1.801277) q[1];
sx q[1];
rz(-2.6621755) q[1];
rz(-pi) q[2];
rz(2.9607356) q[3];
sx q[3];
rz(-1.81362) q[3];
sx q[3];
rz(-2.9202785) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.097229615) q[2];
sx q[2];
rz(-1.2791415) q[2];
sx q[2];
rz(-2.7105159) q[2];
rz(2.9512682) q[3];
sx q[3];
rz(-0.35076916) q[3];
sx q[3];
rz(-0.88461191) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1650319) q[0];
sx q[0];
rz(-2.850552) q[0];
sx q[0];
rz(0.2188368) q[0];
rz(0.95895514) q[1];
sx q[1];
rz(-0.7364277) q[1];
sx q[1];
rz(-1.7494019) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9746124) q[0];
sx q[0];
rz(-1.3777527) q[0];
sx q[0];
rz(-1.6231322) q[0];
x q[1];
rz(-1.4264246) q[2];
sx q[2];
rz(-1.9871658) q[2];
sx q[2];
rz(1.8954111) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.2710378) q[1];
sx q[1];
rz(-0.95083323) q[1];
sx q[1];
rz(1.2024173) q[1];
rz(-pi) q[2];
x q[2];
rz(0.43057118) q[3];
sx q[3];
rz(-2.1408251) q[3];
sx q[3];
rz(1.7209297) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.86768156) q[2];
sx q[2];
rz(-2.3751004) q[2];
sx q[2];
rz(-1.8217746) q[2];
rz(2.6093318) q[3];
sx q[3];
rz(-1.7936734) q[3];
sx q[3];
rz(-1.5490612) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5821447) q[0];
sx q[0];
rz(-1.5912709) q[0];
sx q[0];
rz(-2.7406319) q[0];
rz(2.9136912) q[1];
sx q[1];
rz(-1.0453929) q[1];
sx q[1];
rz(1.4452404) q[1];
rz(-0.41517045) q[2];
sx q[2];
rz(-1.8098469) q[2];
sx q[2];
rz(-1.3630223) q[2];
rz(2.1501272) q[3];
sx q[3];
rz(-1.7064863) q[3];
sx q[3];
rz(1.9319921) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
