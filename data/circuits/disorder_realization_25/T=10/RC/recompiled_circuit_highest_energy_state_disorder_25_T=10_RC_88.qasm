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
rz(-2.146848) q[0];
sx q[0];
rz(-0.58726197) q[0];
sx q[0];
rz(-0.62556148) q[0];
rz(0.62043959) q[1];
sx q[1];
rz(-2.151139) q[1];
sx q[1];
rz(0.78392309) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1358937) q[0];
sx q[0];
rz(-1.5672475) q[0];
sx q[0];
rz(2.7074642) q[0];
rz(0.13909362) q[2];
sx q[2];
rz(-1.4882097) q[2];
sx q[2];
rz(2.6001584) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.16456024) q[1];
sx q[1];
rz(-0.98233084) q[1];
sx q[1];
rz(2.2684437) q[1];
rz(-pi) q[2];
rz(0.089139912) q[3];
sx q[3];
rz(-0.77971346) q[3];
sx q[3];
rz(2.4453795) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.1938532) q[2];
sx q[2];
rz(-1.0449907) q[2];
sx q[2];
rz(2.4742773) q[2];
rz(0.31467485) q[3];
sx q[3];
rz(-0.59942013) q[3];
sx q[3];
rz(0.92710322) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.84745234) q[0];
sx q[0];
rz(-2.2951234) q[0];
sx q[0];
rz(2.8417929) q[0];
rz(-2.1369797) q[1];
sx q[1];
rz(-2.424898) q[1];
sx q[1];
rz(-2.2915548) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3687821) q[0];
sx q[0];
rz(-1.4748992) q[0];
sx q[0];
rz(0.96536388) q[0];
rz(-pi) q[1];
rz(-2.7072435) q[2];
sx q[2];
rz(-2.953989) q[2];
sx q[2];
rz(1.6845195) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.1239224) q[1];
sx q[1];
rz(-0.41819376) q[1];
sx q[1];
rz(1.0355509) q[1];
rz(-1.5845234) q[3];
sx q[3];
rz(-1.1662332) q[3];
sx q[3];
rz(1.3409529) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.5726418) q[2];
sx q[2];
rz(-2.6093542) q[2];
sx q[2];
rz(-2.4786095) q[2];
rz(-0.7524544) q[3];
sx q[3];
rz(-1.3198676) q[3];
sx q[3];
rz(-2.9451356) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3982584) q[0];
sx q[0];
rz(-2.7404116) q[0];
sx q[0];
rz(2.4617526) q[0];
rz(2.4196449) q[1];
sx q[1];
rz(-0.55689055) q[1];
sx q[1];
rz(-0.78071761) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3328898) q[0];
sx q[0];
rz(-1.5631544) q[0];
sx q[0];
rz(1.1701258) q[0];
x q[1];
rz(-1.6569312) q[2];
sx q[2];
rz(-1.5670098) q[2];
sx q[2];
rz(-2.6876793) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.5131849) q[1];
sx q[1];
rz(-0.76090136) q[1];
sx q[1];
rz(0.55879467) q[1];
rz(-1.3296919) q[3];
sx q[3];
rz(-0.74987292) q[3];
sx q[3];
rz(2.1151224) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.9146933) q[2];
sx q[2];
rz(-1.2279899) q[2];
sx q[2];
rz(-2.4719888) q[2];
rz(2.2412444) q[3];
sx q[3];
rz(-2.1293631) q[3];
sx q[3];
rz(-2.3646234) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25903073) q[0];
sx q[0];
rz(-2.7029523) q[0];
sx q[0];
rz(2.2341109) q[0];
rz(0.59440815) q[1];
sx q[1];
rz(-2.2120357) q[1];
sx q[1];
rz(-2.0118735) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.27760273) q[0];
sx q[0];
rz(-0.66501319) q[0];
sx q[0];
rz(1.8825085) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9920861) q[2];
sx q[2];
rz(-2.0975916) q[2];
sx q[2];
rz(3.0164928) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.4817244) q[1];
sx q[1];
rz(-0.65978434) q[1];
sx q[1];
rz(-0.45141141) q[1];
rz(1.8640737) q[3];
sx q[3];
rz(-1.8403111) q[3];
sx q[3];
rz(-2.6285009) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(3.0323459) q[2];
sx q[2];
rz(-1.3999908) q[2];
sx q[2];
rz(0.67231154) q[2];
rz(0.0191056) q[3];
sx q[3];
rz(-2.8002383) q[3];
sx q[3];
rz(-0.13489558) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8371209) q[0];
sx q[0];
rz(-2.3470375) q[0];
sx q[0];
rz(-2.9735612) q[0];
rz(1.9145603) q[1];
sx q[1];
rz(-1.9214168) q[1];
sx q[1];
rz(-0.29774818) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3963602) q[0];
sx q[0];
rz(-1.0567291) q[0];
sx q[0];
rz(-0.6020418) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5644659) q[2];
sx q[2];
rz(-1.071605) q[2];
sx q[2];
rz(1.2291723) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.6444466) q[1];
sx q[1];
rz(-2.9342817) q[1];
sx q[1];
rz(-0.12122341) q[1];
rz(-pi) q[2];
rz(-1.4120123) q[3];
sx q[3];
rz(-1.5701236) q[3];
sx q[3];
rz(1.0345881) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.1905404) q[2];
sx q[2];
rz(-1.0995883) q[2];
sx q[2];
rz(-0.61868787) q[2];
rz(-0.80858532) q[3];
sx q[3];
rz(-0.4777258) q[3];
sx q[3];
rz(1.8019069) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.488778) q[0];
sx q[0];
rz(-1.6365016) q[0];
sx q[0];
rz(0.46522796) q[0];
rz(-0.48509994) q[1];
sx q[1];
rz(-2.7310889) q[1];
sx q[1];
rz(3.0533275) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.34644914) q[0];
sx q[0];
rz(-0.69899054) q[0];
sx q[0];
rz(-0.82620746) q[0];
rz(-pi) q[1];
rz(-0.012862269) q[2];
sx q[2];
rz(-0.84190997) q[2];
sx q[2];
rz(-2.5358605) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.7317071) q[1];
sx q[1];
rz(-2.3165858) q[1];
sx q[1];
rz(-1.7959926) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.29084716) q[3];
sx q[3];
rz(-2.8418782) q[3];
sx q[3];
rz(2.7951473) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.7469067) q[2];
sx q[2];
rz(-2.2862819) q[2];
sx q[2];
rz(-2.6439903) q[2];
rz(0.48007128) q[3];
sx q[3];
rz(-2.2205455) q[3];
sx q[3];
rz(3.0729496) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8467872) q[0];
sx q[0];
rz(-2.6755896) q[0];
sx q[0];
rz(-1.4303327) q[0];
rz(-1.3149186) q[1];
sx q[1];
rz(-2.7480835) q[1];
sx q[1];
rz(1.5865145) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.67051149) q[0];
sx q[0];
rz(-1.1789315) q[0];
sx q[0];
rz(-3.0125822) q[0];
x q[1];
rz(-2.2261593) q[2];
sx q[2];
rz(-0.829773) q[2];
sx q[2];
rz(0.96617736) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.3433229) q[1];
sx q[1];
rz(-0.99162356) q[1];
sx q[1];
rz(0.12076853) q[1];
rz(-pi) q[2];
x q[2];
rz(0.076818941) q[3];
sx q[3];
rz(-2.2008476) q[3];
sx q[3];
rz(2.1253642) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.2985208) q[2];
sx q[2];
rz(-2.4303747) q[2];
sx q[2];
rz(0.92740518) q[2];
rz(1.1585506) q[3];
sx q[3];
rz(-2.4535593) q[3];
sx q[3];
rz(0.096175171) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4226828) q[0];
sx q[0];
rz(-2.2490608) q[0];
sx q[0];
rz(-0.476015) q[0];
rz(-0.99126518) q[1];
sx q[1];
rz(-2.3920993) q[1];
sx q[1];
rz(-3.057726) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.80591969) q[0];
sx q[0];
rz(-1.085338) q[0];
sx q[0];
rz(-2.7341163) q[0];
rz(-pi) q[1];
x q[1];
rz(0.4617347) q[2];
sx q[2];
rz(-1.2424011) q[2];
sx q[2];
rz(-0.4730313) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.8955651) q[1];
sx q[1];
rz(-1.152134) q[1];
sx q[1];
rz(0.093105969) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7937035) q[3];
sx q[3];
rz(-0.37316868) q[3];
sx q[3];
rz(-0.56115967) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.91741651) q[2];
sx q[2];
rz(-2.8755964) q[2];
sx q[2];
rz(0.64845294) q[2];
rz(-0.90483856) q[3];
sx q[3];
rz(-1.4584533) q[3];
sx q[3];
rz(2.8308433) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
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
rz(2.7998578) q[0];
sx q[0];
rz(-0.71001995) q[0];
sx q[0];
rz(-2.572686) q[0];
rz(-0.83373249) q[1];
sx q[1];
rz(-1.8582452) q[1];
sx q[1];
rz(2.1368829) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6093401) q[0];
sx q[0];
rz(-2.1933317) q[0];
sx q[0];
rz(1.3329766) q[0];
x q[1];
rz(0.23022454) q[2];
sx q[2];
rz(-1.0448928) q[2];
sx q[2];
rz(-0.81947015) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.6697996) q[1];
sx q[1];
rz(-2.2679015) q[1];
sx q[1];
rz(-0.094273829) q[1];
rz(-pi) q[2];
rz(-1.4887722) q[3];
sx q[3];
rz(-0.98482212) q[3];
sx q[3];
rz(2.2524633) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.4363165) q[2];
sx q[2];
rz(-0.87203163) q[2];
sx q[2];
rz(-2.7936068) q[2];
rz(-2.3157388) q[3];
sx q[3];
rz(-0.21819849) q[3];
sx q[3];
rz(2.5364449) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0843049) q[0];
sx q[0];
rz(-2.075752) q[0];
sx q[0];
rz(2.9344946) q[0];
rz(0.53889489) q[1];
sx q[1];
rz(-2.5205044) q[1];
sx q[1];
rz(-2.5394687) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.83610632) q[0];
sx q[0];
rz(-2.0417111) q[0];
sx q[0];
rz(-0.54022809) q[0];
x q[1];
rz(-2.5536898) q[2];
sx q[2];
rz(-1.9155115) q[2];
sx q[2];
rz(1.6893139) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.6844897) q[1];
sx q[1];
rz(-2.7278363) q[1];
sx q[1];
rz(-1.7587508) q[1];
rz(-2.3942457) q[3];
sx q[3];
rz(-1.8149733) q[3];
sx q[3];
rz(1.7908733) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.3567317) q[2];
sx q[2];
rz(-0.95180231) q[2];
sx q[2];
rz(1.3447364) q[2];
rz(0.52062672) q[3];
sx q[3];
rz(-2.8713363) q[3];
sx q[3];
rz(2.3394913) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6945334) q[0];
sx q[0];
rz(-0.72493989) q[0];
sx q[0];
rz(-1.3865393) q[0];
rz(0.78320349) q[1];
sx q[1];
rz(-1.422311) q[1];
sx q[1];
rz(1.5625988) q[1];
rz(1.8405909) q[2];
sx q[2];
rz(-1.9023583) q[2];
sx q[2];
rz(-2.6826774) q[2];
rz(-3.093873) q[3];
sx q[3];
rz(-2.9012091) q[3];
sx q[3];
rz(-0.11500959) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
