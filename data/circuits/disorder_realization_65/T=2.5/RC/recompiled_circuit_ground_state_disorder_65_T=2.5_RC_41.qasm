OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.559691) q[0];
sx q[0];
rz(-0.087495916) q[0];
sx q[0];
rz(0.09905941) q[0];
rz(2.8506408) q[1];
sx q[1];
rz(-1.7506316) q[1];
sx q[1];
rz(1.6265534) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0578931) q[0];
sx q[0];
rz(-2.0285995) q[0];
sx q[0];
rz(-2.8843736) q[0];
rz(-pi) q[1];
rz(2.0798105) q[2];
sx q[2];
rz(-1.1689593) q[2];
sx q[2];
rz(-0.33515829) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.6219306) q[1];
sx q[1];
rz(-0.77072137) q[1];
sx q[1];
rz(1.8008485) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8232075) q[3];
sx q[3];
rz(-1.2947467) q[3];
sx q[3];
rz(1.6161671) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.060595162) q[2];
sx q[2];
rz(-3.1298349) q[2];
sx q[2];
rz(-0.11059977) q[2];
rz(-3.1317173) q[3];
sx q[3];
rz(-3.1275833) q[3];
sx q[3];
rz(-0.88147718) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9843543) q[0];
sx q[0];
rz(-0.088033661) q[0];
sx q[0];
rz(-1.3798168) q[0];
rz(0.012160483) q[1];
sx q[1];
rz(-2.1575243) q[1];
sx q[1];
rz(1.5252569) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1943239) q[0];
sx q[0];
rz(-2.8229694) q[0];
sx q[0];
rz(2.9605894) q[0];
rz(-pi) q[1];
rz(-1.3697769) q[2];
sx q[2];
rz(-1.9819385) q[2];
sx q[2];
rz(-2.0697921) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.7213707) q[1];
sx q[1];
rz(-0.18850148) q[1];
sx q[1];
rz(2.7367715) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4472876) q[3];
sx q[3];
rz(-1.12772) q[3];
sx q[3];
rz(-0.62813891) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.3297367) q[2];
sx q[2];
rz(-0.024710329) q[2];
sx q[2];
rz(-3.109566) q[2];
rz(0.61102593) q[3];
sx q[3];
rz(-1.150584) q[3];
sx q[3];
rz(-1.376763) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
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
rz(1.5830314) q[0];
sx q[0];
rz(-0.91201454) q[0];
sx q[0];
rz(-1.4645905) q[0];
rz(-1.6327845) q[1];
sx q[1];
rz(-1.8784411) q[1];
sx q[1];
rz(0.060922932) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2476544) q[0];
sx q[0];
rz(-1.1919855) q[0];
sx q[0];
rz(-1.768186) q[0];
rz(-pi) q[1];
rz(-2.3192984) q[2];
sx q[2];
rz(-2.9884183) q[2];
sx q[2];
rz(0.3016475) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.12184252) q[1];
sx q[1];
rz(-1.3977244) q[1];
sx q[1];
rz(2.9815841) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.39037986) q[3];
sx q[3];
rz(-1.8515311) q[3];
sx q[3];
rz(0.27306199) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.1724725) q[2];
sx q[2];
rz(-3.1388333) q[2];
sx q[2];
rz(-0.74392444) q[2];
rz(-0.38532358) q[3];
sx q[3];
rz(-0.0016366882) q[3];
sx q[3];
rz(-2.2374432) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0124403) q[0];
sx q[0];
rz(-3.0991982) q[0];
sx q[0];
rz(-0.84268919) q[0];
rz(0.95376247) q[1];
sx q[1];
rz(-1.6396435) q[1];
sx q[1];
rz(-1.5488254) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4207537) q[0];
sx q[0];
rz(-0.11891236) q[0];
sx q[0];
rz(-0.17158385) q[0];
rz(2.983613) q[2];
sx q[2];
rz(-1.9648187) q[2];
sx q[2];
rz(-1.3713149) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.4348244) q[1];
sx q[1];
rz(-0.93509669) q[1];
sx q[1];
rz(0.058556538) q[1];
rz(0.78068818) q[3];
sx q[3];
rz(-1.4997291) q[3];
sx q[3];
rz(2.4960757) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(3.039975) q[2];
sx q[2];
rz(-3.1246694) q[2];
sx q[2];
rz(0.42538154) q[2];
rz(1.7915122) q[3];
sx q[3];
rz(-1.5964419) q[3];
sx q[3];
rz(-0.30459705) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72163433) q[0];
sx q[0];
rz(-0.0036792734) q[0];
sx q[0];
rz(-0.71596181) q[0];
rz(1.6250027) q[1];
sx q[1];
rz(-0.45418987) q[1];
sx q[1];
rz(-2.9862278) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6778188) q[0];
sx q[0];
rz(-1.4646104) q[0];
sx q[0];
rz(3.1255577) q[0];
rz(-pi) q[1];
x q[1];
rz(0.22414872) q[2];
sx q[2];
rz(-2.756065) q[2];
sx q[2];
rz(1.7412501) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.6581056) q[1];
sx q[1];
rz(-2.8466821) q[1];
sx q[1];
rz(-0.43842478) q[1];
rz(-2.0955438) q[3];
sx q[3];
rz(-2.397846) q[3];
sx q[3];
rz(2.634457) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.6905489) q[2];
sx q[2];
rz(-2.0205708) q[2];
sx q[2];
rz(1.5345908) q[2];
rz(2.0524041) q[3];
sx q[3];
rz(-0.023622731) q[3];
sx q[3];
rz(1.3634118) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
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
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4264193) q[0];
sx q[0];
rz(-0.024113163) q[0];
sx q[0];
rz(-2.5809848) q[0];
rz(-2.2757065) q[1];
sx q[1];
rz(-0.0090323369) q[1];
sx q[1];
rz(-2.3057356) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7880121) q[0];
sx q[0];
rz(-0.29639527) q[0];
sx q[0];
rz(0.061876492) q[0];
rz(-pi) q[1];
rz(-1.4908815) q[2];
sx q[2];
rz(-0.40849388) q[2];
sx q[2];
rz(-1.6871883) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.33619704) q[1];
sx q[1];
rz(-1.8289036) q[1];
sx q[1];
rz(-1.6805524) q[1];
rz(-pi) q[2];
x q[2];
rz(0.31536169) q[3];
sx q[3];
rz(-1.3360607) q[3];
sx q[3];
rz(-2.0965337) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.5777638) q[2];
sx q[2];
rz(-1.8356297) q[2];
sx q[2];
rz(-1.5888265) q[2];
rz(2.0896437) q[3];
sx q[3];
rz(-0.51783872) q[3];
sx q[3];
rz(-0.53798211) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.8727259) q[0];
sx q[0];
rz(-0.59393847) q[0];
sx q[0];
rz(-1.8178222) q[0];
rz(-2.2920604) q[1];
sx q[1];
rz(-0.001566611) q[1];
sx q[1];
rz(0.82226396) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2182151) q[0];
sx q[0];
rz(-1.4750622) q[0];
sx q[0];
rz(1.7507919) q[0];
rz(0.33820196) q[2];
sx q[2];
rz(-1.6227007) q[2];
sx q[2];
rz(-3.1171892) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.9692053) q[1];
sx q[1];
rz(-1.6403461) q[1];
sx q[1];
rz(-3.1350673) q[1];
x q[2];
rz(2.8760738) q[3];
sx q[3];
rz(-1.3808704) q[3];
sx q[3];
rz(-1.3914212) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.5795341) q[2];
sx q[2];
rz(-1.2349393) q[2];
sx q[2];
rz(-2.298992) q[2];
rz(-0.54251999) q[3];
sx q[3];
rz(-0.018535651) q[3];
sx q[3];
rz(-2.5491469) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(0.056996718) q[0];
sx q[0];
rz(-1.5942986) q[0];
sx q[0];
rz(-2.090276) q[0];
rz(1.7683138) q[1];
sx q[1];
rz(-3.1128502) q[1];
sx q[1];
rz(-1.9017259) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.97530203) q[0];
sx q[0];
rz(-2.8956415) q[0];
sx q[0];
rz(0.16157948) q[0];
rz(-pi) q[1];
rz(-3.1048959) q[2];
sx q[2];
rz(-1.1783491) q[2];
sx q[2];
rz(1.6709171) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.6426601) q[1];
sx q[1];
rz(-2.1524384) q[1];
sx q[1];
rz(-2.8967392) q[1];
rz(-pi) q[2];
rz(-1.3040335) q[3];
sx q[3];
rz(-1.752231) q[3];
sx q[3];
rz(-0.72600466) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.5956868) q[2];
sx q[2];
rz(-3.1133856) q[2];
sx q[2];
rz(-0.14740454) q[2];
rz(-1.5386511) q[3];
sx q[3];
rz(-1.7037062) q[3];
sx q[3];
rz(0.182972) q[3];
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
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1786757) q[0];
sx q[0];
rz(-2.8067532) q[0];
sx q[0];
rz(-1.8033002) q[0];
rz(-1.6637404) q[1];
sx q[1];
rz(-2.0070952) q[1];
sx q[1];
rz(0.17325625) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6918744) q[0];
sx q[0];
rz(-1.0074179) q[0];
sx q[0];
rz(1.059235) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0031172) q[2];
sx q[2];
rz(-0.67174339) q[2];
sx q[2];
rz(-2.9152108) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.4819255) q[1];
sx q[1];
rz(-2.8489981) q[1];
sx q[1];
rz(0.56168075) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.162649) q[3];
sx q[3];
rz(-3.1269073) q[3];
sx q[3];
rz(-0.56416184) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.0080537) q[2];
sx q[2];
rz(-3.1296788) q[2];
sx q[2];
rz(2.1214205) q[2];
rz(3.0607439) q[3];
sx q[3];
rz(-0.0011708502) q[3];
sx q[3];
rz(1.1297869) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5197649) q[0];
sx q[0];
rz(-0.28525678) q[0];
sx q[0];
rz(1.6602302) q[0];
rz(0.18352428) q[1];
sx q[1];
rz(-0.32569519) q[1];
sx q[1];
rz(0.087800177) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7526898) q[0];
sx q[0];
rz(-1.6713033) q[0];
sx q[0];
rz(-2.697562) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4884639) q[2];
sx q[2];
rz(-1.8762686) q[2];
sx q[2];
rz(-0.20037547) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.99276772) q[1];
sx q[1];
rz(-1.449905) q[1];
sx q[1];
rz(-2.8962027) q[1];
rz(-pi) q[2];
rz(-0.3751288) q[3];
sx q[3];
rz(-0.70444706) q[3];
sx q[3];
rz(-0.40788996) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.8537019) q[2];
sx q[2];
rz(-0.012306865) q[2];
sx q[2];
rz(-2.2575209) q[2];
rz(-0.70841241) q[3];
sx q[3];
rz(-3.1413779) q[3];
sx q[3];
rz(0.29402548) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.52083279) q[0];
sx q[0];
rz(-1.5721058) q[0];
sx q[0];
rz(1.5097621) q[0];
rz(-0.028342551) q[1];
sx q[1];
rz(-2.5181073) q[1];
sx q[1];
rz(0.070652031) q[1];
rz(-1.2030884) q[2];
sx q[2];
rz(-1.1496419) q[2];
sx q[2];
rz(-2.9436191) q[2];
rz(-1.9020756) q[3];
sx q[3];
rz(-0.86994949) q[3];
sx q[3];
rz(0.89350064) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
