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
rz(2.202232) q[0];
sx q[0];
rz(-1.1190177) q[0];
sx q[0];
rz(-0.033666704) q[0];
rz(-1.4576003) q[1];
sx q[1];
rz(-1.8063318) q[1];
sx q[1];
rz(1.508498) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2062708) q[0];
sx q[0];
rz(-1.540872) q[0];
sx q[0];
rz(1.4564092) q[0];
rz(-pi) q[1];
rz(-1.8148957) q[2];
sx q[2];
rz(-0.30260049) q[2];
sx q[2];
rz(1.5989593) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.0284217) q[1];
sx q[1];
rz(-1.7152216) q[1];
sx q[1];
rz(2.444205) q[1];
rz(-pi) q[2];
rz(-2.5023828) q[3];
sx q[3];
rz(-2.2230801) q[3];
sx q[3];
rz(-2.2584884) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.17026751) q[2];
sx q[2];
rz(-2.2496932) q[2];
sx q[2];
rz(0.83667052) q[2];
rz(-0.80638805) q[3];
sx q[3];
rz(-2.216335) q[3];
sx q[3];
rz(-2.9562318) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.065501325) q[0];
sx q[0];
rz(-1.0668904) q[0];
sx q[0];
rz(-1.9245032) q[0];
rz(-0.20054664) q[1];
sx q[1];
rz(-0.79981083) q[1];
sx q[1];
rz(-2.5977871) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4763308) q[0];
sx q[0];
rz(-1.3419188) q[0];
sx q[0];
rz(-0.32461597) q[0];
rz(2.5466304) q[2];
sx q[2];
rz(-1.5595072) q[2];
sx q[2];
rz(-1.6807229) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.4098379) q[1];
sx q[1];
rz(-1.2827164) q[1];
sx q[1];
rz(1.9695915) q[1];
rz(-pi) q[2];
rz(-2.2460098) q[3];
sx q[3];
rz(-1.5981658) q[3];
sx q[3];
rz(0.41639183) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.66889888) q[2];
sx q[2];
rz(-0.74287477) q[2];
sx q[2];
rz(1.6061858) q[2];
rz(1.6425447) q[3];
sx q[3];
rz(-0.58321548) q[3];
sx q[3];
rz(-1.8892052) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1010308) q[0];
sx q[0];
rz(-2.5767548) q[0];
sx q[0];
rz(3.132013) q[0];
rz(1.0293695) q[1];
sx q[1];
rz(-1.5734438) q[1];
sx q[1];
rz(1.8052489) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.21214616) q[0];
sx q[0];
rz(-2.1597354) q[0];
sx q[0];
rz(0.94561823) q[0];
rz(-pi) q[1];
rz(-1.0687683) q[2];
sx q[2];
rz(-1.2867376) q[2];
sx q[2];
rz(1.9845923) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.6273986) q[1];
sx q[1];
rz(-1.3954867) q[1];
sx q[1];
rz(1.9216551) q[1];
x q[2];
rz(-1.1747178) q[3];
sx q[3];
rz(-1.4142591) q[3];
sx q[3];
rz(-0.96352947) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(3.0336527) q[2];
sx q[2];
rz(-0.5113217) q[2];
sx q[2];
rz(-1.3134202) q[2];
rz(0.039947346) q[3];
sx q[3];
rz(-2.2291456) q[3];
sx q[3];
rz(1.1279172) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0082598) q[0];
sx q[0];
rz(-1.8513716) q[0];
sx q[0];
rz(-0.3966575) q[0];
rz(0.89504129) q[1];
sx q[1];
rz(-1.4245234) q[1];
sx q[1];
rz(-0.76362124) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1930924) q[0];
sx q[0];
rz(-1.6512799) q[0];
sx q[0];
rz(1.8165239) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0020048) q[2];
sx q[2];
rz(-0.12171037) q[2];
sx q[2];
rz(-1.9892529) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.17814008) q[1];
sx q[1];
rz(-1.8082494) q[1];
sx q[1];
rz(3.1141539) q[1];
rz(-0.29080363) q[3];
sx q[3];
rz(-1.606308) q[3];
sx q[3];
rz(-2.5574656) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.7274196) q[2];
sx q[2];
rz(-1.5199993) q[2];
sx q[2];
rz(2.011389) q[2];
rz(-2.0027509) q[3];
sx q[3];
rz(-1.9500705) q[3];
sx q[3];
rz(-1.1677008) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.30911699) q[0];
sx q[0];
rz(-0.55823767) q[0];
sx q[0];
rz(-0.48156893) q[0];
rz(0.32749495) q[1];
sx q[1];
rz(-1.7227252) q[1];
sx q[1];
rz(0.96424261) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1314581) q[0];
sx q[0];
rz(-1.1842964) q[0];
sx q[0];
rz(2.9054427) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9498001) q[2];
sx q[2];
rz(-0.99770498) q[2];
sx q[2];
rz(-1.4258143) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.22589265) q[1];
sx q[1];
rz(-1.4488765) q[1];
sx q[1];
rz(1.3242766) q[1];
x q[2];
rz(1.7751415) q[3];
sx q[3];
rz(-1.1875523) q[3];
sx q[3];
rz(2.1658773) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.97011906) q[2];
sx q[2];
rz(-2.0233266) q[2];
sx q[2];
rz(1.4946651) q[2];
rz(-1.0129048) q[3];
sx q[3];
rz(-0.92618147) q[3];
sx q[3];
rz(-0.51500285) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
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
rz(-0.2038912) q[0];
sx q[0];
rz(-2.1370115) q[0];
sx q[0];
rz(-1.0585693) q[0];
rz(-0.92264289) q[1];
sx q[1];
rz(-2.1244996) q[1];
sx q[1];
rz(3.0386818) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9227717) q[0];
sx q[0];
rz(-1.1081244) q[0];
sx q[0];
rz(0.86313049) q[0];
rz(-2.5526664) q[2];
sx q[2];
rz(-1.803033) q[2];
sx q[2];
rz(-3.1265772) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.2051112) q[1];
sx q[1];
rz(-2.0983585) q[1];
sx q[1];
rz(2.9283306) q[1];
rz(-1.5636121) q[3];
sx q[3];
rz(-0.71643351) q[3];
sx q[3];
rz(-2.7897657) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.4395478) q[2];
sx q[2];
rz(-2.6921258) q[2];
sx q[2];
rz(-0.31583819) q[2];
rz(0.99509197) q[3];
sx q[3];
rz(-1.8683878) q[3];
sx q[3];
rz(-0.72527138) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.134326) q[0];
sx q[0];
rz(-2.8747929) q[0];
sx q[0];
rz(2.5079978) q[0];
rz(-2.7569356) q[1];
sx q[1];
rz(-1.1793914) q[1];
sx q[1];
rz(0.42919174) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4576526) q[0];
sx q[0];
rz(-0.54163092) q[0];
sx q[0];
rz(-1.0628048) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6056988) q[2];
sx q[2];
rz(-2.2560824) q[2];
sx q[2];
rz(-1.8165782) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.8365578) q[1];
sx q[1];
rz(-1.9227131) q[1];
sx q[1];
rz(2.5266527) q[1];
x q[2];
rz(2.6771992) q[3];
sx q[3];
rz(-2.0804318) q[3];
sx q[3];
rz(-1.9392175) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.089243285) q[2];
sx q[2];
rz(-2.3640335) q[2];
sx q[2];
rz(2.8746129) q[2];
rz(3.0698245) q[3];
sx q[3];
rz(-1.0177344) q[3];
sx q[3];
rz(1.4192386) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-0.6637904) q[0];
sx q[0];
rz(-0.1013805) q[0];
sx q[0];
rz(1.4005533) q[0];
rz(-1.9969223) q[1];
sx q[1];
rz(-1.0619699) q[1];
sx q[1];
rz(2.5544419) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5349197) q[0];
sx q[0];
rz(-0.79100709) q[0];
sx q[0];
rz(-2.530221) q[0];
rz(0.25803395) q[2];
sx q[2];
rz(-0.43083015) q[2];
sx q[2];
rz(-1.8547386) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.352601) q[1];
sx q[1];
rz(-1.7956211) q[1];
sx q[1];
rz(-2.4487752) q[1];
x q[2];
rz(0.31717698) q[3];
sx q[3];
rz(-0.97517386) q[3];
sx q[3];
rz(0.45444876) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.8553541) q[2];
sx q[2];
rz(-1.1796395) q[2];
sx q[2];
rz(0.41909763) q[2];
rz(1.6176443) q[3];
sx q[3];
rz(-0.91846275) q[3];
sx q[3];
rz(-1.5036478) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5986901) q[0];
sx q[0];
rz(-1.221523) q[0];
sx q[0];
rz(-0.30938095) q[0];
rz(-1.3503374) q[1];
sx q[1];
rz(-2.5534936) q[1];
sx q[1];
rz(2.5877171) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2387008) q[0];
sx q[0];
rz(-2.6052778) q[0];
sx q[0];
rz(-0.66030057) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.925681) q[2];
sx q[2];
rz(-0.40635525) q[2];
sx q[2];
rz(2.2333437) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.9123037) q[1];
sx q[1];
rz(-0.69432753) q[1];
sx q[1];
rz(-2.0635384) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.43954) q[3];
sx q[3];
rz(-1.3300196) q[3];
sx q[3];
rz(-2.05621) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.2498593) q[2];
sx q[2];
rz(-1.6102108) q[2];
sx q[2];
rz(-1.7943133) q[2];
rz(2.3740718) q[3];
sx q[3];
rz(-1.9556421) q[3];
sx q[3];
rz(-2.2931113) q[3];
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
rz(0.42221853) q[0];
sx q[0];
rz(-2.8415871) q[0];
sx q[0];
rz(-1.2147709) q[0];
rz(-1.0019373) q[1];
sx q[1];
rz(-1.5129713) q[1];
sx q[1];
rz(-0.92438662) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.11533212) q[0];
sx q[0];
rz(-0.9965653) q[0];
sx q[0];
rz(2.6747245) q[0];
x q[1];
rz(3.1141823) q[2];
sx q[2];
rz(-2.1152585) q[2];
sx q[2];
rz(1.7137063) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.1781285) q[1];
sx q[1];
rz(-2.917082) q[1];
sx q[1];
rz(-0.1046532) q[1];
rz(-pi) q[2];
rz(2.2697023) q[3];
sx q[3];
rz(-2.0465486) q[3];
sx q[3];
rz(-0.94199877) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.7944916) q[2];
sx q[2];
rz(-0.93239409) q[2];
sx q[2];
rz(2.7744228) q[2];
rz(1.1651039) q[3];
sx q[3];
rz(-0.96787435) q[3];
sx q[3];
rz(-0.86618209) q[3];
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
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1871056) q[0];
sx q[0];
rz(-0.73910537) q[0];
sx q[0];
rz(0.11866971) q[0];
rz(2.8754996) q[1];
sx q[1];
rz(-2.2292744) q[1];
sx q[1];
rz(1.6899065) q[1];
rz(2.8546147) q[2];
sx q[2];
rz(-1.3084044) q[2];
sx q[2];
rz(-0.85592782) q[2];
rz(2.1723971) q[3];
sx q[3];
rz(-1.7432913) q[3];
sx q[3];
rz(-0.023719214) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
