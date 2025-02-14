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
rz(0.27920029) q[0];
sx q[0];
rz(1.8912127) q[0];
sx q[0];
rz(9.425107) q[0];
rz(1.1286796) q[1];
sx q[1];
rz(5.1976701) q[1];
sx q[1];
rz(12.357098) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4930388) q[0];
sx q[0];
rz(-1.4161515) q[0];
sx q[0];
rz(1.8039319) q[0];
x q[1];
rz(-1.7369274) q[2];
sx q[2];
rz(-1.0175993) q[2];
sx q[2];
rz(0.90546135) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.5857201) q[1];
sx q[1];
rz(-1.2798847) q[1];
sx q[1];
rz(0.32679273) q[1];
rz(0.23926034) q[3];
sx q[3];
rz(-1.8094899) q[3];
sx q[3];
rz(-0.25588122) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.8802694) q[2];
sx q[2];
rz(-1.7797194) q[2];
sx q[2];
rz(-2.9051991) q[2];
rz(-0.38293019) q[3];
sx q[3];
rz(-2.0536486) q[3];
sx q[3];
rz(1.6898164) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7616538) q[0];
sx q[0];
rz(-2.8528657) q[0];
sx q[0];
rz(0.54288236) q[0];
rz(2.701243) q[1];
sx q[1];
rz(-2.0183225) q[1];
sx q[1];
rz(2.0139258) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2625439) q[0];
sx q[0];
rz(-2.194756) q[0];
sx q[0];
rz(-1.1406374) q[0];
x q[1];
rz(0.93483804) q[2];
sx q[2];
rz(-1.6976154) q[2];
sx q[2];
rz(0.93709125) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.7372514) q[1];
sx q[1];
rz(-1.7187498) q[1];
sx q[1];
rz(-2.3322991) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9502968) q[3];
sx q[3];
rz(-2.2436525) q[3];
sx q[3];
rz(1.6516453) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.27309624) q[2];
sx q[2];
rz(-0.79600969) q[2];
sx q[2];
rz(1.8491171) q[2];
rz(-2.3097307) q[3];
sx q[3];
rz(-1.6658733) q[3];
sx q[3];
rz(-2.7003435) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6865987) q[0];
sx q[0];
rz(-0.43231493) q[0];
sx q[0];
rz(1.5469714) q[0];
rz(3.1349685) q[1];
sx q[1];
rz(-0.5358271) q[1];
sx q[1];
rz(-2.5033902) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2062195) q[0];
sx q[0];
rz(-1.5347592) q[0];
sx q[0];
rz(-1.0956148) q[0];
rz(-pi) q[1];
rz(-2.7225627) q[2];
sx q[2];
rz(-2.3263676) q[2];
sx q[2];
rz(2.4885675) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.51391006) q[1];
sx q[1];
rz(-1.2177151) q[1];
sx q[1];
rz(2.7536235) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.21535899) q[3];
sx q[3];
rz(-2.1073859) q[3];
sx q[3];
rz(-0.99765618) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.39614761) q[2];
sx q[2];
rz(-3.1148995) q[2];
sx q[2];
rz(-0.8566345) q[2];
rz(0.73532295) q[3];
sx q[3];
rz(-1.7263128) q[3];
sx q[3];
rz(1.0480688) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.51744866) q[0];
sx q[0];
rz(-0.051001661) q[0];
sx q[0];
rz(0.93358246) q[0];
rz(-0.75992641) q[1];
sx q[1];
rz(-2.3257207) q[1];
sx q[1];
rz(1.5911721) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4226993) q[0];
sx q[0];
rz(-0.29445364) q[0];
sx q[0];
rz(0.30979777) q[0];
x q[1];
rz(0.43689219) q[2];
sx q[2];
rz(-1.5073188) q[2];
sx q[2];
rz(1.2530099) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.62867579) q[1];
sx q[1];
rz(-2.4077291) q[1];
sx q[1];
rz(0.34725125) q[1];
rz(-2.5752034) q[3];
sx q[3];
rz(-0.29367125) q[3];
sx q[3];
rz(1.3704513) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.7237225) q[2];
sx q[2];
rz(-1.3891862) q[2];
sx q[2];
rz(2.1088481) q[2];
rz(0.51747733) q[3];
sx q[3];
rz(-1.7101219) q[3];
sx q[3];
rz(-1.3051858) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1037769) q[0];
sx q[0];
rz(-3.0396437) q[0];
sx q[0];
rz(1.284449) q[0];
rz(-2.2880554) q[1];
sx q[1];
rz(-1.6910911) q[1];
sx q[1];
rz(-2.6654713) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42890756) q[0];
sx q[0];
rz(-2.0715014) q[0];
sx q[0];
rz(-2.7963573) q[0];
x q[1];
rz(-1.2849152) q[2];
sx q[2];
rz(-1.21884) q[2];
sx q[2];
rz(0.98796036) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.0796226) q[1];
sx q[1];
rz(-1.274612) q[1];
sx q[1];
rz(-2.4825063) q[1];
rz(-1.4912603) q[3];
sx q[3];
rz(-2.400268) q[3];
sx q[3];
rz(1.1363514) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.7242754) q[2];
sx q[2];
rz(-2.3099835) q[2];
sx q[2];
rz(-1.6012491) q[2];
rz(-0.075751461) q[3];
sx q[3];
rz(-1.7509533) q[3];
sx q[3];
rz(0.72995228) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1002355) q[0];
sx q[0];
rz(-2.5900109) q[0];
sx q[0];
rz(-1.3496189) q[0];
rz(1.8655818) q[1];
sx q[1];
rz(-0.7834692) q[1];
sx q[1];
rz(1.0412201) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91105499) q[0];
sx q[0];
rz(-1.7887702) q[0];
sx q[0];
rz(-2.387052) q[0];
rz(1.6068993) q[2];
sx q[2];
rz(-2.3182333) q[2];
sx q[2];
rz(-0.90561101) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.7752616) q[1];
sx q[1];
rz(-1.1995254) q[1];
sx q[1];
rz(-1.8410147) q[1];
x q[2];
rz(-0.75821782) q[3];
sx q[3];
rz(-1.3617523) q[3];
sx q[3];
rz(0.93702873) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.87018806) q[2];
sx q[2];
rz(-1.1689593) q[2];
sx q[2];
rz(2.1427587) q[2];
rz(2.6073604) q[3];
sx q[3];
rz(-1.9637008) q[3];
sx q[3];
rz(-1.6629201) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6967195) q[0];
sx q[0];
rz(-2.1598926) q[0];
sx q[0];
rz(1.5848507) q[0];
rz(-2.059767) q[1];
sx q[1];
rz(-1.4993246) q[1];
sx q[1];
rz(-2.7781442) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2242607) q[0];
sx q[0];
rz(-1.5646569) q[0];
sx q[0];
rz(0.21256769) q[0];
rz(-pi) q[1];
rz(-2.8055023) q[2];
sx q[2];
rz(-0.31340965) q[2];
sx q[2];
rz(1.4614507) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.195943) q[1];
sx q[1];
rz(-2.2684625) q[1];
sx q[1];
rz(2.5673668) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0601252) q[3];
sx q[3];
rz(-2.8301468) q[3];
sx q[3];
rz(-1.2490727) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.58793679) q[2];
sx q[2];
rz(-0.16599545) q[2];
sx q[2];
rz(2.0641649) q[2];
rz(1.5466127) q[3];
sx q[3];
rz(-1.8190106) q[3];
sx q[3];
rz(-0.67302978) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8107635) q[0];
sx q[0];
rz(-1.6442693) q[0];
sx q[0];
rz(-0.28054917) q[0];
rz(2.5827017) q[1];
sx q[1];
rz(-2.2160896) q[1];
sx q[1];
rz(1.2116609) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9292363) q[0];
sx q[0];
rz(-2.2693737) q[0];
sx q[0];
rz(0.86232604) q[0];
x q[1];
rz(-0.49254353) q[2];
sx q[2];
rz(-1.768868) q[2];
sx q[2];
rz(-0.98392526) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.16226823) q[1];
sx q[1];
rz(-0.78020243) q[1];
sx q[1];
rz(3.1092572) q[1];
rz(-1.3727851) q[3];
sx q[3];
rz(-1.565838) q[3];
sx q[3];
rz(-1.229242) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.1739379) q[2];
sx q[2];
rz(-1.329198) q[2];
sx q[2];
rz(2.3363414) q[2];
rz(1.1236745) q[3];
sx q[3];
rz(-1.1403964) q[3];
sx q[3];
rz(2.8203188) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5241908) q[0];
sx q[0];
rz(-2.5179355) q[0];
sx q[0];
rz(-2.6362841) q[0];
rz(-0.71120039) q[1];
sx q[1];
rz(-2.2732747) q[1];
sx q[1];
rz(2.0288846) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.18234564) q[0];
sx q[0];
rz(-1.2768978) q[0];
sx q[0];
rz(2.6698851) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6518902) q[2];
sx q[2];
rz(-1.4633274) q[2];
sx q[2];
rz(-1.7484959) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.6869097) q[1];
sx q[1];
rz(-2.0078604) q[1];
sx q[1];
rz(-0.13961643) q[1];
rz(-pi) q[2];
x q[2];
rz(0.030783665) q[3];
sx q[3];
rz(-1.0153706) q[3];
sx q[3];
rz(-2.8419213) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.6342371) q[2];
sx q[2];
rz(-1.2398182) q[2];
sx q[2];
rz(-2.9061387) q[2];
rz(-2.7638451) q[3];
sx q[3];
rz(-2.7538959) q[3];
sx q[3];
rz(0.050067576) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3560051) q[0];
sx q[0];
rz(-1.2879141) q[0];
sx q[0];
rz(3.0872784) q[0];
rz(2.3771225) q[1];
sx q[1];
rz(-0.46934325) q[1];
sx q[1];
rz(2.4321709) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5886836) q[0];
sx q[0];
rz(-1.3257265) q[0];
sx q[0];
rz(0.10677307) q[0];
x q[1];
rz(-1.7058175) q[2];
sx q[2];
rz(-1.3121989) q[2];
sx q[2];
rz(-1.0279581) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.0065942) q[1];
sx q[1];
rz(-1.6853764) q[1];
sx q[1];
rz(-2.5496818) q[1];
rz(-pi) q[2];
x q[2];
rz(0.34281369) q[3];
sx q[3];
rz(-1.4488932) q[3];
sx q[3];
rz(1.3212622) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.8556122) q[2];
sx q[2];
rz(-1.4748272) q[2];
sx q[2];
rz(3.0484071) q[2];
rz(-2.5423673) q[3];
sx q[3];
rz(-0.56120509) q[3];
sx q[3];
rz(-0.78527251) q[3];
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
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7080606) q[0];
sx q[0];
rz(-0.96994937) q[0];
sx q[0];
rz(2.1920828) q[0];
rz(-0.064619725) q[1];
sx q[1];
rz(-3.0977991) q[1];
sx q[1];
rz(-2.480712) q[1];
rz(-1.4501889) q[2];
sx q[2];
rz(-0.55477886) q[2];
sx q[2];
rz(3.0045493) q[2];
rz(0.6653847) q[3];
sx q[3];
rz(-1.7177842) q[3];
sx q[3];
rz(-0.73242889) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
