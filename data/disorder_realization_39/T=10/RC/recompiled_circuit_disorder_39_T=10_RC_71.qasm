OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.0467779) q[0];
sx q[0];
rz(-1.0682286) q[0];
sx q[0];
rz(-0.46407035) q[0];
rz(-1.1820473) q[1];
sx q[1];
rz(-3.074475) q[1];
sx q[1];
rz(-1.2844515) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.099079236) q[0];
sx q[0];
rz(-0.49603841) q[0];
sx q[0];
rz(-2.2201559) q[0];
rz(0.061878248) q[2];
sx q[2];
rz(-0.49052325) q[2];
sx q[2];
rz(2.2762736) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.1844993) q[1];
sx q[1];
rz(-0.59809369) q[1];
sx q[1];
rz(2.4615272) q[1];
rz(-pi) q[2];
x q[2];
rz(0.10847096) q[3];
sx q[3];
rz(-0.15158187) q[3];
sx q[3];
rz(-0.82419318) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.39711943) q[2];
sx q[2];
rz(-2.2019272) q[2];
sx q[2];
rz(-2.822067) q[2];
rz(-0.5784936) q[3];
sx q[3];
rz(-2.6632023) q[3];
sx q[3];
rz(2.4676676) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4085061) q[0];
sx q[0];
rz(-1.5717614) q[0];
sx q[0];
rz(2.615036) q[0];
rz(-2.5780442) q[1];
sx q[1];
rz(-1.5310409) q[1];
sx q[1];
rz(-2.3449576) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5986886) q[0];
sx q[0];
rz(-1.2447378) q[0];
sx q[0];
rz(0.90755264) q[0];
rz(-pi) q[1];
x q[1];
rz(0.57467069) q[2];
sx q[2];
rz(-1.7386912) q[2];
sx q[2];
rz(1.6811973) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.4371725) q[1];
sx q[1];
rz(-2.17802) q[1];
sx q[1];
rz(1.5864658) q[1];
rz(0.96062406) q[3];
sx q[3];
rz(-1.1517797) q[3];
sx q[3];
rz(-2.744439) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.9600296) q[2];
sx q[2];
rz(-1.7627565) q[2];
sx q[2];
rz(-2.3550418) q[2];
rz(2.6484047) q[3];
sx q[3];
rz(-1.1816053) q[3];
sx q[3];
rz(-2.6942159) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86984533) q[0];
sx q[0];
rz(-0.87690502) q[0];
sx q[0];
rz(1.7012117) q[0];
rz(-0.72021833) q[1];
sx q[1];
rz(-2.5455988) q[1];
sx q[1];
rz(0.70297855) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5539348) q[0];
sx q[0];
rz(-2.8371187) q[0];
sx q[0];
rz(-1.1712043) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5419652) q[2];
sx q[2];
rz(-0.58279524) q[2];
sx q[2];
rz(-1.3404913) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.86537251) q[1];
sx q[1];
rz(-0.11905383) q[1];
sx q[1];
rz(2.7369376) q[1];
x q[2];
rz(-1.8099144) q[3];
sx q[3];
rz(-2.5044887) q[3];
sx q[3];
rz(1.0106196) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.8748223) q[2];
sx q[2];
rz(-1.6569933) q[2];
sx q[2];
rz(2.2300569) q[2];
rz(-2.2276145) q[3];
sx q[3];
rz(-1.4303047) q[3];
sx q[3];
rz(-2.3142557) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5754958) q[0];
sx q[0];
rz(-0.63593447) q[0];
sx q[0];
rz(-2.3655868) q[0];
rz(-1.874118) q[1];
sx q[1];
rz(-1.1088561) q[1];
sx q[1];
rz(-0.56328303) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.253809) q[0];
sx q[0];
rz(-0.94622181) q[0];
sx q[0];
rz(-0.33872351) q[0];
x q[1];
rz(2.7411555) q[2];
sx q[2];
rz(-2.4106328) q[2];
sx q[2];
rz(-1.2618582) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.5717585) q[1];
sx q[1];
rz(-0.27184871) q[1];
sx q[1];
rz(0.024072577) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8068845) q[3];
sx q[3];
rz(-0.74581205) q[3];
sx q[3];
rz(-2.0832182) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.6528066) q[2];
sx q[2];
rz(-2.2061009) q[2];
sx q[2];
rz(0.164786) q[2];
rz(-0.22848836) q[3];
sx q[3];
rz(-2.861664) q[3];
sx q[3];
rz(-0.94648186) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.27424681) q[0];
sx q[0];
rz(-2.2773401) q[0];
sx q[0];
rz(0.85246032) q[0];
rz(2.7903941) q[1];
sx q[1];
rz(-2.5636667) q[1];
sx q[1];
rz(1.16211) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6974555) q[0];
sx q[0];
rz(-1.1778957) q[0];
sx q[0];
rz(-0.59209728) q[0];
x q[1];
rz(0.94050546) q[2];
sx q[2];
rz(-0.61912196) q[2];
sx q[2];
rz(1.8096015) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.3106766) q[1];
sx q[1];
rz(-1.1569996) q[1];
sx q[1];
rz(-2.3526741) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4200053) q[3];
sx q[3];
rz(-1.5820832) q[3];
sx q[3];
rz(1.23502) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.6953485) q[2];
sx q[2];
rz(-2.7480795) q[2];
sx q[2];
rz(-1.2472786) q[2];
rz(0.013109664) q[3];
sx q[3];
rz(-1.40991) q[3];
sx q[3];
rz(-2.9124027) q[3];
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
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2728249) q[0];
sx q[0];
rz(-1.2478963) q[0];
sx q[0];
rz(2.390958) q[0];
rz(1.1095095) q[1];
sx q[1];
rz(-2.7710319) q[1];
sx q[1];
rz(-1.3060588) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7795777) q[0];
sx q[0];
rz(-2.0834196) q[0];
sx q[0];
rz(2.6984452) q[0];
rz(1.7734217) q[2];
sx q[2];
rz(-1.1154004) q[2];
sx q[2];
rz(2.1855598) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.2893147) q[1];
sx q[1];
rz(-1.5631952) q[1];
sx q[1];
rz(-0.72035933) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7197051) q[3];
sx q[3];
rz(-0.9871261) q[3];
sx q[3];
rz(0.35051171) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.7541472) q[2];
sx q[2];
rz(-1.0070609) q[2];
sx q[2];
rz(-0.60069096) q[2];
rz(-1.0026275) q[3];
sx q[3];
rz(-2.6032344) q[3];
sx q[3];
rz(-0.35282648) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7464741) q[0];
sx q[0];
rz(-1.4557319) q[0];
sx q[0];
rz(2.0948998) q[0];
rz(1.5294317) q[1];
sx q[1];
rz(-1.4271095) q[1];
sx q[1];
rz(-0.41710645) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0102651) q[0];
sx q[0];
rz(-1.7908887) q[0];
sx q[0];
rz(-3.1165645) q[0];
x q[1];
rz(-0.76191683) q[2];
sx q[2];
rz(-2.7284405) q[2];
sx q[2];
rz(-2.6921536) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.68361359) q[1];
sx q[1];
rz(-0.90404592) q[1];
sx q[1];
rz(-0.97138202) q[1];
rz(-pi) q[2];
x q[2];
rz(0.38325558) q[3];
sx q[3];
rz(-2.4517422) q[3];
sx q[3];
rz(-0.81724973) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.0059011857) q[2];
sx q[2];
rz(-2.6689172) q[2];
sx q[2];
rz(1.7741514) q[2];
rz(-0.55316365) q[3];
sx q[3];
rz(-1.7382712) q[3];
sx q[3];
rz(2.6161391) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.816514) q[0];
sx q[0];
rz(-1.8186318) q[0];
sx q[0];
rz(1.9518071) q[0];
rz(-1.7143543) q[1];
sx q[1];
rz(-2.863527) q[1];
sx q[1];
rz(-3.0292125) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.053514078) q[0];
sx q[0];
rz(-1.7435762) q[0];
sx q[0];
rz(2.0515576) q[0];
x q[1];
rz(1.4833647) q[2];
sx q[2];
rz(-1.2382675) q[2];
sx q[2];
rz(0.80966262) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.6205794) q[1];
sx q[1];
rz(-0.81110209) q[1];
sx q[1];
rz(2.051342) q[1];
rz(-pi) q[2];
rz(-1.8295248) q[3];
sx q[3];
rz(-1.4530164) q[3];
sx q[3];
rz(0.044101276) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.0671976) q[2];
sx q[2];
rz(-1.2290359) q[2];
sx q[2];
rz(1.7929662) q[2];
rz(-1.2049234) q[3];
sx q[3];
rz(-1.4068853) q[3];
sx q[3];
rz(-2.9437734) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39712054) q[0];
sx q[0];
rz(-1.67698) q[0];
sx q[0];
rz(0.034974139) q[0];
rz(0.84683013) q[1];
sx q[1];
rz(-1.7085107) q[1];
sx q[1];
rz(0.91167489) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.864894) q[0];
sx q[0];
rz(-1.8407341) q[0];
sx q[0];
rz(0.2322659) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9334925) q[2];
sx q[2];
rz(-0.55317438) q[2];
sx q[2];
rz(2.2262239) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.8988077) q[1];
sx q[1];
rz(-0.50711942) q[1];
sx q[1];
rz(0.42906638) q[1];
rz(-pi) q[2];
rz(2.4089912) q[3];
sx q[3];
rz(-2.4676975) q[3];
sx q[3];
rz(0.037308824) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.4336865) q[2];
sx q[2];
rz(-2.6277379) q[2];
sx q[2];
rz(-0.24924499) q[2];
rz(0.76672673) q[3];
sx q[3];
rz(-1.3336351) q[3];
sx q[3];
rz(-0.39961091) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7837759) q[0];
sx q[0];
rz(-2.0932842) q[0];
sx q[0];
rz(0.37208474) q[0];
rz(2.5601939) q[1];
sx q[1];
rz(-2.364368) q[1];
sx q[1];
rz(-1.7262329) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.85522643) q[0];
sx q[0];
rz(-0.33272538) q[0];
sx q[0];
rz(-2.4871728) q[0];
rz(-1.0572817) q[2];
sx q[2];
rz(-1.0143177) q[2];
sx q[2];
rz(-1.8884115) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.6485939) q[1];
sx q[1];
rz(-0.10222888) q[1];
sx q[1];
rz(2.3962254) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7578027) q[3];
sx q[3];
rz(-1.0597611) q[3];
sx q[3];
rz(-0.25110652) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.32594484) q[2];
sx q[2];
rz(-0.14020136) q[2];
sx q[2];
rz(2.9369205) q[2];
rz(1.4137911) q[3];
sx q[3];
rz(-1.9982343) q[3];
sx q[3];
rz(-1.0958825) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.50080147) q[0];
sx q[0];
rz(-1.6727819) q[0];
sx q[0];
rz(0.93760136) q[0];
rz(1.5564556) q[1];
sx q[1];
rz(-1.382348) q[1];
sx q[1];
rz(1.4437645) q[1];
rz(1.3265058) q[2];
sx q[2];
rz(-2.6752224) q[2];
sx q[2];
rz(0.35717076) q[2];
rz(1.4428044) q[3];
sx q[3];
rz(-2.5656869) q[3];
sx q[3];
rz(2.2681469) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
