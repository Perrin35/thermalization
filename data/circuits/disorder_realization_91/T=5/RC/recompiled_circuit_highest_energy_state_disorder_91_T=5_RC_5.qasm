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
rz(0.25536594) q[0];
sx q[0];
rz(-0.74170434) q[0];
sx q[0];
rz(1.1852888) q[0];
rz(1.9396632) q[1];
sx q[1];
rz(-0.97768205) q[1];
sx q[1];
rz(0.77562195) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5016504) q[0];
sx q[0];
rz(-1.3222561) q[0];
sx q[0];
rz(2.2959501) q[0];
rz(-pi) q[1];
rz(-0.94812265) q[2];
sx q[2];
rz(-0.97062696) q[2];
sx q[2];
rz(-2.4400939) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(3.0518912) q[1];
sx q[1];
rz(-0.49124733) q[1];
sx q[1];
rz(-2.6253176) q[1];
x q[2];
rz(-3.1235141) q[3];
sx q[3];
rz(-2.5466777) q[3];
sx q[3];
rz(-0.48479652) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.1620907) q[2];
sx q[2];
rz(-2.4544331) q[2];
sx q[2];
rz(-3.0882351) q[2];
rz(-0.067367628) q[3];
sx q[3];
rz(-2.8991883) q[3];
sx q[3];
rz(-1.7072201) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1783717) q[0];
sx q[0];
rz(-1.0472949) q[0];
sx q[0];
rz(-0.50985992) q[0];
rz(2.0193822) q[1];
sx q[1];
rz(-1.1374067) q[1];
sx q[1];
rz(3.1013464) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3771693) q[0];
sx q[0];
rz(-2.3413045) q[0];
sx q[0];
rz(-1.1616526) q[0];
x q[1];
rz(0.031304403) q[2];
sx q[2];
rz(-2.207216) q[2];
sx q[2];
rz(0.49963504) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.486301) q[1];
sx q[1];
rz(-2.5319063) q[1];
sx q[1];
rz(-0.31368452) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5882309) q[3];
sx q[3];
rz(-2.2477345) q[3];
sx q[3];
rz(1.5172322) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.76790729) q[2];
sx q[2];
rz(-1.7932971) q[2];
sx q[2];
rz(0.56758869) q[2];
rz(-0.23809412) q[3];
sx q[3];
rz(-1.6425491) q[3];
sx q[3];
rz(-0.25180086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0955065) q[0];
sx q[0];
rz(-2.6987785) q[0];
sx q[0];
rz(0.96389043) q[0];
rz(-1.6148199) q[1];
sx q[1];
rz(-1.2232199) q[1];
sx q[1];
rz(0.33178202) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.362865) q[0];
sx q[0];
rz(-2.300867) q[0];
sx q[0];
rz(1.8018468) q[0];
rz(-pi) q[1];
rz(-1.2305698) q[2];
sx q[2];
rz(-1.548049) q[2];
sx q[2];
rz(2.0024254) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.7501237) q[1];
sx q[1];
rz(-1.5771598) q[1];
sx q[1];
rz(-1.6666189) q[1];
rz(2.0121215) q[3];
sx q[3];
rz(-1.6989162) q[3];
sx q[3];
rz(0.81525034) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.31384599) q[2];
sx q[2];
rz(-1.4705986) q[2];
sx q[2];
rz(-3.1249135) q[2];
rz(-1.553675) q[3];
sx q[3];
rz(-1.1976306) q[3];
sx q[3];
rz(-0.29271075) q[3];
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
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68506771) q[0];
sx q[0];
rz(-2.2790788) q[0];
sx q[0];
rz(2.9271794) q[0];
rz(2.9247177) q[1];
sx q[1];
rz(-0.82019359) q[1];
sx q[1];
rz(-2.0490501) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8181145) q[0];
sx q[0];
rz(-1.406455) q[0];
sx q[0];
rz(1.8002609) q[0];
rz(-3.0472146) q[2];
sx q[2];
rz(-2.1571231) q[2];
sx q[2];
rz(0.5273158) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.34784029) q[1];
sx q[1];
rz(-1.5303601) q[1];
sx q[1];
rz(-2.3024125) q[1];
x q[2];
rz(2.2762466) q[3];
sx q[3];
rz(-2.1095037) q[3];
sx q[3];
rz(2.5821843) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.0715926) q[2];
sx q[2];
rz(-2.6061974) q[2];
sx q[2];
rz(0.52616057) q[2];
rz(-1.5782662) q[3];
sx q[3];
rz(-1.9901146) q[3];
sx q[3];
rz(0.17939849) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9569595) q[0];
sx q[0];
rz(-1.7687836) q[0];
sx q[0];
rz(-2.2233295) q[0];
rz(-2.1271162) q[1];
sx q[1];
rz(-2.4160073) q[1];
sx q[1];
rz(-2.0065506) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.43111463) q[0];
sx q[0];
rz(-1.6621792) q[0];
sx q[0];
rz(0.26510948) q[0];
x q[1];
rz(1.1667541) q[2];
sx q[2];
rz(-0.45845505) q[2];
sx q[2];
rz(3.0205883) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.0088494) q[1];
sx q[1];
rz(-0.64450507) q[1];
sx q[1];
rz(-0.82303859) q[1];
rz(-0.26975694) q[3];
sx q[3];
rz(-1.4175102) q[3];
sx q[3];
rz(-2.3333244) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.6404746) q[2];
sx q[2];
rz(-1.2299808) q[2];
sx q[2];
rz(-0.099700363) q[2];
rz(0.68263188) q[3];
sx q[3];
rz(-1.8004902) q[3];
sx q[3];
rz(2.7734267) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8354928) q[0];
sx q[0];
rz(-1.9266799) q[0];
sx q[0];
rz(0.039435506) q[0];
rz(-0.71190747) q[1];
sx q[1];
rz(-2.6075677) q[1];
sx q[1];
rz(1.8338411) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7402504) q[0];
sx q[0];
rz(-1.611705) q[0];
sx q[0];
rz(1.3274819) q[0];
rz(-pi) q[1];
rz(2.4277923) q[2];
sx q[2];
rz(-0.5792743) q[2];
sx q[2];
rz(-0.98503056) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.75406) q[1];
sx q[1];
rz(-2.3002831) q[1];
sx q[1];
rz(2.1977192) q[1];
rz(-pi) q[2];
x q[2];
rz(0.95473632) q[3];
sx q[3];
rz(-2.3619485) q[3];
sx q[3];
rz(0.35552939) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.51440614) q[2];
sx q[2];
rz(-1.4428978) q[2];
sx q[2];
rz(-2.8241039) q[2];
rz(-2.9182331) q[3];
sx q[3];
rz(-0.76627365) q[3];
sx q[3];
rz(1.4661695) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.8873972) q[0];
sx q[0];
rz(-1.5187718) q[0];
sx q[0];
rz(1.9018824) q[0];
rz(0.29048911) q[1];
sx q[1];
rz(-1.5668642) q[1];
sx q[1];
rz(1.3715502) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0199832) q[0];
sx q[0];
rz(-1.0399315) q[0];
sx q[0];
rz(-1.4565827) q[0];
rz(-pi) q[1];
rz(-1.1251584) q[2];
sx q[2];
rz(-1.9547714) q[2];
sx q[2];
rz(1.9572891) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.54372245) q[1];
sx q[1];
rz(-0.67670435) q[1];
sx q[1];
rz(1.1110439) q[1];
rz(-pi) q[2];
rz(1.4441078) q[3];
sx q[3];
rz(-2.3605862) q[3];
sx q[3];
rz(0.26050887) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.394968) q[2];
sx q[2];
rz(-1.1019424) q[2];
sx q[2];
rz(1.9071707) q[2];
rz(1.6055936) q[3];
sx q[3];
rz(-2.3707844) q[3];
sx q[3];
rz(2.3710282) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6813554) q[0];
sx q[0];
rz(-1.3139895) q[0];
sx q[0];
rz(-0.47948691) q[0];
rz(1.208249) q[1];
sx q[1];
rz(-1.6069326) q[1];
sx q[1];
rz(1.23034) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2527178) q[0];
sx q[0];
rz(-1.6299913) q[0];
sx q[0];
rz(2.2857347) q[0];
rz(-pi) q[1];
x q[1];
rz(0.540347) q[2];
sx q[2];
rz(-1.6896938) q[2];
sx q[2];
rz(-0.24248294) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.8248511) q[1];
sx q[1];
rz(-1.4808286) q[1];
sx q[1];
rz(1.852067) q[1];
x q[2];
rz(2.2097384) q[3];
sx q[3];
rz(-0.24792519) q[3];
sx q[3];
rz(2.7558079) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.6464107) q[2];
sx q[2];
rz(-2.9718706) q[2];
sx q[2];
rz(1.0037615) q[2];
rz(1.1209925) q[3];
sx q[3];
rz(-1.6206154) q[3];
sx q[3];
rz(0.77747548) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
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
rz(0.12370387) q[0];
sx q[0];
rz(-0.62010354) q[0];
sx q[0];
rz(-1.7271127) q[0];
rz(1.428528) q[1];
sx q[1];
rz(-2.2500549) q[1];
sx q[1];
rz(-0.052065484) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49567006) q[0];
sx q[0];
rz(-2.033253) q[0];
sx q[0];
rz(0.6301737) q[0];
x q[1];
rz(-2.2744479) q[2];
sx q[2];
rz(-2.646253) q[2];
sx q[2];
rz(0.42967202) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.19329883) q[1];
sx q[1];
rz(-1.5806075) q[1];
sx q[1];
rz(1.7722036) q[1];
x q[2];
rz(-2.089813) q[3];
sx q[3];
rz(-1.4148757) q[3];
sx q[3];
rz(-1.5633769) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.97466126) q[2];
sx q[2];
rz(-2.8278246) q[2];
sx q[2];
rz(-1.3775919) q[2];
rz(0.20197955) q[3];
sx q[3];
rz(-1.4429251) q[3];
sx q[3];
rz(-2.0418237) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5977151) q[0];
sx q[0];
rz(-1.7622204) q[0];
sx q[0];
rz(2.5482225) q[0];
rz(1.6186391) q[1];
sx q[1];
rz(-0.41888371) q[1];
sx q[1];
rz(0.12817344) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3421203) q[0];
sx q[0];
rz(-1.4065259) q[0];
sx q[0];
rz(1.4273433) q[0];
rz(-pi) q[1];
rz(0.32291193) q[2];
sx q[2];
rz(-0.92323136) q[2];
sx q[2];
rz(1.7932359) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.3230774) q[1];
sx q[1];
rz(-1.0836731) q[1];
sx q[1];
rz(0.088177322) q[1];
rz(-pi) q[2];
x q[2];
rz(0.45937755) q[3];
sx q[3];
rz(-0.93139) q[3];
sx q[3];
rz(2.3871445) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.0494978) q[2];
sx q[2];
rz(-2.8785661) q[2];
sx q[2];
rz(0.63801873) q[2];
rz(1.1370031) q[3];
sx q[3];
rz(-0.59320366) q[3];
sx q[3];
rz(3.1021281) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.20208246) q[0];
sx q[0];
rz(-1.7333637) q[0];
sx q[0];
rz(2.4092578) q[0];
rz(-2.9843075) q[1];
sx q[1];
rz(-1.597225) q[1];
sx q[1];
rz(1.0953915) q[1];
rz(-1.1284053) q[2];
sx q[2];
rz(-0.55991713) q[2];
sx q[2];
rz(-1.9922064) q[2];
rz(-0.19098115) q[3];
sx q[3];
rz(-0.92201391) q[3];
sx q[3];
rz(-1.3447522) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
