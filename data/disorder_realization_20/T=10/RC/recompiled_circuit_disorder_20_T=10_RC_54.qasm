OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.5946755) q[0];
sx q[0];
rz(-1.0008873) q[0];
sx q[0];
rz(2.9291908) q[0];
rz(-2.4266333) q[1];
sx q[1];
rz(-0.7874878) q[1];
sx q[1];
rz(-1.2815055) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40385383) q[0];
sx q[0];
rz(-1.3699342) q[0];
sx q[0];
rz(1.2234729) q[0];
x q[1];
rz(1.6002866) q[2];
sx q[2];
rz(-2.2631096) q[2];
sx q[2];
rz(-2.8147547) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.9077397) q[1];
sx q[1];
rz(-0.45957652) q[1];
sx q[1];
rz(1.9563667) q[1];
x q[2];
rz(-2.0251861) q[3];
sx q[3];
rz(-2.9985399) q[3];
sx q[3];
rz(-1.5798626) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.32221258) q[2];
sx q[2];
rz(-3.0266422) q[2];
sx q[2];
rz(-1.5167351) q[2];
rz(2.8939698) q[3];
sx q[3];
rz(-1.473749) q[3];
sx q[3];
rz(-2.4860399) q[3];
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
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6973998) q[0];
sx q[0];
rz(-1.6635165) q[0];
sx q[0];
rz(0.15629388) q[0];
rz(-0.63931757) q[1];
sx q[1];
rz(-2.1289861) q[1];
sx q[1];
rz(1.4888391) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.50870897) q[0];
sx q[0];
rz(-1.7299621) q[0];
sx q[0];
rz(2.5863618) q[0];
rz(-pi) q[1];
x q[1];
rz(0.82350125) q[2];
sx q[2];
rz(-0.89580065) q[2];
sx q[2];
rz(-1.0456955) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.3540346) q[1];
sx q[1];
rz(-2.2984142) q[1];
sx q[1];
rz(2.3393199) q[1];
x q[2];
rz(3.00499) q[3];
sx q[3];
rz(-1.6009637) q[3];
sx q[3];
rz(2.4581916) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.8504101) q[2];
sx q[2];
rz(-0.098878421) q[2];
sx q[2];
rz(2.8448811) q[2];
rz(0.3324278) q[3];
sx q[3];
rz(-1.5184831) q[3];
sx q[3];
rz(-1.9077574) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6333703) q[0];
sx q[0];
rz(-1.8828266) q[0];
sx q[0];
rz(2.6112774) q[0];
rz(-2.5976394) q[1];
sx q[1];
rz(-1.1214316) q[1];
sx q[1];
rz(-1.9525607) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2169164) q[0];
sx q[0];
rz(-0.33736704) q[0];
sx q[0];
rz(-0.49305537) q[0];
rz(-pi) q[1];
rz(-2.9844935) q[2];
sx q[2];
rz(-0.25114775) q[2];
sx q[2];
rz(-2.7709393) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.6310198) q[1];
sx q[1];
rz(-2.2377308) q[1];
sx q[1];
rz(-2.2651947) q[1];
rz(-pi) q[2];
rz(-2.6044106) q[3];
sx q[3];
rz(-2.3677285) q[3];
sx q[3];
rz(2.8230132) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.236078) q[2];
sx q[2];
rz(-1.4730075) q[2];
sx q[2];
rz(1.2515602) q[2];
rz(-0.45423147) q[3];
sx q[3];
rz(-0.43332228) q[3];
sx q[3];
rz(0.35513487) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.77804756) q[0];
sx q[0];
rz(-0.30968928) q[0];
sx q[0];
rz(-1.3409412) q[0];
rz(-0.60607934) q[1];
sx q[1];
rz(-0.93901912) q[1];
sx q[1];
rz(-1.1846503) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9979447) q[0];
sx q[0];
rz(-0.27772433) q[0];
sx q[0];
rz(-2.4585637) q[0];
rz(2.0834288) q[2];
sx q[2];
rz(-1.9789654) q[2];
sx q[2];
rz(0.689091) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.9222316) q[1];
sx q[1];
rz(-0.9833828) q[1];
sx q[1];
rz(0.14031336) q[1];
x q[2];
rz(1.9373478) q[3];
sx q[3];
rz(-1.9703715) q[3];
sx q[3];
rz(-1.2031581) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.9081395) q[2];
sx q[2];
rz(-0.88580695) q[2];
sx q[2];
rz(-2.502029) q[2];
rz(0.8574287) q[3];
sx q[3];
rz(-1.9968417) q[3];
sx q[3];
rz(-2.6517984) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5109167) q[0];
sx q[0];
rz(-0.54172051) q[0];
sx q[0];
rz(2.6211991) q[0];
rz(-3.0386472) q[1];
sx q[1];
rz(-1.0504477) q[1];
sx q[1];
rz(0.72881126) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9279014) q[0];
sx q[0];
rz(-1.7753082) q[0];
sx q[0];
rz(-1.4707028) q[0];
rz(-pi) q[1];
rz(-0.8466709) q[2];
sx q[2];
rz(-1.4482822) q[2];
sx q[2];
rz(-0.12144897) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.55620533) q[1];
sx q[1];
rz(-2.7558748) q[1];
sx q[1];
rz(-1.0005887) q[1];
rz(-pi) q[2];
rz(-0.95727386) q[3];
sx q[3];
rz(-1.1144708) q[3];
sx q[3];
rz(-0.43382713) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.0430498) q[2];
sx q[2];
rz(-2.7816732) q[2];
sx q[2];
rz(2.4318802) q[2];
rz(2.5085311) q[3];
sx q[3];
rz(-1.148162) q[3];
sx q[3];
rz(3.0098797) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1028041) q[0];
sx q[0];
rz(-3.0392201) q[0];
sx q[0];
rz(0.88371712) q[0];
rz(2.2209514) q[1];
sx q[1];
rz(-2.0795152) q[1];
sx q[1];
rz(-2.9439435) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5716498) q[0];
sx q[0];
rz(-2.2075704) q[0];
sx q[0];
rz(-2.3417579) q[0];
rz(-pi) q[1];
rz(3.0350787) q[2];
sx q[2];
rz(-1.4918367) q[2];
sx q[2];
rz(-0.67895141) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.3870961) q[1];
sx q[1];
rz(-1.5281786) q[1];
sx q[1];
rz(-2.6864762) q[1];
rz(-0.71182735) q[3];
sx q[3];
rz(-2.2440352) q[3];
sx q[3];
rz(-0.87513808) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-3.1137696) q[2];
sx q[2];
rz(-1.343507) q[2];
sx q[2];
rz(-3.0659884) q[2];
rz(1.4893701) q[3];
sx q[3];
rz(-0.39947709) q[3];
sx q[3];
rz(-2.3125125) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.41912115) q[0];
sx q[0];
rz(-2.2561181) q[0];
sx q[0];
rz(0.042073123) q[0];
rz(-1.178859) q[1];
sx q[1];
rz(-1.983164) q[1];
sx q[1];
rz(1.1694318) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5116611) q[0];
sx q[0];
rz(-0.17782623) q[0];
sx q[0];
rz(-1.1849665) q[0];
rz(-1.2866227) q[2];
sx q[2];
rz(-1.1948164) q[2];
sx q[2];
rz(0.77501955) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.3010901) q[1];
sx q[1];
rz(-1.2355348) q[1];
sx q[1];
rz(0.5718949) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.38798214) q[3];
sx q[3];
rz(-1.531732) q[3];
sx q[3];
rz(-0.29400533) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.3183257) q[2];
sx q[2];
rz(-1.1151168) q[2];
sx q[2];
rz(-0.014483359) q[2];
rz(-1.621834) q[3];
sx q[3];
rz(-1.2967916) q[3];
sx q[3];
rz(2.5884132) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0018205) q[0];
sx q[0];
rz(-2.791239) q[0];
sx q[0];
rz(2.5792504) q[0];
rz(-1.4315804) q[1];
sx q[1];
rz(-1.3969914) q[1];
sx q[1];
rz(0.45773488) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6477752) q[0];
sx q[0];
rz(-1.4871948) q[0];
sx q[0];
rz(1.093822) q[0];
rz(-pi) q[1];
rz(-2.2675603) q[2];
sx q[2];
rz(-1.1296141) q[2];
sx q[2];
rz(-0.87768427) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.8843958) q[1];
sx q[1];
rz(-0.83409062) q[1];
sx q[1];
rz(-0.019330545) q[1];
rz(-pi) q[2];
rz(1.0192972) q[3];
sx q[3];
rz(-2.5428452) q[3];
sx q[3];
rz(0.37125722) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.3056425) q[2];
sx q[2];
rz(-0.93985158) q[2];
sx q[2];
rz(-2.8590554) q[2];
rz(0.95747581) q[3];
sx q[3];
rz(-1.7493533) q[3];
sx q[3];
rz(-2.6628475) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1829421) q[0];
sx q[0];
rz(-2.4562953) q[0];
sx q[0];
rz(-1.4022934) q[0];
rz(2.4422586) q[1];
sx q[1];
rz(-1.8383086) q[1];
sx q[1];
rz(1.2618077) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48755074) q[0];
sx q[0];
rz(-2.6267509) q[0];
sx q[0];
rz(0.57060711) q[0];
rz(-2.8724573) q[2];
sx q[2];
rz(-0.67101523) q[2];
sx q[2];
rz(-3.105643) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.3952179) q[1];
sx q[1];
rz(-2.4002541) q[1];
sx q[1];
rz(0.37903255) q[1];
x q[2];
rz(-0.52957876) q[3];
sx q[3];
rz(-1.3929318) q[3];
sx q[3];
rz(1.9999268) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.4385779) q[2];
sx q[2];
rz(-1.5038749) q[2];
sx q[2];
rz(-1.4578488) q[2];
rz(-2.4750989) q[3];
sx q[3];
rz(-1.4278744) q[3];
sx q[3];
rz(2.3898747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89501971) q[0];
sx q[0];
rz(-0.81037766) q[0];
sx q[0];
rz(1.1599468) q[0];
rz(3.0874522) q[1];
sx q[1];
rz(-1.663798) q[1];
sx q[1];
rz(2.0711526) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.69660891) q[0];
sx q[0];
rz(-2.3009926) q[0];
sx q[0];
rz(0.22459774) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.1033377) q[2];
sx q[2];
rz(-2.3352211) q[2];
sx q[2];
rz(2.5753266) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.66026238) q[1];
sx q[1];
rz(-0.80033014) q[1];
sx q[1];
rz(-0.34923133) q[1];
rz(2.8937267) q[3];
sx q[3];
rz(-1.7409179) q[3];
sx q[3];
rz(0.47249139) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.3728309) q[2];
sx q[2];
rz(-2.5925345) q[2];
sx q[2];
rz(1.5268415) q[2];
rz(0.57957831) q[3];
sx q[3];
rz(-1.3468346) q[3];
sx q[3];
rz(2.4894864) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51331818) q[0];
sx q[0];
rz(-1.6074629) q[0];
sx q[0];
rz(0.71832023) q[0];
rz(-1.0340446) q[1];
sx q[1];
rz(-0.68987344) q[1];
sx q[1];
rz(-0.72763163) q[1];
rz(1.5288011) q[2];
sx q[2];
rz(-0.73245807) q[2];
sx q[2];
rz(3.0277638) q[2];
rz(-1.0735687) q[3];
sx q[3];
rz(-1.0931001) q[3];
sx q[3];
rz(1.1563374) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];