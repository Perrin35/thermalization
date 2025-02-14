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
rz(-2.8862267) q[0];
sx q[0];
rz(-2.3998883) q[0];
sx q[0];
rz(-1.1852888) q[0];
rz(1.9396632) q[1];
sx q[1];
rz(-0.97768205) q[1];
sx q[1];
rz(-2.3659707) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.14549832) q[0];
sx q[0];
rz(-0.87255635) q[0];
sx q[0];
rz(-2.8146539) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.19347) q[2];
sx q[2];
rz(-2.1709657) q[2];
sx q[2];
rz(-2.4400939) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.089701414) q[1];
sx q[1];
rz(-2.6503453) q[1];
sx q[1];
rz(-2.6253176) q[1];
rz(-pi) q[2];
rz(3.1235141) q[3];
sx q[3];
rz(-2.5466777) q[3];
sx q[3];
rz(0.48479652) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.97950196) q[2];
sx q[2];
rz(-2.4544331) q[2];
sx q[2];
rz(-3.0882351) q[2];
rz(-3.074225) q[3];
sx q[3];
rz(-2.8991883) q[3];
sx q[3];
rz(-1.4343725) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1783717) q[0];
sx q[0];
rz(-2.0942978) q[0];
sx q[0];
rz(0.50985992) q[0];
rz(2.0193822) q[1];
sx q[1];
rz(-1.1374067) q[1];
sx q[1];
rz(-0.040246211) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7644234) q[0];
sx q[0];
rz(-2.3413045) q[0];
sx q[0];
rz(-1.9799401) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.93414238) q[2];
sx q[2];
rz(-1.5959708) q[2];
sx q[2];
rz(-2.0890401) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.65529167) q[1];
sx q[1];
rz(-2.5319063) q[1];
sx q[1];
rz(0.31368452) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.021690856) q[3];
sx q[3];
rz(-2.4644654) q[3];
sx q[3];
rz(1.6521887) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.76790729) q[2];
sx q[2];
rz(-1.3482956) q[2];
sx q[2];
rz(2.574004) q[2];
rz(-2.9034985) q[3];
sx q[3];
rz(-1.6425491) q[3];
sx q[3];
rz(0.25180086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0955065) q[0];
sx q[0];
rz(-2.6987785) q[0];
sx q[0];
rz(2.1777022) q[0];
rz(-1.6148199) q[1];
sx q[1];
rz(-1.9183728) q[1];
sx q[1];
rz(2.8098106) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7787276) q[0];
sx q[0];
rz(-0.84072564) q[0];
sx q[0];
rz(1.3397459) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9110229) q[2];
sx q[2];
rz(-1.5935437) q[2];
sx q[2];
rz(1.1391672) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-3.0283734) q[1];
sx q[1];
rz(-0.096033022) q[1];
sx q[1];
rz(-1.6372096) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0000684) q[3];
sx q[3];
rz(-1.1333395) q[3];
sx q[3];
rz(0.69526068) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.31384599) q[2];
sx q[2];
rz(-1.4705986) q[2];
sx q[2];
rz(3.1249135) q[2];
rz(-1.5879177) q[3];
sx q[3];
rz(-1.1976306) q[3];
sx q[3];
rz(-2.8488819) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4565249) q[0];
sx q[0];
rz(-0.86251384) q[0];
sx q[0];
rz(-0.21441329) q[0];
rz(-2.9247177) q[1];
sx q[1];
rz(-2.3213991) q[1];
sx q[1];
rz(1.0925426) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.93247) q[0];
sx q[0];
rz(-1.7971149) q[0];
sx q[0];
rz(-2.97291) q[0];
x q[1];
rz(-2.1591805) q[2];
sx q[2];
rz(-1.6493754) q[2];
sx q[2];
rz(2.1504398) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.34784029) q[1];
sx q[1];
rz(-1.6112325) q[1];
sx q[1];
rz(-2.3024125) q[1];
rz(0.86534604) q[3];
sx q[3];
rz(-2.1095037) q[3];
sx q[3];
rz(-2.5821843) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.07) q[2];
sx q[2];
rz(-0.53539521) q[2];
sx q[2];
rz(-0.52616057) q[2];
rz(1.5782662) q[3];
sx q[3];
rz(-1.9901146) q[3];
sx q[3];
rz(-0.17939849) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.18463317) q[0];
sx q[0];
rz(-1.7687836) q[0];
sx q[0];
rz(0.9182632) q[0];
rz(-2.1271162) q[1];
sx q[1];
rz(-2.4160073) q[1];
sx q[1];
rz(-2.0065506) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1149108) q[0];
sx q[0];
rz(-1.8347731) q[0];
sx q[0];
rz(-1.6654679) q[0];
rz(1.9967951) q[2];
sx q[2];
rz(-1.3959178) q[2];
sx q[2];
rz(1.3256734) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.1327433) q[1];
sx q[1];
rz(-0.64450507) q[1];
sx q[1];
rz(-0.82303859) q[1];
rz(-1.7297392) q[3];
sx q[3];
rz(-1.30428) q[3];
sx q[3];
rz(-2.4212568) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.6404746) q[2];
sx q[2];
rz(-1.2299808) q[2];
sx q[2];
rz(3.0418923) q[2];
rz(2.4589608) q[3];
sx q[3];
rz(-1.8004902) q[3];
sx q[3];
rz(-2.7734267) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8354928) q[0];
sx q[0];
rz(-1.2149128) q[0];
sx q[0];
rz(0.039435506) q[0];
rz(2.4296852) q[1];
sx q[1];
rz(-0.53402495) q[1];
sx q[1];
rz(-1.8338411) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9619869) q[0];
sx q[0];
rz(-1.8139031) q[0];
sx q[0];
rz(3.0994439) q[0];
rz(-pi) q[1];
rz(2.6823993) q[2];
sx q[2];
rz(-1.9373496) q[2];
sx q[2];
rz(-1.2129779) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.2145526) q[1];
sx q[1];
rz(-0.92260375) q[1];
sx q[1];
rz(-2.560858) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.51898308) q[3];
sx q[3];
rz(-0.95967889) q[3];
sx q[3];
rz(2.713969) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.51440614) q[2];
sx q[2];
rz(-1.4428978) q[2];
sx q[2];
rz(2.8241039) q[2];
rz(0.22335957) q[3];
sx q[3];
rz(-2.375319) q[3];
sx q[3];
rz(1.6754231) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2541955) q[0];
sx q[0];
rz(-1.6228209) q[0];
sx q[0];
rz(-1.2397102) q[0];
rz(0.29048911) q[1];
sx q[1];
rz(-1.5747285) q[1];
sx q[1];
rz(-1.3715502) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1216095) q[0];
sx q[0];
rz(-1.0399315) q[0];
sx q[0];
rz(-1.68501) q[0];
x q[1];
rz(2.0164343) q[2];
sx q[2];
rz(-1.9547714) q[2];
sx q[2];
rz(-1.1843036) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.1095431) q[1];
sx q[1];
rz(-0.97496009) q[1];
sx q[1];
rz(0.34237564) q[1];
rz(-1.6974849) q[3];
sx q[3];
rz(-0.7810065) q[3];
sx q[3];
rz(2.8810838) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.74662465) q[2];
sx q[2];
rz(-2.0396502) q[2];
sx q[2];
rz(-1.234422) q[2];
rz(1.5359991) q[3];
sx q[3];
rz(-2.3707844) q[3];
sx q[3];
rz(0.7705645) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4602373) q[0];
sx q[0];
rz(-1.3139895) q[0];
sx q[0];
rz(2.6621057) q[0];
rz(-1.208249) q[1];
sx q[1];
rz(-1.5346601) q[1];
sx q[1];
rz(1.23034) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.63060627) q[0];
sx q[0];
rz(-0.85737757) q[0];
sx q[0];
rz(0.078321266) q[0];
rz(-pi) q[1];
rz(1.7092136) q[2];
sx q[2];
rz(-2.1069134) q[2];
sx q[2];
rz(1.7422402) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.047428377) q[1];
sx q[1];
rz(-0.29494527) q[1];
sx q[1];
rz(-1.2565681) q[1];
x q[2];
rz(-0.9318542) q[3];
sx q[3];
rz(-2.8936675) q[3];
sx q[3];
rz(0.38578472) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.6464107) q[2];
sx q[2];
rz(-2.9718706) q[2];
sx q[2];
rz(-2.1378311) q[2];
rz(2.0206001) q[3];
sx q[3];
rz(-1.5209773) q[3];
sx q[3];
rz(-2.3641172) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0178888) q[0];
sx q[0];
rz(-2.5214891) q[0];
sx q[0];
rz(1.7271127) q[0];
rz(-1.428528) q[1];
sx q[1];
rz(-0.89153779) q[1];
sx q[1];
rz(3.0895272) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49567006) q[0];
sx q[0];
rz(-1.1083397) q[0];
sx q[0];
rz(-2.511419) q[0];
x q[1];
rz(-1.9615575) q[2];
sx q[2];
rz(-1.8834049) q[2];
sx q[2];
rz(-1.7824204) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.9482938) q[1];
sx q[1];
rz(-1.5609852) q[1];
sx q[1];
rz(-1.7722036) q[1];
rz(-pi) q[2];
rz(2.089813) q[3];
sx q[3];
rz(-1.726717) q[3];
sx q[3];
rz(-1.5633769) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.97466126) q[2];
sx q[2];
rz(-2.8278246) q[2];
sx q[2];
rz(1.3775919) q[2];
rz(0.20197955) q[3];
sx q[3];
rz(-1.4429251) q[3];
sx q[3];
rz(-2.0418237) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5977151) q[0];
sx q[0];
rz(-1.3793722) q[0];
sx q[0];
rz(-2.5482225) q[0];
rz(1.5229535) q[1];
sx q[1];
rz(-0.41888371) q[1];
sx q[1];
rz(-0.12817344) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9365337) q[0];
sx q[0];
rz(-1.4292875) q[0];
sx q[0];
rz(-2.9756484) q[0];
rz(-0.32291193) q[2];
sx q[2];
rz(-2.2183613) q[2];
sx q[2];
rz(-1.3483568) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.81851527) q[1];
sx q[1];
rz(-1.0836731) q[1];
sx q[1];
rz(-0.088177322) q[1];
x q[2];
rz(2.6822151) q[3];
sx q[3];
rz(-0.93139) q[3];
sx q[3];
rz(-2.3871445) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.0920948) q[2];
sx q[2];
rz(-0.2630266) q[2];
sx q[2];
rz(-2.5035739) q[2];
rz(-1.1370031) q[3];
sx q[3];
rz(-2.548389) q[3];
sx q[3];
rz(-0.039464522) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.20208246) q[0];
sx q[0];
rz(-1.7333637) q[0];
sx q[0];
rz(2.4092578) q[0];
rz(0.15728514) q[1];
sx q[1];
rz(-1.597225) q[1];
sx q[1];
rz(1.0953915) q[1];
rz(1.0553817) q[2];
sx q[2];
rz(-1.800174) q[2];
sx q[2];
rz(-0.80309662) q[2];
rz(-0.91315837) q[3];
sx q[3];
rz(-1.4189594) q[3];
sx q[3];
rz(0.34233477) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
