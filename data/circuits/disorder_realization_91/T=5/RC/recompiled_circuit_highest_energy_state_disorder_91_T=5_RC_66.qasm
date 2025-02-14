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
rz(-1.2019295) q[1];
sx q[1];
rz(4.1192747) q[1];
sx q[1];
rz(8.649156) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14549832) q[0];
sx q[0];
rz(-2.2690363) q[0];
sx q[0];
rz(-2.8146539) q[0];
x q[1];
rz(-0.94812265) q[2];
sx q[2];
rz(-2.1709657) q[2];
sx q[2];
rz(2.4400939) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.47992) q[1];
sx q[1];
rz(-1.1480756) q[1];
sx q[1];
rz(-1.3125961) q[1];
rz(-pi) q[2];
rz(-0.01807853) q[3];
sx q[3];
rz(-2.5466777) q[3];
sx q[3];
rz(-2.6567961) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.97950196) q[2];
sx q[2];
rz(-0.68715954) q[2];
sx q[2];
rz(0.053357601) q[2];
rz(0.067367628) q[3];
sx q[3];
rz(-2.8991883) q[3];
sx q[3];
rz(-1.4343725) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1783717) q[0];
sx q[0];
rz(-2.0942978) q[0];
sx q[0];
rz(2.6317327) q[0];
rz(2.0193822) q[1];
sx q[1];
rz(-2.0041859) q[1];
sx q[1];
rz(-3.1013464) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.82032594) q[0];
sx q[0];
rz(-2.2893956) q[0];
sx q[0];
rz(0.38896968) q[0];
rz(-pi) q[1];
rz(-2.2074503) q[2];
sx q[2];
rz(-1.5456219) q[2];
sx q[2];
rz(1.0525525) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.8630835) q[1];
sx q[1];
rz(-0.99480226) q[1];
sx q[1];
rz(1.7830677) q[1];
rz(-pi) q[2];
rz(3.1199018) q[3];
sx q[3];
rz(-2.4644654) q[3];
sx q[3];
rz(1.6521887) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.3736854) q[2];
sx q[2];
rz(-1.7932971) q[2];
sx q[2];
rz(2.574004) q[2];
rz(0.23809412) q[3];
sx q[3];
rz(-1.6425491) q[3];
sx q[3];
rz(0.25180086) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.046086144) q[0];
sx q[0];
rz(-0.44281414) q[0];
sx q[0];
rz(-0.96389043) q[0];
rz(1.5267728) q[1];
sx q[1];
rz(-1.9183728) q[1];
sx q[1];
rz(-0.33178202) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0237494) q[0];
sx q[0];
rz(-0.75928771) q[0];
sx q[0];
rz(-2.8911126) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2305698) q[2];
sx q[2];
rz(-1.5935437) q[2];
sx q[2];
rz(-1.1391672) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.7501237) q[1];
sx q[1];
rz(-1.5771598) q[1];
sx q[1];
rz(1.4749737) q[1];
rz(-pi) q[2];
rz(-1.863722) q[3];
sx q[3];
rz(-0.4583685) q[3];
sx q[3];
rz(-1.019695) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.31384599) q[2];
sx q[2];
rz(-1.670994) q[2];
sx q[2];
rz(-0.016679114) q[2];
rz(1.5879177) q[3];
sx q[3];
rz(-1.943962) q[3];
sx q[3];
rz(-2.8488819) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4565249) q[0];
sx q[0];
rz(-0.86251384) q[0];
sx q[0];
rz(2.9271794) q[0];
rz(2.9247177) q[1];
sx q[1];
rz(-0.82019359) q[1];
sx q[1];
rz(-2.0490501) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2091227) q[0];
sx q[0];
rz(-1.7971149) q[0];
sx q[0];
rz(0.16868261) q[0];
rz(-pi) q[1];
rz(-2.1591805) q[2];
sx q[2];
rz(-1.4922172) q[2];
sx q[2];
rz(0.99115288) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.7937524) q[1];
sx q[1];
rz(-1.5303601) q[1];
sx q[1];
rz(2.3024125) q[1];
x q[2];
rz(-2.3155261) q[3];
sx q[3];
rz(-2.2829901) q[3];
sx q[3];
rz(-0.46921092) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.07) q[2];
sx q[2];
rz(-2.6061974) q[2];
sx q[2];
rz(0.52616057) q[2];
rz(1.5782662) q[3];
sx q[3];
rz(-1.1514781) q[3];
sx q[3];
rz(-2.9621942) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.18463317) q[0];
sx q[0];
rz(-1.3728091) q[0];
sx q[0];
rz(0.9182632) q[0];
rz(-1.0144764) q[1];
sx q[1];
rz(-0.72558534) q[1];
sx q[1];
rz(1.135042) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6776575) q[0];
sx q[0];
rz(-0.2800664) q[0];
sx q[0];
rz(-0.33644648) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9967951) q[2];
sx q[2];
rz(-1.7456749) q[2];
sx q[2];
rz(-1.3256734) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.8682485) q[1];
sx q[1];
rz(-1.1146208) q[1];
sx q[1];
rz(2.6691324) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6162063) q[3];
sx q[3];
rz(-2.8322517) q[3];
sx q[3];
rz(1.8745195) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.6404746) q[2];
sx q[2];
rz(-1.2299808) q[2];
sx q[2];
rz(0.099700363) q[2];
rz(0.68263188) q[3];
sx q[3];
rz(-1.8004902) q[3];
sx q[3];
rz(-0.36816594) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.30609983) q[0];
sx q[0];
rz(-1.2149128) q[0];
sx q[0];
rz(-3.1021571) q[0];
rz(-0.71190747) q[1];
sx q[1];
rz(-0.53402495) q[1];
sx q[1];
rz(-1.8338411) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1354254) q[0];
sx q[0];
rz(-2.8949304) q[0];
sx q[0];
rz(1.7390854) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.4591934) q[2];
sx q[2];
rz(-1.9373496) q[2];
sx q[2];
rz(-1.2129779) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.73346123) q[1];
sx q[1];
rz(-1.1181338) q[1];
sx q[1];
rz(2.3068417) q[1];
x q[2];
rz(2.2496836) q[3];
sx q[3];
rz(-1.1524767) q[3];
sx q[3];
rz(1.4598916) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.6271865) q[2];
sx q[2];
rz(-1.6986948) q[2];
sx q[2];
rz(0.31748873) q[2];
rz(0.22335957) q[3];
sx q[3];
rz(-0.76627365) q[3];
sx q[3];
rz(-1.6754231) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2541955) q[0];
sx q[0];
rz(-1.6228209) q[0];
sx q[0];
rz(-1.9018824) q[0];
rz(-0.29048911) q[1];
sx q[1];
rz(-1.5668642) q[1];
sx q[1];
rz(1.7700425) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6487911) q[0];
sx q[0];
rz(-1.472357) q[0];
sx q[0];
rz(-0.53372328) q[0];
rz(-pi) q[1];
x q[1];
rz(0.42098896) q[2];
sx q[2];
rz(-1.1596646) q[2];
sx q[2];
rz(0.56356424) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.4829258) q[1];
sx q[1];
rz(-1.2892168) q[1];
sx q[1];
rz(0.94689178) q[1];
rz(-pi) q[2];
rz(1.6974849) q[3];
sx q[3];
rz(-2.3605862) q[3];
sx q[3];
rz(2.8810838) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.394968) q[2];
sx q[2];
rz(-2.0396502) q[2];
sx q[2];
rz(1.9071707) q[2];
rz(1.5359991) q[3];
sx q[3];
rz(-2.3707844) q[3];
sx q[3];
rz(0.7705645) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-1.6813554) q[0];
sx q[0];
rz(-1.3139895) q[0];
sx q[0];
rz(0.47948691) q[0];
rz(-1.208249) q[1];
sx q[1];
rz(-1.6069326) q[1];
sx q[1];
rz(1.9112526) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5109864) q[0];
sx q[0];
rz(-2.2842151) q[0];
sx q[0];
rz(-0.078321266) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9134184) q[2];
sx q[2];
rz(-2.58959) q[2];
sx q[2];
rz(1.133103) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.8248511) q[1];
sx q[1];
rz(-1.4808286) q[1];
sx q[1];
rz(-1.852067) q[1];
rz(-2.9917688) q[3];
sx q[3];
rz(-1.3725159) q[3];
sx q[3];
rz(1.0396569) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.49518195) q[2];
sx q[2];
rz(-0.16972204) q[2];
sx q[2];
rz(-1.0037615) q[2];
rz(-1.1209925) q[3];
sx q[3];
rz(-1.6206154) q[3];
sx q[3];
rz(2.3641172) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0178888) q[0];
sx q[0];
rz(-0.62010354) q[0];
sx q[0];
rz(-1.7271127) q[0];
rz(1.7130647) q[1];
sx q[1];
rz(-2.2500549) q[1];
sx q[1];
rz(0.052065484) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3897306) q[0];
sx q[0];
rz(-1.0152752) q[0];
sx q[0];
rz(1.0179505) q[0];
rz(-pi) q[1];
rz(0.33627908) q[2];
sx q[2];
rz(-1.1999201) q[2];
sx q[2];
rz(-2.8039497) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.9482938) q[1];
sx q[1];
rz(-1.5806075) q[1];
sx q[1];
rz(1.3693891) q[1];
x q[2];
rz(-2.9624945) q[3];
sx q[3];
rz(-1.0586959) q[3];
sx q[3];
rz(0.095898151) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.97466126) q[2];
sx q[2];
rz(-0.31376803) q[2];
sx q[2];
rz(1.7640007) q[2];
rz(-2.9396131) q[3];
sx q[3];
rz(-1.6986676) q[3];
sx q[3];
rz(2.0418237) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5977151) q[0];
sx q[0];
rz(-1.3793722) q[0];
sx q[0];
rz(0.5933702) q[0];
rz(1.6186391) q[1];
sx q[1];
rz(-0.41888371) q[1];
sx q[1];
rz(-3.0134192) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9365337) q[0];
sx q[0];
rz(-1.4292875) q[0];
sx q[0];
rz(2.9756484) q[0];
rz(-pi) q[1];
x q[1];
rz(0.89752512) q[2];
sx q[2];
rz(-1.3149259) q[2];
sx q[2];
rz(-3.1183372) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.3230774) q[1];
sx q[1];
rz(-2.0579195) q[1];
sx q[1];
rz(-3.0534153) q[1];
rz(-pi) q[2];
rz(2.1084598) q[3];
sx q[3];
rz(-0.76806107) q[3];
sx q[3];
rz(-1.6949754) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.0920948) q[2];
sx q[2];
rz(-2.8785661) q[2];
sx q[2];
rz(-2.5035739) q[2];
rz(2.0045896) q[3];
sx q[3];
rz(-2.548389) q[3];
sx q[3];
rz(-0.039464522) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.20208246) q[0];
sx q[0];
rz(-1.7333637) q[0];
sx q[0];
rz(2.4092578) q[0];
rz(2.9843075) q[1];
sx q[1];
rz(-1.5443677) q[1];
sx q[1];
rz(-2.0462012) q[1];
rz(-2.8794206) q[2];
sx q[2];
rz(-2.0714348) q[2];
sx q[2];
rz(-2.5019912) q[2];
rz(-2.2284343) q[3];
sx q[3];
rz(-1.7226333) q[3];
sx q[3];
rz(-2.7992579) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
