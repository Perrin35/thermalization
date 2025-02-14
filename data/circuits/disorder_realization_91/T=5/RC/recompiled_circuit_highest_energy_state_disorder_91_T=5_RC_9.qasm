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
rz(-2.1639106) q[1];
sx q[1];
rz(-0.77562195) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6399423) q[0];
sx q[0];
rz(-1.8193366) q[0];
sx q[0];
rz(-2.2959501) q[0];
x q[1];
rz(0.94812265) q[2];
sx q[2];
rz(-0.97062696) q[2];
sx q[2];
rz(2.4400939) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.089701414) q[1];
sx q[1];
rz(-2.6503453) q[1];
sx q[1];
rz(0.516275) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5585639) q[3];
sx q[3];
rz(-2.1656007) q[3];
sx q[3];
rz(-2.6786238) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.97950196) q[2];
sx q[2];
rz(-2.4544331) q[2];
sx q[2];
rz(0.053357601) q[2];
rz(0.067367628) q[3];
sx q[3];
rz(-0.24240436) q[3];
sx q[3];
rz(1.4343725) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1783717) q[0];
sx q[0];
rz(-1.0472949) q[0];
sx q[0];
rz(2.6317327) q[0];
rz(2.0193822) q[1];
sx q[1];
rz(-2.0041859) q[1];
sx q[1];
rz(-3.1013464) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7644234) q[0];
sx q[0];
rz(-0.80028811) q[0];
sx q[0];
rz(-1.1616526) q[0];
rz(-1.5284677) q[2];
sx q[2];
rz(-0.63708239) q[2];
sx q[2];
rz(0.55227623) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.65529167) q[1];
sx q[1];
rz(-2.5319063) q[1];
sx q[1];
rz(2.8279081) q[1];
rz(-pi) q[2];
x q[2];
rz(0.67701234) q[3];
sx q[3];
rz(-1.5843862) q[3];
sx q[3];
rz(-0.06448596) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.76790729) q[2];
sx q[2];
rz(-1.3482956) q[2];
sx q[2];
rz(2.574004) q[2];
rz(-0.23809412) q[3];
sx q[3];
rz(-1.4990436) q[3];
sx q[3];
rz(-2.8897918) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(3.0955065) q[0];
sx q[0];
rz(-2.6987785) q[0];
sx q[0];
rz(-2.1777022) q[0];
rz(1.6148199) q[1];
sx q[1];
rz(-1.2232199) q[1];
sx q[1];
rz(2.8098106) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.362865) q[0];
sx q[0];
rz(-2.300867) q[0];
sx q[0];
rz(-1.8018468) q[0];
rz(-pi) q[1];
rz(-0.024130017) q[2];
sx q[2];
rz(-1.2306613) q[2];
sx q[2];
rz(0.4235776) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.11321924) q[1];
sx q[1];
rz(-3.0455596) q[1];
sx q[1];
rz(-1.5043831) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2778707) q[3];
sx q[3];
rz(-0.4583685) q[3];
sx q[3];
rz(1.019695) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.31384599) q[2];
sx q[2];
rz(-1.4705986) q[2];
sx q[2];
rz(0.016679114) q[2];
rz(-1.5879177) q[3];
sx q[3];
rz(-1.1976306) q[3];
sx q[3];
rz(-2.8488819) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68506771) q[0];
sx q[0];
rz(-2.2790788) q[0];
sx q[0];
rz(-0.21441329) q[0];
rz(0.21687493) q[1];
sx q[1];
rz(-2.3213991) q[1];
sx q[1];
rz(-2.0490501) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8583243) q[0];
sx q[0];
rz(-0.28139925) q[0];
sx q[0];
rz(-0.94080956) q[0];
rz(-pi) q[1];
rz(-0.09437807) q[2];
sx q[2];
rz(-0.98446956) q[2];
sx q[2];
rz(-2.6142769) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.1779611) q[1];
sx q[1];
rz(-0.73252618) q[1];
sx q[1];
rz(1.5103106) q[1];
rz(-pi) q[2];
rz(-0.86534604) q[3];
sx q[3];
rz(-1.0320889) q[3];
sx q[3];
rz(-2.5821843) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.0715926) q[2];
sx q[2];
rz(-2.6061974) q[2];
sx q[2];
rz(-2.6154321) q[2];
rz(1.5782662) q[3];
sx q[3];
rz(-1.1514781) q[3];
sx q[3];
rz(0.17939849) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.18463317) q[0];
sx q[0];
rz(-1.7687836) q[0];
sx q[0];
rz(-2.2233295) q[0];
rz(2.1271162) q[1];
sx q[1];
rz(-2.4160073) q[1];
sx q[1];
rz(2.0065506) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.710478) q[0];
sx q[0];
rz(-1.6621792) q[0];
sx q[0];
rz(2.8764832) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1447976) q[2];
sx q[2];
rz(-1.3959178) q[2];
sx q[2];
rz(1.3256734) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.1327433) q[1];
sx q[1];
rz(-2.4970876) q[1];
sx q[1];
rz(0.82303859) q[1];
rz(-pi) q[2];
rz(1.7297392) q[3];
sx q[3];
rz(-1.30428) q[3];
sx q[3];
rz(2.4212568) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.501118) q[2];
sx q[2];
rz(-1.2299808) q[2];
sx q[2];
rz(-3.0418923) q[2];
rz(-2.4589608) q[3];
sx q[3];
rz(-1.8004902) q[3];
sx q[3];
rz(-0.36816594) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8354928) q[0];
sx q[0];
rz(-1.9266799) q[0];
sx q[0];
rz(3.1021571) q[0];
rz(-2.4296852) q[1];
sx q[1];
rz(-2.6075677) q[1];
sx q[1];
rz(1.3077516) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9619869) q[0];
sx q[0];
rz(-1.8139031) q[0];
sx q[0];
rz(0.042148729) q[0];
rz(-1.9754312) q[2];
sx q[2];
rz(-1.9973997) q[2];
sx q[2];
rz(-0.18243039) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.75406) q[1];
sx q[1];
rz(-2.3002831) q[1];
sx q[1];
rz(-2.1977192) q[1];
rz(-0.51898308) q[3];
sx q[3];
rz(-2.1819138) q[3];
sx q[3];
rz(0.42762363) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.51440614) q[2];
sx q[2];
rz(-1.4428978) q[2];
sx q[2];
rz(0.31748873) q[2];
rz(-0.22335957) q[3];
sx q[3];
rz(-0.76627365) q[3];
sx q[3];
rz(-1.4661695) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2541955) q[0];
sx q[0];
rz(-1.5187718) q[0];
sx q[0];
rz(-1.9018824) q[0];
rz(2.8511035) q[1];
sx q[1];
rz(-1.5747285) q[1];
sx q[1];
rz(1.3715502) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0199832) q[0];
sx q[0];
rz(-1.0399315) q[0];
sx q[0];
rz(-1.68501) q[0];
rz(0.81772501) q[2];
sx q[2];
rz(-2.5618988) q[2];
sx q[2];
rz(-2.8630437) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.54372245) q[1];
sx q[1];
rz(-0.67670435) q[1];
sx q[1];
rz(2.0305487) q[1];
x q[2];
rz(3.0169964) q[3];
sx q[3];
rz(-2.3438896) q[3];
sx q[3];
rz(-0.083051894) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4602373) q[0];
sx q[0];
rz(-1.8276031) q[0];
sx q[0];
rz(-2.6621057) q[0];
rz(1.208249) q[1];
sx q[1];
rz(-1.5346601) q[1];
sx q[1];
rz(-1.23034) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.63060627) q[0];
sx q[0];
rz(-0.85737757) q[0];
sx q[0];
rz(0.078321266) q[0];
rz(1.7092136) q[2];
sx q[2];
rz(-2.1069134) q[2];
sx q[2];
rz(1.7422402) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.8248511) q[1];
sx q[1];
rz(-1.6607641) q[1];
sx q[1];
rz(1.852067) q[1];
rz(0.9318542) q[3];
sx q[3];
rz(-0.24792519) q[3];
sx q[3];
rz(-2.7558079) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.6464107) q[2];
sx q[2];
rz(-0.16972204) q[2];
sx q[2];
rz(1.0037615) q[2];
rz(2.0206001) q[3];
sx q[3];
rz(-1.5209773) q[3];
sx q[3];
rz(-2.3641172) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.12370387) q[0];
sx q[0];
rz(-0.62010354) q[0];
sx q[0];
rz(1.7271127) q[0];
rz(-1.7130647) q[1];
sx q[1];
rz(-2.2500549) q[1];
sx q[1];
rz(3.0895272) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.751862) q[0];
sx q[0];
rz(-1.0152752) q[0];
sx q[0];
rz(1.0179505) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.8671448) q[2];
sx q[2];
rz(-2.646253) q[2];
sx q[2];
rz(2.7119206) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.4255124) q[1];
sx q[1];
rz(-2.9399498) q[1];
sx q[1];
rz(-1.6198026) q[1];
x q[2];
rz(2.089813) q[3];
sx q[3];
rz(-1.4148757) q[3];
sx q[3];
rz(-1.5782158) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.1669314) q[2];
sx q[2];
rz(-2.8278246) q[2];
sx q[2];
rz(1.7640007) q[2];
rz(-0.20197955) q[3];
sx q[3];
rz(-1.4429251) q[3];
sx q[3];
rz(-1.099769) q[3];
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
x q[0];
rz(pi/2) q[0];
x q[1];
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
rz(-1.6186391) q[1];
sx q[1];
rz(-2.7227089) q[1];
sx q[1];
rz(0.12817344) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0656212) q[0];
sx q[0];
rz(-2.9239281) q[0];
sx q[0];
rz(-0.71167167) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9680357) q[2];
sx q[2];
rz(-2.4284869) q[2];
sx q[2];
rz(-1.8548059) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.1364079) q[1];
sx q[1];
rz(-2.6471849) q[1];
sx q[1];
rz(1.4060518) q[1];
rz(2.6822151) q[3];
sx q[3];
rz(-2.2102027) q[3];
sx q[3];
rz(2.3871445) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.0920948) q[2];
sx q[2];
rz(-0.2630266) q[2];
sx q[2];
rz(0.63801873) q[2];
rz(-2.0045896) q[3];
sx q[3];
rz(-2.548389) q[3];
sx q[3];
rz(0.039464522) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20208246) q[0];
sx q[0];
rz(-1.408229) q[0];
sx q[0];
rz(-0.73233488) q[0];
rz(-2.9843075) q[1];
sx q[1];
rz(-1.597225) q[1];
sx q[1];
rz(1.0953915) q[1];
rz(-2.8794206) q[2];
sx q[2];
rz(-2.0714348) q[2];
sx q[2];
rz(-2.5019912) q[2];
rz(2.2284343) q[3];
sx q[3];
rz(-1.4189594) q[3];
sx q[3];
rz(0.34233477) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
