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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8016937) q[0];
sx q[0];
rz(-2.3824132) q[0];
sx q[0];
rz(1.2053427) q[0];
rz(-pi) q[1];
rz(-2.4414659) q[2];
sx q[2];
rz(-1.0686734) q[2];
sx q[2];
rz(1.2545253) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.47992) q[1];
sx q[1];
rz(-1.993517) q[1];
sx q[1];
rz(-1.3125961) q[1];
x q[2];
rz(0.01807853) q[3];
sx q[3];
rz(-0.59491494) q[3];
sx q[3];
rz(-2.6567961) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.97950196) q[2];
sx q[2];
rz(-0.68715954) q[2];
sx q[2];
rz(-0.053357601) q[2];
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
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.963221) q[0];
sx q[0];
rz(-1.0472949) q[0];
sx q[0];
rz(0.50985992) q[0];
rz(-2.0193822) q[1];
sx q[1];
rz(-1.1374067) q[1];
sx q[1];
rz(-3.1013464) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3771693) q[0];
sx q[0];
rz(-0.80028811) q[0];
sx q[0];
rz(1.1616526) q[0];
rz(-3.1102883) q[2];
sx q[2];
rz(-0.93437663) q[2];
sx q[2];
rz(2.6419576) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.65529167) q[1];
sx q[1];
rz(-2.5319063) q[1];
sx q[1];
rz(-2.8279081) q[1];
x q[2];
rz(2.4645803) q[3];
sx q[3];
rz(-1.5572064) q[3];
sx q[3];
rz(3.0771067) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.3736854) q[2];
sx q[2];
rz(-1.7932971) q[2];
sx q[2];
rz(2.574004) q[2];
rz(-0.23809412) q[3];
sx q[3];
rz(-1.6425491) q[3];
sx q[3];
rz(2.8897918) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.046086144) q[0];
sx q[0];
rz(-0.44281414) q[0];
sx q[0];
rz(-2.1777022) q[0];
rz(1.6148199) q[1];
sx q[1];
rz(-1.9183728) q[1];
sx q[1];
rz(0.33178202) q[1];
rz(-pi/2) q[2];
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
rz(1.6388697) q[2];
sx q[2];
rz(-0.34095665) q[2];
sx q[2];
rz(-2.6457977) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.7501237) q[1];
sx q[1];
rz(-1.5771598) q[1];
sx q[1];
rz(-1.4749737) q[1];
x q[2];
rz(-2.0121215) q[3];
sx q[3];
rz(-1.6989162) q[3];
sx q[3];
rz(2.3263423) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.31384599) q[2];
sx q[2];
rz(-1.670994) q[2];
sx q[2];
rz(-0.016679114) q[2];
rz(1.553675) q[3];
sx q[3];
rz(-1.1976306) q[3];
sx q[3];
rz(-2.8488819) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68506771) q[0];
sx q[0];
rz(-2.2790788) q[0];
sx q[0];
rz(-2.9271794) q[0];
rz(0.21687493) q[1];
sx q[1];
rz(-0.82019359) q[1];
sx q[1];
rz(-1.0925426) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8583243) q[0];
sx q[0];
rz(-0.28139925) q[0];
sx q[0];
rz(-0.94080956) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1591805) q[2];
sx q[2];
rz(-1.6493754) q[2];
sx q[2];
rz(-0.99115288) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.8823575) q[1];
sx q[1];
rz(-2.3016788) q[1];
sx q[1];
rz(-0.0543189) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.66555832) q[3];
sx q[3];
rz(-0.98057038) q[3];
sx q[3];
rz(-1.4233103) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.07) q[2];
sx q[2];
rz(-2.6061974) q[2];
sx q[2];
rz(-0.52616057) q[2];
rz(-1.5633265) q[3];
sx q[3];
rz(-1.1514781) q[3];
sx q[3];
rz(-2.9621942) q[3];
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
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.18463317) q[0];
sx q[0];
rz(-1.3728091) q[0];
sx q[0];
rz(-2.2233295) q[0];
rz(-1.0144764) q[1];
sx q[1];
rz(-2.4160073) q[1];
sx q[1];
rz(-1.135042) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1149108) q[0];
sx q[0];
rz(-1.3068195) q[0];
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
sx q[0];
rz(-pi/2) q[0];
rz(2.8682485) q[1];
sx q[1];
rz(-2.0269718) q[1];
sx q[1];
rz(-0.47246023) q[1];
rz(-2.8718357) q[3];
sx q[3];
rz(-1.7240824) q[3];
sx q[3];
rz(0.80826825) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.6404746) q[2];
sx q[2];
rz(-1.9116118) q[2];
sx q[2];
rz(3.0418923) q[2];
rz(-2.4589608) q[3];
sx q[3];
rz(-1.3411025) q[3];
sx q[3];
rz(-2.7734267) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
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
rz(-2.8354928) q[0];
sx q[0];
rz(-1.2149128) q[0];
sx q[0];
rz(0.039435506) q[0];
rz(-2.4296852) q[1];
sx q[1];
rz(-2.6075677) q[1];
sx q[1];
rz(1.3077516) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.17960571) q[0];
sx q[0];
rz(-1.3276895) q[0];
sx q[0];
rz(-3.0994439) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6823993) q[2];
sx q[2];
rz(-1.204243) q[2];
sx q[2];
rz(-1.9286148) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.2145526) q[1];
sx q[1];
rz(-0.92260375) q[1];
sx q[1];
rz(0.58073466) q[1];
x q[2];
rz(0.51898308) q[3];
sx q[3];
rz(-0.95967889) q[3];
sx q[3];
rz(-2.713969) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.51440614) q[2];
sx q[2];
rz(-1.4428978) q[2];
sx q[2];
rz(2.8241039) q[2];
rz(-0.22335957) q[3];
sx q[3];
rz(-2.375319) q[3];
sx q[3];
rz(1.4661695) q[3];
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
rz(-pi/2) q[0];
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
rz(-1.5668642) q[1];
sx q[1];
rz(-1.3715502) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6487911) q[0];
sx q[0];
rz(-1.6692356) q[0];
sx q[0];
rz(-2.6078694) q[0];
rz(-pi) q[1];
rz(-0.42098896) q[2];
sx q[2];
rz(-1.1596646) q[2];
sx q[2];
rz(-0.56356424) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.4829258) q[1];
sx q[1];
rz(-1.2892168) q[1];
sx q[1];
rz(0.94689178) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3477799) q[3];
sx q[3];
rz(-1.6598637) q[3];
sx q[3];
rz(-1.7410914) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.74662465) q[2];
sx q[2];
rz(-2.0396502) q[2];
sx q[2];
rz(-1.234422) q[2];
rz(1.5359991) q[3];
sx q[3];
rz(-0.77080828) q[3];
sx q[3];
rz(-0.7705645) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4602373) q[0];
sx q[0];
rz(-1.8276031) q[0];
sx q[0];
rz(0.47948691) q[0];
rz(-1.208249) q[1];
sx q[1];
rz(-1.6069326) q[1];
sx q[1];
rz(-1.23034) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5109864) q[0];
sx q[0];
rz(-0.85737757) q[0];
sx q[0];
rz(-3.0632714) q[0];
x q[1];
rz(-2.6012457) q[2];
sx q[2];
rz(-1.6896938) q[2];
sx q[2];
rz(2.8991097) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-1.2895256) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2097384) q[3];
sx q[3];
rz(-2.8936675) q[3];
sx q[3];
rz(2.7558079) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.49518195) q[2];
sx q[2];
rz(-0.16972204) q[2];
sx q[2];
rz(2.1378311) q[2];
rz(1.1209925) q[3];
sx q[3];
rz(-1.5209773) q[3];
sx q[3];
rz(2.3641172) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
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
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6459226) q[0];
sx q[0];
rz(-2.033253) q[0];
sx q[0];
rz(-0.6301737) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1800351) q[2];
sx q[2];
rz(-1.8834049) q[2];
sx q[2];
rz(1.7824204) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.19329883) q[1];
sx q[1];
rz(-1.5806075) q[1];
sx q[1];
rz(1.3693891) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8776953) q[3];
sx q[3];
rz(-0.53987316) q[3];
sx q[3];
rz(-0.25800522) q[3];
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
rz(1.099769) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5438775) q[0];
sx q[0];
rz(-1.3793722) q[0];
sx q[0];
rz(-2.5482225) q[0];
rz(-1.5229535) q[1];
sx q[1];
rz(-0.41888371) q[1];
sx q[1];
rz(0.12817344) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20505894) q[0];
sx q[0];
rz(-1.4292875) q[0];
sx q[0];
rz(2.9756484) q[0];
rz(-pi) q[1];
rz(1.9680357) q[2];
sx q[2];
rz(-2.4284869) q[2];
sx q[2];
rz(1.2867868) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.3230774) q[1];
sx q[1];
rz(-1.0836731) q[1];
sx q[1];
rz(-3.0534153) q[1];
rz(-pi) q[2];
rz(-0.87824741) q[3];
sx q[3];
rz(-1.9345648) q[3];
sx q[3];
rz(-2.6122746) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.0920948) q[2];
sx q[2];
rz(-2.8785661) q[2];
sx q[2];
rz(2.5035739) q[2];
rz(-2.0045896) q[3];
sx q[3];
rz(-2.548389) q[3];
sx q[3];
rz(0.039464522) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9395102) q[0];
sx q[0];
rz(-1.408229) q[0];
sx q[0];
rz(-0.73233488) q[0];
rz(-0.15728514) q[1];
sx q[1];
rz(-1.5443677) q[1];
sx q[1];
rz(-2.0462012) q[1];
rz(-1.0553817) q[2];
sx q[2];
rz(-1.3414186) q[2];
sx q[2];
rz(2.338496) q[2];
rz(0.91315837) q[3];
sx q[3];
rz(-1.7226333) q[3];
sx q[3];
rz(-2.7992579) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
