OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.926446) q[0];
sx q[0];
rz(-0.67987052) q[0];
sx q[0];
rz(3.0409066) q[0];
rz(-2.6868532) q[1];
sx q[1];
rz(-2.4603031) q[1];
sx q[1];
rz(1.0256306) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68368841) q[0];
sx q[0];
rz(-1.6207028) q[0];
sx q[0];
rz(1.3752974) q[0];
x q[1];
rz(-2.8335613) q[2];
sx q[2];
rz(-1.6197761) q[2];
sx q[2];
rz(0.9020976) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.6903578) q[1];
sx q[1];
rz(-0.75690311) q[1];
sx q[1];
rz(0.1166269) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0354393) q[3];
sx q[3];
rz(-2.1637754) q[3];
sx q[3];
rz(-1.6361332) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.2166298) q[2];
sx q[2];
rz(-1.7731885) q[2];
sx q[2];
rz(-1.487757) q[2];
rz(-1.9803068) q[3];
sx q[3];
rz(-2.1372644) q[3];
sx q[3];
rz(1.0950834) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0733114) q[0];
sx q[0];
rz(-0.40638766) q[0];
sx q[0];
rz(2.708129) q[0];
rz(-2.0794226) q[1];
sx q[1];
rz(-1.559161) q[1];
sx q[1];
rz(2.4210222) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.90209157) q[0];
sx q[0];
rz(-1.4441617) q[0];
sx q[0];
rz(-0.037221639) q[0];
x q[1];
rz(0.24777221) q[2];
sx q[2];
rz(-1.0508014) q[2];
sx q[2];
rz(1.2733851) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.51746619) q[1];
sx q[1];
rz(-1.4223012) q[1];
sx q[1];
rz(0.35617574) q[1];
rz(-pi) q[2];
rz(2.8470542) q[3];
sx q[3];
rz(-1.1014767) q[3];
sx q[3];
rz(0.32696262) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.7945107) q[2];
sx q[2];
rz(-1.2257267) q[2];
sx q[2];
rz(2.1573055) q[2];
rz(0.85683626) q[3];
sx q[3];
rz(-0.58647668) q[3];
sx q[3];
rz(1.2543359) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1089351) q[0];
sx q[0];
rz(-0.28626838) q[0];
sx q[0];
rz(2.3987067) q[0];
rz(-0.0097097857) q[1];
sx q[1];
rz(-1.566889) q[1];
sx q[1];
rz(2.8579393) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3056639) q[0];
sx q[0];
rz(-1.2244512) q[0];
sx q[0];
rz(-2.5983174) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7314807) q[2];
sx q[2];
rz(-1.3734439) q[2];
sx q[2];
rz(1.6244217) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.4176826) q[1];
sx q[1];
rz(-1.7387448) q[1];
sx q[1];
rz(1.2612832) q[1];
x q[2];
rz(1.795122) q[3];
sx q[3];
rz(-3.0375518) q[3];
sx q[3];
rz(1.6999754) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.11188406) q[2];
sx q[2];
rz(-1.6969029) q[2];
sx q[2];
rz(2.4857944) q[2];
rz(2.9295975) q[3];
sx q[3];
rz(-0.63181221) q[3];
sx q[3];
rz(-1.1727715) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1126605) q[0];
sx q[0];
rz(-0.49429587) q[0];
sx q[0];
rz(-1.5869045) q[0];
rz(-0.2298062) q[1];
sx q[1];
rz(-2.4202012) q[1];
sx q[1];
rz(-1.8386286) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.479871) q[0];
sx q[0];
rz(-0.73307836) q[0];
sx q[0];
rz(-0.69723155) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2755865) q[2];
sx q[2];
rz(-1.3129932) q[2];
sx q[2];
rz(-0.38391963) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.6524383) q[1];
sx q[1];
rz(-1.5492747) q[1];
sx q[1];
rz(-0.52363327) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0082208) q[3];
sx q[3];
rz(-1.8726714) q[3];
sx q[3];
rz(-0.20193298) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.0064696781) q[2];
sx q[2];
rz(-0.80265704) q[2];
sx q[2];
rz(-2.6378677) q[2];
rz(1.6040365) q[3];
sx q[3];
rz(-1.1732912) q[3];
sx q[3];
rz(-1.5776207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.85394323) q[0];
sx q[0];
rz(-2.8048153) q[0];
sx q[0];
rz(-2.7119998) q[0];
rz(0.26643878) q[1];
sx q[1];
rz(-0.8100422) q[1];
sx q[1];
rz(-1.0688759) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7212413) q[0];
sx q[0];
rz(-1.9801894) q[0];
sx q[0];
rz(-2.4024525) q[0];
x q[1];
rz(-0.30854723) q[2];
sx q[2];
rz(-0.82627901) q[2];
sx q[2];
rz(-2.770785) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.0068897) q[1];
sx q[1];
rz(-1.1313032) q[1];
sx q[1];
rz(-1.2864) q[1];
rz(1.3232628) q[3];
sx q[3];
rz(-1.8283852) q[3];
sx q[3];
rz(-3.0158693) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.8916919) q[2];
sx q[2];
rz(-2.4155858) q[2];
sx q[2];
rz(-0.1740087) q[2];
rz(1.8788667) q[3];
sx q[3];
rz(-1.5401253) q[3];
sx q[3];
rz(1.5866097) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0083017666) q[0];
sx q[0];
rz(-2.1228078) q[0];
sx q[0];
rz(0.89023501) q[0];
rz(1.4324191) q[1];
sx q[1];
rz(-2.1232257) q[1];
sx q[1];
rz(-0.11996809) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4910866) q[0];
sx q[0];
rz(-2.1236558) q[0];
sx q[0];
rz(-1.4006459) q[0];
rz(-pi) q[1];
rz(2.7920958) q[2];
sx q[2];
rz(-1.9781402) q[2];
sx q[2];
rz(2.0317589) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.240472) q[1];
sx q[1];
rz(-1.8163306) q[1];
sx q[1];
rz(-0.66048173) q[1];
rz(2.8268485) q[3];
sx q[3];
rz(-0.56296825) q[3];
sx q[3];
rz(1.1551577) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.031781901) q[2];
sx q[2];
rz(-1.6620771) q[2];
sx q[2];
rz(0.39153448) q[2];
rz(-1.7560962) q[3];
sx q[3];
rz(-2.9031495) q[3];
sx q[3];
rz(-2.5029206) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
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
rz(2.2077797) q[0];
sx q[0];
rz(-0.55140984) q[0];
sx q[0];
rz(-1.3054003) q[0];
rz(0.50453672) q[1];
sx q[1];
rz(-1.4089818) q[1];
sx q[1];
rz(0.22444589) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.87418694) q[0];
sx q[0];
rz(-1.5874369) q[0];
sx q[0];
rz(1.5865302) q[0];
rz(-pi) q[1];
rz(0.4928035) q[2];
sx q[2];
rz(-1.7597919) q[2];
sx q[2];
rz(2.1616621) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.0266155) q[1];
sx q[1];
rz(-1.3040286) q[1];
sx q[1];
rz(-1.2342427) q[1];
rz(-0.85783463) q[3];
sx q[3];
rz(-0.47311764) q[3];
sx q[3];
rz(0.83698646) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.18181431) q[2];
sx q[2];
rz(-2.0373127) q[2];
sx q[2];
rz(0.9474729) q[2];
rz(3.0090561) q[3];
sx q[3];
rz(-2.4819076) q[3];
sx q[3];
rz(2.1555677) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2333616) q[0];
sx q[0];
rz(-0.34603226) q[0];
sx q[0];
rz(0.37286266) q[0];
rz(-2.8064959) q[1];
sx q[1];
rz(-1.9357888) q[1];
sx q[1];
rz(2.0705409) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.11075739) q[0];
sx q[0];
rz(-1.5016775) q[0];
sx q[0];
rz(1.5594307) q[0];
rz(-0.33800563) q[2];
sx q[2];
rz(-0.57799423) q[2];
sx q[2];
rz(1.6737311) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.95360699) q[1];
sx q[1];
rz(-1.3380357) q[1];
sx q[1];
rz(-0.25165073) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8926431) q[3];
sx q[3];
rz(-1.6610356) q[3];
sx q[3];
rz(-2.7863996) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.7143453) q[2];
sx q[2];
rz(-1.5421966) q[2];
sx q[2];
rz(2.074312) q[2];
rz(-0.37211564) q[3];
sx q[3];
rz(-0.99388638) q[3];
sx q[3];
rz(2.9414162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9099971) q[0];
sx q[0];
rz(-0.59760004) q[0];
sx q[0];
rz(1.4071314) q[0];
rz(-2.097997) q[1];
sx q[1];
rz(-0.93799543) q[1];
sx q[1];
rz(0.49120206) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8255913) q[0];
sx q[0];
rz(-1.859731) q[0];
sx q[0];
rz(2.2641473) q[0];
rz(-pi) q[1];
rz(-2.3517026) q[2];
sx q[2];
rz(-0.23755506) q[2];
sx q[2];
rz(-1.4804763) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.2031695) q[1];
sx q[1];
rz(-1.7023801) q[1];
sx q[1];
rz(-2.7992976) q[1];
x q[2];
rz(-0.37685412) q[3];
sx q[3];
rz(-1.7389193) q[3];
sx q[3];
rz(0.73790077) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.4444943) q[2];
sx q[2];
rz(-1.631087) q[2];
sx q[2];
rz(-2.0972882) q[2];
rz(2.3704884) q[3];
sx q[3];
rz(-1.6090569) q[3];
sx q[3];
rz(-0.3573629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20030178) q[0];
sx q[0];
rz(-1.8481978) q[0];
sx q[0];
rz(-1.5197165) q[0];
rz(-1.8408076) q[1];
sx q[1];
rz(-1.3011353) q[1];
sx q[1];
rz(-2.4470952) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2570528) q[0];
sx q[0];
rz(-1.5571132) q[0];
sx q[0];
rz(1.3428997) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6171309) q[2];
sx q[2];
rz(-1.6860644) q[2];
sx q[2];
rz(-1.0712717) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.81238231) q[1];
sx q[1];
rz(-2.8073916) q[1];
sx q[1];
rz(-1.5370395) q[1];
x q[2];
rz(1.9574088) q[3];
sx q[3];
rz(-0.85672934) q[3];
sx q[3];
rz(-2.3165645) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.0731611) q[2];
sx q[2];
rz(-1.3621) q[2];
sx q[2];
rz(-2.6279602) q[2];
rz(-2.8237776) q[3];
sx q[3];
rz(-1.603926) q[3];
sx q[3];
rz(2.1224799) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1401055) q[0];
sx q[0];
rz(-1.6248063) q[0];
sx q[0];
rz(-0.90909062) q[0];
rz(-1.9152676) q[1];
sx q[1];
rz(-1.5227804) q[1];
sx q[1];
rz(-2.7979122) q[1];
rz(-2.0991391) q[2];
sx q[2];
rz(-0.94403841) q[2];
sx q[2];
rz(-0.94081139) q[2];
rz(1.0448643) q[3];
sx q[3];
rz(-1.5726225) q[3];
sx q[3];
rz(-2.8468688) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
