OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.2151467) q[0];
sx q[0];
rz(0.67987052) q[0];
sx q[0];
rz(9.3240919) q[0];
rz(-2.6868532) q[1];
sx q[1];
rz(-2.4603031) q[1];
sx q[1];
rz(1.0256306) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2446063) q[0];
sx q[0];
rz(-1.7660487) q[0];
sx q[0];
rz(-3.0907187) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5194015) q[2];
sx q[2];
rz(-1.8784461) q[2];
sx q[2];
rz(0.6531229) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.4512349) q[1];
sx q[1];
rz(-0.75690311) q[1];
sx q[1];
rz(3.0249658) q[1];
rz(-pi) q[2];
rz(-2.1061534) q[3];
sx q[3];
rz(-2.1637754) q[3];
sx q[3];
rz(1.5054594) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.2166298) q[2];
sx q[2];
rz(-1.7731885) q[2];
sx q[2];
rz(1.6538357) q[2];
rz(-1.9803068) q[3];
sx q[3];
rz(-1.0043283) q[3];
sx q[3];
rz(-1.0950834) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0733114) q[0];
sx q[0];
rz(-2.735205) q[0];
sx q[0];
rz(0.43346369) q[0];
rz(1.06217) q[1];
sx q[1];
rz(-1.559161) q[1];
sx q[1];
rz(2.4210222) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.61536371) q[0];
sx q[0];
rz(-3.0096292) q[0];
sx q[0];
rz(1.855164) q[0];
rz(-pi) q[1];
rz(-0.24777221) q[2];
sx q[2];
rz(-2.0907913) q[2];
sx q[2];
rz(1.2733851) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.4318254) q[1];
sx q[1];
rz(-2.7569237) q[1];
sx q[1];
rz(-2.7363214) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0581011) q[3];
sx q[3];
rz(-1.8326899) q[3];
sx q[3];
rz(2.0341121) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.34708193) q[2];
sx q[2];
rz(-1.9158659) q[2];
sx q[2];
rz(2.1573055) q[2];
rz(0.85683626) q[3];
sx q[3];
rz(-0.58647668) q[3];
sx q[3];
rz(-1.8872567) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1089351) q[0];
sx q[0];
rz(-2.8553243) q[0];
sx q[0];
rz(2.3987067) q[0];
rz(-0.0097097857) q[1];
sx q[1];
rz(-1.5747036) q[1];
sx q[1];
rz(0.28365338) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8946202) q[0];
sx q[0];
rz(-0.63480154) q[0];
sx q[0];
rz(2.532124) q[0];
rz(0.67483141) q[2];
sx q[2];
rz(-0.25383224) q[2];
sx q[2];
rz(-0.93364894) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.3481118) q[1];
sx q[1];
rz(-1.2657796) q[1];
sx q[1];
rz(-2.9654337) q[1];
rz(1.4693442) q[3];
sx q[3];
rz(-1.5939004) q[3];
sx q[3];
rz(2.7892609) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-3.0297086) q[2];
sx q[2];
rz(-1.4446898) q[2];
sx q[2];
rz(2.4857944) q[2];
rz(0.2119952) q[3];
sx q[3];
rz(-0.63181221) q[3];
sx q[3];
rz(1.1727715) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.028932171) q[0];
sx q[0];
rz(-2.6472968) q[0];
sx q[0];
rz(1.5869045) q[0];
rz(-2.9117865) q[1];
sx q[1];
rz(-0.72139144) q[1];
sx q[1];
rz(1.302964) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.479871) q[0];
sx q[0];
rz(-2.4085143) q[0];
sx q[0];
rz(-0.69723155) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8660061) q[2];
sx q[2];
rz(-1.8285995) q[2];
sx q[2];
rz(-0.38391963) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.11889549) q[1];
sx q[1];
rz(-2.6175584) q[1];
sx q[1];
rz(0.043020821) q[1];
rz(0.35273055) q[3];
sx q[3];
rz(-2.1051072) q[3];
sx q[3];
rz(-1.9580458) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-3.135123) q[2];
sx q[2];
rz(-0.80265704) q[2];
sx q[2];
rz(-0.50372493) q[2];
rz(1.6040365) q[3];
sx q[3];
rz(-1.1732912) q[3];
sx q[3];
rz(1.5639719) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2876494) q[0];
sx q[0];
rz(-0.33677736) q[0];
sx q[0];
rz(-2.7119998) q[0];
rz(-0.26643878) q[1];
sx q[1];
rz(-0.8100422) q[1];
sx q[1];
rz(1.0688759) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5621786) q[0];
sx q[0];
rz(-0.82579188) q[0];
sx q[0];
rz(0.572227) q[0];
rz(1.2524302) q[2];
sx q[2];
rz(-0.7944383) q[2];
sx q[2];
rz(-0.81046644) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.0068897) q[1];
sx q[1];
rz(-2.0102894) q[1];
sx q[1];
rz(1.8551926) q[1];
x q[2];
rz(-2.8762749) q[3];
sx q[3];
rz(-1.331592) q[3];
sx q[3];
rz(-1.632229) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.2499007) q[2];
sx q[2];
rz(-0.72600681) q[2];
sx q[2];
rz(2.967584) q[2];
rz(1.2627259) q[3];
sx q[3];
rz(-1.6014674) q[3];
sx q[3];
rz(-1.5549829) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
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
rz(-3.1332909) q[0];
sx q[0];
rz(-2.1228078) q[0];
sx q[0];
rz(-2.2513576) q[0];
rz(-1.7091735) q[1];
sx q[1];
rz(-1.018367) q[1];
sx q[1];
rz(-3.0216246) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9719043) q[0];
sx q[0];
rz(-1.426188) q[0];
sx q[0];
rz(2.5822116) q[0];
x q[1];
rz(-1.1402848) q[2];
sx q[2];
rz(-1.2509648) q[2];
sx q[2];
rz(-0.31756535) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.4830345) q[1];
sx q[1];
rz(-2.2081828) q[1];
sx q[1];
rz(-1.8780707) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.54048583) q[3];
sx q[3];
rz(-1.7367762) q[3];
sx q[3];
rz(-2.9946208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.031781901) q[2];
sx q[2];
rz(-1.4795156) q[2];
sx q[2];
rz(0.39153448) q[2];
rz(-1.3854965) q[3];
sx q[3];
rz(-2.9031495) q[3];
sx q[3];
rz(2.5029206) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2077797) q[0];
sx q[0];
rz(-2.5901828) q[0];
sx q[0];
rz(-1.3054003) q[0];
rz(0.50453672) q[1];
sx q[1];
rz(-1.7326109) q[1];
sx q[1];
rz(-0.22444589) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.11672606) q[0];
sx q[0];
rz(-0.022900669) q[0];
sx q[0];
rz(-0.75732996) q[0];
rz(-pi) q[1];
rz(0.38421696) q[2];
sx q[2];
rz(-0.52500341) q[2];
sx q[2];
rz(2.8872761) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.689381) q[1];
sx q[1];
rz(-1.895001) q[1];
sx q[1];
rz(0.28181847) q[1];
rz(2.283758) q[3];
sx q[3];
rz(-0.47311764) q[3];
sx q[3];
rz(-2.3046062) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.9597783) q[2];
sx q[2];
rz(-2.0373127) q[2];
sx q[2];
rz(-0.9474729) q[2];
rz(-0.13253658) q[3];
sx q[3];
rz(-2.4819076) q[3];
sx q[3];
rz(-0.98602492) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2333616) q[0];
sx q[0];
rz(-0.34603226) q[0];
sx q[0];
rz(-0.37286266) q[0];
rz(2.8064959) q[1];
sx q[1];
rz(-1.9357888) q[1];
sx q[1];
rz(-2.0705409) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.11075739) q[0];
sx q[0];
rz(-1.6399151) q[0];
sx q[0];
rz(1.5821619) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.33800563) q[2];
sx q[2];
rz(-2.5635984) q[2];
sx q[2];
rz(1.4678616) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.1879857) q[1];
sx q[1];
rz(-1.803557) q[1];
sx q[1];
rz(2.8899419) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.849412) q[3];
sx q[3];
rz(-0.33383402) q[3];
sx q[3];
rz(0.95164559) q[3];
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
rz(-1.0672807) q[2];
rz(0.37211564) q[3];
sx q[3];
rz(-2.1477063) q[3];
sx q[3];
rz(2.9414162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2315955) q[0];
sx q[0];
rz(-2.5439926) q[0];
sx q[0];
rz(1.7344612) q[0];
rz(1.0435957) q[1];
sx q[1];
rz(-0.93799543) q[1];
sx q[1];
rz(0.49120206) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.585116) q[0];
sx q[0];
rz(-0.74183861) q[0];
sx q[0];
rz(-1.1354574) q[0];
rz(0.16881659) q[2];
sx q[2];
rz(-1.4028609) q[2];
sx q[2];
rz(-0.68532787) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.8206827) q[1];
sx q[1];
rz(-1.9100128) q[1];
sx q[1];
rz(-1.7103819) q[1];
rz(-pi) q[2];
x q[2];
rz(0.43213958) q[3];
sx q[3];
rz(-2.7305805) q[3];
sx q[3];
rz(0.43288818) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.6970984) q[2];
sx q[2];
rz(-1.631087) q[2];
sx q[2];
rz(2.0972882) q[2];
rz(-2.3704884) q[3];
sx q[3];
rz(-1.5325357) q[3];
sx q[3];
rz(2.7842298) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.20030178) q[0];
sx q[0];
rz(-1.2933949) q[0];
sx q[0];
rz(1.6218761) q[0];
rz(1.3007851) q[1];
sx q[1];
rz(-1.3011353) q[1];
sx q[1];
rz(0.69449743) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8845399) q[0];
sx q[0];
rz(-1.5844795) q[0];
sx q[0];
rz(1.798693) q[0];
rz(1.7037665) q[2];
sx q[2];
rz(-1.0501692) q[2];
sx q[2];
rz(-0.43308738) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.3292103) q[1];
sx q[1];
rz(-0.33420104) q[1];
sx q[1];
rz(-1.6045531) q[1];
rz(-pi) q[2];
rz(-1.1841838) q[3];
sx q[3];
rz(-0.85672934) q[3];
sx q[3];
rz(-2.3165645) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.0684315) q[2];
sx q[2];
rz(-1.7794926) q[2];
sx q[2];
rz(2.6279602) q[2];
rz(0.3178151) q[3];
sx q[3];
rz(-1.5376667) q[3];
sx q[3];
rz(1.0191127) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0014872) q[0];
sx q[0];
rz(-1.6248063) q[0];
sx q[0];
rz(-0.90909062) q[0];
rz(1.2263251) q[1];
sx q[1];
rz(-1.5227804) q[1];
sx q[1];
rz(-2.7979122) q[1];
rz(-2.5334755) q[2];
sx q[2];
rz(-0.79610745) q[2];
sx q[2];
rz(2.983762) q[2];
rz(-2.0967284) q[3];
sx q[3];
rz(-1.5726225) q[3];
sx q[3];
rz(-2.8468688) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
