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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5012623) q[0];
sx q[0];
rz(-0.20168951) q[0];
sx q[0];
rz(1.8224688) q[0];
rz(-pi) q[1];
rz(1.6221912) q[2];
sx q[2];
rz(-1.8784461) q[2];
sx q[2];
rz(0.6531229) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-3.1069965) q[1];
sx q[1];
rz(-1.6507848) q[1];
sx q[1];
rz(-0.75350113) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.66459951) q[3];
sx q[3];
rz(-2.0076111) q[3];
sx q[3];
rz(-2.7561881) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.2166298) q[2];
sx q[2];
rz(-1.7731885) q[2];
sx q[2];
rz(-1.487757) q[2];
rz(-1.1612859) q[3];
sx q[3];
rz(-1.0043283) q[3];
sx q[3];
rz(-2.0465093) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0682812) q[0];
sx q[0];
rz(-0.40638766) q[0];
sx q[0];
rz(0.43346369) q[0];
rz(2.0794226) q[1];
sx q[1];
rz(-1.5824317) q[1];
sx q[1];
rz(2.4210222) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.61536371) q[0];
sx q[0];
rz(-0.13196346) q[0];
sx q[0];
rz(1.2864287) q[0];
rz(-0.24777221) q[2];
sx q[2];
rz(-1.0508014) q[2];
sx q[2];
rz(1.8682075) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.4318254) q[1];
sx q[1];
rz(-2.7569237) q[1];
sx q[1];
rz(-0.40527126) q[1];
rz(-pi) q[2];
rz(-2.0907164) q[3];
sx q[3];
rz(-2.5934016) q[3];
sx q[3];
rz(-2.2238126) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.34708193) q[2];
sx q[2];
rz(-1.9158659) q[2];
sx q[2];
rz(-0.98428717) q[2];
rz(-2.2847564) q[3];
sx q[3];
rz(-2.555116) q[3];
sx q[3];
rz(-1.2543359) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1089351) q[0];
sx q[0];
rz(-0.28626838) q[0];
sx q[0];
rz(-0.74288595) q[0];
rz(3.1318829) q[1];
sx q[1];
rz(-1.566889) q[1];
sx q[1];
rz(-0.28365338) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6742636) q[0];
sx q[0];
rz(-1.0630075) q[0];
sx q[0];
rz(-1.9697777) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4667612) q[2];
sx q[2];
rz(-2.8877604) q[2];
sx q[2];
rz(0.93364894) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.8130506) q[1];
sx q[1];
rz(-2.790741) q[1];
sx q[1];
rz(-2.0787129) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3464706) q[3];
sx q[3];
rz(-3.0375518) q[3];
sx q[3];
rz(-1.4416172) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.11188406) q[2];
sx q[2];
rz(-1.4446898) q[2];
sx q[2];
rz(-2.4857944) q[2];
rz(-0.2119952) q[3];
sx q[3];
rz(-2.5097804) q[3];
sx q[3];
rz(1.1727715) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.028932171) q[0];
sx q[0];
rz(-2.6472968) q[0];
sx q[0];
rz(1.5869045) q[0];
rz(0.2298062) q[1];
sx q[1];
rz(-2.4202012) q[1];
sx q[1];
rz(-1.302964) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36520805) q[0];
sx q[0];
rz(-2.1094649) q[0];
sx q[0];
rz(-2.0950277) q[0];
rz(-0.26891687) q[2];
sx q[2];
rz(-1.285621) q[2];
sx q[2];
rz(1.1095003) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.069217056) q[1];
sx q[1];
rz(-2.0942959) q[1];
sx q[1];
rz(1.545946) q[1];
x q[2];
rz(-2.0992364) q[3];
sx q[3];
rz(-0.63066777) q[3];
sx q[3];
rz(-1.809465) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.0064696781) q[2];
sx q[2];
rz(-0.80265704) q[2];
sx q[2];
rz(2.6378677) q[2];
rz(1.5375562) q[3];
sx q[3];
rz(-1.9683014) q[3];
sx q[3];
rz(1.5639719) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2876494) q[0];
sx q[0];
rz(-2.8048153) q[0];
sx q[0];
rz(0.42959282) q[0];
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
rz(-2.5693656) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2524302) q[2];
sx q[2];
rz(-0.7944383) q[2];
sx q[2];
rz(-0.81046644) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.0068897) q[1];
sx q[1];
rz(-2.0102894) q[1];
sx q[1];
rz(1.2864) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8183299) q[3];
sx q[3];
rz(-1.8283852) q[3];
sx q[3];
rz(-3.0158693) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.2499007) q[2];
sx q[2];
rz(-2.4155858) q[2];
sx q[2];
rz(-0.1740087) q[2];
rz(-1.2627259) q[3];
sx q[3];
rz(-1.5401253) q[3];
sx q[3];
rz(-1.5549829) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
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
rz(-pi) q[1];
x q[1];
rz(-2.2416297) q[2];
sx q[2];
rz(-2.6113178) q[2];
sx q[2];
rz(1.8535623) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.90112061) q[1];
sx q[1];
rz(-1.8163306) q[1];
sx q[1];
rz(-2.4811109) q[1];
rz(-pi) q[2];
rz(0.31474416) q[3];
sx q[3];
rz(-2.5786244) q[3];
sx q[3];
rz(1.1551577) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-3.1098108) q[2];
sx q[2];
rz(-1.6620771) q[2];
sx q[2];
rz(2.7500582) q[2];
rz(1.3854965) q[3];
sx q[3];
rz(-2.9031495) q[3];
sx q[3];
rz(-2.5029206) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2077797) q[0];
sx q[0];
rz(-0.55140984) q[0];
sx q[0];
rz(-1.8361924) q[0];
rz(-0.50453672) q[1];
sx q[1];
rz(-1.4089818) q[1];
sx q[1];
rz(2.9171468) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0248666) q[0];
sx q[0];
rz(-3.118692) q[0];
sx q[0];
rz(2.3842627) q[0];
x q[1];
rz(1.7845909) q[2];
sx q[2];
rz(-1.0875306) q[2];
sx q[2];
rz(0.69141206) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.0266155) q[1];
sx q[1];
rz(-1.8375641) q[1];
sx q[1];
rz(-1.90735) q[1];
rz(1.2013632) q[3];
sx q[3];
rz(-1.8734341) q[3];
sx q[3];
rz(1.7517881) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.9597783) q[2];
sx q[2];
rz(-2.0373127) q[2];
sx q[2];
rz(-2.1941198) q[2];
rz(-0.13253658) q[3];
sx q[3];
rz(-2.4819076) q[3];
sx q[3];
rz(2.1555677) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9082311) q[0];
sx q[0];
rz(-2.7955604) q[0];
sx q[0];
rz(-2.76873) q[0];
rz(-2.8064959) q[1];
sx q[1];
rz(-1.9357888) q[1];
sx q[1];
rz(2.0705409) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.459254) q[0];
sx q[0];
rz(-1.5821348) q[0];
sx q[0];
rz(3.0724694) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.803587) q[2];
sx q[2];
rz(-0.57799423) q[2];
sx q[2];
rz(-1.6737311) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.1879857) q[1];
sx q[1];
rz(-1.3380357) q[1];
sx q[1];
rz(2.8899419) q[1];
rz(1.2921807) q[3];
sx q[3];
rz(-0.33383402) q[3];
sx q[3];
rz(-2.1899471) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.7143453) q[2];
sx q[2];
rz(-1.5421966) q[2];
sx q[2];
rz(1.0672807) q[2];
rz(-2.769477) q[3];
sx q[3];
rz(-0.99388638) q[3];
sx q[3];
rz(-2.9414162) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2315955) q[0];
sx q[0];
rz(-0.59760004) q[0];
sx q[0];
rz(-1.4071314) q[0];
rz(-2.097997) q[1];
sx q[1];
rz(-0.93799543) q[1];
sx q[1];
rz(0.49120206) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5564767) q[0];
sx q[0];
rz(-2.399754) q[0];
sx q[0];
rz(-2.0061352) q[0];
rz(2.9727761) q[2];
sx q[2];
rz(-1.7387317) q[2];
sx q[2];
rz(-0.68532787) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.8206827) q[1];
sx q[1];
rz(-1.9100128) q[1];
sx q[1];
rz(-1.7103819) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3902499) q[3];
sx q[3];
rz(-1.9420764) q[3];
sx q[3];
rz(-2.2425687) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.6970984) q[2];
sx q[2];
rz(-1.5105057) q[2];
sx q[2];
rz(-2.0972882) q[2];
rz(-2.3704884) q[3];
sx q[3];
rz(-1.6090569) q[3];
sx q[3];
rz(-2.7842298) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.20030178) q[0];
sx q[0];
rz(-1.8481978) q[0];
sx q[0];
rz(-1.6218761) q[0];
rz(-1.8408076) q[1];
sx q[1];
rz(-1.3011353) q[1];
sx q[1];
rz(-2.4470952) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25481564) q[0];
sx q[0];
rz(-0.22829994) q[0];
sx q[0];
rz(-1.5103025) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4378261) q[2];
sx q[2];
rz(-2.0914234) q[2];
sx q[2];
rz(0.43308738) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.4150691) q[1];
sx q[1];
rz(-1.5818672) q[1];
sx q[1];
rz(1.9048208) q[1];
rz(-pi) q[2];
rz(0.75218529) q[3];
sx q[3];
rz(-1.2818517) q[3];
sx q[3];
rz(-2.1352701) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.0684315) q[2];
sx q[2];
rz(-1.7794926) q[2];
sx q[2];
rz(2.6279602) q[2];
rz(2.8237776) q[3];
sx q[3];
rz(-1.603926) q[3];
sx q[3];
rz(1.0191127) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1401055) q[0];
sx q[0];
rz(-1.6248063) q[0];
sx q[0];
rz(-0.90909062) q[0];
rz(-1.2263251) q[1];
sx q[1];
rz(-1.6188123) q[1];
sx q[1];
rz(0.34368044) q[1];
rz(-0.60811715) q[2];
sx q[2];
rz(-2.3454852) q[2];
sx q[2];
rz(-0.15783068) q[2];
rz(1.5671586) q[3];
sx q[3];
rz(-2.6156577) q[3];
sx q[3];
rz(-1.2729264) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
