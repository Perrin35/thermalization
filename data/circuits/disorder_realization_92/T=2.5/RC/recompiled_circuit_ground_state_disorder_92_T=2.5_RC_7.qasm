OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.2607245) q[0];
sx q[0];
rz(-0.27624929) q[0];
sx q[0];
rz(-2.2978388) q[0];
rz(1.002797) q[1];
sx q[1];
rz(-2.306814) q[1];
sx q[1];
rz(0.53258449) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7850936) q[0];
sx q[0];
rz(-2.4716232) q[0];
sx q[0];
rz(1.3290624) q[0];
x q[1];
rz(1.7108828) q[2];
sx q[2];
rz(-1.0647961) q[2];
sx q[2];
rz(-1.9078474) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.0610032) q[1];
sx q[1];
rz(-0.71032754) q[1];
sx q[1];
rz(1.3377473) q[1];
x q[2];
rz(0.79382503) q[3];
sx q[3];
rz(-0.50758368) q[3];
sx q[3];
rz(0.48963293) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.8325995) q[2];
sx q[2];
rz(-1.7367881) q[2];
sx q[2];
rz(-3.0187606) q[2];
rz(0.041736688) q[3];
sx q[3];
rz(-1.583497) q[3];
sx q[3];
rz(0.48833716) q[3];
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
rz(-1.1500583) q[0];
sx q[0];
rz(-2.8128862) q[0];
sx q[0];
rz(2.026189) q[0];
rz(2.5967122) q[1];
sx q[1];
rz(-1.3976401) q[1];
sx q[1];
rz(0.56745183) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.7780964) q[0];
sx q[0];
rz(-1.6045185) q[0];
sx q[0];
rz(-3.0959227) q[0];
rz(-pi) q[1];
x q[1];
rz(3.1092875) q[2];
sx q[2];
rz(-0.77649335) q[2];
sx q[2];
rz(0.99054289) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.2360252) q[1];
sx q[1];
rz(-1.4705338) q[1];
sx q[1];
rz(2.5699993) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6065234) q[3];
sx q[3];
rz(-2.5949083) q[3];
sx q[3];
rz(-3.1233623) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.51022092) q[2];
sx q[2];
rz(-3.0192182) q[2];
sx q[2];
rz(-3.0369634) q[2];
rz(-0.50466022) q[3];
sx q[3];
rz(-1.3480836) q[3];
sx q[3];
rz(0.73205718) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
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
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1402682) q[0];
sx q[0];
rz(-2.9502385) q[0];
sx q[0];
rz(1.485317) q[0];
rz(-2.5775919) q[1];
sx q[1];
rz(-1.0537078) q[1];
sx q[1];
rz(-0.82411134) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6862515) q[0];
sx q[0];
rz(-1.4061965) q[0];
sx q[0];
rz(1.0353987) q[0];
rz(-2.5512087) q[2];
sx q[2];
rz(-2.2548642) q[2];
sx q[2];
rz(-1.4223398) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.4802758) q[1];
sx q[1];
rz(-1.7719898) q[1];
sx q[1];
rz(0.43432216) q[1];
x q[2];
rz(1.3585659) q[3];
sx q[3];
rz(-1.5023059) q[3];
sx q[3];
rz(0.93671441) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.092502681) q[2];
sx q[2];
rz(-2.8503032) q[2];
sx q[2];
rz(1.606344) q[2];
rz(-2.2384079) q[3];
sx q[3];
rz(-1.817037) q[3];
sx q[3];
rz(0.42455348) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.20151888) q[0];
sx q[0];
rz(-0.98648447) q[0];
sx q[0];
rz(2.7281813) q[0];
rz(2.9460733) q[1];
sx q[1];
rz(-2.1045411) q[1];
sx q[1];
rz(1.4541385) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.80418865) q[0];
sx q[0];
rz(-1.4436663) q[0];
sx q[0];
rz(-0.8781627) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7614582) q[2];
sx q[2];
rz(-1.9205928) q[2];
sx q[2];
rz(-1.4189394) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.6100055) q[1];
sx q[1];
rz(-1.0949228) q[1];
sx q[1];
rz(-2.8240318) q[1];
rz(-pi) q[2];
rz(2.4471483) q[3];
sx q[3];
rz(-1.4141091) q[3];
sx q[3];
rz(-0.82686916) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.89895189) q[2];
sx q[2];
rz(-0.61598888) q[2];
sx q[2];
rz(-1.8757437) q[2];
rz(0.85593456) q[3];
sx q[3];
rz(-1.7549691) q[3];
sx q[3];
rz(-1.9267193) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.6619381) q[0];
sx q[0];
rz(-2.9582773) q[0];
sx q[0];
rz(-1.4187752) q[0];
rz(-1.4310369) q[1];
sx q[1];
rz(-1.6173671) q[1];
sx q[1];
rz(0.72351825) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.60155) q[0];
sx q[0];
rz(-1.7474084) q[0];
sx q[0];
rz(0.58797449) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3695643) q[2];
sx q[2];
rz(-0.92185045) q[2];
sx q[2];
rz(-1.8963501) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.6910839) q[1];
sx q[1];
rz(-1.5683859) q[1];
sx q[1];
rz(2.2775713) q[1];
rz(-1.5635023) q[3];
sx q[3];
rz(-1.2185343) q[3];
sx q[3];
rz(1.9520813) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.4414703) q[2];
sx q[2];
rz(-1.0157601) q[2];
sx q[2];
rz(0.36386841) q[2];
rz(-1.3077334) q[3];
sx q[3];
rz(-1.4738844) q[3];
sx q[3];
rz(1.3522805) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42809197) q[0];
sx q[0];
rz(-2.2245753) q[0];
sx q[0];
rz(-2.2416903) q[0];
rz(-1.046754) q[1];
sx q[1];
rz(-2.7622107) q[1];
sx q[1];
rz(-2.5894763) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13487694) q[0];
sx q[0];
rz(-1.9089886) q[0];
sx q[0];
rz(-1.5408367) q[0];
x q[1];
rz(-2.5234473) q[2];
sx q[2];
rz(-0.80354881) q[2];
sx q[2];
rz(1.8185563) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.1496755) q[1];
sx q[1];
rz(-2.9475355) q[1];
sx q[1];
rz(2.0582576) q[1];
x q[2];
rz(-0.12142809) q[3];
sx q[3];
rz(-2.5478706) q[3];
sx q[3];
rz(-0.93103257) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.8206574) q[2];
sx q[2];
rz(-0.67639095) q[2];
sx q[2];
rz(-2.9429842) q[2];
rz(2.0123539) q[3];
sx q[3];
rz(-1.6695453) q[3];
sx q[3];
rz(2.3614597) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.3405014) q[0];
sx q[0];
rz(-1.3022364) q[0];
sx q[0];
rz(2.4275725) q[0];
rz(-1.9188312) q[1];
sx q[1];
rz(-2.265265) q[1];
sx q[1];
rz(-2.8661935) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31322436) q[0];
sx q[0];
rz(-1.3592915) q[0];
sx q[0];
rz(-1.1660035) q[0];
x q[1];
rz(0.95850079) q[2];
sx q[2];
rz(-1.3103518) q[2];
sx q[2];
rz(2.1672003) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(3.1186953) q[1];
sx q[1];
rz(-2.0400908) q[1];
sx q[1];
rz(-2.1170627) q[1];
rz(-pi) q[2];
rz(-2.4163298) q[3];
sx q[3];
rz(-1.5808839) q[3];
sx q[3];
rz(2.6731104) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.3068646) q[2];
sx q[2];
rz(-1.4707668) q[2];
sx q[2];
rz(1.4677706) q[2];
rz(1.62014) q[3];
sx q[3];
rz(-2.455267) q[3];
sx q[3];
rz(2.2036688) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15243212) q[0];
sx q[0];
rz(-0.91355649) q[0];
sx q[0];
rz(-2.1754225) q[0];
rz(0.78978157) q[1];
sx q[1];
rz(-0.25895324) q[1];
sx q[1];
rz(-2.1678179) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1137266) q[0];
sx q[0];
rz(-1.4513512) q[0];
sx q[0];
rz(1.3238504) q[0];
rz(-pi) q[1];
rz(2.439271) q[2];
sx q[2];
rz(-2.6036545) q[2];
sx q[2];
rz(2.5449716) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.57339797) q[1];
sx q[1];
rz(-2.1713082) q[1];
sx q[1];
rz(-0.37094122) q[1];
rz(-0.059938537) q[3];
sx q[3];
rz(-2.0459243) q[3];
sx q[3];
rz(-2.5157473) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.3070273) q[2];
sx q[2];
rz(-1.4296738) q[2];
sx q[2];
rz(-0.60565051) q[2];
rz(-1.0904795) q[3];
sx q[3];
rz(-2.8993789) q[3];
sx q[3];
rz(3.0428913) q[3];
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
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7379446) q[0];
sx q[0];
rz(-0.048796766) q[0];
sx q[0];
rz(2.5700289) q[0];
rz(2.7538708) q[1];
sx q[1];
rz(-2.3523836) q[1];
sx q[1];
rz(-2.2472084) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.69686517) q[0];
sx q[0];
rz(-2.055671) q[0];
sx q[0];
rz(-2.1567325) q[0];
x q[1];
rz(3.0249075) q[2];
sx q[2];
rz(-1.0067847) q[2];
sx q[2];
rz(2.2224768) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.5928985) q[1];
sx q[1];
rz(-1.1709524) q[1];
sx q[1];
rz(0.26229761) q[1];
rz(-pi) q[2];
rz(-3.0017051) q[3];
sx q[3];
rz(-2.6287946) q[3];
sx q[3];
rz(0.96641216) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.9992708) q[2];
sx q[2];
rz(-2.2597376) q[2];
sx q[2];
rz(2.2692661) q[2];
rz(-3.0669323) q[3];
sx q[3];
rz(-1.9118237) q[3];
sx q[3];
rz(2.5653896) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9999303) q[0];
sx q[0];
rz(-0.44121989) q[0];
sx q[0];
rz(0.73750752) q[0];
rz(-0.69955889) q[1];
sx q[1];
rz(-2.12205) q[1];
sx q[1];
rz(2.2950744) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2071211) q[0];
sx q[0];
rz(-1.1729329) q[0];
sx q[0];
rz(-2.4475696) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4581095) q[2];
sx q[2];
rz(-0.42851617) q[2];
sx q[2];
rz(0.90135114) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.9690035) q[1];
sx q[1];
rz(-2.6903213) q[1];
sx q[1];
rz(-1.8901214) q[1];
rz(-pi) q[2];
rz(-2.722214) q[3];
sx q[3];
rz(-0.43361317) q[3];
sx q[3];
rz(0.24672844) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.6503341) q[2];
sx q[2];
rz(-1.9776055) q[2];
sx q[2];
rz(-2.7872046) q[2];
rz(-0.77704159) q[3];
sx q[3];
rz(-2.5208426) q[3];
sx q[3];
rz(-2.0508155) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1914094) q[0];
sx q[0];
rz(-1.8671028) q[0];
sx q[0];
rz(-2.6381459) q[0];
rz(1.3632111) q[1];
sx q[1];
rz(-2.6111205) q[1];
sx q[1];
rz(3.0725239) q[1];
rz(-1.6949881) q[2];
sx q[2];
rz(-0.71659452) q[2];
sx q[2];
rz(2.334395) q[2];
rz(-2.7991838) q[3];
sx q[3];
rz(-2.0031606) q[3];
sx q[3];
rz(0.99872086) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
