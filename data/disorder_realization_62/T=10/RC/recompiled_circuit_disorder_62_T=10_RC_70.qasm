OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.2405038) q[0];
sx q[0];
rz(3.3189964) q[0];
sx q[0];
rz(11.431974) q[0];
rz(-1.9534684) q[1];
sx q[1];
rz(-1.0367353) q[1];
sx q[1];
rz(0.66361767) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.21270574) q[0];
sx q[0];
rz(-2.0318673) q[0];
sx q[0];
rz(1.5443718) q[0];
rz(-pi) q[1];
rz(2.3637949) q[2];
sx q[2];
rz(-1.3133089) q[2];
sx q[2];
rz(2.4553026) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.8780958) q[1];
sx q[1];
rz(-2.7090008) q[1];
sx q[1];
rz(-2.2357975) q[1];
rz(-1.3271684) q[3];
sx q[3];
rz(-1.4654136) q[3];
sx q[3];
rz(1.5308612) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.87876451) q[2];
sx q[2];
rz(-0.44439134) q[2];
sx q[2];
rz(3.0905241) q[2];
rz(-0.55705327) q[3];
sx q[3];
rz(-2.3414108) q[3];
sx q[3];
rz(-1.5867656) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.58650815) q[0];
sx q[0];
rz(-0.83207911) q[0];
sx q[0];
rz(0.59659514) q[0];
rz(-0.82582981) q[1];
sx q[1];
rz(-1.4412216) q[1];
sx q[1];
rz(-1.2260431) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2873043) q[0];
sx q[0];
rz(-1.1278296) q[0];
sx q[0];
rz(1.3501549) q[0];
x q[1];
rz(2.3833582) q[2];
sx q[2];
rz(-0.97823921) q[2];
sx q[2];
rz(-2.638608) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.5495758) q[1];
sx q[1];
rz(-1.5703652) q[1];
sx q[1];
rz(-1.3509343) q[1];
rz(-1.5854884) q[3];
sx q[3];
rz(-2.5689295) q[3];
sx q[3];
rz(3.1327914) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.79919672) q[2];
sx q[2];
rz(-1.9689955) q[2];
sx q[2];
rz(-0.33102316) q[2];
rz(0.80667574) q[3];
sx q[3];
rz(-1.2253864) q[3];
sx q[3];
rz(1.4276918) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1598635) q[0];
sx q[0];
rz(-1.230343) q[0];
sx q[0];
rz(-1.8925517) q[0];
rz(3.0535835) q[1];
sx q[1];
rz(-1.1227337) q[1];
sx q[1];
rz(-2.1121315) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.35917657) q[0];
sx q[0];
rz(-2.4826907) q[0];
sx q[0];
rz(-1.2078148) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7715893) q[2];
sx q[2];
rz(-0.68379935) q[2];
sx q[2];
rz(-1.5918658) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.0440812) q[1];
sx q[1];
rz(-1.71002) q[1];
sx q[1];
rz(-0.15686762) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0089278) q[3];
sx q[3];
rz(-0.74436114) q[3];
sx q[3];
rz(0.45385195) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.007894667) q[2];
sx q[2];
rz(-1.7241314) q[2];
sx q[2];
rz(-0.50764817) q[2];
rz(1.7525904) q[3];
sx q[3];
rz(-2.8116083) q[3];
sx q[3];
rz(-1.1631789) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65748173) q[0];
sx q[0];
rz(-2.8634475) q[0];
sx q[0];
rz(1.5959651) q[0];
rz(1.0428628) q[1];
sx q[1];
rz(-1.9680126) q[1];
sx q[1];
rz(-1.625659) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0632616) q[0];
sx q[0];
rz(-2.445979) q[0];
sx q[0];
rz(2.9177833) q[0];
x q[1];
rz(1.707294) q[2];
sx q[2];
rz(-2.578306) q[2];
sx q[2];
rz(2.0493281) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.792946) q[1];
sx q[1];
rz(-0.9325087) q[1];
sx q[1];
rz(1.9116886) q[1];
rz(-pi) q[2];
rz(-1.3514148) q[3];
sx q[3];
rz(-0.99321584) q[3];
sx q[3];
rz(2.1882309) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.3778014) q[2];
sx q[2];
rz(-0.8447454) q[2];
sx q[2];
rz(-1.2949004) q[2];
rz(-3.0002248) q[3];
sx q[3];
rz(-0.54261345) q[3];
sx q[3];
rz(-0.98658371) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7871053) q[0];
sx q[0];
rz(-1.1802477) q[0];
sx q[0];
rz(1.6217344) q[0];
rz(-2.5095818) q[1];
sx q[1];
rz(-0.72223392) q[1];
sx q[1];
rz(-2.246726) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8653523) q[0];
sx q[0];
rz(-1.1603174) q[0];
sx q[0];
rz(0.68666896) q[0];
rz(-pi) q[1];
x q[1];
rz(0.72563719) q[2];
sx q[2];
rz(-1.1302395) q[2];
sx q[2];
rz(-2.5051136) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.289031) q[1];
sx q[1];
rz(-0.82419306) q[1];
sx q[1];
rz(-2.186071) q[1];
rz(-pi) q[2];
x q[2];
rz(3.1154247) q[3];
sx q[3];
rz(-0.28655616) q[3];
sx q[3];
rz(-2.9122796) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.52508369) q[2];
sx q[2];
rz(-0.36642763) q[2];
sx q[2];
rz(-1.9990702) q[2];
rz(-0.74674314) q[3];
sx q[3];
rz(-1.8871504) q[3];
sx q[3];
rz(1.07553) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.58364761) q[0];
sx q[0];
rz(-2.8201411) q[0];
sx q[0];
rz(-1.7640132) q[0];
rz(0.47239834) q[1];
sx q[1];
rz(-0.51858416) q[1];
sx q[1];
rz(2.6766052) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2705921) q[0];
sx q[0];
rz(-2.5809079) q[0];
sx q[0];
rz(0.90765783) q[0];
x q[1];
rz(-2.9372413) q[2];
sx q[2];
rz(-2.7933279) q[2];
sx q[2];
rz(1.5262926) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.7558414) q[1];
sx q[1];
rz(-0.63703905) q[1];
sx q[1];
rz(-1.4904651) q[1];
x q[2];
rz(2.432514) q[3];
sx q[3];
rz(-1.7311829) q[3];
sx q[3];
rz(-1.8335613) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.49089367) q[2];
sx q[2];
rz(-1.2630254) q[2];
sx q[2];
rz(-2.6965551) q[2];
rz(0.93368357) q[3];
sx q[3];
rz(-1.4327587) q[3];
sx q[3];
rz(-0.26708189) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.12748195) q[0];
sx q[0];
rz(-1.5662136) q[0];
sx q[0];
rz(-1.7215464) q[0];
rz(0.02380112) q[1];
sx q[1];
rz(-0.61444608) q[1];
sx q[1];
rz(-0.15596095) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6316846) q[0];
sx q[0];
rz(-1.317306) q[0];
sx q[0];
rz(-1.7476728) q[0];
rz(-pi) q[1];
rz(-0.67955741) q[2];
sx q[2];
rz(-2.5153) q[2];
sx q[2];
rz(1.9192413) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.8487726) q[1];
sx q[1];
rz(-0.18310586) q[1];
sx q[1];
rz(-1.8715026) q[1];
rz(-pi) q[2];
rz(-2.6050623) q[3];
sx q[3];
rz(-2.0091669) q[3];
sx q[3];
rz(0.28552548) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.069313958) q[2];
sx q[2];
rz(-2.5645655) q[2];
sx q[2];
rz(2.0689266) q[2];
rz(0.32564751) q[3];
sx q[3];
rz(-1.1762534) q[3];
sx q[3];
rz(-1.6252888) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1919365) q[0];
sx q[0];
rz(-2.7392116) q[0];
sx q[0];
rz(0.33777133) q[0];
rz(-2.0514964) q[1];
sx q[1];
rz(-2.1665159) q[1];
sx q[1];
rz(2.8930194) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.45390689) q[0];
sx q[0];
rz(-1.0605863) q[0];
sx q[0];
rz(-1.8574255) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0806662) q[2];
sx q[2];
rz(-1.14398) q[2];
sx q[2];
rz(2.5271497) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.12411815) q[1];
sx q[1];
rz(-2.0548471) q[1];
sx q[1];
rz(0.89213051) q[1];
x q[2];
rz(-3.0646938) q[3];
sx q[3];
rz(-1.7401164) q[3];
sx q[3];
rz(0.83991915) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.8481855) q[2];
sx q[2];
rz(-0.53331393) q[2];
sx q[2];
rz(0.08671134) q[2];
rz(2.6596206) q[3];
sx q[3];
rz(-0.50589365) q[3];
sx q[3];
rz(1.1530676) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6329353) q[0];
sx q[0];
rz(-2.1338699) q[0];
sx q[0];
rz(-0.28276643) q[0];
rz(2.4400318) q[1];
sx q[1];
rz(-0.82156721) q[1];
sx q[1];
rz(1.3185906) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.628172) q[0];
sx q[0];
rz(-1.4903729) q[0];
sx q[0];
rz(-1.1286939) q[0];
rz(-pi) q[1];
rz(-0.73762383) q[2];
sx q[2];
rz(-0.80543033) q[2];
sx q[2];
rz(-1.1951624) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.2576037) q[1];
sx q[1];
rz(-2.249243) q[1];
sx q[1];
rz(0.18737327) q[1];
rz(0.41140822) q[3];
sx q[3];
rz(-2.3477926) q[3];
sx q[3];
rz(1.0994764) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.7302154) q[2];
sx q[2];
rz(-0.29835478) q[2];
sx q[2];
rz(-2.7424157) q[2];
rz(2.2579851) q[3];
sx q[3];
rz(-1.4727605) q[3];
sx q[3];
rz(1.9201027) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76319641) q[0];
sx q[0];
rz(-2.3662687) q[0];
sx q[0];
rz(-1.5989074) q[0];
rz(-2.0762766) q[1];
sx q[1];
rz(-2.1677446) q[1];
sx q[1];
rz(-1.261196) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6529918) q[0];
sx q[0];
rz(-1.6761259) q[0];
sx q[0];
rz(-0.41098849) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5485498) q[2];
sx q[2];
rz(-0.97264475) q[2];
sx q[2];
rz(-1.0215789) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.65834261) q[1];
sx q[1];
rz(-1.5066506) q[1];
sx q[1];
rz(2.642753) q[1];
x q[2];
rz(-2.946978) q[3];
sx q[3];
rz(-0.6932887) q[3];
sx q[3];
rz(-0.64893901) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.5836872) q[2];
sx q[2];
rz(-2.5186899) q[2];
sx q[2];
rz(-2.5718001) q[2];
rz(1.2184881) q[3];
sx q[3];
rz(-2.4945538) q[3];
sx q[3];
rz(-0.56263721) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31310836) q[0];
sx q[0];
rz(-0.88043558) q[0];
sx q[0];
rz(1.2783929) q[0];
rz(0.60824153) q[1];
sx q[1];
rz(-2.6651762) q[1];
sx q[1];
rz(-2.6574635) q[1];
rz(-0.026253168) q[2];
sx q[2];
rz(-0.99925169) q[2];
sx q[2];
rz(3.0796438) q[2];
rz(-1.431987) q[3];
sx q[3];
rz(-0.98416735) q[3];
sx q[3];
rz(-2.9744801) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
