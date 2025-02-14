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
rz(0.0079732815) q[0];
sx q[0];
rz(-0.43626269) q[0];
sx q[0];
rz(2.1313957) q[0];
rz(0.17207347) q[1];
sx q[1];
rz(-2.3354524) q[1];
sx q[1];
rz(2.1176718) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14903325) q[0];
sx q[0];
rz(-1.5776538) q[0];
sx q[0];
rz(-2.0586781) q[0];
rz(0.51895835) q[2];
sx q[2];
rz(-0.4441688) q[2];
sx q[2];
rz(-2.3553395) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.0144452) q[1];
sx q[1];
rz(-1.5902441) q[1];
sx q[1];
rz(1.4102742) q[1];
rz(2.8918582) q[3];
sx q[3];
rz(-1.6563376) q[3];
sx q[3];
rz(1.6289323) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.94850928) q[2];
sx q[2];
rz(-2.812959) q[2];
sx q[2];
rz(2.9060717) q[2];
rz(-2.113302) q[3];
sx q[3];
rz(-1.2582658) q[3];
sx q[3];
rz(1.1576912) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.704772) q[0];
sx q[0];
rz(-1.3587767) q[0];
sx q[0];
rz(-0.92192465) q[0];
rz(3.0251265) q[1];
sx q[1];
rz(-1.7200229) q[1];
sx q[1];
rz(0.88752735) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.015892643) q[0];
sx q[0];
rz(-1.3845446) q[0];
sx q[0];
rz(0.13293477) q[0];
rz(1.5121801) q[2];
sx q[2];
rz(-1.2208856) q[2];
sx q[2];
rz(1.0998124) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.077558806) q[1];
sx q[1];
rz(-2.56018) q[1];
sx q[1];
rz(-0.20661776) q[1];
x q[2];
rz(-1.4027897) q[3];
sx q[3];
rz(-0.6989494) q[3];
sx q[3];
rz(0.93937554) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.55887115) q[2];
sx q[2];
rz(-2.116394) q[2];
sx q[2];
rz(-2.9026046) q[2];
rz(-2.41921) q[3];
sx q[3];
rz(-0.84010774) q[3];
sx q[3];
rz(1.1649789) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31579414) q[0];
sx q[0];
rz(-0.9676942) q[0];
sx q[0];
rz(0.78829366) q[0];
rz(1.7514508) q[1];
sx q[1];
rz(-0.43594053) q[1];
sx q[1];
rz(-1.699126) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9787131) q[0];
sx q[0];
rz(-1.4294693) q[0];
sx q[0];
rz(-3.0706329) q[0];
rz(-pi) q[1];
rz(-2.2295775) q[2];
sx q[2];
rz(-2.680353) q[2];
sx q[2];
rz(1.3818503) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.07787598) q[1];
sx q[1];
rz(-2.2635679) q[1];
sx q[1];
rz(0.31097842) q[1];
x q[2];
rz(1.0679108) q[3];
sx q[3];
rz(-1.7278132) q[3];
sx q[3];
rz(-0.63105052) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.4060789) q[2];
sx q[2];
rz(-1.9903851) q[2];
sx q[2];
rz(0.25685143) q[2];
rz(0.86366051) q[3];
sx q[3];
rz(-2.05859) q[3];
sx q[3];
rz(-0.5736205) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0470806) q[0];
sx q[0];
rz(-0.032945078) q[0];
sx q[0];
rz(-1.6324014) q[0];
rz(-2.8250561) q[1];
sx q[1];
rz(-1.5959975) q[1];
sx q[1];
rz(-1.0155771) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2779396) q[0];
sx q[0];
rz(-2.4653593) q[0];
sx q[0];
rz(2.3689224) q[0];
rz(1.847083) q[2];
sx q[2];
rz(-2.3450548) q[2];
sx q[2];
rz(-0.73819064) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.16962651) q[1];
sx q[1];
rz(-1.3463273) q[1];
sx q[1];
rz(-1.1826993) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9731941) q[3];
sx q[3];
rz(-0.83074289) q[3];
sx q[3];
rz(-0.5357045) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.8429026) q[2];
sx q[2];
rz(-2.0899253) q[2];
sx q[2];
rz(-1.451937) q[2];
rz(1.2339833) q[3];
sx q[3];
rz(-2.0472287) q[3];
sx q[3];
rz(2.2598677) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
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
rz(-1.7538309) q[0];
sx q[0];
rz(-1.8330638) q[0];
sx q[0];
rz(-2.0552788) q[0];
rz(-1.4007252) q[1];
sx q[1];
rz(-2.3298405) q[1];
sx q[1];
rz(0.0033671826) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8683748) q[0];
sx q[0];
rz(-1.5766607) q[0];
sx q[0];
rz(2.547079) q[0];
rz(-pi) q[1];
rz(2.8068266) q[2];
sx q[2];
rz(-2.1328983) q[2];
sx q[2];
rz(2.7239852) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.6033404) q[1];
sx q[1];
rz(-1.0639186) q[1];
sx q[1];
rz(-1.6350621) q[1];
rz(-pi) q[2];
rz(-0.60008757) q[3];
sx q[3];
rz(-1.5494124) q[3];
sx q[3];
rz(0.96259538) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.9747666) q[2];
sx q[2];
rz(-0.70009118) q[2];
sx q[2];
rz(2.8066446) q[2];
rz(-0.29609984) q[3];
sx q[3];
rz(-2.1985603) q[3];
sx q[3];
rz(0.12224841) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.629338) q[0];
sx q[0];
rz(-0.4991931) q[0];
sx q[0];
rz(-0.41803023) q[0];
rz(-0.66608518) q[1];
sx q[1];
rz(-1.2434375) q[1];
sx q[1];
rz(-0.24806771) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9739415) q[0];
sx q[0];
rz(-2.535916) q[0];
sx q[0];
rz(-2.1781237) q[0];
rz(-pi) q[1];
rz(-0.054912189) q[2];
sx q[2];
rz(-1.1217199) q[2];
sx q[2];
rz(-0.95740805) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.0679761) q[1];
sx q[1];
rz(-2.4779137) q[1];
sx q[1];
rz(-0.90866088) q[1];
rz(-pi) q[2];
x q[2];
rz(0.94132386) q[3];
sx q[3];
rz(-1.7793312) q[3];
sx q[3];
rz(-1.0104313) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.4416113) q[2];
sx q[2];
rz(-1.8786083) q[2];
sx q[2];
rz(1.9786394) q[2];
rz(0.59398389) q[3];
sx q[3];
rz(-1.6006399) q[3];
sx q[3];
rz(2.1514413) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.537259) q[0];
sx q[0];
rz(-2.8626677) q[0];
sx q[0];
rz(0.064924031) q[0];
rz(-2.7760432) q[1];
sx q[1];
rz(-1.7100916) q[1];
sx q[1];
rz(-2.4117267) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6014272) q[0];
sx q[0];
rz(-1.4608723) q[0];
sx q[0];
rz(-0.92782995) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7634298) q[2];
sx q[2];
rz(-0.81708497) q[2];
sx q[2];
rz(-2.0157004) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.0522451) q[1];
sx q[1];
rz(-2.5762853) q[1];
sx q[1];
rz(-1.7491399) q[1];
x q[2];
rz(-2.6486126) q[3];
sx q[3];
rz(-0.74876174) q[3];
sx q[3];
rz(-1.6822588) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.2713752) q[2];
sx q[2];
rz(-1.7729746) q[2];
sx q[2];
rz(2.0591056) q[2];
rz(1.6535951) q[3];
sx q[3];
rz(-2.7797785) q[3];
sx q[3];
rz(0.35058072) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38700405) q[0];
sx q[0];
rz(-2.258774) q[0];
sx q[0];
rz(-2.8104695) q[0];
rz(1.6638727) q[1];
sx q[1];
rz(-0.96550566) q[1];
sx q[1];
rz(0.32041916) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3382771) q[0];
sx q[0];
rz(-1.471963) q[0];
sx q[0];
rz(-1.5966948) q[0];
rz(-0.1542145) q[2];
sx q[2];
rz(-1.4262098) q[2];
sx q[2];
rz(0.54942451) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.3255405) q[1];
sx q[1];
rz(-2.7150702) q[1];
sx q[1];
rz(3.1040621) q[1];
rz(-pi) q[2];
rz(3.0842264) q[3];
sx q[3];
rz(-2.2542977) q[3];
sx q[3];
rz(0.43903109) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.42935818) q[2];
sx q[2];
rz(-1.7244491) q[2];
sx q[2];
rz(-0.36163914) q[2];
rz(-2.2278191) q[3];
sx q[3];
rz(-0.29699609) q[3];
sx q[3];
rz(0.44338068) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2509895) q[0];
sx q[0];
rz(-2.3539703) q[0];
sx q[0];
rz(-3.0255764) q[0];
rz(-1.3277671) q[1];
sx q[1];
rz(-1.9783744) q[1];
sx q[1];
rz(1.9414925) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3108771) q[0];
sx q[0];
rz(-2.9492464) q[0];
sx q[0];
rz(3.0430692) q[0];
x q[1];
rz(1.7723254) q[2];
sx q[2];
rz(-0.69567615) q[2];
sx q[2];
rz(-1.8303378) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.6873106) q[1];
sx q[1];
rz(-0.56182278) q[1];
sx q[1];
rz(3.1141206) q[1];
rz(-pi) q[2];
rz(-1.8103792) q[3];
sx q[3];
rz(-1.7322043) q[3];
sx q[3];
rz(2.7454698) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.8414595) q[2];
sx q[2];
rz(-1.7068784) q[2];
sx q[2];
rz(0.33162281) q[2];
rz(2.8578109) q[3];
sx q[3];
rz(-2.0399358) q[3];
sx q[3];
rz(2.7200429) q[3];
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
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.83061853) q[0];
sx q[0];
rz(-2.0299439) q[0];
sx q[0];
rz(-2.9750138) q[0];
rz(-0.84959787) q[1];
sx q[1];
rz(-2.6764937) q[1];
sx q[1];
rz(3.0679682) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66921418) q[0];
sx q[0];
rz(-1.3166691) q[0];
sx q[0];
rz(-2.4082802) q[0];
rz(2.2715015) q[2];
sx q[2];
rz(-0.67591643) q[2];
sx q[2];
rz(1.6411587) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.6839809) q[1];
sx q[1];
rz(-1.906412) q[1];
sx q[1];
rz(2.6102351) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0898548) q[3];
sx q[3];
rz(-2.1409799) q[3];
sx q[3];
rz(0.4057623) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.2354043) q[2];
sx q[2];
rz(-0.81934682) q[2];
sx q[2];
rz(-2.4833615) q[2];
rz(1.341358) q[3];
sx q[3];
rz(-0.96786371) q[3];
sx q[3];
rz(-0.85097504) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3157208) q[0];
sx q[0];
rz(-0.82875874) q[0];
sx q[0];
rz(2.2525633) q[0];
rz(-0.52147621) q[1];
sx q[1];
rz(-1.5670525) q[1];
sx q[1];
rz(-1.570931) q[1];
rz(-0.53311382) q[2];
sx q[2];
rz(-0.86400142) q[2];
sx q[2];
rz(0.17178847) q[2];
rz(-2.9110649) q[3];
sx q[3];
rz(-2.5067825) q[3];
sx q[3];
rz(1.4092177) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
