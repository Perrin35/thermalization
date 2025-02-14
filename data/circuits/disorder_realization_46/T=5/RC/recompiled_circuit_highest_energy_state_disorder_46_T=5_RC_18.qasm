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
rz(-1.3396076) q[0];
sx q[0];
rz(-0.34914246) q[0];
sx q[0];
rz(2.1139297) q[0];
rz(3.1066306) q[1];
sx q[1];
rz(-2.1244013) q[1];
sx q[1];
rz(1.4453759) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.18939173) q[0];
sx q[0];
rz(-0.9993573) q[0];
sx q[0];
rz(-1.9791358) q[0];
x q[1];
rz(1.6048172) q[2];
sx q[2];
rz(-1.7482702) q[2];
sx q[2];
rz(-0.79171514) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.4852745) q[1];
sx q[1];
rz(-2.7492011) q[1];
sx q[1];
rz(-2.5618784) q[1];
x q[2];
rz(1.2054005) q[3];
sx q[3];
rz(-2.2789848) q[3];
sx q[3];
rz(-2.303249) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.33862904) q[2];
sx q[2];
rz(-0.45946071) q[2];
sx q[2];
rz(2.6098693) q[2];
rz(1.3945329) q[3];
sx q[3];
rz(-1.2732384) q[3];
sx q[3];
rz(-1.9664221) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26674536) q[0];
sx q[0];
rz(-1.6371472) q[0];
sx q[0];
rz(-2.3496085) q[0];
rz(2.2199471) q[1];
sx q[1];
rz(-1.6519203) q[1];
sx q[1];
rz(-1.1080866) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8130694) q[0];
sx q[0];
rz(-1.4591154) q[0];
sx q[0];
rz(0.96891667) q[0];
rz(-pi) q[1];
x q[1];
rz(0.048594995) q[2];
sx q[2];
rz(-2.193795) q[2];
sx q[2];
rz(0.6374661) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.0818357) q[1];
sx q[1];
rz(-1.949777) q[1];
sx q[1];
rz(0.80767085) q[1];
rz(-1.8872126) q[3];
sx q[3];
rz(-1.5174447) q[3];
sx q[3];
rz(2.2406468) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.7532928) q[2];
sx q[2];
rz(-2.1123977) q[2];
sx q[2];
rz(-1.252906) q[2];
rz(-2.5168354) q[3];
sx q[3];
rz(-1.2478991) q[3];
sx q[3];
rz(0.46318444) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6817634) q[0];
sx q[0];
rz(-0.1777996) q[0];
sx q[0];
rz(-2.8681927) q[0];
rz(-3.0699442) q[1];
sx q[1];
rz(-1.1337846) q[1];
sx q[1];
rz(-1.515548) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0833383) q[0];
sx q[0];
rz(-1.951955) q[0];
sx q[0];
rz(2.5883771) q[0];
rz(-2.5070174) q[2];
sx q[2];
rz(-2.5800642) q[2];
sx q[2];
rz(-0.47088366) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.3699941) q[1];
sx q[1];
rz(-2.150476) q[1];
sx q[1];
rz(1.3808151) q[1];
x q[2];
rz(1.203556) q[3];
sx q[3];
rz(-0.62588309) q[3];
sx q[3];
rz(0.82928951) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.7800954) q[2];
sx q[2];
rz(-0.7146892) q[2];
sx q[2];
rz(-1.7717465) q[2];
rz(-1.0684377) q[3];
sx q[3];
rz(-1.7504102) q[3];
sx q[3];
rz(-0.95503241) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.84900981) q[0];
sx q[0];
rz(-2.7136901) q[0];
sx q[0];
rz(-0.12246116) q[0];
rz(3.124584) q[1];
sx q[1];
rz(-1.6288792) q[1];
sx q[1];
rz(-2.2564127) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2816078) q[0];
sx q[0];
rz(-2.6329319) q[0];
sx q[0];
rz(1.2889494) q[0];
x q[1];
rz(1.1954284) q[2];
sx q[2];
rz(-2.3705707) q[2];
sx q[2];
rz(0.32685243) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.7834691) q[1];
sx q[1];
rz(-2.0253638) q[1];
sx q[1];
rz(0.44444167) q[1];
rz(-1.3766798) q[3];
sx q[3];
rz(-2.2633865) q[3];
sx q[3];
rz(-0.10759456) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.433832) q[2];
sx q[2];
rz(-1.8261352) q[2];
sx q[2];
rz(0.07746499) q[2];
rz(-0.42099434) q[3];
sx q[3];
rz(-1.1869895) q[3];
sx q[3];
rz(2.8903956) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0635327) q[0];
sx q[0];
rz(-3.1227626) q[0];
sx q[0];
rz(-2.4596762) q[0];
rz(-2.6497427) q[1];
sx q[1];
rz(-1.0268772) q[1];
sx q[1];
rz(1.9815725) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9043961) q[0];
sx q[0];
rz(-2.2771796) q[0];
sx q[0];
rz(-0.24817962) q[0];
x q[1];
rz(2.5712058) q[2];
sx q[2];
rz(-1.5231992) q[2];
sx q[2];
rz(-2.334495) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.8012538) q[1];
sx q[1];
rz(-1.7427708) q[1];
sx q[1];
rz(0.94229631) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1669768) q[3];
sx q[3];
rz(-0.49562956) q[3];
sx q[3];
rz(2.0295124) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.68449768) q[2];
sx q[2];
rz(-0.86898494) q[2];
sx q[2];
rz(-2.1333466) q[2];
rz(1.5021987) q[3];
sx q[3];
rz(-2.0366171) q[3];
sx q[3];
rz(0.26386279) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0428001) q[0];
sx q[0];
rz(-1.5473939) q[0];
sx q[0];
rz(0.40476009) q[0];
rz(-2.3233991) q[1];
sx q[1];
rz(-1.300756) q[1];
sx q[1];
rz(-2.4249446) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9500026) q[0];
sx q[0];
rz(-1.5609589) q[0];
sx q[0];
rz(-2.8623926) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0990853) q[2];
sx q[2];
rz(-0.87954885) q[2];
sx q[2];
rz(-1.7886358) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.08581768) q[1];
sx q[1];
rz(-1.1113864) q[1];
sx q[1];
rz(-1.2797848) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6022609) q[3];
sx q[3];
rz(-1.2507273) q[3];
sx q[3];
rz(2.0254997) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.8019668) q[2];
sx q[2];
rz(-0.64029396) q[2];
sx q[2];
rz(-0.65529811) q[2];
rz(-0.35704923) q[3];
sx q[3];
rz(-1.7506295) q[3];
sx q[3];
rz(0.51957875) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.201467) q[0];
sx q[0];
rz(-0.31620142) q[0];
sx q[0];
rz(2.1215718) q[0];
rz(-0.54234281) q[1];
sx q[1];
rz(-1.1154563) q[1];
sx q[1];
rz(1.7592336) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.64847222) q[0];
sx q[0];
rz(-2.1628404) q[0];
sx q[0];
rz(2.4061794) q[0];
rz(-pi) q[1];
rz(0.57047259) q[2];
sx q[2];
rz(-0.24952023) q[2];
sx q[2];
rz(0.0970627) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.51782896) q[1];
sx q[1];
rz(-1.6490855) q[1];
sx q[1];
rz(1.109156) q[1];
x q[2];
rz(-1.5730619) q[3];
sx q[3];
rz(-1.6544154) q[3];
sx q[3];
rz(-0.22606848) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.0688087) q[2];
sx q[2];
rz(-1.5280318) q[2];
sx q[2];
rz(-0.35045785) q[2];
rz(2.6751878) q[3];
sx q[3];
rz(-2.2240708) q[3];
sx q[3];
rz(0.94356999) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35895178) q[0];
sx q[0];
rz(-2.751001) q[0];
sx q[0];
rz(-0.95640957) q[0];
rz(1.6920754) q[1];
sx q[1];
rz(-0.78069514) q[1];
sx q[1];
rz(-0.67109674) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3365375) q[0];
sx q[0];
rz(-1.4292246) q[0];
sx q[0];
rz(-0.086351589) q[0];
x q[1];
rz(-2.8196149) q[2];
sx q[2];
rz(-0.59315943) q[2];
sx q[2];
rz(-0.71074394) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.8099212) q[1];
sx q[1];
rz(-1.6345895) q[1];
sx q[1];
rz(1.2170736) q[1];
rz(-pi) q[2];
rz(2.8550034) q[3];
sx q[3];
rz(-0.72460382) q[3];
sx q[3];
rz(2.2514908) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.87390071) q[2];
sx q[2];
rz(-2.716422) q[2];
sx q[2];
rz(-2.0153866) q[2];
rz(-0.33024427) q[3];
sx q[3];
rz(-1.7405905) q[3];
sx q[3];
rz(-0.11708524) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0432334) q[0];
sx q[0];
rz(-2.5625304) q[0];
sx q[0];
rz(0.057057127) q[0];
rz(3.0077899) q[1];
sx q[1];
rz(-1.0452784) q[1];
sx q[1];
rz(-2.4287976) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0781435) q[0];
sx q[0];
rz(-2.8196555) q[0];
sx q[0];
rz(-2.9418895) q[0];
rz(1.1634689) q[2];
sx q[2];
rz(-1.9453269) q[2];
sx q[2];
rz(0.56331149) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.9660898) q[1];
sx q[1];
rz(-2.0398629) q[1];
sx q[1];
rz(-2.0211981) q[1];
rz(-pi) q[2];
rz(1.4737878) q[3];
sx q[3];
rz(-0.69655124) q[3];
sx q[3];
rz(2.4318621) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.5130676) q[2];
sx q[2];
rz(-1.8551989) q[2];
sx q[2];
rz(-3.0600424) q[2];
rz(1.9214123) q[3];
sx q[3];
rz(-0.27118513) q[3];
sx q[3];
rz(-2.0512106) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14209014) q[0];
sx q[0];
rz(-1.5220078) q[0];
sx q[0];
rz(1.581544) q[0];
rz(-2.0060495) q[1];
sx q[1];
rz(-1.138843) q[1];
sx q[1];
rz(-1.1467689) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1449908) q[0];
sx q[0];
rz(-1.0787316) q[0];
sx q[0];
rz(-2.4053615) q[0];
x q[1];
rz(-0.76374526) q[2];
sx q[2];
rz(-2.1830705) q[2];
sx q[2];
rz(-0.51240048) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.4301871) q[1];
sx q[1];
rz(-1.3608977) q[1];
sx q[1];
rz(-1.6564356) q[1];
x q[2];
rz(-2.1358397) q[3];
sx q[3];
rz(-1.4970137) q[3];
sx q[3];
rz(-2.0948054) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.59209383) q[2];
sx q[2];
rz(-1.3201951) q[2];
sx q[2];
rz(-0.36592323) q[2];
rz(0.38771114) q[3];
sx q[3];
rz(-2.0230899) q[3];
sx q[3];
rz(-1.4661192) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.53030071) q[0];
sx q[0];
rz(-2.3983751) q[0];
sx q[0];
rz(-0.30839738) q[0];
rz(0.74956924) q[1];
sx q[1];
rz(-2.0105965) q[1];
sx q[1];
rz(-1.2731193) q[1];
rz(0.91084008) q[2];
sx q[2];
rz(-1.2963933) q[2];
sx q[2];
rz(-2.8252841) q[2];
rz(-1.9182792) q[3];
sx q[3];
rz(-2.1741304) q[3];
sx q[3];
rz(2.771467) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
