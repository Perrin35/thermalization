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
rz(0.53606755) q[0];
sx q[0];
rz(7.4458691) q[0];
sx q[0];
rz(9.8557628) q[0];
rz(0.35297901) q[1];
sx q[1];
rz(-1.9184435) q[1];
sx q[1];
rz(0.42981848) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0551222) q[0];
sx q[0];
rz(-1.4807379) q[0];
sx q[0];
rz(-0.49205972) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5360898) q[2];
sx q[2];
rz(-2.7089951) q[2];
sx q[2];
rz(-2.2207867) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.7215772) q[1];
sx q[1];
rz(-0.41999451) q[1];
sx q[1];
rz(2.6073522) q[1];
x q[2];
rz(1.1898891) q[3];
sx q[3];
rz(-2.0365247) q[3];
sx q[3];
rz(-1.7349752) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.063804403) q[2];
sx q[2];
rz(-2.6142075) q[2];
sx q[2];
rz(-1.0112313) q[2];
rz(-2.1003335) q[3];
sx q[3];
rz(-0.51308378) q[3];
sx q[3];
rz(2.7019555) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0332727) q[0];
sx q[0];
rz(-0.99650538) q[0];
sx q[0];
rz(2.394115) q[0];
rz(0.29391089) q[1];
sx q[1];
rz(-1.1312609) q[1];
sx q[1];
rz(2.8864313) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0145484) q[0];
sx q[0];
rz(-0.60227312) q[0];
sx q[0];
rz(1.2748347) q[0];
x q[1];
rz(2.601971) q[2];
sx q[2];
rz(-0.35225062) q[2];
sx q[2];
rz(0.43708153) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.16431674) q[1];
sx q[1];
rz(-2.212388) q[1];
sx q[1];
rz(-0.56820984) q[1];
rz(-pi) q[2];
rz(0.94631715) q[3];
sx q[3];
rz(-1.7229652) q[3];
sx q[3];
rz(0.63380235) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-3.1110288) q[2];
sx q[2];
rz(-0.054628987) q[2];
sx q[2];
rz(-0.89152208) q[2];
rz(2.8920178) q[3];
sx q[3];
rz(-2.2587743) q[3];
sx q[3];
rz(-2.7454624) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4075277) q[0];
sx q[0];
rz(-2.0052818) q[0];
sx q[0];
rz(-0.044128142) q[0];
rz(-1.8924425) q[1];
sx q[1];
rz(-0.42611486) q[1];
sx q[1];
rz(-2.6557907) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0853538) q[0];
sx q[0];
rz(-1.5512497) q[0];
sx q[0];
rz(-2.9162558) q[0];
rz(-pi) q[1];
rz(3.0032512) q[2];
sx q[2];
rz(-1.3855583) q[2];
sx q[2];
rz(-2.0092162) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.9618157) q[1];
sx q[1];
rz(-2.0992341) q[1];
sx q[1];
rz(-2.4584199) q[1];
x q[2];
rz(0.39704571) q[3];
sx q[3];
rz(-1.1685123) q[3];
sx q[3];
rz(2.8232529) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.32758632) q[2];
sx q[2];
rz(-2.0531211) q[2];
sx q[2];
rz(2.4468454) q[2];
rz(0.57573777) q[3];
sx q[3];
rz(-1.271194) q[3];
sx q[3];
rz(-1.4076788) q[3];
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
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.95076743) q[0];
sx q[0];
rz(-2.8468843) q[0];
sx q[0];
rz(2.8488391) q[0];
rz(0.5689019) q[1];
sx q[1];
rz(-2.3141373) q[1];
sx q[1];
rz(0.84691602) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0305188) q[0];
sx q[0];
rz(-1.9785709) q[0];
sx q[0];
rz(3.0570583) q[0];
rz(-0.072242468) q[2];
sx q[2];
rz(-2.2390167) q[2];
sx q[2];
rz(2.7492858) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.12754) q[1];
sx q[1];
rz(-1.7636646) q[1];
sx q[1];
rz(-2.8991153) q[1];
rz(-1.9027223) q[3];
sx q[3];
rz(-0.20461018) q[3];
sx q[3];
rz(2.0505035) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.34951052) q[2];
sx q[2];
rz(-1.8831207) q[2];
sx q[2];
rz(-0.2087896) q[2];
rz(2.4236692) q[3];
sx q[3];
rz(-0.28306285) q[3];
sx q[3];
rz(-0.94055241) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1013041) q[0];
sx q[0];
rz(-0.52241075) q[0];
sx q[0];
rz(2.8029602) q[0];
rz(-0.70308095) q[1];
sx q[1];
rz(-2.4026726) q[1];
sx q[1];
rz(2.3032761) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.188952) q[0];
sx q[0];
rz(-1.8243655) q[0];
sx q[0];
rz(0.44373657) q[0];
rz(-2.9432326) q[2];
sx q[2];
rz(-2.3075929) q[2];
sx q[2];
rz(1.5892346) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.014480249) q[1];
sx q[1];
rz(-2.4408484) q[1];
sx q[1];
rz(1.5190009) q[1];
rz(0.57788604) q[3];
sx q[3];
rz(-0.69212038) q[3];
sx q[3];
rz(-2.1731165) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.4579953) q[2];
sx q[2];
rz(-2.0532875) q[2];
sx q[2];
rz(-2.1387157) q[2];
rz(2.9966127) q[3];
sx q[3];
rz(-3.0478015) q[3];
sx q[3];
rz(0.063044757) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.018464) q[0];
sx q[0];
rz(-0.58582425) q[0];
sx q[0];
rz(0.96328324) q[0];
rz(-2.9604498) q[1];
sx q[1];
rz(-2.2434442) q[1];
sx q[1];
rz(1.2394261) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.081112208) q[0];
sx q[0];
rz(-1.1976997) q[0];
sx q[0];
rz(-0.65848559) q[0];
rz(-pi) q[1];
rz(-0.507429) q[2];
sx q[2];
rz(-2.6118738) q[2];
sx q[2];
rz(-1.5786174) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.0964862) q[1];
sx q[1];
rz(-0.54212057) q[1];
sx q[1];
rz(0.78150918) q[1];
rz(-pi) q[2];
rz(0.12902) q[3];
sx q[3];
rz(-1.4203398) q[3];
sx q[3];
rz(-0.64808382) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.88508254) q[2];
sx q[2];
rz(-2.2233267) q[2];
sx q[2];
rz(0.93290848) q[2];
rz(-1.4126623) q[3];
sx q[3];
rz(-1.4224854) q[3];
sx q[3];
rz(0.2093813) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4466208) q[0];
sx q[0];
rz(-0.3681204) q[0];
sx q[0];
rz(-0.80286655) q[0];
rz(2.39957) q[1];
sx q[1];
rz(-1.4890393) q[1];
sx q[1];
rz(-0.41025695) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0502346) q[0];
sx q[0];
rz(-1.0438393) q[0];
sx q[0];
rz(-1.1731253) q[0];
rz(-3.018961) q[2];
sx q[2];
rz(-0.76803128) q[2];
sx q[2];
rz(-1.0379858) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.81933016) q[1];
sx q[1];
rz(-2.9633489) q[1];
sx q[1];
rz(1.1398386) q[1];
x q[2];
rz(1.0954082) q[3];
sx q[3];
rz(-0.96964004) q[3];
sx q[3];
rz(2.9132995) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.2731169) q[2];
sx q[2];
rz(-1.5831524) q[2];
sx q[2];
rz(-2.9203019) q[2];
rz(-2.7224329) q[3];
sx q[3];
rz(-0.48210382) q[3];
sx q[3];
rz(-0.47590772) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9047423) q[0];
sx q[0];
rz(-2.1433647) q[0];
sx q[0];
rz(-1.6118443) q[0];
rz(1.3329685) q[1];
sx q[1];
rz(-0.95104853) q[1];
sx q[1];
rz(-1.6882247) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7573563) q[0];
sx q[0];
rz(-1.0125375) q[0];
sx q[0];
rz(-2.4008958) q[0];
x q[1];
rz(-2.9506893) q[2];
sx q[2];
rz(-1.4530572) q[2];
sx q[2];
rz(-2.8536053) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.087880922) q[1];
sx q[1];
rz(-2.2715927) q[1];
sx q[1];
rz(-1.7675964) q[1];
rz(2.0734208) q[3];
sx q[3];
rz(-0.31538559) q[3];
sx q[3];
rz(-2.634545) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.86614418) q[2];
sx q[2];
rz(-2.8871705) q[2];
sx q[2];
rz(-0.010738372) q[2];
rz(-0.3206611) q[3];
sx q[3];
rz(-1.2977707) q[3];
sx q[3];
rz(0.58967725) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7084932) q[0];
sx q[0];
rz(-0.95516094) q[0];
sx q[0];
rz(0.47472111) q[0];
rz(0.2933329) q[1];
sx q[1];
rz(-1.1055929) q[1];
sx q[1];
rz(-1.4795823) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0728481) q[0];
sx q[0];
rz(-1.6075875) q[0];
sx q[0];
rz(-1.6108914) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8443498) q[2];
sx q[2];
rz(-2.4435776) q[2];
sx q[2];
rz(2.6797803) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.868728) q[1];
sx q[1];
rz(-2.0365148) q[1];
sx q[1];
rz(-1.8293423) q[1];
rz(-pi) q[2];
x q[2];
rz(1.817486) q[3];
sx q[3];
rz(-1.7462915) q[3];
sx q[3];
rz(-1.3934203) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.044363) q[2];
sx q[2];
rz(-1.2791415) q[2];
sx q[2];
rz(-0.43107671) q[2];
rz(-0.19032446) q[3];
sx q[3];
rz(-0.35076916) q[3];
sx q[3];
rz(-0.88461191) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.97656074) q[0];
sx q[0];
rz(-2.850552) q[0];
sx q[0];
rz(2.9227559) q[0];
rz(2.1826375) q[1];
sx q[1];
rz(-2.405165) q[1];
sx q[1];
rz(1.3921907) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7277273) q[0];
sx q[0];
rz(-1.6221592) q[0];
sx q[0];
rz(0.19330175) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8270589) q[2];
sx q[2];
rz(-2.7022903) q[2];
sx q[2];
rz(-0.9010992) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.47962418) q[1];
sx q[1];
rz(-1.8682518) q[1];
sx q[1];
rz(-0.65315078) q[1];
rz(-pi) q[2];
rz(0.99361657) q[3];
sx q[3];
rz(-2.4419067) q[3];
sx q[3];
rz(2.4260739) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.2739111) q[2];
sx q[2];
rz(-0.76649222) q[2];
sx q[2];
rz(-1.8217746) q[2];
rz(2.6093318) q[3];
sx q[3];
rz(-1.7936734) q[3];
sx q[3];
rz(1.5925315) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5821447) q[0];
sx q[0];
rz(-1.5912709) q[0];
sx q[0];
rz(-2.7406319) q[0];
rz(-2.9136912) q[1];
sx q[1];
rz(-2.0961998) q[1];
sx q[1];
rz(-1.6963522) q[1];
rz(-2.7264222) q[2];
sx q[2];
rz(-1.3317458) q[2];
sx q[2];
rz(1.7785704) q[2];
rz(0.99146546) q[3];
sx q[3];
rz(-1.4351063) q[3];
sx q[3];
rz(-1.2096005) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
