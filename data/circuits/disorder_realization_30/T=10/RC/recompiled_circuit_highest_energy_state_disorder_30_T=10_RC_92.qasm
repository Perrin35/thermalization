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
rz(-2.2089145) q[0];
sx q[0];
rz(-0.39251602) q[0];
sx q[0];
rz(1.0970595) q[0];
rz(-1.8707844) q[1];
sx q[1];
rz(4.3391736) q[1];
sx q[1];
rz(14.270562) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28406806) q[0];
sx q[0];
rz(-1.2204613) q[0];
sx q[0];
rz(-2.461238) q[0];
rz(-pi) q[1];
rz(-1.497606) q[2];
sx q[2];
rz(-1.3393667) q[2];
sx q[2];
rz(2.481593) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.7667707) q[1];
sx q[1];
rz(-1.8954191) q[1];
sx q[1];
rz(2.2036821) q[1];
rz(3.1278856) q[3];
sx q[3];
rz(-0.64798149) q[3];
sx q[3];
rz(1.614384) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.2854332) q[2];
sx q[2];
rz(-0.23468748) q[2];
sx q[2];
rz(2.6402546) q[2];
rz(-1.8173789) q[3];
sx q[3];
rz(-2.5960077) q[3];
sx q[3];
rz(1.6747624) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.45601869) q[0];
sx q[0];
rz(-2.7360003) q[0];
sx q[0];
rz(-1.6803886) q[0];
rz(0.17924084) q[1];
sx q[1];
rz(-2.3724809) q[1];
sx q[1];
rz(-0.63757149) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9623827) q[0];
sx q[0];
rz(-1.3342627) q[0];
sx q[0];
rz(-1.80577) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4316971) q[2];
sx q[2];
rz(-2.8219328) q[2];
sx q[2];
rz(0.57397288) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.3481787) q[1];
sx q[1];
rz(-1.2169588) q[1];
sx q[1];
rz(-2.2458312) q[1];
x q[2];
rz(-2.6735503) q[3];
sx q[3];
rz(-2.4222825) q[3];
sx q[3];
rz(-1.5507789) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.8403988) q[2];
sx q[2];
rz(-2.0638778) q[2];
sx q[2];
rz(-0.044142874) q[2];
rz(0.61331493) q[3];
sx q[3];
rz(-1.3215348) q[3];
sx q[3];
rz(0.80673748) q[3];
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
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.052385656) q[0];
sx q[0];
rz(-3.0969924) q[0];
sx q[0];
rz(-1.7983623) q[0];
rz(1.4187468) q[1];
sx q[1];
rz(-2.2371465) q[1];
sx q[1];
rz(-2.338063) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.270954) q[0];
sx q[0];
rz(-0.89193908) q[0];
sx q[0];
rz(1.4951597) q[0];
rz(-pi) q[1];
rz(0.64068303) q[2];
sx q[2];
rz(-2.4523458) q[2];
sx q[2];
rz(0.32459637) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.61345873) q[1];
sx q[1];
rz(-2.8430004) q[1];
sx q[1];
rz(-1.9419844) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.66763287) q[3];
sx q[3];
rz(-2.5270695) q[3];
sx q[3];
rz(-0.4947084) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.58328763) q[2];
sx q[2];
rz(-0.95197695) q[2];
sx q[2];
rz(-1.0179016) q[2];
rz(1.2244276) q[3];
sx q[3];
rz(-1.3320351) q[3];
sx q[3];
rz(-1.1289736) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99854904) q[0];
sx q[0];
rz(-2.3426549) q[0];
sx q[0];
rz(-0.96027389) q[0];
rz(2.9472561) q[1];
sx q[1];
rz(-1.4413709) q[1];
sx q[1];
rz(-3.1365373) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0329629) q[0];
sx q[0];
rz(-0.89187183) q[0];
sx q[0];
rz(2.4358948) q[0];
x q[1];
rz(-2.9580367) q[2];
sx q[2];
rz(-1.608299) q[2];
sx q[2];
rz(-1.3958193) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.3761085) q[1];
sx q[1];
rz(-1.6836201) q[1];
sx q[1];
rz(2.9214431) q[1];
rz(-2.9937135) q[3];
sx q[3];
rz(-2.085146) q[3];
sx q[3];
rz(1.0837931) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.2527577) q[2];
sx q[2];
rz(-1.5965261) q[2];
sx q[2];
rz(1.4463536) q[2];
rz(1.3085922) q[3];
sx q[3];
rz(-1.3819709) q[3];
sx q[3];
rz(2.5368209) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62110353) q[0];
sx q[0];
rz(-2.4172754) q[0];
sx q[0];
rz(-0.58897585) q[0];
rz(-3.1336054) q[1];
sx q[1];
rz(-2.0365448) q[1];
sx q[1];
rz(1.084682) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5863335) q[0];
sx q[0];
rz(-1.8249092) q[0];
sx q[0];
rz(1.3293367) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.1118664) q[2];
sx q[2];
rz(-2.8619302) q[2];
sx q[2];
rz(-3.0936808) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.1639203) q[1];
sx q[1];
rz(-2.2798988) q[1];
sx q[1];
rz(-2.6129938) q[1];
rz(-0.9081697) q[3];
sx q[3];
rz(-1.623933) q[3];
sx q[3];
rz(-0.65180627) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7320554) q[2];
sx q[2];
rz(-2.1246702) q[2];
sx q[2];
rz(-1.3746877) q[2];
rz(-3.0459259) q[3];
sx q[3];
rz(-1.5547662) q[3];
sx q[3];
rz(0.63053757) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-1.10504) q[0];
sx q[0];
rz(-0.51346546) q[0];
sx q[0];
rz(1.8010944) q[0];
rz(-2.4758677) q[1];
sx q[1];
rz(-1.7981073) q[1];
sx q[1];
rz(-0.7241157) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.501717) q[0];
sx q[0];
rz(-2.6611009) q[0];
sx q[0];
rz(2.2697422) q[0];
rz(-pi) q[1];
rz(2.1445455) q[2];
sx q[2];
rz(-1.0862175) q[2];
sx q[2];
rz(-0.73198971) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.49289068) q[1];
sx q[1];
rz(-2.6600983) q[1];
sx q[1];
rz(0.29701155) q[1];
rz(-pi) q[2];
x q[2];
rz(0.0047896623) q[3];
sx q[3];
rz(-1.4744722) q[3];
sx q[3];
rz(-0.19118689) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.93806997) q[2];
sx q[2];
rz(-1.4504434) q[2];
sx q[2];
rz(-3.0435437) q[2];
rz(2.8010662) q[3];
sx q[3];
rz(-0.90172684) q[3];
sx q[3];
rz(2.9448729) q[3];
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
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.53511867) q[0];
sx q[0];
rz(-0.70347324) q[0];
sx q[0];
rz(0.055543609) q[0];
rz(2.7659888) q[1];
sx q[1];
rz(-2.3665078) q[1];
sx q[1];
rz(1.6966049) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8880896) q[0];
sx q[0];
rz(-1.6763335) q[0];
sx q[0];
rz(0.77368931) q[0];
rz(1.7710502) q[2];
sx q[2];
rz(-0.25875388) q[2];
sx q[2];
rz(-3.0508397) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.58126175) q[1];
sx q[1];
rz(-2.7637097) q[1];
sx q[1];
rz(-0.037200971) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.7020601) q[3];
sx q[3];
rz(-2.5581048) q[3];
sx q[3];
rz(-1.4017915) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-3.0982509) q[2];
sx q[2];
rz(-1.8378374) q[2];
sx q[2];
rz(-1.8275758) q[2];
rz(0.029953778) q[3];
sx q[3];
rz(-1.5104537) q[3];
sx q[3];
rz(0.31908527) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.39716649) q[0];
sx q[0];
rz(-2.501896) q[0];
sx q[0];
rz(-1.2514914) q[0];
rz(-2.331612) q[1];
sx q[1];
rz(-0.73695838) q[1];
sx q[1];
rz(-1.4553778) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90235119) q[0];
sx q[0];
rz(-0.7860113) q[0];
sx q[0];
rz(2.8186574) q[0];
rz(-pi) q[1];
rz(-1.6258442) q[2];
sx q[2];
rz(-2.0924797) q[2];
sx q[2];
rz(1.5417527) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.4474214) q[1];
sx q[1];
rz(-0.93915597) q[1];
sx q[1];
rz(-0.76404913) q[1];
rz(-pi) q[2];
rz(0.87883605) q[3];
sx q[3];
rz(-0.86711001) q[3];
sx q[3];
rz(-1.4181942) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.2927085) q[2];
sx q[2];
rz(-1.7851189) q[2];
sx q[2];
rz(-1.9902309) q[2];
rz(-0.014160841) q[3];
sx q[3];
rz(-1.3580946) q[3];
sx q[3];
rz(-1.2667228) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8779811) q[0];
sx q[0];
rz(-0.71313715) q[0];
sx q[0];
rz(2.1942595) q[0];
rz(0.5961279) q[1];
sx q[1];
rz(-1.4214186) q[1];
sx q[1];
rz(1.5790342) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.57347711) q[0];
sx q[0];
rz(-1.7422424) q[0];
sx q[0];
rz(0.72400064) q[0];
rz(-pi) q[1];
rz(1.7315699) q[2];
sx q[2];
rz(-1.1005745) q[2];
sx q[2];
rz(1.1755113) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.5489997) q[1];
sx q[1];
rz(-2.3301417) q[1];
sx q[1];
rz(0.92632067) q[1];
rz(2.4189255) q[3];
sx q[3];
rz(-0.30084601) q[3];
sx q[3];
rz(2.2900555) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.8922213) q[2];
sx q[2];
rz(-2.5669079) q[2];
sx q[2];
rz(-0.48386827) q[2];
rz(-2.5441235) q[3];
sx q[3];
rz(-1.4474844) q[3];
sx q[3];
rz(-2.5634403) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4018965) q[0];
sx q[0];
rz(-1.551349) q[0];
sx q[0];
rz(0.36648146) q[0];
rz(1.81555) q[1];
sx q[1];
rz(-2.1791024) q[1];
sx q[1];
rz(2.0749626) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.90049998) q[0];
sx q[0];
rz(-2.0632732) q[0];
sx q[0];
rz(-0.89333138) q[0];
x q[1];
rz(-2.5353955) q[2];
sx q[2];
rz(-1.6087039) q[2];
sx q[2];
rz(2.4647922) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.047160867) q[1];
sx q[1];
rz(-2.4991813) q[1];
sx q[1];
rz(-1.98114) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0935547) q[3];
sx q[3];
rz(-1.541409) q[3];
sx q[3];
rz(0.060197006) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.2610953) q[2];
sx q[2];
rz(-2.3636621) q[2];
sx q[2];
rz(1.0658537) q[2];
rz(-2.8737658) q[3];
sx q[3];
rz(-1.4466176) q[3];
sx q[3];
rz(0.48521391) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5814701) q[0];
sx q[0];
rz(-1.0030092) q[0];
sx q[0];
rz(0.21448294) q[0];
rz(-0.32196925) q[1];
sx q[1];
rz(-1.8245158) q[1];
sx q[1];
rz(1.239924) q[1];
rz(1.5628846) q[2];
sx q[2];
rz(-1.2894165) q[2];
sx q[2];
rz(-3.0647562) q[2];
rz(-0.81372502) q[3];
sx q[3];
rz(-2.57434) q[3];
sx q[3];
rz(2.5634585) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
