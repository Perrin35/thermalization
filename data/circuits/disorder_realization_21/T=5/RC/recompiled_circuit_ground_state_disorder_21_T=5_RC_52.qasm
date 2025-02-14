OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.8747099) q[0];
sx q[0];
rz(4.2031718) q[0];
sx q[0];
rz(9.940552) q[0];
rz(-2.2401659) q[1];
sx q[1];
rz(-2.9401448) q[1];
sx q[1];
rz(-2.1405061) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2767746) q[0];
sx q[0];
rz(-1.0207812) q[0];
sx q[0];
rz(0.986542) q[0];
x q[1];
rz(-1.564525) q[2];
sx q[2];
rz(-1.1561818) q[2];
sx q[2];
rz(-0.27801499) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.6382003) q[1];
sx q[1];
rz(-0.93825785) q[1];
sx q[1];
rz(-1.0122385) q[1];
rz(2.4523425) q[3];
sx q[3];
rz(-0.19016506) q[3];
sx q[3];
rz(-0.81091698) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.72656816) q[2];
sx q[2];
rz(-2.6308306) q[2];
sx q[2];
rz(1.6815574) q[2];
rz(-2.7728752) q[3];
sx q[3];
rz(-1.5548778) q[3];
sx q[3];
rz(-1.4812428) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6129795) q[0];
sx q[0];
rz(-3.0280085) q[0];
sx q[0];
rz(-3.0533277) q[0];
rz(-2.2093692) q[1];
sx q[1];
rz(-1.1270707) q[1];
sx q[1];
rz(-0.50311911) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1189592) q[0];
sx q[0];
rz(-1.467407) q[0];
sx q[0];
rz(-2.6539283) q[0];
rz(-2.2406949) q[2];
sx q[2];
rz(-1.42282) q[2];
sx q[2];
rz(-1.0690546) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.4331665) q[1];
sx q[1];
rz(-1.6787156) q[1];
sx q[1];
rz(2.7923461) q[1];
rz(-pi) q[2];
x q[2];
rz(0.52946584) q[3];
sx q[3];
rz(-0.90116167) q[3];
sx q[3];
rz(-0.96082276) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.8253537) q[2];
sx q[2];
rz(-2.1672858) q[2];
sx q[2];
rz(0.45315722) q[2];
rz(0.00099269021) q[3];
sx q[3];
rz(-1.5627292) q[3];
sx q[3];
rz(2.4003975) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.93333018) q[0];
sx q[0];
rz(-0.18284155) q[0];
sx q[0];
rz(-2.6728447) q[0];
rz(0.58724976) q[1];
sx q[1];
rz(-1.8533555) q[1];
sx q[1];
rz(-2.8900878) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6856988) q[0];
sx q[0];
rz(-2.532428) q[0];
sx q[0];
rz(2.3745911) q[0];
rz(2.7526546) q[2];
sx q[2];
rz(-1.374482) q[2];
sx q[2];
rz(1.0599355) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.7461485) q[1];
sx q[1];
rz(-1.4471635) q[1];
sx q[1];
rz(-3.0828031) q[1];
rz(0.44106828) q[3];
sx q[3];
rz(-1.2510643) q[3];
sx q[3];
rz(-2.5351304) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.544439) q[2];
sx q[2];
rz(-2.4033098) q[2];
sx q[2];
rz(-0.047133751) q[2];
rz(-0.86826396) q[3];
sx q[3];
rz(-1.2406113) q[3];
sx q[3];
rz(-2.6422083) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.78859162) q[0];
sx q[0];
rz(-0.40591875) q[0];
sx q[0];
rz(1.8274008) q[0];
rz(0.81002533) q[1];
sx q[1];
rz(-1.6255197) q[1];
sx q[1];
rz(-2.8257418) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0421214) q[0];
sx q[0];
rz(-2.1232623) q[0];
sx q[0];
rz(-0.61037678) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.48607488) q[2];
sx q[2];
rz(-2.7856084) q[2];
sx q[2];
rz(2.4497368) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.12055732) q[1];
sx q[1];
rz(-1.8623973) q[1];
sx q[1];
rz(-2.3910644) q[1];
x q[2];
rz(2.9258435) q[3];
sx q[3];
rz(-2.5255754) q[3];
sx q[3];
rz(2.8293138) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.4874068) q[2];
sx q[2];
rz(-0.6260637) q[2];
sx q[2];
rz(-0.89228863) q[2];
rz(-1.6728632) q[3];
sx q[3];
rz(-1.059831) q[3];
sx q[3];
rz(1.0631801) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47657087) q[0];
sx q[0];
rz(-2.5044818) q[0];
sx q[0];
rz(-0.88660216) q[0];
rz(2.2992112) q[1];
sx q[1];
rz(-1.3096755) q[1];
sx q[1];
rz(-2.0319891) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.78326016) q[0];
sx q[0];
rz(-1.9985804) q[0];
sx q[0];
rz(0.88029998) q[0];
x q[1];
rz(-0.72317041) q[2];
sx q[2];
rz(-1.4338655) q[2];
sx q[2];
rz(-1.3201158) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.230727) q[1];
sx q[1];
rz(-0.80692569) q[1];
sx q[1];
rz(-0.74300503) q[1];
rz(-pi) q[2];
rz(2.4585633) q[3];
sx q[3];
rz(-1.08324) q[3];
sx q[3];
rz(1.4157996) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.2716918) q[2];
sx q[2];
rz(-2.7503408) q[2];
sx q[2];
rz(0.58763495) q[2];
rz(2.9247126) q[3];
sx q[3];
rz(-1.2809332) q[3];
sx q[3];
rz(1.7873526) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0103535) q[0];
sx q[0];
rz(-0.47061798) q[0];
sx q[0];
rz(-1.3838029) q[0];
rz(-0.17419392) q[1];
sx q[1];
rz(-2.0400338) q[1];
sx q[1];
rz(-1.9536288) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39858946) q[0];
sx q[0];
rz(-1.849353) q[0];
sx q[0];
rz(-0.085531959) q[0];
x q[1];
rz(-1.401004) q[2];
sx q[2];
rz(-1.0708969) q[2];
sx q[2];
rz(-0.57725805) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.9932258) q[1];
sx q[1];
rz(-1.6794551) q[1];
sx q[1];
rz(0.0309561) q[1];
rz(1.2634041) q[3];
sx q[3];
rz(-0.83019054) q[3];
sx q[3];
rz(2.8678558) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.69372988) q[2];
sx q[2];
rz(-1.2749981) q[2];
sx q[2];
rz(-2.0873783) q[2];
rz(2.5888455) q[3];
sx q[3];
rz(-0.25944969) q[3];
sx q[3];
rz(1.118008) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0746821) q[0];
sx q[0];
rz(-0.56147611) q[0];
sx q[0];
rz(3.1232324) q[0];
rz(-1.8439937) q[1];
sx q[1];
rz(-1.3341787) q[1];
sx q[1];
rz(1.1150572) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.499704) q[0];
sx q[0];
rz(-1.5532237) q[0];
sx q[0];
rz(-3.0123424) q[0];
rz(-pi) q[1];
rz(-2.1586559) q[2];
sx q[2];
rz(-0.52861428) q[2];
sx q[2];
rz(-1.7129218) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.12488481) q[1];
sx q[1];
rz(-2.36824) q[1];
sx q[1];
rz(0.70752899) q[1];
rz(-2.6874198) q[3];
sx q[3];
rz(-2.0273367) q[3];
sx q[3];
rz(-0.87107108) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.4516419) q[2];
sx q[2];
rz(-0.74863282) q[2];
sx q[2];
rz(-0.79317036) q[2];
rz(-2.4801109) q[3];
sx q[3];
rz(-1.7659148) q[3];
sx q[3];
rz(0.64600265) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.20872214) q[0];
sx q[0];
rz(-2.0168309) q[0];
sx q[0];
rz(0.17054184) q[0];
rz(1.0199245) q[1];
sx q[1];
rz(-2.7248757) q[1];
sx q[1];
rz(-0.65933093) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.19457283) q[0];
sx q[0];
rz(-1.6259369) q[0];
sx q[0];
rz(1.5558262) q[0];
rz(0.77240981) q[2];
sx q[2];
rz(-0.67467022) q[2];
sx q[2];
rz(-1.2090558) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.7119321) q[1];
sx q[1];
rz(-1.871456) q[1];
sx q[1];
rz(0.64325227) q[1];
rz(2.1390247) q[3];
sx q[3];
rz(-0.89876491) q[3];
sx q[3];
rz(-0.51154414) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.3063804) q[2];
sx q[2];
rz(-2.7472718) q[2];
sx q[2];
rz(-0.96768704) q[2];
rz(-0.59099284) q[3];
sx q[3];
rz(-1.7715745) q[3];
sx q[3];
rz(1.7041357) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10385926) q[0];
sx q[0];
rz(-0.93872207) q[0];
sx q[0];
rz(-0.18873225) q[0];
rz(2.0813148) q[1];
sx q[1];
rz(-0.66245586) q[1];
sx q[1];
rz(2.2417384) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8215733) q[0];
sx q[0];
rz(-2.061141) q[0];
sx q[0];
rz(2.1925395) q[0];
rz(-pi) q[1];
rz(-1.5678039) q[2];
sx q[2];
rz(-0.72751497) q[2];
sx q[2];
rz(-1.8513917) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.65727216) q[1];
sx q[1];
rz(-1.5725699) q[1];
sx q[1];
rz(0.65048762) q[1];
rz(-2.9505355) q[3];
sx q[3];
rz(-2.2828498) q[3];
sx q[3];
rz(-2.4993757) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.4047644) q[2];
sx q[2];
rz(-1.9177723) q[2];
sx q[2];
rz(-1.9723816) q[2];
rz(-2.2660008) q[3];
sx q[3];
rz(-2.4647522) q[3];
sx q[3];
rz(-0.08918795) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(-1.635023) q[0];
sx q[0];
rz(-1.4590141) q[0];
sx q[0];
rz(-0.42612472) q[0];
rz(1.2395073) q[1];
sx q[1];
rz(-1.0303717) q[1];
sx q[1];
rz(0.32136163) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6792435) q[0];
sx q[0];
rz(-2.5658742) q[0];
sx q[0];
rz(-1.4690983) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.05409492) q[2];
sx q[2];
rz(-2.0129497) q[2];
sx q[2];
rz(-1.4070321) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.4453955) q[1];
sx q[1];
rz(-1.2330258) q[1];
sx q[1];
rz(0.15134751) q[1];
rz(-pi) q[2];
rz(-1.0635183) q[3];
sx q[3];
rz(-0.92417704) q[3];
sx q[3];
rz(-2.3648928) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.1539803) q[2];
sx q[2];
rz(-0.93903792) q[2];
sx q[2];
rz(-0.067848094) q[2];
rz(0.96505729) q[3];
sx q[3];
rz(-0.32876757) q[3];
sx q[3];
rz(-1.4062101) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6168552) q[0];
sx q[0];
rz(-2.6829834) q[0];
sx q[0];
rz(0.99880698) q[0];
rz(-2.2899992) q[1];
sx q[1];
rz(-1.0706182) q[1];
sx q[1];
rz(0.65382438) q[1];
rz(0.55642301) q[2];
sx q[2];
rz(-1.6243373) q[2];
sx q[2];
rz(-2.7173964) q[2];
rz(2.4236267) q[3];
sx q[3];
rz(-1.5893275) q[3];
sx q[3];
rz(0.5940819) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
