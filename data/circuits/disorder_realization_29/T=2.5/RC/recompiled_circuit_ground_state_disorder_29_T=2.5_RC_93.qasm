OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.27958265) q[0];
sx q[0];
rz(-2.6065338) q[0];
sx q[0];
rz(2.9076599) q[0];
rz(0.84199953) q[1];
sx q[1];
rz(4.5555727) q[1];
sx q[1];
rz(11.04784) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5667701) q[0];
sx q[0];
rz(-1.260551) q[0];
sx q[0];
rz(-0.53939087) q[0];
rz(-2.0056304) q[2];
sx q[2];
rz(-1.7644492) q[2];
sx q[2];
rz(2.6676712) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.3564811) q[1];
sx q[1];
rz(-2.3143396) q[1];
sx q[1];
rz(0.77036304) q[1];
rz(0.99125864) q[3];
sx q[3];
rz(-0.75014866) q[3];
sx q[3];
rz(2.0969491) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.5412377) q[2];
sx q[2];
rz(-1.4404094) q[2];
sx q[2];
rz(0.8055299) q[2];
rz(-0.39189288) q[3];
sx q[3];
rz(-1.3561748) q[3];
sx q[3];
rz(-2.384757) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.58981744) q[0];
sx q[0];
rz(-0.68142319) q[0];
sx q[0];
rz(-0.43847325) q[0];
rz(-1.8665727) q[1];
sx q[1];
rz(-1.6908815) q[1];
sx q[1];
rz(-3.0335887) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9697125) q[0];
sx q[0];
rz(-1.6816655) q[0];
sx q[0];
rz(0.025048704) q[0];
rz(-0.82307016) q[2];
sx q[2];
rz(-1.0727173) q[2];
sx q[2];
rz(0.66253875) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.4181136) q[1];
sx q[1];
rz(-1.3176085) q[1];
sx q[1];
rz(-0.0865571) q[1];
x q[2];
rz(-2.9505902) q[3];
sx q[3];
rz(-1.0878948) q[3];
sx q[3];
rz(-1.0109028) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.7763623) q[2];
sx q[2];
rz(-2.6680816) q[2];
sx q[2];
rz(3.0253809) q[2];
rz(2.727437) q[3];
sx q[3];
rz(-2.0678346) q[3];
sx q[3];
rz(-2.3381086) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6666343) q[0];
sx q[0];
rz(-1.8284429) q[0];
sx q[0];
rz(0.83589244) q[0];
rz(1.6235141) q[1];
sx q[1];
rz(-1.6267136) q[1];
sx q[1];
rz(-1.2380884) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.19003632) q[0];
sx q[0];
rz(-1.0685896) q[0];
sx q[0];
rz(-0.42597187) q[0];
x q[1];
rz(-1.8169077) q[2];
sx q[2];
rz(-1.9045343) q[2];
sx q[2];
rz(-0.98881236) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.9596531) q[1];
sx q[1];
rz(-1.6283855) q[1];
sx q[1];
rz(-1.4532386) q[1];
rz(-pi) q[2];
rz(-1.2830623) q[3];
sx q[3];
rz(-1.2997174) q[3];
sx q[3];
rz(-2.6379741) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.8572924) q[2];
sx q[2];
rz(-2.9967873) q[2];
sx q[2];
rz(-2.4227552) q[2];
rz(-0.58088628) q[3];
sx q[3];
rz(-1.7234756) q[3];
sx q[3];
rz(-2.412793) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6179287) q[0];
sx q[0];
rz(-1.5476462) q[0];
sx q[0];
rz(-2.0148328) q[0];
rz(1.2976546) q[1];
sx q[1];
rz(-1.490386) q[1];
sx q[1];
rz(1.0430956) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.299351) q[0];
sx q[0];
rz(-1.8411921) q[0];
sx q[0];
rz(0.031456703) q[0];
rz(-pi) q[1];
x q[1];
rz(0.2506855) q[2];
sx q[2];
rz(-2.0935241) q[2];
sx q[2];
rz(-1.2831068) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.14341893) q[1];
sx q[1];
rz(-0.34356782) q[1];
sx q[1];
rz(-2.0433389) q[1];
rz(-3.0303589) q[3];
sx q[3];
rz(-2.4832442) q[3];
sx q[3];
rz(1.8452132) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.2780693) q[2];
sx q[2];
rz(-1.3845283) q[2];
sx q[2];
rz(1.0370022) q[2];
rz(-2.1841124) q[3];
sx q[3];
rz(-1.0629531) q[3];
sx q[3];
rz(-0.83468848) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1400414) q[0];
sx q[0];
rz(-2.934444) q[0];
sx q[0];
rz(3.0539883) q[0];
rz(0.43830782) q[1];
sx q[1];
rz(-2.0594845) q[1];
sx q[1];
rz(-1.8278488) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2976332) q[0];
sx q[0];
rz(-0.54415138) q[0];
sx q[0];
rz(2.0354969) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2095954) q[2];
sx q[2];
rz(-3.0362077) q[2];
sx q[2];
rz(3.1230645) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.96927724) q[1];
sx q[1];
rz(-1.7427497) q[1];
sx q[1];
rz(0.090222619) q[1];
rz(-pi) q[2];
x q[2];
rz(0.99788061) q[3];
sx q[3];
rz(-0.30227236) q[3];
sx q[3];
rz(-2.5469766) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.6614512) q[2];
sx q[2];
rz(-1.2378614) q[2];
sx q[2];
rz(0.56337774) q[2];
rz(1.2207458) q[3];
sx q[3];
rz(-1.7505587) q[3];
sx q[3];
rz(2.5105072) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0193943) q[0];
sx q[0];
rz(-0.65827426) q[0];
sx q[0];
rz(-3.1296545) q[0];
rz(0.38966933) q[1];
sx q[1];
rz(-0.7020815) q[1];
sx q[1];
rz(0.24519244) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6547873) q[0];
sx q[0];
rz(-1.1410895) q[0];
sx q[0];
rz(-0.91581099) q[0];
rz(-pi) q[1];
rz(-0.25845627) q[2];
sx q[2];
rz(-1.788013) q[2];
sx q[2];
rz(0.9612135) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.4575544) q[1];
sx q[1];
rz(-1.1361309) q[1];
sx q[1];
rz(2.5358389) q[1];
rz(-2.526029) q[3];
sx q[3];
rz(-0.48071733) q[3];
sx q[3];
rz(-1.2445039) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.7616854) q[2];
sx q[2];
rz(-1.8319538) q[2];
sx q[2];
rz(-1.3690108) q[2];
rz(-0.53064972) q[3];
sx q[3];
rz(-0.68415087) q[3];
sx q[3];
rz(1.0859038) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.948792) q[0];
sx q[0];
rz(-1.8185607) q[0];
sx q[0];
rz(-2.8444994) q[0];
rz(1.6311215) q[1];
sx q[1];
rz(-2.511697) q[1];
sx q[1];
rz(-0.43103257) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4925028) q[0];
sx q[0];
rz(-1.2823055) q[0];
sx q[0];
rz(-1.7899075) q[0];
x q[1];
rz(-2.9107793) q[2];
sx q[2];
rz(-1.5195519) q[2];
sx q[2];
rz(1.9153999) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.7707211) q[1];
sx q[1];
rz(-1.2464735) q[1];
sx q[1];
rz(1.0350219) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0484974) q[3];
sx q[3];
rz(-2.5601031) q[3];
sx q[3];
rz(0.026611004) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3057574) q[2];
sx q[2];
rz(-0.77360669) q[2];
sx q[2];
rz(-2.8744899) q[2];
rz(0.032020656) q[3];
sx q[3];
rz(-1.1705541) q[3];
sx q[3];
rz(3.035868) q[3];
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
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.06374643) q[0];
sx q[0];
rz(-1.8505322) q[0];
sx q[0];
rz(-0.63999501) q[0];
rz(-0.85482875) q[1];
sx q[1];
rz(-1.4850441) q[1];
sx q[1];
rz(-1.7291501) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8518191) q[0];
sx q[0];
rz(-0.69727325) q[0];
sx q[0];
rz(-2.4757705) q[0];
rz(-pi) q[1];
rz(-1.0785901) q[2];
sx q[2];
rz(-1.5255431) q[2];
sx q[2];
rz(2.947888) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.392726) q[1];
sx q[1];
rz(-2.7905474) q[1];
sx q[1];
rz(2.8944394) q[1];
rz(-pi) q[2];
rz(2.7667609) q[3];
sx q[3];
rz(-1.171798) q[3];
sx q[3];
rz(2.3712036) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.60014805) q[2];
sx q[2];
rz(-1.7057799) q[2];
sx q[2];
rz(2.7195462) q[2];
rz(2.6680434) q[3];
sx q[3];
rz(-2.519042) q[3];
sx q[3];
rz(-2.9464909) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7710829) q[0];
sx q[0];
rz(-0.93733731) q[0];
sx q[0];
rz(-1.3516082) q[0];
rz(-0.31556684) q[1];
sx q[1];
rz(-1.4548929) q[1];
sx q[1];
rz(-1.9884761) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0457927) q[0];
sx q[0];
rz(-1.6216767) q[0];
sx q[0];
rz(1.6232383) q[0];
x q[1];
rz(-3.0598203) q[2];
sx q[2];
rz(-0.53682971) q[2];
sx q[2];
rz(2.5935136) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.8421346) q[1];
sx q[1];
rz(-2.149035) q[1];
sx q[1];
rz(-1.4186354) q[1];
x q[2];
rz(1.159129) q[3];
sx q[3];
rz(-1.921706) q[3];
sx q[3];
rz(1.0699492) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.6844668) q[2];
sx q[2];
rz(-1.9292597) q[2];
sx q[2];
rz(-0.37377629) q[2];
rz(1.9155546) q[3];
sx q[3];
rz(-0.91807476) q[3];
sx q[3];
rz(1.3396243) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1821197) q[0];
sx q[0];
rz(-2.4612893) q[0];
sx q[0];
rz(-1.3207588) q[0];
rz(0.93718115) q[1];
sx q[1];
rz(-1.8056185) q[1];
sx q[1];
rz(-2.6722867) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.73399765) q[0];
sx q[0];
rz(-3.0014801) q[0];
sx q[0];
rz(2.1489643) q[0];
rz(-pi) q[1];
rz(0.19005084) q[2];
sx q[2];
rz(-0.99744883) q[2];
sx q[2];
rz(0.1694451) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.1782068) q[1];
sx q[1];
rz(-1.4904516) q[1];
sx q[1];
rz(-2.9909913) q[1];
x q[2];
rz(-0.68113459) q[3];
sx q[3];
rz(-0.38060846) q[3];
sx q[3];
rz(1.8073624) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.4361973) q[2];
sx q[2];
rz(-0.86279482) q[2];
sx q[2];
rz(1.2104642) q[2];
rz(0.56810275) q[3];
sx q[3];
rz(-1.7567239) q[3];
sx q[3];
rz(-2.1657522) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9217011) q[0];
sx q[0];
rz(-0.85048631) q[0];
sx q[0];
rz(1.462107) q[0];
rz(2.2484491) q[1];
sx q[1];
rz(-1.3124663) q[1];
sx q[1];
rz(1.5193473) q[1];
rz(1.2814796) q[2];
sx q[2];
rz(-1.791496) q[2];
sx q[2];
rz(1.6481177) q[2];
rz(2.1956069) q[3];
sx q[3];
rz(-0.66910268) q[3];
sx q[3];
rz(0.55629081) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
