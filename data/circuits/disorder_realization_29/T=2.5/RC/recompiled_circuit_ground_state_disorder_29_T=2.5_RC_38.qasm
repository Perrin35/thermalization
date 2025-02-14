OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.86201) q[0];
sx q[0];
rz(-0.53505889) q[0];
sx q[0];
rz(0.23393272) q[0];
rz(-2.2995931) q[1];
sx q[1];
rz(-1.41398) q[1];
sx q[1];
rz(-1.6230621) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5748225) q[0];
sx q[0];
rz(-1.260551) q[0];
sx q[0];
rz(-0.53939087) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0056304) q[2];
sx q[2];
rz(-1.3771435) q[2];
sx q[2];
rz(0.47392148) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.7747353) q[1];
sx q[1];
rz(-2.1090057) q[1];
sx q[1];
rz(-0.66267207) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.47187658) q[3];
sx q[3];
rz(-2.1778244) q[3];
sx q[3];
rz(2.8267953) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.5412377) q[2];
sx q[2];
rz(-1.7011832) q[2];
sx q[2];
rz(2.3360628) q[2];
rz(-0.39189288) q[3];
sx q[3];
rz(-1.7854179) q[3];
sx q[3];
rz(2.384757) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.58981744) q[0];
sx q[0];
rz(-2.4601695) q[0];
sx q[0];
rz(-0.43847325) q[0];
rz(1.8665727) q[1];
sx q[1];
rz(-1.6908815) q[1];
sx q[1];
rz(3.0335887) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1718801) q[0];
sx q[0];
rz(-1.4599271) q[0];
sx q[0];
rz(0.025048704) q[0];
rz(2.2453868) q[2];
sx q[2];
rz(-2.2707078) q[2];
sx q[2];
rz(0.43255478) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.17441544) q[1];
sx q[1];
rz(-1.4870054) q[1];
sx q[1];
rz(-1.3166974) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9181973) q[3];
sx q[3];
rz(-2.6250771) q[3];
sx q[3];
rz(-1.7361189) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.7763623) q[2];
sx q[2];
rz(-2.6680816) q[2];
sx q[2];
rz(3.0253809) q[2];
rz(-2.727437) q[3];
sx q[3];
rz(-1.073758) q[3];
sx q[3];
rz(0.80348408) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-0.47495833) q[0];
sx q[0];
rz(-1.3131498) q[0];
sx q[0];
rz(-2.3057002) q[0];
rz(1.6235141) q[1];
sx q[1];
rz(-1.6267136) q[1];
sx q[1];
rz(1.9035043) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9515563) q[0];
sx q[0];
rz(-2.073003) q[0];
sx q[0];
rz(0.42597187) q[0];
rz(-1.8169077) q[2];
sx q[2];
rz(-1.2370584) q[2];
sx q[2];
rz(0.98881236) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.7459384) q[1];
sx q[1];
rz(-1.4534344) q[1];
sx q[1];
rz(0.057988564) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8594871) q[3];
sx q[3];
rz(-1.2938525) q[3];
sx q[3];
rz(-1.9953342) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.8572924) q[2];
sx q[2];
rz(-2.9967873) q[2];
sx q[2];
rz(-0.71883744) q[2];
rz(0.58088628) q[3];
sx q[3];
rz(-1.418117) q[3];
sx q[3];
rz(-2.412793) q[3];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6179287) q[0];
sx q[0];
rz(-1.5476462) q[0];
sx q[0];
rz(1.1267598) q[0];
rz(1.2976546) q[1];
sx q[1];
rz(-1.6512066) q[1];
sx q[1];
rz(-1.0430956) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8422416) q[0];
sx q[0];
rz(-1.3004005) q[0];
sx q[0];
rz(-3.110136) q[0];
rz(-pi) q[1];
rz(0.2506855) q[2];
sx q[2];
rz(-1.0480685) q[2];
sx q[2];
rz(1.2831068) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.9981737) q[1];
sx q[1];
rz(-2.7980248) q[1];
sx q[1];
rz(-2.0433389) q[1];
rz(2.486242) q[3];
sx q[3];
rz(-1.5028302) q[3];
sx q[3];
rz(-2.7790536) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.8635233) q[2];
sx q[2];
rz(-1.7570644) q[2];
sx q[2];
rz(1.0370022) q[2];
rz(0.95748025) q[3];
sx q[3];
rz(-2.0786395) q[3];
sx q[3];
rz(0.83468848) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1400414) q[0];
sx q[0];
rz(-0.20714864) q[0];
sx q[0];
rz(3.0539883) q[0];
rz(0.43830782) q[1];
sx q[1];
rz(-2.0594845) q[1];
sx q[1];
rz(-1.8278488) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8439595) q[0];
sx q[0];
rz(-0.54415138) q[0];
sx q[0];
rz(2.0354969) q[0];
rz(-1.9319973) q[2];
sx q[2];
rz(-0.10538498) q[2];
sx q[2];
rz(-3.1230645) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.6586761) q[1];
sx q[1];
rz(-2.9476142) q[1];
sx q[1];
rz(1.0922171) q[1];
rz(-pi) q[2];
rz(-0.99788061) q[3];
sx q[3];
rz(-0.30227236) q[3];
sx q[3];
rz(2.5469766) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.4801415) q[2];
sx q[2];
rz(-1.2378614) q[2];
sx q[2];
rz(2.5782149) q[2];
rz(-1.2207458) q[3];
sx q[3];
rz(-1.7505587) q[3];
sx q[3];
rz(-2.5105072) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0193943) q[0];
sx q[0];
rz(-2.4833184) q[0];
sx q[0];
rz(0.011938183) q[0];
rz(2.7519233) q[1];
sx q[1];
rz(-2.4395112) q[1];
sx q[1];
rz(-2.8964002) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5545776) q[0];
sx q[0];
rz(-2.3759807) q[0];
sx q[0];
rz(0.9258201) q[0];
x q[1];
rz(1.7952314) q[2];
sx q[2];
rz(-1.8230453) q[2];
sx q[2];
rz(-0.66649619) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.9710412) q[1];
sx q[1];
rz(-1.028113) q[1];
sx q[1];
rz(-2.0849063) q[1];
x q[2];
rz(-2.7390476) q[3];
sx q[3];
rz(-1.8410826) q[3];
sx q[3];
rz(-2.9077934) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.37990722) q[2];
sx q[2];
rz(-1.8319538) q[2];
sx q[2];
rz(1.7725819) q[2];
rz(0.53064972) q[3];
sx q[3];
rz(-0.68415087) q[3];
sx q[3];
rz(2.0556889) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1928007) q[0];
sx q[0];
rz(-1.3230319) q[0];
sx q[0];
rz(-2.8444994) q[0];
rz(-1.6311215) q[1];
sx q[1];
rz(-0.62989569) q[1];
sx q[1];
rz(2.7105601) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4925028) q[0];
sx q[0];
rz(-1.8592872) q[0];
sx q[0];
rz(-1.3516851) q[0];
x q[1];
rz(1.5181584) q[2];
sx q[2];
rz(-1.3402914) q[2];
sx q[2];
rz(0.33256691) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.7073133) q[1];
sx q[1];
rz(-2.523604) q[1];
sx q[1];
rz(2.1531578) q[1];
rz(-pi) q[2];
rz(2.8481315) q[3];
sx q[3];
rz(-2.0803422) q[3];
sx q[3];
rz(-2.5603385) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.83583528) q[2];
sx q[2];
rz(-2.367986) q[2];
sx q[2];
rz(0.26710278) q[2];
rz(-3.109572) q[3];
sx q[3];
rz(-1.1705541) q[3];
sx q[3];
rz(-0.10572461) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0778462) q[0];
sx q[0];
rz(-1.2910605) q[0];
sx q[0];
rz(0.63999501) q[0];
rz(0.85482875) q[1];
sx q[1];
rz(-1.4850441) q[1];
sx q[1];
rz(-1.4124426) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0873358) q[0];
sx q[0];
rz(-1.0414413) q[0];
sx q[0];
rz(-2.0482778) q[0];
x q[1];
rz(-3.090254) q[2];
sx q[2];
rz(-1.079139) q[2];
sx q[2];
rz(-1.3528388) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.1302274) q[1];
sx q[1];
rz(-1.2308569) q[1];
sx q[1];
rz(1.66015) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.85573393) q[3];
sx q[3];
rz(-2.6011356) q[3];
sx q[3];
rz(-0.021322535) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.60014805) q[2];
sx q[2];
rz(-1.4358127) q[2];
sx q[2];
rz(2.7195462) q[2];
rz(0.47354928) q[3];
sx q[3];
rz(-2.519042) q[3];
sx q[3];
rz(2.9464909) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37050978) q[0];
sx q[0];
rz(-2.2042553) q[0];
sx q[0];
rz(-1.3516082) q[0];
rz(-0.31556684) q[1];
sx q[1];
rz(-1.4548929) q[1];
sx q[1];
rz(1.1531166) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2446072) q[0];
sx q[0];
rz(-3.0685406) q[0];
sx q[0];
rz(2.3417419) q[0];
rz(-pi) q[1];
x q[1];
rz(0.081772371) q[2];
sx q[2];
rz(-0.53682971) q[2];
sx q[2];
rz(-0.54807907) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.35495423) q[1];
sx q[1];
rz(-1.6980722) q[1];
sx q[1];
rz(-0.58357012) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.38001506) q[3];
sx q[3];
rz(-1.1855864) q[3];
sx q[3];
rz(-0.64982254) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.6844668) q[2];
sx q[2];
rz(-1.212333) q[2];
sx q[2];
rz(-2.7678164) q[2];
rz(-1.2260381) q[3];
sx q[3];
rz(-2.2235179) q[3];
sx q[3];
rz(1.8019684) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.95947295) q[0];
sx q[0];
rz(-2.4612893) q[0];
sx q[0];
rz(-1.8208338) q[0];
rz(-0.93718115) q[1];
sx q[1];
rz(-1.8056185) q[1];
sx q[1];
rz(2.6722867) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26312882) q[0];
sx q[0];
rz(-1.4944021) q[0];
sx q[0];
rz(1.4532277) q[0];
rz(-pi) q[1];
rz(1.2861757) q[2];
sx q[2];
rz(-2.5409343) q[2];
sx q[2];
rz(0.17135581) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.7368245) q[1];
sx q[1];
rz(-1.4206845) q[1];
sx q[1];
rz(1.4895358) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8402249) q[3];
sx q[3];
rz(-1.8068988) q[3];
sx q[3];
rz(2.2599041) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.70539537) q[2];
sx q[2];
rz(-0.86279482) q[2];
sx q[2];
rz(1.2104642) q[2];
rz(-0.56810275) q[3];
sx q[3];
rz(-1.7567239) q[3];
sx q[3];
rz(-0.97584045) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2198915) q[0];
sx q[0];
rz(-0.85048631) q[0];
sx q[0];
rz(1.462107) q[0];
rz(0.89314356) q[1];
sx q[1];
rz(-1.8291263) q[1];
sx q[1];
rz(-1.6222454) q[1];
rz(-0.90441119) q[2];
sx q[2];
rz(-0.36199649) q[2];
sx q[2];
rz(0.71142759) q[2];
rz(2.1410971) q[3];
sx q[3];
rz(-1.199493) q[3];
sx q[3];
rz(-1.5293157) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
