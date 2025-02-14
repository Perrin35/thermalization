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
rz(5.7481264) q[0];
sx q[0];
rz(9.6587107) q[0];
rz(0.84199953) q[1];
sx q[1];
rz(4.5555727) q[1];
sx q[1];
rz(11.04784) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.18478014) q[0];
sx q[0];
rz(-2.0818458) q[0];
sx q[0];
rz(1.9283811) q[0];
rz(-pi) q[1];
rz(-1.1351003) q[2];
sx q[2];
rz(-2.6681136) q[2];
sx q[2];
rz(1.651929) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.9649992) q[1];
sx q[1];
rz(-1.0142583) q[1];
sx q[1];
rz(-2.2189369) q[1];
rz(-pi) q[2];
rz(-2.6697161) q[3];
sx q[3];
rz(-2.1778244) q[3];
sx q[3];
rz(0.31479731) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.6003549) q[2];
sx q[2];
rz(-1.4404094) q[2];
sx q[2];
rz(0.8055299) q[2];
rz(0.39189288) q[3];
sx q[3];
rz(-1.3561748) q[3];
sx q[3];
rz(-0.75683561) q[3];
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
rz(-pi/2) q[0];
x q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5517752) q[0];
sx q[0];
rz(-0.68142319) q[0];
sx q[0];
rz(2.7031194) q[0];
rz(-1.8665727) q[1];
sx q[1];
rz(-1.4507111) q[1];
sx q[1];
rz(-0.10800392) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.94919449) q[0];
sx q[0];
rz(-3.0279403) q[0];
sx q[0];
rz(1.7920919) q[0];
rz(-2.503452) q[2];
sx q[2];
rz(-2.2109988) q[2];
sx q[2];
rz(-1.3253044) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.17441544) q[1];
sx q[1];
rz(-1.6545873) q[1];
sx q[1];
rz(-1.3166974) q[1];
rz(-pi) q[2];
rz(-1.0803079) q[3];
sx q[3];
rz(-1.401859) q[3];
sx q[3];
rz(-0.47034697) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.3652304) q[2];
sx q[2];
rz(-0.47351101) q[2];
sx q[2];
rz(-3.0253809) q[2];
rz(-2.727437) q[3];
sx q[3];
rz(-1.073758) q[3];
sx q[3];
rz(-2.3381086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47495833) q[0];
sx q[0];
rz(-1.8284429) q[0];
sx q[0];
rz(0.83589244) q[0];
rz(1.5180786) q[1];
sx q[1];
rz(-1.6267136) q[1];
sx q[1];
rz(1.2380884) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.19003632) q[0];
sx q[0];
rz(-2.073003) q[0];
sx q[0];
rz(2.7156208) q[0];
rz(-pi) q[1];
rz(-1.8169077) q[2];
sx q[2];
rz(-1.2370584) q[2];
sx q[2];
rz(-2.1527803) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.39565428) q[1];
sx q[1];
rz(-1.6881583) q[1];
sx q[1];
rz(-3.0836041) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8594871) q[3];
sx q[3];
rz(-1.8477401) q[3];
sx q[3];
rz(-1.9953342) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.28430024) q[2];
sx q[2];
rz(-0.14480536) q[2];
sx q[2];
rz(-2.4227552) q[2];
rz(0.58088628) q[3];
sx q[3];
rz(-1.7234756) q[3];
sx q[3];
rz(2.412793) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
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
rz(-pi) q[1];
x q[1];
sx q[2];
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
rz(2.0984971) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.299351) q[0];
sx q[0];
rz(-1.8411921) q[0];
sx q[0];
rz(-0.031456703) q[0];
x q[1];
rz(1.1642493) q[2];
sx q[2];
rz(-2.5669328) q[2];
sx q[2];
rz(2.3324147) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.9981737) q[1];
sx q[1];
rz(-2.7980248) q[1];
sx q[1];
rz(1.0982537) q[1];
x q[2];
rz(-0.11123379) q[3];
sx q[3];
rz(-2.4832442) q[3];
sx q[3];
rz(1.2963795) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.8635233) q[2];
sx q[2];
rz(-1.3845283) q[2];
sx q[2];
rz(1.0370022) q[2];
rz(-0.95748025) q[3];
sx q[3];
rz(-2.0786395) q[3];
sx q[3];
rz(2.3069042) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
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
rz(1.3137438) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4632516) q[0];
sx q[0];
rz(-1.3366564) q[0];
sx q[0];
rz(2.0666422) q[0];
x q[1];
rz(1.2095954) q[2];
sx q[2];
rz(-3.0362077) q[2];
sx q[2];
rz(-0.018528136) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.524595) q[1];
sx q[1];
rz(-1.6596848) q[1];
sx q[1];
rz(1.3981546) q[1];
x q[2];
rz(-2.143712) q[3];
sx q[3];
rz(-2.8393203) q[3];
sx q[3];
rz(-0.59461601) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.6614512) q[2];
sx q[2];
rz(-1.2378614) q[2];
sx q[2];
rz(2.5782149) q[2];
rz(1.9208469) q[3];
sx q[3];
rz(-1.7505587) q[3];
sx q[3];
rz(-2.5105072) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0193943) q[0];
sx q[0];
rz(-0.65827426) q[0];
sx q[0];
rz(3.1296545) q[0];
rz(-2.7519233) q[1];
sx q[1];
rz(-0.7020815) q[1];
sx q[1];
rz(0.24519244) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48680533) q[0];
sx q[0];
rz(-2.0005032) q[0];
sx q[0];
rz(2.2257817) q[0];
x q[1];
rz(1.3463613) q[2];
sx q[2];
rz(-1.8230453) q[2];
sx q[2];
rz(0.66649619) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.17055146) q[1];
sx q[1];
rz(-2.1134796) q[1];
sx q[1];
rz(-2.0849063) q[1];
rz(-pi) q[2];
rz(-2.7390476) q[3];
sx q[3];
rz(-1.8410826) q[3];
sx q[3];
rz(0.2337993) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7616854) q[2];
sx q[2];
rz(-1.8319538) q[2];
sx q[2];
rz(1.7725819) q[2];
rz(-2.6109429) q[3];
sx q[3];
rz(-2.4574418) q[3];
sx q[3];
rz(1.0859038) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.948792) q[0];
sx q[0];
rz(-1.3230319) q[0];
sx q[0];
rz(0.2970933) q[0];
rz(-1.6311215) q[1];
sx q[1];
rz(-0.62989569) q[1];
sx q[1];
rz(2.7105601) q[1];
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
rz(-pi) q[1];
rz(-1.6234342) q[2];
sx q[2];
rz(-1.8013012) q[2];
sx q[2];
rz(2.8090257) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.4342793) q[1];
sx q[1];
rz(-0.61798862) q[1];
sx q[1];
rz(0.98843482) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0424352) q[3];
sx q[3];
rz(-1.3155121) q[3];
sx q[3];
rz(-2.0057037) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.3057574) q[2];
sx q[2];
rz(-0.77360669) q[2];
sx q[2];
rz(-2.8744899) q[2];
rz(-0.032020656) q[3];
sx q[3];
rz(-1.9710385) q[3];
sx q[3];
rz(-0.10572461) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.06374643) q[0];
sx q[0];
rz(-1.8505322) q[0];
sx q[0];
rz(2.5015976) q[0];
rz(-0.85482875) q[1];
sx q[1];
rz(-1.4850441) q[1];
sx q[1];
rz(1.4124426) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0542568) q[0];
sx q[0];
rz(-1.0414413) q[0];
sx q[0];
rz(-2.0482778) q[0];
rz(-1.0785901) q[2];
sx q[2];
rz(-1.5255431) q[2];
sx q[2];
rz(-0.19370463) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.58929491) q[1];
sx q[1];
rz(-1.4865685) q[1];
sx q[1];
rz(0.34119795) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.37483172) q[3];
sx q[3];
rz(-1.9697947) q[3];
sx q[3];
rz(0.77038902) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.5414446) q[2];
sx q[2];
rz(-1.4358127) q[2];
sx q[2];
rz(0.42204648) q[2];
rz(0.47354928) q[3];
sx q[3];
rz(-0.62255064) q[3];
sx q[3];
rz(0.1951018) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37050978) q[0];
sx q[0];
rz(-2.2042553) q[0];
sx q[0];
rz(1.7899845) q[0];
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
rz(0.095799965) q[0];
sx q[0];
rz(-1.519916) q[0];
sx q[0];
rz(1.5183543) q[0];
rz(-pi) q[1];
rz(2.6062327) q[2];
sx q[2];
rz(-1.5290112) q[2];
sx q[2];
rz(2.1891862) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.2994581) q[1];
sx q[1];
rz(-0.99255764) q[1];
sx q[1];
rz(-1.7229573) q[1];
x q[2];
rz(1.159129) q[3];
sx q[3];
rz(-1.921706) q[3];
sx q[3];
rz(1.0699492) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.45712581) q[2];
sx q[2];
rz(-1.9292597) q[2];
sx q[2];
rz(0.37377629) q[2];
rz(-1.2260381) q[3];
sx q[3];
rz(-2.2235179) q[3];
sx q[3];
rz(-1.3396243) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
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
rz(0.95947295) q[0];
sx q[0];
rz(-2.4612893) q[0];
sx q[0];
rz(-1.8208338) q[0];
rz(0.93718115) q[1];
sx q[1];
rz(-1.3359741) q[1];
sx q[1];
rz(-0.46930596) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.824911) q[0];
sx q[0];
rz(-1.6880205) q[0];
sx q[0];
rz(-0.076923142) q[0];
rz(-0.19005084) q[2];
sx q[2];
rz(-2.1441438) q[2];
sx q[2];
rz(0.1694451) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(3.0476956) q[1];
sx q[1];
rz(-2.9710431) q[1];
sx q[1];
rz(-2.6490414) q[1];
x q[2];
rz(0.30136775) q[3];
sx q[3];
rz(-1.8068988) q[3];
sx q[3];
rz(-2.2599041) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.70539537) q[2];
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
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9217011) q[0];
sx q[0];
rz(-2.2911063) q[0];
sx q[0];
rz(-1.6794857) q[0];
rz(-0.89314356) q[1];
sx q[1];
rz(-1.3124663) q[1];
sx q[1];
rz(1.5193473) q[1];
rz(0.90441119) q[2];
sx q[2];
rz(-2.7795962) q[2];
sx q[2];
rz(-2.4301651) q[2];
rz(-0.4332581) q[3];
sx q[3];
rz(-1.0436202) q[3];
sx q[3];
rz(2.9544261) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
