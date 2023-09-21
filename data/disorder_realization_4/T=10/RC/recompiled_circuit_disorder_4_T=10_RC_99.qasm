OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.4274347) q[0];
sx q[0];
rz(2.7058869) q[0];
sx q[0];
rz(11.640179) q[0];
rz(-1.1821094) q[1];
sx q[1];
rz(3.8745772) q[1];
sx q[1];
rz(12.193845) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6719138) q[0];
sx q[0];
rz(-2.1436474) q[0];
sx q[0];
rz(1.8125305) q[0];
x q[1];
rz(2.9973642) q[2];
sx q[2];
rz(-2.5055474) q[2];
sx q[2];
rz(2.6682105) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.8718611) q[1];
sx q[1];
rz(-1.6850123) q[1];
sx q[1];
rz(-1.3807339) q[1];
rz(-pi) q[2];
rz(-0.11301179) q[3];
sx q[3];
rz(-1.7859308) q[3];
sx q[3];
rz(-0.61258951) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.42840502) q[2];
sx q[2];
rz(-1.5593854) q[2];
sx q[2];
rz(2.2170846) q[2];
rz(-1.472578) q[3];
sx q[3];
rz(-1.893483) q[3];
sx q[3];
rz(-1.4991466) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.18838841) q[0];
sx q[0];
rz(-3.0792455) q[0];
sx q[0];
rz(1.6054608) q[0];
rz(-2.9470782) q[1];
sx q[1];
rz(-1.3214) q[1];
sx q[1];
rz(0.054873437) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0458192) q[0];
sx q[0];
rz(-1.4574483) q[0];
sx q[0];
rz(-0.093009526) q[0];
rz(-pi) q[1];
rz(2.4404281) q[2];
sx q[2];
rz(-2.3904281) q[2];
sx q[2];
rz(-0.0093815087) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.730719) q[1];
sx q[1];
rz(-1.5748071) q[1];
sx q[1];
rz(1.9133168) q[1];
x q[2];
rz(-2.8781434) q[3];
sx q[3];
rz(-0.69034319) q[3];
sx q[3];
rz(-0.024397959) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.43869552) q[2];
sx q[2];
rz(-2.5800811) q[2];
sx q[2];
rz(-1.6595718) q[2];
rz(0.99059087) q[3];
sx q[3];
rz(-1.4086658) q[3];
sx q[3];
rz(1.7747442) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7822587) q[0];
sx q[0];
rz(-0.5193091) q[0];
sx q[0];
rz(-0.3749795) q[0];
rz(2.8938876) q[1];
sx q[1];
rz(-1.9539555) q[1];
sx q[1];
rz(-1.8050271) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.61993116) q[0];
sx q[0];
rz(-1.3494028) q[0];
sx q[0];
rz(0.059376052) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5252684) q[2];
sx q[2];
rz(-0.68683544) q[2];
sx q[2];
rz(0.64885215) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.3805483) q[1];
sx q[1];
rz(-0.79954631) q[1];
sx q[1];
rz(1.9995081) q[1];
rz(-pi) q[2];
rz(-2.7258354) q[3];
sx q[3];
rz(-0.81167479) q[3];
sx q[3];
rz(0.31879253) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.1027362) q[2];
sx q[2];
rz(-2.3114624) q[2];
sx q[2];
rz(-2.1170199) q[2];
rz(1.864795) q[3];
sx q[3];
rz(-1.5921311) q[3];
sx q[3];
rz(-1.5156486) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17722002) q[0];
sx q[0];
rz(-1.6765046) q[0];
sx q[0];
rz(-3.1080416) q[0];
rz(1.230348) q[1];
sx q[1];
rz(-0.76554811) q[1];
sx q[1];
rz(-1.4039325) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0008154) q[0];
sx q[0];
rz(-0.40110943) q[0];
sx q[0];
rz(-2.6252069) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3061182) q[2];
sx q[2];
rz(-0.88047853) q[2];
sx q[2];
rz(-0.26963216) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.2851706) q[1];
sx q[1];
rz(-0.31177786) q[1];
sx q[1];
rz(-1.9429893) q[1];
rz(-pi) q[2];
rz(0.085539354) q[3];
sx q[3];
rz(-0.88118689) q[3];
sx q[3];
rz(0.52348189) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.6491062) q[2];
sx q[2];
rz(-2.2968569) q[2];
sx q[2];
rz(2.9283294) q[2];
rz(-0.02031859) q[3];
sx q[3];
rz(-1.820194) q[3];
sx q[3];
rz(0.4030574) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3605109) q[0];
sx q[0];
rz(-1.1163982) q[0];
sx q[0];
rz(-2.1257341) q[0];
rz(-1.6332743) q[1];
sx q[1];
rz(-0.95817482) q[1];
sx q[1];
rz(-3.1075081) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7406697) q[0];
sx q[0];
rz(-1.4105083) q[0];
sx q[0];
rz(2.7435061) q[0];
rz(0.57406117) q[2];
sx q[2];
rz(-1.7957557) q[2];
sx q[2];
rz(0.056919295) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.9220229) q[1];
sx q[1];
rz(-1.2864188) q[1];
sx q[1];
rz(3.1344096) q[1];
rz(-pi) q[2];
rz(2.6328153) q[3];
sx q[3];
rz(-2.1366589) q[3];
sx q[3];
rz(2.1600427) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.74026996) q[2];
sx q[2];
rz(-0.57381845) q[2];
sx q[2];
rz(-1.1711228) q[2];
rz(-1.8088388) q[3];
sx q[3];
rz(-1.7141902) q[3];
sx q[3];
rz(0.12935054) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68833441) q[0];
sx q[0];
rz(-1.1880705) q[0];
sx q[0];
rz(2.7057498) q[0];
rz(-2.0360937) q[1];
sx q[1];
rz(-1.2970129) q[1];
sx q[1];
rz(1.9357392) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5642731) q[0];
sx q[0];
rz(-2.3015129) q[0];
sx q[0];
rz(-2.8217836) q[0];
rz(-pi) q[1];
rz(-0.016024307) q[2];
sx q[2];
rz(-1.6659701) q[2];
sx q[2];
rz(1.8297838) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.6916794) q[1];
sx q[1];
rz(-2.0813585) q[1];
sx q[1];
rz(1.540586) q[1];
rz(-1.4924963) q[3];
sx q[3];
rz(-2.7408528) q[3];
sx q[3];
rz(-2.0164255) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.2445406) q[2];
sx q[2];
rz(-0.62770939) q[2];
sx q[2];
rz(-0.051606027) q[2];
rz(-2.7339325) q[3];
sx q[3];
rz(-1.9577273) q[3];
sx q[3];
rz(0.44529644) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5200941) q[0];
sx q[0];
rz(-2.5243653) q[0];
sx q[0];
rz(0.67333418) q[0];
rz(-2.3576221) q[1];
sx q[1];
rz(-1.7242804) q[1];
sx q[1];
rz(0.53692445) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4250701) q[0];
sx q[0];
rz(-1.636583) q[0];
sx q[0];
rz(-1.487088) q[0];
rz(-pi) q[1];
x q[1];
rz(1.827048) q[2];
sx q[2];
rz(-1.7125687) q[2];
sx q[2];
rz(1.2830551) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.010579022) q[1];
sx q[1];
rz(-0.77242935) q[1];
sx q[1];
rz(-1.7060682) q[1];
rz(-pi) q[2];
rz(1.4173996) q[3];
sx q[3];
rz(-2.1563357) q[3];
sx q[3];
rz(-3.1068387) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.9817104) q[2];
sx q[2];
rz(-1.2572224) q[2];
sx q[2];
rz(2.9505777) q[2];
rz(-0.28132176) q[3];
sx q[3];
rz(-1.9502935) q[3];
sx q[3];
rz(-0.34240001) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1087588) q[0];
sx q[0];
rz(-2.7712951) q[0];
sx q[0];
rz(2.7539745) q[0];
rz(-0.11501137) q[1];
sx q[1];
rz(-1.4287881) q[1];
sx q[1];
rz(0.33755916) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.176929) q[0];
sx q[0];
rz(-1.7690072) q[0];
sx q[0];
rz(-2.6508509) q[0];
rz(-pi) q[1];
rz(-0.93162025) q[2];
sx q[2];
rz(-1.5134303) q[2];
sx q[2];
rz(-0.66040874) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.74886125) q[1];
sx q[1];
rz(-1.0034605) q[1];
sx q[1];
rz(0.15565236) q[1];
rz(1.5070595) q[3];
sx q[3];
rz(-0.81634854) q[3];
sx q[3];
rz(0.07894978) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.18743029) q[2];
sx q[2];
rz(-0.53415161) q[2];
sx q[2];
rz(1.2672651) q[2];
rz(2.3665442) q[3];
sx q[3];
rz(-1.5379484) q[3];
sx q[3];
rz(0.78139853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9545814) q[0];
sx q[0];
rz(-1.0382074) q[0];
sx q[0];
rz(0.22098456) q[0];
rz(2.2194608) q[1];
sx q[1];
rz(-1.8716967) q[1];
sx q[1];
rz(2.1386713) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10119469) q[0];
sx q[0];
rz(-1.4324491) q[0];
sx q[0];
rz(0.25427108) q[0];
x q[1];
rz(2.6384764) q[2];
sx q[2];
rz(-2.0965577) q[2];
sx q[2];
rz(-2.2803277) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.3468614) q[1];
sx q[1];
rz(-2.0309629) q[1];
sx q[1];
rz(2.5882583) q[1];
x q[2];
rz(-2.2569611) q[3];
sx q[3];
rz(-1.7836708) q[3];
sx q[3];
rz(0.5205982) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.4724491) q[2];
sx q[2];
rz(-2.3995235) q[2];
sx q[2];
rz(0.15850244) q[2];
rz(0.45378271) q[3];
sx q[3];
rz(-2.3588389) q[3];
sx q[3];
rz(-2.2973072) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6145988) q[0];
sx q[0];
rz(-0.90899962) q[0];
sx q[0];
rz(2.9504543) q[0];
rz(2.846431) q[1];
sx q[1];
rz(-2.252153) q[1];
sx q[1];
rz(0.89231649) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44360456) q[0];
sx q[0];
rz(-2.3058878) q[0];
sx q[0];
rz(2.9025335) q[0];
rz(2.2116682) q[2];
sx q[2];
rz(-1.495549) q[2];
sx q[2];
rz(1.6844695) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.72496966) q[1];
sx q[1];
rz(-1.6367216) q[1];
sx q[1];
rz(0.23857393) q[1];
rz(2.7097706) q[3];
sx q[3];
rz(-1.338827) q[3];
sx q[3];
rz(0.98917978) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.7717379) q[2];
sx q[2];
rz(-2.3042046) q[2];
sx q[2];
rz(2.2820293) q[2];
rz(-1.2236979) q[3];
sx q[3];
rz(-0.92646354) q[3];
sx q[3];
rz(-2.3971476) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
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
rz(1.1214462) q[0];
sx q[0];
rz(-0.93562026) q[0];
sx q[0];
rz(1.390441) q[0];
rz(1.6620811) q[1];
sx q[1];
rz(-2.6585487) q[1];
sx q[1];
rz(-1.9225635) q[1];
rz(3.1142278) q[2];
sx q[2];
rz(-1.27925) q[2];
sx q[2];
rz(-1.7640511) q[2];
rz(-1.6395232) q[3];
sx q[3];
rz(-1.5309661) q[3];
sx q[3];
rz(-1.0504709) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
