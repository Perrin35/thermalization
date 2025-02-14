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
rz(-1.4827363) q[0];
sx q[0];
rz(-2.1602477) q[0];
sx q[0];
rz(2.0438097) q[0];
rz(-0.78805796) q[1];
sx q[1];
rz(-2.0907953) q[1];
sx q[1];
rz(-0.0069590574) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1084676) q[0];
sx q[0];
rz(-1.5146717) q[0];
sx q[0];
rz(2.2214487) q[0];
x q[1];
rz(-0.9240146) q[2];
sx q[2];
rz(-1.7597464) q[2];
sx q[2];
rz(-1.762977) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.6485868) q[1];
sx q[1];
rz(-1.8064335) q[1];
sx q[1];
rz(-0.034502397) q[1];
rz(-0.29421774) q[3];
sx q[3];
rz(-2.7876884) q[3];
sx q[3];
rz(0.82415165) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.2385345) q[2];
sx q[2];
rz(-1.3471341) q[2];
sx q[2];
rz(-1.6713589) q[2];
rz(3.0155731) q[3];
sx q[3];
rz(-2.6991548) q[3];
sx q[3];
rz(2.457705) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.11494342) q[0];
sx q[0];
rz(-2.6925955) q[0];
sx q[0];
rz(0.63585109) q[0];
rz(0.77330971) q[1];
sx q[1];
rz(-2.4644303) q[1];
sx q[1];
rz(-1.4453452) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6757386) q[0];
sx q[0];
rz(-2.2096118) q[0];
sx q[0];
rz(-2.5771228) q[0];
rz(-2.0436297) q[2];
sx q[2];
rz(-0.84542984) q[2];
sx q[2];
rz(2.5918317) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.6856001) q[1];
sx q[1];
rz(-1.7280792) q[1];
sx q[1];
rz(2.434751) q[1];
rz(-0.068447114) q[3];
sx q[3];
rz(-1.4395778) q[3];
sx q[3];
rz(0.4438627) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.9952952) q[2];
sx q[2];
rz(-1.7701021) q[2];
sx q[2];
rz(0.13776097) q[2];
rz(2.6855101) q[3];
sx q[3];
rz(-2.5465953) q[3];
sx q[3];
rz(-2.6661787) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59738961) q[0];
sx q[0];
rz(-1.5578288) q[0];
sx q[0];
rz(2.5624045) q[0];
rz(0.10313615) q[1];
sx q[1];
rz(-1.8201647) q[1];
sx q[1];
rz(-0.78027049) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0294757) q[0];
sx q[0];
rz(-1.5519841) q[0];
sx q[0];
rz(-0.095186724) q[0];
rz(-1.4109072) q[2];
sx q[2];
rz(-2.2457651) q[2];
sx q[2];
rz(-0.28197786) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.5835598) q[1];
sx q[1];
rz(-1.7906902) q[1];
sx q[1];
rz(-1.9755441) q[1];
x q[2];
rz(1.8336201) q[3];
sx q[3];
rz(-0.71280957) q[3];
sx q[3];
rz(-1.3808911) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.66037336) q[2];
sx q[2];
rz(-1.6796835) q[2];
sx q[2];
rz(-0.090864651) q[2];
rz(-1.4800492) q[3];
sx q[3];
rz(-1.816498) q[3];
sx q[3];
rz(0.81502325) q[3];
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
rz(-pi) q[3];
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
rz(-0.12260967) q[0];
sx q[0];
rz(-0.18773395) q[0];
sx q[0];
rz(-0.51938272) q[0];
rz(-1.1795562) q[1];
sx q[1];
rz(-0.68125454) q[1];
sx q[1];
rz(-0.048351668) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3607442) q[0];
sx q[0];
rz(-2.0791302) q[0];
sx q[0];
rz(-2.2717649) q[0];
rz(-pi) q[1];
rz(0.91421336) q[2];
sx q[2];
rz(-0.58744741) q[2];
sx q[2];
rz(2.812831) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.1241403) q[1];
sx q[1];
rz(-1.0832979) q[1];
sx q[1];
rz(-1.5392787) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8818733) q[3];
sx q[3];
rz(-0.81213299) q[3];
sx q[3];
rz(-0.5838697) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.4590596) q[2];
sx q[2];
rz(-2.8352663) q[2];
sx q[2];
rz(1.691386) q[2];
rz(-2.0356483) q[3];
sx q[3];
rz(-1.9927497) q[3];
sx q[3];
rz(-0.86627427) q[3];
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
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3806216) q[0];
sx q[0];
rz(-2.3411317) q[0];
sx q[0];
rz(-2.3642484) q[0];
rz(2.6457973) q[1];
sx q[1];
rz(-2.7340041) q[1];
sx q[1];
rz(0.19283238) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4254166) q[0];
sx q[0];
rz(-1.8766073) q[0];
sx q[0];
rz(-0.56434631) q[0];
rz(-pi) q[1];
rz(-1.8106789) q[2];
sx q[2];
rz(-2.3238306) q[2];
sx q[2];
rz(-2.8155934) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.7954171) q[1];
sx q[1];
rz(-0.3095135) q[1];
sx q[1];
rz(-1.8619821) q[1];
x q[2];
rz(-1.2019346) q[3];
sx q[3];
rz(-2.7141124) q[3];
sx q[3];
rz(-0.99586801) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.58482802) q[2];
sx q[2];
rz(-2.3636221) q[2];
sx q[2];
rz(-2.3053816) q[2];
rz(-0.95544514) q[3];
sx q[3];
rz(-1.7183869) q[3];
sx q[3];
rz(-2.6098765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8504976) q[0];
sx q[0];
rz(-1.9170772) q[0];
sx q[0];
rz(0.033705458) q[0];
rz(0.91824245) q[1];
sx q[1];
rz(-2.4372209) q[1];
sx q[1];
rz(-0.69923002) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.79141533) q[0];
sx q[0];
rz(-2.3473513) q[0];
sx q[0];
rz(2.3167531) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8644833) q[2];
sx q[2];
rz(-2.94706) q[2];
sx q[2];
rz(-2.3152318) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.0306176) q[1];
sx q[1];
rz(-0.4781107) q[1];
sx q[1];
rz(0.39076372) q[1];
rz(-pi) q[2];
rz(0.27856234) q[3];
sx q[3];
rz(-1.8980366) q[3];
sx q[3];
rz(0.95765169) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.08192) q[2];
sx q[2];
rz(-0.6531738) q[2];
sx q[2];
rz(-0.095452249) q[2];
rz(-1.0517612) q[3];
sx q[3];
rz(-1.8973408) q[3];
sx q[3];
rz(-0.89422798) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28602257) q[0];
sx q[0];
rz(-0.48208553) q[0];
sx q[0];
rz(1.1676189) q[0];
rz(-2.7883912) q[1];
sx q[1];
rz(-0.51274061) q[1];
sx q[1];
rz(0.28164992) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.024806) q[0];
sx q[0];
rz(-1.6580947) q[0];
sx q[0];
rz(-1.0952107) q[0];
rz(-pi) q[1];
rz(-3.0007576) q[2];
sx q[2];
rz(-0.65208921) q[2];
sx q[2];
rz(-1.9186794) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.49127625) q[1];
sx q[1];
rz(-1.4230844) q[1];
sx q[1];
rz(3.0056535) q[1];
rz(-pi) q[2];
rz(-2.2719215) q[3];
sx q[3];
rz(-2.4574349) q[3];
sx q[3];
rz(-1.9811833) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.99792751) q[2];
sx q[2];
rz(-1.2219656) q[2];
sx q[2];
rz(1.3227051) q[2];
rz(-1.5023242) q[3];
sx q[3];
rz(-2.0570677) q[3];
sx q[3];
rz(3.0433906) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.99763501) q[0];
sx q[0];
rz(-1.2459545) q[0];
sx q[0];
rz(2.8676046) q[0];
rz(-2.5471089) q[1];
sx q[1];
rz(-1.7285873) q[1];
sx q[1];
rz(-0.53057539) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0230867) q[0];
sx q[0];
rz(-3.102031) q[0];
sx q[0];
rz(-3.0450185) q[0];
rz(-pi) q[1];
rz(-1.4945901) q[2];
sx q[2];
rz(-1.0193902) q[2];
sx q[2];
rz(1.4879256) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.6484687) q[1];
sx q[1];
rz(-1.1965355) q[1];
sx q[1];
rz(-1.3918274) q[1];
rz(0.85803558) q[3];
sx q[3];
rz(-1.6948366) q[3];
sx q[3];
rz(0.8197561) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.24255594) q[2];
sx q[2];
rz(-1.9344923) q[2];
sx q[2];
rz(0.25699082) q[2];
rz(-1.9422003) q[3];
sx q[3];
rz(-1.896984) q[3];
sx q[3];
rz(-1.6283584) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5523819) q[0];
sx q[0];
rz(-0.95443812) q[0];
sx q[0];
rz(-0.5300262) q[0];
rz(1.6389305) q[1];
sx q[1];
rz(-1.9866147) q[1];
sx q[1];
rz(-2.2030305) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5650944) q[0];
sx q[0];
rz(-0.78701729) q[0];
sx q[0];
rz(2.4502343) q[0];
x q[1];
rz(-2.9979894) q[2];
sx q[2];
rz(-1.68949) q[2];
sx q[2];
rz(-0.76583662) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(3.0238513) q[1];
sx q[1];
rz(-1.7833901) q[1];
sx q[1];
rz(2.0076058) q[1];
rz(1.7750562) q[3];
sx q[3];
rz(-0.59511772) q[3];
sx q[3];
rz(0.68347733) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.75371257) q[2];
sx q[2];
rz(-2.2215999) q[2];
sx q[2];
rz(-2.06125) q[2];
rz(-2.3267817) q[3];
sx q[3];
rz(-1.1956513) q[3];
sx q[3];
rz(-1.4174392) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.998488) q[0];
sx q[0];
rz(-1.4838706) q[0];
sx q[0];
rz(2.0527573) q[0];
rz(0.067527436) q[1];
sx q[1];
rz(-1.2779526) q[1];
sx q[1];
rz(-0.75928226) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6367594) q[0];
sx q[0];
rz(-2.8593316) q[0];
sx q[0];
rz(2.0105848) q[0];
x q[1];
rz(-0.22589182) q[2];
sx q[2];
rz(-1.4498269) q[2];
sx q[2];
rz(-2.8793983) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.2893709) q[1];
sx q[1];
rz(-1.3557649) q[1];
sx q[1];
rz(2.380886) q[1];
x q[2];
rz(-1.3410232) q[3];
sx q[3];
rz(-1.9073351) q[3];
sx q[3];
rz(-2.0391109) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.8074983) q[2];
sx q[2];
rz(-2.0457902) q[2];
sx q[2];
rz(0.57541543) q[2];
rz(-2.5790162) q[3];
sx q[3];
rz(-0.96499363) q[3];
sx q[3];
rz(1.3728728) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
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
rz(0.060487735) q[0];
sx q[0];
rz(-1.2860379) q[0];
sx q[0];
rz(2.2338569) q[0];
rz(-1.0870712) q[1];
sx q[1];
rz(-2.0094951) q[1];
sx q[1];
rz(3.1241945) q[1];
rz(2.9870177) q[2];
sx q[2];
rz(-1.4308725) q[2];
sx q[2];
rz(1.9785656) q[2];
rz(-1.9418874) q[3];
sx q[3];
rz(-2.1368847) q[3];
sx q[3];
rz(0.18659244) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
