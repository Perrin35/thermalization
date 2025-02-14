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
rz(2.8132791) q[0];
sx q[0];
rz(-2.0067196) q[0];
sx q[0];
rz(2.5757134) q[0];
rz(0.24428754) q[1];
sx q[1];
rz(-1.3574418) q[1];
sx q[1];
rz(0.53425962) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8847647) q[0];
sx q[0];
rz(-2.3951911) q[0];
sx q[0];
rz(1.9749667) q[0];
x q[1];
rz(-1.2674321) q[2];
sx q[2];
rz(-2.144226) q[2];
sx q[2];
rz(-1.0614741) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.6837343) q[1];
sx q[1];
rz(-2.8401655) q[1];
sx q[1];
rz(3.0617759) q[1];
rz(-pi) q[2];
x q[2];
rz(1.322106) q[3];
sx q[3];
rz(-2.0157031) q[3];
sx q[3];
rz(-1.1550732) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.013022097) q[2];
sx q[2];
rz(-0.2505005) q[2];
sx q[2];
rz(2.8924083) q[2];
rz(1.3439882) q[3];
sx q[3];
rz(-1.5430217) q[3];
sx q[3];
rz(0.070076076) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
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
rz(-0.4942112) q[0];
sx q[0];
rz(-1.9152315) q[0];
sx q[0];
rz(-0.98980728) q[0];
rz(0.79065943) q[1];
sx q[1];
rz(-0.76690563) q[1];
sx q[1];
rz(-0.31165037) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3450736) q[0];
sx q[0];
rz(-0.88788827) q[0];
sx q[0];
rz(0.89181283) q[0];
rz(-pi) q[1];
rz(3.0333854) q[2];
sx q[2];
rz(-1.756486) q[2];
sx q[2];
rz(-0.78150392) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-3.0428704) q[1];
sx q[1];
rz(-1.8095892) q[1];
sx q[1];
rz(2.9643994) q[1];
rz(-pi) q[2];
rz(2.8007224) q[3];
sx q[3];
rz(-0.98353993) q[3];
sx q[3];
rz(1.0140918) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.0011757294) q[2];
sx q[2];
rz(-1.1602465) q[2];
sx q[2];
rz(-2.1902093) q[2];
rz(2.9133255) q[3];
sx q[3];
rz(-0.97672668) q[3];
sx q[3];
rz(1.2736646) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9839639) q[0];
sx q[0];
rz(-1.552859) q[0];
sx q[0];
rz(-0.11446318) q[0];
rz(-2.139367) q[1];
sx q[1];
rz(-2.7148235) q[1];
sx q[1];
rz(0.1618298) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9660258) q[0];
sx q[0];
rz(-1.1822261) q[0];
sx q[0];
rz(-0.92312621) q[0];
x q[1];
rz(-0.21185565) q[2];
sx q[2];
rz(-2.7447034) q[2];
sx q[2];
rz(2.3174163) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.5787631) q[1];
sx q[1];
rz(-2.17318) q[1];
sx q[1];
rz(-1.8627326) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7416218) q[3];
sx q[3];
rz(-0.84163991) q[3];
sx q[3];
rz(-0.10310452) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.74730045) q[2];
sx q[2];
rz(-0.44555274) q[2];
sx q[2];
rz(0.14075819) q[2];
rz(2.5898139) q[3];
sx q[3];
rz(-1.8675624) q[3];
sx q[3];
rz(-2.3902635) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
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
rz(-0.20525876) q[0];
sx q[0];
rz(-0.2660428) q[0];
sx q[0];
rz(0.67489135) q[0];
rz(0.15448013) q[1];
sx q[1];
rz(-1.7325502) q[1];
sx q[1];
rz(-2.6836269) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5602011) q[0];
sx q[0];
rz(-1.0835674) q[0];
sx q[0];
rz(1.1283406) q[0];
rz(-pi) q[1];
x q[1];
rz(0.9651009) q[2];
sx q[2];
rz(-2.6332246) q[2];
sx q[2];
rz(1.8338211) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.030589015) q[1];
sx q[1];
rz(-1.9408424) q[1];
sx q[1];
rz(-1.7683328) q[1];
x q[2];
rz(-1.4639888) q[3];
sx q[3];
rz(-1.065101) q[3];
sx q[3];
rz(-0.13733521) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(3.0080879) q[2];
sx q[2];
rz(-0.68024457) q[2];
sx q[2];
rz(-0.31944719) q[2];
rz(-1.3301298) q[3];
sx q[3];
rz(-2.1250686) q[3];
sx q[3];
rz(0.5602079) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9820246) q[0];
sx q[0];
rz(-1.1363131) q[0];
sx q[0];
rz(1.2215479) q[0];
rz(0.23280652) q[1];
sx q[1];
rz(-2.5722952) q[1];
sx q[1];
rz(-0.98308841) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1637835) q[0];
sx q[0];
rz(-2.0931232) q[0];
sx q[0];
rz(-2.1420826) q[0];
rz(-3.0276322) q[2];
sx q[2];
rz(-0.95319213) q[2];
sx q[2];
rz(-1.8007939) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.32199643) q[1];
sx q[1];
rz(-2.9038972) q[1];
sx q[1];
rz(2.4758215) q[1];
rz(-pi) q[2];
x q[2];
rz(0.1666937) q[3];
sx q[3];
rz(-1.987926) q[3];
sx q[3];
rz(2.6670109) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.7588707) q[2];
sx q[2];
rz(-1.4843586) q[2];
sx q[2];
rz(0.37528428) q[2];
rz(-2.0659857) q[3];
sx q[3];
rz(-0.38638249) q[3];
sx q[3];
rz(-2.6066499) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(1.7406834) q[0];
sx q[0];
rz(-1.704819) q[0];
sx q[0];
rz(1.4373454) q[0];
rz(3.0628693) q[1];
sx q[1];
rz(-1.3767786) q[1];
sx q[1];
rz(-2.5934503) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5054436) q[0];
sx q[0];
rz(-2.4432805) q[0];
sx q[0];
rz(0.33238073) q[0];
x q[1];
rz(2.5896543) q[2];
sx q[2];
rz(-1.1480398) q[2];
sx q[2];
rz(1.1814337) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.2386475) q[1];
sx q[1];
rz(-0.79552197) q[1];
sx q[1];
rz(-2.5832157) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.94977149) q[3];
sx q[3];
rz(-2.2249537) q[3];
sx q[3];
rz(2.014924) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.4569725) q[2];
sx q[2];
rz(-0.62070864) q[2];
sx q[2];
rz(2.4662628) q[2];
rz(-1.5843377) q[3];
sx q[3];
rz(-1.8341589) q[3];
sx q[3];
rz(-2.5244782) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(-2.5976582) q[0];
sx q[0];
rz(-1.958853) q[0];
sx q[0];
rz(0.4959929) q[0];
rz(1.687382) q[1];
sx q[1];
rz(-2.0861237) q[1];
sx q[1];
rz(-2.0248263) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6750609) q[0];
sx q[0];
rz(-1.2197615) q[0];
sx q[0];
rz(-2.768633) q[0];
x q[1];
rz(1.1363813) q[2];
sx q[2];
rz(-2.5602753) q[2];
sx q[2];
rz(0.28806799) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.36193725) q[1];
sx q[1];
rz(-1.2314267) q[1];
sx q[1];
rz(1.2324059) q[1];
rz(0.12567606) q[3];
sx q[3];
rz(-0.45914859) q[3];
sx q[3];
rz(-0.3885551) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.83069673) q[2];
sx q[2];
rz(-2.1572957) q[2];
sx q[2];
rz(-2.7867219) q[2];
rz(-2.3762459) q[3];
sx q[3];
rz(-0.86330515) q[3];
sx q[3];
rz(3.0520458) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4527721) q[0];
sx q[0];
rz(-1.6625762) q[0];
sx q[0];
rz(-2.3759957) q[0];
rz(0.87144026) q[1];
sx q[1];
rz(-2.3759418) q[1];
sx q[1];
rz(-1.8188459) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8894371) q[0];
sx q[0];
rz(-1.8312635) q[0];
sx q[0];
rz(-1.5171264) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9051612) q[2];
sx q[2];
rz(-2.2364103) q[2];
sx q[2];
rz(0.98275358) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.2603399) q[1];
sx q[1];
rz(-1.5872721) q[1];
sx q[1];
rz(0.50838281) q[1];
x q[2];
rz(-2.3115308) q[3];
sx q[3];
rz(-2.8403211) q[3];
sx q[3];
rz(-0.64727441) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.65779296) q[2];
sx q[2];
rz(-1.0427661) q[2];
sx q[2];
rz(2.9198809) q[2];
rz(-2.0486369) q[3];
sx q[3];
rz(-1.4128128) q[3];
sx q[3];
rz(0.67813412) q[3];
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
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.073864989) q[0];
sx q[0];
rz(-1.2465957) q[0];
sx q[0];
rz(-2.8644323) q[0];
rz(2.3932338) q[1];
sx q[1];
rz(-2.6942418) q[1];
sx q[1];
rz(-1.1433196) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4410981) q[0];
sx q[0];
rz(-1.7713393) q[0];
sx q[0];
rz(-1.2593377) q[0];
x q[1];
rz(-1.3242993) q[2];
sx q[2];
rz(-0.72028226) q[2];
sx q[2];
rz(-1.5697073) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.34830293) q[1];
sx q[1];
rz(-1.7277246) q[1];
sx q[1];
rz(1.2497529) q[1];
rz(-pi) q[2];
rz(-2.7797749) q[3];
sx q[3];
rz(-1.731195) q[3];
sx q[3];
rz(-0.058506207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.2271759) q[2];
sx q[2];
rz(-1.02355) q[2];
sx q[2];
rz(0.048967036) q[2];
rz(-2.07552) q[3];
sx q[3];
rz(-0.62838307) q[3];
sx q[3];
rz(-2.9858885) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.12298909) q[0];
sx q[0];
rz(-1.6139231) q[0];
sx q[0];
rz(0.019512026) q[0];
rz(0.77876577) q[1];
sx q[1];
rz(-2.4486783) q[1];
sx q[1];
rz(1.5787517) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.21694788) q[0];
sx q[0];
rz(-2.9836982) q[0];
sx q[0];
rz(0.96952207) q[0];
x q[1];
rz(-1.8291924) q[2];
sx q[2];
rz(-0.62587291) q[2];
sx q[2];
rz(1.3385119) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.07933087) q[1];
sx q[1];
rz(-2.375646) q[1];
sx q[1];
rz(-2.2205611) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4061798) q[3];
sx q[3];
rz(-2.4009628) q[3];
sx q[3];
rz(3.0082138) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.94893) q[2];
sx q[2];
rz(-0.86343416) q[2];
sx q[2];
rz(2.9760402) q[2];
rz(0.60756573) q[3];
sx q[3];
rz(-1.2847885) q[3];
sx q[3];
rz(-2.6732388) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7623357) q[0];
sx q[0];
rz(-2.6597334) q[0];
sx q[0];
rz(-0.045482176) q[0];
rz(1.5070076) q[1];
sx q[1];
rz(-1.4611117) q[1];
sx q[1];
rz(-1.5932105) q[1];
rz(2.3661939) q[2];
sx q[2];
rz(-1.4682894) q[2];
sx q[2];
rz(-0.52649047) q[2];
rz(-2.8589917) q[3];
sx q[3];
rz(-2.5167214) q[3];
sx q[3];
rz(-1.4759397) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
