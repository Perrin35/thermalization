OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.4047591) q[0];
sx q[0];
rz(-1.7801378) q[0];
sx q[0];
rz(-1.7629495) q[0];
rz(-0.8575851) q[1];
sx q[1];
rz(-1.4839988) q[1];
sx q[1];
rz(-2.690697) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2338976) q[0];
sx q[0];
rz(-1.6595708) q[0];
sx q[0];
rz(1.5560454) q[0];
x q[1];
rz(3.0300006) q[2];
sx q[2];
rz(-1.0942232) q[2];
sx q[2];
rz(-3.1389719) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.7887468) q[1];
sx q[1];
rz(-2.5874918) q[1];
sx q[1];
rz(2.7098141) q[1];
rz(2.2545635) q[3];
sx q[3];
rz(-1.4966045) q[3];
sx q[3];
rz(0.85899734) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.9840055) q[2];
sx q[2];
rz(-1.459534) q[2];
sx q[2];
rz(-2.297304) q[2];
rz(0.44101161) q[3];
sx q[3];
rz(-0.35566548) q[3];
sx q[3];
rz(-0.60602337) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
rz(2.5490897) q[0];
sx q[0];
rz(-1.2298158) q[0];
sx q[0];
rz(-0.26309183) q[0];
rz(2.198055) q[1];
sx q[1];
rz(-2.5448006) q[1];
sx q[1];
rz(-1.9553604) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0378368) q[0];
sx q[0];
rz(-3.0825966) q[0];
sx q[0];
rz(2.8179413) q[0];
rz(1.6329174) q[2];
sx q[2];
rz(-1.3717692) q[2];
sx q[2];
rz(1.8387427) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.8801404) q[1];
sx q[1];
rz(-1.4457236) q[1];
sx q[1];
rz(-2.2373881) q[1];
rz(-2.9119592) q[3];
sx q[3];
rz(-2.4335055) q[3];
sx q[3];
rz(-1.9830444) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.1295604) q[2];
sx q[2];
rz(-1.0027145) q[2];
sx q[2];
rz(-1.9821232) q[2];
rz(-2.7705079) q[3];
sx q[3];
rz(-1.6371195) q[3];
sx q[3];
rz(-2.8306567) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.7611258) q[0];
sx q[0];
rz(-1.1300056) q[0];
sx q[0];
rz(-0.80672112) q[0];
rz(-0.21356788) q[1];
sx q[1];
rz(-0.49626207) q[1];
sx q[1];
rz(0.82021964) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0637174) q[0];
sx q[0];
rz(-0.88839196) q[0];
sx q[0];
rz(-0.14494411) q[0];
x q[1];
rz(1.0772466) q[2];
sx q[2];
rz(-1.3771257) q[2];
sx q[2];
rz(0.94656241) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.080938235) q[1];
sx q[1];
rz(-1.9229691) q[1];
sx q[1];
rz(2.0209795) q[1];
rz(-pi) q[2];
rz(-0.1180325) q[3];
sx q[3];
rz(-1.1088088) q[3];
sx q[3];
rz(3.1363917) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.31072581) q[2];
sx q[2];
rz(-1.6409637) q[2];
sx q[2];
rz(-0.93079981) q[2];
rz(-2.9860949) q[3];
sx q[3];
rz(-1.6379387) q[3];
sx q[3];
rz(0.29155198) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.61313066) q[0];
sx q[0];
rz(-0.72137946) q[0];
sx q[0];
rz(-0.91127515) q[0];
rz(-2.7032734) q[1];
sx q[1];
rz(-1.8194018) q[1];
sx q[1];
rz(-1.320425) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.60162773) q[0];
sx q[0];
rz(-2.7318582) q[0];
sx q[0];
rz(1.5967303) q[0];
rz(-2.7337012) q[2];
sx q[2];
rz(-0.48569277) q[2];
sx q[2];
rz(2.8913468) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.98852324) q[1];
sx q[1];
rz(-1.3077902) q[1];
sx q[1];
rz(2.4371229) q[1];
rz(-pi) q[2];
rz(-2.5370595) q[3];
sx q[3];
rz(-0.88450888) q[3];
sx q[3];
rz(1.3674919) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.1057672) q[2];
sx q[2];
rz(-2.2129009) q[2];
sx q[2];
rz(-0.34238112) q[2];
rz(2.9648182) q[3];
sx q[3];
rz(-2.7084559) q[3];
sx q[3];
rz(1.140973) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
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
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.3115561) q[0];
sx q[0];
rz(-0.72768584) q[0];
sx q[0];
rz(-0.86529055) q[0];
rz(1.226549) q[1];
sx q[1];
rz(-0.98926917) q[1];
sx q[1];
rz(1.3006166) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0375992) q[0];
sx q[0];
rz(-0.35498699) q[0];
sx q[0];
rz(-0.07261891) q[0];
rz(0.36655764) q[2];
sx q[2];
rz(-1.6379116) q[2];
sx q[2];
rz(1.3973438) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(3.0504426) q[1];
sx q[1];
rz(-1.0346518) q[1];
sx q[1];
rz(1.9552783) q[1];
rz(-pi) q[2];
rz(-2.6404068) q[3];
sx q[3];
rz(-0.29705829) q[3];
sx q[3];
rz(0.81851573) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.1317923) q[2];
sx q[2];
rz(-1.9921781) q[2];
sx q[2];
rz(2.664393) q[2];
rz(0.19208433) q[3];
sx q[3];
rz(-1.6936857) q[3];
sx q[3];
rz(2.208476) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3451097) q[0];
sx q[0];
rz(-0.61426291) q[0];
sx q[0];
rz(-0.011750301) q[0];
rz(-2.5911962) q[1];
sx q[1];
rz(-1.3563211) q[1];
sx q[1];
rz(-1.5531497) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5396744) q[0];
sx q[0];
rz(-1.2201021) q[0];
sx q[0];
rz(-2.7838216) q[0];
rz(-pi) q[1];
rz(1.2049098) q[2];
sx q[2];
rz(-1.591218) q[2];
sx q[2];
rz(0.4991971) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.8223871) q[1];
sx q[1];
rz(-0.86383312) q[1];
sx q[1];
rz(-1.6824526) q[1];
x q[2];
rz(0.86798571) q[3];
sx q[3];
rz(-1.6276974) q[3];
sx q[3];
rz(1.1662607) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.50756303) q[2];
sx q[2];
rz(-0.66036779) q[2];
sx q[2];
rz(1.2825512) q[2];
rz(-1.7717308) q[3];
sx q[3];
rz(-1.7462574) q[3];
sx q[3];
rz(-1.1184568) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5605374) q[0];
sx q[0];
rz(-0.16796172) q[0];
sx q[0];
rz(0.67725956) q[0];
rz(-0.15180763) q[1];
sx q[1];
rz(-1.3744524) q[1];
sx q[1];
rz(-2.1645434) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86977406) q[0];
sx q[0];
rz(-2.5442113) q[0];
sx q[0];
rz(3.0074044) q[0];
x q[1];
rz(-3.0095519) q[2];
sx q[2];
rz(-1.570895) q[2];
sx q[2];
rz(1.4591109) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.44927412) q[1];
sx q[1];
rz(-2.793503) q[1];
sx q[1];
rz(1.3432137) q[1];
rz(0.75120039) q[3];
sx q[3];
rz(-0.45414543) q[3];
sx q[3];
rz(2.5129012) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.7523505) q[2];
sx q[2];
rz(-0.82169473) q[2];
sx q[2];
rz(-1.0127257) q[2];
rz(-1.1879454) q[3];
sx q[3];
rz(-2.0690737) q[3];
sx q[3];
rz(-0.48721203) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0174719) q[0];
sx q[0];
rz(-0.033360632) q[0];
sx q[0];
rz(-0.69865984) q[0];
rz(2.0195122) q[1];
sx q[1];
rz(-2.2955003) q[1];
sx q[1];
rz(-1.8922071) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65803618) q[0];
sx q[0];
rz(-1.617096) q[0];
sx q[0];
rz(-0.21090837) q[0];
rz(2.0935358) q[2];
sx q[2];
rz(-1.9881696) q[2];
sx q[2];
rz(-1.0440895) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.3932712) q[1];
sx q[1];
rz(-1.5763348) q[1];
sx q[1];
rz(0.022151532) q[1];
rz(-pi) q[2];
rz(2.8532393) q[3];
sx q[3];
rz(-0.95803146) q[3];
sx q[3];
rz(1.8307277) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.32968783) q[2];
sx q[2];
rz(-2.3554282) q[2];
sx q[2];
rz(-1.9630986) q[2];
rz(1.4568436) q[3];
sx q[3];
rz(-1.0624351) q[3];
sx q[3];
rz(-2.7594574) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6417398) q[0];
sx q[0];
rz(-1.7620182) q[0];
sx q[0];
rz(1.2930124) q[0];
rz(1.7199843) q[1];
sx q[1];
rz(-1.0363818) q[1];
sx q[1];
rz(2.5440149) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15411988) q[0];
sx q[0];
rz(-0.14053908) q[0];
sx q[0];
rz(-1.2921635) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0963247) q[2];
sx q[2];
rz(-2.0647486) q[2];
sx q[2];
rz(-1.8906821) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.5850726) q[1];
sx q[1];
rz(-0.44610281) q[1];
sx q[1];
rz(-2.3162566) q[1];
rz(-pi) q[2];
rz(0.58595539) q[3];
sx q[3];
rz(-1.6349941) q[3];
sx q[3];
rz(-0.5639329) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.9188345) q[2];
sx q[2];
rz(-1.6901878) q[2];
sx q[2];
rz(1.9082327) q[2];
rz(0.90138609) q[3];
sx q[3];
rz(-0.12005761) q[3];
sx q[3];
rz(1.4982769) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0937061) q[0];
sx q[0];
rz(-2.369635) q[0];
sx q[0];
rz(0.023660252) q[0];
rz(-2.1854782) q[1];
sx q[1];
rz(-1.8319943) q[1];
sx q[1];
rz(2.4694209) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4155054) q[0];
sx q[0];
rz(-3.130075) q[0];
sx q[0];
rz(2.0477717) q[0];
x q[1];
rz(2.6300738) q[2];
sx q[2];
rz(-1.5626972) q[2];
sx q[2];
rz(2.9715003) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.9384267) q[1];
sx q[1];
rz(-0.90183479) q[1];
sx q[1];
rz(1.1346362) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8176114) q[3];
sx q[3];
rz(-2.6880662) q[3];
sx q[3];
rz(1.4676263) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.6293634) q[2];
sx q[2];
rz(-1.2269292) q[2];
sx q[2];
rz(2.771634) q[2];
rz(-1.6379179) q[3];
sx q[3];
rz(-0.88589293) q[3];
sx q[3];
rz(1.9406208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5794012) q[0];
sx q[0];
rz(-0.36407064) q[0];
sx q[0];
rz(-1.9343485) q[0];
rz(0.72369408) q[1];
sx q[1];
rz(-0.98725286) q[1];
sx q[1];
rz(-0.90686803) q[1];
rz(1.205668) q[2];
sx q[2];
rz(-0.42386133) q[2];
sx q[2];
rz(-0.78122666) q[2];
rz(-2.1914992) q[3];
sx q[3];
rz(-0.48662574) q[3];
sx q[3];
rz(1.1403198) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];