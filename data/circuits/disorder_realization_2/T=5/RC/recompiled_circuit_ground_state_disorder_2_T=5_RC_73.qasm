OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.53192294) q[0];
sx q[0];
rz(-1.5250396) q[0];
sx q[0];
rz(1.2582231) q[0];
rz(7.8426709) q[1];
sx q[1];
rz(5.7962228) q[1];
sx q[1];
rz(16.143057) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5814911) q[0];
sx q[0];
rz(-2.1061497) q[0];
sx q[0];
rz(1.0161607) q[0];
x q[1];
rz(2.6039106) q[2];
sx q[2];
rz(-1.8694365) q[2];
sx q[2];
rz(2.5387105) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.62578979) q[1];
sx q[1];
rz(-2.8608659) q[1];
sx q[1];
rz(1.2899542) q[1];
x q[2];
rz(-0.70617999) q[3];
sx q[3];
rz(-0.67863388) q[3];
sx q[3];
rz(2.5265001) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.9293999) q[2];
sx q[2];
rz(-1.8636899) q[2];
sx q[2];
rz(-1.0431935) q[2];
rz(2.7405401) q[3];
sx q[3];
rz(-1.4570823) q[3];
sx q[3];
rz(2.5636087) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.22724085) q[0];
sx q[0];
rz(-0.14509097) q[0];
sx q[0];
rz(2.5703854) q[0];
rz(-0.97194833) q[1];
sx q[1];
rz(-2.4010039) q[1];
sx q[1];
rz(-2.0170225) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1069431) q[0];
sx q[0];
rz(-0.9616344) q[0];
sx q[0];
rz(-0.27584313) q[0];
rz(-pi) q[1];
rz(0.6508462) q[2];
sx q[2];
rz(-1.2799601) q[2];
sx q[2];
rz(-0.42465045) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.1536869) q[1];
sx q[1];
rz(-0.73927727) q[1];
sx q[1];
rz(1.4312782) q[1];
rz(-pi) q[2];
rz(1.589682) q[3];
sx q[3];
rz(-1.369993) q[3];
sx q[3];
rz(-2.7310179) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.22175114) q[2];
sx q[2];
rz(-1.5849179) q[2];
sx q[2];
rz(0.43608967) q[2];
rz(0.74887216) q[3];
sx q[3];
rz(-0.80513969) q[3];
sx q[3];
rz(-2.3124636) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1083168) q[0];
sx q[0];
rz(-1.746614) q[0];
sx q[0];
rz(0.088031553) q[0];
rz(2.2875359) q[1];
sx q[1];
rz(-0.66103926) q[1];
sx q[1];
rz(-2.0565654) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89660925) q[0];
sx q[0];
rz(-1.6379781) q[0];
sx q[0];
rz(-1.843473) q[0];
rz(-pi) q[1];
rz(-1.3543812) q[2];
sx q[2];
rz(-1.7209574) q[2];
sx q[2];
rz(0.79126287) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.4702556) q[1];
sx q[1];
rz(-0.47154676) q[1];
sx q[1];
rz(-1.2045453) q[1];
x q[2];
rz(1.4973699) q[3];
sx q[3];
rz(-1.182397) q[3];
sx q[3];
rz(-1.8457613) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.6701086) q[2];
sx q[2];
rz(-1.4205616) q[2];
sx q[2];
rz(-2.3954083) q[2];
rz(-2.9255195) q[3];
sx q[3];
rz(-1.2553299) q[3];
sx q[3];
rz(2.9358673) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3716607) q[0];
sx q[0];
rz(-0.07337229) q[0];
sx q[0];
rz(1.7727456) q[0];
rz(0.61028496) q[1];
sx q[1];
rz(-1.3657602) q[1];
sx q[1];
rz(-0.38017166) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.73715009) q[0];
sx q[0];
rz(-1.2594873) q[0];
sx q[0];
rz(2.690517) q[0];
rz(2.751558) q[2];
sx q[2];
rz(-0.48214285) q[2];
sx q[2];
rz(2.6576561) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.026555) q[1];
sx q[1];
rz(-1.4241744) q[1];
sx q[1];
rz(-1.0653593) q[1];
x q[2];
rz(-0.065297619) q[3];
sx q[3];
rz(-1.9717798) q[3];
sx q[3];
rz(-0.15006615) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.88468203) q[2];
sx q[2];
rz(-2.1941049) q[2];
sx q[2];
rz(2.103629) q[2];
rz(0.3642309) q[3];
sx q[3];
rz(-0.7887775) q[3];
sx q[3];
rz(-2.4558333) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6895741) q[0];
sx q[0];
rz(-1.7114534) q[0];
sx q[0];
rz(0.83056393) q[0];
rz(-0.05016249) q[1];
sx q[1];
rz(-0.12718931) q[1];
sx q[1];
rz(-0.62279472) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2723214) q[0];
sx q[0];
rz(-1.6487445) q[0];
sx q[0];
rz(0.21792953) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7217614) q[2];
sx q[2];
rz(-1.3908885) q[2];
sx q[2];
rz(-2.3617488) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.2166136) q[1];
sx q[1];
rz(-2.2053524) q[1];
sx q[1];
rz(-0.44559176) q[1];
rz(-3.1386802) q[3];
sx q[3];
rz(-0.022449819) q[3];
sx q[3];
rz(-2.8685304) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.097816618) q[2];
sx q[2];
rz(-0.4349097) q[2];
sx q[2];
rz(-0.80769509) q[2];
rz(-1.5378599) q[3];
sx q[3];
rz(-1.5065498) q[3];
sx q[3];
rz(2.9174771) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1767126) q[0];
sx q[0];
rz(-1.8105312) q[0];
sx q[0];
rz(-1.4759195) q[0];
rz(1.9077979) q[1];
sx q[1];
rz(-1.7314311) q[1];
sx q[1];
rz(0.15353157) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3672402) q[0];
sx q[0];
rz(-1.4611547) q[0];
sx q[0];
rz(1.5745274) q[0];
x q[1];
rz(0.84014272) q[2];
sx q[2];
rz(-2.5135927) q[2];
sx q[2];
rz(-0.59115228) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(3.081806) q[1];
sx q[1];
rz(-2.012552) q[1];
sx q[1];
rz(-0.55952832) q[1];
rz(1.4755523) q[3];
sx q[3];
rz(-0.93425865) q[3];
sx q[3];
rz(1.7004636) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.70225707) q[2];
sx q[2];
rz(-0.80436891) q[2];
sx q[2];
rz(-2.5852933) q[2];
rz(-0.016911658) q[3];
sx q[3];
rz(-2.1664679) q[3];
sx q[3];
rz(-0.66439605) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.21338129) q[0];
sx q[0];
rz(-3.0688372) q[0];
sx q[0];
rz(0.67759204) q[0];
rz(1.3202336) q[1];
sx q[1];
rz(-2.0537387) q[1];
sx q[1];
rz(-0.56603146) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4844263) q[0];
sx q[0];
rz(-0.56492794) q[0];
sx q[0];
rz(0.66175263) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6230311) q[2];
sx q[2];
rz(-1.5855316) q[2];
sx q[2];
rz(-2.7434204) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.5929451) q[1];
sx q[1];
rz(-0.8742395) q[1];
sx q[1];
rz(2.9440483) q[1];
rz(1.5291847) q[3];
sx q[3];
rz(-2.0609612) q[3];
sx q[3];
rz(-3.0050638) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.90938202) q[2];
sx q[2];
rz(-2.7558694) q[2];
sx q[2];
rz(-2.5717226) q[2];
rz(-1.0424403) q[3];
sx q[3];
rz(-1.9614377) q[3];
sx q[3];
rz(-1.0038143) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37782272) q[0];
sx q[0];
rz(-0.42709392) q[0];
sx q[0];
rz(1.0231934) q[0];
rz(0.91776735) q[1];
sx q[1];
rz(-0.75420165) q[1];
sx q[1];
rz(-0.86407152) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7650493) q[0];
sx q[0];
rz(-1.4481539) q[0];
sx q[0];
rz(-1.8160233) q[0];
rz(-pi) q[1];
x q[1];
rz(2.239924) q[2];
sx q[2];
rz(-2.5732917) q[2];
sx q[2];
rz(-1.2327884) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(3.1004677) q[1];
sx q[1];
rz(-2.0534614) q[1];
sx q[1];
rz(0.80759279) q[1];
rz(2.731063) q[3];
sx q[3];
rz(-2.3709112) q[3];
sx q[3];
rz(0.16134027) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.0920022) q[2];
sx q[2];
rz(-1.2668173) q[2];
sx q[2];
rz(-2.3769296) q[2];
rz(2.2406254) q[3];
sx q[3];
rz(-2.4668906) q[3];
sx q[3];
rz(2.0425792) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6242999) q[0];
sx q[0];
rz(-0.39532548) q[0];
sx q[0];
rz(-0.98315352) q[0];
rz(-0.82955018) q[1];
sx q[1];
rz(-1.7770276) q[1];
sx q[1];
rz(-0.57873908) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1252223) q[0];
sx q[0];
rz(-1.9645321) q[0];
sx q[0];
rz(-3.1278174) q[0];
rz(-pi) q[1];
rz(2.3624973) q[2];
sx q[2];
rz(-0.70487521) q[2];
sx q[2];
rz(-1.5388956) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.14565258) q[1];
sx q[1];
rz(-2.4697587) q[1];
sx q[1];
rz(-2.2493659) q[1];
rz(-pi) q[2];
rz(1.6074568) q[3];
sx q[3];
rz(-1.9158773) q[3];
sx q[3];
rz(2.7861129) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.4911554) q[2];
sx q[2];
rz(-1.8530242) q[2];
sx q[2];
rz(0.11162652) q[2];
rz(2.1857183) q[3];
sx q[3];
rz(-2.1501232) q[3];
sx q[3];
rz(-1.8808232) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.88804203) q[0];
sx q[0];
rz(-1.0567867) q[0];
sx q[0];
rz(0.014130935) q[0];
rz(1.1125394) q[1];
sx q[1];
rz(-2.2281149) q[1];
sx q[1];
rz(-1.5412451) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.74792101) q[0];
sx q[0];
rz(-1.4427516) q[0];
sx q[0];
rz(-2.9927954) q[0];
rz(-pi) q[1];
rz(-2.2152973) q[2];
sx q[2];
rz(-1.6716812) q[2];
sx q[2];
rz(-0.48357329) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.66110509) q[1];
sx q[1];
rz(-0.36927642) q[1];
sx q[1];
rz(-0.1657471) q[1];
x q[2];
rz(2.9716206) q[3];
sx q[3];
rz(-1.7331496) q[3];
sx q[3];
rz(-1.270783) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.9666226) q[2];
sx q[2];
rz(-1.2302159) q[2];
sx q[2];
rz(0.54565412) q[2];
rz(-0.086006554) q[3];
sx q[3];
rz(-1.514148) q[3];
sx q[3];
rz(-2.8470794) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.94588146) q[0];
sx q[0];
rz(-2.0989037) q[0];
sx q[0];
rz(1.209191) q[0];
rz(0.60278268) q[1];
sx q[1];
rz(-2.5310015) q[1];
sx q[1];
rz(-2.3878154) q[1];
rz(1.8281585) q[2];
sx q[2];
rz(-2.4429697) q[2];
sx q[2];
rz(2.2945678) q[2];
rz(1.6794593) q[3];
sx q[3];
rz(-2.1908219) q[3];
sx q[3];
rz(0.11809668) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
