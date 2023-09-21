OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.73683357) q[0];
sx q[0];
rz(4.9217304) q[0];
sx q[0];
rz(11.187727) q[0];
rz(2.2840075) q[1];
sx q[1];
rz(-1.6575939) q[1];
sx q[1];
rz(-0.4508957) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3987797) q[0];
sx q[0];
rz(-3.0516041) q[0];
sx q[0];
rz(-2.9773657) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0916753) q[2];
sx q[2];
rz(-1.6699104) q[2];
sx q[2];
rz(-1.5168158) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.7887468) q[1];
sx q[1];
rz(-2.5874918) q[1];
sx q[1];
rz(0.43177859) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6879184) q[3];
sx q[3];
rz(-2.4544567) q[3];
sx q[3];
rz(-0.62108921) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.1575872) q[2];
sx q[2];
rz(-1.6820587) q[2];
sx q[2];
rz(0.84428865) q[2];
rz(0.44101161) q[3];
sx q[3];
rz(-2.7859272) q[3];
sx q[3];
rz(0.60602337) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59250295) q[0];
sx q[0];
rz(-1.2298158) q[0];
sx q[0];
rz(0.26309183) q[0];
rz(2.198055) q[1];
sx q[1];
rz(-2.5448006) q[1];
sx q[1];
rz(1.1862322) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1037558) q[0];
sx q[0];
rz(-0.058996011) q[0];
sx q[0];
rz(-0.32365139) q[0];
rz(2.942191) q[2];
sx q[2];
rz(-1.5099031) q[2];
sx q[2];
rz(-2.8859438) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.26145229) q[1];
sx q[1];
rz(-1.6958691) q[1];
sx q[1];
rz(2.2373881) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7632742) q[3];
sx q[3];
rz(-2.2566183) q[3];
sx q[3];
rz(-1.457085) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.0120323) q[2];
sx q[2];
rz(-2.1388781) q[2];
sx q[2];
rz(-1.1594695) q[2];
rz(-0.37108478) q[3];
sx q[3];
rz(-1.5044731) q[3];
sx q[3];
rz(-2.8306567) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7611258) q[0];
sx q[0];
rz(-1.1300056) q[0];
sx q[0];
rz(0.80672112) q[0];
rz(-0.21356788) q[1];
sx q[1];
rz(-0.49626207) q[1];
sx q[1];
rz(0.82021964) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0778753) q[0];
sx q[0];
rz(-2.2532007) q[0];
sx q[0];
rz(2.9966485) q[0];
rz(-pi) q[1];
rz(-2.9224612) q[2];
sx q[2];
rz(-1.0872772) q[2];
sx q[2];
rz(-0.52106524) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.1097475) q[1];
sx q[1];
rz(-0.56402962) q[1];
sx q[1];
rz(-0.86947039) q[1];
rz(-pi) q[2];
rz(1.8030274) q[3];
sx q[3];
rz(-0.47577061) q[3];
sx q[3];
rz(0.25482086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.31072581) q[2];
sx q[2];
rz(-1.5006289) q[2];
sx q[2];
rz(-2.2107928) q[2];
rz(-2.9860949) q[3];
sx q[3];
rz(-1.5036539) q[3];
sx q[3];
rz(2.8500407) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.61313066) q[0];
sx q[0];
rz(-2.4202132) q[0];
sx q[0];
rz(2.2303175) q[0];
rz(2.7032734) q[1];
sx q[1];
rz(-1.8194018) q[1];
sx q[1];
rz(-1.8211676) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.94538044) q[0];
sx q[0];
rz(-1.5811265) q[0];
sx q[0];
rz(1.1611847) q[0];
rz(-2.690372) q[2];
sx q[2];
rz(-1.7570474) q[2];
sx q[2];
rz(-0.95552432) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.3418386) q[1];
sx q[1];
rz(-2.2463887) q[1];
sx q[1];
rz(-1.9104596) q[1];
x q[2];
rz(0.78762357) q[3];
sx q[3];
rz(-2.0260603) q[3];
sx q[3];
rz(0.61592197) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.1057672) q[2];
sx q[2];
rz(-2.2129009) q[2];
sx q[2];
rz(0.34238112) q[2];
rz(-2.9648182) q[3];
sx q[3];
rz(-2.7084559) q[3];
sx q[3];
rz(2.0006196) q[3];
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
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8300366) q[0];
sx q[0];
rz(-2.4139068) q[0];
sx q[0];
rz(2.2763021) q[0];
rz(1.9150437) q[1];
sx q[1];
rz(-2.1523235) q[1];
sx q[1];
rz(-1.8409761) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1814225) q[0];
sx q[0];
rz(-1.2167861) q[0];
sx q[0];
rz(-1.5439073) q[0];
rz(-pi) q[1];
rz(1.6426716) q[2];
sx q[2];
rz(-1.9364898) q[2];
sx q[2];
rz(-2.9938811) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.76104858) q[1];
sx q[1];
rz(-0.64861464) q[1];
sx q[1];
rz(0.56306871) q[1];
x q[2];
rz(2.8793094) q[3];
sx q[3];
rz(-1.4296921) q[3];
sx q[3];
rz(-1.9067681) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.1317923) q[2];
sx q[2];
rz(-1.9921781) q[2];
sx q[2];
rz(-2.664393) q[2];
rz(2.9495083) q[3];
sx q[3];
rz(-1.447907) q[3];
sx q[3];
rz(-0.93311667) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(2.3451097) q[0];
sx q[0];
rz(-0.61426291) q[0];
sx q[0];
rz(3.1298424) q[0];
rz(-0.55039644) q[1];
sx q[1];
rz(-1.7852716) q[1];
sx q[1];
rz(-1.5531497) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5396744) q[0];
sx q[0];
rz(-1.9214905) q[0];
sx q[0];
rz(2.7838216) q[0];
rz(-pi) q[1];
x q[1];
rz(3.1197238) q[2];
sx q[2];
rz(-1.936603) q[2];
sx q[2];
rz(1.0794229) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.1482684) q[1];
sx q[1];
rz(-0.71422186) q[1];
sx q[1];
rz(0.12970129) q[1];
rz(-0.86798571) q[3];
sx q[3];
rz(-1.5138953) q[3];
sx q[3];
rz(-1.9753319) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.50756303) q[2];
sx q[2];
rz(-2.4812249) q[2];
sx q[2];
rz(1.8590415) q[2];
rz(1.3698618) q[3];
sx q[3];
rz(-1.7462574) q[3];
sx q[3];
rz(2.0231358) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.58105528) q[0];
sx q[0];
rz(-0.16796172) q[0];
sx q[0];
rz(-2.4643331) q[0];
rz(2.989785) q[1];
sx q[1];
rz(-1.3744524) q[1];
sx q[1];
rz(-2.1645434) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2718186) q[0];
sx q[0];
rz(-2.5442113) q[0];
sx q[0];
rz(-0.13418829) q[0];
x q[1];
rz(-1.5706967) q[2];
sx q[2];
rz(-1.4387555) q[2];
sx q[2];
rz(3.0299203) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.2077142) q[1];
sx q[1];
rz(-1.9095451) q[1];
sx q[1];
rz(-3.0599041) q[1];
x q[2];
rz(2.3903923) q[3];
sx q[3];
rz(-2.6874472) q[3];
sx q[3];
rz(-0.62869149) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.3892422) q[2];
sx q[2];
rz(-0.82169473) q[2];
sx q[2];
rz(-2.1288669) q[2];
rz(-1.1879454) q[3];
sx q[3];
rz(-2.0690737) q[3];
sx q[3];
rz(-0.48721203) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
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
rz(-2.1241207) q[0];
sx q[0];
rz(-0.033360632) q[0];
sx q[0];
rz(-0.69865984) q[0];
rz(-1.1220804) q[1];
sx q[1];
rz(-2.2955003) q[1];
sx q[1];
rz(-1.8922071) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65803618) q[0];
sx q[0];
rz(-1.617096) q[0];
sx q[0];
rz(2.9306843) q[0];
rz(-pi) q[1];
rz(2.6685733) q[2];
sx q[2];
rz(-2.0447391) q[2];
sx q[2];
rz(-2.3854286) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.5640806) q[1];
sx q[1];
rz(-3.1187594) q[1];
sx q[1];
rz(-2.8965685) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1864248) q[3];
sx q[3];
rz(-0.66925183) q[3];
sx q[3];
rz(0.8347019) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.8119048) q[2];
sx q[2];
rz(-0.78616443) q[2];
sx q[2];
rz(1.9630986) q[2];
rz(-1.4568436) q[3];
sx q[3];
rz(-1.0624351) q[3];
sx q[3];
rz(-0.38213521) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6417398) q[0];
sx q[0];
rz(-1.3795744) q[0];
sx q[0];
rz(-1.2930124) q[0];
rz(1.7199843) q[1];
sx q[1];
rz(-1.0363818) q[1];
sx q[1];
rz(2.5440149) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4353838) q[0];
sx q[0];
rz(-1.4357114) q[0];
sx q[0];
rz(-3.1027017) q[0];
rz(-pi) q[1];
rz(0.75002807) q[2];
sx q[2];
rz(-0.70493297) q[2];
sx q[2];
rz(-0.36545576) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.4328879) q[1];
sx q[1];
rz(-1.27379) q[1];
sx q[1];
rz(1.2328641) q[1];
rz(-pi) q[2];
rz(-2.5556373) q[3];
sx q[3];
rz(-1.6349941) q[3];
sx q[3];
rz(2.5776598) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.22275816) q[2];
sx q[2];
rz(-1.6901878) q[2];
sx q[2];
rz(-1.9082327) q[2];
rz(0.90138609) q[3];
sx q[3];
rz(-3.021535) q[3];
sx q[3];
rz(1.6433158) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0937061) q[0];
sx q[0];
rz(-0.77195764) q[0];
sx q[0];
rz(-0.023660252) q[0];
rz(2.1854782) q[1];
sx q[1];
rz(-1.8319943) q[1];
sx q[1];
rz(0.67217174) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2490847) q[0];
sx q[0];
rz(-1.5810284) q[0];
sx q[0];
rz(-0.0052878629) q[0];
x q[1];
rz(-0.51151885) q[2];
sx q[2];
rz(-1.5788955) q[2];
sx q[2];
rz(0.17009232) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.086239554) q[1];
sx q[1];
rz(-1.2330016) q[1];
sx q[1];
rz(-2.4243381) q[1];
rz(-pi) q[2];
x q[2];
rz(0.43283312) q[3];
sx q[3];
rz(-1.7107309) q[3];
sx q[3];
rz(0.39633745) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.51222926) q[2];
sx q[2];
rz(-1.2269292) q[2];
sx q[2];
rz(-0.36995861) q[2];
rz(1.6379179) q[3];
sx q[3];
rz(-2.2556997) q[3];
sx q[3];
rz(-1.2009719) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5621915) q[0];
sx q[0];
rz(-2.777522) q[0];
sx q[0];
rz(1.2072442) q[0];
rz(-0.72369408) q[1];
sx q[1];
rz(-2.1543398) q[1];
sx q[1];
rz(2.2347246) q[1];
rz(-1.1719218) q[2];
sx q[2];
rz(-1.7181859) q[2];
sx q[2];
rz(-2.0167375) q[2];
rz(-1.1643812) q[3];
sx q[3];
rz(-1.8462528) q[3];
sx q[3];
rz(0.13312199) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
