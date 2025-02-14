OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.7527591) q[0];
sx q[0];
rz(-1.4491117) q[0];
sx q[0];
rz(-1.6911072) q[0];
rz(0.9736355) q[1];
sx q[1];
rz(-1.7042301) q[1];
sx q[1];
rz(-0.91926423) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.60828078) q[0];
sx q[0];
rz(-1.601378) q[0];
sx q[0];
rz(1.5879059) q[0];
rz(-pi) q[1];
rz(0.64627846) q[2];
sx q[2];
rz(-0.15028223) q[2];
sx q[2];
rz(-2.8719547) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.8074933) q[1];
sx q[1];
rz(-1.4791282) q[1];
sx q[1];
rz(2.2794072) q[1];
x q[2];
rz(1.1000865) q[3];
sx q[3];
rz(-1.19095) q[3];
sx q[3];
rz(-2.0947411) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.3322525) q[2];
sx q[2];
rz(-1.6925749) q[2];
sx q[2];
rz(-1.7285041) q[2];
rz(-2.938802) q[3];
sx q[3];
rz(-1.7594124) q[3];
sx q[3];
rz(0.10281674) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24421144) q[0];
sx q[0];
rz(-1.0845217) q[0];
sx q[0];
rz(-0.71075034) q[0];
rz(-2.3731025) q[1];
sx q[1];
rz(-1.0667421) q[1];
sx q[1];
rz(-2.1220727) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7026414) q[0];
sx q[0];
rz(-2.0176689) q[0];
sx q[0];
rz(-1.6146445) q[0];
rz(-2.0020194) q[2];
sx q[2];
rz(-1.6399472) q[2];
sx q[2];
rz(-0.56996843) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.4162035) q[1];
sx q[1];
rz(-2.7692911) q[1];
sx q[1];
rz(-1.6561693) q[1];
rz(-1.016576) q[3];
sx q[3];
rz(-1.6416405) q[3];
sx q[3];
rz(-1.6984102) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.3780313) q[2];
sx q[2];
rz(-1.8969994) q[2];
sx q[2];
rz(1.2223318) q[2];
rz(-1.9289121) q[3];
sx q[3];
rz(-1.0168394) q[3];
sx q[3];
rz(-0.71162629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.057230495) q[0];
sx q[0];
rz(-0.98452345) q[0];
sx q[0];
rz(1.2179751) q[0];
rz(2.6257264) q[1];
sx q[1];
rz(-2.5936544) q[1];
sx q[1];
rz(2.2191494) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2005916) q[0];
sx q[0];
rz(-0.070644826) q[0];
sx q[0];
rz(2.5289422) q[0];
rz(2.416196) q[2];
sx q[2];
rz(-1.9624406) q[2];
sx q[2];
rz(1.2466516) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.45097414) q[1];
sx q[1];
rz(-0.96267525) q[1];
sx q[1];
rz(-2.4324424) q[1];
rz(-pi) q[2];
rz(0.46307989) q[3];
sx q[3];
rz(-0.36470098) q[3];
sx q[3];
rz(-1.1187584) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.1232274) q[2];
sx q[2];
rz(-1.0845228) q[2];
sx q[2];
rz(1.8017192) q[2];
rz(2.5943622) q[3];
sx q[3];
rz(-0.42357835) q[3];
sx q[3];
rz(-2.4303998) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.055534) q[0];
sx q[0];
rz(-2.9965897) q[0];
sx q[0];
rz(-0.15922971) q[0];
rz(3.1314462) q[1];
sx q[1];
rz(-2.1187014) q[1];
sx q[1];
rz(-0.25746447) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2672294) q[0];
sx q[0];
rz(-0.82634514) q[0];
sx q[0];
rz(-2.6494725) q[0];
rz(-pi) q[1];
rz(3.0765216) q[2];
sx q[2];
rz(-0.79539585) q[2];
sx q[2];
rz(0.8000904) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.3336948) q[1];
sx q[1];
rz(-2.195524) q[1];
sx q[1];
rz(-1.6935392) q[1];
rz(-pi) q[2];
rz(0.61473989) q[3];
sx q[3];
rz(-0.90723824) q[3];
sx q[3];
rz(3.1082982) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.7633535) q[2];
sx q[2];
rz(-1.8469609) q[2];
sx q[2];
rz(-2.6687458) q[2];
rz(-0.7044479) q[3];
sx q[3];
rz(-1.7770146) q[3];
sx q[3];
rz(2.4956467) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5064297) q[0];
sx q[0];
rz(-1.1744873) q[0];
sx q[0];
rz(3.0761062) q[0];
rz(-0.40924117) q[1];
sx q[1];
rz(-2.0114653) q[1];
sx q[1];
rz(-1.4261036) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1714892) q[0];
sx q[0];
rz(-1.7894151) q[0];
sx q[0];
rz(-1.6660369) q[0];
x q[1];
rz(0.21435301) q[2];
sx q[2];
rz(-0.39696908) q[2];
sx q[2];
rz(3.1052542) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.1115426) q[1];
sx q[1];
rz(-0.36007127) q[1];
sx q[1];
rz(-2.4637725) q[1];
rz(-2.9687443) q[3];
sx q[3];
rz(-1.2842872) q[3];
sx q[3];
rz(-2.2423173) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.19491974) q[2];
sx q[2];
rz(-1.0871004) q[2];
sx q[2];
rz(-1.3067513) q[2];
rz(-1.170018) q[3];
sx q[3];
rz(-0.40478671) q[3];
sx q[3];
rz(-0.10725966) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65866798) q[0];
sx q[0];
rz(-1.985745) q[0];
sx q[0];
rz(0.40147993) q[0];
rz(-0.67277706) q[1];
sx q[1];
rz(-0.79877001) q[1];
sx q[1];
rz(2.2875517) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14388785) q[0];
sx q[0];
rz(-1.524462) q[0];
sx q[0];
rz(1.7884939) q[0];
rz(-pi) q[1];
rz(1.4208904) q[2];
sx q[2];
rz(-1.3873563) q[2];
sx q[2];
rz(-2.4114087) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.2931765) q[1];
sx q[1];
rz(-2.5418315) q[1];
sx q[1];
rz(2.9648215) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5402732) q[3];
sx q[3];
rz(-2.4375705) q[3];
sx q[3];
rz(0.46425925) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.65471571) q[2];
sx q[2];
rz(-0.87149039) q[2];
sx q[2];
rz(-1.1775449) q[2];
rz(-2.5901637) q[3];
sx q[3];
rz(-1.5827554) q[3];
sx q[3];
rz(0.83561713) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42171445) q[0];
sx q[0];
rz(-2.8604909) q[0];
sx q[0];
rz(1.9792492) q[0];
rz(1.989919) q[1];
sx q[1];
rz(-1.3538227) q[1];
sx q[1];
rz(2.2231359) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7699444) q[0];
sx q[0];
rz(-1.838883) q[0];
sx q[0];
rz(1.9140585) q[0];
rz(0.8753885) q[2];
sx q[2];
rz(-1.3686485) q[2];
sx q[2];
rz(3.1234158) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.20275252) q[1];
sx q[1];
rz(-1.563464) q[1];
sx q[1];
rz(-1.5775024) q[1];
rz(-1.2270532) q[3];
sx q[3];
rz(-1.4101646) q[3];
sx q[3];
rz(-2.4317222) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-3.1071757) q[2];
sx q[2];
rz(-1.7566661) q[2];
sx q[2];
rz(2.6386063) q[2];
rz(0.77477396) q[3];
sx q[3];
rz(-0.46190244) q[3];
sx q[3];
rz(1.3431965) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69304943) q[0];
sx q[0];
rz(-0.22668426) q[0];
sx q[0];
rz(-1.5013129) q[0];
rz(2.8003108) q[1];
sx q[1];
rz(-1.4309859) q[1];
sx q[1];
rz(-1.1192082) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6014746) q[0];
sx q[0];
rz(-2.1261423) q[0];
sx q[0];
rz(2.0106273) q[0];
x q[1];
rz(-2.0430492) q[2];
sx q[2];
rz(-2.611428) q[2];
sx q[2];
rz(1.7165754) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.71257617) q[1];
sx q[1];
rz(-1.6916995) q[1];
sx q[1];
rz(-0.6625794) q[1];
rz(1.4358892) q[3];
sx q[3];
rz(-1.9622318) q[3];
sx q[3];
rz(-1.6329444) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.1559653) q[2];
sx q[2];
rz(-0.95981821) q[2];
sx q[2];
rz(-0.36925527) q[2];
rz(2.8333832) q[3];
sx q[3];
rz(-2.1741368) q[3];
sx q[3];
rz(1.5306028) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8592598) q[0];
sx q[0];
rz(-1.3066602) q[0];
sx q[0];
rz(0.34307137) q[0];
rz(-1.169091) q[1];
sx q[1];
rz(-1.7770551) q[1];
sx q[1];
rz(1.7128568) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0051826) q[0];
sx q[0];
rz(-2.6652968) q[0];
sx q[0];
rz(0.57203697) q[0];
rz(-pi) q[1];
rz(2.3071204) q[2];
sx q[2];
rz(-2.1205466) q[2];
sx q[2];
rz(0.85981926) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.6862168) q[1];
sx q[1];
rz(-1.9460524) q[1];
sx q[1];
rz(-3.0757559) q[1];
rz(-pi) q[2];
rz(0.47023021) q[3];
sx q[3];
rz(-0.93358002) q[3];
sx q[3];
rz(-1.5643844) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.2009361) q[2];
sx q[2];
rz(-0.20874615) q[2];
sx q[2];
rz(1.7500056) q[2];
rz(2.7348147) q[3];
sx q[3];
rz(-1.4429561) q[3];
sx q[3];
rz(-2.2627635) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10130356) q[0];
sx q[0];
rz(-2.5387796) q[0];
sx q[0];
rz(-1.3579177) q[0];
rz(-1.9039924) q[1];
sx q[1];
rz(-2.1289181) q[1];
sx q[1];
rz(0.79992574) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42521503) q[0];
sx q[0];
rz(-0.82934643) q[0];
sx q[0];
rz(0.31405507) q[0];
rz(-pi) q[1];
rz(-1.8736035) q[2];
sx q[2];
rz(-1.0410415) q[2];
sx q[2];
rz(0.15632665) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.0837896) q[1];
sx q[1];
rz(-0.68843319) q[1];
sx q[1];
rz(-2.8958984) q[1];
x q[2];
rz(1.0855476) q[3];
sx q[3];
rz(-1.7929412) q[3];
sx q[3];
rz(2.9911161) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.37821975) q[2];
sx q[2];
rz(-1.7174481) q[2];
sx q[2];
rz(0.073089449) q[2];
rz(-0.25589219) q[3];
sx q[3];
rz(-2.3292694) q[3];
sx q[3];
rz(-1.6528486) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8681317) q[0];
sx q[0];
rz(-1.4556226) q[0];
sx q[0];
rz(1.8709394) q[0];
rz(-1.3399667) q[1];
sx q[1];
rz(-1.7908962) q[1];
sx q[1];
rz(-2.4790196) q[1];
rz(-0.70941464) q[2];
sx q[2];
rz(-1.569289) q[2];
sx q[2];
rz(-2.8808165) q[2];
rz(-1.003391) q[3];
sx q[3];
rz(-1.0091253) q[3];
sx q[3];
rz(-0.71970018) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
