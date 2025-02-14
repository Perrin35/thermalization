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
rz(-0.69775692) q[0];
sx q[0];
rz(-1.193576) q[0];
sx q[0];
rz(-0.20456631) q[0];
rz(-0.76072955) q[1];
sx q[1];
rz(-1.3890356) q[1];
sx q[1];
rz(1.3995481) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7273063) q[0];
sx q[0];
rz(-2.901863) q[0];
sx q[0];
rz(1.1778465) q[0];
rz(1.1209773) q[2];
sx q[2];
rz(-1.9391141) q[2];
sx q[2];
rz(0.30003795) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.5198361) q[1];
sx q[1];
rz(-0.43070983) q[1];
sx q[1];
rz(-0.21060305) q[1];
x q[2];
rz(0.30866039) q[3];
sx q[3];
rz(-2.4081569) q[3];
sx q[3];
rz(-1.6876458) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.7064887) q[2];
sx q[2];
rz(-0.032328345) q[2];
sx q[2];
rz(2.6522563) q[2];
rz(2.1183744) q[3];
sx q[3];
rz(-3.1232941) q[3];
sx q[3];
rz(2.0160915) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0959051) q[0];
sx q[0];
rz(-2.4850595) q[0];
sx q[0];
rz(-2.3387961) q[0];
rz(-0.071391694) q[1];
sx q[1];
rz(-2.8749021) q[1];
sx q[1];
rz(0.057770483) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.18081576) q[0];
sx q[0];
rz(-1.0313927) q[0];
sx q[0];
rz(0.30110724) q[0];
x q[1];
rz(-1.6768912) q[2];
sx q[2];
rz(-1.2857784) q[2];
sx q[2];
rz(0.34044701) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.3398185) q[1];
sx q[1];
rz(-1.8027824) q[1];
sx q[1];
rz(-2.0663227) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.14795078) q[3];
sx q[3];
rz(-0.86455621) q[3];
sx q[3];
rz(2.1220292) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.1959261) q[2];
sx q[2];
rz(-2.0464996) q[2];
sx q[2];
rz(-1.8460974) q[2];
rz(-0.96674353) q[3];
sx q[3];
rz(-2.3714378) q[3];
sx q[3];
rz(-2.3841592) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1118689) q[0];
sx q[0];
rz(-1.3628549) q[0];
sx q[0];
rz(-1.7080074) q[0];
rz(0.068610527) q[1];
sx q[1];
rz(-1.5674633) q[1];
sx q[1];
rz(0.57919085) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9703909) q[0];
sx q[0];
rz(-0.10848898) q[0];
sx q[0];
rz(-1.0068588) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5689108) q[2];
sx q[2];
rz(-2.0465188) q[2];
sx q[2];
rz(1.204551) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.75603682) q[1];
sx q[1];
rz(-1.3446523) q[1];
sx q[1];
rz(-1.6075587) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0924358) q[3];
sx q[3];
rz(-1.1598827) q[3];
sx q[3];
rz(2.876407) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.8706943) q[2];
sx q[2];
rz(-1.8879994) q[2];
sx q[2];
rz(-2.9581621) q[2];
rz(0.80426788) q[3];
sx q[3];
rz(-0.99884123) q[3];
sx q[3];
rz(2.6075294) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9428228) q[0];
sx q[0];
rz(-0.11480055) q[0];
sx q[0];
rz(-2.542069) q[0];
rz(0.41246688) q[1];
sx q[1];
rz(-3.1209374) q[1];
sx q[1];
rz(-2.1309158) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28579564) q[0];
sx q[0];
rz(-1.2210346) q[0];
sx q[0];
rz(-1.4361022) q[0];
rz(-3.0907107) q[2];
sx q[2];
rz(-2.2023099) q[2];
sx q[2];
rz(0.33755195) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.49941555) q[1];
sx q[1];
rz(-0.66153971) q[1];
sx q[1];
rz(-2.0209577) q[1];
rz(-2.385731) q[3];
sx q[3];
rz(-1.469765) q[3];
sx q[3];
rz(1.3974691) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.7708873) q[2];
sx q[2];
rz(-0.34660307) q[2];
sx q[2];
rz(0.31751219) q[2];
rz(-0.62234771) q[3];
sx q[3];
rz(-2.205866) q[3];
sx q[3];
rz(-0.58421016) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36251003) q[0];
sx q[0];
rz(-2.2267987) q[0];
sx q[0];
rz(-1.0269748) q[0];
rz(-0.55038553) q[1];
sx q[1];
rz(-3.0774979) q[1];
sx q[1];
rz(1.9245573) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37792045) q[0];
sx q[0];
rz(-2.2224713) q[0];
sx q[0];
rz(2.4125189) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6486859) q[2];
sx q[2];
rz(-2.1078034) q[2];
sx q[2];
rz(0.80918771) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.1149619) q[1];
sx q[1];
rz(-1.3991881) q[1];
sx q[1];
rz(2.2725355) q[1];
rz(-pi) q[2];
rz(-0.86037614) q[3];
sx q[3];
rz(-1.851322) q[3];
sx q[3];
rz(2.4376412) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.43651849) q[2];
sx q[2];
rz(-1.3631835) q[2];
sx q[2];
rz(-0.77318937) q[2];
rz(-3.0298722) q[3];
sx q[3];
rz(-1.819928) q[3];
sx q[3];
rz(-0.93938655) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6620827) q[0];
sx q[0];
rz(-0.22308068) q[0];
sx q[0];
rz(2.7146085) q[0];
rz(-0.93049479) q[1];
sx q[1];
rz(-0.016914802) q[1];
sx q[1];
rz(0.46447909) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1982959) q[0];
sx q[0];
rz(-1.7664599) q[0];
sx q[0];
rz(-0.29384675) q[0];
rz(-pi) q[1];
x q[1];
rz(0.33916766) q[2];
sx q[2];
rz(-0.99771032) q[2];
sx q[2];
rz(-3.1011875) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.528094) q[1];
sx q[1];
rz(-1.534675) q[1];
sx q[1];
rz(-0.36904676) q[1];
x q[2];
rz(-3.1103856) q[3];
sx q[3];
rz(-1.9949119) q[3];
sx q[3];
rz(-1.1348789) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.6102607) q[2];
sx q[2];
rz(-1.320763) q[2];
sx q[2];
rz(0.28826928) q[2];
rz(2.098295) q[3];
sx q[3];
rz(-2.5259924) q[3];
sx q[3];
rz(2.402795) q[3];
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
x q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6601324) q[0];
sx q[0];
rz(-1.2876502) q[0];
sx q[0];
rz(-2.3840391) q[0];
rz(3.0689012) q[1];
sx q[1];
rz(-3.1158267) q[1];
sx q[1];
rz(-0.048197897) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2295727) q[0];
sx q[0];
rz(-0.98668146) q[0];
sx q[0];
rz(0.035365625) q[0];
rz(-pi) q[1];
x q[1];
rz(2.169007) q[2];
sx q[2];
rz(-1.3763104) q[2];
sx q[2];
rz(2.011428) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.1485702) q[1];
sx q[1];
rz(-2.3766368) q[1];
sx q[1];
rz(2.2123442) q[1];
x q[2];
rz(2.0406575) q[3];
sx q[3];
rz(-1.5239232) q[3];
sx q[3];
rz(0.27134233) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.4797719) q[2];
sx q[2];
rz(-1.6419819) q[2];
sx q[2];
rz(-0.069395937) q[2];
rz(-1.5540468) q[3];
sx q[3];
rz(-2.3543251) q[3];
sx q[3];
rz(2.9230996) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3875535) q[0];
sx q[0];
rz(-2.0856922) q[0];
sx q[0];
rz(1.7283424) q[0];
rz(2.3130401) q[1];
sx q[1];
rz(-3.0998402) q[1];
sx q[1];
rz(0.54263306) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4601645) q[0];
sx q[0];
rz(-1.6896473) q[0];
sx q[0];
rz(-0.14560856) q[0];
x q[1];
rz(-1.2130402) q[2];
sx q[2];
rz(-1.925549) q[2];
sx q[2];
rz(-2.8194129) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.16281505) q[1];
sx q[1];
rz(-3.0138624) q[1];
sx q[1];
rz(0.43600299) q[1];
rz(-pi) q[2];
rz(1.1922791) q[3];
sx q[3];
rz(-1.9782441) q[3];
sx q[3];
rz(0.70127869) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.98216206) q[2];
sx q[2];
rz(-2.7320778) q[2];
sx q[2];
rz(0.25992599) q[2];
rz(-2.121117) q[3];
sx q[3];
rz(-2.8838938) q[3];
sx q[3];
rz(0.79403383) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6771773) q[0];
sx q[0];
rz(-0.18615119) q[0];
sx q[0];
rz(1.4790685) q[0];
rz(1.6089449) q[1];
sx q[1];
rz(-2.1023991) q[1];
sx q[1];
rz(0.74302465) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4585658) q[0];
sx q[0];
rz(-1.4560149) q[0];
sx q[0];
rz(-2.0857138) q[0];
x q[1];
rz(0.3826999) q[2];
sx q[2];
rz(-1.0596794) q[2];
sx q[2];
rz(-1.9267043) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.0155269) q[1];
sx q[1];
rz(-0.085745009) q[1];
sx q[1];
rz(2.2669906) q[1];
rz(-pi) q[2];
rz(0.10492341) q[3];
sx q[3];
rz(-1.5399974) q[3];
sx q[3];
rz(-1.2237751) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.6741901) q[2];
sx q[2];
rz(-0.82618606) q[2];
sx q[2];
rz(-2.3361333) q[2];
rz(-1.6921267) q[3];
sx q[3];
rz(-1.9129246) q[3];
sx q[3];
rz(-0.8031556) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2529124) q[0];
sx q[0];
rz(-0.54179931) q[0];
sx q[0];
rz(2.3895277) q[0];
rz(-1.9883142) q[1];
sx q[1];
rz(-2.2572932) q[1];
sx q[1];
rz(2.8582252) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3842938) q[0];
sx q[0];
rz(-1.6597676) q[0];
sx q[0];
rz(2.6020537) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.14788515) q[2];
sx q[2];
rz(-2.3171632) q[2];
sx q[2];
rz(1.7764719) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.4979889) q[1];
sx q[1];
rz(-2.170553) q[1];
sx q[1];
rz(1.1685755) q[1];
x q[2];
rz(2.649077) q[3];
sx q[3];
rz(-2.4253143) q[3];
sx q[3];
rz(3.1066343) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.2470384) q[2];
sx q[2];
rz(-0.082823195) q[2];
sx q[2];
rz(1.7029597) q[2];
rz(-2.8476207) q[3];
sx q[3];
rz(-3.1271264) q[3];
sx q[3];
rz(-2.1132052) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.033584874) q[0];
sx q[0];
rz(-1.4172194) q[0];
sx q[0];
rz(-1.5237756) q[0];
rz(0.53957466) q[1];
sx q[1];
rz(-0.78782606) q[1];
sx q[1];
rz(0.16001564) q[1];
rz(-2.9931184) q[2];
sx q[2];
rz(-1.9273026) q[2];
sx q[2];
rz(-2.9782563) q[2];
rz(-2.7593437) q[3];
sx q[3];
rz(-1.3180238) q[3];
sx q[3];
rz(-2.870306) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
