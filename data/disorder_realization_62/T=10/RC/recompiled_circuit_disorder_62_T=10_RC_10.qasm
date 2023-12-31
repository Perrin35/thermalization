OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.2405038) q[0];
sx q[0];
rz(-2.9641889) q[0];
sx q[0];
rz(2.0071964) q[0];
rz(-1.9534684) q[1];
sx q[1];
rz(-1.0367353) q[1];
sx q[1];
rz(-2.477975) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7952607) q[0];
sx q[0];
rz(-1.5471317) q[0];
sx q[0];
rz(-0.46121009) q[0];
rz(-pi) q[1];
rz(1.2167591) q[2];
sx q[2];
rz(-2.3166222) q[2];
sx q[2];
rz(1.1302469) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.2634969) q[1];
sx q[1];
rz(-0.43259183) q[1];
sx q[1];
rz(2.2357975) q[1];
x q[2];
rz(0.10856467) q[3];
sx q[3];
rz(-1.8130455) q[3];
sx q[3];
rz(3.0755175) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.2628281) q[2];
sx q[2];
rz(-0.44439134) q[2];
sx q[2];
rz(0.051068548) q[2];
rz(-0.55705327) q[3];
sx q[3];
rz(-2.3414108) q[3];
sx q[3];
rz(1.5548271) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.58650815) q[0];
sx q[0];
rz(-2.3095135) q[0];
sx q[0];
rz(-0.59659514) q[0];
rz(-2.3157628) q[1];
sx q[1];
rz(-1.4412216) q[1];
sx q[1];
rz(-1.9155496) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3292424) q[0];
sx q[0];
rz(-1.3717522) q[0];
sx q[0];
rz(-2.6890486) q[0];
rz(-pi) q[1];
x q[1];
rz(0.77483564) q[2];
sx q[2];
rz(-2.2171387) q[2];
sx q[2];
rz(-1.6006084) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.021124161) q[1];
sx q[1];
rz(-1.7906584) q[1];
sx q[1];
rz(3.1411509) q[1];
rz(-1.5854884) q[3];
sx q[3];
rz(-0.5726632) q[3];
sx q[3];
rz(-3.1327914) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.79919672) q[2];
sx q[2];
rz(-1.9689955) q[2];
sx q[2];
rz(0.33102316) q[2];
rz(-0.80667574) q[3];
sx q[3];
rz(-1.2253864) q[3];
sx q[3];
rz(1.7139009) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1598635) q[0];
sx q[0];
rz(-1.9112497) q[0];
sx q[0];
rz(-1.8925517) q[0];
rz(3.0535835) q[1];
sx q[1];
rz(-2.0188589) q[1];
sx q[1];
rz(2.1121315) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2217076) q[0];
sx q[0];
rz(-1.789933) q[0];
sx q[0];
rz(0.94421454) q[0];
rz(-pi) q[1];
rz(0.37000334) q[2];
sx q[2];
rz(-0.68379935) q[2];
sx q[2];
rz(-1.5918658) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.50476915) q[1];
sx q[1];
rz(-1.7261337) q[1];
sx q[1];
rz(1.7117281) q[1];
x q[2];
rz(0.13266487) q[3];
sx q[3];
rz(-2.3972315) q[3];
sx q[3];
rz(-0.45385195) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.007894667) q[2];
sx q[2];
rz(-1.4174613) q[2];
sx q[2];
rz(2.6339445) q[2];
rz(1.3890022) q[3];
sx q[3];
rz(-2.8116083) q[3];
sx q[3];
rz(-1.9784137) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4841109) q[0];
sx q[0];
rz(-0.27814516) q[0];
sx q[0];
rz(1.5959651) q[0];
rz(1.0428628) q[1];
sx q[1];
rz(-1.9680126) q[1];
sx q[1];
rz(1.5159336) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.078331) q[0];
sx q[0];
rz(-2.445979) q[0];
sx q[0];
rz(-0.22380933) q[0];
rz(-pi) q[1];
rz(-0.085725351) q[2];
sx q[2];
rz(-1.013373) q[2];
sx q[2];
rz(-2.2103708) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.3486466) q[1];
sx q[1];
rz(-0.9325087) q[1];
sx q[1];
rz(-1.2299041) q[1];
rz(-0.58873119) q[3];
sx q[3];
rz(-1.3874467) q[3];
sx q[3];
rz(-2.6453032) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.3778014) q[2];
sx q[2];
rz(-0.8447454) q[2];
sx q[2];
rz(1.2949004) q[2];
rz(0.14136782) q[3];
sx q[3];
rz(-2.5989792) q[3];
sx q[3];
rz(0.98658371) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(2.7871053) q[0];
sx q[0];
rz(-1.1802477) q[0];
sx q[0];
rz(1.6217344) q[0];
rz(2.5095818) q[1];
sx q[1];
rz(-2.4193587) q[1];
sx q[1];
rz(0.89486665) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.97840727) q[0];
sx q[0];
rz(-0.95046959) q[0];
sx q[0];
rz(1.0582256) q[0];
rz(-pi) q[1];
rz(-2.1331482) q[2];
sx q[2];
rz(-0.92698669) q[2];
sx q[2];
rz(1.84562) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.85256165) q[1];
sx q[1];
rz(-0.82419306) q[1];
sx q[1];
rz(2.186071) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5785061) q[3];
sx q[3];
rz(-1.2843411) q[3];
sx q[3];
rz(-0.25659284) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.52508369) q[2];
sx q[2];
rz(-2.775165) q[2];
sx q[2];
rz(-1.9990702) q[2];
rz(0.74674314) q[3];
sx q[3];
rz(-1.8871504) q[3];
sx q[3];
rz(2.0660627) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.58364761) q[0];
sx q[0];
rz(-2.8201411) q[0];
sx q[0];
rz(1.7640132) q[0];
rz(-0.47239834) q[1];
sx q[1];
rz(-0.51858416) q[1];
sx q[1];
rz(-2.6766052) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.87100055) q[0];
sx q[0];
rz(-2.5809079) q[0];
sx q[0];
rz(-2.2339348) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.497252) q[2];
sx q[2];
rz(-1.2300756) q[2];
sx q[2];
rz(1.7433012) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.38575129) q[1];
sx q[1];
rz(-2.5045536) q[1];
sx q[1];
rz(1.6511276) q[1];
x q[2];
rz(-2.898071) q[3];
sx q[3];
rz(-2.4176819) q[3];
sx q[3];
rz(0.078725423) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.49089367) q[2];
sx q[2];
rz(-1.2630254) q[2];
sx q[2];
rz(2.6965551) q[2];
rz(-0.93368357) q[3];
sx q[3];
rz(-1.7088339) q[3];
sx q[3];
rz(2.8745108) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0141107) q[0];
sx q[0];
rz(-1.575379) q[0];
sx q[0];
rz(-1.4200462) q[0];
rz(3.1177915) q[1];
sx q[1];
rz(-0.61444608) q[1];
sx q[1];
rz(0.15596095) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6316846) q[0];
sx q[0];
rz(-1.8242867) q[0];
sx q[0];
rz(-1.3939199) q[0];
x q[1];
rz(-1.1440802) q[2];
sx q[2];
rz(-2.0442171) q[2];
sx q[2];
rz(1.13525) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.1289039) q[1];
sx q[1];
rz(-1.3959937) q[1];
sx q[1];
rz(0.054794475) q[1];
rz(-pi) q[2];
rz(0.74219269) q[3];
sx q[3];
rz(-2.4626197) q[3];
sx q[3];
rz(-1.9051444) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-3.0722787) q[2];
sx q[2];
rz(-0.57702714) q[2];
sx q[2];
rz(-2.0689266) q[2];
rz(-2.8159451) q[3];
sx q[3];
rz(-1.9653392) q[3];
sx q[3];
rz(-1.5163039) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
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
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9496562) q[0];
sx q[0];
rz(-2.7392116) q[0];
sx q[0];
rz(-2.8038213) q[0];
rz(1.0900963) q[1];
sx q[1];
rz(-2.1665159) q[1];
sx q[1];
rz(-0.24857323) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2598495) q[0];
sx q[0];
rz(-1.8200841) q[0];
sx q[0];
rz(-0.52815204) q[0];
rz(-pi) q[1];
rz(1.0609264) q[2];
sx q[2];
rz(-1.14398) q[2];
sx q[2];
rz(-2.5271497) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.3358811) q[1];
sx q[1];
rz(-0.98166775) q[1];
sx q[1];
rz(-0.59405234) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9931273) q[3];
sx q[3];
rz(-0.18581192) q[3];
sx q[3];
rz(-1.8728135) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.2934072) q[2];
sx q[2];
rz(-0.53331393) q[2];
sx q[2];
rz(3.0548813) q[2];
rz(2.6596206) q[3];
sx q[3];
rz(-2.635699) q[3];
sx q[3];
rz(-1.1530676) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6329353) q[0];
sx q[0];
rz(-1.0077227) q[0];
sx q[0];
rz(-2.8588262) q[0];
rz(0.70156082) q[1];
sx q[1];
rz(-0.82156721) q[1];
sx q[1];
rz(1.823002) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9160999) q[0];
sx q[0];
rz(-2.6927104) q[0];
sx q[0];
rz(1.3845969) q[0];
rz(-pi) q[1];
rz(-2.4852072) q[2];
sx q[2];
rz(-2.0771386) q[2];
sx q[2];
rz(0.93751794) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.56837592) q[1];
sx q[1];
rz(-1.7163367) q[1];
sx q[1];
rz(-0.88370609) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.41140822) q[3];
sx q[3];
rz(-0.79380006) q[3];
sx q[3];
rz(-2.0421162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.41137722) q[2];
sx q[2];
rz(-2.8432379) q[2];
sx q[2];
rz(-2.7424157) q[2];
rz(-0.88360751) q[3];
sx q[3];
rz(-1.4727605) q[3];
sx q[3];
rz(1.9201027) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3783962) q[0];
sx q[0];
rz(-2.3662687) q[0];
sx q[0];
rz(1.5426853) q[0];
rz(2.0762766) q[1];
sx q[1];
rz(-2.1677446) q[1];
sx q[1];
rz(1.261196) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3189145) q[0];
sx q[0];
rz(-0.42352522) q[0];
sx q[0];
rz(-0.25869297) q[0];
rz(-pi) q[1];
rz(-0.032632685) q[2];
sx q[2];
rz(-2.543078) q[2];
sx q[2];
rz(1.061071) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.0295769) q[1];
sx q[1];
rz(-2.638991) q[1];
sx q[1];
rz(3.0081248) q[1];
rz(-pi) q[2];
rz(2.946978) q[3];
sx q[3];
rz(-2.448304) q[3];
sx q[3];
rz(2.4926536) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.5836872) q[2];
sx q[2];
rz(-2.5186899) q[2];
sx q[2];
rz(-2.5718001) q[2];
rz(-1.9231046) q[3];
sx q[3];
rz(-0.64703882) q[3];
sx q[3];
rz(-2.5789554) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31310836) q[0];
sx q[0];
rz(-2.2611571) q[0];
sx q[0];
rz(-1.8631998) q[0];
rz(2.5333511) q[1];
sx q[1];
rz(-0.47641644) q[1];
sx q[1];
rz(0.48412916) q[1];
rz(0.99909487) q[2];
sx q[2];
rz(-1.5928762) q[2];
sx q[2];
rz(1.5230509) q[2];
rz(0.20523397) q[3];
sx q[3];
rz(-2.5406465) q[3];
sx q[3];
rz(3.061486) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
