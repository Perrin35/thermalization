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
rz(-0.18345565) q[0];
sx q[0];
rz(-0.7440716) q[0];
sx q[0];
rz(-1.0519354) q[0];
rz(0.19368859) q[1];
sx q[1];
rz(-0.63313484) q[1];
sx q[1];
rz(0.57300353) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5405226) q[0];
sx q[0];
rz(-0.94694505) q[0];
sx q[0];
rz(-3.0649351) q[0];
rz(2.6296205) q[2];
sx q[2];
rz(-0.7535156) q[2];
sx q[2];
rz(1.7640424) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.51542789) q[1];
sx q[1];
rz(-2.0671131) q[1];
sx q[1];
rz(1.4908916) q[1];
rz(-pi) q[2];
rz(-1.4680193) q[3];
sx q[3];
rz(-2.2647018) q[3];
sx q[3];
rz(-0.68460195) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.28213349) q[2];
sx q[2];
rz(-2.8318475) q[2];
sx q[2];
rz(-1.0563043) q[2];
rz(-3.1193962) q[3];
sx q[3];
rz(-2.3789417) q[3];
sx q[3];
rz(1.7215884) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9098772) q[0];
sx q[0];
rz(-0.48119369) q[0];
sx q[0];
rz(0.43014446) q[0];
rz(0.12769708) q[1];
sx q[1];
rz(-1.1558497) q[1];
sx q[1];
rz(-1.7040303) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8970222) q[0];
sx q[0];
rz(-1.26957) q[0];
sx q[0];
rz(1.8409077) q[0];
rz(-pi) q[1];
rz(-1.6280074) q[2];
sx q[2];
rz(-1.621641) q[2];
sx q[2];
rz(1.0982996) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.3900657) q[1];
sx q[1];
rz(-1.8987149) q[1];
sx q[1];
rz(-2.4310914) q[1];
x q[2];
rz(2.6019179) q[3];
sx q[3];
rz(-1.7575348) q[3];
sx q[3];
rz(-0.58247551) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.11640707) q[2];
sx q[2];
rz(-1.6802639) q[2];
sx q[2];
rz(-0.44542584) q[2];
rz(0.40924254) q[3];
sx q[3];
rz(-1.0990812) q[3];
sx q[3];
rz(-0.87944952) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5889848) q[0];
sx q[0];
rz(-0.092985066) q[0];
sx q[0];
rz(2.4267922) q[0];
rz(-1.0572761) q[1];
sx q[1];
rz(-0.42066586) q[1];
sx q[1];
rz(0.32726273) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6598073) q[0];
sx q[0];
rz(-2.5707173) q[0];
sx q[0];
rz(-1.2587738) q[0];
x q[1];
rz(-0.94560854) q[2];
sx q[2];
rz(-2.473712) q[2];
sx q[2];
rz(0.57648522) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.63673692) q[1];
sx q[1];
rz(-1.550192) q[1];
sx q[1];
rz(-1.4192753) q[1];
rz(-pi) q[2];
rz(2.8873575) q[3];
sx q[3];
rz(-1.2336858) q[3];
sx q[3];
rz(1.997662) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.0032234) q[2];
sx q[2];
rz(-2.1681163) q[2];
sx q[2];
rz(1.6376015) q[2];
rz(1.2184294) q[3];
sx q[3];
rz(-0.4156433) q[3];
sx q[3];
rz(2.159582) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0729436) q[0];
sx q[0];
rz(-1.2462085) q[0];
sx q[0];
rz(0.15596998) q[0];
rz(-2.7929557) q[1];
sx q[1];
rz(-0.60395423) q[1];
sx q[1];
rz(-1.9085931) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.975978) q[0];
sx q[0];
rz(-2.8535643) q[0];
sx q[0];
rz(-2.0095129) q[0];
rz(-pi) q[1];
rz(0.92278752) q[2];
sx q[2];
rz(-1.7402667) q[2];
sx q[2];
rz(-3.0344935) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.11248842) q[1];
sx q[1];
rz(-1.2891411) q[1];
sx q[1];
rz(-2.2206578) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.057179515) q[3];
sx q[3];
rz(-0.25186497) q[3];
sx q[3];
rz(0.72309031) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.6984581) q[2];
sx q[2];
rz(-0.036018697) q[2];
sx q[2];
rz(-1.6301463) q[2];
rz(1.1394507) q[3];
sx q[3];
rz(-1.3316493) q[3];
sx q[3];
rz(0.99036923) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.776942) q[0];
sx q[0];
rz(-0.41780892) q[0];
sx q[0];
rz(-1.1530217) q[0];
rz(2.5979089) q[1];
sx q[1];
rz(-1.2666356) q[1];
sx q[1];
rz(-2.8797454) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76806289) q[0];
sx q[0];
rz(-0.090001194) q[0];
sx q[0];
rz(1.9530208) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.98723094) q[2];
sx q[2];
rz(-1.2886184) q[2];
sx q[2];
rz(0.63167494) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.10298577) q[1];
sx q[1];
rz(-2.1521882) q[1];
sx q[1];
rz(1.3036779) q[1];
rz(-pi) q[2];
rz(-0.20528593) q[3];
sx q[3];
rz(-0.79278273) q[3];
sx q[3];
rz(2.1403811) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.61949817) q[2];
sx q[2];
rz(-0.60569373) q[2];
sx q[2];
rz(-1.7363133) q[2];
rz(-2.3960579) q[3];
sx q[3];
rz(-1.8757952) q[3];
sx q[3];
rz(-2.1982927) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
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
rz(-3.140542) q[0];
sx q[0];
rz(-3.0241835) q[0];
sx q[0];
rz(-2.7897799) q[0];
rz(-2.70128) q[1];
sx q[1];
rz(-1.7761209) q[1];
sx q[1];
rz(2.9702759) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.139411) q[0];
sx q[0];
rz(-1.8137964) q[0];
sx q[0];
rz(-2.4619106) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2662884) q[2];
sx q[2];
rz(-0.84104462) q[2];
sx q[2];
rz(0.27308057) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.76207668) q[1];
sx q[1];
rz(-0.87240309) q[1];
sx q[1];
rz(-1.5435436) q[1];
rz(-pi) q[2];
rz(-0.8796575) q[3];
sx q[3];
rz(-2.1591957) q[3];
sx q[3];
rz(-0.0042875687) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-3.0773086) q[2];
sx q[2];
rz(-1.7694387) q[2];
sx q[2];
rz(-2.1997931) q[2];
rz(1.9866379) q[3];
sx q[3];
rz(-0.20808163) q[3];
sx q[3];
rz(1.340516) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(0.44523859) q[0];
sx q[0];
rz(-1.1055163) q[0];
sx q[0];
rz(-0.0048986991) q[0];
rz(-0.44662961) q[1];
sx q[1];
rz(-2.547867) q[1];
sx q[1];
rz(0.68797025) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3197927) q[0];
sx q[0];
rz(-1.3988136) q[0];
sx q[0];
rz(0.71806192) q[0];
rz(-pi) q[1];
rz(0.61882682) q[2];
sx q[2];
rz(-1.8558673) q[2];
sx q[2];
rz(0.10933354) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.1157284) q[1];
sx q[1];
rz(-1.1394355) q[1];
sx q[1];
rz(0.23484767) q[1];
rz(-0.87422411) q[3];
sx q[3];
rz(-2.1316574) q[3];
sx q[3];
rz(2.0340597) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.530431) q[2];
sx q[2];
rz(-0.40803424) q[2];
sx q[2];
rz(0.92998695) q[2];
rz(1.2215349) q[3];
sx q[3];
rz(-2.0285716) q[3];
sx q[3];
rz(-2.5482224) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4596443) q[0];
sx q[0];
rz(-1.3197897) q[0];
sx q[0];
rz(0.090106877) q[0];
rz(1.4247165) q[1];
sx q[1];
rz(-0.63980353) q[1];
sx q[1];
rz(0.43509126) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1662707) q[0];
sx q[0];
rz(-1.4251801) q[0];
sx q[0];
rz(-2.6377489) q[0];
rz(1.3899562) q[2];
sx q[2];
rz(-2.2308084) q[2];
sx q[2];
rz(-1.3446128) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.43109801) q[1];
sx q[1];
rz(-0.63023797) q[1];
sx q[1];
rz(1.789854) q[1];
x q[2];
rz(1.5699205) q[3];
sx q[3];
rz(-1.1464351) q[3];
sx q[3];
rz(2.3800338) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.4789751) q[2];
sx q[2];
rz(-2.0765897) q[2];
sx q[2];
rz(-0.62620658) q[2];
rz(0.19691697) q[3];
sx q[3];
rz(-0.99072376) q[3];
sx q[3];
rz(-1.7530493) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.58186746) q[0];
sx q[0];
rz(-1.6904866) q[0];
sx q[0];
rz(0.054280601) q[0];
rz(-0.42298969) q[1];
sx q[1];
rz(-1.3754247) q[1];
sx q[1];
rz(1.3351701) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2847808) q[0];
sx q[0];
rz(-2.3216212) q[0];
sx q[0];
rz(0.7818082) q[0];
rz(-pi) q[1];
rz(-2.1077431) q[2];
sx q[2];
rz(-0.62990153) q[2];
sx q[2];
rz(0.065187188) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.685251) q[1];
sx q[1];
rz(-1.047118) q[1];
sx q[1];
rz(-2.1728553) q[1];
rz(-pi) q[2];
rz(2.0306251) q[3];
sx q[3];
rz(-2.9879192) q[3];
sx q[3];
rz(-0.75017649) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.52428952) q[2];
sx q[2];
rz(-1.6733988) q[2];
sx q[2];
rz(-0.35150251) q[2];
rz(1.7678123) q[3];
sx q[3];
rz(-0.51353729) q[3];
sx q[3];
rz(2.8721749) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3897301) q[0];
sx q[0];
rz(-0.26079145) q[0];
sx q[0];
rz(0.75916284) q[0];
rz(1.0836541) q[1];
sx q[1];
rz(-1.6108797) q[1];
sx q[1];
rz(1.0677451) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40475527) q[0];
sx q[0];
rz(-2.1421297) q[0];
sx q[0];
rz(0.57147567) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.847151) q[2];
sx q[2];
rz(-1.3047555) q[2];
sx q[2];
rz(-2.5584115) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.2764923) q[1];
sx q[1];
rz(-2.4363764) q[1];
sx q[1];
rz(-2.221285) q[1];
x q[2];
rz(0.87369793) q[3];
sx q[3];
rz(-1.3109968) q[3];
sx q[3];
rz(-1.8099305) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.4296253) q[2];
sx q[2];
rz(-1.2099313) q[2];
sx q[2];
rz(-2.9313226) q[2];
rz(2.353904) q[3];
sx q[3];
rz(-1.7588047) q[3];
sx q[3];
rz(2.6113966) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6982211) q[0];
sx q[0];
rz(-0.5589232) q[0];
sx q[0];
rz(1.9236175) q[0];
rz(-1.632985) q[1];
sx q[1];
rz(-1.2651545) q[1];
sx q[1];
rz(-1.9427585) q[1];
rz(2.0532578) q[2];
sx q[2];
rz(-2.3934622) q[2];
sx q[2];
rz(-0.29665034) q[2];
rz(-2.9859424) q[3];
sx q[3];
rz(-1.2990824) q[3];
sx q[3];
rz(-1.6344447) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
