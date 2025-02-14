OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.3396575) q[0];
sx q[0];
rz(-0.8987838) q[0];
sx q[0];
rz(1.3026097) q[0];
rz(0.12401914) q[1];
sx q[1];
rz(4.8766512) q[1];
sx q[1];
rz(7.7968346) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2500656) q[0];
sx q[0];
rz(-2.6882049) q[0];
sx q[0];
rz(-0.84693036) q[0];
rz(-pi) q[1];
x q[1];
rz(0.75998016) q[2];
sx q[2];
rz(-0.70356762) q[2];
sx q[2];
rz(-0.68389308) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.096949654) q[1];
sx q[1];
rz(-1.6253119) q[1];
sx q[1];
rz(1.9115555) q[1];
rz(2.7645993) q[3];
sx q[3];
rz(-0.97929472) q[3];
sx q[3];
rz(-1.5821004) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.0836432) q[2];
sx q[2];
rz(-1.3336072) q[2];
sx q[2];
rz(0.38457125) q[2];
rz(0.40258506) q[3];
sx q[3];
rz(-1.4339002) q[3];
sx q[3];
rz(-0.80818498) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.758051) q[0];
sx q[0];
rz(-1.533968) q[0];
sx q[0];
rz(2.5480399) q[0];
rz(1.8202164) q[1];
sx q[1];
rz(-1.079419) q[1];
sx q[1];
rz(3.1022601) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.93146053) q[0];
sx q[0];
rz(-1.1173741) q[0];
sx q[0];
rz(2.7260145) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4565384) q[2];
sx q[2];
rz(-0.93747497) q[2];
sx q[2];
rz(-2.3267724) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.63168282) q[1];
sx q[1];
rz(-1.0541774) q[1];
sx q[1];
rz(1.5138018) q[1];
x q[2];
rz(-1.6013299) q[3];
sx q[3];
rz(-1.021592) q[3];
sx q[3];
rz(-2.1516678) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.84052691) q[2];
sx q[2];
rz(-1.7037921) q[2];
sx q[2];
rz(-2.5145516) q[2];
rz(2.7035233) q[3];
sx q[3];
rz(-1.6431243) q[3];
sx q[3];
rz(-2.0048678) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42191926) q[0];
sx q[0];
rz(-1.0872343) q[0];
sx q[0];
rz(1.3713974) q[0];
rz(-0.60945359) q[1];
sx q[1];
rz(-0.89027673) q[1];
sx q[1];
rz(-1.9166463) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6951675) q[0];
sx q[0];
rz(-1.7621627) q[0];
sx q[0];
rz(-1.6686977) q[0];
rz(-pi) q[1];
rz(-2.0535499) q[2];
sx q[2];
rz(-0.73763371) q[2];
sx q[2];
rz(2.9575153) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.48771) q[1];
sx q[1];
rz(-1.9284029) q[1];
sx q[1];
rz(-0.34642152) q[1];
rz(-pi) q[2];
rz(0.21594343) q[3];
sx q[3];
rz(-2.2547997) q[3];
sx q[3];
rz(2.0404599) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.8372535) q[2];
sx q[2];
rz(-2.2746268) q[2];
sx q[2];
rz(2.2875817) q[2];
rz(-2.617406) q[3];
sx q[3];
rz(-1.5093404) q[3];
sx q[3];
rz(2.1715651) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9147341) q[0];
sx q[0];
rz(-0.15758841) q[0];
sx q[0];
rz(2.3478813) q[0];
rz(-0.0018250068) q[1];
sx q[1];
rz(-0.55627126) q[1];
sx q[1];
rz(2.3307641) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2398259) q[0];
sx q[0];
rz(-1.3890653) q[0];
sx q[0];
rz(2.316409) q[0];
rz(-pi) q[1];
rz(-0.43528683) q[2];
sx q[2];
rz(-0.94418028) q[2];
sx q[2];
rz(2.1666908) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.8673265) q[1];
sx q[1];
rz(-1.9073448) q[1];
sx q[1];
rz(0.95086581) q[1];
rz(-pi) q[2];
rz(-1.2620737) q[3];
sx q[3];
rz(-2.6938872) q[3];
sx q[3];
rz(-1.6310504) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.68373716) q[2];
sx q[2];
rz(-0.33317864) q[2];
sx q[2];
rz(1.6443171) q[2];
rz(1.3582683) q[3];
sx q[3];
rz(-2.0446916) q[3];
sx q[3];
rz(1.3076521) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
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
rz(-1.7233647) q[0];
sx q[0];
rz(-1.6809373) q[0];
sx q[0];
rz(-0.25890589) q[0];
rz(0.12340165) q[1];
sx q[1];
rz(-0.63107189) q[1];
sx q[1];
rz(2.6554328) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0927057) q[0];
sx q[0];
rz(-0.98711038) q[0];
sx q[0];
rz(2.5539469) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7068549) q[2];
sx q[2];
rz(-1.1216838) q[2];
sx q[2];
rz(0.27829188) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.9100719) q[1];
sx q[1];
rz(-1.9329482) q[1];
sx q[1];
rz(-0.2015743) q[1];
rz(-pi) q[2];
rz(-2.2712689) q[3];
sx q[3];
rz(-1.1253256) q[3];
sx q[3];
rz(1.1401759) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.90998021) q[2];
sx q[2];
rz(-2.6284802) q[2];
sx q[2];
rz(1.7463589) q[2];
rz(0.95476556) q[3];
sx q[3];
rz(-1.3937817) q[3];
sx q[3];
rz(2.098293) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.93916494) q[0];
sx q[0];
rz(-1.9068149) q[0];
sx q[0];
rz(-2.3811316) q[0];
rz(-2.3072534) q[1];
sx q[1];
rz(-1.1015247) q[1];
sx q[1];
rz(-1.1044501) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6532947) q[0];
sx q[0];
rz(-2.1135931) q[0];
sx q[0];
rz(-2.8858868) q[0];
x q[1];
rz(0.34453407) q[2];
sx q[2];
rz(-1.1921368) q[2];
sx q[2];
rz(-1.6123503) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.84127562) q[1];
sx q[1];
rz(-1.0659443) q[1];
sx q[1];
rz(-2.3077479) q[1];
rz(-pi) q[2];
rz(-2.0310294) q[3];
sx q[3];
rz(-0.62556872) q[3];
sx q[3];
rz(-2.9265508) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.96222782) q[2];
sx q[2];
rz(-0.70415512) q[2];
sx q[2];
rz(3.1000225) q[2];
rz(-2.9465594) q[3];
sx q[3];
rz(-2.8988367) q[3];
sx q[3];
rz(-2.354505) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.70260173) q[0];
sx q[0];
rz(-2.2881303) q[0];
sx q[0];
rz(-2.3882197) q[0];
rz(-1.0607177) q[1];
sx q[1];
rz(-2.3371688) q[1];
sx q[1];
rz(2.9341968) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3531137) q[0];
sx q[0];
rz(-0.99035701) q[0];
sx q[0];
rz(-2.3222695) q[0];
x q[1];
rz(-1.8600402) q[2];
sx q[2];
rz(-1.1221894) q[2];
sx q[2];
rz(-1.9236652) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.8589391) q[1];
sx q[1];
rz(-0.6047073) q[1];
sx q[1];
rz(-2.7207123) q[1];
rz(-pi) q[2];
x q[2];
rz(3.1025299) q[3];
sx q[3];
rz(-1.563058) q[3];
sx q[3];
rz(1.4259065) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.9994026) q[2];
sx q[2];
rz(-1.157607) q[2];
sx q[2];
rz(-3.0863777) q[2];
rz(0.84290543) q[3];
sx q[3];
rz(-2.3838145) q[3];
sx q[3];
rz(-2.260476) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.371599) q[0];
sx q[0];
rz(-2.433625) q[0];
sx q[0];
rz(1.9644894) q[0];
rz(2.2616995) q[1];
sx q[1];
rz(-2.8113007) q[1];
sx q[1];
rz(3.0807307) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0165981) q[0];
sx q[0];
rz(-1.6654764) q[0];
sx q[0];
rz(-1.4461317) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6598236) q[2];
sx q[2];
rz(-1.2854115) q[2];
sx q[2];
rz(2.9586359) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.60778516) q[1];
sx q[1];
rz(-2.5609511) q[1];
sx q[1];
rz(2.9279079) q[1];
x q[2];
rz(-3.0677972) q[3];
sx q[3];
rz(-1.3072469) q[3];
sx q[3];
rz(2.4228754) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.047478288) q[2];
sx q[2];
rz(-0.55314174) q[2];
sx q[2];
rz(1.4609569) q[2];
rz(0.85584062) q[3];
sx q[3];
rz(-2.1688921) q[3];
sx q[3];
rz(2.262825) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(2.0395373) q[0];
sx q[0];
rz(-2.2886031) q[0];
sx q[0];
rz(-0.0041740388) q[0];
rz(-1.8904842) q[1];
sx q[1];
rz(-0.1336385) q[1];
sx q[1];
rz(2.4600696) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.174356) q[0];
sx q[0];
rz(-0.78034329) q[0];
sx q[0];
rz(1.1514949) q[0];
x q[1];
rz(2.2987492) q[2];
sx q[2];
rz(-2.3327391) q[2];
sx q[2];
rz(-3.0191772) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.2118069) q[1];
sx q[1];
rz(-2.8007467) q[1];
sx q[1];
rz(-0.58675672) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2596016) q[3];
sx q[3];
rz(-0.83899812) q[3];
sx q[3];
rz(1.9669718) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.556813) q[2];
sx q[2];
rz(-0.4759554) q[2];
sx q[2];
rz(-1.3947831) q[2];
rz(2.7252588) q[3];
sx q[3];
rz(-1.2952015) q[3];
sx q[3];
rz(-0.39286119) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.75547472) q[0];
sx q[0];
rz(-2.8135354) q[0];
sx q[0];
rz(1.0020483) q[0];
rz(0.26456061) q[1];
sx q[1];
rz(-1.8160276) q[1];
sx q[1];
rz(-1.6511668) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0984983) q[0];
sx q[0];
rz(-1.7014456) q[0];
sx q[0];
rz(0.54389531) q[0];
rz(-pi) q[1];
rz(2.3333346) q[2];
sx q[2];
rz(-1.6811451) q[2];
sx q[2];
rz(-2.7791952) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.7539053) q[1];
sx q[1];
rz(-1.4567944) q[1];
sx q[1];
rz(3.0579159) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0496554) q[3];
sx q[3];
rz(-1.6306393) q[3];
sx q[3];
rz(-0.53093225) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.1313021) q[2];
sx q[2];
rz(-1.4069858) q[2];
sx q[2];
rz(-1.4170125) q[2];
rz(2.8085282) q[3];
sx q[3];
rz(-0.76996961) q[3];
sx q[3];
rz(-1.3781579) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3604816) q[0];
sx q[0];
rz(-1.8753373) q[0];
sx q[0];
rz(-1.2486096) q[0];
rz(-0.026451182) q[1];
sx q[1];
rz(-1.5298264) q[1];
sx q[1];
rz(-1.5419921) q[1];
rz(1.1556861) q[2];
sx q[2];
rz(-2.2397556) q[2];
sx q[2];
rz(-0.02469183) q[2];
rz(-0.61458434) q[3];
sx q[3];
rz(-1.3086223) q[3];
sx q[3];
rz(2.3400459) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
