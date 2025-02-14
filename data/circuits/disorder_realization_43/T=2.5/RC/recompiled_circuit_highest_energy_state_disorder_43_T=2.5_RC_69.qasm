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
rz(0.046959538) q[0];
sx q[0];
rz(1.5927915) q[0];
sx q[0];
rz(11.172063) q[0];
rz(-0.46051639) q[1];
sx q[1];
rz(4.4098201) q[1];
sx q[1];
rz(8.9486651) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.61834419) q[0];
sx q[0];
rz(-1.1914413) q[0];
sx q[0];
rz(0.98486395) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3018487) q[2];
sx q[2];
rz(-1.7978872) q[2];
sx q[2];
rz(-0.64193945) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.3232267) q[1];
sx q[1];
rz(-1.4225679) q[1];
sx q[1];
rz(-1.5882951) q[1];
x q[2];
rz(-1.5067817) q[3];
sx q[3];
rz(-1.5435757) q[3];
sx q[3];
rz(0.79827362) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.8146879) q[2];
sx q[2];
rz(-1.9196332) q[2];
sx q[2];
rz(1.8587221) q[2];
rz(0.13947105) q[3];
sx q[3];
rz(-1.3160416) q[3];
sx q[3];
rz(-2.4563346) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
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
rz(-2.4271456) q[0];
sx q[0];
rz(-2.7834263) q[0];
sx q[0];
rz(-2.8156679) q[0];
rz(-0.70949316) q[1];
sx q[1];
rz(-0.26115099) q[1];
sx q[1];
rz(0.66169468) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.046592043) q[0];
sx q[0];
rz(-1.5612771) q[0];
sx q[0];
rz(3.1341888) q[0];
rz(-pi) q[1];
rz(2.2019655) q[2];
sx q[2];
rz(-1.9847615) q[2];
sx q[2];
rz(-2.8394073) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.0256309) q[1];
sx q[1];
rz(-0.44646663) q[1];
sx q[1];
rz(2.129351) q[1];
rz(-0.98425166) q[3];
sx q[3];
rz(-2.1681166) q[3];
sx q[3];
rz(2.190499) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.1208531) q[2];
sx q[2];
rz(-1.8031305) q[2];
sx q[2];
rz(-2.8815114) q[2];
rz(0.86380473) q[3];
sx q[3];
rz(-1.443202) q[3];
sx q[3];
rz(-1.9728194) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.23673713) q[0];
sx q[0];
rz(-0.78176347) q[0];
sx q[0];
rz(2.8535063) q[0];
rz(-1.9347363) q[1];
sx q[1];
rz(-2.0286045) q[1];
sx q[1];
rz(0.77817121) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.659063) q[0];
sx q[0];
rz(-0.89119833) q[0];
sx q[0];
rz(-2.8975455) q[0];
rz(-3.0765509) q[2];
sx q[2];
rz(-2.3823526) q[2];
sx q[2];
rz(-1.0763847) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.6699804) q[1];
sx q[1];
rz(-1.1410574) q[1];
sx q[1];
rz(2.5915036) q[1];
rz(-pi) q[2];
rz(-2.0708848) q[3];
sx q[3];
rz(-0.71083655) q[3];
sx q[3];
rz(2.5612166) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.8975767) q[2];
sx q[2];
rz(-2.0147169) q[2];
sx q[2];
rz(-0.24813949) q[2];
rz(1.0330307) q[3];
sx q[3];
rz(-1.4890198) q[3];
sx q[3];
rz(1.2584794) q[3];
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
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7018062) q[0];
sx q[0];
rz(-0.93455625) q[0];
sx q[0];
rz(0.55950657) q[0];
rz(1.7350381) q[1];
sx q[1];
rz(-2.1969257) q[1];
sx q[1];
rz(-1.14538) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9713719) q[0];
sx q[0];
rz(-1.6165501) q[0];
sx q[0];
rz(2.6849823) q[0];
rz(-pi) q[1];
x q[1];
rz(0.68762806) q[2];
sx q[2];
rz(-2.3634644) q[2];
sx q[2];
rz(-1.8541921) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.0068906764) q[1];
sx q[1];
rz(-0.81738735) q[1];
sx q[1];
rz(2.7322053) q[1];
rz(0.01500742) q[3];
sx q[3];
rz(-2.1405485) q[3];
sx q[3];
rz(0.80591312) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.8765325) q[2];
sx q[2];
rz(-2.127485) q[2];
sx q[2];
rz(1.8939135) q[2];
rz(-1.0809336) q[3];
sx q[3];
rz(-1.7526046) q[3];
sx q[3];
rz(-1.8814253) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0697698) q[0];
sx q[0];
rz(-0.48957303) q[0];
sx q[0];
rz(-2.3967632) q[0];
rz(1.3818332) q[1];
sx q[1];
rz(-1.1188743) q[1];
sx q[1];
rz(0.098085731) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25135219) q[0];
sx q[0];
rz(-0.47516631) q[0];
sx q[0];
rz(0.71936468) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6977243) q[2];
sx q[2];
rz(-1.265268) q[2];
sx q[2];
rz(1.533184) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.262935) q[1];
sx q[1];
rz(-2.2644357) q[1];
sx q[1];
rz(-1.5946112) q[1];
rz(-pi) q[2];
rz(2.0518584) q[3];
sx q[3];
rz(-1.6599382) q[3];
sx q[3];
rz(-1.3463904) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.44094008) q[2];
sx q[2];
rz(-2.7589189) q[2];
sx q[2];
rz(1.0901701) q[2];
rz(0.32554659) q[3];
sx q[3];
rz(-1.5109589) q[3];
sx q[3];
rz(-0.24698273) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.22661041) q[0];
sx q[0];
rz(-2.60422) q[0];
sx q[0];
rz(1.2363303) q[0];
rz(-0.63713282) q[1];
sx q[1];
rz(-0.56012154) q[1];
sx q[1];
rz(0.93322745) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7819311) q[0];
sx q[0];
rz(-1.0190411) q[0];
sx q[0];
rz(-0.42858584) q[0];
rz(-pi) q[1];
x q[1];
rz(0.18913375) q[2];
sx q[2];
rz(-1.1580843) q[2];
sx q[2];
rz(-2.7688062) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.56343397) q[1];
sx q[1];
rz(-1.2799885) q[1];
sx q[1];
rz(-1.6573424) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0108443) q[3];
sx q[3];
rz(-2.5820316) q[3];
sx q[3];
rz(2.6236629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.25202641) q[2];
sx q[2];
rz(-1.7539975) q[2];
sx q[2];
rz(1.265556) q[2];
rz(1.7172074) q[3];
sx q[3];
rz(-0.5286743) q[3];
sx q[3];
rz(2.2831447) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3739361) q[0];
sx q[0];
rz(-0.38145426) q[0];
sx q[0];
rz(-0.23442991) q[0];
rz(2.7081721) q[1];
sx q[1];
rz(-2.0442918) q[1];
sx q[1];
rz(2.0030599) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.075468242) q[0];
sx q[0];
rz(-0.99695092) q[0];
sx q[0];
rz(-2.2908151) q[0];
rz(-pi) q[1];
rz(-0.01485558) q[2];
sx q[2];
rz(-0.49070036) q[2];
sx q[2];
rz(-2.4717836) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.1416152) q[1];
sx q[1];
rz(-2.6385698) q[1];
sx q[1];
rz(-1.5161689) q[1];
rz(-pi) q[2];
rz(-2.906231) q[3];
sx q[3];
rz(-2.4409264) q[3];
sx q[3];
rz(1.2149902) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.66702691) q[2];
sx q[2];
rz(-2.462025) q[2];
sx q[2];
rz(2.1843074) q[2];
rz(2.3840617) q[3];
sx q[3];
rz(-1.4742955) q[3];
sx q[3];
rz(-0.1434513) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2738508) q[0];
sx q[0];
rz(-1.1026646) q[0];
sx q[0];
rz(-2.962501) q[0];
rz(0.20009072) q[1];
sx q[1];
rz(-0.6691907) q[1];
sx q[1];
rz(2.3650513) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6497407) q[0];
sx q[0];
rz(-2.092768) q[0];
sx q[0];
rz(-2.9328547) q[0];
x q[1];
rz(2.8440823) q[2];
sx q[2];
rz(-2.1665467) q[2];
sx q[2];
rz(-1.7864625) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.84785968) q[1];
sx q[1];
rz(-1.9256448) q[1];
sx q[1];
rz(-1.4925641) q[1];
rz(-pi) q[2];
rz(3.0902571) q[3];
sx q[3];
rz(-1.8881137) q[3];
sx q[3];
rz(2.7961344) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.17579235) q[2];
sx q[2];
rz(-1.0719904) q[2];
sx q[2];
rz(1.779186) q[2];
rz(1.8386748) q[3];
sx q[3];
rz(-1.4785654) q[3];
sx q[3];
rz(2.1030857) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68504828) q[0];
sx q[0];
rz(-1.9467204) q[0];
sx q[0];
rz(-3.0273279) q[0];
rz(-1.9396797) q[1];
sx q[1];
rz(-2.6004531) q[1];
sx q[1];
rz(-3.0466383) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1916434) q[0];
sx q[0];
rz(-1.7472634) q[0];
sx q[0];
rz(2.1032189) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0509148) q[2];
sx q[2];
rz(-0.78453817) q[2];
sx q[2];
rz(2.926411) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.11647955) q[1];
sx q[1];
rz(-1.0214865) q[1];
sx q[1];
rz(1.504937) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0614653) q[3];
sx q[3];
rz(-1.2042787) q[3];
sx q[3];
rz(1.7297922) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.8408884) q[2];
sx q[2];
rz(-1.3450832) q[2];
sx q[2];
rz(-2.4697206) q[2];
rz(-2.4554409) q[3];
sx q[3];
rz(-0.66303623) q[3];
sx q[3];
rz(-1.9239976) q[3];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1662064) q[0];
sx q[0];
rz(-2.9500742) q[0];
sx q[0];
rz(1.6084877) q[0];
rz(0.32471049) q[1];
sx q[1];
rz(-1.277781) q[1];
sx q[1];
rz(0.3130354) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.19713773) q[0];
sx q[0];
rz(-0.99992434) q[0];
sx q[0];
rz(-0.095143347) q[0];
rz(-1.8386449) q[2];
sx q[2];
rz(-0.22780475) q[2];
sx q[2];
rz(-1.1033664) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.6021784) q[1];
sx q[1];
rz(-1.0262607) q[1];
sx q[1];
rz(2.3474752) q[1];
rz(-0.62449091) q[3];
sx q[3];
rz(-1.1905128) q[3];
sx q[3];
rz(2.6759669) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.2031871) q[2];
sx q[2];
rz(-1.9878191) q[2];
sx q[2];
rz(-1.3333092) q[2];
rz(-0.61776727) q[3];
sx q[3];
rz(-0.94527644) q[3];
sx q[3];
rz(-1.1355404) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7137322) q[0];
sx q[0];
rz(-1.8129616) q[0];
sx q[0];
rz(0.34995361) q[0];
rz(2.2433544) q[1];
sx q[1];
rz(-1.0844834) q[1];
sx q[1];
rz(-0.35379298) q[1];
rz(-2.3980065) q[2];
sx q[2];
rz(-1.0844885) q[2];
sx q[2];
rz(1.1906638) q[2];
rz(1.5469503) q[3];
sx q[3];
rz(-2.3004354) q[3];
sx q[3];
rz(-2.7247747) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
