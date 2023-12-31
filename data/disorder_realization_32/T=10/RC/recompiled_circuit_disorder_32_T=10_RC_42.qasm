OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.10652868) q[0];
sx q[0];
rz(-1.0892692) q[0];
sx q[0];
rz(-2.9805592) q[0];
rz(1.610202) q[1];
sx q[1];
rz(2.664497) q[1];
sx q[1];
rz(8.9283979) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.30713233) q[0];
sx q[0];
rz(-1.4559329) q[0];
sx q[0];
rz(-2.1644724) q[0];
x q[1];
rz(0.52486323) q[2];
sx q[2];
rz(-0.3877936) q[2];
sx q[2];
rz(0.85688574) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(3.1094184) q[1];
sx q[1];
rz(-0.86583455) q[1];
sx q[1];
rz(-1.6900307) q[1];
rz(-1.6707889) q[3];
sx q[3];
rz(-1.9766207) q[3];
sx q[3];
rz(2.0506746) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.51497841) q[2];
sx q[2];
rz(-0.86003059) q[2];
sx q[2];
rz(0.28764763) q[2];
rz(-1.7488165) q[3];
sx q[3];
rz(-2.1578622) q[3];
sx q[3];
rz(1.1028517) q[3];
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
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2709133) q[0];
sx q[0];
rz(-2.3864855) q[0];
sx q[0];
rz(2.5640008) q[0];
rz(-1.6060991) q[1];
sx q[1];
rz(-1.1178733) q[1];
sx q[1];
rz(-1.577852) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0751511) q[0];
sx q[0];
rz(-1.6834714) q[0];
sx q[0];
rz(0.22002797) q[0];
rz(-pi) q[1];
rz(0.73194506) q[2];
sx q[2];
rz(-0.39790301) q[2];
sx q[2];
rz(2.0741472) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.3040846) q[1];
sx q[1];
rz(-1.6581074) q[1];
sx q[1];
rz(-2.0259894) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3489477) q[3];
sx q[3];
rz(-2.9229197) q[3];
sx q[3];
rz(-1.6817026) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-3.0343798) q[2];
sx q[2];
rz(-0.97637525) q[2];
sx q[2];
rz(-1.0822901) q[2];
rz(1.9783431) q[3];
sx q[3];
rz(-1.9415559) q[3];
sx q[3];
rz(2.5015586) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0796233) q[0];
sx q[0];
rz(-2.847023) q[0];
sx q[0];
rz(0.18297718) q[0];
rz(0.043116365) q[1];
sx q[1];
rz(-0.95298195) q[1];
sx q[1];
rz(1.9972237) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2079175) q[0];
sx q[0];
rz(-1.1753923) q[0];
sx q[0];
rz(-2.140722) q[0];
rz(-pi) q[1];
rz(-0.47282131) q[2];
sx q[2];
rz(-2.678215) q[2];
sx q[2];
rz(2.9544427) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.80013393) q[1];
sx q[1];
rz(-1.1143436) q[1];
sx q[1];
rz(-1.9564499) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8539682) q[3];
sx q[3];
rz(-1.5737185) q[3];
sx q[3];
rz(-0.50141108) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.1220876) q[2];
sx q[2];
rz(-1.959789) q[2];
sx q[2];
rz(-2.2369177) q[2];
rz(2.4217862) q[3];
sx q[3];
rz(-0.42680877) q[3];
sx q[3];
rz(0.75511801) q[3];
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
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6782137) q[0];
sx q[0];
rz(-2.1713874) q[0];
sx q[0];
rz(-2.4257207) q[0];
rz(2.8158358) q[1];
sx q[1];
rz(-2.082086) q[1];
sx q[1];
rz(0.99266565) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2199812) q[0];
sx q[0];
rz(-2.5489759) q[0];
sx q[0];
rz(2.6556334) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2797212) q[2];
sx q[2];
rz(-2.3700691) q[2];
sx q[2];
rz(0.18415235) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.0419473) q[1];
sx q[1];
rz(-1.4775659) q[1];
sx q[1];
rz(-2.9779742) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.643814) q[3];
sx q[3];
rz(-2.0911502) q[3];
sx q[3];
rz(1.2342412) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.0585534) q[2];
sx q[2];
rz(-0.71648592) q[2];
sx q[2];
rz(1.7129718) q[2];
rz(2.2955017) q[3];
sx q[3];
rz(-1.6276136) q[3];
sx q[3];
rz(-1.9184453) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5457526) q[0];
sx q[0];
rz(-1.7359474) q[0];
sx q[0];
rz(-2.3748421) q[0];
rz(1.9013566) q[1];
sx q[1];
rz(-2.2940472) q[1];
sx q[1];
rz(0.3516745) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36067097) q[0];
sx q[0];
rz(-0.70510266) q[0];
sx q[0];
rz(2.6893925) q[0];
rz(-pi) q[1];
x q[1];
rz(0.50260966) q[2];
sx q[2];
rz(-0.99684274) q[2];
sx q[2];
rz(0.9290907) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.000948) q[1];
sx q[1];
rz(-1.5854156) q[1];
sx q[1];
rz(-0.39477243) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1122123) q[3];
sx q[3];
rz(-0.57313985) q[3];
sx q[3];
rz(1.0669607) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.9436283) q[2];
sx q[2];
rz(-1.4732271) q[2];
sx q[2];
rz(-2.6082805) q[2];
rz(0.42896459) q[3];
sx q[3];
rz(-0.52754378) q[3];
sx q[3];
rz(-1.126948) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11739843) q[0];
sx q[0];
rz(-1.0929996) q[0];
sx q[0];
rz(-2.5999516) q[0];
rz(2.533124) q[1];
sx q[1];
rz(-0.21921961) q[1];
sx q[1];
rz(-1.0995964) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8655411) q[0];
sx q[0];
rz(-1.2035032) q[0];
sx q[0];
rz(2.0588576) q[0];
x q[1];
rz(-2.1018283) q[2];
sx q[2];
rz(-1.4073939) q[2];
sx q[2];
rz(0.49887564) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.9601599) q[1];
sx q[1];
rz(-0.90248855) q[1];
sx q[1];
rz(-1.798435) q[1];
rz(-pi) q[2];
rz(-1.2332637) q[3];
sx q[3];
rz(-0.81157717) q[3];
sx q[3];
rz(-2.8525762) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.6301443) q[2];
sx q[2];
rz(-1.6161329) q[2];
sx q[2];
rz(-0.28298322) q[2];
rz(-2.4354368) q[3];
sx q[3];
rz(-1.5420087) q[3];
sx q[3];
rz(0.63505665) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(-3.0866078) q[0];
sx q[0];
rz(-2.1540756) q[0];
sx q[0];
rz(-1.6116066) q[0];
rz(-2.4967172) q[1];
sx q[1];
rz(-0.49352831) q[1];
sx q[1];
rz(-2.9842916) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.776942) q[0];
sx q[0];
rz(-1.4770664) q[0];
sx q[0];
rz(-2.4073497) q[0];
rz(-1.332875) q[2];
sx q[2];
rz(-0.750713) q[2];
sx q[2];
rz(-2.6964292) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.7423319) q[1];
sx q[1];
rz(-0.6472339) q[1];
sx q[1];
rz(-1.2433692) q[1];
x q[2];
rz(1.8623452) q[3];
sx q[3];
rz(-0.96372094) q[3];
sx q[3];
rz(2.8214422) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.9592231) q[2];
sx q[2];
rz(-2.5490641) q[2];
sx q[2];
rz(2.3106993) q[2];
rz(0.043878555) q[3];
sx q[3];
rz(-1.9079804) q[3];
sx q[3];
rz(2.9390745) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6458994) q[0];
sx q[0];
rz(-2.2024787) q[0];
sx q[0];
rz(-3.1307401) q[0];
rz(-0.27663484) q[1];
sx q[1];
rz(-2.3697772) q[1];
sx q[1];
rz(2.8483134) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5921708) q[0];
sx q[0];
rz(-2.5991419) q[0];
sx q[0];
rz(-0.011499238) q[0];
rz(-pi) q[1];
rz(-2.9759334) q[2];
sx q[2];
rz(-0.93369166) q[2];
sx q[2];
rz(0.13656244) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.7640904) q[1];
sx q[1];
rz(-1.1524156) q[1];
sx q[1];
rz(-0.80470316) q[1];
rz(-pi) q[2];
rz(-1.9803489) q[3];
sx q[3];
rz(-2.000994) q[3];
sx q[3];
rz(2.6524515) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(3.098112) q[2];
sx q[2];
rz(-0.094823368) q[2];
sx q[2];
rz(-3.0528255) q[2];
rz(-0.25012112) q[3];
sx q[3];
rz(-1.1238031) q[3];
sx q[3];
rz(0.5789825) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0042689) q[0];
sx q[0];
rz(-2.9172638) q[0];
sx q[0];
rz(1.0639169) q[0];
rz(2.6990199) q[1];
sx q[1];
rz(-1.4241709) q[1];
sx q[1];
rz(-1.7907422) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3063072) q[0];
sx q[0];
rz(-1.6500236) q[0];
sx q[0];
rz(1.6552734) q[0];
rz(-pi) q[1];
rz(-1.4883792) q[2];
sx q[2];
rz(-0.91225183) q[2];
sx q[2];
rz(1.8694307) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.0693448) q[1];
sx q[1];
rz(-0.43466553) q[1];
sx q[1];
rz(-1.5321295) q[1];
x q[2];
rz(-1.7926932) q[3];
sx q[3];
rz(-1.4699234) q[3];
sx q[3];
rz(1.6290806) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.110048) q[2];
sx q[2];
rz(-1.2908547) q[2];
sx q[2];
rz(-2.6943977) q[2];
rz(1.3859008) q[3];
sx q[3];
rz(-0.30062506) q[3];
sx q[3];
rz(0.51469222) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31496012) q[0];
sx q[0];
rz(-1.5216014) q[0];
sx q[0];
rz(0.78053027) q[0];
rz(2.7138117) q[1];
sx q[1];
rz(-0.3573187) q[1];
sx q[1];
rz(-2.9719877) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1039625) q[0];
sx q[0];
rz(-2.5314405) q[0];
sx q[0];
rz(-0.9141586) q[0];
rz(-pi) q[1];
rz(-2.7197003) q[2];
sx q[2];
rz(-1.5932398) q[2];
sx q[2];
rz(0.23082146) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.6360213) q[1];
sx q[1];
rz(-2.1850056) q[1];
sx q[1];
rz(-3.0843656) q[1];
x q[2];
rz(-1.2492368) q[3];
sx q[3];
rz(-1.7064629) q[3];
sx q[3];
rz(-0.78007573) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.85049373) q[2];
sx q[2];
rz(-1.9591745) q[2];
sx q[2];
rz(-0.69236857) q[2];
rz(2.678357) q[3];
sx q[3];
rz(-1.4866456) q[3];
sx q[3];
rz(1.108981) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2904084) q[0];
sx q[0];
rz(-1.52607) q[0];
sx q[0];
rz(-2.4551256) q[0];
rz(-0.51207536) q[1];
sx q[1];
rz(-2.8072186) q[1];
sx q[1];
rz(1.0288815) q[1];
rz(0.26906536) q[2];
sx q[2];
rz(-1.8620327) q[2];
sx q[2];
rz(-3.0897683) q[2];
rz(0.11733304) q[3];
sx q[3];
rz(-1.3640079) q[3];
sx q[3];
rz(-0.97466536) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
