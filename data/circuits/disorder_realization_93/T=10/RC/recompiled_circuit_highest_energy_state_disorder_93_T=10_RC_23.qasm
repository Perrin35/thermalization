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
rz(-2.2095069) q[0];
sx q[0];
rz(-1.0340438) q[0];
sx q[0];
rz(-2.3873868) q[0];
rz(1.7475313) q[1];
sx q[1];
rz(-0.59352195) q[1];
sx q[1];
rz(-1.4374179) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0707929) q[0];
sx q[0];
rz(-1.5120929) q[0];
sx q[0];
rz(0.12639166) q[0];
rz(-pi) q[1];
rz(-0.038226323) q[2];
sx q[2];
rz(-1.9598205) q[2];
sx q[2];
rz(-1.807275) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.61753786) q[1];
sx q[1];
rz(-2.1270942) q[1];
sx q[1];
rz(-1.9281045) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2487605) q[3];
sx q[3];
rz(-0.40169558) q[3];
sx q[3];
rz(-0.63689771) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.6393911) q[2];
sx q[2];
rz(-2.907967) q[2];
sx q[2];
rz(-2.8761253) q[2];
rz(1.4207077) q[3];
sx q[3];
rz(-1.1666433) q[3];
sx q[3];
rz(-1.407912) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0625286) q[0];
sx q[0];
rz(-1.9363576) q[0];
sx q[0];
rz(-2.089654) q[0];
rz(-0.99986347) q[1];
sx q[1];
rz(-1.544416) q[1];
sx q[1];
rz(-2.3725407) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6175943) q[0];
sx q[0];
rz(-0.72223653) q[0];
sx q[0];
rz(-2.2666032) q[0];
rz(2.5639655) q[2];
sx q[2];
rz(-1.5211454) q[2];
sx q[2];
rz(0.28896633) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.96892541) q[1];
sx q[1];
rz(-1.7255867) q[1];
sx q[1];
rz(1.118578) q[1];
rz(-2.6399986) q[3];
sx q[3];
rz(-0.47712773) q[3];
sx q[3];
rz(-0.17233822) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.179788) q[2];
sx q[2];
rz(-0.804681) q[2];
sx q[2];
rz(-1.5847607) q[2];
rz(-0.70476091) q[3];
sx q[3];
rz(-0.2125936) q[3];
sx q[3];
rz(-2.052665) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8162808) q[0];
sx q[0];
rz(-0.75958696) q[0];
sx q[0];
rz(-2.3749206) q[0];
rz(-0.39241544) q[1];
sx q[1];
rz(-0.86154834) q[1];
sx q[1];
rz(-0.79889417) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1846599) q[0];
sx q[0];
rz(-2.6574832) q[0];
sx q[0];
rz(1.5501322) q[0];
rz(-pi) q[1];
rz(-1.4806108) q[2];
sx q[2];
rz(-1.8205368) q[2];
sx q[2];
rz(-1.1426999) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.4370712) q[1];
sx q[1];
rz(-2.066631) q[1];
sx q[1];
rz(-0.35192546) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9641482) q[3];
sx q[3];
rz(-1.465762) q[3];
sx q[3];
rz(-0.80897409) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.40156349) q[2];
sx q[2];
rz(-0.26569685) q[2];
sx q[2];
rz(-3.1166039) q[2];
rz(1.7449024) q[3];
sx q[3];
rz(-1.9091505) q[3];
sx q[3];
rz(-0.11740824) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5806737) q[0];
sx q[0];
rz(-0.54045254) q[0];
sx q[0];
rz(-0.29676357) q[0];
rz(-0.98776039) q[1];
sx q[1];
rz(-1.8901653) q[1];
sx q[1];
rz(2.3938649) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4982939) q[0];
sx q[0];
rz(-0.81658634) q[0];
sx q[0];
rz(0.30188456) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.50570391) q[2];
sx q[2];
rz(-0.78185588) q[2];
sx q[2];
rz(-3.0145558) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.26655134) q[1];
sx q[1];
rz(-1.5593441) q[1];
sx q[1];
rz(-2.7391866) q[1];
x q[2];
rz(-0.9763262) q[3];
sx q[3];
rz(-2.1552999) q[3];
sx q[3];
rz(0.85660958) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.1838386) q[2];
sx q[2];
rz(-0.87205333) q[2];
sx q[2];
rz(-0.5198861) q[2];
rz(-0.94684354) q[3];
sx q[3];
rz(-1.1917453) q[3];
sx q[3];
rz(-0.051518353) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2648322) q[0];
sx q[0];
rz(-2.2417289) q[0];
sx q[0];
rz(-1.9448036) q[0];
rz(-1.6668677) q[1];
sx q[1];
rz(-2.3366172) q[1];
sx q[1];
rz(0.29874745) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7525879) q[0];
sx q[0];
rz(-2.2225131) q[0];
sx q[0];
rz(1.1260004) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.79738209) q[2];
sx q[2];
rz(-1.935531) q[2];
sx q[2];
rz(1.9759751) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.1890002) q[1];
sx q[1];
rz(-0.77356427) q[1];
sx q[1];
rz(-0.1632593) q[1];
rz(-pi) q[2];
x q[2];
rz(0.87033724) q[3];
sx q[3];
rz(-1.8856109) q[3];
sx q[3];
rz(-0.39488068) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.8220736) q[2];
sx q[2];
rz(-2.3910429) q[2];
sx q[2];
rz(3.1263515) q[2];
rz(3.0139253) q[3];
sx q[3];
rz(-2.7805507) q[3];
sx q[3];
rz(0.84967363) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3540038) q[0];
sx q[0];
rz(-2.6241527) q[0];
sx q[0];
rz(-0.87376839) q[0];
rz(-1.2242225) q[1];
sx q[1];
rz(-2.1189549) q[1];
sx q[1];
rz(-1.2194941) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2827493) q[0];
sx q[0];
rz(-1.97012) q[0];
sx q[0];
rz(2.0612677) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5415433) q[2];
sx q[2];
rz(-0.82227899) q[2];
sx q[2];
rz(0.5136958) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.5883792) q[1];
sx q[1];
rz(-1.3822275) q[1];
sx q[1];
rz(1.7836573) q[1];
rz(-pi) q[2];
rz(2.5490646) q[3];
sx q[3];
rz(-1.9059407) q[3];
sx q[3];
rz(-0.0026897653) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.445861) q[2];
sx q[2];
rz(-1.4706688) q[2];
sx q[2];
rz(0.50103465) q[2];
rz(-1.3387574) q[3];
sx q[3];
rz(-0.41156358) q[3];
sx q[3];
rz(-1.1494273) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
x q[3];
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
rz(-2.5434791) q[0];
sx q[0];
rz(-1.1971594) q[0];
sx q[0];
rz(-0.57394779) q[0];
rz(2.5698938) q[1];
sx q[1];
rz(-2.1015002) q[1];
sx q[1];
rz(-1.5843102) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.27423647) q[0];
sx q[0];
rz(-2.353243) q[0];
sx q[0];
rz(2.1091745) q[0];
rz(-pi) q[1];
x q[1];
rz(0.36841285) q[2];
sx q[2];
rz(-0.94295687) q[2];
sx q[2];
rz(-1.6729015) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.3347335) q[1];
sx q[1];
rz(-0.99913952) q[1];
sx q[1];
rz(2.566659) q[1];
x q[2];
rz(-1.9628223) q[3];
sx q[3];
rz(-1.2268664) q[3];
sx q[3];
rz(1.0052538) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.456936) q[2];
sx q[2];
rz(-1.1494278) q[2];
sx q[2];
rz(1.1535025) q[2];
rz(1.614511) q[3];
sx q[3];
rz(-1.0228415) q[3];
sx q[3];
rz(1.3326741) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.033878) q[0];
sx q[0];
rz(-0.77924538) q[0];
sx q[0];
rz(-0.23767924) q[0];
rz(2.4029845) q[1];
sx q[1];
rz(-1.9269201) q[1];
sx q[1];
rz(-0.019066378) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4486074) q[0];
sx q[0];
rz(-3.0685022) q[0];
sx q[0];
rz(-1.4574217) q[0];
rz(2.3735061) q[2];
sx q[2];
rz(-1.0407018) q[2];
sx q[2];
rz(1.4370611) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.993666) q[1];
sx q[1];
rz(-1.4892409) q[1];
sx q[1];
rz(-3.102735) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4138149) q[3];
sx q[3];
rz(-0.78633868) q[3];
sx q[3];
rz(-2.0008519) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.0154524) q[2];
sx q[2];
rz(-1.3498053) q[2];
sx q[2];
rz(0.17507412) q[2];
rz(-1.763688) q[3];
sx q[3];
rz(-2.521069) q[3];
sx q[3];
rz(-3.0231045) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.058067583) q[0];
sx q[0];
rz(-1.3578992) q[0];
sx q[0];
rz(-0.51496664) q[0];
rz(-2.1452451) q[1];
sx q[1];
rz(-1.0641655) q[1];
sx q[1];
rz(-1.4774342) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.9290976) q[0];
sx q[0];
rz(-0.97508865) q[0];
sx q[0];
rz(0.86796772) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6040049) q[2];
sx q[2];
rz(-2.1589212) q[2];
sx q[2];
rz(0.7809446) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.2075726) q[1];
sx q[1];
rz(-2.0193384) q[1];
sx q[1];
rz(1.0058844) q[1];
rz(-pi) q[2];
rz(2.1636073) q[3];
sx q[3];
rz(-1.3372652) q[3];
sx q[3];
rz(-1.9625036) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.065206334) q[2];
sx q[2];
rz(-1.4896769) q[2];
sx q[2];
rz(-0.15929407) q[2];
rz(1.4984891) q[3];
sx q[3];
rz(-2.4703333) q[3];
sx q[3];
rz(1.8026277) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0945187) q[0];
sx q[0];
rz(-0.044059489) q[0];
sx q[0];
rz(2.2799168) q[0];
rz(-1.2791951) q[1];
sx q[1];
rz(-0.28293124) q[1];
sx q[1];
rz(-0.18712015) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8409654) q[0];
sx q[0];
rz(-1.8442796) q[0];
sx q[0];
rz(-1.3329559) q[0];
rz(-1.8420024) q[2];
sx q[2];
rz(-1.091963) q[2];
sx q[2];
rz(-1.2287272) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.71003733) q[1];
sx q[1];
rz(-1.8373243) q[1];
sx q[1];
rz(-2.4026067) q[1];
rz(1.8871103) q[3];
sx q[3];
rz(-1.1939002) q[3];
sx q[3];
rz(-2.047647) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.6992135) q[2];
sx q[2];
rz(-2.1670161) q[2];
sx q[2];
rz(-1.3628192) q[2];
rz(-1.5022701) q[3];
sx q[3];
rz(-2.5876741) q[3];
sx q[3];
rz(-1.7207918) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5300071) q[0];
sx q[0];
rz(-1.8066318) q[0];
sx q[0];
rz(-3.1351177) q[0];
rz(-0.50637983) q[1];
sx q[1];
rz(-2.9492999) q[1];
sx q[1];
rz(0.86023387) q[1];
rz(-2.996247) q[2];
sx q[2];
rz(-1.1517364) q[2];
sx q[2];
rz(0.20530063) q[2];
rz(-2.4615749) q[3];
sx q[3];
rz(-0.33878762) q[3];
sx q[3];
rz(1.9382077) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
