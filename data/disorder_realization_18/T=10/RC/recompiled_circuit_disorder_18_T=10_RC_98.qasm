OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.3192531) q[0];
sx q[0];
rz(-2.1847794) q[0];
sx q[0];
rz(1.7083038) q[0];
rz(-0.23694555) q[1];
sx q[1];
rz(-1.911093) q[1];
sx q[1];
rz(1.9967611) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5231397) q[0];
sx q[0];
rz(-1.3292399) q[0];
sx q[0];
rz(3.0735344) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.98923367) q[2];
sx q[2];
rz(-0.49058149) q[2];
sx q[2];
rz(-1.4197592) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.45935985) q[1];
sx q[1];
rz(-1.1182922) q[1];
sx q[1];
rz(-0.020999055) q[1];
rz(-pi) q[2];
rz(-2.6712816) q[3];
sx q[3];
rz(-0.59578005) q[3];
sx q[3];
rz(-2.6252928) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(3.0248489) q[2];
sx q[2];
rz(-1.3329093) q[2];
sx q[2];
rz(-1.9072745) q[2];
rz(1.0553137) q[3];
sx q[3];
rz(-1.0588667) q[3];
sx q[3];
rz(-1.9887259) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.18594436) q[0];
sx q[0];
rz(-1.1971104) q[0];
sx q[0];
rz(-2.3117075) q[0];
rz(-2.8886967) q[1];
sx q[1];
rz(-2.3650832) q[1];
sx q[1];
rz(0.84709644) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9053099) q[0];
sx q[0];
rz(-0.61203996) q[0];
sx q[0];
rz(2.406714) q[0];
x q[1];
rz(-0.77913021) q[2];
sx q[2];
rz(-1.3574294) q[2];
sx q[2];
rz(1.9386292) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.9287195) q[1];
sx q[1];
rz(-1.6391014) q[1];
sx q[1];
rz(-0.54547711) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2291525) q[3];
sx q[3];
rz(-0.53386253) q[3];
sx q[3];
rz(-2.562059) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.0236686) q[2];
sx q[2];
rz(-0.37332049) q[2];
sx q[2];
rz(-0.71933293) q[2];
rz(-1.2891278) q[3];
sx q[3];
rz(-0.23580655) q[3];
sx q[3];
rz(1.4891362) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
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
rz(-1.2704724) q[0];
sx q[0];
rz(-1.8914762) q[0];
sx q[0];
rz(2.5701994) q[0];
rz(2.1748623) q[1];
sx q[1];
rz(-1.2582018) q[1];
sx q[1];
rz(0.5273231) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8427211) q[0];
sx q[0];
rz(-1.5639595) q[0];
sx q[0];
rz(-0.30928916) q[0];
rz(-pi) q[1];
x q[1];
rz(0.6730404) q[2];
sx q[2];
rz(-1.3125784) q[2];
sx q[2];
rz(-2.8845127) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.61422435) q[1];
sx q[1];
rz(-0.80263153) q[1];
sx q[1];
rz(-0.36348344) q[1];
x q[2];
rz(1.2445883) q[3];
sx q[3];
rz(-2.2766114) q[3];
sx q[3];
rz(2.0446442) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.8335235) q[2];
sx q[2];
rz(-1.4493194) q[2];
sx q[2];
rz(2.4148338) q[2];
rz(0.99749342) q[3];
sx q[3];
rz(-0.62524978) q[3];
sx q[3];
rz(1.5184901) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8961287) q[0];
sx q[0];
rz(-1.4828869) q[0];
sx q[0];
rz(1.0035275) q[0];
rz(3.1009122) q[1];
sx q[1];
rz(-1.129312) q[1];
sx q[1];
rz(-0.8262659) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4429312) q[0];
sx q[0];
rz(-0.94416617) q[0];
sx q[0];
rz(2.4311964) q[0];
rz(-pi) q[1];
x q[1];
rz(1.516953) q[2];
sx q[2];
rz(-0.87140897) q[2];
sx q[2];
rz(-2.2211423) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.0127718) q[1];
sx q[1];
rz(-1.9073309) q[1];
sx q[1];
rz(-3.0613042) q[1];
rz(-2.6310001) q[3];
sx q[3];
rz(-0.69805745) q[3];
sx q[3];
rz(-1.5654246) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.56759175) q[2];
sx q[2];
rz(-1.2446128) q[2];
sx q[2];
rz(-0.91439247) q[2];
rz(-3.0363723) q[3];
sx q[3];
rz(-1.6243694) q[3];
sx q[3];
rz(-2.5523394) q[3];
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
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2545664) q[0];
sx q[0];
rz(-1.6396602) q[0];
sx q[0];
rz(1.1337093) q[0];
rz(3.1359613) q[1];
sx q[1];
rz(-1.7040323) q[1];
sx q[1];
rz(2.5240135) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3591946) q[0];
sx q[0];
rz(-3.087128) q[0];
sx q[0];
rz(-0.99476238) q[0];
rz(2.4992906) q[2];
sx q[2];
rz(-1.2080492) q[2];
sx q[2];
rz(1.158266) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.51296556) q[1];
sx q[1];
rz(-1.8497397) q[1];
sx q[1];
rz(-0.97475027) q[1];
rz(-pi) q[2];
rz(0.17554749) q[3];
sx q[3];
rz(-1.5581589) q[3];
sx q[3];
rz(-2.4059084) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.3474943) q[2];
sx q[2];
rz(-1.2604875) q[2];
sx q[2];
rz(-0.23362544) q[2];
rz(2.1485093) q[3];
sx q[3];
rz(-2.1916094) q[3];
sx q[3];
rz(-2.5036507) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.96930209) q[0];
sx q[0];
rz(-0.06047051) q[0];
sx q[0];
rz(-0.0897952) q[0];
rz(-0.88109294) q[1];
sx q[1];
rz(-2.6125364) q[1];
sx q[1];
rz(0.18009137) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.8343617) q[0];
sx q[0];
rz(-2.0198856) q[0];
sx q[0];
rz(2.2063072) q[0];
rz(-pi) q[1];
rz(-1.8413576) q[2];
sx q[2];
rz(-2.5205043) q[2];
sx q[2];
rz(1.6294711) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.084701531) q[1];
sx q[1];
rz(-1.7329721) q[1];
sx q[1];
rz(0.2918891) q[1];
rz(-3.1230934) q[3];
sx q[3];
rz(-1.0780932) q[3];
sx q[3];
rz(2.8871356) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.883541) q[2];
sx q[2];
rz(-1.5774612) q[2];
sx q[2];
rz(-1.1759261) q[2];
rz(0.073143395) q[3];
sx q[3];
rz(-0.72435838) q[3];
sx q[3];
rz(2.6426962) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.183855) q[0];
sx q[0];
rz(-2.4860005) q[0];
sx q[0];
rz(1.7234329) q[0];
rz(-1.1649959) q[1];
sx q[1];
rz(-0.55924758) q[1];
sx q[1];
rz(-0.040239008) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7559164) q[0];
sx q[0];
rz(-1.3323116) q[0];
sx q[0];
rz(3.0576865) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4637429) q[2];
sx q[2];
rz(-0.40823001) q[2];
sx q[2];
rz(1.6183311) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.7525967) q[1];
sx q[1];
rz(-2.1850359) q[1];
sx q[1];
rz(2.9700301) q[1];
rz(-1.8764898) q[3];
sx q[3];
rz(-1.8743268) q[3];
sx q[3];
rz(-2.2895253) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.13850257) q[2];
sx q[2];
rz(-2.3562727) q[2];
sx q[2];
rz(0.17383943) q[2];
rz(-1.396817) q[3];
sx q[3];
rz(-1.6766179) q[3];
sx q[3];
rz(2.2824536) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.023507) q[0];
sx q[0];
rz(-0.21723391) q[0];
sx q[0];
rz(-1.7370976) q[0];
rz(-1.2069758) q[1];
sx q[1];
rz(-1.4852306) q[1];
sx q[1];
rz(1.6361902) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0597499) q[0];
sx q[0];
rz(-1.2568226) q[0];
sx q[0];
rz(-0.44072515) q[0];
rz(-pi) q[1];
rz(2.5152399) q[2];
sx q[2];
rz(-0.92712958) q[2];
sx q[2];
rz(2.6477637) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.9951524) q[1];
sx q[1];
rz(-0.69680981) q[1];
sx q[1];
rz(3.0909096) q[1];
x q[2];
rz(-2.1721341) q[3];
sx q[3];
rz(-1.2974032) q[3];
sx q[3];
rz(-2.6858342) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.3796842) q[2];
sx q[2];
rz(-0.48019335) q[2];
sx q[2];
rz(-2.9910679) q[2];
rz(-1.5395509) q[3];
sx q[3];
rz(-1.1952885) q[3];
sx q[3];
rz(0.62301821) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4432916) q[0];
sx q[0];
rz(-1.1346096) q[0];
sx q[0];
rz(0.64569965) q[0];
rz(-2.6422016) q[1];
sx q[1];
rz(-1.7130518) q[1];
sx q[1];
rz(-1.2649149) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9547573) q[0];
sx q[0];
rz(-1.1823913) q[0];
sx q[0];
rz(0.79083058) q[0];
rz(-pi) q[1];
x q[1];
rz(1.132498) q[2];
sx q[2];
rz(-2.0813137) q[2];
sx q[2];
rz(-2.8709708) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.1688655) q[1];
sx q[1];
rz(-1.7164413) q[1];
sx q[1];
rz(-1.5345854) q[1];
x q[2];
rz(-1.7725138) q[3];
sx q[3];
rz(-1.9864051) q[3];
sx q[3];
rz(-2.6430074) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.50679961) q[2];
sx q[2];
rz(-0.96968499) q[2];
sx q[2];
rz(2.6169422) q[2];
rz(-2.5850885) q[3];
sx q[3];
rz(-1.6289214) q[3];
sx q[3];
rz(2.7867253) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5234579) q[0];
sx q[0];
rz(-2.3045461) q[0];
sx q[0];
rz(0.24542228) q[0];
rz(2.5852809) q[1];
sx q[1];
rz(-1.6106771) q[1];
sx q[1];
rz(1.6773178) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.9058162) q[0];
sx q[0];
rz(-0.24484867) q[0];
sx q[0];
rz(-0.4723627) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9015058) q[2];
sx q[2];
rz(-0.61243528) q[2];
sx q[2];
rz(-1.9233179) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.6357395) q[1];
sx q[1];
rz(-0.74843279) q[1];
sx q[1];
rz(-1.9745419) q[1];
rz(-pi) q[2];
rz(3.0562923) q[3];
sx q[3];
rz(-1.395426) q[3];
sx q[3];
rz(-2.8773957) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.8993373) q[2];
sx q[2];
rz(-1.9031886) q[2];
sx q[2];
rz(-2.396092) q[2];
rz(1.7177104) q[3];
sx q[3];
rz(-2.8427377) q[3];
sx q[3];
rz(2.0132813) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2765008) q[0];
sx q[0];
rz(-1.6456974) q[0];
sx q[0];
rz(0.67281848) q[0];
rz(1.2561692) q[1];
sx q[1];
rz(-0.80614631) q[1];
sx q[1];
rz(2.0731906) q[1];
rz(2.4297498) q[2];
sx q[2];
rz(-1.1887475) q[2];
sx q[2];
rz(2.9954994) q[2];
rz(0.54995723) q[3];
sx q[3];
rz(-2.4933542) q[3];
sx q[3];
rz(-1.8174432) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
