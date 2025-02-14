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
rz(-1.363938) q[0];
sx q[0];
rz(3.5600297) q[0];
sx q[0];
rz(10.726396) q[0];
rz(-2.9786181) q[1];
sx q[1];
rz(-1.6956704) q[1];
sx q[1];
rz(-0.016782848) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.421464) q[0];
sx q[0];
rz(-1.6454433) q[0];
sx q[0];
rz(1.9326841) q[0];
rz(-0.016525908) q[2];
sx q[2];
rz(-1.3912829) q[2];
sx q[2];
rz(-1.5915888) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.3455194) q[1];
sx q[1];
rz(-1.5658448) q[1];
sx q[1];
rz(-1.5758745) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9892606) q[3];
sx q[3];
rz(-1.6248584) q[3];
sx q[3];
rz(-0.033933725) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.5825384) q[2];
sx q[2];
rz(-2.4808919) q[2];
sx q[2];
rz(-1.5622697) q[2];
rz(-2.9301379) q[3];
sx q[3];
rz(-3.1410757) q[3];
sx q[3];
rz(2.9805984) q[3];
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
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6120537) q[0];
sx q[0];
rz(-0.27612975) q[0];
sx q[0];
rz(-1.8147234) q[0];
rz(0.57948411) q[1];
sx q[1];
rz(-3.1377628) q[1];
sx q[1];
rz(-0.63900596) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9352183) q[0];
sx q[0];
rz(-0.91418302) q[0];
sx q[0];
rz(-2.4103863) q[0];
x q[1];
rz(-0.1370378) q[2];
sx q[2];
rz(-0.1220905) q[2];
sx q[2];
rz(-0.11954319) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.745805) q[1];
sx q[1];
rz(-1.5823872) q[1];
sx q[1];
rz(-3.1245924) q[1];
rz(-2.3611446) q[3];
sx q[3];
rz(-1.5089581) q[3];
sx q[3];
rz(-2.0831747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.4352033) q[2];
sx q[2];
rz(-3.0053164) q[2];
sx q[2];
rz(-1.5380247) q[2];
rz(-1.5609353) q[3];
sx q[3];
rz(-3.1272562) q[3];
sx q[3];
rz(-3.1106136) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2494217) q[0];
sx q[0];
rz(-0.51512655) q[0];
sx q[0];
rz(0.38145915) q[0];
rz(2.4341266) q[1];
sx q[1];
rz(-0.019376945) q[1];
sx q[1];
rz(1.1245419) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8302001) q[0];
sx q[0];
rz(-1.8233577) q[0];
sx q[0];
rz(2.8930032) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4558305) q[2];
sx q[2];
rz(-1.5449776) q[2];
sx q[2];
rz(1.7242277) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.687011) q[1];
sx q[1];
rz(-3.0774766) q[1];
sx q[1];
rz(-2.9518576) q[1];
rz(-pi) q[2];
rz(-0.55508681) q[3];
sx q[3];
rz(-2.0197809) q[3];
sx q[3];
rz(2.4628696) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.4418929) q[2];
sx q[2];
rz(-0.012233891) q[2];
sx q[2];
rz(-3.0921248) q[2];
rz(-2.5345645) q[3];
sx q[3];
rz(-3.1403465) q[3];
sx q[3];
rz(-1.2073257) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1599051) q[0];
sx q[0];
rz(-2.9732381) q[0];
sx q[0];
rz(0.017177563) q[0];
rz(0.29302868) q[1];
sx q[1];
rz(-0.79048645) q[1];
sx q[1];
rz(-1.5944098) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3395383) q[0];
sx q[0];
rz(-1.0419163) q[0];
sx q[0];
rz(0.89978973) q[0];
x q[1];
rz(0.36357673) q[2];
sx q[2];
rz(-0.57206094) q[2];
sx q[2];
rz(0.043302082) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.35369998) q[1];
sx q[1];
rz(-0.13772923) q[1];
sx q[1];
rz(2.6929123) q[1];
x q[2];
rz(-3.1286448) q[3];
sx q[3];
rz(-2.4955813) q[3];
sx q[3];
rz(2.4821441) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.8881417) q[2];
sx q[2];
rz(-0.47051045) q[2];
sx q[2];
rz(0.68224254) q[2];
rz(3.0864129) q[3];
sx q[3];
rz(-0.0076871593) q[3];
sx q[3];
rz(1.3050219) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3799915) q[0];
sx q[0];
rz(-0.43552265) q[0];
sx q[0];
rz(0.4250266) q[0];
rz(-1.60166) q[1];
sx q[1];
rz(-2.6585177) q[1];
sx q[1];
rz(2.3262598) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3638301) q[0];
sx q[0];
rz(-1.6319931) q[0];
sx q[0];
rz(1.5817002) q[0];
x q[1];
rz(-1.4311976) q[2];
sx q[2];
rz(-3.1181212) q[2];
sx q[2];
rz(-1.0271629) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.968135) q[1];
sx q[1];
rz(-1.6575529) q[1];
sx q[1];
rz(3.0220093) q[1];
x q[2];
rz(0.32834239) q[3];
sx q[3];
rz(-0.36188618) q[3];
sx q[3];
rz(-0.34235024) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.0067979) q[2];
sx q[2];
rz(-3.1289913) q[2];
sx q[2];
rz(1.4726144) q[2];
rz(2.2115479) q[3];
sx q[3];
rz(-0.01447066) q[3];
sx q[3];
rz(-0.84021935) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(-pi/2) q[3];
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
rz(1.8188266) q[0];
sx q[0];
rz(-0.023094026) q[0];
sx q[0];
rz(-1.4267138) q[0];
rz(2.4115883) q[1];
sx q[1];
rz(-0.58861029) q[1];
sx q[1];
rz(-1.1013365) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5663213) q[0];
sx q[0];
rz(-2.3744517) q[0];
sx q[0];
rz(-1.9643892) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9726102) q[2];
sx q[2];
rz(-1.3436396) q[2];
sx q[2];
rz(-0.81534602) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.28336542) q[1];
sx q[1];
rz(-1.6943329) q[1];
sx q[1];
rz(-3.0459411) q[1];
rz(-pi) q[2];
rz(-1.6953129) q[3];
sx q[3];
rz(-1.6027556) q[3];
sx q[3];
rz(-2.4886142) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.4250028) q[2];
sx q[2];
rz(-0.060881946) q[2];
sx q[2];
rz(-1.3072183) q[2];
rz(2.763125) q[3];
sx q[3];
rz(-3.1186447) q[3];
sx q[3];
rz(-0.65346658) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8398447) q[0];
sx q[0];
rz(-1.2675588) q[0];
sx q[0];
rz(-0.92754716) q[0];
rz(1.357366) q[1];
sx q[1];
rz(-0.83186847) q[1];
sx q[1];
rz(1.6102788) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6283693) q[0];
sx q[0];
rz(-1.5484865) q[0];
sx q[0];
rz(-1.5324027) q[0];
rz(-pi) q[1];
rz(-0.39674098) q[2];
sx q[2];
rz(-1.3441685) q[2];
sx q[2];
rz(0.47725783) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.5591874) q[1];
sx q[1];
rz(-1.6998763) q[1];
sx q[1];
rz(-3.1388793) q[1];
rz(-pi) q[2];
rz(-2.82656) q[3];
sx q[3];
rz(-2.0612392) q[3];
sx q[3];
rz(-0.88446188) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.7808468) q[2];
sx q[2];
rz(-3.1372034) q[2];
sx q[2];
rz(-1.8878262) q[2];
rz(0.69416657) q[3];
sx q[3];
rz(-0.73918754) q[3];
sx q[3];
rz(-2.9085801) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
rz(-2.6041782) q[0];
sx q[0];
rz(-1.0056714) q[0];
sx q[0];
rz(-1.0373212) q[0];
rz(-1.6088156) q[1];
sx q[1];
rz(-0.2205801) q[1];
sx q[1];
rz(-1.467009) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5418422) q[0];
sx q[0];
rz(-1.6381725) q[0];
sx q[0];
rz(1.3753433) q[0];
rz(-pi) q[1];
x q[1];
rz(0.27811082) q[2];
sx q[2];
rz(-1.5928972) q[2];
sx q[2];
rz(1.2985566) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.92865366) q[1];
sx q[1];
rz(-0.0012081971) q[1];
sx q[1];
rz(1.1986284) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6653196) q[3];
sx q[3];
rz(-2.0243939) q[3];
sx q[3];
rz(-1.9840553) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(3.1303225) q[2];
sx q[2];
rz(-0.20877561) q[2];
sx q[2];
rz(-3.0976963) q[2];
rz(-2.6210426) q[3];
sx q[3];
rz(-0.0046516727) q[3];
sx q[3];
rz(1.4588149) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3540102) q[0];
sx q[0];
rz(-0.0024604877) q[0];
sx q[0];
rz(1.822923) q[0];
rz(-1.7240546) q[1];
sx q[1];
rz(-2.8520165) q[1];
sx q[1];
rz(-1.5444548) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.8108687) q[0];
sx q[0];
rz(-1.6858584) q[0];
sx q[0];
rz(0.48145357) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.6575281) q[2];
sx q[2];
rz(-1.3597466) q[2];
sx q[2];
rz(-1.576265) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.456157) q[1];
sx q[1];
rz(-2.7815869) q[1];
sx q[1];
rz(2.092157) q[1];
rz(-2.4539976) q[3];
sx q[3];
rz(-1.15117) q[3];
sx q[3];
rz(-0.18292038) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.80457193) q[2];
sx q[2];
rz(-1.2995517) q[2];
sx q[2];
rz(-2.9447832) q[2];
rz(-1.9491516) q[3];
sx q[3];
rz(-0.20659031) q[3];
sx q[3];
rz(0.20096745) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.30109677) q[0];
sx q[0];
rz(-1.7844642) q[0];
sx q[0];
rz(-1.1916196) q[0];
rz(-1.5246897) q[1];
sx q[1];
rz(-2.4947417) q[1];
sx q[1];
rz(1.5764538) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3197446) q[0];
sx q[0];
rz(-0.34043202) q[0];
sx q[0];
rz(-2.6186187) q[0];
x q[1];
rz(-2.0247518) q[2];
sx q[2];
rz(-1.364214) q[2];
sx q[2];
rz(-1.1517186) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.79542613) q[1];
sx q[1];
rz(-0.0010716575) q[1];
sx q[1];
rz(2.6571855) q[1];
x q[2];
rz(-3.0281046) q[3];
sx q[3];
rz(-1.7120512) q[3];
sx q[3];
rz(-3.0761513) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.9578751) q[2];
sx q[2];
rz(-0.58791939) q[2];
sx q[2];
rz(1.6831762) q[2];
rz(-0.030473907) q[3];
sx q[3];
rz(-3.1320429) q[3];
sx q[3];
rz(2.9392346) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4771582) q[0];
sx q[0];
rz(-1.3344593) q[0];
sx q[0];
rz(1.6819171) q[0];
rz(-1.5741813) q[1];
sx q[1];
rz(-1.3290783) q[1];
sx q[1];
rz(-3.0507416) q[1];
rz(3.136619) q[2];
sx q[2];
rz(-1.4952352) q[2];
sx q[2];
rz(0.16005439) q[2];
rz(2.1021765) q[3];
sx q[3];
rz(-2.3134091) q[3];
sx q[3];
rz(-0.29534657) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
