OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.0032089631) q[0];
sx q[0];
rz(-0.15455833) q[0];
sx q[0];
rz(0.69252339) q[0];
rz(1.9321631) q[1];
sx q[1];
rz(-1.2485319) q[1];
sx q[1];
rz(-1.385153) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6801075) q[0];
sx q[0];
rz(-0.84851096) q[0];
sx q[0];
rz(1.1397584) q[0];
rz(-pi) q[1];
rz(-2.7117549) q[2];
sx q[2];
rz(-2.5463383) q[2];
sx q[2];
rz(-2.0855479) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.79917158) q[1];
sx q[1];
rz(-0.28563269) q[1];
sx q[1];
rz(-0.61113961) q[1];
rz(-pi) q[2];
rz(2.4888943) q[3];
sx q[3];
rz(-1.1678809) q[3];
sx q[3];
rz(-0.17890113) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.8866855) q[2];
sx q[2];
rz(-2.343785) q[2];
sx q[2];
rz(2.936426) q[2];
rz(-2.3702879) q[3];
sx q[3];
rz(-2.3588534) q[3];
sx q[3];
rz(-2.0390959) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7339864) q[0];
sx q[0];
rz(-2.3953231) q[0];
sx q[0];
rz(-2.6876887) q[0];
rz(2.1167963) q[1];
sx q[1];
rz(-0.4075993) q[1];
sx q[1];
rz(-1.9143547) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5204029) q[0];
sx q[0];
rz(-1.495201) q[0];
sx q[0];
rz(-1.441342) q[0];
rz(-2.0298376) q[2];
sx q[2];
rz(-1.0978205) q[2];
sx q[2];
rz(0.74726653) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.8933567) q[1];
sx q[1];
rz(-1.3125988) q[1];
sx q[1];
rz(-0.30171079) q[1];
rz(-pi) q[2];
x q[2];
rz(0.20083986) q[3];
sx q[3];
rz(-1.4174995) q[3];
sx q[3];
rz(-2.4917345) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-3.1002905) q[2];
sx q[2];
rz(-1.9561448) q[2];
sx q[2];
rz(0.56742898) q[2];
rz(-0.36519095) q[3];
sx q[3];
rz(-1.4130211) q[3];
sx q[3];
rz(0.96810961) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.658618) q[0];
sx q[0];
rz(-2.5768319) q[0];
sx q[0];
rz(-2.2429402) q[0];
rz(-2.1458416) q[1];
sx q[1];
rz(-1.5581222) q[1];
sx q[1];
rz(-2.8083037) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.916327) q[0];
sx q[0];
rz(-0.81575459) q[0];
sx q[0];
rz(-2.4832721) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7477481) q[2];
sx q[2];
rz(-1.6672009) q[2];
sx q[2];
rz(-2.5071438) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.4132727) q[1];
sx q[1];
rz(-1.4188758) q[1];
sx q[1];
rz(0.9539714) q[1];
x q[2];
rz(2.1786147) q[3];
sx q[3];
rz(-0.64219785) q[3];
sx q[3];
rz(-1.9574788) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.4553392) q[2];
sx q[2];
rz(-1.3505961) q[2];
sx q[2];
rz(1.1509482) q[2];
rz(-2.3006556) q[3];
sx q[3];
rz(-2.0261814) q[3];
sx q[3];
rz(1.1545198) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6999917) q[0];
sx q[0];
rz(-1.561152) q[0];
sx q[0];
rz(2.4568795) q[0];
rz(1.0355863) q[1];
sx q[1];
rz(-2.6338449) q[1];
sx q[1];
rz(1.9365786) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.372615) q[0];
sx q[0];
rz(-0.87991558) q[0];
sx q[0];
rz(0.0432424) q[0];
rz(-pi) q[1];
rz(-1.3022468) q[2];
sx q[2];
rz(-2.4797202) q[2];
sx q[2];
rz(2.4115987) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.11338621) q[1];
sx q[1];
rz(-2.3909667) q[1];
sx q[1];
rz(-1.1651462) q[1];
rz(-3.0699176) q[3];
sx q[3];
rz(-2.274401) q[3];
sx q[3];
rz(3.049831) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.24923199) q[2];
sx q[2];
rz(-1.426733) q[2];
sx q[2];
rz(0.37115804) q[2];
rz(1.4012198) q[3];
sx q[3];
rz(-0.6597844) q[3];
sx q[3];
rz(-1.1192809) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.054759653) q[0];
sx q[0];
rz(-2.355447) q[0];
sx q[0];
rz(0.13312419) q[0];
rz(-0.99331028) q[1];
sx q[1];
rz(-1.3860093) q[1];
sx q[1];
rz(0.55508074) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.13984891) q[0];
sx q[0];
rz(-1.2554597) q[0];
sx q[0];
rz(0.01339162) q[0];
rz(1.1558756) q[2];
sx q[2];
rz(-1.5083815) q[2];
sx q[2];
rz(1.8780564) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.75025573) q[1];
sx q[1];
rz(-1.5515944) q[1];
sx q[1];
rz(2.0160497) q[1];
rz(-0.57226945) q[3];
sx q[3];
rz(-0.096147691) q[3];
sx q[3];
rz(-2.8807246) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.30620265) q[2];
sx q[2];
rz(-2.1357048) q[2];
sx q[2];
rz(0.13892697) q[2];
rz(-0.94240087) q[3];
sx q[3];
rz(-1.5709632) q[3];
sx q[3];
rz(2.5901103) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5979364) q[0];
sx q[0];
rz(-0.56977001) q[0];
sx q[0];
rz(0.57975769) q[0];
rz(0.12750164) q[1];
sx q[1];
rz(-1.9516877) q[1];
sx q[1];
rz(-1.5396083) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.72551661) q[0];
sx q[0];
rz(-2.2726739) q[0];
sx q[0];
rz(0.26728018) q[0];
rz(-pi) q[1];
rz(1.6164262) q[2];
sx q[2];
rz(-1.7012193) q[2];
sx q[2];
rz(-1.6861196) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.786799) q[1];
sx q[1];
rz(-1.8537632) q[1];
sx q[1];
rz(-1.1250886) q[1];
x q[2];
rz(1.2195915) q[3];
sx q[3];
rz(-0.71577365) q[3];
sx q[3];
rz(0.90482611) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.5876028) q[2];
sx q[2];
rz(-2.8911399) q[2];
sx q[2];
rz(2.8721151) q[2];
rz(-2.907471) q[3];
sx q[3];
rz(-2.6168489) q[3];
sx q[3];
rz(-3.0814734) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44678974) q[0];
sx q[0];
rz(-2.5233874) q[0];
sx q[0];
rz(-0.68429464) q[0];
rz(-0.11958312) q[1];
sx q[1];
rz(-1.2922829) q[1];
sx q[1];
rz(-0.51876846) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9414026) q[0];
sx q[0];
rz(-1.7046283) q[0];
sx q[0];
rz(-0.83394136) q[0];
rz(-1.7905551) q[2];
sx q[2];
rz(-2.3292654) q[2];
sx q[2];
rz(-2.1213558) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.24592933) q[1];
sx q[1];
rz(-0.83918011) q[1];
sx q[1];
rz(0.64114665) q[1];
x q[2];
rz(-1.5090844) q[3];
sx q[3];
rz(-0.40137526) q[3];
sx q[3];
rz(2.9369831) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.2542904) q[2];
sx q[2];
rz(-1.7742949) q[2];
sx q[2];
rz(2.5781412) q[2];
rz(3.0900132) q[3];
sx q[3];
rz(-2.0023465) q[3];
sx q[3];
rz(0.95190597) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
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
rz(0.96034399) q[0];
sx q[0];
rz(-2.612817) q[0];
sx q[0];
rz(-1.3990336) q[0];
rz(0.78701204) q[1];
sx q[1];
rz(-2.0128638) q[1];
sx q[1];
rz(2.3972437) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.23119152) q[0];
sx q[0];
rz(-1.4230799) q[0];
sx q[0];
rz(1.6258679) q[0];
x q[1];
rz(2.0853945) q[2];
sx q[2];
rz(-1.7040164) q[2];
sx q[2];
rz(2.104987) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.12606049) q[1];
sx q[1];
rz(-1.0711526) q[1];
sx q[1];
rz(-1.8263032) q[1];
x q[2];
rz(-3.0934422) q[3];
sx q[3];
rz(-2.6865494) q[3];
sx q[3];
rz(0.20850785) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.4259592) q[2];
sx q[2];
rz(-1.7461494) q[2];
sx q[2];
rz(0.60950935) q[2];
rz(0.65731796) q[3];
sx q[3];
rz(-2.4980563) q[3];
sx q[3];
rz(-2.8801584) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1881926) q[0];
sx q[0];
rz(-3.0472026) q[0];
sx q[0];
rz(-1.5040065) q[0];
rz(1.9001182) q[1];
sx q[1];
rz(-1.1499317) q[1];
sx q[1];
rz(2.3666568) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.06668815) q[0];
sx q[0];
rz(-2.3432891) q[0];
sx q[0];
rz(2.128128) q[0];
rz(-pi) q[1];
x q[1];
rz(0.1562642) q[2];
sx q[2];
rz(-2.0093577) q[2];
sx q[2];
rz(1.4807448) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.1800268) q[1];
sx q[1];
rz(-1.119009) q[1];
sx q[1];
rz(1.3614484) q[1];
rz(-1.7899412) q[3];
sx q[3];
rz(-1.431576) q[3];
sx q[3];
rz(-2.5185891) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.9541786) q[2];
sx q[2];
rz(-2.9118907) q[2];
sx q[2];
rz(2.9837218) q[2];
rz(-1.212451) q[3];
sx q[3];
rz(-2.082943) q[3];
sx q[3];
rz(1.2780179) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.071844) q[0];
sx q[0];
rz(-0.97244111) q[0];
sx q[0];
rz(0.21433314) q[0];
rz(-0.65746039) q[1];
sx q[1];
rz(-0.22413707) q[1];
sx q[1];
rz(1.0459895) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2803672) q[0];
sx q[0];
rz(-1.4645637) q[0];
sx q[0];
rz(-2.7809814) q[0];
rz(-2.5506053) q[2];
sx q[2];
rz(-1.5948442) q[2];
sx q[2];
rz(-1.9901333) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.18098772) q[1];
sx q[1];
rz(-1.029403) q[1];
sx q[1];
rz(-2.7156746) q[1];
rz(-1.5442185) q[3];
sx q[3];
rz(-1.1655032) q[3];
sx q[3];
rz(1.6534896) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.41632286) q[2];
sx q[2];
rz(-0.060083397) q[2];
sx q[2];
rz(-1.0894758) q[2];
rz(1.5661092) q[3];
sx q[3];
rz(-1.243306) q[3];
sx q[3];
rz(0.65264788) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2789223) q[0];
sx q[0];
rz(-0.60386064) q[0];
sx q[0];
rz(0.84558564) q[0];
rz(-1.5325585) q[1];
sx q[1];
rz(-1.4668203) q[1];
sx q[1];
rz(-1.1046881) q[1];
rz(-2.7325148) q[2];
sx q[2];
rz(-1.2862051) q[2];
sx q[2];
rz(-0.4551879) q[2];
rz(-1.5076751) q[3];
sx q[3];
rz(-2.1199385) q[3];
sx q[3];
rz(-1.4169823) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];