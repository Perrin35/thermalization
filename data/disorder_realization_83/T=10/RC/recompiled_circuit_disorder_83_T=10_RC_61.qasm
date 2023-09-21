OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.5542334) q[0];
sx q[0];
rz(4.1242546) q[0];
sx q[0];
rz(10.186515) q[0];
rz(0.76454437) q[1];
sx q[1];
rz(4.2188797) q[1];
sx q[1];
rz(10.168434) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8964891) q[0];
sx q[0];
rz(-1.6365543) q[0];
sx q[0];
rz(1.8114281) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.6526297) q[2];
sx q[2];
rz(-1.4748117) q[2];
sx q[2];
rz(-1.3053615) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.88120645) q[1];
sx q[1];
rz(-1.8477866) q[1];
sx q[1];
rz(-1.9828412) q[1];
x q[2];
rz(-1.2067912) q[3];
sx q[3];
rz(-1.5036821) q[3];
sx q[3];
rz(1.8458927) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.4347697) q[2];
sx q[2];
rz(-2.7396024) q[2];
sx q[2];
rz(-0.10786954) q[2];
rz(-2.9989631) q[3];
sx q[3];
rz(-1.7248036) q[3];
sx q[3];
rz(0.67255783) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0649439) q[0];
sx q[0];
rz(-2.3363484) q[0];
sx q[0];
rz(2.9108677) q[0];
rz(1.327286) q[1];
sx q[1];
rz(-0.67115152) q[1];
sx q[1];
rz(3.1006295) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.165867) q[0];
sx q[0];
rz(-1.3250933) q[0];
sx q[0];
rz(2.9742572) q[0];
rz(-1.2626921) q[2];
sx q[2];
rz(-1.8381881) q[2];
sx q[2];
rz(-1.845713) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.4119271) q[1];
sx q[1];
rz(-1.1150868) q[1];
sx q[1];
rz(-2.5018442) q[1];
x q[2];
rz(-1.1702234) q[3];
sx q[3];
rz(-1.922732) q[3];
sx q[3];
rz(2.7924201) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.9791947) q[2];
sx q[2];
rz(-1.4731864) q[2];
sx q[2];
rz(-2.5734148) q[2];
rz(0.5125106) q[3];
sx q[3];
rz(-0.5957225) q[3];
sx q[3];
rz(-0.56186831) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7115962) q[0];
sx q[0];
rz(-1.5189518) q[0];
sx q[0];
rz(0.45021737) q[0];
rz(1.846116) q[1];
sx q[1];
rz(-1.1643012) q[1];
sx q[1];
rz(0.67726642) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.30401858) q[0];
sx q[0];
rz(-1.7390334) q[0];
sx q[0];
rz(-1.5181245) q[0];
rz(1.5718939) q[2];
sx q[2];
rz(-1.564431) q[2];
sx q[2];
rz(0.15417834) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.22563572) q[1];
sx q[1];
rz(-1.3924053) q[1];
sx q[1];
rz(3.008328) q[1];
rz(-pi) q[2];
rz(2.4694091) q[3];
sx q[3];
rz(-0.58478343) q[3];
sx q[3];
rz(-0.72090805) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.9230817) q[2];
sx q[2];
rz(-1.9873025) q[2];
sx q[2];
rz(-3.0947321) q[2];
rz(2.3299407) q[3];
sx q[3];
rz(-0.67957687) q[3];
sx q[3];
rz(2.9714382) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.99455225) q[0];
sx q[0];
rz(-3.0299598) q[0];
sx q[0];
rz(0.051483367) q[0];
rz(-0.4908081) q[1];
sx q[1];
rz(-0.97508109) q[1];
sx q[1];
rz(1.1725918) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6068891) q[0];
sx q[0];
rz(-2.4692315) q[0];
sx q[0];
rz(2.9260203) q[0];
rz(0.02247359) q[2];
sx q[2];
rz(-0.60612504) q[2];
sx q[2];
rz(-1.7200574) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.1662894) q[1];
sx q[1];
rz(-1.2180093) q[1];
sx q[1];
rz(0.3198448) q[1];
rz(-pi) q[2];
rz(-1.4106393) q[3];
sx q[3];
rz(-0.43678108) q[3];
sx q[3];
rz(0.76524759) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.88901687) q[2];
sx q[2];
rz(-2.2106407) q[2];
sx q[2];
rz(2.6546997) q[2];
rz(1.9481109) q[3];
sx q[3];
rz(-0.32961696) q[3];
sx q[3];
rz(-0.025432767) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.837773) q[0];
sx q[0];
rz(-2.1210414) q[0];
sx q[0];
rz(1.3254962) q[0];
rz(-0.59108132) q[1];
sx q[1];
rz(-2.4609844) q[1];
sx q[1];
rz(2.4868734) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7703169) q[0];
sx q[0];
rz(-1.9195942) q[0];
sx q[0];
rz(0.47229292) q[0];
x q[1];
rz(-0.40277092) q[2];
sx q[2];
rz(-2.6908814) q[2];
sx q[2];
rz(-0.41926256) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.7655189) q[1];
sx q[1];
rz(-2.3507833) q[1];
sx q[1];
rz(-2.5942624) q[1];
x q[2];
rz(-2.914364) q[3];
sx q[3];
rz(-2.7471746) q[3];
sx q[3];
rz(-1.7928746) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.4513662) q[2];
sx q[2];
rz(-1.2008685) q[2];
sx q[2];
rz(2.5704685) q[2];
rz(0.58978224) q[3];
sx q[3];
rz(-0.46001205) q[3];
sx q[3];
rz(0.067972876) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1968483) q[0];
sx q[0];
rz(-1.1905043) q[0];
sx q[0];
rz(-2.8425472) q[0];
rz(1.3202745) q[1];
sx q[1];
rz(-2.8838005) q[1];
sx q[1];
rz(-1.4978131) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2346674) q[0];
sx q[0];
rz(-1.5186131) q[0];
sx q[0];
rz(2.232057) q[0];
x q[1];
rz(-2.3408893) q[2];
sx q[2];
rz(-0.75468894) q[2];
sx q[2];
rz(-0.18037361) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.3469997) q[1];
sx q[1];
rz(-1.1113249) q[1];
sx q[1];
rz(0.32230349) q[1];
x q[2];
rz(-1.3214345) q[3];
sx q[3];
rz(-1.8421838) q[3];
sx q[3];
rz(2.4621778) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.6254639) q[2];
sx q[2];
rz(-1.7136145) q[2];
sx q[2];
rz(-3.1141172) q[2];
rz(-0.52250683) q[3];
sx q[3];
rz(-2.3404739) q[3];
sx q[3];
rz(2.3251422) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7291173) q[0];
sx q[0];
rz(-2.0232047) q[0];
sx q[0];
rz(0.11418848) q[0];
rz(2.1633637) q[1];
sx q[1];
rz(-0.54662919) q[1];
sx q[1];
rz(0.79089975) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.233477) q[0];
sx q[0];
rz(-0.95196264) q[0];
sx q[0];
rz(3.0254362) q[0];
rz(-pi) q[1];
rz(1.9060471) q[2];
sx q[2];
rz(-2.2860048) q[2];
sx q[2];
rz(0.4244948) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.0930209) q[1];
sx q[1];
rz(-1.934811) q[1];
sx q[1];
rz(-0.084333468) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7780276) q[3];
sx q[3];
rz(-2.225038) q[3];
sx q[3];
rz(3.0594861) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.2699282) q[2];
sx q[2];
rz(-2.0285138) q[2];
sx q[2];
rz(-2.2765735) q[2];
rz(-0.67251742) q[3];
sx q[3];
rz(-1.1285684) q[3];
sx q[3];
rz(0.095656693) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6487811) q[0];
sx q[0];
rz(-2.6219941) q[0];
sx q[0];
rz(-2.6742324) q[0];
rz(2.6043747) q[1];
sx q[1];
rz(-0.98058128) q[1];
sx q[1];
rz(2.8875202) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.92914903) q[0];
sx q[0];
rz(-1.7663167) q[0];
sx q[0];
rz(1.6385727) q[0];
rz(-pi) q[1];
x q[1];
rz(3.1363856) q[2];
sx q[2];
rz(-1.569283) q[2];
sx q[2];
rz(3.0031799) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.64995631) q[1];
sx q[1];
rz(-1.2600113) q[1];
sx q[1];
rz(1.7040571) q[1];
rz(-pi) q[2];
rz(-0.59928008) q[3];
sx q[3];
rz(-1.5044971) q[3];
sx q[3];
rz(1.146572) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.49999985) q[2];
sx q[2];
rz(-2.3708512) q[2];
sx q[2];
rz(2.8059778) q[2];
rz(2.8149758) q[3];
sx q[3];
rz(-0.89544046) q[3];
sx q[3];
rz(1.4183104) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
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
rz(-0.083387233) q[0];
sx q[0];
rz(-0.21033062) q[0];
sx q[0];
rz(1.0539508) q[0];
rz(-2.9846233) q[1];
sx q[1];
rz(-1.7228246) q[1];
sx q[1];
rz(0.98186791) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1218402) q[0];
sx q[0];
rz(-1.4593908) q[0];
sx q[0];
rz(-0.41282546) q[0];
x q[1];
rz(-1.0018714) q[2];
sx q[2];
rz(-1.47942) q[2];
sx q[2];
rz(2.5779974) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.30584221) q[1];
sx q[1];
rz(-1.7562477) q[1];
sx q[1];
rz(-0.058538392) q[1];
rz(2.4307067) q[3];
sx q[3];
rz(-1.7611739) q[3];
sx q[3];
rz(0.058660942) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.92423576) q[2];
sx q[2];
rz(-2.0841667) q[2];
sx q[2];
rz(-2.1441933) q[2];
rz(-2.635397) q[3];
sx q[3];
rz(-2.180407) q[3];
sx q[3];
rz(0.034328073) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85703325) q[0];
sx q[0];
rz(-0.71145809) q[0];
sx q[0];
rz(2.6249264) q[0];
rz(0.38756469) q[1];
sx q[1];
rz(-1.060408) q[1];
sx q[1];
rz(0.26836747) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3049406) q[0];
sx q[0];
rz(-2.1450844) q[0];
sx q[0];
rz(0.71026295) q[0];
x q[1];
rz(-2.7029413) q[2];
sx q[2];
rz(-1.029656) q[2];
sx q[2];
rz(1.3676608) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.8473709) q[1];
sx q[1];
rz(-0.98681824) q[1];
sx q[1];
rz(1.2482615) q[1];
rz(-pi) q[2];
rz(1.1770505) q[3];
sx q[3];
rz(-0.6150133) q[3];
sx q[3];
rz(-2.7494489) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.2051852) q[2];
sx q[2];
rz(-2.3582017) q[2];
sx q[2];
rz(0.56383413) q[2];
rz(-1.9514203) q[3];
sx q[3];
rz(-2.1824013) q[3];
sx q[3];
rz(-0.64175516) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.664809) q[0];
sx q[0];
rz(-1.8905147) q[0];
sx q[0];
rz(2.074194) q[0];
rz(-1.8019567) q[1];
sx q[1];
rz(-1.4383153) q[1];
sx q[1];
rz(-1.7972606) q[1];
rz(1.2673169) q[2];
sx q[2];
rz(-2.0694642) q[2];
sx q[2];
rz(2.4533761) q[2];
rz(-1.3700804) q[3];
sx q[3];
rz(-1.2643391) q[3];
sx q[3];
rz(-1.5872965) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
