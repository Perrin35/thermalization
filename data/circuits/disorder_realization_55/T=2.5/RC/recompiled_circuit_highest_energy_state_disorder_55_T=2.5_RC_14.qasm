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
rz(0.16297451) q[1];
sx q[1];
rz(-1.4459223) q[1];
sx q[1];
rz(0.016782848) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.34384511) q[0];
sx q[0];
rz(-0.36917403) q[0];
sx q[0];
rz(1.7789715) q[0];
x q[1];
rz(-1.4799825) q[2];
sx q[2];
rz(-0.1802643) q[2];
sx q[2];
rz(-1.4992876) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.3455194) q[1];
sx q[1];
rz(-1.5658448) q[1];
sx q[1];
rz(1.5657182) q[1];
rz(-1.5161022) q[3];
sx q[3];
rz(-1.7229041) q[3];
sx q[3];
rz(1.5285672) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.5590543) q[2];
sx q[2];
rz(-2.4808919) q[2];
sx q[2];
rz(1.5622697) q[2];
rz(-2.9301379) q[3];
sx q[3];
rz(-0.00051694218) q[3];
sx q[3];
rz(-2.9805984) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.529539) q[0];
sx q[0];
rz(-0.27612975) q[0];
sx q[0];
rz(-1.3268693) q[0];
rz(-2.5621085) q[1];
sx q[1];
rz(-0.0038298413) q[1];
sx q[1];
rz(0.63900596) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20637437) q[0];
sx q[0];
rz(-2.2274096) q[0];
sx q[0];
rz(-0.73120631) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5875568) q[2];
sx q[2];
rz(-1.4498561) q[2];
sx q[2];
rz(3.1230833) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.7182448) q[1];
sx q[1];
rz(-3.1210174) q[1];
sx q[1];
rz(-2.5431376) q[1];
x q[2];
rz(-2.3611446) q[3];
sx q[3];
rz(-1.5089581) q[3];
sx q[3];
rz(1.058418) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.4352033) q[2];
sx q[2];
rz(-3.0053164) q[2];
sx q[2];
rz(-1.5380247) q[2];
rz(1.5806574) q[3];
sx q[3];
rz(-3.1272562) q[3];
sx q[3];
rz(-3.1106136) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8921709) q[0];
sx q[0];
rz(-0.51512655) q[0];
sx q[0];
rz(-2.7601335) q[0];
rz(2.4341266) q[1];
sx q[1];
rz(-3.1222157) q[1];
sx q[1];
rz(-1.1245419) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8302001) q[0];
sx q[0];
rz(-1.8233577) q[0];
sx q[0];
rz(0.24858944) q[0];
rz(-pi) q[1];
rz(1.3493645) q[2];
sx q[2];
rz(-3.0237758) q[2];
sx q[2];
rz(0.066514579) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.8360236) q[1];
sx q[1];
rz(-1.5828805) q[1];
sx q[1];
rz(-0.062968465) q[1];
rz(-pi) q[2];
rz(0.74060969) q[3];
sx q[3];
rz(-0.69878529) q[3];
sx q[3];
rz(-1.6388338) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.6996998) q[2];
sx q[2];
rz(-3.1293588) q[2];
sx q[2];
rz(-3.0921248) q[2];
rz(0.60702819) q[3];
sx q[3];
rz(-0.0012461239) q[3];
sx q[3];
rz(-1.9342669) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1599051) q[0];
sx q[0];
rz(-2.9732381) q[0];
sx q[0];
rz(-3.1244151) q[0];
rz(-0.29302868) q[1];
sx q[1];
rz(-2.3511062) q[1];
sx q[1];
rz(1.5471829) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7538504) q[0];
sx q[0];
rz(-1.004129) q[0];
sx q[0];
rz(-2.500534) q[0];
rz(-2.7780159) q[2];
sx q[2];
rz(-2.5695317) q[2];
sx q[2];
rz(-0.043302082) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.35369998) q[1];
sx q[1];
rz(-3.0038634) q[1];
sx q[1];
rz(0.44868033) q[1];
rz(-pi) q[2];
rz(1.5805575) q[3];
sx q[3];
rz(-0.92484821) q[3];
sx q[3];
rz(-0.6432337) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.8881417) q[2];
sx q[2];
rz(-0.47051045) q[2];
sx q[2];
rz(-2.4593501) q[2];
rz(-0.055179723) q[3];
sx q[3];
rz(-0.0076871593) q[3];
sx q[3];
rz(-1.8365708) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76160112) q[0];
sx q[0];
rz(-2.70607) q[0];
sx q[0];
rz(-0.4250266) q[0];
rz(1.60166) q[1];
sx q[1];
rz(-0.48307499) q[1];
sx q[1];
rz(2.3262598) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20629932) q[0];
sx q[0];
rz(-1.5816798) q[0];
sx q[0];
rz(-3.0803922) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5475531) q[2];
sx q[2];
rz(-1.574062) q[2];
sx q[2];
rz(-0.68319418) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.968135) q[1];
sx q[1];
rz(-1.4840398) q[1];
sx q[1];
rz(-0.11958338) q[1];
x q[2];
rz(-0.32834239) q[3];
sx q[3];
rz(-2.7797065) q[3];
sx q[3];
rz(2.7992424) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.0067979) q[2];
sx q[2];
rz(-3.1289913) q[2];
sx q[2];
rz(-1.4726144) q[2];
rz(2.2115479) q[3];
sx q[3];
rz(-3.127122) q[3];
sx q[3];
rz(0.84021935) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3227661) q[0];
sx q[0];
rz(-3.1184986) q[0];
sx q[0];
rz(-1.4267138) q[0];
rz(2.4115883) q[1];
sx q[1];
rz(-2.5529824) q[1];
sx q[1];
rz(-2.0402562) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5663213) q[0];
sx q[0];
rz(-2.3744517) q[0];
sx q[0];
rz(1.9643892) q[0];
rz(-2.1997994) q[2];
sx q[2];
rz(-0.28224361) q[2];
sx q[2];
rz(1.6779225) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.8423398) q[1];
sx q[1];
rz(-1.6657167) q[1];
sx q[1];
rz(1.4466982) q[1];
rz(1.6953129) q[3];
sx q[3];
rz(-1.6027556) q[3];
sx q[3];
rz(2.4886142) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.4250028) q[2];
sx q[2];
rz(-3.0807107) q[2];
sx q[2];
rz(-1.3072183) q[2];
rz(-2.763125) q[3];
sx q[3];
rz(-3.1186447) q[3];
sx q[3];
rz(-2.4881261) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3017479) q[0];
sx q[0];
rz(-1.8740338) q[0];
sx q[0];
rz(-2.2140455) q[0];
rz(-1.357366) q[1];
sx q[1];
rz(-0.83186847) q[1];
sx q[1];
rz(1.5313139) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6283693) q[0];
sx q[0];
rz(-1.5931061) q[0];
sx q[0];
rz(1.5324027) q[0];
rz(-2.6035735) q[2];
sx q[2];
rz(-2.6876861) q[2];
sx q[2];
rz(1.5557784) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.011259638) q[1];
sx q[1];
rz(-1.5734871) q[1];
sx q[1];
rz(1.6998768) q[1];
x q[2];
rz(-2.0824984) q[3];
sx q[3];
rz(-1.8476433) q[3];
sx q[3];
rz(0.83864318) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.36074582) q[2];
sx q[2];
rz(-0.0043892269) q[2];
sx q[2];
rz(1.2537664) q[2];
rz(-2.4474261) q[3];
sx q[3];
rz(-0.73918754) q[3];
sx q[3];
rz(0.23301253) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6041782) q[0];
sx q[0];
rz(-2.1359213) q[0];
sx q[0];
rz(-1.0373212) q[0];
rz(-1.5327771) q[1];
sx q[1];
rz(-2.9210126) q[1];
sx q[1];
rz(-1.467009) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.64325414) q[0];
sx q[0];
rz(-0.20659978) q[0];
sx q[0];
rz(1.905196) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5478126) q[2];
sx q[2];
rz(-1.2927552) q[2];
sx q[2];
rz(2.8756623) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.3008219) q[1];
sx q[1];
rz(-1.5719218) q[1];
sx q[1];
rz(-0.00043934396) q[1];
rz(-1.476273) q[3];
sx q[3];
rz(-2.0243939) q[3];
sx q[3];
rz(-1.9840553) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.011270114) q[2];
sx q[2];
rz(-0.20877561) q[2];
sx q[2];
rz(3.0976963) q[2];
rz(-2.6210426) q[3];
sx q[3];
rz(-0.0046516727) q[3];
sx q[3];
rz(-1.6827778) q[3];
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
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7875824) q[0];
sx q[0];
rz(-3.1391322) q[0];
sx q[0];
rz(1.3186697) q[0];
rz(1.4175381) q[1];
sx q[1];
rz(-0.28957614) q[1];
sx q[1];
rz(-1.5971378) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.97623551) q[0];
sx q[0];
rz(-0.49396038) q[0];
sx q[0];
rz(-0.24458347) q[0];
rz(-pi) q[1];
rz(2.4840646) q[2];
sx q[2];
rz(-1.3597466) q[2];
sx q[2];
rz(1.5653277) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.456157) q[1];
sx q[1];
rz(-2.7815869) q[1];
sx q[1];
rz(2.092157) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.68759509) q[3];
sx q[3];
rz(-1.15117) q[3];
sx q[3];
rz(-2.9586723) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.3370207) q[2];
sx q[2];
rz(-1.8420409) q[2];
sx q[2];
rz(-0.19680944) q[2];
rz(1.9491516) q[3];
sx q[3];
rz(-2.9350023) q[3];
sx q[3];
rz(0.20096745) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8404959) q[0];
sx q[0];
rz(-1.7844642) q[0];
sx q[0];
rz(-1.949973) q[0];
rz(-1.5246897) q[1];
sx q[1];
rz(-2.4947417) q[1];
sx q[1];
rz(1.5764538) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25116205) q[0];
sx q[0];
rz(-1.7383456) q[0];
sx q[0];
rz(-2.8438389) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9124969) q[2];
sx q[2];
rz(-2.014403) q[2];
sx q[2];
rz(2.6227621) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.79542613) q[1];
sx q[1];
rz(-0.0010716575) q[1];
sx q[1];
rz(-2.6571855) q[1];
rz(-pi) q[2];
x q[2];
rz(0.11348806) q[3];
sx q[3];
rz(-1.7120512) q[3];
sx q[3];
rz(-3.0761513) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.18371753) q[2];
sx q[2];
rz(-2.5536733) q[2];
sx q[2];
rz(1.6831762) q[2];
rz(-0.030473907) q[3];
sx q[3];
rz(-0.009549791) q[3];
sx q[3];
rz(0.20235801) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4771582) q[0];
sx q[0];
rz(-1.8071334) q[0];
sx q[0];
rz(-1.4596756) q[0];
rz(-1.5674113) q[1];
sx q[1];
rz(-1.8125143) q[1];
sx q[1];
rz(0.090851091) q[1];
rz(1.5051928) q[2];
sx q[2];
rz(-3.0658683) q[2];
sx q[2];
rz(0.22584596) q[2];
rz(-2.6371757) q[3];
sx q[3];
rz(-0.88263369) q[3];
sx q[3];
rz(0.42019444) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
