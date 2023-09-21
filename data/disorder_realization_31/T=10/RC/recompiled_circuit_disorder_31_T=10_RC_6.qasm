OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.6150317) q[0];
sx q[0];
rz(2.5685413) q[0];
sx q[0];
rz(11.723784) q[0];
rz(2.1057582) q[1];
sx q[1];
rz(8.3254568) q[1];
sx q[1];
rz(7.96666) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5320839) q[0];
sx q[0];
rz(-2.3713787) q[0];
sx q[0];
rz(-0.30755933) q[0];
x q[1];
rz(0.61383944) q[2];
sx q[2];
rz(-1.5868574) q[2];
sx q[2];
rz(1.2889372) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.7588501) q[1];
sx q[1];
rz(-1.9470125) q[1];
sx q[1];
rz(-1.7099027) q[1];
rz(-2.5472766) q[3];
sx q[3];
rz(-1.672097) q[3];
sx q[3];
rz(3.1241824) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.6136916) q[2];
sx q[2];
rz(-2.1353022) q[2];
sx q[2];
rz(-0.17949417) q[2];
rz(1.9159296) q[3];
sx q[3];
rz(-1.7951199) q[3];
sx q[3];
rz(2.3195482) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
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
rz(-2.3935788) q[0];
sx q[0];
rz(-0.8809692) q[0];
sx q[0];
rz(0.32546145) q[0];
rz(-1.7851967) q[1];
sx q[1];
rz(-2.0929095) q[1];
sx q[1];
rz(-1.9869841) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.017529537) q[0];
sx q[0];
rz(-1.5520099) q[0];
sx q[0];
rz(1.5536874) q[0];
x q[1];
rz(-2.7484659) q[2];
sx q[2];
rz(-0.98193411) q[2];
sx q[2];
rz(1.4002422) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.5408389) q[1];
sx q[1];
rz(-2.3327017) q[1];
sx q[1];
rz(-1.6879338) q[1];
rz(-2.1784337) q[3];
sx q[3];
rz(-2.8249486) q[3];
sx q[3];
rz(-0.39099993) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.6894199) q[2];
sx q[2];
rz(-1.8916811) q[2];
sx q[2];
rz(-0.88341218) q[2];
rz(-2.6702821) q[3];
sx q[3];
rz(-1.703197) q[3];
sx q[3];
rz(0.78770351) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31323355) q[0];
sx q[0];
rz(-1.4947083) q[0];
sx q[0];
rz(-1.6261684) q[0];
rz(-0.60107636) q[1];
sx q[1];
rz(-2.5939012) q[1];
sx q[1];
rz(2.0498958) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1356782) q[0];
sx q[0];
rz(-1.0251097) q[0];
sx q[0];
rz(-0.90555993) q[0];
rz(-2.0902363) q[2];
sx q[2];
rz(-0.72548496) q[2];
sx q[2];
rz(1.6348334) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.4746689) q[1];
sx q[1];
rz(-0.83819929) q[1];
sx q[1];
rz(-1.0679507) q[1];
x q[2];
rz(-2.6767119) q[3];
sx q[3];
rz(-2.998623) q[3];
sx q[3];
rz(-0.79899341) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.320257) q[2];
sx q[2];
rz(-2.6358423) q[2];
sx q[2];
rz(2.2606405) q[2];
rz(-1.3736003) q[3];
sx q[3];
rz(-1.526984) q[3];
sx q[3];
rz(-1.0176456) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3110733) q[0];
sx q[0];
rz(-1.7493462) q[0];
sx q[0];
rz(-0.4367035) q[0];
rz(-0.23315915) q[1];
sx q[1];
rz(-1.2522839) q[1];
sx q[1];
rz(-2.8312347) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6757641) q[0];
sx q[0];
rz(-2.7390263) q[0];
sx q[0];
rz(-0.34253828) q[0];
rz(2.4565059) q[2];
sx q[2];
rz(-1.6685467) q[2];
sx q[2];
rz(-0.098066559) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.0536641) q[1];
sx q[1];
rz(-1.2128608) q[1];
sx q[1];
rz(-0.019836516) q[1];
rz(3.1060018) q[3];
sx q[3];
rz(-1.4471874) q[3];
sx q[3];
rz(-1.3328758) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.13005304) q[2];
sx q[2];
rz(-2.4235642) q[2];
sx q[2];
rz(1.0774353) q[2];
rz(-3.0854026) q[3];
sx q[3];
rz(-2.5037933) q[3];
sx q[3];
rz(-1.594054) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.189165) q[0];
sx q[0];
rz(-2.0963033) q[0];
sx q[0];
rz(-2.8919343) q[0];
rz(-1.5769618) q[1];
sx q[1];
rz(-2.3639634) q[1];
sx q[1];
rz(0.87019428) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91709671) q[0];
sx q[0];
rz(-1.774569) q[0];
sx q[0];
rz(-2.8208371) q[0];
rz(1.7587897) q[2];
sx q[2];
rz(-0.88289875) q[2];
sx q[2];
rz(-0.57304136) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.40904564) q[1];
sx q[1];
rz(-2.2605719) q[1];
sx q[1];
rz(1.9802666) q[1];
rz(-1.0088483) q[3];
sx q[3];
rz(-1.4268488) q[3];
sx q[3];
rz(2.2465003) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.27328086) q[2];
sx q[2];
rz(-1.8153278) q[2];
sx q[2];
rz(2.4678521) q[2];
rz(-0.30361787) q[3];
sx q[3];
rz(-1.9165336) q[3];
sx q[3];
rz(-1.822086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7917787) q[0];
sx q[0];
rz(-0.93739167) q[0];
sx q[0];
rz(-2.8836024) q[0];
rz(2.7164283) q[1];
sx q[1];
rz(-2.185967) q[1];
sx q[1];
rz(1.4917096) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8404322) q[0];
sx q[0];
rz(-1.4727117) q[0];
sx q[0];
rz(-0.43099404) q[0];
rz(1.9906524) q[2];
sx q[2];
rz(-0.93517762) q[2];
sx q[2];
rz(0.90204001) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.3867492) q[1];
sx q[1];
rz(-2.1347087) q[1];
sx q[1];
rz(-0.90653231) q[1];
x q[2];
rz(1.2061938) q[3];
sx q[3];
rz(-1.4758665) q[3];
sx q[3];
rz(-1.7942384) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.012718) q[2];
sx q[2];
rz(-0.96367633) q[2];
sx q[2];
rz(3.0498665) q[2];
rz(2.2979459) q[3];
sx q[3];
rz(-0.97976145) q[3];
sx q[3];
rz(2.2475524) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3180852) q[0];
sx q[0];
rz(-1.2473236) q[0];
sx q[0];
rz(-2.7303625) q[0];
rz(2.2757018) q[1];
sx q[1];
rz(-0.31232467) q[1];
sx q[1];
rz(0.033989865) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.54995178) q[0];
sx q[0];
rz(-0.80054987) q[0];
sx q[0];
rz(-0.22649015) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2074239) q[2];
sx q[2];
rz(-1.120943) q[2];
sx q[2];
rz(-2.5069782) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.9719203) q[1];
sx q[1];
rz(-1.5202513) q[1];
sx q[1];
rz(-2.2488942) q[1];
x q[2];
rz(2.9339318) q[3];
sx q[3];
rz(-2.5163979) q[3];
sx q[3];
rz(3.0561662) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.5380481) q[2];
sx q[2];
rz(-0.59843439) q[2];
sx q[2];
rz(-0.87654385) q[2];
rz(0.34902469) q[3];
sx q[3];
rz(-1.9411496) q[3];
sx q[3];
rz(-0.14311895) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2440764) q[0];
sx q[0];
rz(-1.709047) q[0];
sx q[0];
rz(0.39392719) q[0];
rz(0.36755964) q[1];
sx q[1];
rz(-1.7575248) q[1];
sx q[1];
rz(-1.4454909) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1720393) q[0];
sx q[0];
rz(-1.5656099) q[0];
sx q[0];
rz(2.0041549) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1922853) q[2];
sx q[2];
rz(-1.0770123) q[2];
sx q[2];
rz(-2.2311503) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.1576924) q[1];
sx q[1];
rz(-1.0219814) q[1];
sx q[1];
rz(1.7556612) q[1];
rz(-pi) q[2];
rz(2.3950855) q[3];
sx q[3];
rz(-2.3254447) q[3];
sx q[3];
rz(1.5796803) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.8470856) q[2];
sx q[2];
rz(-0.89035788) q[2];
sx q[2];
rz(-2.7344446) q[2];
rz(1.5173222) q[3];
sx q[3];
rz(-1.9842691) q[3];
sx q[3];
rz(2.8919162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.80609926) q[0];
sx q[0];
rz(-2.6265916) q[0];
sx q[0];
rz(-1.2517713) q[0];
rz(2.4720526) q[1];
sx q[1];
rz(-1.1839097) q[1];
sx q[1];
rz(2.8318185) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.6471841) q[0];
sx q[0];
rz(-1.697288) q[0];
sx q[0];
rz(2.6148318) q[0];
x q[1];
rz(-1.6236213) q[2];
sx q[2];
rz(-2.9644358) q[2];
sx q[2];
rz(2.0668541) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.579638) q[1];
sx q[1];
rz(-1.6288174) q[1];
sx q[1];
rz(1.3647563) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8272607) q[3];
sx q[3];
rz(-0.98955742) q[3];
sx q[3];
rz(0.47282156) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.70242515) q[2];
sx q[2];
rz(-2.4283786) q[2];
sx q[2];
rz(-1.2072198) q[2];
rz(-1.0369982) q[3];
sx q[3];
rz(-1.8959277) q[3];
sx q[3];
rz(-2.4859378) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0913775) q[0];
sx q[0];
rz(-1.8176879) q[0];
sx q[0];
rz(-1.2058831) q[0];
rz(2.5559015) q[1];
sx q[1];
rz(-1.0810477) q[1];
sx q[1];
rz(1.4996128) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8382032) q[0];
sx q[0];
rz(-1.2801542) q[0];
sx q[0];
rz(-3.035726) q[0];
rz(-pi) q[1];
x q[1];
rz(1.502938) q[2];
sx q[2];
rz(-0.80791622) q[2];
sx q[2];
rz(0.18005904) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.75913945) q[1];
sx q[1];
rz(-1.1083974) q[1];
sx q[1];
rz(-2.9049302) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3979982) q[3];
sx q[3];
rz(-1.5732592) q[3];
sx q[3];
rz(-0.60009225) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.2575834) q[2];
sx q[2];
rz(-1.3486226) q[2];
sx q[2];
rz(-2.1949027) q[2];
rz(-0.36869129) q[3];
sx q[3];
rz(-1.5654516) q[3];
sx q[3];
rz(-0.45599109) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7832227) q[0];
sx q[0];
rz(-1.1483648) q[0];
sx q[0];
rz(-0.4233465) q[0];
rz(-3.070667) q[1];
sx q[1];
rz(-1.6880886) q[1];
sx q[1];
rz(-0.26500519) q[1];
rz(0.82540853) q[2];
sx q[2];
rz(-0.63612973) q[2];
sx q[2];
rz(-0.3542184) q[2];
rz(2.4214217) q[3];
sx q[3];
rz(-1.1838893) q[3];
sx q[3];
rz(0.67223687) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
