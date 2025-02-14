OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.1333753) q[0];
sx q[0];
rz(-1.8741338) q[0];
sx q[0];
rz(-3.128669) q[0];
rz(-2.456993) q[1];
sx q[1];
rz(-0.79889387) q[1];
sx q[1];
rz(2.0838783) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7803181) q[0];
sx q[0];
rz(-2.4343505) q[0];
sx q[0];
rz(0.35298423) q[0];
rz(2.0796258) q[2];
sx q[2];
rz(-0.86039174) q[2];
sx q[2];
rz(1.8343385) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.8243421) q[1];
sx q[1];
rz(-1.2095272) q[1];
sx q[1];
rz(-2.3995705) q[1];
rz(-1.9625447) q[3];
sx q[3];
rz(-1.6101675) q[3];
sx q[3];
rz(1.5848499) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.58618033) q[2];
sx q[2];
rz(-1.5427898) q[2];
sx q[2];
rz(3.0482698) q[2];
rz(0.020545067) q[3];
sx q[3];
rz(-2.8829657) q[3];
sx q[3];
rz(-1.3625905) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.071863197) q[0];
sx q[0];
rz(-1.7427895) q[0];
sx q[0];
rz(-0.82759696) q[0];
rz(0.96356511) q[1];
sx q[1];
rz(-1.5255442) q[1];
sx q[1];
rz(-2.4049984) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6927836) q[0];
sx q[0];
rz(-2.8570456) q[0];
sx q[0];
rz(1.297516) q[0];
rz(0.57666333) q[2];
sx q[2];
rz(-2.0915589) q[2];
sx q[2];
rz(-1.7885838) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.5399038) q[1];
sx q[1];
rz(-1.772953) q[1];
sx q[1];
rz(-1.8784932) q[1];
rz(1.5932036) q[3];
sx q[3];
rz(-1.7051892) q[3];
sx q[3];
rz(-0.16772863) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.5654512) q[2];
sx q[2];
rz(-0.57500035) q[2];
sx q[2];
rz(-0.74419332) q[2];
rz(-0.60892504) q[3];
sx q[3];
rz(-2.3600793) q[3];
sx q[3];
rz(-1.5003381) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3610483) q[0];
sx q[0];
rz(-2.2760976) q[0];
sx q[0];
rz(1.3737099) q[0];
rz(-2.3420077) q[1];
sx q[1];
rz(-2.1345963) q[1];
sx q[1];
rz(1.132157) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.06057418) q[0];
sx q[0];
rz(-1.3558398) q[0];
sx q[0];
rz(-0.092553986) q[0];
rz(-0.27081174) q[2];
sx q[2];
rz(-0.083148227) q[2];
sx q[2];
rz(-1.6390334) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.4375648) q[1];
sx q[1];
rz(-1.4933741) q[1];
sx q[1];
rz(0.67005007) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6145094) q[3];
sx q[3];
rz(-0.96230405) q[3];
sx q[3];
rz(-0.19245806) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.75480294) q[2];
sx q[2];
rz(-0.58480442) q[2];
sx q[2];
rz(-1.7542138) q[2];
rz(0.34902188) q[3];
sx q[3];
rz(-1.6882221) q[3];
sx q[3];
rz(0.52946985) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.741852) q[0];
sx q[0];
rz(-2.1735503) q[0];
sx q[0];
rz(-2.7834748) q[0];
rz(2.070836) q[1];
sx q[1];
rz(-2.4929969) q[1];
sx q[1];
rz(1.4395641) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1308129) q[0];
sx q[0];
rz(-2.3284916) q[0];
sx q[0];
rz(1.2050425) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1522572) q[2];
sx q[2];
rz(-2.3609567) q[2];
sx q[2];
rz(2.6065741) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.52709748) q[1];
sx q[1];
rz(-2.1367837) q[1];
sx q[1];
rz(-0.66733349) q[1];
x q[2];
rz(-2.8791588) q[3];
sx q[3];
rz(-0.78426266) q[3];
sx q[3];
rz(-1.5465496) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.7122571) q[2];
sx q[2];
rz(-1.2306932) q[2];
sx q[2];
rz(0.62475359) q[2];
rz(1.6338927) q[3];
sx q[3];
rz(-0.75993901) q[3];
sx q[3];
rz(-2.0190575) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
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
rz(2.7655012) q[0];
sx q[0];
rz(-2.2608345) q[0];
sx q[0];
rz(1.1464024) q[0];
rz(1.435185) q[1];
sx q[1];
rz(-2.621666) q[1];
sx q[1];
rz(0.82040876) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.7026742) q[0];
sx q[0];
rz(-1.6587388) q[0];
sx q[0];
rz(-1.8431208) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9632513) q[2];
sx q[2];
rz(-2.2296411) q[2];
sx q[2];
rz(-0.73769157) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.0443521) q[1];
sx q[1];
rz(-0.28438452) q[1];
sx q[1];
rz(-2.3769955) q[1];
rz(-pi) q[2];
rz(-0.58973226) q[3];
sx q[3];
rz(-1.5431649) q[3];
sx q[3];
rz(-0.10527912) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.72941214) q[2];
sx q[2];
rz(-1.0554487) q[2];
sx q[2];
rz(-2.6030276) q[2];
rz(-1.4696848) q[3];
sx q[3];
rz(-1.4476176) q[3];
sx q[3];
rz(-3.0830429) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3013714) q[0];
sx q[0];
rz(-3.0026307) q[0];
sx q[0];
rz(3.050991) q[0];
rz(2.7719356) q[1];
sx q[1];
rz(-1.5934817) q[1];
sx q[1];
rz(0.37818092) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2205296) q[0];
sx q[0];
rz(-2.605633) q[0];
sx q[0];
rz(-0.45951636) q[0];
x q[1];
rz(-1.5289076) q[2];
sx q[2];
rz(-1.48078) q[2];
sx q[2];
rz(2.5652094) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.0269811) q[1];
sx q[1];
rz(-2.0617194) q[1];
sx q[1];
rz(1.4308962) q[1];
x q[2];
rz(1.1797253) q[3];
sx q[3];
rz(-2.9357233) q[3];
sx q[3];
rz(0.34032527) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.6211264) q[2];
sx q[2];
rz(-1.4755321) q[2];
sx q[2];
rz(-3.0401518) q[2];
rz(1.8488047) q[3];
sx q[3];
rz(-2.3088876) q[3];
sx q[3];
rz(1.4147991) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(1.8349649) q[0];
sx q[0];
rz(-1.2766301) q[0];
sx q[0];
rz(2.2391338) q[0];
rz(0.2392256) q[1];
sx q[1];
rz(-1.3419515) q[1];
sx q[1];
rz(-2.7801524) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7976686) q[0];
sx q[0];
rz(-1.6822681) q[0];
sx q[0];
rz(-1.072207) q[0];
rz(1.6267516) q[2];
sx q[2];
rz(-1.4437041) q[2];
sx q[2];
rz(-0.54886078) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.8351961) q[1];
sx q[1];
rz(-1.3521242) q[1];
sx q[1];
rz(1.1347215) q[1];
rz(-1.7291728) q[3];
sx q[3];
rz(-1.2353727) q[3];
sx q[3];
rz(0.80663825) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.07301894) q[2];
sx q[2];
rz(-1.34015) q[2];
sx q[2];
rz(-0.87693357) q[2];
rz(-2.1584611) q[3];
sx q[3];
rz(-1.5638331) q[3];
sx q[3];
rz(2.1985998) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8125732) q[0];
sx q[0];
rz(-0.1952157) q[0];
sx q[0];
rz(-2.3028497) q[0];
rz(2.4328531) q[1];
sx q[1];
rz(-2.510431) q[1];
sx q[1];
rz(-0.90604025) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2748286) q[0];
sx q[0];
rz(-1.5918918) q[0];
sx q[0];
rz(2.4167908) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9259813) q[2];
sx q[2];
rz(-1.2279358) q[2];
sx q[2];
rz(-0.19419369) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.69605322) q[1];
sx q[1];
rz(-1.0372682) q[1];
sx q[1];
rz(-0.75977709) q[1];
rz(-2.1791547) q[3];
sx q[3];
rz(-1.3624853) q[3];
sx q[3];
rz(2.359451) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.9771542) q[2];
sx q[2];
rz(-0.70455569) q[2];
sx q[2];
rz(-1.1158811) q[2];
rz(-2.5130533) q[3];
sx q[3];
rz(-1.081859) q[3];
sx q[3];
rz(-0.61454296) q[3];
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
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.81166613) q[0];
sx q[0];
rz(-2.3217432) q[0];
sx q[0];
rz(2.9803168) q[0];
rz(2.2992086) q[1];
sx q[1];
rz(-1.4295652) q[1];
sx q[1];
rz(-0.17002034) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5710167) q[0];
sx q[0];
rz(-1.5890772) q[0];
sx q[0];
rz(1.5789501) q[0];
x q[1];
rz(2.7964994) q[2];
sx q[2];
rz(-0.4640013) q[2];
sx q[2];
rz(0.74079266) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.0898897) q[1];
sx q[1];
rz(-2.0031628) q[1];
sx q[1];
rz(2.927455) q[1];
rz(-pi) q[2];
rz(1.5000072) q[3];
sx q[3];
rz(-0.38372358) q[3];
sx q[3];
rz(-0.078670382) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.1539479) q[2];
sx q[2];
rz(-1.5445856) q[2];
sx q[2];
rz(-1.6551931) q[2];
rz(-0.060128309) q[3];
sx q[3];
rz(-0.82258737) q[3];
sx q[3];
rz(1.4415584) q[3];
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
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.63115591) q[0];
sx q[0];
rz(-3.1102409) q[0];
sx q[0];
rz(1.7106868) q[0];
rz(2.778964) q[1];
sx q[1];
rz(-1.6094094) q[1];
sx q[1];
rz(-1.3154202) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5595488) q[0];
sx q[0];
rz(-1.8616195) q[0];
sx q[0];
rz(1.2469588) q[0];
x q[1];
rz(-3.0193826) q[2];
sx q[2];
rz(-1.9463681) q[2];
sx q[2];
rz(-1.5236601) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.9818287) q[1];
sx q[1];
rz(-1.8366645) q[1];
sx q[1];
rz(-2.5362316) q[1];
rz(-pi) q[2];
rz(2.682456) q[3];
sx q[3];
rz(-0.92184421) q[3];
sx q[3];
rz(-0.54455633) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.057283904) q[2];
sx q[2];
rz(-1.5949275) q[2];
sx q[2];
rz(0.16051897) q[2];
rz(-0.48216835) q[3];
sx q[3];
rz(-2.4653258) q[3];
sx q[3];
rz(2.4619861) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
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
rz(3.1228444) q[0];
sx q[0];
rz(-1.3958805) q[0];
sx q[0];
rz(-0.91176283) q[0];
rz(1.1625166) q[1];
sx q[1];
rz(-1.569842) q[1];
sx q[1];
rz(-0.40649489) q[1];
rz(-3.054259) q[2];
sx q[2];
rz(-1.6651911) q[2];
sx q[2];
rz(-1.390425) q[2];
rz(2.3707262) q[3];
sx q[3];
rz(-1.1419747) q[3];
sx q[3];
rz(2.9621073) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
