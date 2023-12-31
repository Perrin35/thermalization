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
rz(-1.5313907) q[1];
sx q[1];
rz(-2.664497) q[1];
sx q[1];
rz(2.6452126) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3408605) q[0];
sx q[0];
rz(-0.98156089) q[0];
sx q[0];
rz(0.13829921) q[0];
x q[1];
rz(-1.7726937) q[2];
sx q[2];
rz(-1.9041833) q[2];
sx q[2];
rz(-0.29793973) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.5254933) q[1];
sx q[1];
rz(-1.4800737) q[1];
sx q[1];
rz(2.4331122) q[1];
rz(-2.7339488) q[3];
sx q[3];
rz(-1.4789494) q[3];
sx q[3];
rz(2.7013005) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.51497841) q[2];
sx q[2];
rz(-0.86003059) q[2];
sx q[2];
rz(-2.853945) q[2];
rz(-1.3927762) q[3];
sx q[3];
rz(-0.98373047) q[3];
sx q[3];
rz(-2.0387409) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-2.2709133) q[0];
sx q[0];
rz(-2.3864855) q[0];
sx q[0];
rz(-0.57759181) q[0];
rz(1.6060991) q[1];
sx q[1];
rz(-1.1178733) q[1];
sx q[1];
rz(-1.5637406) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0751511) q[0];
sx q[0];
rz(-1.6834714) q[0];
sx q[0];
rz(-0.22002797) q[0];
rz(-pi) q[1];
rz(-2.4096476) q[2];
sx q[2];
rz(-0.39790301) q[2];
sx q[2];
rz(2.0741472) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.8375081) q[1];
sx q[1];
rz(-1.4834852) q[1];
sx q[1];
rz(1.1156032) q[1];
rz(-pi) q[2];
x q[2];
rz(1.413826) q[3];
sx q[3];
rz(-1.7236712) q[3];
sx q[3];
rz(-0.65519858) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.10721283) q[2];
sx q[2];
rz(-2.1652174) q[2];
sx q[2];
rz(-1.0822901) q[2];
rz(-1.9783431) q[3];
sx q[3];
rz(-1.2000368) q[3];
sx q[3];
rz(2.5015586) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
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
rz(-1.0796233) q[0];
sx q[0];
rz(-2.847023) q[0];
sx q[0];
rz(0.18297718) q[0];
rz(0.043116365) q[1];
sx q[1];
rz(-2.1886107) q[1];
sx q[1];
rz(1.144369) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9632918) q[0];
sx q[0];
rz(-2.460647) q[0];
sx q[0];
rz(-2.2292024) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6687713) q[2];
sx q[2];
rz(-2.678215) q[2];
sx q[2];
rz(-2.9544427) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.3414587) q[1];
sx q[1];
rz(-2.027249) q[1];
sx q[1];
rz(1.1851428) q[1];
rz(-3.1312917) q[3];
sx q[3];
rz(-0.28763887) q[3];
sx q[3];
rz(-2.0623296) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.019505067) q[2];
sx q[2];
rz(-1.959789) q[2];
sx q[2];
rz(-2.2369177) q[2];
rz(-0.71980643) q[3];
sx q[3];
rz(-0.42680877) q[3];
sx q[3];
rz(-2.3864746) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
x q[3];
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
rz(-1.4633789) q[0];
sx q[0];
rz(-2.1713874) q[0];
sx q[0];
rz(0.71587193) q[0];
rz(-2.8158358) q[1];
sx q[1];
rz(-1.0595067) q[1];
sx q[1];
rz(0.99266565) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9216114) q[0];
sx q[0];
rz(-2.5489759) q[0];
sx q[0];
rz(-0.48595925) q[0];
rz(-0.9348346) q[2];
sx q[2];
rz(-2.0419428) q[2];
sx q[2];
rz(2.3061894) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.48651931) q[1];
sx q[1];
rz(-1.7336978) q[1];
sx q[1];
rz(-1.6652813) q[1];
rz(2.2654815) q[3];
sx q[3];
rz(-0.70385859) q[3];
sx q[3];
rz(2.0640404) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.0830393) q[2];
sx q[2];
rz(-0.71648592) q[2];
sx q[2];
rz(-1.7129718) q[2];
rz(-0.84609091) q[3];
sx q[3];
rz(-1.5139791) q[3];
sx q[3];
rz(-1.2231474) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59584004) q[0];
sx q[0];
rz(-1.7359474) q[0];
sx q[0];
rz(-0.76675057) q[0];
rz(-1.2402361) q[1];
sx q[1];
rz(-2.2940472) q[1];
sx q[1];
rz(0.3516745) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36067097) q[0];
sx q[0];
rz(-2.43649) q[0];
sx q[0];
rz(2.6893925) q[0];
rz(-pi) q[1];
x q[1];
rz(2.211116) q[2];
sx q[2];
rz(-0.74379197) q[2];
sx q[2];
rz(-0.1375246) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.1406447) q[1];
sx q[1];
rz(-1.5854156) q[1];
sx q[1];
rz(0.39477243) q[1];
rz(-2.8204927) q[3];
sx q[3];
rz(-2.0541111) q[3];
sx q[3];
rz(-1.6881642) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.19796431) q[2];
sx q[2];
rz(-1.6683656) q[2];
sx q[2];
rz(0.53331214) q[2];
rz(-0.42896459) q[3];
sx q[3];
rz(-2.6140489) q[3];
sx q[3];
rz(-1.126948) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0241942) q[0];
sx q[0];
rz(-2.0485931) q[0];
sx q[0];
rz(0.54164106) q[0];
rz(-2.533124) q[1];
sx q[1];
rz(-0.21921961) q[1];
sx q[1];
rz(-2.0419962) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0352286) q[0];
sx q[0];
rz(-1.117825) q[0];
sx q[0];
rz(0.41082541) q[0];
rz(2.1018283) q[2];
sx q[2];
rz(-1.4073939) q[2];
sx q[2];
rz(2.642717) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.82367831) q[1];
sx q[1];
rz(-2.4412529) q[1];
sx q[1];
rz(0.27842303) q[1];
x q[2];
rz(-0.78824708) q[3];
sx q[3];
rz(-1.3282093) q[3];
sx q[3];
rz(-1.5188252) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.6301443) q[2];
sx q[2];
rz(-1.5254598) q[2];
sx q[2];
rz(0.28298322) q[2];
rz(0.7061559) q[3];
sx q[3];
rz(-1.599584) q[3];
sx q[3];
rz(2.506536) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0866078) q[0];
sx q[0];
rz(-0.98751706) q[0];
sx q[0];
rz(-1.6116066) q[0];
rz(2.4967172) q[1];
sx q[1];
rz(-0.49352831) q[1];
sx q[1];
rz(-0.15730102) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10282117) q[0];
sx q[0];
rz(-0.73909315) q[0];
sx q[0];
rz(3.0022013) q[0];
rz(-pi) q[1];
x q[1];
rz(0.21643164) q[2];
sx q[2];
rz(-2.2955403) q[2];
sx q[2];
rz(0.12491465) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(3.0485222) q[1];
sx q[1];
rz(-1.7659566) q[1];
sx q[1];
rz(-2.1919769) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5141684) q[3];
sx q[3];
rz(-1.3324696) q[3];
sx q[3];
rz(2.0605007) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.1823696) q[2];
sx q[2];
rz(-2.5490641) q[2];
sx q[2];
rz(-2.3106993) q[2];
rz(3.0977141) q[3];
sx q[3];
rz(-1.9079804) q[3];
sx q[3];
rz(-2.9390745) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6458994) q[0];
sx q[0];
rz(-0.93911397) q[0];
sx q[0];
rz(3.1307401) q[0];
rz(-2.8649578) q[1];
sx q[1];
rz(-0.77181548) q[1];
sx q[1];
rz(2.8483134) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.54942185) q[0];
sx q[0];
rz(-2.5991419) q[0];
sx q[0];
rz(0.011499238) q[0];
rz(-pi) q[1];
rz(0.92708712) q[2];
sx q[2];
rz(-1.4378528) q[2];
sx q[2];
rz(-1.5333652) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.3775023) q[1];
sx q[1];
rz(-1.1524156) q[1];
sx q[1];
rz(-0.80470316) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6777612) q[3];
sx q[3];
rz(-1.2004735) q[3];
sx q[3];
rz(0.90255373) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.04348065) q[2];
sx q[2];
rz(-3.0467693) q[2];
sx q[2];
rz(0.088767178) q[2];
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
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1373238) q[0];
sx q[0];
rz(-0.22432888) q[0];
sx q[0];
rz(1.0639169) q[0];
rz(0.44257277) q[1];
sx q[1];
rz(-1.4241709) q[1];
sx q[1];
rz(1.7907422) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1254743) q[0];
sx q[0];
rz(-0.11575143) q[0];
sx q[0];
rz(2.3257757) q[0];
rz(-1.6532134) q[2];
sx q[2];
rz(-0.91225183) q[2];
sx q[2];
rz(1.272162) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.0693448) q[1];
sx q[1];
rz(-2.7069271) q[1];
sx q[1];
rz(1.5321295) q[1];
rz(2.001858) q[3];
sx q[3];
rz(-2.8981879) q[3];
sx q[3];
rz(-0.36153015) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.0315447) q[2];
sx q[2];
rz(-1.2908547) q[2];
sx q[2];
rz(-2.6943977) q[2];
rz(-1.7556919) q[3];
sx q[3];
rz(-2.8409676) q[3];
sx q[3];
rz(-0.51469222) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8266325) q[0];
sx q[0];
rz(-1.6199912) q[0];
sx q[0];
rz(0.78053027) q[0];
rz(0.42778095) q[1];
sx q[1];
rz(-0.3573187) q[1];
sx q[1];
rz(2.9719877) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0376301) q[0];
sx q[0];
rz(-0.61015218) q[0];
sx q[0];
rz(-0.9141586) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5461966) q[2];
sx q[2];
rz(-1.149017) q[2];
sx q[2];
rz(1.3299024) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.5055713) q[1];
sx q[1];
rz(-2.1850056) q[1];
sx q[1];
rz(-0.057227055) q[1];
rz(-2.9986936) q[3];
sx q[3];
rz(-1.2522962) q[3];
sx q[3];
rz(-2.3058476) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.85049373) q[2];
sx q[2];
rz(-1.1824181) q[2];
sx q[2];
rz(0.69236857) q[2];
rz(2.678357) q[3];
sx q[3];
rz(-1.654947) q[3];
sx q[3];
rz(2.0326116) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2904084) q[0];
sx q[0];
rz(-1.6155227) q[0];
sx q[0];
rz(0.68646705) q[0];
rz(-2.6295173) q[1];
sx q[1];
rz(-0.33437406) q[1];
sx q[1];
rz(-2.1127111) q[1];
rz(-2.8725273) q[2];
sx q[2];
rz(-1.8620327) q[2];
sx q[2];
rz(-3.0897683) q[2];
rz(-1.3626171) q[3];
sx q[3];
rz(-1.6856185) q[3];
sx q[3];
rz(0.57193397) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
