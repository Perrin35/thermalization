OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.6168851) q[0];
sx q[0];
rz(-0.73862139) q[0];
sx q[0];
rz(0.11697669) q[0];
rz(0.95247954) q[1];
sx q[1];
rz(7.7533819) q[1];
sx q[1];
rz(8.5511313) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0067209816) q[0];
sx q[0];
rz(-2.3804579) q[0];
sx q[0];
rz(1.8675682) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.9175542) q[2];
sx q[2];
rz(-0.54346701) q[2];
sx q[2];
rz(-2.562491) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.80109875) q[1];
sx q[1];
rz(-1.5135161) q[1];
sx q[1];
rz(-1.0216453) q[1];
rz(-1.9756894) q[3];
sx q[3];
rz(-1.3171853) q[3];
sx q[3];
rz(3.1133127) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.17254193) q[2];
sx q[2];
rz(-1.8092864) q[2];
sx q[2];
rz(-1.2782028) q[2];
rz(-1.605002) q[3];
sx q[3];
rz(-1.2383818) q[3];
sx q[3];
rz(-0.31755519) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69739598) q[0];
sx q[0];
rz(-2.9228656) q[0];
sx q[0];
rz(2.2636407) q[0];
rz(-0.76389337) q[1];
sx q[1];
rz(-1.0822252) q[1];
sx q[1];
rz(2.0108932) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5754841) q[0];
sx q[0];
rz(-0.16666767) q[0];
sx q[0];
rz(-1.9420293) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.133059) q[2];
sx q[2];
rz(-1.4302876) q[2];
sx q[2];
rz(2.2214691) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.89040989) q[1];
sx q[1];
rz(-1.5580388) q[1];
sx q[1];
rz(0.81752543) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.14826397) q[3];
sx q[3];
rz(-2.5103323) q[3];
sx q[3];
rz(0.12693044) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.6905602) q[2];
sx q[2];
rz(-1.6310383) q[2];
sx q[2];
rz(2.4859599) q[2];
rz(1.9111274) q[3];
sx q[3];
rz(-2.1793607) q[3];
sx q[3];
rz(2.3088096) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3967628) q[0];
sx q[0];
rz(-1.9704882) q[0];
sx q[0];
rz(2.6387264) q[0];
rz(2.9961131) q[1];
sx q[1];
rz(-1.611404) q[1];
sx q[1];
rz(1.1605638) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7877951) q[0];
sx q[0];
rz(-0.43913865) q[0];
sx q[0];
rz(-0.77180688) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.34973379) q[2];
sx q[2];
rz(-1.4473414) q[2];
sx q[2];
rz(0.48179212) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.8050011) q[1];
sx q[1];
rz(-2.1865926) q[1];
sx q[1];
rz(-2.6686882) q[1];
x q[2];
rz(-0.3480556) q[3];
sx q[3];
rz(-1.8249434) q[3];
sx q[3];
rz(-2.4055907) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.9460556) q[2];
sx q[2];
rz(-2.1058407) q[2];
sx q[2];
rz(2.3717086) q[2];
rz(-2.5464673) q[3];
sx q[3];
rz(-2.6502521) q[3];
sx q[3];
rz(-2.6814521) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6339517) q[0];
sx q[0];
rz(-2.0946298) q[0];
sx q[0];
rz(-1.7858343) q[0];
rz(-0.62375623) q[1];
sx q[1];
rz(-1.535894) q[1];
sx q[1];
rz(2.6502868) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76440109) q[0];
sx q[0];
rz(-2.1024414) q[0];
sx q[0];
rz(1.3926943) q[0];
rz(-pi) q[1];
rz(1.7801966) q[2];
sx q[2];
rz(-1.2970957) q[2];
sx q[2];
rz(-1.4865246) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.2196893) q[1];
sx q[1];
rz(-1.8796726) q[1];
sx q[1];
rz(0.83654735) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0912498) q[3];
sx q[3];
rz(-2.1698275) q[3];
sx q[3];
rz(1.1831621) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.19646159) q[2];
sx q[2];
rz(-1.6170231) q[2];
sx q[2];
rz(2.5932236) q[2];
rz(1.8114926) q[3];
sx q[3];
rz(-2.8434704) q[3];
sx q[3];
rz(0.33838457) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2746628) q[0];
sx q[0];
rz(-0.45329705) q[0];
sx q[0];
rz(1.0484265) q[0];
rz(0.10920814) q[1];
sx q[1];
rz(-1.5501153) q[1];
sx q[1];
rz(1.106326) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7503742) q[0];
sx q[0];
rz(-2.8958909) q[0];
sx q[0];
rz(-1.1121763) q[0];
x q[1];
rz(0.84594251) q[2];
sx q[2];
rz(-1.3157613) q[2];
sx q[2];
rz(-0.77470335) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.26004654) q[1];
sx q[1];
rz(-1.4235919) q[1];
sx q[1];
rz(-1.8974426) q[1];
rz(-pi) q[2];
x q[2];
rz(2.901793) q[3];
sx q[3];
rz(-1.2811525) q[3];
sx q[3];
rz(1.0847131) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3722374) q[2];
sx q[2];
rz(-1.2662338) q[2];
sx q[2];
rz(0.20277578) q[2];
rz(2.3619704) q[3];
sx q[3];
rz(-2.0519665) q[3];
sx q[3];
rz(2.7705418) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5146273) q[0];
sx q[0];
rz(-3.096464) q[0];
sx q[0];
rz(-2.2623862) q[0];
rz(-2.5458287) q[1];
sx q[1];
rz(-0.87409449) q[1];
sx q[1];
rz(-1.2557868) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7065053) q[0];
sx q[0];
rz(-1.5179978) q[0];
sx q[0];
rz(-3.0493124) q[0];
x q[1];
rz(-2.6379616) q[2];
sx q[2];
rz(-2.5066895) q[2];
sx q[2];
rz(-2.0201928) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.75460282) q[1];
sx q[1];
rz(-0.55596724) q[1];
sx q[1];
rz(-2.703859) q[1];
rz(-pi) q[2];
rz(-1.8196657) q[3];
sx q[3];
rz(-0.96967317) q[3];
sx q[3];
rz(2.7271276) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.6229728) q[2];
sx q[2];
rz(-2.8039248) q[2];
sx q[2];
rz(-1.6597623) q[2];
rz(1.4431813) q[3];
sx q[3];
rz(-1.3866164) q[3];
sx q[3];
rz(1.1186918) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5959394) q[0];
sx q[0];
rz(-2.312199) q[0];
sx q[0];
rz(0.19350061) q[0];
rz(-0.045348383) q[1];
sx q[1];
rz(-1.1646611) q[1];
sx q[1];
rz(-2.4605816) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4880002) q[0];
sx q[0];
rz(-2.9569599) q[0];
sx q[0];
rz(-0.49166481) q[0];
rz(-0.25410507) q[2];
sx q[2];
rz(-2.0850344) q[2];
sx q[2];
rz(2.6387173) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.633887) q[1];
sx q[1];
rz(-1.5649598) q[1];
sx q[1];
rz(-1.5836442) q[1];
rz(-2.7549823) q[3];
sx q[3];
rz(-1.4550687) q[3];
sx q[3];
rz(-0.9936617) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.27651522) q[2];
sx q[2];
rz(-0.64537185) q[2];
sx q[2];
rz(1.3481677) q[2];
rz(1.3639785) q[3];
sx q[3];
rz(-1.8230702) q[3];
sx q[3];
rz(1.9510423) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6337223) q[0];
sx q[0];
rz(-0.81955925) q[0];
sx q[0];
rz(2.4622397) q[0];
rz(-0.23718111) q[1];
sx q[1];
rz(-2.2229384) q[1];
sx q[1];
rz(2.0749345) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6045881) q[0];
sx q[0];
rz(-0.51031546) q[0];
sx q[0];
rz(-3.133581) q[0];
rz(-2.1751498) q[2];
sx q[2];
rz(-1.4928515) q[2];
sx q[2];
rz(-1.9029531) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.28344) q[1];
sx q[1];
rz(-1.4464708) q[1];
sx q[1];
rz(-1.0470301) q[1];
rz(-pi) q[2];
rz(0.32853029) q[3];
sx q[3];
rz(-1.9764998) q[3];
sx q[3];
rz(-0.87268396) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.25881585) q[2];
sx q[2];
rz(-1.9100185) q[2];
sx q[2];
rz(-2.7302177) q[2];
rz(-1.2808293) q[3];
sx q[3];
rz(-0.6461834) q[3];
sx q[3];
rz(-1.0584186) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48788747) q[0];
sx q[0];
rz(-1.3478841) q[0];
sx q[0];
rz(0.93061289) q[0];
rz(-2.6764684) q[1];
sx q[1];
rz(-0.5828751) q[1];
sx q[1];
rz(1.4177657) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.055061888) q[0];
sx q[0];
rz(-1.3320005) q[0];
sx q[0];
rz(-2.5361674) q[0];
x q[1];
rz(-0.64196824) q[2];
sx q[2];
rz(-0.99364181) q[2];
sx q[2];
rz(-0.3954487) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.90365961) q[1];
sx q[1];
rz(-1.2301908) q[1];
sx q[1];
rz(-1.8800859) q[1];
x q[2];
rz(1.3950138) q[3];
sx q[3];
rz(-0.58628288) q[3];
sx q[3];
rz(2.0678584) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.174939) q[2];
sx q[2];
rz(-0.83028364) q[2];
sx q[2];
rz(1.9261599) q[2];
rz(-2.4231966) q[3];
sx q[3];
rz(-1.4677488) q[3];
sx q[3];
rz(-2.4992656) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4833118) q[0];
sx q[0];
rz(-0.42228666) q[0];
sx q[0];
rz(-1.8102113) q[0];
rz(-2.4347958) q[1];
sx q[1];
rz(-1.4247318) q[1];
sx q[1];
rz(1.7165064) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9909716) q[0];
sx q[0];
rz(-2.2036834) q[0];
sx q[0];
rz(-1.0101684) q[0];
x q[1];
rz(1.6548884) q[2];
sx q[2];
rz(-2.483209) q[2];
sx q[2];
rz(0.7158196) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.6687732) q[1];
sx q[1];
rz(-2.4309078) q[1];
sx q[1];
rz(0.70699923) q[1];
rz(-2.6985907) q[3];
sx q[3];
rz(-2.2878929) q[3];
sx q[3];
rz(-2.9505961) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.66296545) q[2];
sx q[2];
rz(-0.68296432) q[2];
sx q[2];
rz(1.9662201) q[2];
rz(0.023224467) q[3];
sx q[3];
rz(-2.569779) q[3];
sx q[3];
rz(2.8779252) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5671253) q[0];
sx q[0];
rz(-0.99526417) q[0];
sx q[0];
rz(-1.9010726) q[0];
rz(0.65470882) q[1];
sx q[1];
rz(-1.8617102) q[1];
sx q[1];
rz(-1.1690232) q[1];
rz(-3.0289247) q[2];
sx q[2];
rz(-1.7074399) q[2];
sx q[2];
rz(1.8497333) q[2];
rz(0.8620048) q[3];
sx q[3];
rz(-1.0222407) q[3];
sx q[3];
rz(0.23289451) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
