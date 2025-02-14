OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.92276031) q[0];
sx q[0];
rz(-3.0781167) q[0];
sx q[0];
rz(1.3702962) q[0];
rz(0.81075794) q[1];
sx q[1];
rz(1.6914565) q[1];
sx q[1];
rz(9.20426) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17973404) q[0];
sx q[0];
rz(-0.73886907) q[0];
sx q[0];
rz(-1.671285) q[0];
x q[1];
rz(0.96376586) q[2];
sx q[2];
rz(-2.6015601) q[2];
sx q[2];
rz(2.2390197) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-3.0444572) q[1];
sx q[1];
rz(-2.9191764) q[1];
sx q[1];
rz(3.0916721) q[1];
rz(-pi) q[2];
rz(-2.070921) q[3];
sx q[3];
rz(-0.4005188) q[3];
sx q[3];
rz(-2.5409215) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-3.0567131) q[2];
sx q[2];
rz(-0.006338174) q[2];
sx q[2];
rz(2.32178) q[2];
rz(3.1208755) q[3];
sx q[3];
rz(-2.8262704) q[3];
sx q[3];
rz(-2.0899541) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3725975) q[0];
sx q[0];
rz(-0.1854493) q[0];
sx q[0];
rz(2.4107667) q[0];
rz(1.7164879) q[1];
sx q[1];
rz(-2.4162879) q[1];
sx q[1];
rz(-2.6740429) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.088652164) q[0];
sx q[0];
rz(-3.0346617) q[0];
sx q[0];
rz(-1.8631975) q[0];
x q[1];
rz(-0.40981648) q[2];
sx q[2];
rz(-1.5404319) q[2];
sx q[2];
rz(-1.8259468) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.4380355) q[1];
sx q[1];
rz(-2.272637) q[1];
sx q[1];
rz(1.7662394) q[1];
rz(1.2456513) q[3];
sx q[3];
rz(-2.3811901) q[3];
sx q[3];
rz(-1.6154621) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.1656437) q[2];
sx q[2];
rz(-3.1303945) q[2];
sx q[2];
rz(-0.33640081) q[2];
rz(-2.8339556) q[3];
sx q[3];
rz(-0.010713723) q[3];
sx q[3];
rz(-0.78905869) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.12863185) q[0];
sx q[0];
rz(-2.9427981) q[0];
sx q[0];
rz(0.13191731) q[0];
rz(1.7031274) q[1];
sx q[1];
rz(-1.7273644) q[1];
sx q[1];
rz(-1.4836813) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0795201) q[0];
sx q[0];
rz(-1.8490845) q[0];
sx q[0];
rz(0.62836439) q[0];
x q[1];
rz(3.1171666) q[2];
sx q[2];
rz(-1.5533864) q[2];
sx q[2];
rz(2.8581564) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.7805788) q[1];
sx q[1];
rz(-3.0254361) q[1];
sx q[1];
rz(1.1067944) q[1];
rz(-pi) q[2];
x q[2];
rz(0.086016969) q[3];
sx q[3];
rz(-1.6413851) q[3];
sx q[3];
rz(-0.020078192) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.5710473) q[2];
sx q[2];
rz(-1.3344301) q[2];
sx q[2];
rz(1.4578311) q[2];
rz(2.3633862) q[3];
sx q[3];
rz(-0.0056548803) q[3];
sx q[3];
rz(-2.9290504) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0197765) q[0];
sx q[0];
rz(-1.7541405) q[0];
sx q[0];
rz(-0.29086599) q[0];
rz(-1.5485171) q[1];
sx q[1];
rz(-0.80268186) q[1];
sx q[1];
rz(-0.026550857) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.3817859) q[0];
sx q[0];
rz(-1.2273754) q[0];
sx q[0];
rz(-1.1939826) q[0];
rz(-pi) q[1];
x q[1];
rz(0.10992421) q[2];
sx q[2];
rz(-0.4532686) q[2];
sx q[2];
rz(-2.088415) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.0504405) q[1];
sx q[1];
rz(-2.181297) q[1];
sx q[1];
rz(-0.97368413) q[1];
x q[2];
rz(3.13843) q[3];
sx q[3];
rz(-1.589865) q[3];
sx q[3];
rz(-2.9657488) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.2619541) q[2];
sx q[2];
rz(-3.1231572) q[2];
sx q[2];
rz(-2.7035942) q[2];
rz(-1.1079463) q[3];
sx q[3];
rz(-0.017939311) q[3];
sx q[3];
rz(-1.7855135) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9172025) q[0];
sx q[0];
rz(-2.6438535) q[0];
sx q[0];
rz(-2.9955731) q[0];
rz(-0.13305013) q[1];
sx q[1];
rz(-1.0597884) q[1];
sx q[1];
rz(-2.6859247) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.21996343) q[0];
sx q[0];
rz(-1.5507924) q[0];
sx q[0];
rz(-3.0476493) q[0];
rz(-pi) q[1];
x q[1];
rz(0.59203573) q[2];
sx q[2];
rz(-1.8285654) q[2];
sx q[2];
rz(1.7532312) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.4670438) q[1];
sx q[1];
rz(-1.4321126) q[1];
sx q[1];
rz(-1.3402746) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.066448575) q[3];
sx q[3];
rz(-0.22804582) q[3];
sx q[3];
rz(-1.4888637) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.2666152) q[2];
sx q[2];
rz(-0.019860331) q[2];
sx q[2];
rz(1.2561426) q[2];
rz(0.52136326) q[3];
sx q[3];
rz(-0.25786906) q[3];
sx q[3];
rz(2.7969587) q[3];
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
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.94619036) q[0];
sx q[0];
rz(-2.7552216) q[0];
sx q[0];
rz(0.56889164) q[0];
rz(2.9733114) q[1];
sx q[1];
rz(-1.556309) q[1];
sx q[1];
rz(3.1313484) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.83547384) q[0];
sx q[0];
rz(-2.9635307) q[0];
sx q[0];
rz(2.8639069) q[0];
x q[1];
rz(-0.019316761) q[2];
sx q[2];
rz(-1.3903119) q[2];
sx q[2];
rz(1.3767124) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.3553487) q[1];
sx q[1];
rz(-3.0867483) q[1];
sx q[1];
rz(2.2516903) q[1];
rz(-0.63416173) q[3];
sx q[3];
rz(-1.595257) q[3];
sx q[3];
rz(-2.0195877) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.16113082) q[2];
sx q[2];
rz(-2.9176596) q[2];
sx q[2];
rz(1.0978318) q[2];
rz(2.791642) q[3];
sx q[3];
rz(-1.8643458) q[3];
sx q[3];
rz(-2.2866975) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2646645) q[0];
sx q[0];
rz(-2.9783037) q[0];
sx q[0];
rz(2.7352585) q[0];
rz(-1.8304652) q[1];
sx q[1];
rz(-2.6268112) q[1];
sx q[1];
rz(-3.0313671) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0144008) q[0];
sx q[0];
rz(-1.6699635) q[0];
sx q[0];
rz(-1.5983461) q[0];
rz(-pi) q[1];
rz(3.1378205) q[2];
sx q[2];
rz(-1.5654148) q[2];
sx q[2];
rz(-2.7710235) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.4496042) q[1];
sx q[1];
rz(-1.1432198) q[1];
sx q[1];
rz(-0.009593462) q[1];
rz(0.73506395) q[3];
sx q[3];
rz(-1.0197612) q[3];
sx q[3];
rz(-2.6363274) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.8608287) q[2];
sx q[2];
rz(-0.0088366652) q[2];
sx q[2];
rz(2.3542118) q[2];
rz(2.8619838) q[3];
sx q[3];
rz(-0.044737261) q[3];
sx q[3];
rz(-2.4217822) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1396007) q[0];
sx q[0];
rz(-0.75374341) q[0];
sx q[0];
rz(1.4308223) q[0];
rz(2.9665973) q[1];
sx q[1];
rz(-2.4038959) q[1];
sx q[1];
rz(-1.4281248) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5663877) q[0];
sx q[0];
rz(-3.1194127) q[0];
sx q[0];
rz(1.7087858) q[0];
rz(-pi) q[1];
rz(3.1375132) q[2];
sx q[2];
rz(-1.5634087) q[2];
sx q[2];
rz(1.0122062) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.1654351) q[1];
sx q[1];
rz(-1.4095029) q[1];
sx q[1];
rz(-1.2978202) q[1];
x q[2];
rz(-0.66760701) q[3];
sx q[3];
rz(-0.35304754) q[3];
sx q[3];
rz(-1.8065344) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.3565107) q[2];
sx q[2];
rz(-2.3729237) q[2];
sx q[2];
rz(1.53995) q[2];
rz(0.29940638) q[3];
sx q[3];
rz(-1.6573903) q[3];
sx q[3];
rz(0.33777344) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6701819) q[0];
sx q[0];
rz(-0.010936745) q[0];
sx q[0];
rz(-0.48728824) q[0];
rz(-1.7022853) q[1];
sx q[1];
rz(-0.7936365) q[1];
sx q[1];
rz(0.38584858) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13284616) q[0];
sx q[0];
rz(-1.1782968) q[0];
sx q[0];
rz(1.5061492) q[0];
rz(-pi) q[1];
rz(-0.11211498) q[2];
sx q[2];
rz(-1.9939853) q[2];
sx q[2];
rz(1.8433169) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.5380521) q[1];
sx q[1];
rz(-1.5376989) q[1];
sx q[1];
rz(-0.69149986) q[1];
rz(-pi) q[2];
rz(3.039592) q[3];
sx q[3];
rz(-0.33683646) q[3];
sx q[3];
rz(-0.61397314) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.7172598) q[2];
sx q[2];
rz(-3.11185) q[2];
sx q[2];
rz(-0.40561238) q[2];
rz(-0.82617104) q[3];
sx q[3];
rz(-0.043353733) q[3];
sx q[3];
rz(0.9894754) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9721603) q[0];
sx q[0];
rz(-0.5652453) q[0];
sx q[0];
rz(1.3954847) q[0];
rz(1.6204429) q[1];
sx q[1];
rz(-0.65137678) q[1];
sx q[1];
rz(1.7924538) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1698818) q[0];
sx q[0];
rz(-0.75877178) q[0];
sx q[0];
rz(1.2618598) q[0];
rz(1.8113891) q[2];
sx q[2];
rz(-1.7363915) q[2];
sx q[2];
rz(-0.22364932) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.3685128) q[1];
sx q[1];
rz(-3.0798612) q[1];
sx q[1];
rz(1.3365937) q[1];
rz(-pi) q[2];
rz(-1.5557846) q[3];
sx q[3];
rz(-1.6009357) q[3];
sx q[3];
rz(0.82885259) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.8162083) q[2];
sx q[2];
rz(-3.1367229) q[2];
sx q[2];
rz(1.4332786) q[2];
rz(-1.8520744) q[3];
sx q[3];
rz(-3.0472445) q[3];
sx q[3];
rz(-1.2963699) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0636487) q[0];
sx q[0];
rz(-1.0366806) q[0];
sx q[0];
rz(0.59564577) q[0];
rz(3.0972277) q[1];
sx q[1];
rz(-2.8652419) q[1];
sx q[1];
rz(-1.2038632) q[1];
rz(2.0478838) q[2];
sx q[2];
rz(-1.3774584) q[2];
sx q[2];
rz(1.8130345) q[2];
rz(1.1607004) q[3];
sx q[3];
rz(-1.6312508) q[3];
sx q[3];
rz(-1.4283258) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
