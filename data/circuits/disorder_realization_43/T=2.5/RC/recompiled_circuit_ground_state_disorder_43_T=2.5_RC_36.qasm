OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.7795774) q[0];
sx q[0];
rz(-0.22688046) q[0];
sx q[0];
rz(1.7664631) q[0];
rz(3.0472164) q[1];
sx q[1];
rz(-0.91369319) q[1];
sx q[1];
rz(0.28092608) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8178806) q[0];
sx q[0];
rz(-1.4120988) q[0];
sx q[0];
rz(-2.236332) q[0];
rz(1.2035349) q[2];
sx q[2];
rz(-0.31031552) q[2];
sx q[2];
rz(-0.65904891) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.0620131) q[1];
sx q[1];
rz(-1.8342736) q[1];
sx q[1];
rz(1.1398214) q[1];
x q[2];
rz(-0.57262086) q[3];
sx q[3];
rz(-1.2702281) q[3];
sx q[3];
rz(2.4550748) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.0278339) q[2];
sx q[2];
rz(-1.7366624) q[2];
sx q[2];
rz(0.11715451) q[2];
rz(0.33200085) q[3];
sx q[3];
rz(-0.74688512) q[3];
sx q[3];
rz(-0.78506708) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4431045) q[0];
sx q[0];
rz(-0.76258689) q[0];
sx q[0];
rz(-0.37345988) q[0];
rz(2.7442878) q[1];
sx q[1];
rz(-1.1445069) q[1];
sx q[1];
rz(0.73748803) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0126368) q[0];
sx q[0];
rz(-1.0397433) q[0];
sx q[0];
rz(0.80967243) q[0];
x q[1];
rz(0.55464427) q[2];
sx q[2];
rz(-2.0337542) q[2];
sx q[2];
rz(0.92051586) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.0651689) q[1];
sx q[1];
rz(-2.5701414) q[1];
sx q[1];
rz(2.854611) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2190388) q[3];
sx q[3];
rz(-2.2969807) q[3];
sx q[3];
rz(1.5016457) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.071659) q[2];
sx q[2];
rz(-1.5295014) q[2];
sx q[2];
rz(2.7640479) q[2];
rz(0.30970445) q[3];
sx q[3];
rz(-2.0524502) q[3];
sx q[3];
rz(-1.496605) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6250703) q[0];
sx q[0];
rz(-2.3023038) q[0];
sx q[0];
rz(-1.2155493) q[0];
rz(-1.0141605) q[1];
sx q[1];
rz(-1.4233669) q[1];
sx q[1];
rz(-2.1713712) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1078915) q[0];
sx q[0];
rz(-1.530307) q[0];
sx q[0];
rz(2.070436) q[0];
rz(-pi) q[1];
rz(-0.90486352) q[2];
sx q[2];
rz(-2.3653125) q[2];
sx q[2];
rz(1.5011476) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.6128224) q[1];
sx q[1];
rz(-0.38685683) q[1];
sx q[1];
rz(-1.0861138) q[1];
rz(-pi) q[2];
x q[2];
rz(0.35200624) q[3];
sx q[3];
rz(-1.4624634) q[3];
sx q[3];
rz(-1.4198522) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.0936475) q[2];
sx q[2];
rz(-2.2990172) q[2];
sx q[2];
rz(0.45316163) q[2];
rz(-1.2293182) q[3];
sx q[3];
rz(-1.8639576) q[3];
sx q[3];
rz(2.9696828) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1969084) q[0];
sx q[0];
rz(-2.4561645) q[0];
sx q[0];
rz(0.36566439) q[0];
rz(-2.2293495) q[1];
sx q[1];
rz(-1.8625926) q[1];
sx q[1];
rz(2.7361187) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1039818) q[0];
sx q[0];
rz(-1.77586) q[0];
sx q[0];
rz(0.12139856) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.916476) q[2];
sx q[2];
rz(-0.76864132) q[2];
sx q[2];
rz(-2.3942238) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.7995389) q[1];
sx q[1];
rz(-2.5785682) q[1];
sx q[1];
rz(-2.8232226) q[1];
rz(-0.63946569) q[3];
sx q[3];
rz(-1.3124976) q[3];
sx q[3];
rz(2.4909508) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.071986467) q[2];
sx q[2];
rz(-1.4402086) q[2];
sx q[2];
rz(-1.5988505) q[2];
rz(-2.767848) q[3];
sx q[3];
rz(-1.475324) q[3];
sx q[3];
rz(-1.4603978) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54295802) q[0];
sx q[0];
rz(-2.3668508) q[0];
sx q[0];
rz(0.69611088) q[0];
rz(-1.8822582) q[1];
sx q[1];
rz(-1.101661) q[1];
sx q[1];
rz(2.7764244) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8697978) q[0];
sx q[0];
rz(-2.4569017) q[0];
sx q[0];
rz(-2.4224046) q[0];
x q[1];
rz(2.8234473) q[2];
sx q[2];
rz(-1.1765624) q[2];
sx q[2];
rz(2.1925558) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.9112253) q[1];
sx q[1];
rz(-2.5960698) q[1];
sx q[1];
rz(0.0026774252) q[1];
rz(-pi) q[2];
rz(0.62778421) q[3];
sx q[3];
rz(-1.9115698) q[3];
sx q[3];
rz(-2.9140811) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.68088561) q[2];
sx q[2];
rz(-1.978771) q[2];
sx q[2];
rz(2.9648901) q[2];
rz(1.7548615) q[3];
sx q[3];
rz(-1.0597798) q[3];
sx q[3];
rz(0.37469125) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
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
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.52727592) q[0];
sx q[0];
rz(-0.85955954) q[0];
sx q[0];
rz(1.543462) q[0];
rz(1.8371001) q[1];
sx q[1];
rz(-1.0544798) q[1];
sx q[1];
rz(-2.3354882) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4325754) q[0];
sx q[0];
rz(-2.3516549) q[0];
sx q[0];
rz(-0.33613251) q[0];
rz(-0.51050425) q[2];
sx q[2];
rz(-2.1091828) q[2];
sx q[2];
rz(-1.5427854) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.42000586) q[1];
sx q[1];
rz(-1.3665821) q[1];
sx q[1];
rz(3.0908683) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0982355) q[3];
sx q[3];
rz(-1.0602078) q[3];
sx q[3];
rz(-3.1289738) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.1690037) q[2];
sx q[2];
rz(-0.4946332) q[2];
sx q[2];
rz(-3.0779823) q[2];
rz(-2.5566067) q[3];
sx q[3];
rz(-1.4002742) q[3];
sx q[3];
rz(1.5144279) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1503898) q[0];
sx q[0];
rz(-1.7504033) q[0];
sx q[0];
rz(-2.4124131) q[0];
rz(-1.179262) q[1];
sx q[1];
rz(-2.3919892) q[1];
sx q[1];
rz(2.8470305) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.23239947) q[0];
sx q[0];
rz(-2.5721484) q[0];
sx q[0];
rz(-2.6528986) q[0];
x q[1];
rz(1.137758) q[2];
sx q[2];
rz(-1.6707509) q[2];
sx q[2];
rz(-1.86509) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.2475491) q[1];
sx q[1];
rz(-1.6881384) q[1];
sx q[1];
rz(2.3168922) q[1];
rz(-2.758243) q[3];
sx q[3];
rz(-2.2046996) q[3];
sx q[3];
rz(-2.1103566) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.4906759) q[2];
sx q[2];
rz(-1.7016405) q[2];
sx q[2];
rz(2.4093464) q[2];
rz(0.8693153) q[3];
sx q[3];
rz(-1.9711875) q[3];
sx q[3];
rz(-1.4066345) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-1.231584) q[0];
sx q[0];
rz(-1.5861479) q[0];
sx q[0];
rz(-0.79767942) q[0];
rz(-0.32294598) q[1];
sx q[1];
rz(-0.99635092) q[1];
sx q[1];
rz(1.4083883) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.080264576) q[0];
sx q[0];
rz(-1.7751851) q[0];
sx q[0];
rz(0.56046446) q[0];
x q[1];
rz(2.9041821) q[2];
sx q[2];
rz(-2.0185197) q[2];
sx q[2];
rz(-2.6765649) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.3288142) q[1];
sx q[1];
rz(-1.8056889) q[1];
sx q[1];
rz(0.77130227) q[1];
x q[2];
rz(1.3878294) q[3];
sx q[3];
rz(-1.2550233) q[3];
sx q[3];
rz(-2.0980199) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.8352167) q[2];
sx q[2];
rz(-1.110346) q[2];
sx q[2];
rz(-2.9442673) q[2];
rz(-0.80162445) q[3];
sx q[3];
rz(-0.97951952) q[3];
sx q[3];
rz(-2.7200123) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1590969) q[0];
sx q[0];
rz(-1.6668586) q[0];
sx q[0];
rz(2.2267447) q[0];
rz(-0.21028701) q[1];
sx q[1];
rz(-0.57840127) q[1];
sx q[1];
rz(-2.4519144) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.294201) q[0];
sx q[0];
rz(-2.8970844) q[0];
sx q[0];
rz(-0.39347009) q[0];
rz(-pi) q[1];
rz(-1.44697) q[2];
sx q[2];
rz(-3.0231907) q[2];
sx q[2];
rz(-1.3294544) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.0090078) q[1];
sx q[1];
rz(-1.0632391) q[1];
sx q[1];
rz(-2.1193488) q[1];
rz(-pi) q[2];
rz(-0.34387572) q[3];
sx q[3];
rz(-1.0707885) q[3];
sx q[3];
rz(-0.22658928) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.0476394) q[2];
sx q[2];
rz(-0.84838212) q[2];
sx q[2];
rz(0.32996714) q[2];
rz(-1.0432358) q[3];
sx q[3];
rz(-1.7236575) q[3];
sx q[3];
rz(0.51658336) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1107776) q[0];
sx q[0];
rz(-1.3309706) q[0];
sx q[0];
rz(-2.0844841) q[0];
rz(2.8874176) q[1];
sx q[1];
rz(-1.9994241) q[1];
sx q[1];
rz(2.4370297) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0443665) q[0];
sx q[0];
rz(-1.2597879) q[0];
sx q[0];
rz(-1.0754271) q[0];
rz(-pi) q[1];
rz(-1.754934) q[2];
sx q[2];
rz(-1.6312851) q[2];
sx q[2];
rz(-2.9764701) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.89082644) q[1];
sx q[1];
rz(-1.1264404) q[1];
sx q[1];
rz(-0.19877427) q[1];
rz(-0.20829717) q[3];
sx q[3];
rz(-1.7476255) q[3];
sx q[3];
rz(3.1272776) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.5083984) q[2];
sx q[2];
rz(-1.9693547) q[2];
sx q[2];
rz(0.50676662) q[2];
rz(0.1861598) q[3];
sx q[3];
rz(-0.23701826) q[3];
sx q[3];
rz(-0.53562927) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-1.7645466) q[0];
sx q[0];
rz(-1.4519539) q[0];
sx q[0];
rz(-2.0013381) q[0];
rz(2.2346732) q[1];
sx q[1];
rz(-0.74105558) q[1];
sx q[1];
rz(1.8190609) q[1];
rz(2.604031) q[2];
sx q[2];
rz(-0.96502177) q[2];
sx q[2];
rz(1.6966664) q[2];
rz(2.4558057) q[3];
sx q[3];
rz(-2.7332173) q[3];
sx q[3];
rz(0.12242534) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
