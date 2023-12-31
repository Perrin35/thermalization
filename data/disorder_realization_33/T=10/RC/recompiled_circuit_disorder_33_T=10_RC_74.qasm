OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.6361976) q[0];
sx q[0];
rz(-0.27591053) q[0];
sx q[0];
rz(1.3077868) q[0];
rz(1.1360599) q[1];
sx q[1];
rz(-0.93568957) q[1];
sx q[1];
rz(1.5703262) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6558134) q[0];
sx q[0];
rz(-1.2278779) q[0];
sx q[0];
rz(-1.9302084) q[0];
rz(-pi) q[1];
rz(0.70648944) q[2];
sx q[2];
rz(-0.90105614) q[2];
sx q[2];
rz(2.0073839) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.624144) q[1];
sx q[1];
rz(-1.2341208) q[1];
sx q[1];
rz(1.8271853) q[1];
rz(-1.0899815) q[3];
sx q[3];
rz(-0.40502031) q[3];
sx q[3];
rz(-0.68457505) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.2661665) q[2];
sx q[2];
rz(-2.8484919) q[2];
sx q[2];
rz(-2.0092633) q[2];
rz(1.4663565) q[3];
sx q[3];
rz(-1.8050067) q[3];
sx q[3];
rz(1.0124538) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
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
rz(-2.9448626) q[0];
sx q[0];
rz(-2.9319627) q[0];
sx q[0];
rz(0.18584132) q[0];
rz(-2.5813685) q[1];
sx q[1];
rz(-1.2954243) q[1];
sx q[1];
rz(2.9247608) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6235979) q[0];
sx q[0];
rz(-1.2835842) q[0];
sx q[0];
rz(0.90201305) q[0];
x q[1];
rz(-2.2611513) q[2];
sx q[2];
rz(-2.464622) q[2];
sx q[2];
rz(-2.0073236) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.9358881) q[1];
sx q[1];
rz(-1.7824031) q[1];
sx q[1];
rz(2.2536224) q[1];
rz(-pi) q[2];
rz(2.300755) q[3];
sx q[3];
rz(-0.84078046) q[3];
sx q[3];
rz(-0.18883146) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.310114) q[2];
sx q[2];
rz(-0.82565132) q[2];
sx q[2];
rz(1.2878093) q[2];
rz(-0.76256049) q[3];
sx q[3];
rz(-1.972714) q[3];
sx q[3];
rz(0.30502239) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6771616) q[0];
sx q[0];
rz(-2.7966249) q[0];
sx q[0];
rz(2.537354) q[0];
rz(-1.3263946) q[1];
sx q[1];
rz(-1.3605958) q[1];
sx q[1];
rz(2.2089829) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5843825) q[0];
sx q[0];
rz(-2.819448) q[0];
sx q[0];
rz(1.4003217) q[0];
rz(-0.18205299) q[2];
sx q[2];
rz(-1.7943873) q[2];
sx q[2];
rz(-1.6857266) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.4438666) q[1];
sx q[1];
rz(-1.6892471) q[1];
sx q[1];
rz(-2.5812134) q[1];
x q[2];
rz(-0.44585769) q[3];
sx q[3];
rz(-1.0599531) q[3];
sx q[3];
rz(1.0872935) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.147826) q[2];
sx q[2];
rz(-1.0819165) q[2];
sx q[2];
rz(-2.0489342) q[2];
rz(-0.5422194) q[3];
sx q[3];
rz(-1.0850302) q[3];
sx q[3];
rz(-0.96737635) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3595235) q[0];
sx q[0];
rz(-0.096465915) q[0];
sx q[0];
rz(0.50022593) q[0];
rz(-0.80530986) q[1];
sx q[1];
rz(-1.9814682) q[1];
sx q[1];
rz(-1.4979699) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26101199) q[0];
sx q[0];
rz(-1.8646761) q[0];
sx q[0];
rz(2.0902993) q[0];
rz(1.8565606) q[2];
sx q[2];
rz(-2.9432202) q[2];
sx q[2];
rz(2.6464268) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.30724635) q[1];
sx q[1];
rz(-1.8028959) q[1];
sx q[1];
rz(2.8744065) q[1];
rz(-pi) q[2];
rz(-2.6763776) q[3];
sx q[3];
rz(-1.9808931) q[3];
sx q[3];
rz(1.0614392) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.74636373) q[2];
sx q[2];
rz(-0.56240288) q[2];
sx q[2];
rz(2.4397819) q[2];
rz(0.83135215) q[3];
sx q[3];
rz(-2.1777007) q[3];
sx q[3];
rz(-0.62197661) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9005301) q[0];
sx q[0];
rz(-0.59589544) q[0];
sx q[0];
rz(-2.3262614) q[0];
rz(-1.5218081) q[1];
sx q[1];
rz(-2.3074469) q[1];
sx q[1];
rz(-2.0933847) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5532903) q[0];
sx q[0];
rz(-2.7748845) q[0];
sx q[0];
rz(-0.81272965) q[0];
rz(-pi) q[1];
rz(-1.5830718) q[2];
sx q[2];
rz(-0.91931146) q[2];
sx q[2];
rz(0.94142454) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.8909) q[1];
sx q[1];
rz(-0.39784583) q[1];
sx q[1];
rz(-2.489151) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7311677) q[3];
sx q[3];
rz(-0.99567185) q[3];
sx q[3];
rz(-0.23469532) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.6158225) q[2];
sx q[2];
rz(-0.56695357) q[2];
sx q[2];
rz(-2.0416416) q[2];
rz(-0.82529092) q[3];
sx q[3];
rz(-1.0422948) q[3];
sx q[3];
rz(2.2560789) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0734171) q[0];
sx q[0];
rz(-0.59403479) q[0];
sx q[0];
rz(-2.2391879) q[0];
rz(-2.1249318) q[1];
sx q[1];
rz(-2.0817751) q[1];
sx q[1];
rz(0.12983233) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8311365) q[0];
sx q[0];
rz(-1.2015011) q[0];
sx q[0];
rz(-0.62311689) q[0];
x q[1];
rz(2.1075222) q[2];
sx q[2];
rz(-2.495129) q[2];
sx q[2];
rz(-1.5922286) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-3.0778724) q[1];
sx q[1];
rz(-1.93396) q[1];
sx q[1];
rz(0.6086463) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5927605) q[3];
sx q[3];
rz(-1.4026814) q[3];
sx q[3];
rz(-0.43743922) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.8292024) q[2];
sx q[2];
rz(-0.94909334) q[2];
sx q[2];
rz(2.9373346) q[2];
rz(1.2060818) q[3];
sx q[3];
rz(-1.5217425) q[3];
sx q[3];
rz(2.9061785) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4181353) q[0];
sx q[0];
rz(-1.3292987) q[0];
sx q[0];
rz(1.4468505) q[0];
rz(-1.8824668) q[1];
sx q[1];
rz(-0.99021688) q[1];
sx q[1];
rz(-2.4553305) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9335564) q[0];
sx q[0];
rz(-0.66970034) q[0];
sx q[0];
rz(-2.3963388) q[0];
rz(-pi) q[1];
rz(-2.1967728) q[2];
sx q[2];
rz(-0.84569028) q[2];
sx q[2];
rz(-0.041989728) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.2346748) q[1];
sx q[1];
rz(-1.7559116) q[1];
sx q[1];
rz(-3.0659552) q[1];
rz(-pi) q[2];
rz(1.8920184) q[3];
sx q[3];
rz(-1.9739082) q[3];
sx q[3];
rz(-1.5461127) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.4454322) q[2];
sx q[2];
rz(-1.3779209) q[2];
sx q[2];
rz(-0.0017722842) q[2];
rz(2.5799675) q[3];
sx q[3];
rz(-2.2300945) q[3];
sx q[3];
rz(-1.5047489) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6034265) q[0];
sx q[0];
rz(-2.4551233) q[0];
sx q[0];
rz(1.4461393) q[0];
rz(2.360545) q[1];
sx q[1];
rz(-1.8361517) q[1];
sx q[1];
rz(-1.6400281) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38599309) q[0];
sx q[0];
rz(-2.0634683) q[0];
sx q[0];
rz(-0.022105769) q[0];
x q[1];
rz(-1.5090452) q[2];
sx q[2];
rz(-1.5442863) q[2];
sx q[2];
rz(0.74552958) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.6540263) q[1];
sx q[1];
rz(-1.8974202) q[1];
sx q[1];
rz(1.4766272) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.6286962) q[3];
sx q[3];
rz(-2.0572212) q[3];
sx q[3];
rz(1.6823671) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.35187307) q[2];
sx q[2];
rz(-1.7799653) q[2];
sx q[2];
rz(-1.3191351) q[2];
rz(1.2119279) q[3];
sx q[3];
rz(-1.2865678) q[3];
sx q[3];
rz(0.31931988) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.33655745) q[0];
sx q[0];
rz(-2.5890077) q[0];
sx q[0];
rz(1.9375027) q[0];
rz(0.38326344) q[1];
sx q[1];
rz(-0.52572322) q[1];
sx q[1];
rz(0.35167545) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.29031819) q[0];
sx q[0];
rz(-1.1997249) q[0];
sx q[0];
rz(-2.0593658) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8750538) q[2];
sx q[2];
rz(-2.2396302) q[2];
sx q[2];
rz(-0.63389102) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.5896776) q[1];
sx q[1];
rz(-2.3304686) q[1];
sx q[1];
rz(-3.0685436) q[1];
x q[2];
rz(-1.6230691) q[3];
sx q[3];
rz(-1.6802603) q[3];
sx q[3];
rz(2.0930406) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.3433156) q[2];
sx q[2];
rz(-1.1077935) q[2];
sx q[2];
rz(1.8593672) q[2];
rz(1.4964237) q[3];
sx q[3];
rz(-1.6069501) q[3];
sx q[3];
rz(1.055868) q[3];
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
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4984109) q[0];
sx q[0];
rz(-1.2675985) q[0];
sx q[0];
rz(2.9472651) q[0];
rz(2.1037897) q[1];
sx q[1];
rz(-2.5732645) q[1];
sx q[1];
rz(2.1077572) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72458306) q[0];
sx q[0];
rz(-1.4405182) q[0];
sx q[0];
rz(0.91086046) q[0];
x q[1];
rz(-0.98722234) q[2];
sx q[2];
rz(-2.3505031) q[2];
sx q[2];
rz(-0.9466048) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.49627134) q[1];
sx q[1];
rz(-1.6722582) q[1];
sx q[1];
rz(-2.1437777) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8909573) q[3];
sx q[3];
rz(-0.72892979) q[3];
sx q[3];
rz(0.38754101) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.0795435) q[2];
sx q[2];
rz(-2.1958308) q[2];
sx q[2];
rz(-0.6357843) q[2];
rz(-2.87129) q[3];
sx q[3];
rz(-0.79939866) q[3];
sx q[3];
rz(-1.5283782) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4476267) q[0];
sx q[0];
rz(-1.8287369) q[0];
sx q[0];
rz(1.0736314) q[0];
rz(1.4355961) q[1];
sx q[1];
rz(-1.5789079) q[1];
sx q[1];
rz(0.78067738) q[1];
rz(2.1736017) q[2];
sx q[2];
rz(-1.597076) q[2];
sx q[2];
rz(1.1521641) q[2];
rz(0.070449645) q[3];
sx q[3];
rz(-1.9000713) q[3];
sx q[3];
rz(-2.6303359) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
