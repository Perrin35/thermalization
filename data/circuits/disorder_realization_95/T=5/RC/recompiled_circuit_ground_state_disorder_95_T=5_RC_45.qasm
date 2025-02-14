OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.65741444) q[0];
sx q[0];
rz(-1.8622458) q[0];
sx q[0];
rz(2.4980463) q[0];
rz(-2.9930288) q[1];
sx q[1];
rz(-0.64639503) q[1];
sx q[1];
rz(0.57599154) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51645378) q[0];
sx q[0];
rz(-1.0363434) q[0];
sx q[0];
rz(2.778976) q[0];
rz(-pi) q[1];
rz(-1.2744489) q[2];
sx q[2];
rz(-1.7360032) q[2];
sx q[2];
rz(-2.8458386) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.6580087) q[1];
sx q[1];
rz(-1.600126) q[1];
sx q[1];
rz(1.6701677) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7891145) q[3];
sx q[3];
rz(-2.3872833) q[3];
sx q[3];
rz(-1.393817) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.6286455) q[2];
sx q[2];
rz(-0.60428667) q[2];
sx q[2];
rz(0.71913546) q[2];
rz(2.5288845) q[3];
sx q[3];
rz(-2.6477974) q[3];
sx q[3];
rz(1.4601716) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4145819) q[0];
sx q[0];
rz(-0.56986037) q[0];
sx q[0];
rz(1.3563096) q[0];
rz(-0.59056979) q[1];
sx q[1];
rz(-1.5547724) q[1];
sx q[1];
rz(0.86033598) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7588494) q[0];
sx q[0];
rz(-2.4573987) q[0];
sx q[0];
rz(-1.280715) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7458631) q[2];
sx q[2];
rz(-1.0222728) q[2];
sx q[2];
rz(-2.0344987) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.0257657) q[1];
sx q[1];
rz(-1.888926) q[1];
sx q[1];
rz(-0.33237388) q[1];
x q[2];
rz(3.0290696) q[3];
sx q[3];
rz(-0.72970245) q[3];
sx q[3];
rz(-2.2000748) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.7343973) q[2];
sx q[2];
rz(-2.037214) q[2];
sx q[2];
rz(0.26367903) q[2];
rz(-3.0401491) q[3];
sx q[3];
rz(-2.0439309) q[3];
sx q[3];
rz(1.4640079) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2270671) q[0];
sx q[0];
rz(-2.3093746) q[0];
sx q[0];
rz(1.9050003) q[0];
rz(0.49691686) q[1];
sx q[1];
rz(-2.2221815) q[1];
sx q[1];
rz(0.20981728) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1001756) q[0];
sx q[0];
rz(-0.13209535) q[0];
sx q[0];
rz(0.14762525) q[0];
rz(-pi) q[1];
rz(-3.0044697) q[2];
sx q[2];
rz(-1.1242766) q[2];
sx q[2];
rz(1.5369122) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.57269872) q[1];
sx q[1];
rz(-2.0742744) q[1];
sx q[1];
rz(3.0996432) q[1];
x q[2];
rz(-2.9852226) q[3];
sx q[3];
rz(-0.84330446) q[3];
sx q[3];
rz(2.7263977) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.0158374) q[2];
sx q[2];
rz(-2.0070984) q[2];
sx q[2];
rz(-2.8718359) q[2];
rz(-2.6848327) q[3];
sx q[3];
rz(-0.41958198) q[3];
sx q[3];
rz(2.8381798) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1064827) q[0];
sx q[0];
rz(-0.39316097) q[0];
sx q[0];
rz(-3.0495354) q[0];
rz(2.5606142) q[1];
sx q[1];
rz(-2.5216504) q[1];
sx q[1];
rz(0.43749014) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5212183) q[0];
sx q[0];
rz(-2.552444) q[0];
sx q[0];
rz(-1.8626088) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0552513) q[2];
sx q[2];
rz(-0.24366442) q[2];
sx q[2];
rz(2.551384) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.535634) q[1];
sx q[1];
rz(-1.8211963) q[1];
sx q[1];
rz(1.1634924) q[1];
rz(0.32467036) q[3];
sx q[3];
rz(-0.89327795) q[3];
sx q[3];
rz(0.087740104) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.1979394) q[2];
sx q[2];
rz(-1.7655756) q[2];
sx q[2];
rz(0.1680689) q[2];
rz(-2.4456612) q[3];
sx q[3];
rz(-1.9441425) q[3];
sx q[3];
rz(1.7456938) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.54774061) q[0];
sx q[0];
rz(-2.934179) q[0];
sx q[0];
rz(-0.70163027) q[0];
rz(-0.54948366) q[1];
sx q[1];
rz(-2.0618849) q[1];
sx q[1];
rz(-2.2061677) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0129378) q[0];
sx q[0];
rz(-1.8194981) q[0];
sx q[0];
rz(2.0123345) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0183724) q[2];
sx q[2];
rz(-1.5180523) q[2];
sx q[2];
rz(0.067600943) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.6657609) q[1];
sx q[1];
rz(-2.40959) q[1];
sx q[1];
rz(-0.93722645) q[1];
rz(-pi) q[2];
rz(-0.071604972) q[3];
sx q[3];
rz(-1.5061028) q[3];
sx q[3];
rz(-0.66401362) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.51297411) q[2];
sx q[2];
rz(-1.8206035) q[2];
sx q[2];
rz(-0.67250195) q[2];
rz(2.3682112) q[3];
sx q[3];
rz(-2.0956109) q[3];
sx q[3];
rz(2.8092303) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(0.32909876) q[0];
sx q[0];
rz(-0.90730539) q[0];
sx q[0];
rz(-1.4688274) q[0];
rz(2.2724197) q[1];
sx q[1];
rz(-0.69907993) q[1];
sx q[1];
rz(-2.8156143) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7237968) q[0];
sx q[0];
rz(-1.4362596) q[0];
sx q[0];
rz(0.30067443) q[0];
rz(-pi) q[1];
rz(2.5766854) q[2];
sx q[2];
rz(-2.0015026) q[2];
sx q[2];
rz(-2.4912094) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.80854844) q[1];
sx q[1];
rz(-0.99892975) q[1];
sx q[1];
rz(0.99355917) q[1];
x q[2];
rz(-3.0748653) q[3];
sx q[3];
rz(-1.1432093) q[3];
sx q[3];
rz(0.59837435) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.088923067) q[2];
sx q[2];
rz(-1.5648877) q[2];
sx q[2];
rz(-0.79724533) q[2];
rz(-0.097660216) q[3];
sx q[3];
rz(-0.74334136) q[3];
sx q[3];
rz(-1.9184448) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.866211) q[0];
sx q[0];
rz(-0.99358639) q[0];
sx q[0];
rz(0.98268253) q[0];
rz(-1.4879701) q[1];
sx q[1];
rz(-0.5653615) q[1];
sx q[1];
rz(0.32378325) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.57162962) q[0];
sx q[0];
rz(-1.4870207) q[0];
sx q[0];
rz(-0.28599583) q[0];
rz(-pi) q[1];
x q[1];
rz(0.95877842) q[2];
sx q[2];
rz(-1.806815) q[2];
sx q[2];
rz(0.70597202) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.1947136) q[1];
sx q[1];
rz(-2.4733798) q[1];
sx q[1];
rz(0.18385033) q[1];
rz(-pi) q[2];
rz(2.5681132) q[3];
sx q[3];
rz(-2.5635168) q[3];
sx q[3];
rz(0.40595192) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.4299778) q[2];
sx q[2];
rz(-0.51402503) q[2];
sx q[2];
rz(3.1058822) q[2];
rz(-2.4014373) q[3];
sx q[3];
rz(-1.454708) q[3];
sx q[3];
rz(-1.1470225) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3586054) q[0];
sx q[0];
rz(-0.37707034) q[0];
sx q[0];
rz(-0.20088917) q[0];
rz(-0.56328493) q[1];
sx q[1];
rz(-1.3689684) q[1];
sx q[1];
rz(2.7422781) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2173846) q[0];
sx q[0];
rz(-1.426412) q[0];
sx q[0];
rz(-1.7213137) q[0];
rz(-pi) q[1];
rz(-3.1004347) q[2];
sx q[2];
rz(-1.090547) q[2];
sx q[2];
rz(0.89370773) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.3914403) q[1];
sx q[1];
rz(-1.2160157) q[1];
sx q[1];
rz(-1.1417692) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3763555) q[3];
sx q[3];
rz(-1.2685734) q[3];
sx q[3];
rz(2.0207495) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.1839972) q[2];
sx q[2];
rz(-0.60773578) q[2];
sx q[2];
rz(0.11873928) q[2];
rz(0.27058288) q[3];
sx q[3];
rz(-2.1515473) q[3];
sx q[3];
rz(2.9873007) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1110693) q[0];
sx q[0];
rz(-1.1708165) q[0];
sx q[0];
rz(0.4554553) q[0];
rz(-2.9083374) q[1];
sx q[1];
rz(-0.74359727) q[1];
sx q[1];
rz(0.0049678405) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.57039047) q[0];
sx q[0];
rz(-1.5367147) q[0];
sx q[0];
rz(0.081460074) q[0];
rz(-pi) q[1];
rz(-0.67989477) q[2];
sx q[2];
rz(-0.48326421) q[2];
sx q[2];
rz(-2.6269185) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.3525654) q[1];
sx q[1];
rz(-2.0754756) q[1];
sx q[1];
rz(0.53225174) q[1];
rz(0.016214215) q[3];
sx q[3];
rz(-1.1450279) q[3];
sx q[3];
rz(-0.92416422) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.7838955) q[2];
sx q[2];
rz(-1.8507439) q[2];
sx q[2];
rz(2.2515187) q[2];
rz(-0.87738532) q[3];
sx q[3];
rz(-0.93534094) q[3];
sx q[3];
rz(-2.4963511) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.18168618) q[0];
sx q[0];
rz(-3.063995) q[0];
sx q[0];
rz(-1.0542057) q[0];
rz(-0.28610817) q[1];
sx q[1];
rz(-1.0461297) q[1];
sx q[1];
rz(1.2623164) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6968203) q[0];
sx q[0];
rz(-0.69545924) q[0];
sx q[0];
rz(-2.6878342) q[0];
rz(-pi) q[1];
rz(-1.3081864) q[2];
sx q[2];
rz(-1.3725026) q[2];
sx q[2];
rz(-1.8526371) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.4642548) q[1];
sx q[1];
rz(-2.6644775) q[1];
sx q[1];
rz(0.40661637) q[1];
x q[2];
rz(1.0159053) q[3];
sx q[3];
rz(-1.6421179) q[3];
sx q[3];
rz(-1.9147022) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.012764843) q[2];
sx q[2];
rz(-0.78170693) q[2];
sx q[2];
rz(-0.26415602) q[2];
rz(-0.24010298) q[3];
sx q[3];
rz(-2.0930347) q[3];
sx q[3];
rz(2.8255294) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24513292) q[0];
sx q[0];
rz(-2.3871853) q[0];
sx q[0];
rz(2.2375794) q[0];
rz(2.8527507) q[1];
sx q[1];
rz(-1.387351) q[1];
sx q[1];
rz(3.0659061) q[1];
rz(0.31927989) q[2];
sx q[2];
rz(-1.8055331) q[2];
sx q[2];
rz(-1.8100678) q[2];
rz(3.0700923) q[3];
sx q[3];
rz(-1.5300453) q[3];
sx q[3];
rz(1.9634043) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
