OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.77914971) q[0];
sx q[0];
rz(-0.40894142) q[0];
sx q[0];
rz(-0.9932819) q[0];
rz(-0.14708695) q[1];
sx q[1];
rz(-0.99647254) q[1];
sx q[1];
rz(1.4176523) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.111747) q[0];
sx q[0];
rz(-0.8518712) q[0];
sx q[0];
rz(-2.0540038) q[0];
x q[1];
rz(0.63739802) q[2];
sx q[2];
rz(-2.3294724) q[2];
sx q[2];
rz(0.82773436) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.1163535) q[1];
sx q[1];
rz(-0.95121801) q[1];
sx q[1];
rz(-1.2337633) q[1];
rz(-pi) q[2];
rz(1.3325813) q[3];
sx q[3];
rz(-3.0034667) q[3];
sx q[3];
rz(-2.3254243) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.19109569) q[2];
sx q[2];
rz(-1.3206626) q[2];
sx q[2];
rz(1.0547868) q[2];
rz(0.47109207) q[3];
sx q[3];
rz(-1.4548917) q[3];
sx q[3];
rz(0.86827046) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7489814) q[0];
sx q[0];
rz(-1.6845717) q[0];
sx q[0];
rz(2.9066322) q[0];
rz(1.8670392) q[1];
sx q[1];
rz(-1.8320558) q[1];
sx q[1];
rz(2.4539006) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0318555) q[0];
sx q[0];
rz(-2.8068455) q[0];
sx q[0];
rz(-2.9358923) q[0];
x q[1];
rz(0.48414064) q[2];
sx q[2];
rz(-2.2233367) q[2];
sx q[2];
rz(2.2732041) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.5968273) q[1];
sx q[1];
rz(-3.1278962) q[1];
sx q[1];
rz(-1.9910452) q[1];
rz(-pi) q[2];
rz(0.092472381) q[3];
sx q[3];
rz(-2.1824129) q[3];
sx q[3];
rz(-0.05449748) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.71622744) q[2];
sx q[2];
rz(-1.8391515) q[2];
sx q[2];
rz(0.20935527) q[2];
rz(2.1685205) q[3];
sx q[3];
rz(-0.89997411) q[3];
sx q[3];
rz(-2.5488241) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8415602) q[0];
sx q[0];
rz(-0.40599269) q[0];
sx q[0];
rz(5/(11*pi)) q[0];
rz(-1.9035829) q[1];
sx q[1];
rz(-2.369945) q[1];
sx q[1];
rz(3.0498116) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.84052575) q[0];
sx q[0];
rz(-0.65945461) q[0];
sx q[0];
rz(-2.1546138) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1917159) q[2];
sx q[2];
rz(-0.11813049) q[2];
sx q[2];
rz(-2.0565513) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.8615177) q[1];
sx q[1];
rz(-0.28805486) q[1];
sx q[1];
rz(-0.83863284) q[1];
rz(-pi) q[2];
rz(-0.067236891) q[3];
sx q[3];
rz(-1.1394258) q[3];
sx q[3];
rz(2.9500285) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.7564275) q[2];
sx q[2];
rz(-1.6331208) q[2];
sx q[2];
rz(2.7623994) q[2];
rz(0.82342974) q[3];
sx q[3];
rz(-0.40612602) q[3];
sx q[3];
rz(-2.8520975) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32651153) q[0];
sx q[0];
rz(-0.87711763) q[0];
sx q[0];
rz(-0.99863482) q[0];
rz(0.48209349) q[1];
sx q[1];
rz(-0.95078743) q[1];
sx q[1];
rz(-1.3179717) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40967746) q[0];
sx q[0];
rz(-1.240726) q[0];
sx q[0];
rz(-2.5797352) q[0];
rz(2.5561294) q[2];
sx q[2];
rz(-1.3635067) q[2];
sx q[2];
rz(0.49400615) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.4596583) q[1];
sx q[1];
rz(-1.0409969) q[1];
sx q[1];
rz(1.3270686) q[1];
x q[2];
rz(0.79341268) q[3];
sx q[3];
rz(-1.9559091) q[3];
sx q[3];
rz(-2.957893) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.92129293) q[2];
sx q[2];
rz(-2.2396294) q[2];
sx q[2];
rz(0.48545066) q[2];
rz(2.7576533) q[3];
sx q[3];
rz(-1.5555614) q[3];
sx q[3];
rz(-2.3575822) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.081472814) q[0];
sx q[0];
rz(-3.0373242) q[0];
sx q[0];
rz(-3.1299348) q[0];
rz(-0.13721379) q[1];
sx q[1];
rz(-1.5882746) q[1];
sx q[1];
rz(-1.3020017) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.046776) q[0];
sx q[0];
rz(-1.4621549) q[0];
sx q[0];
rz(-2.957445) q[0];
rz(-pi) q[1];
rz(0.88413357) q[2];
sx q[2];
rz(-0.72956402) q[2];
sx q[2];
rz(-1.3183678) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.72921371) q[1];
sx q[1];
rz(-1.4625616) q[1];
sx q[1];
rz(1.257105) q[1];
rz(1.9309031) q[3];
sx q[3];
rz(-1.2146287) q[3];
sx q[3];
rz(2.3115013) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.4181218) q[2];
sx q[2];
rz(-1.8401517) q[2];
sx q[2];
rz(-2.8483025) q[2];
rz(-0.063118525) q[3];
sx q[3];
rz(-1.9189574) q[3];
sx q[3];
rz(-2.2466834) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0224595) q[0];
sx q[0];
rz(-0.28894153) q[0];
sx q[0];
rz(0.038507842) q[0];
rz(-1.0844082) q[1];
sx q[1];
rz(-0.42250982) q[1];
sx q[1];
rz(-2.1748621) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.777323) q[0];
sx q[0];
rz(-1.7614953) q[0];
sx q[0];
rz(0.1994812) q[0];
x q[1];
rz(-3.128563) q[2];
sx q[2];
rz(-1.3037062) q[2];
sx q[2];
rz(-0.14946562) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.8340862) q[1];
sx q[1];
rz(-0.53552578) q[1];
sx q[1];
rz(2.2811625) q[1];
x q[2];
rz(2.2214063) q[3];
sx q[3];
rz(-1.3704066) q[3];
sx q[3];
rz(-1.5707627) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.66095573) q[2];
sx q[2];
rz(-0.2350685) q[2];
sx q[2];
rz(0.98141518) q[2];
rz(1.4977945) q[3];
sx q[3];
rz(-1.8794329) q[3];
sx q[3];
rz(1.5952544) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3129231) q[0];
sx q[0];
rz(-0.68083119) q[0];
sx q[0];
rz(0.31461) q[0];
rz(2.9529052) q[1];
sx q[1];
rz(-0.43470946) q[1];
sx q[1];
rz(-0.46868086) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9148736) q[0];
sx q[0];
rz(-1.813199) q[0];
sx q[0];
rz(1.4366494) q[0];
rz(-pi) q[1];
rz(1.5596703) q[2];
sx q[2];
rz(-2.9738148) q[2];
sx q[2];
rz(-2.0027225) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.6275121) q[1];
sx q[1];
rz(-0.8610763) q[1];
sx q[1];
rz(-2.7372317) q[1];
x q[2];
rz(-2.5545337) q[3];
sx q[3];
rz(-0.54577845) q[3];
sx q[3];
rz(3.0445638) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.9357052) q[2];
sx q[2];
rz(-0.84172717) q[2];
sx q[2];
rz(2.1708798) q[2];
rz(0.5361706) q[3];
sx q[3];
rz(-1.2090809) q[3];
sx q[3];
rz(-0.71603388) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69990528) q[0];
sx q[0];
rz(-2.5869885) q[0];
sx q[0];
rz(-1.4458789) q[0];
rz(1.0384167) q[1];
sx q[1];
rz(-2.0380135) q[1];
sx q[1];
rz(3.0605002) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4202284) q[0];
sx q[0];
rz(-2.2330771) q[0];
sx q[0];
rz(1.9382947) q[0];
rz(-pi) q[1];
rz(-2.7620188) q[2];
sx q[2];
rz(-0.69131572) q[2];
sx q[2];
rz(2.8065681) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.1782185) q[1];
sx q[1];
rz(-1.6859173) q[1];
sx q[1];
rz(1.8646152) q[1];
rz(-pi) q[2];
x q[2];
rz(0.57402667) q[3];
sx q[3];
rz(-1.9732631) q[3];
sx q[3];
rz(2.5695679) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.1991835) q[2];
sx q[2];
rz(-0.77745357) q[2];
sx q[2];
rz(2.8988163) q[2];
rz(1.5593922) q[3];
sx q[3];
rz(-2.715761) q[3];
sx q[3];
rz(1.2529713) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2389857) q[0];
sx q[0];
rz(-2.1098397) q[0];
sx q[0];
rz(1.2549988) q[0];
rz(1.3764508) q[1];
sx q[1];
rz(-2.1170728) q[1];
sx q[1];
rz(-0.66744101) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8534626) q[0];
sx q[0];
rz(-1.2391744) q[0];
sx q[0];
rz(-2.1037071) q[0];
x q[1];
rz(2.6450577) q[2];
sx q[2];
rz(-1.5333954) q[2];
sx q[2];
rz(2.2384245) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.63782843) q[1];
sx q[1];
rz(-1.003028) q[1];
sx q[1];
rz(2.7953887) q[1];
rz(-0.38102229) q[3];
sx q[3];
rz(-1.2667315) q[3];
sx q[3];
rz(-1.1883433) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.61975512) q[2];
sx q[2];
rz(-0.73885584) q[2];
sx q[2];
rz(0.45836207) q[2];
rz(1.6058263) q[3];
sx q[3];
rz(-2.063844) q[3];
sx q[3];
rz(-2.2569236) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.89852029) q[0];
sx q[0];
rz(-2.5694818) q[0];
sx q[0];
rz(-2.7761053) q[0];
rz(0.99994031) q[1];
sx q[1];
rz(-1.4391856) q[1];
sx q[1];
rz(1.684729) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4936876) q[0];
sx q[0];
rz(-2.5468605) q[0];
sx q[0];
rz(-1.2059284) q[0];
rz(-1.1816977) q[2];
sx q[2];
rz(-0.28841296) q[2];
sx q[2];
rz(-2.4450977) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.3927276) q[1];
sx q[1];
rz(-0.76887776) q[1];
sx q[1];
rz(2.8906872) q[1];
rz(-pi) q[2];
rz(2.445128) q[3];
sx q[3];
rz(-1.9380155) q[3];
sx q[3];
rz(-2.5270278) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.140787) q[2];
sx q[2];
rz(-1.8203338) q[2];
sx q[2];
rz(2.136039) q[2];
rz(2.4908861) q[3];
sx q[3];
rz(-1.7020099) q[3];
sx q[3];
rz(0.85056359) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0477796) q[0];
sx q[0];
rz(-2.35738) q[0];
sx q[0];
rz(2.1398075) q[0];
rz(1.7306937) q[1];
sx q[1];
rz(-1.7041364) q[1];
sx q[1];
rz(-2.0028353) q[1];
rz(-1.2391113) q[2];
sx q[2];
rz(-0.3467691) q[2];
sx q[2];
rz(-2.3820387) q[2];
rz(-0.10436124) q[3];
sx q[3];
rz(-1.4781393) q[3];
sx q[3];
rz(0.27235336) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
