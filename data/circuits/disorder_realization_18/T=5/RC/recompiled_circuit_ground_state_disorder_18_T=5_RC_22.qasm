OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.1825778) q[0];
sx q[0];
rz(-0.20354095) q[0];
sx q[0];
rz(1.8939053) q[0];
rz(-1.7893451) q[1];
sx q[1];
rz(-2.4440553) q[1];
sx q[1];
rz(-2.4434659) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9331935) q[0];
sx q[0];
rz(-2.0740447) q[0];
sx q[0];
rz(-1.1017975) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.013555992) q[2];
sx q[2];
rz(-0.45326158) q[2];
sx q[2];
rz(2.1238985) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.9351077) q[1];
sx q[1];
rz(-0.84151959) q[1];
sx q[1];
rz(-0.40928276) q[1];
x q[2];
rz(1.8258445) q[3];
sx q[3];
rz(-1.5849893) q[3];
sx q[3];
rz(0.04081459) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.84833604) q[2];
sx q[2];
rz(-1.9170599) q[2];
sx q[2];
rz(-0.87948925) q[2];
rz(0.73421156) q[3];
sx q[3];
rz(-2.0101533) q[3];
sx q[3];
rz(2.3122299) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.64049983) q[0];
sx q[0];
rz(-1.1859897) q[0];
sx q[0];
rz(-2.1511141) q[0];
rz(0.99539202) q[1];
sx q[1];
rz(-1.4442911) q[1];
sx q[1];
rz(-0.46517864) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.746691) q[0];
sx q[0];
rz(-1.9865661) q[0];
sx q[0];
rz(-1.3711934) q[0];
x q[1];
rz(-2.1976333) q[2];
sx q[2];
rz(-1.2612245) q[2];
sx q[2];
rz(2.7322526) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.3639646) q[1];
sx q[1];
rz(-0.87715699) q[1];
sx q[1];
rz(-0.65366807) q[1];
rz(-1.3690794) q[3];
sx q[3];
rz(-1.8276102) q[3];
sx q[3];
rz(-1.6981924) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.62749949) q[2];
sx q[2];
rz(-1.2343312) q[2];
sx q[2];
rz(0.17847432) q[2];
rz(-0.56525362) q[3];
sx q[3];
rz(-1.7491128) q[3];
sx q[3];
rz(2.5843411) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31767118) q[0];
sx q[0];
rz(-1.3452106) q[0];
sx q[0];
rz(-2.4406261) q[0];
rz(-2.7340381) q[1];
sx q[1];
rz(-1.2113672) q[1];
sx q[1];
rz(-1.0754546) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.70301818) q[0];
sx q[0];
rz(-2.10963) q[0];
sx q[0];
rz(-1.1807022) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8604061) q[2];
sx q[2];
rz(-2.4649005) q[2];
sx q[2];
rz(1.8829568) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.13710626) q[1];
sx q[1];
rz(-1.9705233) q[1];
sx q[1];
rz(-0.14182472) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.386679) q[3];
sx q[3];
rz(-2.872481) q[3];
sx q[3];
rz(1.6412013) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.0740697) q[2];
sx q[2];
rz(-1.8068376) q[2];
sx q[2];
rz(-1.4074116) q[2];
rz(-0.47809005) q[3];
sx q[3];
rz(-2.7985088) q[3];
sx q[3];
rz(2.7966255) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.57426977) q[0];
sx q[0];
rz(-1.3936717) q[0];
sx q[0];
rz(-1.0465013) q[0];
rz(-2.8872755) q[1];
sx q[1];
rz(-2.3260702) q[1];
sx q[1];
rz(-2.8291124) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.69092551) q[0];
sx q[0];
rz(-1.1617799) q[0];
sx q[0];
rz(-0.14139195) q[0];
x q[1];
rz(-2.4813969) q[2];
sx q[2];
rz(-1.1822299) q[2];
sx q[2];
rz(1.6903433) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.99287886) q[1];
sx q[1];
rz(-1.8853469) q[1];
sx q[1];
rz(1.3158362) q[1];
rz(-pi) q[2];
rz(-1.0279827) q[3];
sx q[3];
rz(-2.5590754) q[3];
sx q[3];
rz(-1.5912671) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.3089402) q[2];
sx q[2];
rz(-0.71844429) q[2];
sx q[2];
rz(-2.9529875) q[2];
rz(0.23379937) q[3];
sx q[3];
rz(-2.2457687) q[3];
sx q[3];
rz(2.5578267) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4466062) q[0];
sx q[0];
rz(-1.3124895) q[0];
sx q[0];
rz(-3.1332698) q[0];
rz(-2.7024929) q[1];
sx q[1];
rz(-1.3689901) q[1];
sx q[1];
rz(1.8956553) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9326235) q[0];
sx q[0];
rz(-2.0216938) q[0];
sx q[0];
rz(1.3572925) q[0];
rz(1.0143004) q[2];
sx q[2];
rz(-0.27703962) q[2];
sx q[2];
rz(3.1173681) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.51366561) q[1];
sx q[1];
rz(-0.49066077) q[1];
sx q[1];
rz(2.2918244) q[1];
x q[2];
rz(-1.9381136) q[3];
sx q[3];
rz(-2.4532336) q[3];
sx q[3];
rz(1.7662774) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.40450725) q[2];
sx q[2];
rz(-3.0948907) q[2];
sx q[2];
rz(1.5076293) q[2];
rz(-2.5493933) q[3];
sx q[3];
rz(-1.3785572) q[3];
sx q[3];
rz(-0.73970214) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.020788766) q[0];
sx q[0];
rz(-1.4494267) q[0];
sx q[0];
rz(2.8756496) q[0];
rz(-1.2159011) q[1];
sx q[1];
rz(-2.3561056) q[1];
sx q[1];
rz(-1.9507834) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.94534369) q[0];
sx q[0];
rz(-2.4085345) q[0];
sx q[0];
rz(-1.5668012) q[0];
rz(-pi) q[1];
rz(3.1053084) q[2];
sx q[2];
rz(-2.1232833) q[2];
sx q[2];
rz(-0.82922574) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.3490847) q[1];
sx q[1];
rz(-0.81696063) q[1];
sx q[1];
rz(1.5813827) q[1];
rz(-pi) q[2];
rz(-0.27482185) q[3];
sx q[3];
rz(-1.6259818) q[3];
sx q[3];
rz(-1.382687) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.8094981) q[2];
sx q[2];
rz(-1.7569434) q[2];
sx q[2];
rz(0.14796999) q[2];
rz(-2.68908) q[3];
sx q[3];
rz(-2.0544402) q[3];
sx q[3];
rz(-0.40542671) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0026523503) q[0];
sx q[0];
rz(-1.2661221) q[0];
sx q[0];
rz(-1.7953605) q[0];
rz(-0.0070618709) q[1];
sx q[1];
rz(-2.4738753) q[1];
sx q[1];
rz(2.6341569) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.19120309) q[0];
sx q[0];
rz(-1.0325214) q[0];
sx q[0];
rz(0.61417587) q[0];
rz(-pi) q[1];
rz(-0.15501898) q[2];
sx q[2];
rz(-1.2298202) q[2];
sx q[2];
rz(2.2237958) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.7120613) q[1];
sx q[1];
rz(-1.1768394) q[1];
sx q[1];
rz(1.560657) q[1];
rz(1.8431013) q[3];
sx q[3];
rz(-1.57988) q[3];
sx q[3];
rz(-2.7177725) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.4355882) q[2];
sx q[2];
rz(-2.0861349) q[2];
sx q[2];
rz(2.9662507) q[2];
rz(-0.38241479) q[3];
sx q[3];
rz(-1.9491111) q[3];
sx q[3];
rz(0.46149883) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14313702) q[0];
sx q[0];
rz(-2.1303506) q[0];
sx q[0];
rz(-1.829041) q[0];
rz(-0.10874272) q[1];
sx q[1];
rz(-2.5332632) q[1];
sx q[1];
rz(2.463602) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0339542) q[0];
sx q[0];
rz(-1.0773412) q[0];
sx q[0];
rz(-3.1163868) q[0];
rz(-pi) q[1];
x q[1];
rz(0.021804734) q[2];
sx q[2];
rz(-2.3924689) q[2];
sx q[2];
rz(-2.5072012) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.40174088) q[1];
sx q[1];
rz(-1.051552) q[1];
sx q[1];
rz(-0.3335764) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.548382) q[3];
sx q[3];
rz(-2.3311391) q[3];
sx q[3];
rz(2.5814501) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.2737736) q[2];
sx q[2];
rz(-1.5742233) q[2];
sx q[2];
rz(1.1567953) q[2];
rz(-2.6371238) q[3];
sx q[3];
rz(-2.1095095) q[3];
sx q[3];
rz(-2.5696136) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2061763) q[0];
sx q[0];
rz(-1.529006) q[0];
sx q[0];
rz(-2.5231498) q[0];
rz(0.84856021) q[1];
sx q[1];
rz(-2.6237374) q[1];
sx q[1];
rz(2.3458164) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9002118) q[0];
sx q[0];
rz(-2.6707895) q[0];
sx q[0];
rz(-1.5174675) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3168542) q[2];
sx q[2];
rz(-1.6055664) q[2];
sx q[2];
rz(-2.8512213) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.8164715) q[1];
sx q[1];
rz(-2.7181648) q[1];
sx q[1];
rz(0.85514833) q[1];
rz(-pi) q[2];
rz(-0.14988092) q[3];
sx q[3];
rz(-1.1544041) q[3];
sx q[3];
rz(-0.021377953) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.985118) q[2];
sx q[2];
rz(-1.5072344) q[2];
sx q[2];
rz(-1.1953243) q[2];
rz(-0.34504238) q[3];
sx q[3];
rz(-2.2796567) q[3];
sx q[3];
rz(2.2752458) q[3];
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
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7095551) q[0];
sx q[0];
rz(-0.23831743) q[0];
sx q[0];
rz(-0.077614345) q[0];
rz(1.2331102) q[1];
sx q[1];
rz(-1.0401007) q[1];
sx q[1];
rz(2.3495823) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4302942) q[0];
sx q[0];
rz(-0.92204198) q[0];
sx q[0];
rz(-2.4374967) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9819005) q[2];
sx q[2];
rz(-1.4174685) q[2];
sx q[2];
rz(0.13736573) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.35490552) q[1];
sx q[1];
rz(-2.6514158) q[1];
sx q[1];
rz(1.9016983) q[1];
x q[2];
rz(-2.3378793) q[3];
sx q[3];
rz(-1.4856385) q[3];
sx q[3];
rz(-0.61060315) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.21976694) q[2];
sx q[2];
rz(-2.2302088) q[2];
sx q[2];
rz(-1.0731953) q[2];
rz(-2.8940767) q[3];
sx q[3];
rz(-1.211834) q[3];
sx q[3];
rz(-0.016935067) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76219227) q[0];
sx q[0];
rz(-1.301855) q[0];
sx q[0];
rz(-0.72574885) q[0];
rz(-0.87860592) q[1];
sx q[1];
rz(-0.85522063) q[1];
sx q[1];
rz(-1.9488889) q[1];
rz(2.939777) q[2];
sx q[2];
rz(-1.0952264) q[2];
sx q[2];
rz(0.14080382) q[2];
rz(-1.5069275) q[3];
sx q[3];
rz(-0.72790925) q[3];
sx q[3];
rz(2.6016605) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
