OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.80823094) q[0];
sx q[0];
rz(-2.1699177) q[0];
sx q[0];
rz(1.4810286) q[0];
rz(-0.12227585) q[1];
sx q[1];
rz(3.0552157) q[1];
sx q[1];
rz(9.4019158) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.45863736) q[0];
sx q[0];
rz(-1.6874755) q[0];
sx q[0];
rz(1.7200243) q[0];
x q[1];
rz(2.4542698) q[2];
sx q[2];
rz(-0.63185531) q[2];
sx q[2];
rz(-0.99422115) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.37576807) q[1];
sx q[1];
rz(-1.4210912) q[1];
sx q[1];
rz(1.9181812) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1894737) q[3];
sx q[3];
rz(-1.3631571) q[3];
sx q[3];
rz(1.729897) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.4188529) q[2];
sx q[2];
rz(-0.8272233) q[2];
sx q[2];
rz(-0.76618761) q[2];
rz(0.17151672) q[3];
sx q[3];
rz(-0.73232108) q[3];
sx q[3];
rz(-0.58656251) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3294753) q[0];
sx q[0];
rz(-0.43918878) q[0];
sx q[0];
rz(0.086659327) q[0];
rz(-0.86241972) q[1];
sx q[1];
rz(-2.5984867) q[1];
sx q[1];
rz(3.0564953) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.800195) q[0];
sx q[0];
rz(-1.5872103) q[0];
sx q[0];
rz(-0.62231681) q[0];
rz(-pi) q[1];
x q[1];
rz(0.9688103) q[2];
sx q[2];
rz(-1.3170871) q[2];
sx q[2];
rz(1.3041376) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.6669586) q[1];
sx q[1];
rz(-2.9999795) q[1];
sx q[1];
rz(-2.2255564) q[1];
rz(-0.99382932) q[3];
sx q[3];
rz(-2.5964063) q[3];
sx q[3];
rz(-0.59515566) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.8791085) q[2];
sx q[2];
rz(-1.5510677) q[2];
sx q[2];
rz(1.6764199) q[2];
rz(-1.1761459) q[3];
sx q[3];
rz(-2.2361103) q[3];
sx q[3];
rz(-0.24648497) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48433205) q[0];
sx q[0];
rz(-3.0069139) q[0];
sx q[0];
rz(0.9056257) q[0];
rz(2.8088645) q[1];
sx q[1];
rz(-2.2867124) q[1];
sx q[1];
rz(-1.3844301) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8850088) q[0];
sx q[0];
rz(-2.6589767) q[0];
sx q[0];
rz(-0.17871876) q[0];
rz(-pi) q[1];
rz(-2.7483814) q[2];
sx q[2];
rz(-1.097743) q[2];
sx q[2];
rz(-0.87393239) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.3789931) q[1];
sx q[1];
rz(-1.5484719) q[1];
sx q[1];
rz(-0.75548817) q[1];
rz(-pi) q[2];
rz(-1.0043886) q[3];
sx q[3];
rz(-1.1096138) q[3];
sx q[3];
rz(0.084463488) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.24421346) q[2];
sx q[2];
rz(-2.801557) q[2];
sx q[2];
rz(-1.3282233) q[2];
rz(1.7437079) q[3];
sx q[3];
rz(-0.54026794) q[3];
sx q[3];
rz(-3.0298997) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.06552799) q[0];
sx q[0];
rz(-0.17834839) q[0];
sx q[0];
rz(-0.67101014) q[0];
rz(-1.7975851) q[1];
sx q[1];
rz(-1.1616511) q[1];
sx q[1];
rz(2.9342594) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3610791) q[0];
sx q[0];
rz(-1.5618389) q[0];
sx q[0];
rz(1.5293967) q[0];
rz(-pi) q[1];
rz(3.1111654) q[2];
sx q[2];
rz(-1.3969587) q[2];
sx q[2];
rz(-2.0174842) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.3934717) q[1];
sx q[1];
rz(-1.8078184) q[1];
sx q[1];
rz(-2.9841828) q[1];
x q[2];
rz(-1.585235) q[3];
sx q[3];
rz(-1.0591905) q[3];
sx q[3];
rz(0.59359854) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.8308782) q[2];
sx q[2];
rz(-2.0534616) q[2];
sx q[2];
rz(-2.3332398) q[2];
rz(0.80777848) q[3];
sx q[3];
rz(-0.72179663) q[3];
sx q[3];
rz(-2.9479153) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59072524) q[0];
sx q[0];
rz(-2.4243675) q[0];
sx q[0];
rz(-0.89548683) q[0];
rz(-2.8764309) q[1];
sx q[1];
rz(-0.83509713) q[1];
sx q[1];
rz(2.6079544) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69913188) q[0];
sx q[0];
rz(-1.2958382) q[0];
sx q[0];
rz(1.2752227) q[0];
rz(-pi) q[1];
rz(-0.76782121) q[2];
sx q[2];
rz(-1.7349093) q[2];
sx q[2];
rz(0.62190157) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.4644949) q[1];
sx q[1];
rz(-0.2650731) q[1];
sx q[1];
rz(-0.58880834) q[1];
x q[2];
rz(-0.20547262) q[3];
sx q[3];
rz(-2.1901202) q[3];
sx q[3];
rz(1.6831786) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.2395893) q[2];
sx q[2];
rz(-1.5782372) q[2];
sx q[2];
rz(2.5732102) q[2];
rz(1.2166294) q[3];
sx q[3];
rz(-2.670848) q[3];
sx q[3];
rz(-0.34354982) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6333273) q[0];
sx q[0];
rz(-2.2631622) q[0];
sx q[0];
rz(-1.9653962) q[0];
rz(1.5559224) q[1];
sx q[1];
rz(-0.95247477) q[1];
sx q[1];
rz(1.0046545) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3096136) q[0];
sx q[0];
rz(-1.2011315) q[0];
sx q[0];
rz(1.2020822) q[0];
rz(2.7447694) q[2];
sx q[2];
rz(-1.6453711) q[2];
sx q[2];
rz(2.4900988) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-3.1130044) q[1];
sx q[1];
rz(-0.60739809) q[1];
sx q[1];
rz(-0.44087704) q[1];
x q[2];
rz(-1.1569571) q[3];
sx q[3];
rz(-1.5817747) q[3];
sx q[3];
rz(1.4423086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.4040318) q[2];
sx q[2];
rz(-2.1514251) q[2];
sx q[2];
rz(-0.44815865) q[2];
rz(-0.29799497) q[3];
sx q[3];
rz(-2.4500676) q[3];
sx q[3];
rz(-0.34860778) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3649243) q[0];
sx q[0];
rz(-1.8402599) q[0];
sx q[0];
rz(-2.4834494) q[0];
rz(-1.2843885) q[1];
sx q[1];
rz(-2.6874266) q[1];
sx q[1];
rz(2.4694494) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2015398) q[0];
sx q[0];
rz(-1.7367898) q[0];
sx q[0];
rz(0.0075017651) q[0];
rz(-1.1291814) q[2];
sx q[2];
rz(-1.9600944) q[2];
sx q[2];
rz(-1.0600952) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.3232279) q[1];
sx q[1];
rz(-2.0548901) q[1];
sx q[1];
rz(1.9869884) q[1];
x q[2];
rz(2.4403205) q[3];
sx q[3];
rz(-1.4605195) q[3];
sx q[3];
rz(2.9993204) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.60823524) q[2];
sx q[2];
rz(-0.63295263) q[2];
sx q[2];
rz(0.39880025) q[2];
rz(-0.82459015) q[3];
sx q[3];
rz(-1.7493068) q[3];
sx q[3];
rz(2.5430172) q[3];
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
sx q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7010715) q[0];
sx q[0];
rz(-2.7109787) q[0];
sx q[0];
rz(2.9833802) q[0];
rz(-2.054706) q[1];
sx q[1];
rz(-2.465076) q[1];
sx q[1];
rz(2.3349082) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8950302) q[0];
sx q[0];
rz(-0.44996214) q[0];
sx q[0];
rz(-1.8438086) q[0];
rz(-1.6786472) q[2];
sx q[2];
rz(-1.7753522) q[2];
sx q[2];
rz(0.66609913) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.1572691) q[1];
sx q[1];
rz(-1.8586129) q[1];
sx q[1];
rz(2.4411574) q[1];
rz(-pi) q[2];
x q[2];
rz(0.032012352) q[3];
sx q[3];
rz(-1.2700998) q[3];
sx q[3];
rz(0.34919448) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.5640101) q[2];
sx q[2];
rz(-1.6489886) q[2];
sx q[2];
rz(-2.1419443) q[2];
rz(-3.0380761) q[3];
sx q[3];
rz(-3.0045356) q[3];
sx q[3];
rz(-2.0172393) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.12061159) q[0];
sx q[0];
rz(-2.3636901) q[0];
sx q[0];
rz(-2.4840684) q[0];
rz(-2.855037) q[1];
sx q[1];
rz(-0.92266881) q[1];
sx q[1];
rz(-0.34067571) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1239615) q[0];
sx q[0];
rz(-2.258856) q[0];
sx q[0];
rz(0.4171564) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.22206837) q[2];
sx q[2];
rz(-2.1246315) q[2];
sx q[2];
rz(1.008322) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.77010158) q[1];
sx q[1];
rz(-2.1226468) q[1];
sx q[1];
rz(2.9279207) q[1];
rz(-pi) q[2];
rz(-2.4386028) q[3];
sx q[3];
rz(-0.25745108) q[3];
sx q[3];
rz(1.5233745) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.3866773) q[2];
sx q[2];
rz(-1.200054) q[2];
sx q[2];
rz(-2.6436451) q[2];
rz(2.8399816) q[3];
sx q[3];
rz(-0.40062723) q[3];
sx q[3];
rz(0.51914674) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.72934812) q[0];
sx q[0];
rz(-0.059818581) q[0];
sx q[0];
rz(2.8630032) q[0];
rz(0.57922286) q[1];
sx q[1];
rz(-0.93943739) q[1];
sx q[1];
rz(3.0648807) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.840832) q[0];
sx q[0];
rz(-1.4474086) q[0];
sx q[0];
rz(-1.2542017) q[0];
rz(1.3086583) q[2];
sx q[2];
rz(-0.58912504) q[2];
sx q[2];
rz(-0.72074705) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.0753239) q[1];
sx q[1];
rz(-1.2053718) q[1];
sx q[1];
rz(1.5484023) q[1];
x q[2];
rz(-2.6807908) q[3];
sx q[3];
rz(-0.19690234) q[3];
sx q[3];
rz(2.868696) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.77438337) q[2];
sx q[2];
rz(-0.70720208) q[2];
sx q[2];
rz(0.45551604) q[2];
rz(-2.7101743) q[3];
sx q[3];
rz(-0.12532561) q[3];
sx q[3];
rz(-0.063522696) q[3];
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
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1800304) q[0];
sx q[0];
rz(-1.7407692) q[0];
sx q[0];
rz(-2.2041007) q[0];
rz(-0.58615276) q[1];
sx q[1];
rz(-1.6393607) q[1];
sx q[1];
rz(1.6929109) q[1];
rz(2.8557599) q[2];
sx q[2];
rz(-2.0252068) q[2];
sx q[2];
rz(0.037638738) q[2];
rz(2.4424845) q[3];
sx q[3];
rz(-0.53634488) q[3];
sx q[3];
rz(2.8041822) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
