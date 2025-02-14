OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.32970348) q[0];
sx q[0];
rz(-2.831037) q[0];
sx q[0];
rz(1.9429053) q[0];
rz(-2.0074453) q[1];
sx q[1];
rz(-0.79678798) q[1];
sx q[1];
rz(0.30153433) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3350007) q[0];
sx q[0];
rz(-1.3863239) q[0];
sx q[0];
rz(-2.9770699) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9542175) q[2];
sx q[2];
rz(-1.2760386) q[2];
sx q[2];
rz(-2.8127138) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.63362279) q[1];
sx q[1];
rz(-0.75532856) q[1];
sx q[1];
rz(1.669674) q[1];
x q[2];
rz(1.1490378) q[3];
sx q[3];
rz(-0.43614498) q[3];
sx q[3];
rz(-2.8835322) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.6887168) q[2];
sx q[2];
rz(-1.0998925) q[2];
sx q[2];
rz(-2.0067046) q[2];
rz(-0.12198837) q[3];
sx q[3];
rz(-0.9361836) q[3];
sx q[3];
rz(0.056099135) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4695456) q[0];
sx q[0];
rz(-2.8944954) q[0];
sx q[0];
rz(-3.1392414) q[0];
rz(-2.7408842) q[1];
sx q[1];
rz(-1.7727163) q[1];
sx q[1];
rz(-2.2854663) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5120289) q[0];
sx q[0];
rz(-2.1855547) q[0];
sx q[0];
rz(0.31350664) q[0];
rz(-pi) q[1];
rz(-2.8555774) q[2];
sx q[2];
rz(-1.0723503) q[2];
sx q[2];
rz(1.1666544) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.5634671) q[1];
sx q[1];
rz(-1.5208066) q[1];
sx q[1];
rz(-2.0983138) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.54243954) q[3];
sx q[3];
rz(-0.63839165) q[3];
sx q[3];
rz(2.8505489) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.037228435) q[2];
sx q[2];
rz(-1.2040441) q[2];
sx q[2];
rz(-2.6715703) q[2];
rz(-0.59703279) q[3];
sx q[3];
rz(-3.0266914) q[3];
sx q[3];
rz(1.443583) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3700579) q[0];
sx q[0];
rz(-0.4929339) q[0];
sx q[0];
rz(-2.9135627) q[0];
rz(0.44949624) q[1];
sx q[1];
rz(-1.0003961) q[1];
sx q[1];
rz(-2.1953348) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3920604) q[0];
sx q[0];
rz(-2.6596053) q[0];
sx q[0];
rz(-2.1493692) q[0];
x q[1];
rz(0.94257943) q[2];
sx q[2];
rz(-1.0120856) q[2];
sx q[2];
rz(-3.1329346) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.98639224) q[1];
sx q[1];
rz(-0.3861392) q[1];
sx q[1];
rz(-0.50215747) q[1];
x q[2];
rz(2.0722996) q[3];
sx q[3];
rz(-0.32078241) q[3];
sx q[3];
rz(-0.71170413) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.2366644) q[2];
sx q[2];
rz(-1.6764574) q[2];
sx q[2];
rz(1.4795335) q[2];
rz(-2.9605401) q[3];
sx q[3];
rz(-0.79020399) q[3];
sx q[3];
rz(-0.886206) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(0.36412305) q[0];
sx q[0];
rz(-1.5793707) q[0];
sx q[0];
rz(2.2430578) q[0];
rz(0.50615519) q[1];
sx q[1];
rz(-1.8471142) q[1];
sx q[1];
rz(-1.6815965) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0943992) q[0];
sx q[0];
rz(-1.5051821) q[0];
sx q[0];
rz(-0.53334566) q[0];
rz(-pi) q[1];
rz(-1.1300259) q[2];
sx q[2];
rz(-0.87723644) q[2];
sx q[2];
rz(2.9206985) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.0467028) q[1];
sx q[1];
rz(-3.0071435) q[1];
sx q[1];
rz(-2.1220757) q[1];
rz(-pi) q[2];
rz(-1.4921032) q[3];
sx q[3];
rz(-0.65641145) q[3];
sx q[3];
rz(-2.479913) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.8755181) q[2];
sx q[2];
rz(-0.47553277) q[2];
sx q[2];
rz(-1.5024705) q[2];
rz(-1.255704) q[3];
sx q[3];
rz(-1.7399104) q[3];
sx q[3];
rz(-3.0953395) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9083549) q[0];
sx q[0];
rz(-2.4112356) q[0];
sx q[0];
rz(2.9610942) q[0];
rz(-1.7620697) q[1];
sx q[1];
rz(-2.5738398) q[1];
sx q[1];
rz(-2.0464121) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3334072) q[0];
sx q[0];
rz(-1.6327724) q[0];
sx q[0];
rz(-3.005885) q[0];
x q[1];
rz(3.0600431) q[2];
sx q[2];
rz(-2.3429541) q[2];
sx q[2];
rz(2.7872686) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.74873456) q[1];
sx q[1];
rz(-1.431141) q[1];
sx q[1];
rz(-0.77844324) q[1];
rz(-pi) q[2];
rz(1.7322695) q[3];
sx q[3];
rz(-2.4465909) q[3];
sx q[3];
rz(1.7308066) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.77202648) q[2];
sx q[2];
rz(-0.84731805) q[2];
sx q[2];
rz(1.0687211) q[2];
rz(-2.7145743) q[3];
sx q[3];
rz(-2.2618099) q[3];
sx q[3];
rz(1.8900185) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1123947) q[0];
sx q[0];
rz(-2.2258832) q[0];
sx q[0];
rz(0.29119626) q[0];
rz(-1.683782) q[1];
sx q[1];
rz(-1.0271415) q[1];
sx q[1];
rz(0.88868946) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5425974) q[0];
sx q[0];
rz(-2.3610736) q[0];
sx q[0];
rz(0.88135834) q[0];
rz(-pi) q[1];
rz(2.3843147) q[2];
sx q[2];
rz(-1.5520397) q[2];
sx q[2];
rz(0.64964408) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.86119628) q[1];
sx q[1];
rz(-0.79842192) q[1];
sx q[1];
rz(-2.3905928) q[1];
rz(1.8848041) q[3];
sx q[3];
rz(-1.2824138) q[3];
sx q[3];
rz(-2.7183804) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.7500744) q[2];
sx q[2];
rz(-1.3041648) q[2];
sx q[2];
rz(0.94775003) q[2];
rz(2.9446972) q[3];
sx q[3];
rz(-2.9010549) q[3];
sx q[3];
rz(-0.31015629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
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
rz(2.8295558) q[0];
sx q[0];
rz(-2.0289679) q[0];
sx q[0];
rz(2.6716676) q[0];
rz(-1.4620818) q[1];
sx q[1];
rz(-1.3009678) q[1];
sx q[1];
rz(-1.7376815) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1184153) q[0];
sx q[0];
rz(-1.0141393) q[0];
sx q[0];
rz(-2.7283637) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.092488) q[2];
sx q[2];
rz(-1.7564536) q[2];
sx q[2];
rz(-2.3379679) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.1199855) q[1];
sx q[1];
rz(-2.48263) q[1];
sx q[1];
rz(-0.95665415) q[1];
rz(-pi) q[2];
x q[2];
rz(0.49755712) q[3];
sx q[3];
rz(-2.4973739) q[3];
sx q[3];
rz(-1.8807172) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.5002084) q[2];
sx q[2];
rz(-0.83157867) q[2];
sx q[2];
rz(0.10759648) q[2];
rz(-2.522116) q[3];
sx q[3];
rz(-1.8457125) q[3];
sx q[3];
rz(-1.3639601) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7463995) q[0];
sx q[0];
rz(-1.0522319) q[0];
sx q[0];
rz(-0.25303823) q[0];
rz(-0.088317618) q[1];
sx q[1];
rz(-0.87366906) q[1];
sx q[1];
rz(-2.1926682) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.8924361) q[0];
sx q[0];
rz(-1.2741486) q[0];
sx q[0];
rz(0.15693024) q[0];
rz(-pi) q[1];
rz(-1.9789655) q[2];
sx q[2];
rz(-1.7498831) q[2];
sx q[2];
rz(0.31291747) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.9946549) q[1];
sx q[1];
rz(-2.6587464) q[1];
sx q[1];
rz(1.6864683) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.003951) q[3];
sx q[3];
rz(-1.1205871) q[3];
sx q[3];
rz(2.4382044) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.39763149) q[2];
sx q[2];
rz(-0.88699114) q[2];
sx q[2];
rz(-0.31069791) q[2];
rz(0.24142309) q[3];
sx q[3];
rz(-1.2025236) q[3];
sx q[3];
rz(2.0588622) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0017515) q[0];
sx q[0];
rz(-0.40626353) q[0];
sx q[0];
rz(1.9061506) q[0];
rz(2.6605117) q[1];
sx q[1];
rz(-2.1886539) q[1];
sx q[1];
rz(1.4280041) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7803376) q[0];
sx q[0];
rz(-1.5819823) q[0];
sx q[0];
rz(0.7594603) q[0];
x q[1];
rz(-2.253654) q[2];
sx q[2];
rz(-0.35338923) q[2];
sx q[2];
rz(1.8977752) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.2216827) q[1];
sx q[1];
rz(-2.7836728) q[1];
sx q[1];
rz(-2.329925) q[1];
rz(-2.7747267) q[3];
sx q[3];
rz(-0.75682029) q[3];
sx q[3];
rz(-3.0928825) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.1741751) q[2];
sx q[2];
rz(-2.316541) q[2];
sx q[2];
rz(-0.027912557) q[2];
rz(-0.86999718) q[3];
sx q[3];
rz(-0.78012192) q[3];
sx q[3];
rz(-1.9546668) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
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
rz(0.025539909) q[0];
sx q[0];
rz(-1.5486131) q[0];
sx q[0];
rz(-1.0593587) q[0];
rz(-1.7268044) q[1];
sx q[1];
rz(-1.6635514) q[1];
sx q[1];
rz(-0.47763225) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.85738289) q[0];
sx q[0];
rz(-0.83440953) q[0];
sx q[0];
rz(1.7718023) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.576409) q[2];
sx q[2];
rz(-2.8057348) q[2];
sx q[2];
rz(0.72959057) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.3463626) q[1];
sx q[1];
rz(-1.6246965) q[1];
sx q[1];
rz(1.6832215) q[1];
rz(-0.60250256) q[3];
sx q[3];
rz(-2.2244766) q[3];
sx q[3];
rz(2.9127171) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.16964218) q[2];
sx q[2];
rz(-2.1596238) q[2];
sx q[2];
rz(-0.36894813) q[2];
rz(-1.87489) q[3];
sx q[3];
rz(-0.72487512) q[3];
sx q[3];
rz(0.96316159) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0030768) q[0];
sx q[0];
rz(-1.2110447) q[0];
sx q[0];
rz(-1.3883653) q[0];
rz(-0.44454642) q[1];
sx q[1];
rz(-0.81041705) q[1];
sx q[1];
rz(3.0141426) q[1];
rz(0.62194555) q[2];
sx q[2];
rz(-2.6126886) q[2];
sx q[2];
rz(-2.9168203) q[2];
rz(-2.6740549) q[3];
sx q[3];
rz(-1.553592) q[3];
sx q[3];
rz(-1.6910465) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
