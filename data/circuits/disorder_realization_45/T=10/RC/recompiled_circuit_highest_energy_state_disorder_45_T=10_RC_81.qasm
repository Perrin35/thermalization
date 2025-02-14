OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.9396012) q[0];
sx q[0];
rz(-0.48299462) q[0];
sx q[0];
rz(1.6311837) q[0];
rz(3.0022439) q[1];
sx q[1];
rz(-2.5820093) q[1];
sx q[1];
rz(0.72996563) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.18099526) q[0];
sx q[0];
rz(-0.9708403) q[0];
sx q[0];
rz(0.17790312) q[0];
x q[1];
rz(0.32663871) q[2];
sx q[2];
rz(-0.62033287) q[2];
sx q[2];
rz(-2.6503944) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(3.0180359) q[1];
sx q[1];
rz(-2.5606321) q[1];
sx q[1];
rz(3.0835129) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2830986) q[3];
sx q[3];
rz(-1.7904864) q[3];
sx q[3];
rz(-0.9385329) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.0240747) q[2];
sx q[2];
rz(-0.2710318) q[2];
sx q[2];
rz(0.81895858) q[2];
rz(0.0062746127) q[3];
sx q[3];
rz(-1.2333906) q[3];
sx q[3];
rz(-0.98244572) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5821424) q[0];
sx q[0];
rz(-1.5511976) q[0];
sx q[0];
rz(0.5994125) q[0];
rz(-0.86743152) q[1];
sx q[1];
rz(-1.0299094) q[1];
sx q[1];
rz(-0.74554602) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.065577995) q[0];
sx q[0];
rz(-1.6981372) q[0];
sx q[0];
rz(-0.29177427) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3815161) q[2];
sx q[2];
rz(-1.809263) q[2];
sx q[2];
rz(0.60205215) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.9245431) q[1];
sx q[1];
rz(-2.1991208) q[1];
sx q[1];
rz(-1.7679067) q[1];
rz(-2.638444) q[3];
sx q[3];
rz(-2.1981578) q[3];
sx q[3];
rz(-0.37093758) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.6856689) q[2];
sx q[2];
rz(-0.29752877) q[2];
sx q[2];
rz(-1.4073184) q[2];
rz(2.2972441) q[3];
sx q[3];
rz(-0.83365369) q[3];
sx q[3];
rz(-1.0367397) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5582964) q[0];
sx q[0];
rz(-1.1953657) q[0];
sx q[0];
rz(-1.012828) q[0];
rz(-2.5335675) q[1];
sx q[1];
rz(-1.5826179) q[1];
sx q[1];
rz(-1.2695405) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0657261) q[0];
sx q[0];
rz(-1.1320496) q[0];
sx q[0];
rz(-0.42515305) q[0];
rz(-pi) q[1];
x q[1];
rz(3.059194) q[2];
sx q[2];
rz(-0.90148704) q[2];
sx q[2];
rz(1.5562039) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.041965466) q[1];
sx q[1];
rz(-1.7668952) q[1];
sx q[1];
rz(-0.72814299) q[1];
rz(1.2392912) q[3];
sx q[3];
rz(-1.1024144) q[3];
sx q[3];
rz(1.4862332) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.515392) q[2];
sx q[2];
rz(-1.0582558) q[2];
sx q[2];
rz(-1.8499648) q[2];
rz(-3.1332704) q[3];
sx q[3];
rz(-0.95299995) q[3];
sx q[3];
rz(-2.9279809) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13852791) q[0];
sx q[0];
rz(-0.55532885) q[0];
sx q[0];
rz(1.7648765) q[0];
rz(1.6142913) q[1];
sx q[1];
rz(-1.9381899) q[1];
sx q[1];
rz(-2.6208904) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0390065) q[0];
sx q[0];
rz(-1.163536) q[0];
sx q[0];
rz(-2.7558221) q[0];
rz(-2.1200373) q[2];
sx q[2];
rz(-0.62473544) q[2];
sx q[2];
rz(-2.3083589) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.6286875) q[1];
sx q[1];
rz(-2.4081552) q[1];
sx q[1];
rz(-0.67104407) q[1];
x q[2];
rz(-2.2124452) q[3];
sx q[3];
rz(-2.3722674) q[3];
sx q[3];
rz(-1.3206467) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.6167831) q[2];
sx q[2];
rz(-0.7889792) q[2];
sx q[2];
rz(-1.7899803) q[2];
rz(1.0194408) q[3];
sx q[3];
rz(-2.5717058) q[3];
sx q[3];
rz(1.1714237) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5919507) q[0];
sx q[0];
rz(-1.6787981) q[0];
sx q[0];
rz(-0.098966448) q[0];
rz(0.70676604) q[1];
sx q[1];
rz(-0.87044972) q[1];
sx q[1];
rz(-2.6720572) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5564559) q[0];
sx q[0];
rz(-3.0017108) q[0];
sx q[0];
rz(-1.3993457) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7623873) q[2];
sx q[2];
rz(-2.0742886) q[2];
sx q[2];
rz(2.6852754) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.3217056) q[1];
sx q[1];
rz(-2.4155136) q[1];
sx q[1];
rz(-0.5989845) q[1];
rz(-2.585586) q[3];
sx q[3];
rz(-1.4463498) q[3];
sx q[3];
rz(-0.57226244) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.88064757) q[2];
sx q[2];
rz(-1.3346883) q[2];
sx q[2];
rz(-1.3060695) q[2];
rz(2.7169054) q[3];
sx q[3];
rz(-1.3771907) q[3];
sx q[3];
rz(-2.2556321) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.11923085) q[0];
sx q[0];
rz(-0.62218085) q[0];
sx q[0];
rz(-1.3223883) q[0];
rz(-1.127683) q[1];
sx q[1];
rz(-1.5195945) q[1];
sx q[1];
rz(-1.7599531) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.19162543) q[0];
sx q[0];
rz(-2.3268496) q[0];
sx q[0];
rz(0.11884584) q[0];
x q[1];
rz(0.99389771) q[2];
sx q[2];
rz(-2.0303147) q[2];
sx q[2];
rz(1.7651059) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.1028916) q[1];
sx q[1];
rz(-2.3159317) q[1];
sx q[1];
rz(1.3172512) q[1];
x q[2];
rz(-1.6448037) q[3];
sx q[3];
rz(-1.8714219) q[3];
sx q[3];
rz(-0.12412503) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.7606925) q[2];
sx q[2];
rz(-1.6221294) q[2];
sx q[2];
rz(2.6603928) q[2];
rz(0.5091269) q[3];
sx q[3];
rz(-2.9042518) q[3];
sx q[3];
rz(2.7595162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5803489) q[0];
sx q[0];
rz(-0.48155293) q[0];
sx q[0];
rz(-0.089381889) q[0];
rz(-2.0629758) q[1];
sx q[1];
rz(-1.4559454) q[1];
sx q[1];
rz(0.99348974) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.29879323) q[0];
sx q[0];
rz(-1.5325108) q[0];
sx q[0];
rz(-2.3596615) q[0];
rz(-pi) q[1];
rz(2.2682796) q[2];
sx q[2];
rz(-0.9160348) q[2];
sx q[2];
rz(0.26604929) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.9536994) q[1];
sx q[1];
rz(-2.7479798) q[1];
sx q[1];
rz(1.1792437) q[1];
rz(-pi) q[2];
rz(2.5105623) q[3];
sx q[3];
rz(-2.0027805) q[3];
sx q[3];
rz(-0.61842266) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.79954687) q[2];
sx q[2];
rz(-0.61508238) q[2];
sx q[2];
rz(-2.7395524) q[2];
rz(-1.0953995) q[3];
sx q[3];
rz(-1.2145372) q[3];
sx q[3];
rz(-0.57797617) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(-0.47356975) q[0];
sx q[0];
rz(-0.80900017) q[0];
sx q[0];
rz(1.8079669) q[0];
rz(-2.2857621) q[1];
sx q[1];
rz(-1.7355093) q[1];
sx q[1];
rz(0.91748253) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.30675754) q[0];
sx q[0];
rz(-1.4687755) q[0];
sx q[0];
rz(-1.7708885) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.11882527) q[2];
sx q[2];
rz(-2.1911494) q[2];
sx q[2];
rz(0.22681776) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.6010103) q[1];
sx q[1];
rz(-1.9153908) q[1];
sx q[1];
rz(0.2356727) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0621655) q[3];
sx q[3];
rz(-0.98782238) q[3];
sx q[3];
rz(0.82750083) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.63460073) q[2];
sx q[2];
rz(-1.736182) q[2];
sx q[2];
rz(1.290192) q[2];
rz(-0.90977943) q[3];
sx q[3];
rz(-1.4434283) q[3];
sx q[3];
rz(0.040987404) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
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
rz(0.23583394) q[0];
sx q[0];
rz(-1.1474778) q[0];
sx q[0];
rz(-0.27715096) q[0];
rz(1.2241036) q[1];
sx q[1];
rz(-1.5382907) q[1];
sx q[1];
rz(1.3714429) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.087589892) q[0];
sx q[0];
rz(-1.3773019) q[0];
sx q[0];
rz(-0.67275472) q[0];
rz(-pi) q[1];
rz(0.54746898) q[2];
sx q[2];
rz(-1.1927422) q[2];
sx q[2];
rz(-2.7965656) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.38199785) q[1];
sx q[1];
rz(-1.7773668) q[1];
sx q[1];
rz(2.2269511) q[1];
rz(2.9670466) q[3];
sx q[3];
rz(-2.2851599) q[3];
sx q[3];
rz(0.81596953) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.93854967) q[2];
sx q[2];
rz(-2.2800192) q[2];
sx q[2];
rz(2.4533563) q[2];
rz(0.31442434) q[3];
sx q[3];
rz(-0.36948547) q[3];
sx q[3];
rz(1.2615874) q[3];
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
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7670583) q[0];
sx q[0];
rz(-2.2735167) q[0];
sx q[0];
rz(0.57149291) q[0];
rz(-2.4608965) q[1];
sx q[1];
rz(-1.9711767) q[1];
sx q[1];
rz(0.64819711) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2637973) q[0];
sx q[0];
rz(-1.5476942) q[0];
sx q[0];
rz(-0.013851555) q[0];
rz(-pi) q[1];
rz(0.86096455) q[2];
sx q[2];
rz(-2.2662244) q[2];
sx q[2];
rz(-2.8784213) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.1614589) q[1];
sx q[1];
rz(-2.4422283) q[1];
sx q[1];
rz(0.034566391) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.51697124) q[3];
sx q[3];
rz(-1.7423769) q[3];
sx q[3];
rz(-1.3763381) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.8615243) q[2];
sx q[2];
rz(-1.0736977) q[2];
sx q[2];
rz(1.466922) q[2];
rz(-0.17659771) q[3];
sx q[3];
rz(-0.64591518) q[3];
sx q[3];
rz(-1.8471898) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8661154) q[0];
sx q[0];
rz(-1.7452411) q[0];
sx q[0];
rz(-1.2757975) q[0];
rz(-2.7571309) q[1];
sx q[1];
rz(-1.5059595) q[1];
sx q[1];
rz(-0.62677871) q[1];
rz(-2.7914417) q[2];
sx q[2];
rz(-1.2731009) q[2];
sx q[2];
rz(0.44558744) q[2];
rz(2.3220358) q[3];
sx q[3];
rz(-0.71376505) q[3];
sx q[3];
rz(-2.5160088) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
