OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.3396575) q[0];
sx q[0];
rz(-0.8987838) q[0];
sx q[0];
rz(1.3026097) q[0];
rz(-3.0175735) q[1];
sx q[1];
rz(-1.7350585) q[1];
sx q[1];
rz(1.6279434) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0078137579) q[0];
sx q[0];
rz(-1.2764751) q[0];
sx q[0];
rz(1.9208287) q[0];
rz(-2.0997203) q[2];
sx q[2];
rz(-2.0588819) q[2];
sx q[2];
rz(1.5785335) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.6484272) q[1];
sx q[1];
rz(-1.9110288) q[1];
sx q[1];
rz(-3.0837584) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7645993) q[3];
sx q[3];
rz(-0.97929472) q[3];
sx q[3];
rz(1.5821004) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.0836432) q[2];
sx q[2];
rz(-1.8079855) q[2];
sx q[2];
rz(-2.7570214) q[2];
rz(-0.40258506) q[3];
sx q[3];
rz(-1.4339002) q[3];
sx q[3];
rz(-2.3334077) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3835417) q[0];
sx q[0];
rz(-1.533968) q[0];
sx q[0];
rz(0.59355271) q[0];
rz(-1.3213762) q[1];
sx q[1];
rz(-1.079419) q[1];
sx q[1];
rz(3.1022601) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6932019) q[0];
sx q[0];
rz(-1.9422008) q[0];
sx q[0];
rz(2.0601963) q[0];
x q[1];
rz(-2.987542) q[2];
sx q[2];
rz(-0.64214911) q[2];
sx q[2];
rz(-2.1352445) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.5099098) q[1];
sx q[1];
rz(-1.0541774) q[1];
sx q[1];
rz(-1.6277908) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.049841471) q[3];
sx q[3];
rz(-0.54996544) q[3];
sx q[3];
rz(-1.0483688) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.3010657) q[2];
sx q[2];
rz(-1.4378005) q[2];
sx q[2];
rz(0.62704101) q[2];
rz(-2.7035233) q[3];
sx q[3];
rz(-1.4984683) q[3];
sx q[3];
rz(-2.0048678) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
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
rz(0.42191926) q[0];
sx q[0];
rz(-1.0872343) q[0];
sx q[0];
rz(1.3713974) q[0];
rz(-2.5321391) q[1];
sx q[1];
rz(-2.2513159) q[1];
sx q[1];
rz(1.2249464) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6951675) q[0];
sx q[0];
rz(-1.7621627) q[0];
sx q[0];
rz(1.4728949) q[0];
rz(-pi) q[1];
rz(2.742381) q[2];
sx q[2];
rz(-2.2089094) q[2];
sx q[2];
rz(-2.7094584) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.6869192) q[1];
sx q[1];
rz(-0.49263601) q[1];
sx q[1];
rz(0.83322816) q[1];
x q[2];
rz(0.21594343) q[3];
sx q[3];
rz(-0.88679291) q[3];
sx q[3];
rz(-2.0404599) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.3043392) q[2];
sx q[2];
rz(-2.2746268) q[2];
sx q[2];
rz(2.2875817) q[2];
rz(2.617406) q[3];
sx q[3];
rz(-1.5093404) q[3];
sx q[3];
rz(-2.1715651) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.22685856) q[0];
sx q[0];
rz(-0.15758841) q[0];
sx q[0];
rz(-2.3478813) q[0];
rz(0.0018250068) q[1];
sx q[1];
rz(-2.5853214) q[1];
sx q[1];
rz(2.3307641) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9017667) q[0];
sx q[0];
rz(-1.7525273) q[0];
sx q[0];
rz(-0.82518362) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0433719) q[2];
sx q[2];
rz(-0.7458936) q[2];
sx q[2];
rz(2.8371642) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.27426611) q[1];
sx q[1];
rz(-1.9073448) q[1];
sx q[1];
rz(-0.95086581) q[1];
rz(-pi) q[2];
x q[2];
rz(1.879519) q[3];
sx q[3];
rz(-0.44770542) q[3];
sx q[3];
rz(1.6310504) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.4578555) q[2];
sx q[2];
rz(-0.33317864) q[2];
sx q[2];
rz(-1.4972756) q[2];
rz(-1.3582683) q[3];
sx q[3];
rz(-2.0446916) q[3];
sx q[3];
rz(1.8339405) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7233647) q[0];
sx q[0];
rz(-1.6809373) q[0];
sx q[0];
rz(-0.25890589) q[0];
rz(0.12340165) q[1];
sx q[1];
rz(-0.63107189) q[1];
sx q[1];
rz(2.6554328) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.21307316) q[0];
sx q[0];
rz(-2.3384636) q[0];
sx q[0];
rz(-2.2691239) q[0];
rz(-pi) q[1];
rz(2.6888383) q[2];
sx q[2];
rz(-1.448302) q[2];
sx q[2];
rz(1.233135) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.26700324) q[1];
sx q[1];
rz(-1.7591333) q[1];
sx q[1];
rz(1.9397771) q[1];
rz(-pi) q[2];
x q[2];
rz(0.93323243) q[3];
sx q[3];
rz(-0.8094111) q[3];
sx q[3];
rz(-0.903086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.90998021) q[2];
sx q[2];
rz(-2.6284802) q[2];
sx q[2];
rz(1.7463589) q[2];
rz(0.95476556) q[3];
sx q[3];
rz(-1.3937817) q[3];
sx q[3];
rz(-1.0432997) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2024277) q[0];
sx q[0];
rz(-1.9068149) q[0];
sx q[0];
rz(-2.3811316) q[0];
rz(2.3072534) q[1];
sx q[1];
rz(-2.0400679) q[1];
sx q[1];
rz(2.0371425) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0898665) q[0];
sx q[0];
rz(-1.3524922) q[0];
sx q[0];
rz(-1.0132683) q[0];
rz(-2.274674) q[2];
sx q[2];
rz(-0.50627497) q[2];
sx q[2];
rz(-2.3828854) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.23984662) q[1];
sx q[1];
rz(-0.86565986) q[1];
sx q[1];
rz(2.2589931) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1105632) q[3];
sx q[3];
rz(-0.62556872) q[3];
sx q[3];
rz(-0.21504185) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.1793648) q[2];
sx q[2];
rz(-0.70415512) q[2];
sx q[2];
rz(-0.041570138) q[2];
rz(2.9465594) q[3];
sx q[3];
rz(-2.8988367) q[3];
sx q[3];
rz(2.354505) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4389909) q[0];
sx q[0];
rz(-0.85346237) q[0];
sx q[0];
rz(-0.75337291) q[0];
rz(-1.0607177) q[1];
sx q[1];
rz(-0.80442387) q[1];
sx q[1];
rz(-2.9341968) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3531137) q[0];
sx q[0];
rz(-2.1512356) q[0];
sx q[0];
rz(-0.81932318) q[0];
rz(-pi) q[1];
rz(-1.8600402) q[2];
sx q[2];
rz(-1.1221894) q[2];
sx q[2];
rz(-1.9236652) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.3606229) q[1];
sx q[1];
rz(-2.1163673) q[1];
sx q[1];
rz(1.8459794) q[1];
rz(-pi) q[2];
rz(2.9459729) q[3];
sx q[3];
rz(-3.1017711) q[3];
sx q[3];
rz(-2.8012343) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.14219001) q[2];
sx q[2];
rz(-1.9839857) q[2];
sx q[2];
rz(-0.055214971) q[2];
rz(-0.84290543) q[3];
sx q[3];
rz(-2.3838145) q[3];
sx q[3];
rz(2.260476) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.371599) q[0];
sx q[0];
rz(-2.433625) q[0];
sx q[0];
rz(-1.1771033) q[0];
rz(-0.87989315) q[1];
sx q[1];
rz(-0.33029193) q[1];
sx q[1];
rz(-3.0807307) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0165981) q[0];
sx q[0];
rz(-1.6654764) q[0];
sx q[0];
rz(1.4461317) q[0];
rz(-pi) q[1];
rz(0.481769) q[2];
sx q[2];
rz(-1.8561811) q[2];
sx q[2];
rz(0.18295675) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.9990998) q[1];
sx q[1];
rz(-1.4542034) q[1];
sx q[1];
rz(-0.57013843) q[1];
x q[2];
rz(1.8350321) q[3];
sx q[3];
rz(-1.4995534) q[3];
sx q[3];
rz(-0.83282214) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-3.0941144) q[2];
sx q[2];
rz(-2.5884509) q[2];
sx q[2];
rz(1.6806357) q[2];
rz(2.285752) q[3];
sx q[3];
rz(-0.97270054) q[3];
sx q[3];
rz(-0.87876764) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1020553) q[0];
sx q[0];
rz(-0.85298959) q[0];
sx q[0];
rz(0.0041740388) q[0];
rz(-1.8904842) q[1];
sx q[1];
rz(-0.1336385) q[1];
sx q[1];
rz(-0.68152308) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.96723667) q[0];
sx q[0];
rz(-0.78034329) q[0];
sx q[0];
rz(-1.1514949) q[0];
rz(-pi) q[1];
rz(-0.84284346) q[2];
sx q[2];
rz(-2.3327391) q[2];
sx q[2];
rz(0.12241546) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.3423667) q[1];
sx q[1];
rz(-1.384642) q[1];
sx q[1];
rz(-0.28719814) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2596016) q[3];
sx q[3];
rz(-2.3025945) q[3];
sx q[3];
rz(-1.9669718) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.556813) q[2];
sx q[2];
rz(-0.4759554) q[2];
sx q[2];
rz(-1.7468096) q[2];
rz(0.41633385) q[3];
sx q[3];
rz(-1.8463912) q[3];
sx q[3];
rz(2.7487315) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3861179) q[0];
sx q[0];
rz(-0.32805726) q[0];
sx q[0];
rz(2.1395444) q[0];
rz(2.877032) q[1];
sx q[1];
rz(-1.8160276) q[1];
sx q[1];
rz(-1.4904259) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68448418) q[0];
sx q[0];
rz(-0.55782813) q[0];
sx q[0];
rz(-0.24866636) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4117496) q[2];
sx q[2];
rz(-2.3727086) q[2];
sx q[2];
rz(1.8184219) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.7539053) q[1];
sx q[1];
rz(-1.6847982) q[1];
sx q[1];
rz(-3.0579159) q[1];
x q[2];
rz(-0.57820676) q[3];
sx q[3];
rz(-0.1096519) q[3];
sx q[3];
rz(-2.6771817) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(3.1313021) q[2];
sx q[2];
rz(-1.4069858) q[2];
sx q[2];
rz(1.7245801) q[2];
rz(-0.33306444) q[3];
sx q[3];
rz(-0.76996961) q[3];
sx q[3];
rz(1.7634348) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3604816) q[0];
sx q[0];
rz(-1.8753373) q[0];
sx q[0];
rz(-1.2486096) q[0];
rz(3.1151415) q[1];
sx q[1];
rz(-1.5298264) q[1];
sx q[1];
rz(-1.5419921) q[1];
rz(-0.47172038) q[2];
sx q[2];
rz(-0.77009554) q[2];
sx q[2];
rz(-0.64252616) q[2];
rz(-0.61458434) q[3];
sx q[3];
rz(-1.3086223) q[3];
sx q[3];
rz(2.3400459) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
