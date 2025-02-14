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
rz(3.003886) q[0];
sx q[0];
rz(-2.6258111) q[0];
sx q[0];
rz(-0.04096026) q[0];
rz(-2.3907258) q[1];
sx q[1];
rz(-2.6584396) q[1];
sx q[1];
rz(2.9274489) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.53759585) q[0];
sx q[0];
rz(-1.6600837) q[0];
sx q[0];
rz(0.16601899) q[0];
rz(1.7688355) q[2];
sx q[2];
rz(-1.4841586) q[2];
sx q[2];
rz(-1.600304) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.5109066) q[1];
sx q[1];
rz(-1.7244461) q[1];
sx q[1];
rz(-0.84801482) q[1];
rz(-pi) q[2];
rz(-0.80506206) q[3];
sx q[3];
rz(-0.65316155) q[3];
sx q[3];
rz(-1.4555318) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.21138771) q[2];
sx q[2];
rz(-1.3896421) q[2];
sx q[2];
rz(0.99838057) q[2];
rz(-1.8006648) q[3];
sx q[3];
rz(-1.0705592) q[3];
sx q[3];
rz(-2.4836331) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9407161) q[0];
sx q[0];
rz(-2.2287892) q[0];
sx q[0];
rz(0.87769133) q[0];
rz(0.25582036) q[1];
sx q[1];
rz(-2.1820549) q[1];
sx q[1];
rz(2.6928601) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.259346) q[0];
sx q[0];
rz(-1.2417698) q[0];
sx q[0];
rz(1.6297608) q[0];
rz(-3.0340292) q[2];
sx q[2];
rz(-1.7834181) q[2];
sx q[2];
rz(-0.68356252) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.4325368) q[1];
sx q[1];
rz(-0.87608713) q[1];
sx q[1];
rz(-1.1251775) q[1];
rz(-pi) q[2];
x q[2];
rz(0.99440672) q[3];
sx q[3];
rz(-2.4197516) q[3];
sx q[3];
rz(0.54366684) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.9933219) q[2];
sx q[2];
rz(-1.693087) q[2];
sx q[2];
rz(-1.5427962) q[2];
rz(-2.350542) q[3];
sx q[3];
rz(-1.6796716) q[3];
sx q[3];
rz(-0.79234523) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.032884447) q[0];
sx q[0];
rz(-2.5486163) q[0];
sx q[0];
rz(-1.2805043) q[0];
rz(2.9325824) q[1];
sx q[1];
rz(-2.0031395) q[1];
sx q[1];
rz(-1.4271522) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66433231) q[0];
sx q[0];
rz(-2.0149079) q[0];
sx q[0];
rz(0.61927619) q[0];
x q[1];
rz(1.4032321) q[2];
sx q[2];
rz(-1.6779876) q[2];
sx q[2];
rz(-0.78533781) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.64066891) q[1];
sx q[1];
rz(-0.66447778) q[1];
sx q[1];
rz(-2.6946103) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9455037) q[3];
sx q[3];
rz(-0.59299378) q[3];
sx q[3];
rz(-0.65341572) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.7872494) q[2];
sx q[2];
rz(-2.3086583) q[2];
sx q[2];
rz(2.1429817) q[2];
rz(0.35501114) q[3];
sx q[3];
rz(-1.7576926) q[3];
sx q[3];
rz(1.2242873) q[3];
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
sx q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3806279) q[0];
sx q[0];
rz(-2.9686718) q[0];
sx q[0];
rz(0.91249102) q[0];
rz(1.6418705) q[1];
sx q[1];
rz(-1.8955756) q[1];
sx q[1];
rz(-1.6778827) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.283962) q[0];
sx q[0];
rz(-1.2799834) q[0];
sx q[0];
rz(1.6785851) q[0];
rz(-pi) q[1];
rz(-1.6557515) q[2];
sx q[2];
rz(-0.96829295) q[2];
sx q[2];
rz(2.9311644) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.1596502) q[1];
sx q[1];
rz(-1.8022418) q[1];
sx q[1];
rz(0.39868928) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2578418) q[3];
sx q[3];
rz(-0.74713641) q[3];
sx q[3];
rz(-1.136387) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.0616167) q[2];
sx q[2];
rz(-1.7023664) q[2];
sx q[2];
rz(1.7495135) q[2];
rz(-2.020284) q[3];
sx q[3];
rz(-2.0544572) q[3];
sx q[3];
rz(-0.054718941) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1295117) q[0];
sx q[0];
rz(-1.5665781) q[0];
sx q[0];
rz(-0.47796252) q[0];
rz(-1.4683918) q[1];
sx q[1];
rz(-1.5696399) q[1];
sx q[1];
rz(1.063331) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37603518) q[0];
sx q[0];
rz(-2.4004694) q[0];
sx q[0];
rz(0.82668178) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9929286) q[2];
sx q[2];
rz(-2.1227266) q[2];
sx q[2];
rz(-0.27503289) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.0374839) q[1];
sx q[1];
rz(-1.9870305) q[1];
sx q[1];
rz(-2.3690696) q[1];
x q[2];
rz(1.8724779) q[3];
sx q[3];
rz(-1.2815803) q[3];
sx q[3];
rz(2.7986643) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.8308476) q[2];
sx q[2];
rz(-1.9644535) q[2];
sx q[2];
rz(1.8929405) q[2];
rz(-0.8904852) q[3];
sx q[3];
rz(-1.9200446) q[3];
sx q[3];
rz(0.42731592) q[3];
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
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2763727) q[0];
sx q[0];
rz(-0.38580147) q[0];
sx q[0];
rz(1.9355829) q[0];
rz(1.457006) q[1];
sx q[1];
rz(-2.147069) q[1];
sx q[1];
rz(1.0135894) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.550279) q[0];
sx q[0];
rz(-0.065803615) q[0];
sx q[0];
rz(-2.4273511) q[0];
rz(2.0034143) q[2];
sx q[2];
rz(-1.0241531) q[2];
sx q[2];
rz(-2.0631323) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.3675843) q[1];
sx q[1];
rz(-2.1220397) q[1];
sx q[1];
rz(1.8249976) q[1];
x q[2];
rz(-1.417755) q[3];
sx q[3];
rz(-0.42649999) q[3];
sx q[3];
rz(2.7271326) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.69063416) q[2];
sx q[2];
rz(-2.3834507) q[2];
sx q[2];
rz(3.0699733) q[2];
rz(0.13876638) q[3];
sx q[3];
rz(-1.7641188) q[3];
sx q[3];
rz(0.57851401) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11107681) q[0];
sx q[0];
rz(-1.0826305) q[0];
sx q[0];
rz(-1.6336596) q[0];
rz(1.1208447) q[1];
sx q[1];
rz(-1.0034674) q[1];
sx q[1];
rz(2.9337163) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5198398) q[0];
sx q[0];
rz(-0.43448453) q[0];
sx q[0];
rz(-2.1293917) q[0];
x q[1];
rz(-0.25983622) q[2];
sx q[2];
rz(-1.4186315) q[2];
sx q[2];
rz(2.2155188) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.4411197) q[1];
sx q[1];
rz(-1.9090096) q[1];
sx q[1];
rz(0.21041056) q[1];
rz(-1.3399501) q[3];
sx q[3];
rz(-1.5192598) q[3];
sx q[3];
rz(-1.7973943) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.8230744) q[2];
sx q[2];
rz(-1.069331) q[2];
sx q[2];
rz(-11/(7*pi)) q[2];
rz(-2.9253166) q[3];
sx q[3];
rz(-1.6612771) q[3];
sx q[3];
rz(1.1094619) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.08865393) q[0];
sx q[0];
rz(-1.5733938) q[0];
sx q[0];
rz(1.9807568) q[0];
rz(2.898518) q[1];
sx q[1];
rz(-1.4832152) q[1];
sx q[1];
rz(-1.9879139) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6223654) q[0];
sx q[0];
rz(-1.7891684) q[0];
sx q[0];
rz(1.9879935) q[0];
x q[1];
rz(-2.6188208) q[2];
sx q[2];
rz(-0.47975329) q[2];
sx q[2];
rz(0.72296491) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.8550314) q[1];
sx q[1];
rz(-1.5530545) q[1];
sx q[1];
rz(1.6231389) q[1];
rz(1.6457993) q[3];
sx q[3];
rz(-0.60216367) q[3];
sx q[3];
rz(-1.259089) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.2005177) q[2];
sx q[2];
rz(-0.87824559) q[2];
sx q[2];
rz(-0.78594691) q[2];
rz(-0.032622967) q[3];
sx q[3];
rz(-2.2226108) q[3];
sx q[3];
rz(-2.6526764) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0068479) q[0];
sx q[0];
rz(-2.5866046) q[0];
sx q[0];
rz(-0.81163374) q[0];
rz(-0.45400485) q[1];
sx q[1];
rz(-1.3521399) q[1];
sx q[1];
rz(-0.46617359) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3562856) q[0];
sx q[0];
rz(-1.6652096) q[0];
sx q[0];
rz(-2.055302) q[0];
x q[1];
rz(1.097686) q[2];
sx q[2];
rz(-1.1856831) q[2];
sx q[2];
rz(-1.3375741) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.55250231) q[1];
sx q[1];
rz(-1.9297761) q[1];
sx q[1];
rz(-0.41429452) q[1];
rz(-pi) q[2];
rz(-2.2884263) q[3];
sx q[3];
rz(-2.2090342) q[3];
sx q[3];
rz(-0.555942) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0635282) q[2];
sx q[2];
rz(-1.9643276) q[2];
sx q[2];
rz(1.9130116) q[2];
rz(2.581253) q[3];
sx q[3];
rz(-1.8599583) q[3];
sx q[3];
rz(1.3250215) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.233577) q[0];
sx q[0];
rz(-0.19291872) q[0];
sx q[0];
rz(2.8237901) q[0];
rz(0.69215149) q[1];
sx q[1];
rz(-1.0617804) q[1];
sx q[1];
rz(0.94921509) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5957221) q[0];
sx q[0];
rz(-0.4341653) q[0];
sx q[0];
rz(1.6950785) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8525392) q[2];
sx q[2];
rz(-2.2069283) q[2];
sx q[2];
rz(-0.67109443) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.91864955) q[1];
sx q[1];
rz(-2.6368196) q[1];
sx q[1];
rz(0.10279067) q[1];
x q[2];
rz(1.0383706) q[3];
sx q[3];
rz(-3.0255613) q[3];
sx q[3];
rz(0.86933245) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-3.074241) q[2];
sx q[2];
rz(-0.94506741) q[2];
sx q[2];
rz(2.2828339) q[2];
rz(2.7680715) q[3];
sx q[3];
rz(-1.7645323) q[3];
sx q[3];
rz(-0.080373272) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9141948) q[0];
sx q[0];
rz(-0.62794958) q[0];
sx q[0];
rz(-3.0218883) q[0];
rz(1.1788728) q[1];
sx q[1];
rz(-0.96794712) q[1];
sx q[1];
rz(1.6603574) q[1];
rz(-2.1356077) q[2];
sx q[2];
rz(-1.0848252) q[2];
sx q[2];
rz(1.8262524) q[2];
rz(-1.4734442) q[3];
sx q[3];
rz(-0.90383263) q[3];
sx q[3];
rz(-2.2624349) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
