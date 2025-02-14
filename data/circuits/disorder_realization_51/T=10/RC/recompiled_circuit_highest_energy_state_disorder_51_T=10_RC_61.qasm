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
rz(0.32938862) q[0];
sx q[0];
rz(-2.5100799) q[0];
sx q[0];
rz(0.00087498571) q[0];
rz(-0.64088351) q[1];
sx q[1];
rz(-1.0008608) q[1];
sx q[1];
rz(0.34520087) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76287718) q[0];
sx q[0];
rz(-1.595531) q[0];
sx q[0];
rz(1.6616761) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4232741) q[2];
sx q[2];
rz(-1.194724) q[2];
sx q[2];
rz(3.0126115) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.7172437) q[1];
sx q[1];
rz(-2.1645438) q[1];
sx q[1];
rz(2.6301066) q[1];
x q[2];
rz(-1.6451938) q[3];
sx q[3];
rz(-1.5870023) q[3];
sx q[3];
rz(-2.5287573) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.9125646) q[2];
sx q[2];
rz(-1.8077069) q[2];
sx q[2];
rz(-0.20797569) q[2];
rz(0.29911706) q[3];
sx q[3];
rz(-2.5695473) q[3];
sx q[3];
rz(-2.0007029) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3308554) q[0];
sx q[0];
rz(-2.9006697) q[0];
sx q[0];
rz(2.148707) q[0];
rz(1.8244686) q[1];
sx q[1];
rz(-0.31610745) q[1];
sx q[1];
rz(0.77450007) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3178783) q[0];
sx q[0];
rz(-1.5817989) q[0];
sx q[0];
rz(1.7490134) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6235145) q[2];
sx q[2];
rz(-0.87730125) q[2];
sx q[2];
rz(-0.77502807) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.6933813) q[1];
sx q[1];
rz(-0.44084545) q[1];
sx q[1];
rz(-0.54852672) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1011613) q[3];
sx q[3];
rz(-1.7723933) q[3];
sx q[3];
rz(0.93702173) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.896686) q[2];
sx q[2];
rz(-1.8550355) q[2];
sx q[2];
rz(-1.0150821) q[2];
rz(0.24909881) q[3];
sx q[3];
rz(-2.2738012) q[3];
sx q[3];
rz(-0.36756137) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.86431137) q[0];
sx q[0];
rz(-0.69121498) q[0];
sx q[0];
rz(2.9840898) q[0];
rz(2.1306254) q[1];
sx q[1];
rz(-0.33282655) q[1];
sx q[1];
rz(-3.0027622) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9897744) q[0];
sx q[0];
rz(-2.646137) q[0];
sx q[0];
rz(-0.83924967) q[0];
rz(-pi) q[1];
rz(-1.9857668) q[2];
sx q[2];
rz(-2.2425644) q[2];
sx q[2];
rz(-0.49723724) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.87698679) q[1];
sx q[1];
rz(-2.786676) q[1];
sx q[1];
rz(-0.97477977) q[1];
x q[2];
rz(-2.5957727) q[3];
sx q[3];
rz(-2.8132317) q[3];
sx q[3];
rz(-2.2203746) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.6802754) q[2];
sx q[2];
rz(-0.91247827) q[2];
sx q[2];
rz(2.7834564) q[2];
rz(-0.67874587) q[3];
sx q[3];
rz(-0.94921422) q[3];
sx q[3];
rz(2.1024735) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0963652) q[0];
sx q[0];
rz(-2.1684833) q[0];
sx q[0];
rz(-2.8367693) q[0];
rz(0.40799704) q[1];
sx q[1];
rz(-1.4381189) q[1];
sx q[1];
rz(-2.2136484) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1243195) q[0];
sx q[0];
rz(-1.7649179) q[0];
sx q[0];
rz(-3.0108736) q[0];
rz(-pi) q[1];
rz(-0.81176968) q[2];
sx q[2];
rz(-2.3583745) q[2];
sx q[2];
rz(-2.4350172) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.9548861) q[1];
sx q[1];
rz(-2.2518603) q[1];
sx q[1];
rz(-1.7658556) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4797012) q[3];
sx q[3];
rz(-2.3958979) q[3];
sx q[3];
rz(0.5058561) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.1912332) q[2];
sx q[2];
rz(-0.57098907) q[2];
sx q[2];
rz(-2.9534269) q[2];
rz(-1.865271) q[3];
sx q[3];
rz(-1.8264344) q[3];
sx q[3];
rz(1.7108542) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4417878) q[0];
sx q[0];
rz(-2.7975174) q[0];
sx q[0];
rz(0.019388327) q[0];
rz(-2.0806606) q[1];
sx q[1];
rz(-1.414199) q[1];
sx q[1];
rz(-1.9020938) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4864522) q[0];
sx q[0];
rz(-1.0024655) q[0];
sx q[0];
rz(0.33354946) q[0];
rz(2.2197228) q[2];
sx q[2];
rz(-2.1989294) q[2];
sx q[2];
rz(-2.3482196) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.6790948) q[1];
sx q[1];
rz(-2.0789008) q[1];
sx q[1];
rz(-1.8930045) q[1];
x q[2];
rz(0.89770384) q[3];
sx q[3];
rz(-0.60808676) q[3];
sx q[3];
rz(-1.2230108) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.0197319) q[2];
sx q[2];
rz(-1.3132361) q[2];
sx q[2];
rz(1.8149553) q[2];
rz(-2.895368) q[3];
sx q[3];
rz(-2.0387869) q[3];
sx q[3];
rz(-0.8518014) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5354079) q[0];
sx q[0];
rz(-2.1562205) q[0];
sx q[0];
rz(0.73915172) q[0];
rz(-2.7958561) q[1];
sx q[1];
rz(-1.3290936) q[1];
sx q[1];
rz(2.2703222) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5369878) q[0];
sx q[0];
rz(-2.3571627) q[0];
sx q[0];
rz(-3.0434199) q[0];
rz(0.61001444) q[2];
sx q[2];
rz(-0.84976174) q[2];
sx q[2];
rz(-2.339156) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.024684357) q[1];
sx q[1];
rz(-1.643812) q[1];
sx q[1];
rz(-1.1817929) q[1];
x q[2];
rz(-2.5302251) q[3];
sx q[3];
rz(-1.1725559) q[3];
sx q[3];
rz(0.37722019) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.31200108) q[2];
sx q[2];
rz(-0.63903725) q[2];
sx q[2];
rz(-0.87257067) q[2];
rz(-1.5729337) q[3];
sx q[3];
rz(-2.5376153) q[3];
sx q[3];
rz(0.93305552) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0912112) q[0];
sx q[0];
rz(-2.8822883) q[0];
sx q[0];
rz(-2.3799489) q[0];
rz(0.42539445) q[1];
sx q[1];
rz(-2.0014706) q[1];
sx q[1];
rz(-0.92686191) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8445209) q[0];
sx q[0];
rz(-1.8064587) q[0];
sx q[0];
rz(1.8051487) q[0];
rz(-pi) q[1];
x q[1];
rz(0.71304597) q[2];
sx q[2];
rz(-2.2879061) q[2];
sx q[2];
rz(0.46592679) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.95185223) q[1];
sx q[1];
rz(-1.9323255) q[1];
sx q[1];
rz(-2.9041163) q[1];
rz(-1.5038483) q[3];
sx q[3];
rz(-2.5173325) q[3];
sx q[3];
rz(-0.74893307) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.1130134) q[2];
sx q[2];
rz(-0.69602746) q[2];
sx q[2];
rz(-0.87116233) q[2];
rz(-2.8431559) q[3];
sx q[3];
rz(-1.8961743) q[3];
sx q[3];
rz(-1.3309853) q[3];
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
rz(-pi) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7262064) q[0];
sx q[0];
rz(-2.6415249) q[0];
sx q[0];
rz(0.4183847) q[0];
rz(2.8437974) q[1];
sx q[1];
rz(-1.9458385) q[1];
sx q[1];
rz(0.13430886) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90420049) q[0];
sx q[0];
rz(-0.89590329) q[0];
sx q[0];
rz(-0.86752059) q[0];
rz(-2.9555126) q[2];
sx q[2];
rz(-1.1717516) q[2];
sx q[2];
rz(-3.0732791) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.2063576) q[1];
sx q[1];
rz(-2.9138953) q[1];
sx q[1];
rz(-2.5050194) q[1];
rz(-pi) q[2];
x q[2];
rz(0.28388309) q[3];
sx q[3];
rz(-2.0238658) q[3];
sx q[3];
rz(2.4456152) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.7316651) q[2];
sx q[2];
rz(-1.8525367) q[2];
sx q[2];
rz(0.64201391) q[2];
rz(2.0783453) q[3];
sx q[3];
rz(-2.4898873) q[3];
sx q[3];
rz(1.0412019) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6006271) q[0];
sx q[0];
rz(-0.0890812) q[0];
sx q[0];
rz(0.11216057) q[0];
rz(1.8070096) q[1];
sx q[1];
rz(-2.5490675) q[1];
sx q[1];
rz(-1.0293915) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4698668) q[0];
sx q[0];
rz(-2.436536) q[0];
sx q[0];
rz(1.5211578) q[0];
rz(1.6104944) q[2];
sx q[2];
rz(-1.8744812) q[2];
sx q[2];
rz(0.025321753) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.7299906) q[1];
sx q[1];
rz(-1.8008324) q[1];
sx q[1];
rz(3.0964628) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4791895) q[3];
sx q[3];
rz(-2.0370968) q[3];
sx q[3];
rz(0.3084076) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.42334291) q[2];
sx q[2];
rz(-0.074904718) q[2];
sx q[2];
rz(0.038012803) q[2];
rz(0.612261) q[3];
sx q[3];
rz(-2.2050048) q[3];
sx q[3];
rz(-0.14122252) q[3];
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
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8925979) q[0];
sx q[0];
rz(-1.4709512) q[0];
sx q[0];
rz(-2.7227962) q[0];
rz(-1.2576125) q[1];
sx q[1];
rz(-1.6323171) q[1];
sx q[1];
rz(0.14258252) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.52209848) q[0];
sx q[0];
rz(-1.5428679) q[0];
sx q[0];
rz(2.9832207) q[0];
rz(-1.293574) q[2];
sx q[2];
rz(-2.1810594) q[2];
sx q[2];
rz(-1.1278314) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.8293162) q[1];
sx q[1];
rz(-1.8354776) q[1];
sx q[1];
rz(-2.1618202) q[1];
rz(-pi) q[2];
rz(-2.0311277) q[3];
sx q[3];
rz(-1.2602196) q[3];
sx q[3];
rz(0.47720695) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.7964145) q[2];
sx q[2];
rz(-3.0609481) q[2];
sx q[2];
rz(0.26819116) q[2];
rz(1.7372519) q[3];
sx q[3];
rz(-1.0920478) q[3];
sx q[3];
rz(-1.0442737) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0161229) q[0];
sx q[0];
rz(-2.7653427) q[0];
sx q[0];
rz(0.34633037) q[0];
rz(-0.12915962) q[1];
sx q[1];
rz(-1.8864514) q[1];
sx q[1];
rz(1.4056978) q[1];
rz(-0.62025537) q[2];
sx q[2];
rz(-2.5778985) q[2];
sx q[2];
rz(2.3621205) q[2];
rz(-1.2760194) q[3];
sx q[3];
rz(-1.634292) q[3];
sx q[3];
rz(2.5802142) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
