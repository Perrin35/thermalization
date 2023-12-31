OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.6991601) q[0];
sx q[0];
rz(-1.7572829) q[0];
sx q[0];
rz(1.260489) q[0];
rz(2.1029544) q[1];
sx q[1];
rz(-1.3488052) q[1];
sx q[1];
rz(0.92372149) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7605654) q[0];
sx q[0];
rz(-1.1228704) q[0];
sx q[0];
rz(-0.3737803) q[0];
rz(-1.2878296) q[2];
sx q[2];
rz(-1.4322865) q[2];
sx q[2];
rz(-1.1793009) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-3.0337861) q[1];
sx q[1];
rz(-0.68192712) q[1];
sx q[1];
rz(-0.71180196) q[1];
x q[2];
rz(2.5039623) q[3];
sx q[3];
rz(-0.59083592) q[3];
sx q[3];
rz(1.6319815) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.91360056) q[2];
sx q[2];
rz(-1.8929409) q[2];
sx q[2];
rz(-0.16201924) q[2];
rz(-2.2062733) q[3];
sx q[3];
rz(-0.98615065) q[3];
sx q[3];
rz(-2.4285765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0682003) q[0];
sx q[0];
rz(-2.91495) q[0];
sx q[0];
rz(-1.9447928) q[0];
rz(2.4616922) q[1];
sx q[1];
rz(-2.6459243) q[1];
sx q[1];
rz(-1.686036) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6360977) q[0];
sx q[0];
rz(-2.539145) q[0];
sx q[0];
rz(-1.2732182) q[0];
rz(-0.20632867) q[2];
sx q[2];
rz(-0.38197877) q[2];
sx q[2];
rz(-1.0155592) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.1728954) q[1];
sx q[1];
rz(-2.6186133) q[1];
sx q[1];
rz(-1.5437267) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0807642) q[3];
sx q[3];
rz(-2.7292477) q[3];
sx q[3];
rz(-2.6044248) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.42852795) q[2];
sx q[2];
rz(-1.45168) q[2];
sx q[2];
rz(-1.7896174) q[2];
rz(-2.9591566) q[3];
sx q[3];
rz(-2.1648516) q[3];
sx q[3];
rz(0.3119719) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3882554) q[0];
sx q[0];
rz(-0.68080807) q[0];
sx q[0];
rz(0.80048168) q[0];
rz(-0.02877409) q[1];
sx q[1];
rz(-2.0859699) q[1];
sx q[1];
rz(1.9690537) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.81329936) q[0];
sx q[0];
rz(-1.0070224) q[0];
sx q[0];
rz(-1.203712) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3861888) q[2];
sx q[2];
rz(-1.8131776) q[2];
sx q[2];
rz(-2.541045) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.064627083) q[1];
sx q[1];
rz(-2.5588227) q[1];
sx q[1];
rz(-1.1175734) q[1];
rz(0.058733744) q[3];
sx q[3];
rz(-1.8211094) q[3];
sx q[3];
rz(2.5854923) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.0744434) q[2];
sx q[2];
rz(-1.477244) q[2];
sx q[2];
rz(-0.91119901) q[2];
rz(-0.95101142) q[3];
sx q[3];
rz(-0.8042897) q[3];
sx q[3];
rz(-2.2495911) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7610385) q[0];
sx q[0];
rz(-0.13042139) q[0];
sx q[0];
rz(-0.12810853) q[0];
rz(0.076106636) q[1];
sx q[1];
rz(-1.2144621) q[1];
sx q[1];
rz(-2.6180843) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2142221) q[0];
sx q[0];
rz(-2.4282051) q[0];
sx q[0];
rz(0.58332304) q[0];
rz(-1.6558311) q[2];
sx q[2];
rz(-1.2418613) q[2];
sx q[2];
rz(-2.0563682) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.88672968) q[1];
sx q[1];
rz(-0.88725315) q[1];
sx q[1];
rz(2.8121594) q[1];
rz(-pi) q[2];
rz(-0.97949667) q[3];
sx q[3];
rz(-1.6846091) q[3];
sx q[3];
rz(0.98609656) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.6161502) q[2];
sx q[2];
rz(-1.5972861) q[2];
sx q[2];
rz(-2.5775487) q[2];
rz(-0.28856746) q[3];
sx q[3];
rz(-2.7189062) q[3];
sx q[3];
rz(2.585876) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-0.48150912) q[0];
sx q[0];
rz(-2.4531589) q[0];
sx q[0];
rz(-1.4915285) q[0];
rz(2.2619757) q[1];
sx q[1];
rz(-1.2477701) q[1];
sx q[1];
rz(0.99194828) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.49682) q[0];
sx q[0];
rz(-2.9368375) q[0];
sx q[0];
rz(0.55069189) q[0];
x q[1];
rz(-1.5422103) q[2];
sx q[2];
rz(-0.30297849) q[2];
sx q[2];
rz(2.6945393) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.93047749) q[1];
sx q[1];
rz(-2.8301297) q[1];
sx q[1];
rz(0.13179563) q[1];
rz(3.0117412) q[3];
sx q[3];
rz(-0.81749812) q[3];
sx q[3];
rz(-2.0714456) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.0118959) q[2];
sx q[2];
rz(-2.7719438) q[2];
sx q[2];
rz(2.8707855) q[2];
rz(2.9233542) q[3];
sx q[3];
rz(-1.3202347) q[3];
sx q[3];
rz(0.22578421) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4145684) q[0];
sx q[0];
rz(-0.71838656) q[0];
sx q[0];
rz(-1.3487934) q[0];
rz(-2.7596966) q[1];
sx q[1];
rz(-2.8254639) q[1];
sx q[1];
rz(-1.4250925) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2857367) q[0];
sx q[0];
rz(-0.55186134) q[0];
sx q[0];
rz(-2.7049271) q[0];
rz(-pi) q[1];
x q[1];
rz(1.233333) q[2];
sx q[2];
rz(-2.87185) q[2];
sx q[2];
rz(-2.595682) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.7879494) q[1];
sx q[1];
rz(-1.1218346) q[1];
sx q[1];
rz(-3.06762) q[1];
x q[2];
rz(0.53374966) q[3];
sx q[3];
rz(-2.1042049) q[3];
sx q[3];
rz(0.48983869) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.0075334) q[2];
sx q[2];
rz(-0.39742658) q[2];
sx q[2];
rz(2.5777204) q[2];
rz(-0.18051906) q[3];
sx q[3];
rz(-1.518395) q[3];
sx q[3];
rz(-2.738651) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6329704) q[0];
sx q[0];
rz(-2.9794725) q[0];
sx q[0];
rz(2.7222743) q[0];
rz(1.5527027) q[1];
sx q[1];
rz(-1.2607375) q[1];
sx q[1];
rz(-2.3197876) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9200631) q[0];
sx q[0];
rz(-0.95525817) q[0];
sx q[0];
rz(-0.91381844) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6452438) q[2];
sx q[2];
rz(-0.90663547) q[2];
sx q[2];
rz(-1.6830483) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.070548363) q[1];
sx q[1];
rz(-0.76862915) q[1];
sx q[1];
rz(0.32960906) q[1];
x q[2];
rz(2.5542198) q[3];
sx q[3];
rz(-1.8837187) q[3];
sx q[3];
rz(1.3164933) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.3372779) q[2];
sx q[2];
rz(-0.75411212) q[2];
sx q[2];
rz(-2.896893) q[2];
rz(-3.0120567) q[3];
sx q[3];
rz(-1.9774388) q[3];
sx q[3];
rz(1.5130419) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.725175) q[0];
sx q[0];
rz(-0.019151909) q[0];
sx q[0];
rz(-2.3186671) q[0];
rz(-2.8322463) q[1];
sx q[1];
rz(-1.3920709) q[1];
sx q[1];
rz(-1.8364505) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37874052) q[0];
sx q[0];
rz(-1.2438602) q[0];
sx q[0];
rz(-2.5255894) q[0];
rz(0.42758503) q[2];
sx q[2];
rz(-2.3901849) q[2];
sx q[2];
rz(2.4497355) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.0481938) q[1];
sx q[1];
rz(-1.9462002) q[1];
sx q[1];
rz(2.4673389) q[1];
rz(-pi) q[2];
rz(-3.1180624) q[3];
sx q[3];
rz(-1.9100034) q[3];
sx q[3];
rz(-0.48285218) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.1398853) q[2];
sx q[2];
rz(-1.3858162) q[2];
sx q[2];
rz(1.4902327) q[2];
rz(1.0772609) q[3];
sx q[3];
rz(-2.1765985) q[3];
sx q[3];
rz(3.0100477) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7548783) q[0];
sx q[0];
rz(-1.2972378) q[0];
sx q[0];
rz(2.819678) q[0];
rz(-1.6053258) q[1];
sx q[1];
rz(-1.9202817) q[1];
sx q[1];
rz(2.4386491) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.87899938) q[0];
sx q[0];
rz(-1.2125373) q[0];
sx q[0];
rz(0.4155638) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0800955) q[2];
sx q[2];
rz(-1.5361538) q[2];
sx q[2];
rz(-2.408037) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.326509) q[1];
sx q[1];
rz(-1.0352967) q[1];
sx q[1];
rz(2.9498847) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1861107) q[3];
sx q[3];
rz(-1.8468879) q[3];
sx q[3];
rz(2.626782) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.90074173) q[2];
sx q[2];
rz(-2.8618331) q[2];
sx q[2];
rz(1.3396324) q[2];
rz(0.30570269) q[3];
sx q[3];
rz(-1.327508) q[3];
sx q[3];
rz(1.3302749) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.16383485) q[0];
sx q[0];
rz(-0.75755388) q[0];
sx q[0];
rz(-1.2257858) q[0];
rz(0.90351358) q[1];
sx q[1];
rz(-0.61360306) q[1];
sx q[1];
rz(0.46863619) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7077431) q[0];
sx q[0];
rz(-1.9540457) q[0];
sx q[0];
rz(-1.7953403) q[0];
x q[1];
rz(-2.2655728) q[2];
sx q[2];
rz(-2.6555736) q[2];
sx q[2];
rz(1.8234058) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.66099) q[1];
sx q[1];
rz(-0.27462474) q[1];
sx q[1];
rz(-1.8250699) q[1];
x q[2];
rz(-2.0300794) q[3];
sx q[3];
rz(-1.9643524) q[3];
sx q[3];
rz(0.72898385) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.48352155) q[2];
sx q[2];
rz(-1.2926241) q[2];
sx q[2];
rz(-1.998385) q[2];
rz(0.11463595) q[3];
sx q[3];
rz(-2.1879523) q[3];
sx q[3];
rz(1.6121929) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6464012) q[0];
sx q[0];
rz(-1.1947182) q[0];
sx q[0];
rz(2.4583046) q[0];
rz(-0.62190965) q[1];
sx q[1];
rz(-1.6786631) q[1];
sx q[1];
rz(2.8181029) q[1];
rz(0.8016349) q[2];
sx q[2];
rz(-0.87052204) q[2];
sx q[2];
rz(-1.082765) q[2];
rz(3.0631089) q[3];
sx q[3];
rz(-0.92072903) q[3];
sx q[3];
rz(-2.846684) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
