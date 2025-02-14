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
rz(3.0930003) q[0];
sx q[0];
rz(-1.4253923) q[0];
sx q[0];
rz(0.26560321) q[0];
rz(-0.44759294) q[1];
sx q[1];
rz(4.6456479) q[1];
sx q[1];
rz(9.3001443) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4701047) q[0];
sx q[0];
rz(-2.3186145) q[0];
sx q[0];
rz(-2.1822071) q[0];
x q[1];
rz(1.2536178) q[2];
sx q[2];
rz(-2.7493959) q[2];
sx q[2];
rz(-2.9756096) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.0693486) q[1];
sx q[1];
rz(-0.26990971) q[1];
sx q[1];
rz(-2.4000485) q[1];
rz(0.83104317) q[3];
sx q[3];
rz(-0.22749113) q[3];
sx q[3];
rz(-0.25615197) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.0894185) q[2];
sx q[2];
rz(-2.3201421) q[2];
sx q[2];
rz(2.2649412) q[2];
rz(-2.9540201) q[3];
sx q[3];
rz(-2.1552174) q[3];
sx q[3];
rz(-0.11876373) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.12980421) q[0];
sx q[0];
rz(-3.0846444) q[0];
sx q[0];
rz(-0.38527986) q[0];
rz(0.41703364) q[1];
sx q[1];
rz(-0.72154355) q[1];
sx q[1];
rz(-1.0122274) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4436662) q[0];
sx q[0];
rz(-1.3529643) q[0];
sx q[0];
rz(0.99866726) q[0];
rz(-2.5778722) q[2];
sx q[2];
rz(-2.4032058) q[2];
sx q[2];
rz(2.9132879) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.8292147) q[1];
sx q[1];
rz(-0.80655231) q[1];
sx q[1];
rz(-0.4017425) q[1];
rz(-pi) q[2];
rz(-0.46383922) q[3];
sx q[3];
rz(-0.94479221) q[3];
sx q[3];
rz(-0.30445619) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.9932844) q[2];
sx q[2];
rz(-2.7316284) q[2];
sx q[2];
rz(-2.6170464) q[2];
rz(-2.2927393) q[3];
sx q[3];
rz(-2.4436804) q[3];
sx q[3];
rz(-2.2822288) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9808905) q[0];
sx q[0];
rz(-2.7301259) q[0];
sx q[0];
rz(0.32660642) q[0];
rz(-1.3702565) q[1];
sx q[1];
rz(-1.0561008) q[1];
sx q[1];
rz(0.036570963) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9865788) q[0];
sx q[0];
rz(-2.554092) q[0];
sx q[0];
rz(-0.67836887) q[0];
rz(-0.49679773) q[2];
sx q[2];
rz(-0.76947509) q[2];
sx q[2];
rz(2.0747607) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.4160683) q[1];
sx q[1];
rz(-1.5276794) q[1];
sx q[1];
rz(1.3836224) q[1];
rz(-pi) q[2];
rz(-1.8518515) q[3];
sx q[3];
rz(-1.7450168) q[3];
sx q[3];
rz(-0.094467316) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.44564351) q[2];
sx q[2];
rz(-0.30569884) q[2];
sx q[2];
rz(-1.9574399) q[2];
rz(-0.46237692) q[3];
sx q[3];
rz(-1.0824883) q[3];
sx q[3];
rz(2.6760127) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9584123) q[0];
sx q[0];
rz(-1.4875702) q[0];
sx q[0];
rz(2.2271449) q[0];
rz(0.86638266) q[1];
sx q[1];
rz(-1.2976846) q[1];
sx q[1];
rz(2.8241209) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1007967) q[0];
sx q[0];
rz(-1.6799752) q[0];
sx q[0];
rz(-0.37495701) q[0];
rz(-pi) q[1];
rz(0.16902618) q[2];
sx q[2];
rz(-1.9155972) q[2];
sx q[2];
rz(-1.930869) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.9383512) q[1];
sx q[1];
rz(-2.8305191) q[1];
sx q[1];
rz(-2.0955906) q[1];
rz(1.5741979) q[3];
sx q[3];
rz(-2.7614584) q[3];
sx q[3];
rz(-1.8148607) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.094396599) q[2];
sx q[2];
rz(-2.6698343) q[2];
sx q[2];
rz(-0.45486927) q[2];
rz(1.1826078) q[3];
sx q[3];
rz(-0.22795658) q[3];
sx q[3];
rz(-3.0003149) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8126548) q[0];
sx q[0];
rz(-1.8998572) q[0];
sx q[0];
rz(-3.090233) q[0];
rz(-2.8142169) q[1];
sx q[1];
rz(-0.86762571) q[1];
sx q[1];
rz(-3.0816269) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.77207309) q[0];
sx q[0];
rz(-0.74893337) q[0];
sx q[0];
rz(-1.0706484) q[0];
x q[1];
rz(-0.48624235) q[2];
sx q[2];
rz(-0.71388968) q[2];
sx q[2];
rz(-1.3522268) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.9717945) q[1];
sx q[1];
rz(-0.64193875) q[1];
sx q[1];
rz(-2.1084014) q[1];
x q[2];
rz(-2.9465605) q[3];
sx q[3];
rz(-2.2761506) q[3];
sx q[3];
rz(0.56567276) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.7472234) q[2];
sx q[2];
rz(-0.9129492) q[2];
sx q[2];
rz(-0.24820122) q[2];
rz(-0.40211755) q[3];
sx q[3];
rz(-0.44200236) q[3];
sx q[3];
rz(2.2615652) q[3];
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
rz(-pi) q[0];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.0070655951) q[0];
sx q[0];
rz(-1.8777254) q[0];
sx q[0];
rz(-1.7279351) q[0];
rz(-1.2029348) q[1];
sx q[1];
rz(-2.2212432) q[1];
sx q[1];
rz(0.10341067) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0132842) q[0];
sx q[0];
rz(-0.09311267) q[0];
sx q[0];
rz(0.72269209) q[0];
rz(-2.8273583) q[2];
sx q[2];
rz(-1.7467919) q[2];
sx q[2];
rz(0.90547784) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.9050871) q[1];
sx q[1];
rz(-0.12612045) q[1];
sx q[1];
rz(-2.990475) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8887022) q[3];
sx q[3];
rz(-1.9441981) q[3];
sx q[3];
rz(-1.6643844) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.6695413) q[2];
sx q[2];
rz(-2.5504888) q[2];
sx q[2];
rz(2.9284787) q[2];
rz(1.6040246) q[3];
sx q[3];
rz(-3.0209916) q[3];
sx q[3];
rz(-0.11046256) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9977426) q[0];
sx q[0];
rz(-0.85875964) q[0];
sx q[0];
rz(0.48458883) q[0];
rz(0.46802014) q[1];
sx q[1];
rz(-1.3222398) q[1];
sx q[1];
rz(0.13672926) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3829684) q[0];
sx q[0];
rz(-0.45238414) q[0];
sx q[0];
rz(-2.1470619) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9152476) q[2];
sx q[2];
rz(-1.9392804) q[2];
sx q[2];
rz(-1.1158352) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.9715106) q[1];
sx q[1];
rz(-2.2822343) q[1];
sx q[1];
rz(0.17632874) q[1];
x q[2];
rz(2.7220821) q[3];
sx q[3];
rz(-1.4624573) q[3];
sx q[3];
rz(-0.51589032) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.7318657) q[2];
sx q[2];
rz(-0.12696433) q[2];
sx q[2];
rz(-0.48180386) q[2];
rz(1.2305772) q[3];
sx q[3];
rz(-1.1426208) q[3];
sx q[3];
rz(-0.13550152) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.19732533) q[0];
sx q[0];
rz(-0.33536401) q[0];
sx q[0];
rz(-2.770597) q[0];
rz(-0.15759298) q[1];
sx q[1];
rz(-1.7540365) q[1];
sx q[1];
rz(-0.88034672) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.77773422) q[0];
sx q[0];
rz(-0.4798747) q[0];
sx q[0];
rz(-1.4196447) q[0];
rz(2.9463459) q[2];
sx q[2];
rz(-0.71788997) q[2];
sx q[2];
rz(0.860329) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.62746) q[1];
sx q[1];
rz(-1.4480906) q[1];
sx q[1];
rz(3.1377327) q[1];
x q[2];
rz(0.28435376) q[3];
sx q[3];
rz(-1.2907249) q[3];
sx q[3];
rz(-1.3765311) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.2788435) q[2];
sx q[2];
rz(-2.4767488) q[2];
sx q[2];
rz(2.1319907) q[2];
rz(0.75262117) q[3];
sx q[3];
rz(-1.8372583) q[3];
sx q[3];
rz(-0.93839222) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2272334) q[0];
sx q[0];
rz(-3.0000913) q[0];
sx q[0];
rz(-3.09521) q[0];
rz(-1.6498238) q[1];
sx q[1];
rz(-1.7581419) q[1];
sx q[1];
rz(0.47328624) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2641122) q[0];
sx q[0];
rz(-0.22019596) q[0];
sx q[0];
rz(1.0408677) q[0];
rz(-pi) q[1];
rz(2.0436486) q[2];
sx q[2];
rz(-1.2769498) q[2];
sx q[2];
rz(-2.1515255) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.5492268) q[1];
sx q[1];
rz(-1.8344546) q[1];
sx q[1];
rz(0.63559611) q[1];
rz(-pi) q[2];
rz(-1.6026845) q[3];
sx q[3];
rz(-1.569935) q[3];
sx q[3];
rz(-0.79142785) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.16372323) q[2];
sx q[2];
rz(-0.59000868) q[2];
sx q[2];
rz(0.02656492) q[2];
rz(2.6350392) q[3];
sx q[3];
rz(-0.81304628) q[3];
sx q[3];
rz(-2.626239) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.16515054) q[0];
sx q[0];
rz(-0.9557752) q[0];
sx q[0];
rz(2.9859848) q[0];
rz(2.7071629) q[1];
sx q[1];
rz(-2.5948718) q[1];
sx q[1];
rz(2.853552) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.018943448) q[0];
sx q[0];
rz(-2.3068006) q[0];
sx q[0];
rz(0.74170603) q[0];
x q[1];
rz(-3.1278438) q[2];
sx q[2];
rz(-1.4176157) q[2];
sx q[2];
rz(1.8977697) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.1571785) q[1];
sx q[1];
rz(-2.0044185) q[1];
sx q[1];
rz(0.84624802) q[1];
rz(-pi) q[2];
x q[2];
rz(0.94567255) q[3];
sx q[3];
rz(-1.9145619) q[3];
sx q[3];
rz(-1.6812066) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.41409278) q[2];
sx q[2];
rz(-2.6850271) q[2];
sx q[2];
rz(0.045819316) q[2];
rz(-3.0231061) q[3];
sx q[3];
rz(-0.89209569) q[3];
sx q[3];
rz(-2.4011325) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4988149) q[0];
sx q[0];
rz(-1.6294263) q[0];
sx q[0];
rz(-1.908041) q[0];
rz(-1.9998101) q[1];
sx q[1];
rz(-0.70527609) q[1];
sx q[1];
rz(-1.3377778) q[1];
rz(0.98182044) q[2];
sx q[2];
rz(-2.1661988) q[2];
sx q[2];
rz(1.7924776) q[2];
rz(1.6330887) q[3];
sx q[3];
rz(-1.6332165) q[3];
sx q[3];
rz(-1.6467057) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
