OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.8060057) q[0];
sx q[0];
rz(-0.94526362) q[0];
sx q[0];
rz(2.6160016) q[0];
rz(0.2431915) q[1];
sx q[1];
rz(-1.9089729) q[1];
sx q[1];
rz(0.90484172) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5072767) q[0];
sx q[0];
rz(-1.1391369) q[0];
sx q[0];
rz(-1.0213724) q[0];
rz(-pi) q[1];
rz(2.4279847) q[2];
sx q[2];
rz(-0.2954233) q[2];
sx q[2];
rz(1.3908536) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.875784) q[1];
sx q[1];
rz(-1.335653) q[1];
sx q[1];
rz(-1.1633412) q[1];
x q[2];
rz(1.0337898) q[3];
sx q[3];
rz(-2.783943) q[3];
sx q[3];
rz(1.3057749) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.20377542) q[2];
sx q[2];
rz(-1.4062466) q[2];
sx q[2];
rz(-0.096244372) q[2];
rz(1.0359267) q[3];
sx q[3];
rz(-2.7544498) q[3];
sx q[3];
rz(0.15371418) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44089833) q[0];
sx q[0];
rz(-0.39114025) q[0];
sx q[0];
rz(-0.76517117) q[0];
rz(-1.8493429) q[1];
sx q[1];
rz(-2.6563829) q[1];
sx q[1];
rz(-0.66295019) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0239319) q[0];
sx q[0];
rz(-1.6370602) q[0];
sx q[0];
rz(-3.0793385) q[0];
rz(-pi) q[1];
x q[1];
rz(0.40132482) q[2];
sx q[2];
rz(-2.1839645) q[2];
sx q[2];
rz(-2.3344628) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.7141124) q[1];
sx q[1];
rz(-2.6926827) q[1];
sx q[1];
rz(-2.950579) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6547673) q[3];
sx q[3];
rz(-1.4273705) q[3];
sx q[3];
rz(-1.5407345) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-3.092676) q[2];
sx q[2];
rz(-2.0930223) q[2];
sx q[2];
rz(-2.8125787) q[2];
rz(2.4760903) q[3];
sx q[3];
rz(-0.21829675) q[3];
sx q[3];
rz(-1.3177419) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5927785) q[0];
sx q[0];
rz(-2.9187262) q[0];
sx q[0];
rz(-2.9192525) q[0];
rz(-1.0173343) q[1];
sx q[1];
rz(-2.4203114) q[1];
sx q[1];
rz(-0.51868784) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0761622) q[0];
sx q[0];
rz(-0.87834529) q[0];
sx q[0];
rz(2.6718219) q[0];
x q[1];
rz(1.6665886) q[2];
sx q[2];
rz(-1.3057858) q[2];
sx q[2];
rz(1.4694627) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.31836244) q[1];
sx q[1];
rz(-1.3981817) q[1];
sx q[1];
rz(2.361972) q[1];
x q[2];
rz(3.1268901) q[3];
sx q[3];
rz(-3.0869752) q[3];
sx q[3];
rz(-0.86044776) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.8609994) q[2];
sx q[2];
rz(-2.7797647) q[2];
sx q[2];
rz(3.0730491) q[2];
rz(-0.60244256) q[3];
sx q[3];
rz(-0.76255637) q[3];
sx q[3];
rz(0.13901916) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4908726) q[0];
sx q[0];
rz(-2.3644709) q[0];
sx q[0];
rz(2.9673476) q[0];
rz(-2.6113367) q[1];
sx q[1];
rz(-1.4825876) q[1];
sx q[1];
rz(-0.51309103) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1888694) q[0];
sx q[0];
rz(-0.9523305) q[0];
sx q[0];
rz(1.4737211) q[0];
x q[1];
rz(2.1999173) q[2];
sx q[2];
rz(-2.2721014) q[2];
sx q[2];
rz(-0.47052449) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.2345703) q[1];
sx q[1];
rz(-1.7090194) q[1];
sx q[1];
rz(0.45781086) q[1];
rz(1.3644049) q[3];
sx q[3];
rz(-1.6009814) q[3];
sx q[3];
rz(-2.9914732) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.6667368) q[2];
sx q[2];
rz(-2.7754144) q[2];
sx q[2];
rz(0.22988698) q[2];
rz(-0.41904467) q[3];
sx q[3];
rz(-1.34904) q[3];
sx q[3];
rz(-0.45927799) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.6722365) q[0];
sx q[0];
rz(-0.75119632) q[0];
sx q[0];
rz(-0.67681926) q[0];
rz(-2.6485486) q[1];
sx q[1];
rz(-0.9489916) q[1];
sx q[1];
rz(-0.61606032) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4426684) q[0];
sx q[0];
rz(-1.6232423) q[0];
sx q[0];
rz(3.1057182) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7763543) q[2];
sx q[2];
rz(-1.4589981) q[2];
sx q[2];
rz(-2.0888084) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.9628145) q[1];
sx q[1];
rz(-1.5652834) q[1];
sx q[1];
rz(2.5639736) q[1];
x q[2];
rz(0.89415278) q[3];
sx q[3];
rz(-1.3291306) q[3];
sx q[3];
rz(1.579293) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.4074576) q[2];
sx q[2];
rz(-2.5871758) q[2];
sx q[2];
rz(0.25536728) q[2];
rz(-1.6051965) q[3];
sx q[3];
rz(-2.1891749) q[3];
sx q[3];
rz(-0.77409625) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7191294) q[0];
sx q[0];
rz(-2.1753949) q[0];
sx q[0];
rz(0.47250026) q[0];
rz(-2.6155112) q[1];
sx q[1];
rz(-2.9317347) q[1];
sx q[1];
rz(-0.88476673) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5621983) q[0];
sx q[0];
rz(-1.644855) q[0];
sx q[0];
rz(1.9522569) q[0];
x q[1];
rz(-3.0753291) q[2];
sx q[2];
rz(-2.931086) q[2];
sx q[2];
rz(1.9312242) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.5342418) q[1];
sx q[1];
rz(-2.631175) q[1];
sx q[1];
rz(2.1900602) q[1];
rz(2.1045477) q[3];
sx q[3];
rz(-1.783672) q[3];
sx q[3];
rz(-2.6900929) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.3970967) q[2];
sx q[2];
rz(-2.3366191) q[2];
sx q[2];
rz(-2.8302144) q[2];
rz(1.3686251) q[3];
sx q[3];
rz(-2.6840648) q[3];
sx q[3];
rz(-0.5293203) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.41806528) q[0];
sx q[0];
rz(-2.8630246) q[0];
sx q[0];
rz(-3.080522) q[0];
rz(0.04018499) q[1];
sx q[1];
rz(-1.9804852) q[1];
sx q[1];
rz(0.73289245) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7948579) q[0];
sx q[0];
rz(-1.1746527) q[0];
sx q[0];
rz(0.57647716) q[0];
x q[1];
rz(1.0211208) q[2];
sx q[2];
rz(-0.78133821) q[2];
sx q[2];
rz(2.6598425) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.9961173) q[1];
sx q[1];
rz(-0.62578326) q[1];
sx q[1];
rz(1.5819509) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2039127) q[3];
sx q[3];
rz(-0.11493472) q[3];
sx q[3];
rz(0.68655187) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.0733033) q[2];
sx q[2];
rz(-1.1969593) q[2];
sx q[2];
rz(2.627009) q[2];
rz(-1.9040646) q[3];
sx q[3];
rz(-2.5585744) q[3];
sx q[3];
rz(-0.54491836) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
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
rz(0.056203689) q[0];
sx q[0];
rz(-1.6563002) q[0];
sx q[0];
rz(-2.432166) q[0];
rz(-1.6363232) q[1];
sx q[1];
rz(-2.067833) q[1];
sx q[1];
rz(-2.8628796) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7473135) q[0];
sx q[0];
rz(-0.25082591) q[0];
sx q[0];
rz(-2.5670811) q[0];
x q[1];
rz(0.3785554) q[2];
sx q[2];
rz(-2.5494529) q[2];
sx q[2];
rz(0.90422599) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.0540788) q[1];
sx q[1];
rz(-1.8808108) q[1];
sx q[1];
rz(2.9127496) q[1];
rz(1.2392427) q[3];
sx q[3];
rz(-2.9058876) q[3];
sx q[3];
rz(1.9086259) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.8401106) q[2];
sx q[2];
rz(-1.8436517) q[2];
sx q[2];
rz(3.0855132) q[2];
rz(2.2864443) q[3];
sx q[3];
rz(-0.44848281) q[3];
sx q[3];
rz(0.40518951) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.99496019) q[0];
sx q[0];
rz(-0.17630795) q[0];
sx q[0];
rz(0.96889281) q[0];
rz(0.47337198) q[1];
sx q[1];
rz(-2.3469766) q[1];
sx q[1];
rz(1.999058) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8965217) q[0];
sx q[0];
rz(-2.8327496) q[0];
sx q[0];
rz(1.2678498) q[0];
rz(3.1387024) q[2];
sx q[2];
rz(-2.4382466) q[2];
sx q[2];
rz(-3.0386915) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-3.0513623) q[1];
sx q[1];
rz(-1.6962595) q[1];
sx q[1];
rz(-2.5999703) q[1];
rz(-pi) q[2];
rz(-1.3167131) q[3];
sx q[3];
rz(-2.4554606) q[3];
sx q[3];
rz(0.38476598) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.2450976) q[2];
sx q[2];
rz(-1.2197887) q[2];
sx q[2];
rz(1.0207821) q[2];
rz(-2.8178689) q[3];
sx q[3];
rz(-0.75298572) q[3];
sx q[3];
rz(-0.011172115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1616515) q[0];
sx q[0];
rz(-0.027898235) q[0];
sx q[0];
rz(0.7014057) q[0];
rz(0.91570634) q[1];
sx q[1];
rz(-1.0083102) q[1];
sx q[1];
rz(1.9030301) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.21047132) q[0];
sx q[0];
rz(-2.5744372) q[0];
sx q[0];
rz(1.9803067) q[0];
x q[1];
rz(2.3886016) q[2];
sx q[2];
rz(-1.5333311) q[2];
sx q[2];
rz(2.3245036) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.30675948) q[1];
sx q[1];
rz(-0.66654897) q[1];
sx q[1];
rz(-1.761761) q[1];
rz(-0.978312) q[3];
sx q[3];
rz(-1.6654135) q[3];
sx q[3];
rz(0.30944165) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.68676585) q[2];
sx q[2];
rz(-2.0451615) q[2];
sx q[2];
rz(-3.0977541) q[2];
rz(1.94058) q[3];
sx q[3];
rz(-0.73533708) q[3];
sx q[3];
rz(-2.1380077) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.6702406) q[0];
sx q[0];
rz(-2.4155407) q[0];
sx q[0];
rz(1.7759905) q[0];
rz(-0.025370601) q[1];
sx q[1];
rz(-1.3062968) q[1];
sx q[1];
rz(1.2702373) q[1];
rz(1.17169) q[2];
sx q[2];
rz(-1.8532955) q[2];
sx q[2];
rz(2.4886139) q[2];
rz(2.3425441) q[3];
sx q[3];
rz(-1.4481164) q[3];
sx q[3];
rz(-1.8275402) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
