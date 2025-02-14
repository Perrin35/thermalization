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
rz(0.41843721) q[0];
sx q[0];
rz(-0.96324459) q[0];
sx q[0];
rz(0.20382398) q[0];
rz(1.0526429) q[1];
sx q[1];
rz(3.9349603) q[1];
sx q[1];
rz(10.959672) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2685854) q[0];
sx q[0];
rz(-2.632336) q[0];
sx q[0];
rz(1.0727706) q[0];
rz(-pi) q[1];
rz(-1.9594749) q[2];
sx q[2];
rz(-0.76180327) q[2];
sx q[2];
rz(-2.3488059) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.878256) q[1];
sx q[1];
rz(-1.9227826) q[1];
sx q[1];
rz(0.89603591) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9382557) q[3];
sx q[3];
rz(-1.7032188) q[3];
sx q[3];
rz(1.2156957) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.4172998) q[2];
sx q[2];
rz(-2.5827926) q[2];
sx q[2];
rz(2.144045) q[2];
rz(-2.3857462) q[3];
sx q[3];
rz(-1.7852802) q[3];
sx q[3];
rz(2.1006179) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.34870979) q[0];
sx q[0];
rz(-0.097276874) q[0];
sx q[0];
rz(1.8541699) q[0];
rz(0.48311326) q[1];
sx q[1];
rz(-2.099642) q[1];
sx q[1];
rz(1.2189254) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.78271237) q[0];
sx q[0];
rz(-0.25280372) q[0];
sx q[0];
rz(2.2532399) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8701209) q[2];
sx q[2];
rz(-0.69005943) q[2];
sx q[2];
rz(0.59218237) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.7658577) q[1];
sx q[1];
rz(-1.2925783) q[1];
sx q[1];
rz(0.15032676) q[1];
rz(-pi) q[2];
rz(-0.87872259) q[3];
sx q[3];
rz(-1.2936397) q[3];
sx q[3];
rz(-0.73776484) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.7634742) q[2];
sx q[2];
rz(-0.066630445) q[2];
sx q[2];
rz(-0.62180579) q[2];
rz(-2.5032737) q[3];
sx q[3];
rz(-2.392605) q[3];
sx q[3];
rz(2.2973255) q[3];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8353552) q[0];
sx q[0];
rz(-2.4995646) q[0];
sx q[0];
rz(-0.21632347) q[0];
rz(-0.59858876) q[1];
sx q[1];
rz(-1.3052156) q[1];
sx q[1];
rz(-2.6187706) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4958333) q[0];
sx q[0];
rz(-1.2649987) q[0];
sx q[0];
rz(-1.9093139) q[0];
rz(1.7274569) q[2];
sx q[2];
rz(-1.9576503) q[2];
sx q[2];
rz(1.1413107) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.26979687) q[1];
sx q[1];
rz(-1.3482454) q[1];
sx q[1];
rz(2.6128164) q[1];
rz(-0.82562311) q[3];
sx q[3];
rz(-0.43681991) q[3];
sx q[3];
rz(-2.2568011) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.181695) q[2];
sx q[2];
rz(-1.2744224) q[2];
sx q[2];
rz(1.6240906) q[2];
rz(2.6767139) q[3];
sx q[3];
rz(-0.93022323) q[3];
sx q[3];
rz(-1.6200861) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14438039) q[0];
sx q[0];
rz(-2.7535487) q[0];
sx q[0];
rz(-1.602518) q[0];
rz(2.1082361) q[1];
sx q[1];
rz(-2.3732503) q[1];
sx q[1];
rz(1.296952) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.715623) q[0];
sx q[0];
rz(-0.85943009) q[0];
sx q[0];
rz(-2.864896) q[0];
rz(-3.0791666) q[2];
sx q[2];
rz(-0.87597825) q[2];
sx q[2];
rz(0.42499396) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.0694643) q[1];
sx q[1];
rz(-1.9481244) q[1];
sx q[1];
rz(1.1262141) q[1];
x q[2];
rz(2.6565032) q[3];
sx q[3];
rz(-2.5311433) q[3];
sx q[3];
rz(-1.4097253) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.3492655) q[2];
sx q[2];
rz(-2.5040864) q[2];
sx q[2];
rz(-1.0456592) q[2];
rz(0.21444923) q[3];
sx q[3];
rz(-0.72434536) q[3];
sx q[3];
rz(1.1469871) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8056718) q[0];
sx q[0];
rz(-1.2821953) q[0];
sx q[0];
rz(-0.51934284) q[0];
rz(2.1995811) q[1];
sx q[1];
rz(-0.28128925) q[1];
sx q[1];
rz(1.5325783) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62270861) q[0];
sx q[0];
rz(-1.2075338) q[0];
sx q[0];
rz(0.69464442) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4759859) q[2];
sx q[2];
rz(-2.9133248) q[2];
sx q[2];
rz(-2.1537202) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.457442) q[1];
sx q[1];
rz(-1.5760211) q[1];
sx q[1];
rz(-1.0702151) q[1];
x q[2];
rz(-1.4761837) q[3];
sx q[3];
rz(-1.4615371) q[3];
sx q[3];
rz(3.008568) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.4970826) q[2];
sx q[2];
rz(-0.77748674) q[2];
sx q[2];
rz(-0.28437781) q[2];
rz(0.55822462) q[3];
sx q[3];
rz(-1.7773209) q[3];
sx q[3];
rz(-0.24793454) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
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
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.072642) q[0];
sx q[0];
rz(-0.41256368) q[0];
sx q[0];
rz(-0.64055842) q[0];
rz(-1.5156281) q[1];
sx q[1];
rz(-2.0104355) q[1];
sx q[1];
rz(-2.9277149) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.86588496) q[0];
sx q[0];
rz(-1.2684221) q[0];
sx q[0];
rz(-2.9024966) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3967314) q[2];
sx q[2];
rz(-1.8554885) q[2];
sx q[2];
rz(0.25279754) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.0391024) q[1];
sx q[1];
rz(-0.63345816) q[1];
sx q[1];
rz(-1.6879115) q[1];
x q[2];
rz(0.071482166) q[3];
sx q[3];
rz(-1.051828) q[3];
sx q[3];
rz(-1.9148358) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.79036698) q[2];
sx q[2];
rz(-2.6595317) q[2];
sx q[2];
rz(-0.17769979) q[2];
rz(1.7443582) q[3];
sx q[3];
rz(-1.1763562) q[3];
sx q[3];
rz(2.7241404) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.099667065) q[0];
sx q[0];
rz(-0.61138994) q[0];
sx q[0];
rz(1.1055111) q[0];
rz(0.28327495) q[1];
sx q[1];
rz(-2.6073644) q[1];
sx q[1];
rz(-1.6800605) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.71587976) q[0];
sx q[0];
rz(-1.0277896) q[0];
sx q[0];
rz(-2.3433861) q[0];
x q[1];
rz(-0.93223946) q[2];
sx q[2];
rz(-1.9912212) q[2];
sx q[2];
rz(-1.0300385) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-3.0431598) q[1];
sx q[1];
rz(-0.57505703) q[1];
sx q[1];
rz(-2.9270323) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.57638611) q[3];
sx q[3];
rz(-1.295305) q[3];
sx q[3];
rz(-2.798693) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.0561515) q[2];
sx q[2];
rz(-1.5458919) q[2];
sx q[2];
rz(0.20480569) q[2];
rz(2.0917995) q[3];
sx q[3];
rz(-0.27725163) q[3];
sx q[3];
rz(2.7753196) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7751223) q[0];
sx q[0];
rz(-2.6752052) q[0];
sx q[0];
rz(-0.033576641) q[0];
rz(1.6665943) q[1];
sx q[1];
rz(-0.74445236) q[1];
sx q[1];
rz(-1.9974476) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0559674) q[0];
sx q[0];
rz(-2.6542943) q[0];
sx q[0];
rz(3.0874599) q[0];
rz(-pi) q[1];
rz(-1.894925) q[2];
sx q[2];
rz(-1.6648544) q[2];
sx q[2];
rz(-2.6878074) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.5976286) q[1];
sx q[1];
rz(-1.9457327) q[1];
sx q[1];
rz(0.72245325) q[1];
rz(-pi) q[2];
rz(1.9162634) q[3];
sx q[3];
rz(-1.8356712) q[3];
sx q[3];
rz(-0.027907413) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.7555776) q[2];
sx q[2];
rz(-0.80499804) q[2];
sx q[2];
rz(1.8899274) q[2];
rz(1.7898611) q[3];
sx q[3];
rz(-0.91910619) q[3];
sx q[3];
rz(2.1613817) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1397454) q[0];
sx q[0];
rz(-2.503105) q[0];
sx q[0];
rz(-1.2563323) q[0];
rz(0.095257692) q[1];
sx q[1];
rz(-2.2525411) q[1];
sx q[1];
rz(-2.6108066) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9448816) q[0];
sx q[0];
rz(-1.8183578) q[0];
sx q[0];
rz(2.976368) q[0];
x q[1];
rz(-3.0246434) q[2];
sx q[2];
rz(-2.8215234) q[2];
sx q[2];
rz(2.9695599) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.8618882) q[1];
sx q[1];
rz(-1.8100396) q[1];
sx q[1];
rz(1.7025392) q[1];
rz(-pi) q[2];
x q[2];
rz(0.44638951) q[3];
sx q[3];
rz(-1.4671031) q[3];
sx q[3];
rz(-0.54561347) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.329616) q[2];
sx q[2];
rz(-1.5831542) q[2];
sx q[2];
rz(-0.54083332) q[2];
rz(0.81816188) q[3];
sx q[3];
rz(-1.6988138) q[3];
sx q[3];
rz(-2.5087859) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
sx q[3];
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
rz(-2.5166017) q[0];
sx q[0];
rz(-1.3636959) q[0];
sx q[0];
rz(2.7428108) q[0];
rz(-1.7575691) q[1];
sx q[1];
rz(-2.3868491) q[1];
sx q[1];
rz(0.049364518) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2731664) q[0];
sx q[0];
rz(-3.0169562) q[0];
sx q[0];
rz(-1.0602555) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9513543) q[2];
sx q[2];
rz(-2.8597288) q[2];
sx q[2];
rz(-2.6388002) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.4409742) q[1];
sx q[1];
rz(-1.5262394) q[1];
sx q[1];
rz(-0.99127165) q[1];
rz(2.8667198) q[3];
sx q[3];
rz(-2.4913553) q[3];
sx q[3];
rz(-0.45675983) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.12386879) q[2];
sx q[2];
rz(-1.8327291) q[2];
sx q[2];
rz(1.6309942) q[2];
rz(0.92783582) q[3];
sx q[3];
rz(-1.3414914) q[3];
sx q[3];
rz(-1.9063037) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8265726) q[0];
sx q[0];
rz(-2.6980504) q[0];
sx q[0];
rz(2.0499688) q[0];
rz(-2.0768968) q[1];
sx q[1];
rz(-2.6152492) q[1];
sx q[1];
rz(-1.1660887) q[1];
rz(2.0666368) q[2];
sx q[2];
rz(-2.0315764) q[2];
sx q[2];
rz(-2.365664) q[2];
rz(-2.06729) q[3];
sx q[3];
rz(-0.60582325) q[3];
sx q[3];
rz(-0.7614991) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
