OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.7251627) q[0];
sx q[0];
rz(0.13983146) q[0];
sx q[0];
rz(10.034372) q[0];
rz(0.66863376) q[1];
sx q[1];
rz(-2.2761087) q[1];
sx q[1];
rz(3.0545711) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4092769) q[0];
sx q[0];
rz(-1.1999994) q[0];
sx q[0];
rz(-0.91863527) q[0];
rz(-pi) q[1];
rz(-0.1331698) q[2];
sx q[2];
rz(-2.3699017) q[2];
sx q[2];
rz(2.0193677) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.062292369) q[1];
sx q[1];
rz(-1.0907409) q[1];
sx q[1];
rz(-2.9855707) q[1];
rz(3.0356785) q[3];
sx q[3];
rz(-2.7273791) q[3];
sx q[3];
rz(0.29990444) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.13880754) q[2];
sx q[2];
rz(-1.3802718) q[2];
sx q[2];
rz(2.7677317) q[2];
rz(-2.8047681) q[3];
sx q[3];
rz(-1.5954433) q[3];
sx q[3];
rz(0.22836223) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8841298) q[0];
sx q[0];
rz(-0.3733491) q[0];
sx q[0];
rz(-1.9467547) q[0];
rz(-3.0589814) q[1];
sx q[1];
rz(-1.1673085) q[1];
sx q[1];
rz(3.1412178) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9276792) q[0];
sx q[0];
rz(-1.1535026) q[0];
sx q[0];
rz(0.37950619) q[0];
x q[1];
rz(-2.9412494) q[2];
sx q[2];
rz(-1.3140972) q[2];
sx q[2];
rz(2.8368907) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.5217168) q[1];
sx q[1];
rz(-2.4881425) q[1];
sx q[1];
rz(0.39342777) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5156636) q[3];
sx q[3];
rz(-2.030636) q[3];
sx q[3];
rz(2.9210747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.3327545) q[2];
sx q[2];
rz(-0.19503441) q[2];
sx q[2];
rz(-0.17671281) q[2];
rz(-2.3475032) q[3];
sx q[3];
rz(-0.77787557) q[3];
sx q[3];
rz(-0.55364048) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90213838) q[0];
sx q[0];
rz(-0.70394009) q[0];
sx q[0];
rz(-0.40476558) q[0];
rz(1.8602712) q[1];
sx q[1];
rz(-0.57360137) q[1];
sx q[1];
rz(1.8331029) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1256516) q[0];
sx q[0];
rz(-1.1106655) q[0];
sx q[0];
rz(2.5815651) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8070418) q[2];
sx q[2];
rz(-2.0958732) q[2];
sx q[2];
rz(-0.46307785) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.96924671) q[1];
sx q[1];
rz(-1.3174651) q[1];
sx q[1];
rz(2.2971056) q[1];
rz(-pi) q[2];
rz(-1.6609459) q[3];
sx q[3];
rz(-1.1336859) q[3];
sx q[3];
rz(-2.3428832) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.219316) q[2];
sx q[2];
rz(-2.6185161) q[2];
sx q[2];
rz(-0.02040872) q[2];
rz(2.0698047) q[3];
sx q[3];
rz(-2.0139549) q[3];
sx q[3];
rz(2.4782457) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14169176) q[0];
sx q[0];
rz(-0.90943709) q[0];
sx q[0];
rz(0.91598696) q[0];
rz(2.6782716) q[1];
sx q[1];
rz(-2.0756192) q[1];
sx q[1];
rz(1.0571009) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4823285) q[0];
sx q[0];
rz(-1.5963012) q[0];
sx q[0];
rz(2.7042424) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5713423) q[2];
sx q[2];
rz(-1.0409365) q[2];
sx q[2];
rz(2.3392764) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-3.108236) q[1];
sx q[1];
rz(-1.557096) q[1];
sx q[1];
rz(1.800888) q[1];
rz(2.3841342) q[3];
sx q[3];
rz(-1.7372903) q[3];
sx q[3];
rz(-1.7050626) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.2477734) q[2];
sx q[2];
rz(-2.5017068) q[2];
sx q[2];
rz(-0.34269732) q[2];
rz(1.4438859) q[3];
sx q[3];
rz(-2.2659437) q[3];
sx q[3];
rz(0.7152344) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7856359) q[0];
sx q[0];
rz(-2.5781093) q[0];
sx q[0];
rz(-3.0551531) q[0];
rz(1.7516288) q[1];
sx q[1];
rz(-0.93170634) q[1];
sx q[1];
rz(-2.6729029) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.3062895) q[0];
sx q[0];
rz(-1.9625184) q[0];
sx q[0];
rz(-0.83175559) q[0];
x q[1];
rz(0.30829633) q[2];
sx q[2];
rz(-1.397965) q[2];
sx q[2];
rz(1.0002713) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.75697452) q[1];
sx q[1];
rz(-2.9006835) q[1];
sx q[1];
rz(0.5730281) q[1];
x q[2];
rz(2.3272446) q[3];
sx q[3];
rz(-1.6296248) q[3];
sx q[3];
rz(-1.3853663) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.9995352) q[2];
sx q[2];
rz(-2.2613328) q[2];
sx q[2];
rz(-2.2772677) q[2];
rz(-2.9243829) q[3];
sx q[3];
rz(-0.61621284) q[3];
sx q[3];
rz(2.8488081) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42751673) q[0];
sx q[0];
rz(-1.4070516) q[0];
sx q[0];
rz(-2.9445904) q[0];
rz(1.3621832) q[1];
sx q[1];
rz(-1.4809337) q[1];
sx q[1];
rz(0.33624712) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.67664458) q[0];
sx q[0];
rz(-1.1332382) q[0];
sx q[0];
rz(2.4462571) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5480812) q[2];
sx q[2];
rz(-2.0715908) q[2];
sx q[2];
rz(1.3154495) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.5832526) q[1];
sx q[1];
rz(-1.6181769) q[1];
sx q[1];
rz(-2.5517795) q[1];
rz(-pi) q[2];
rz(-1.7429966) q[3];
sx q[3];
rz(-0.88168722) q[3];
sx q[3];
rz(-2.8840051) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.53529915) q[2];
sx q[2];
rz(-1.4761816) q[2];
sx q[2];
rz(2.2952648) q[2];
rz(-1.9019295) q[3];
sx q[3];
rz(-1.9097493) q[3];
sx q[3];
rz(3.1404176) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7476615) q[0];
sx q[0];
rz(-2.1037536) q[0];
sx q[0];
rz(2.8826662) q[0];
rz(-1.7954284) q[1];
sx q[1];
rz(-1.7566453) q[1];
sx q[1];
rz(-1.0940201) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0518814) q[0];
sx q[0];
rz(-1.2298755) q[0];
sx q[0];
rz(1.107035) q[0];
rz(-1.6254243) q[2];
sx q[2];
rz(-2.9187181) q[2];
sx q[2];
rz(-2.5958027) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.5309696) q[1];
sx q[1];
rz(-1.5429284) q[1];
sx q[1];
rz(1.665297) q[1];
rz(-pi) q[2];
x q[2];
rz(0.21906994) q[3];
sx q[3];
rz(-3.0411358) q[3];
sx q[3];
rz(-0.62517525) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.1743494) q[2];
sx q[2];
rz(-1.7847585) q[2];
sx q[2];
rz(2.8209177) q[2];
rz(-2.705412) q[3];
sx q[3];
rz(-0.47302055) q[3];
sx q[3];
rz(-0.44803739) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86673474) q[0];
sx q[0];
rz(-2.6651356) q[0];
sx q[0];
rz(2.136769) q[0];
rz(-2.5514305) q[1];
sx q[1];
rz(-2.2361123) q[1];
sx q[1];
rz(0.80387962) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3748976) q[0];
sx q[0];
rz(-0.98537966) q[0];
sx q[0];
rz(-1.6181519) q[0];
rz(-pi) q[1];
rz(2.8477746) q[2];
sx q[2];
rz(-2.567798) q[2];
sx q[2];
rz(2.3125355) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.5775634) q[1];
sx q[1];
rz(-1.659435) q[1];
sx q[1];
rz(-0.010239756) q[1];
x q[2];
rz(1.6611093) q[3];
sx q[3];
rz(-2.2323425) q[3];
sx q[3];
rz(-2.6722398) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.33264318) q[2];
sx q[2];
rz(-2.2968764) q[2];
sx q[2];
rz(-1.0104898) q[2];
rz(-2.5381952) q[3];
sx q[3];
rz(-1.5248652) q[3];
sx q[3];
rz(2.5914014) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6036966) q[0];
sx q[0];
rz(-1.864707) q[0];
sx q[0];
rz(2.7340775) q[0];
rz(-0.28911668) q[1];
sx q[1];
rz(-1.1228077) q[1];
sx q[1];
rz(0.75072748) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5068881) q[0];
sx q[0];
rz(-1.8543715) q[0];
sx q[0];
rz(2.8371235) q[0];
x q[1];
rz(2.8620371) q[2];
sx q[2];
rz(-1.4923555) q[2];
sx q[2];
rz(0.58141764) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.1474485) q[1];
sx q[1];
rz(-1.7732883) q[1];
sx q[1];
rz(-1.8840428) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.30131807) q[3];
sx q[3];
rz(-1.5958061) q[3];
sx q[3];
rz(0.1015639) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.0687381) q[2];
sx q[2];
rz(-0.72454238) q[2];
sx q[2];
rz(-1.9753974) q[2];
rz(-1.3646305) q[3];
sx q[3];
rz(-2.3659673) q[3];
sx q[3];
rz(0.021818074) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5831379) q[0];
sx q[0];
rz(-2.3151509) q[0];
sx q[0];
rz(1.7512084) q[0];
rz(-0.33069262) q[1];
sx q[1];
rz(-0.76534098) q[1];
sx q[1];
rz(1.4601382) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6351785) q[0];
sx q[0];
rz(-1.9902475) q[0];
sx q[0];
rz(2.123453) q[0];
x q[1];
rz(1.5286469) q[2];
sx q[2];
rz(-2.6372006) q[2];
sx q[2];
rz(1.1500037) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.3223443) q[1];
sx q[1];
rz(-1.7548314) q[1];
sx q[1];
rz(-0.40477246) q[1];
rz(-0.71698935) q[3];
sx q[3];
rz(-0.90962142) q[3];
sx q[3];
rz(2.9644074) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.7918487) q[2];
sx q[2];
rz(-1.8169553) q[2];
sx q[2];
rz(-1.4617408) q[2];
rz(-1.1200303) q[3];
sx q[3];
rz(-0.62168613) q[3];
sx q[3];
rz(2.6859443) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.6256975) q[0];
sx q[0];
rz(-1.5213756) q[0];
sx q[0];
rz(-1.9949927) q[0];
rz(-1.3810146) q[1];
sx q[1];
rz(-1.8534503) q[1];
sx q[1];
rz(1.9402515) q[1];
rz(-2.3537221) q[2];
sx q[2];
rz(-1.8265767) q[2];
sx q[2];
rz(1.5192601) q[2];
rz(-1.9588884) q[3];
sx q[3];
rz(-2.2065065) q[3];
sx q[3];
rz(1.828215) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
