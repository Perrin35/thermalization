OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.0410864) q[0];
sx q[0];
rz(2.9334928) q[0];
sx q[0];
rz(9.8280332) q[0];
rz(2.9085605) q[1];
sx q[1];
rz(-1.7014528) q[1];
sx q[1];
rz(-0.22388248) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1772645) q[0];
sx q[0];
rz(-1.7814264) q[0];
sx q[0];
rz(2.5151279) q[0];
x q[1];
rz(2.2751132) q[2];
sx q[2];
rz(-1.472855) q[2];
sx q[2];
rz(-0.78664727) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.76624291) q[1];
sx q[1];
rz(-2.2170984) q[1];
sx q[1];
rz(0.28278858) q[1];
x q[2];
rz(-2.3305064) q[3];
sx q[3];
rz(-2.1032277) q[3];
sx q[3];
rz(0.48871751) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.69748059) q[2];
sx q[2];
rz(-0.31284249) q[2];
sx q[2];
rz(-3.1057788) q[2];
rz(2.4114285) q[3];
sx q[3];
rz(-1.1532447) q[3];
sx q[3];
rz(0.20619503) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8908454) q[0];
sx q[0];
rz(-2.146281) q[0];
sx q[0];
rz(-0.19749755) q[0];
rz(0.47710553) q[1];
sx q[1];
rz(-0.66480607) q[1];
sx q[1];
rz(-2.3359931) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.8957414) q[0];
sx q[0];
rz(-0.68994683) q[0];
sx q[0];
rz(2.8474381) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.637403) q[2];
sx q[2];
rz(-1.9506467) q[2];
sx q[2];
rz(2.0646281) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.9979084) q[1];
sx q[1];
rz(-0.98947462) q[1];
sx q[1];
rz(-2.3475926) q[1];
x q[2];
rz(0.74975852) q[3];
sx q[3];
rz(-0.20515144) q[3];
sx q[3];
rz(1.4914025) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.3104559) q[2];
sx q[2];
rz(-1.4996935) q[2];
sx q[2];
rz(-1.7384701) q[2];
rz(2.5095615) q[3];
sx q[3];
rz(-2.8681614) q[3];
sx q[3];
rz(3.0145751) q[3];
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
rz(-pi) q[0];
sx q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3056575) q[0];
sx q[0];
rz(-2.5129565) q[0];
sx q[0];
rz(-2.054731) q[0];
rz(2.944259) q[1];
sx q[1];
rz(-1.5060164) q[1];
sx q[1];
rz(-0.33263439) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.63276382) q[0];
sx q[0];
rz(-2.6869505) q[0];
sx q[0];
rz(-1.2484545) q[0];
x q[1];
rz(1.5015177) q[2];
sx q[2];
rz(-2.4054619) q[2];
sx q[2];
rz(-0.87321216) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.3217993) q[1];
sx q[1];
rz(-1.7463356) q[1];
sx q[1];
rz(-1.790253) q[1];
rz(-pi) q[2];
rz(-0.72686355) q[3];
sx q[3];
rz(-1.7439902) q[3];
sx q[3];
rz(1.3722591) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.7670333) q[2];
sx q[2];
rz(-1.550753) q[2];
sx q[2];
rz(1.4683051) q[2];
rz(-0.0091920216) q[3];
sx q[3];
rz(-2.6720948) q[3];
sx q[3];
rz(-0.79206842) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.511635) q[0];
sx q[0];
rz(-2.503643) q[0];
sx q[0];
rz(-1.6988423) q[0];
rz(-1.0244145) q[1];
sx q[1];
rz(-1.9099648) q[1];
sx q[1];
rz(2.3146497) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31931811) q[0];
sx q[0];
rz(-1.7546904) q[0];
sx q[0];
rz(0.075550373) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2587772) q[2];
sx q[2];
rz(-1.6903631) q[2];
sx q[2];
rz(-1.9012698) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.1093028) q[1];
sx q[1];
rz(-0.48498165) q[1];
sx q[1];
rz(0.46320199) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.8884646) q[3];
sx q[3];
rz(-1.6649705) q[3];
sx q[3];
rz(-1.3741796) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.8168762) q[2];
sx q[2];
rz(-2.3574895) q[2];
sx q[2];
rz(-1.2775705) q[2];
rz(2.4534524) q[3];
sx q[3];
rz(-2.1161931) q[3];
sx q[3];
rz(-2.1791606) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0012896) q[0];
sx q[0];
rz(-1.7004509) q[0];
sx q[0];
rz(-2.9241614) q[0];
rz(-2.8561719) q[1];
sx q[1];
rz(-0.60573429) q[1];
sx q[1];
rz(-2.6434456) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.065154) q[0];
sx q[0];
rz(-0.74238837) q[0];
sx q[0];
rz(-2.0641293) q[0];
x q[1];
rz(0.26136036) q[2];
sx q[2];
rz(-2.0730505) q[2];
sx q[2];
rz(0.72797352) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.0951257) q[1];
sx q[1];
rz(-0.53494638) q[1];
sx q[1];
rz(1.9880268) q[1];
rz(-pi) q[2];
rz(1.3181456) q[3];
sx q[3];
rz(-1.9985191) q[3];
sx q[3];
rz(-1.3633756) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.6238326) q[2];
sx q[2];
rz(-1.0972923) q[2];
sx q[2];
rz(-2.9916054) q[2];
rz(3.127626) q[3];
sx q[3];
rz(-1.59168) q[3];
sx q[3];
rz(-2.6278031) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49195313) q[0];
sx q[0];
rz(-2.5998901) q[0];
sx q[0];
rz(-0.81277043) q[0];
rz(-1.4179519) q[1];
sx q[1];
rz(-1.6287454) q[1];
sx q[1];
rz(1.0636122) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15062885) q[0];
sx q[0];
rz(-1.0526987) q[0];
sx q[0];
rz(1.9512779) q[0];
rz(-pi) q[1];
rz(3.0399357) q[2];
sx q[2];
rz(-0.16451193) q[2];
sx q[2];
rz(-3.0992584) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.2460766) q[1];
sx q[1];
rz(-2.5427596) q[1];
sx q[1];
rz(-0.020222874) q[1];
x q[2];
rz(-0.1891047) q[3];
sx q[3];
rz(-1.8207887) q[3];
sx q[3];
rz(2.7562906) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.873988) q[2];
sx q[2];
rz(-0.559811) q[2];
sx q[2];
rz(-1.416729) q[2];
rz(-0.060976107) q[3];
sx q[3];
rz(-0.91106001) q[3];
sx q[3];
rz(1.4643668) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7199719) q[0];
sx q[0];
rz(-0.82141972) q[0];
sx q[0];
rz(-3.022497) q[0];
rz(-2.6630317) q[1];
sx q[1];
rz(-2.3275972) q[1];
sx q[1];
rz(-2.4518769) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0499303) q[0];
sx q[0];
rz(-2.1436467) q[0];
sx q[0];
rz(0.86688231) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1999667) q[2];
sx q[2];
rz(-2.0304108) q[2];
sx q[2];
rz(0.72061611) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.8552095) q[1];
sx q[1];
rz(-0.70806009) q[1];
sx q[1];
rz(-1.5016012) q[1];
x q[2];
rz(2.2834217) q[3];
sx q[3];
rz(-1.6827876) q[3];
sx q[3];
rz(1.6046765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.0488284) q[2];
sx q[2];
rz(-0.060828716) q[2];
sx q[2];
rz(-2.0737341) q[2];
rz(2.974406) q[3];
sx q[3];
rz(-1.3719631) q[3];
sx q[3];
rz(-3.0021477) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6086513) q[0];
sx q[0];
rz(-0.81357384) q[0];
sx q[0];
rz(-0.90840489) q[0];
rz(0.99336973) q[1];
sx q[1];
rz(-0.16255957) q[1];
sx q[1];
rz(-0.87432528) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4607166) q[0];
sx q[0];
rz(-1.5534288) q[0];
sx q[0];
rz(-0.13219035) q[0];
x q[1];
rz(0.019432391) q[2];
sx q[2];
rz(-1.8284869) q[2];
sx q[2];
rz(-1.2033705) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.5914375) q[1];
sx q[1];
rz(-0.14796013) q[1];
sx q[1];
rz(1.1569389) q[1];
rz(-pi) q[2];
rz(-2.1217974) q[3];
sx q[3];
rz(-1.8809109) q[3];
sx q[3];
rz(-1.1834902) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.61451644) q[2];
sx q[2];
rz(-2.9824342) q[2];
sx q[2];
rz(-0.85421526) q[2];
rz(0.45520511) q[3];
sx q[3];
rz(-2.5762317) q[3];
sx q[3];
rz(2.8860886) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47827569) q[0];
sx q[0];
rz(-1.9238967) q[0];
sx q[0];
rz(-1.664337) q[0];
rz(-1.5746337) q[1];
sx q[1];
rz(-0.68395558) q[1];
sx q[1];
rz(-2.4050567) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2701243) q[0];
sx q[0];
rz(-2.0379319) q[0];
sx q[0];
rz(0.26345912) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4669424) q[2];
sx q[2];
rz(-1.9496634) q[2];
sx q[2];
rz(1.1844289) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.354035) q[1];
sx q[1];
rz(-1.327146) q[1];
sx q[1];
rz(0.32275782) q[1];
x q[2];
rz(-0.97379721) q[3];
sx q[3];
rz(-1.6000119) q[3];
sx q[3];
rz(2.4643498) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.7443202) q[2];
sx q[2];
rz(-0.87817764) q[2];
sx q[2];
rz(1.5273904) q[2];
rz(2.7496036) q[3];
sx q[3];
rz(-1.1992998) q[3];
sx q[3];
rz(-0.049662445) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.117332) q[0];
sx q[0];
rz(-1.4735104) q[0];
sx q[0];
rz(2.9397553) q[0];
rz(-2.9182538) q[1];
sx q[1];
rz(-2.6170862) q[1];
sx q[1];
rz(2.7490659) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.72001624) q[0];
sx q[0];
rz(-0.47552738) q[0];
sx q[0];
rz(-1.3386734) q[0];
x q[1];
rz(-0.052458737) q[2];
sx q[2];
rz(-0.83485583) q[2];
sx q[2];
rz(0.8142161) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.1087677) q[1];
sx q[1];
rz(-0.6995753) q[1];
sx q[1];
rz(2.2684069) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6018312) q[3];
sx q[3];
rz(-2.5934412) q[3];
sx q[3];
rz(0.87566151) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.4203804) q[2];
sx q[2];
rz(-0.98782295) q[2];
sx q[2];
rz(-2.5610899) q[2];
rz(-0.37603363) q[3];
sx q[3];
rz(-0.97085634) q[3];
sx q[3];
rz(1.5413126) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9911983) q[0];
sx q[0];
rz(-1.5491485) q[0];
sx q[0];
rz(-1.3843672) q[0];
rz(-1.4906384) q[1];
sx q[1];
rz(-1.8883659) q[1];
sx q[1];
rz(-2.2813003) q[1];
rz(-2.5096624) q[2];
sx q[2];
rz(-2.2298553) q[2];
sx q[2];
rz(-2.3134445) q[2];
rz(3.0720465) q[3];
sx q[3];
rz(-2.8290717) q[3];
sx q[3];
rz(1.1858218) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
