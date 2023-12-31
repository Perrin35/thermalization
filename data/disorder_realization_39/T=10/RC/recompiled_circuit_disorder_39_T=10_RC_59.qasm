OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.0467779) q[0];
sx q[0];
rz(-1.0682286) q[0];
sx q[0];
rz(2.6775223) q[0];
rz(-1.1820473) q[1];
sx q[1];
rz(-3.074475) q[1];
sx q[1];
rz(-1.2844515) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5285437) q[0];
sx q[0];
rz(-1.9595946) q[0];
sx q[0];
rz(-0.31624985) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5377827) q[2];
sx q[2];
rz(-1.0812949) q[2];
sx q[2];
rz(-2.2061493) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.9570934) q[1];
sx q[1];
rz(-2.543499) q[1];
sx q[1];
rz(-2.4615272) q[1];
x q[2];
rz(-1.5873317) q[3];
sx q[3];
rz(-1.7214805) q[3];
sx q[3];
rz(-2.2076803) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.39711943) q[2];
sx q[2];
rz(-0.93966547) q[2];
sx q[2];
rz(2.822067) q[2];
rz(0.5784936) q[3];
sx q[3];
rz(-2.6632023) q[3];
sx q[3];
rz(-2.4676676) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4085061) q[0];
sx q[0];
rz(-1.5717614) q[0];
sx q[0];
rz(-2.615036) q[0];
rz(0.56354848) q[1];
sx q[1];
rz(-1.5310409) q[1];
sx q[1];
rz(0.79663509) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5986886) q[0];
sx q[0];
rz(-1.2447378) q[0];
sx q[0];
rz(-0.90755264) q[0];
rz(-pi) q[1];
rz(-0.30226207) q[2];
sx q[2];
rz(-2.5455591) q[2];
sx q[2];
rz(-2.9994534) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.7044201) q[1];
sx q[1];
rz(-0.96357268) q[1];
sx q[1];
rz(1.5551268) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.49780952) q[3];
sx q[3];
rz(-1.0199162) q[3];
sx q[3];
rz(0.8964955) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.9600296) q[2];
sx q[2];
rz(-1.7627565) q[2];
sx q[2];
rz(2.3550418) q[2];
rz(-2.6484047) q[3];
sx q[3];
rz(-1.9599873) q[3];
sx q[3];
rz(0.44737679) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2717473) q[0];
sx q[0];
rz(-2.2646876) q[0];
sx q[0];
rz(-1.440381) q[0];
rz(-2.4213743) q[1];
sx q[1];
rz(-0.59599382) q[1];
sx q[1];
rz(-2.4386141) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9706791) q[0];
sx q[0];
rz(-1.2909856) q[0];
sx q[0];
rz(0.12165102) q[0];
rz(-pi) q[1];
rz(1.2146644) q[2];
sx q[2];
rz(-2.0424358) q[2];
sx q[2];
rz(2.4871662) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.30333334) q[1];
sx q[1];
rz(-1.5240182) q[1];
sx q[1];
rz(0.10951885) q[1];
rz(-2.968077) q[3];
sx q[3];
rz(-2.1870038) q[3];
sx q[3];
rz(0.71615744) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.8748223) q[2];
sx q[2];
rz(-1.6569933) q[2];
sx q[2];
rz(-2.2300569) q[2];
rz(-2.2276145) q[3];
sx q[3];
rz(-1.4303047) q[3];
sx q[3];
rz(-2.3142557) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5660969) q[0];
sx q[0];
rz(-2.5056582) q[0];
sx q[0];
rz(2.3655868) q[0];
rz(-1.874118) q[1];
sx q[1];
rz(-1.1088561) q[1];
sx q[1];
rz(-0.56328303) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8877836) q[0];
sx q[0];
rz(-0.94622181) q[0];
sx q[0];
rz(-2.8028691) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2345384) q[2];
sx q[2];
rz(-2.2328937) q[2];
sx q[2];
rz(1.3627571) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.024151) q[1];
sx q[1];
rz(-1.5643331) q[1];
sx q[1];
rz(0.27177377) q[1];
x q[2];
rz(-2.8068845) q[3];
sx q[3];
rz(-0.74581205) q[3];
sx q[3];
rz(1.0583744) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.48878601) q[2];
sx q[2];
rz(-0.93549171) q[2];
sx q[2];
rz(-0.164786) q[2];
rz(2.9131043) q[3];
sx q[3];
rz(-2.861664) q[3];
sx q[3];
rz(2.1951108) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.27424681) q[0];
sx q[0];
rz(-2.2773401) q[0];
sx q[0];
rz(-0.85246032) q[0];
rz(0.35119855) q[1];
sx q[1];
rz(-2.5636667) q[1];
sx q[1];
rz(-1.16211) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1253818) q[0];
sx q[0];
rz(-2.1124766) q[0];
sx q[0];
rz(-2.0340232) q[0];
rz(-0.39761333) q[2];
sx q[2];
rz(-1.0828472) q[2];
sx q[2];
rz(-2.5401126) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.3106766) q[1];
sx q[1];
rz(-1.984593) q[1];
sx q[1];
rz(-2.3526741) q[1];
rz(3.1245072) q[3];
sx q[3];
rz(-0.72165976) q[3];
sx q[3];
rz(2.7929896) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.44624415) q[2];
sx q[2];
rz(-2.7480795) q[2];
sx q[2];
rz(-1.2472786) q[2];
rz(0.013109664) q[3];
sx q[3];
rz(-1.7316827) q[3];
sx q[3];
rz(2.9124027) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86876774) q[0];
sx q[0];
rz(-1.2478963) q[0];
sx q[0];
rz(2.390958) q[0];
rz(-1.1095095) q[1];
sx q[1];
rz(-2.7710319) q[1];
sx q[1];
rz(1.3060588) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36201492) q[0];
sx q[0];
rz(-2.0834196) q[0];
sx q[0];
rz(-2.6984452) q[0];
x q[1];
rz(2.6779804) q[2];
sx q[2];
rz(-1.7525275) q[2];
sx q[2];
rz(-0.52464991) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.8522779) q[1];
sx q[1];
rz(-1.5783974) q[1];
sx q[1];
rz(2.4212333) q[1];
x q[2];
rz(-1.7197051) q[3];
sx q[3];
rz(-0.9871261) q[3];
sx q[3];
rz(0.35051171) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.3874454) q[2];
sx q[2];
rz(-1.0070609) q[2];
sx q[2];
rz(0.60069096) q[2];
rz(2.1389652) q[3];
sx q[3];
rz(-0.53835821) q[3];
sx q[3];
rz(-2.7887662) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7464741) q[0];
sx q[0];
rz(-1.6858608) q[0];
sx q[0];
rz(-2.0948998) q[0];
rz(1.5294317) q[1];
sx q[1];
rz(-1.7144831) q[1];
sx q[1];
rz(0.41710645) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55506599) q[0];
sx q[0];
rz(-1.546372) q[0];
sx q[0];
rz(-1.3506372) q[0];
rz(-pi) q[1];
rz(2.8344526) q[2];
sx q[2];
rz(-1.8516314) q[2];
sx q[2];
rz(1.3020696) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.48737803) q[1];
sx q[1];
rz(-1.1114792) q[1];
sx q[1];
rz(-0.76141255) q[1];
rz(-2.4884175) q[3];
sx q[3];
rz(-1.8110868) q[3];
sx q[3];
rz(-2.0865292) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.0059011857) q[2];
sx q[2];
rz(-0.4726755) q[2];
sx q[2];
rz(1.7741514) q[2];
rz(-0.55316365) q[3];
sx q[3];
rz(-1.4033214) q[3];
sx q[3];
rz(0.52545351) q[3];
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
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32507867) q[0];
sx q[0];
rz(-1.3229609) q[0];
sx q[0];
rz(-1.1897855) q[0];
rz(1.4272383) q[1];
sx q[1];
rz(-2.863527) q[1];
sx q[1];
rz(0.11238012) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.053514078) q[0];
sx q[0];
rz(-1.7435762) q[0];
sx q[0];
rz(2.0515576) q[0];
rz(1.6582279) q[2];
sx q[2];
rz(-1.9033252) q[2];
sx q[2];
rz(-2.33193) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.29490678) q[1];
sx q[1];
rz(-1.2290188) q[1];
sx q[1];
rz(0.81975598) q[1];
x q[2];
rz(1.1376082) q[3];
sx q[3];
rz(-0.28372753) q[3];
sx q[3];
rz(-1.9445436) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.0671976) q[2];
sx q[2];
rz(-1.9125568) q[2];
sx q[2];
rz(-1.3486264) q[2];
rz(1.9366692) q[3];
sx q[3];
rz(-1.7347074) q[3];
sx q[3];
rz(2.9437734) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.39712054) q[0];
sx q[0];
rz(-1.67698) q[0];
sx q[0];
rz(3.1066185) q[0];
rz(-2.2947625) q[1];
sx q[1];
rz(-1.433082) q[1];
sx q[1];
rz(2.2299178) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.002279) q[0];
sx q[0];
rz(-2.7873435) q[0];
sx q[0];
rz(0.87689633) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4439092) q[2];
sx q[2];
rz(-1.03089) q[2];
sx q[2];
rz(-1.9829696) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.1940143) q[1];
sx q[1];
rz(-1.36735) q[1];
sx q[1];
rz(-2.6737763) q[1];
rz(-0.53578844) q[3];
sx q[3];
rz(-1.1402604) q[3];
sx q[3];
rz(2.146194) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.4336865) q[2];
sx q[2];
rz(-2.6277379) q[2];
sx q[2];
rz(-0.24924499) q[2];
rz(2.3748659) q[3];
sx q[3];
rz(-1.3336351) q[3];
sx q[3];
rz(-2.7419817) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3578167) q[0];
sx q[0];
rz(-1.0483085) q[0];
sx q[0];
rz(0.37208474) q[0];
rz(-2.5601939) q[1];
sx q[1];
rz(-2.364368) q[1];
sx q[1];
rz(-1.4153597) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6045195) q[0];
sx q[0];
rz(-1.3086638) q[0];
sx q[0];
rz(1.7781236) q[0];
rz(-2.4731589) q[2];
sx q[2];
rz(-0.738315) q[2];
sx q[2];
rz(2.7065606) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.33503867) q[1];
sx q[1];
rz(-1.6400669) q[1];
sx q[1];
rz(-3.0663504) q[1];
rz(-pi) q[2];
rz(-0.38378999) q[3];
sx q[3];
rz(-1.0597611) q[3];
sx q[3];
rz(-2.8904861) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.8156478) q[2];
sx q[2];
rz(-3.0013913) q[2];
sx q[2];
rz(2.9369205) q[2];
rz(-1.4137911) q[3];
sx q[3];
rz(-1.9982343) q[3];
sx q[3];
rz(1.0958825) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(0.50080147) q[0];
sx q[0];
rz(-1.4688107) q[0];
sx q[0];
rz(-2.2039913) q[0];
rz(1.5564556) q[1];
sx q[1];
rz(-1.382348) q[1];
sx q[1];
rz(1.4437645) q[1];
rz(3.0204308) q[2];
sx q[2];
rz(-1.1193174) q[2];
sx q[2];
rz(0.08502273) q[2];
rz(-0.08269357) q[3];
sx q[3];
rz(-1.0001928) q[3];
sx q[3];
rz(2.115888) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
