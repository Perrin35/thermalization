OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.0948148) q[0];
sx q[0];
rz(4.2098213) q[0];
sx q[0];
rz(9.8888483) q[0];
rz(-1.1820473) q[1];
sx q[1];
rz(-3.074475) q[1];
sx q[1];
rz(1.8571412) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0604295) q[0];
sx q[0];
rz(-1.2788749) q[0];
sx q[0];
rz(-1.1638327) q[0];
rz(-1.60381) q[2];
sx q[2];
rz(-1.0812949) q[2];
sx q[2];
rz(2.2061493) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.9570934) q[1];
sx q[1];
rz(-2.543499) q[1];
sx q[1];
rz(-2.4615272) q[1];
x q[2];
rz(1.5873317) q[3];
sx q[3];
rz(-1.4201122) q[3];
sx q[3];
rz(0.93391234) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.7444732) q[2];
sx q[2];
rz(-2.2019272) q[2];
sx q[2];
rz(2.822067) q[2];
rz(-2.5630991) q[3];
sx q[3];
rz(-2.6632023) q[3];
sx q[3];
rz(-2.4676676) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4085061) q[0];
sx q[0];
rz(-1.5717614) q[0];
sx q[0];
rz(2.615036) q[0];
rz(2.5780442) q[1];
sx q[1];
rz(-1.5310409) q[1];
sx q[1];
rz(2.3449576) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2731199) q[0];
sx q[0];
rz(-0.94808775) q[0];
sx q[0];
rz(2.7362583) q[0];
rz(-0.57467069) q[2];
sx q[2];
rz(-1.7386912) q[2];
sx q[2];
rz(1.4603953) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.7318774) q[1];
sx q[1];
rz(-0.60740031) q[1];
sx q[1];
rz(-0.022547988) q[1];
rz(-2.6437831) q[3];
sx q[3];
rz(-2.1216765) q[3];
sx q[3];
rz(0.8964955) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.18156302) q[2];
sx q[2];
rz(-1.3788362) q[2];
sx q[2];
rz(-0.78655085) q[2];
rz(-2.6484047) q[3];
sx q[3];
rz(-1.9599873) q[3];
sx q[3];
rz(0.44737679) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2717473) q[0];
sx q[0];
rz(-0.87690502) q[0];
sx q[0];
rz(-1.7012117) q[0];
rz(0.72021833) q[1];
sx q[1];
rz(-2.5455988) q[1];
sx q[1];
rz(2.4386141) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7754606) q[0];
sx q[0];
rz(-1.4538987) q[0];
sx q[0];
rz(1.2890105) q[0];
x q[1];
rz(-1.9269283) q[2];
sx q[2];
rz(-2.0424358) q[2];
sx q[2];
rz(-0.65442649) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.30333334) q[1];
sx q[1];
rz(-1.5240182) q[1];
sx q[1];
rz(0.10951885) q[1];
rz(-pi) q[2];
rz(2.968077) q[3];
sx q[3];
rz(-0.95458889) q[3];
sx q[3];
rz(0.71615744) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.8748223) q[2];
sx q[2];
rz(-1.6569933) q[2];
sx q[2];
rz(-2.2300569) q[2];
rz(0.91397816) q[3];
sx q[3];
rz(-1.4303047) q[3];
sx q[3];
rz(-2.3142557) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5754958) q[0];
sx q[0];
rz(-0.63593447) q[0];
sx q[0];
rz(2.3655868) q[0];
rz(1.874118) q[1];
sx q[1];
rz(-1.1088561) q[1];
sx q[1];
rz(-2.5783096) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8877836) q[0];
sx q[0];
rz(-0.94622181) q[0];
sx q[0];
rz(2.8028691) q[0];
rz(-pi) q[1];
rz(-0.40043719) q[2];
sx q[2];
rz(-2.4106328) q[2];
sx q[2];
rz(-1.2618582) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.54484425) q[1];
sx q[1];
rz(-1.2990284) q[1];
sx q[1];
rz(-1.5775058) q[1];
rz(2.8068845) q[3];
sx q[3];
rz(-2.3957806) q[3];
sx q[3];
rz(-2.0832182) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.6528066) q[2];
sx q[2];
rz(-2.2061009) q[2];
sx q[2];
rz(-0.164786) q[2];
rz(-0.22848836) q[3];
sx q[3];
rz(-2.861664) q[3];
sx q[3];
rz(-0.94648186) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.27424681) q[0];
sx q[0];
rz(-2.2773401) q[0];
sx q[0];
rz(-0.85246032) q[0];
rz(-2.7903941) q[1];
sx q[1];
rz(-0.57792592) q[1];
sx q[1];
rz(-1.9794827) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4974385) q[0];
sx q[0];
rz(-0.69735202) q[0];
sx q[0];
rz(0.63875142) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0931307) q[2];
sx q[2];
rz(-1.9198717) q[2];
sx q[2];
rz(-0.77490865) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.8309161) q[1];
sx q[1];
rz(-1.1569996) q[1];
sx q[1];
rz(2.3526741) q[1];
rz(1.5557628) q[3];
sx q[3];
rz(-2.2923277) q[3];
sx q[3];
rz(2.815747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.6953485) q[2];
sx q[2];
rz(-2.7480795) q[2];
sx q[2];
rz(1.8943141) q[2];
rz(-3.128483) q[3];
sx q[3];
rz(-1.7316827) q[3];
sx q[3];
rz(2.9124027) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86876774) q[0];
sx q[0];
rz(-1.8936963) q[0];
sx q[0];
rz(0.75063467) q[0];
rz(-2.0320832) q[1];
sx q[1];
rz(-0.37056071) q[1];
sx q[1];
rz(1.3060588) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7040946) q[0];
sx q[0];
rz(-1.1878345) q[0];
sx q[0];
rz(1.0136481) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6779804) q[2];
sx q[2];
rz(-1.7525275) q[2];
sx q[2];
rz(-0.52464991) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.2893147) q[1];
sx q[1];
rz(-1.5783974) q[1];
sx q[1];
rz(-0.72035933) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7197051) q[3];
sx q[3];
rz(-0.9871261) q[3];
sx q[3];
rz(2.7910809) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.3874454) q[2];
sx q[2];
rz(-2.1345317) q[2];
sx q[2];
rz(0.60069096) q[2];
rz(-2.1389652) q[3];
sx q[3];
rz(-2.6032344) q[3];
sx q[3];
rz(0.35282648) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7464741) q[0];
sx q[0];
rz(-1.4557319) q[0];
sx q[0];
rz(1.0466928) q[0];
rz(1.5294317) q[1];
sx q[1];
rz(-1.7144831) q[1];
sx q[1];
rz(-2.7244862) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1244303) q[0];
sx q[0];
rz(-0.22148795) q[0];
sx q[0];
rz(-1.4593967) q[0];
rz(2.3796758) q[2];
sx q[2];
rz(-0.41315213) q[2];
sx q[2];
rz(-0.44943902) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.6542146) q[1];
sx q[1];
rz(-1.1114792) q[1];
sx q[1];
rz(-2.3801801) q[1];
x q[2];
rz(-2.7583371) q[3];
sx q[3];
rz(-0.68985046) q[3];
sx q[3];
rz(0.81724973) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.1356915) q[2];
sx q[2];
rz(-2.6689172) q[2];
sx q[2];
rz(-1.3674412) q[2];
rz(-2.588429) q[3];
sx q[3];
rz(-1.4033214) q[3];
sx q[3];
rz(-0.52545351) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.816514) q[0];
sx q[0];
rz(-1.3229609) q[0];
sx q[0];
rz(1.1897855) q[0];
rz(-1.4272383) q[1];
sx q[1];
rz(-0.27806565) q[1];
sx q[1];
rz(0.11238012) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.534879) q[0];
sx q[0];
rz(-1.0977854) q[0];
sx q[0];
rz(-0.19434778) q[0];
rz(-pi) q[1];
rz(-2.8939395) q[2];
sx q[2];
rz(-2.7981749) q[2];
sx q[2];
rz(2.0695956) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.97265128) q[1];
sx q[1];
rz(-2.2691138) q[1];
sx q[1];
rz(0.45291839) q[1];
rz(-1.3120679) q[3];
sx q[3];
rz(-1.6885763) q[3];
sx q[3];
rz(-3.0974914) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.0743951) q[2];
sx q[2];
rz(-1.9125568) q[2];
sx q[2];
rz(1.7929662) q[2];
rz(-1.2049234) q[3];
sx q[3];
rz(-1.7347074) q[3];
sx q[3];
rz(-0.19781923) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.39712054) q[0];
sx q[0];
rz(-1.4646126) q[0];
sx q[0];
rz(-3.1066185) q[0];
rz(-2.2947625) q[1];
sx q[1];
rz(-1.433082) q[1];
sx q[1];
rz(-0.91167489) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.23110403) q[0];
sx q[0];
rz(-1.3470874) q[0];
sx q[0];
rz(-1.847812) q[0];
x q[1];
rz(0.54347221) q[2];
sx q[2];
rz(-1.6795571) q[2];
sx q[2];
rz(-0.47765884) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.1940143) q[1];
sx q[1];
rz(-1.36735) q[1];
sx q[1];
rz(0.46781637) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6058042) q[3];
sx q[3];
rz(-1.1402604) q[3];
sx q[3];
rz(-2.146194) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.70790616) q[2];
sx q[2];
rz(-2.6277379) q[2];
sx q[2];
rz(2.8923477) q[2];
rz(-2.3748659) q[3];
sx q[3];
rz(-1.8079575) q[3];
sx q[3];
rz(-2.7419817) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3578167) q[0];
sx q[0];
rz(-2.0932842) q[0];
sx q[0];
rz(2.7695079) q[0];
rz(0.58139873) q[1];
sx q[1];
rz(-0.77722469) q[1];
sx q[1];
rz(-1.7262329) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2863662) q[0];
sx q[0];
rz(-0.33272538) q[0];
sx q[0];
rz(-0.65441982) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5214228) q[2];
sx q[2];
rz(-1.1405186) q[2];
sx q[2];
rz(-0.60713965) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.2409754) q[1];
sx q[1];
rz(-1.6458578) q[1];
sx q[1];
rz(-1.5013298) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1146718) q[3];
sx q[3];
rz(-1.2380935) q[3];
sx q[3];
rz(-2.016891) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.32594484) q[2];
sx q[2];
rz(-3.0013913) q[2];
sx q[2];
rz(0.20467219) q[2];
rz(1.4137911) q[3];
sx q[3];
rz(-1.1433583) q[3];
sx q[3];
rz(1.0958825) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6407912) q[0];
sx q[0];
rz(-1.4688107) q[0];
sx q[0];
rz(-2.2039913) q[0];
rz(1.5564556) q[1];
sx q[1];
rz(-1.382348) q[1];
sx q[1];
rz(1.4437645) q[1];
rz(-3.0204308) q[2];
sx q[2];
rz(-2.0222752) q[2];
sx q[2];
rz(-3.0565699) q[2];
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
