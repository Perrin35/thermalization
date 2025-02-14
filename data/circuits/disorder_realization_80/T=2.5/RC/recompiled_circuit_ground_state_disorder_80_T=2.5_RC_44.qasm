OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.03117938) q[0];
sx q[0];
rz(-0.92209446) q[0];
sx q[0];
rz(2.3715012) q[0];
rz(-2.4251179) q[1];
sx q[1];
rz(-0.94243503) q[1];
sx q[1];
rz(-0.64185774) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7866369) q[0];
sx q[0];
rz(-1.4572976) q[0];
sx q[0];
rz(0.20239511) q[0];
rz(2.4840706) q[2];
sx q[2];
rz(-0.0018334963) q[2];
sx q[2];
rz(-0.73500234) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.227046) q[1];
sx q[1];
rz(-1.4150591) q[1];
sx q[1];
rz(-0.36392938) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.558488) q[3];
sx q[3];
rz(-2.2197086) q[3];
sx q[3];
rz(2.976604) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.98051071) q[2];
sx q[2];
rz(-2.0445721) q[2];
sx q[2];
rz(2.3316627) q[2];
rz(0.03446456) q[3];
sx q[3];
rz(-0.66027111) q[3];
sx q[3];
rz(-0.1035498) q[3];
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
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.63475364) q[0];
sx q[0];
rz(-2.6320808) q[0];
sx q[0];
rz(2.4822045) q[0];
rz(1.7022645) q[1];
sx q[1];
rz(-1.5123475) q[1];
sx q[1];
rz(-2.6606681) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.07671) q[0];
sx q[0];
rz(-2.0079029) q[0];
sx q[0];
rz(-2.9124898) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4056817) q[2];
sx q[2];
rz(-0.99756587) q[2];
sx q[2];
rz(-0.57511759) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.9316599) q[1];
sx q[1];
rz(-1.6012962) q[1];
sx q[1];
rz(1.362544) q[1];
rz(2.0388076) q[3];
sx q[3];
rz(-2.0986522) q[3];
sx q[3];
rz(-2.6754987) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.7646358) q[2];
sx q[2];
rz(-2.3048293) q[2];
sx q[2];
rz(0.62977201) q[2];
rz(-1.9624286) q[3];
sx q[3];
rz(-2.4383014) q[3];
sx q[3];
rz(2.0298957) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.69029194) q[0];
sx q[0];
rz(-3.1307104) q[0];
sx q[0];
rz(0.96845281) q[0];
rz(0.14006607) q[1];
sx q[1];
rz(-1.3532956) q[1];
sx q[1];
rz(-2.5586939) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66630689) q[0];
sx q[0];
rz(-0.17284849) q[0];
sx q[0];
rz(-0.26474196) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3034322) q[2];
sx q[2];
rz(-1.7485022) q[2];
sx q[2];
rz(-3.1237683) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.1000468) q[1];
sx q[1];
rz(-0.94521451) q[1];
sx q[1];
rz(0.85432521) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.85286136) q[3];
sx q[3];
rz(-2.7282469) q[3];
sx q[3];
rz(2.1242382) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.53050238) q[2];
sx q[2];
rz(-0.077433057) q[2];
sx q[2];
rz(-2.9275242) q[2];
rz(2.8392082) q[3];
sx q[3];
rz(-0.74367911) q[3];
sx q[3];
rz(0.42588699) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9532303) q[0];
sx q[0];
rz(-2.132405) q[0];
sx q[0];
rz(2.0971712) q[0];
rz(-0.87772477) q[1];
sx q[1];
rz(-1.5889771) q[1];
sx q[1];
rz(-0.16876076) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2917738) q[0];
sx q[0];
rz(-0.55233228) q[0];
sx q[0];
rz(1.1293579) q[0];
rz(1.243882) q[2];
sx q[2];
rz(-1.765328) q[2];
sx q[2];
rz(3.1049797) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.5280594) q[1];
sx q[1];
rz(-0.643104) q[1];
sx q[1];
rz(2.7763441) q[1];
rz(-2.2764858) q[3];
sx q[3];
rz(-1.5388515) q[3];
sx q[3];
rz(0.77199304) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.7793444) q[2];
sx q[2];
rz(-1.8579973) q[2];
sx q[2];
rz(1.0603504) q[2];
rz(-2.849071) q[3];
sx q[3];
rz(-0.7258324) q[3];
sx q[3];
rz(-0.084107548) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(-0.58920687) q[0];
sx q[0];
rz(-0.64738208) q[0];
sx q[0];
rz(0.86828434) q[0];
rz(-0.6262511) q[1];
sx q[1];
rz(-1.7528844) q[1];
sx q[1];
rz(-1.3583604) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.043644301) q[0];
sx q[0];
rz(-1.4820423) q[0];
sx q[0];
rz(0.28907816) q[0];
x q[1];
rz(0.38154885) q[2];
sx q[2];
rz(-1.7722436) q[2];
sx q[2];
rz(2.3214981) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.080278668) q[1];
sx q[1];
rz(-1.3484203) q[1];
sx q[1];
rz(0.23134065) q[1];
rz(-pi) q[2];
rz(0.56760889) q[3];
sx q[3];
rz(-2.4205812) q[3];
sx q[3];
rz(0.22338671) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.2061578) q[2];
sx q[2];
rz(-2.3249966) q[2];
sx q[2];
rz(0.52331501) q[2];
rz(0.37007904) q[3];
sx q[3];
rz(-2.3612634) q[3];
sx q[3];
rz(1.3453329) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10941457) q[0];
sx q[0];
rz(-2.9258756) q[0];
sx q[0];
rz(-0.35414645) q[0];
rz(2.1996563) q[1];
sx q[1];
rz(-1.4553921) q[1];
sx q[1];
rz(1.3166434) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.38625941) q[0];
sx q[0];
rz(-1.0585896) q[0];
sx q[0];
rz(-1.451158) q[0];
rz(-pi) q[1];
rz(-1.5572335) q[2];
sx q[2];
rz(-2.1199653) q[2];
sx q[2];
rz(1.4209335) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.2159913) q[1];
sx q[1];
rz(-2.9570438) q[1];
sx q[1];
rz(0.33949236) q[1];
rz(2.7735071) q[3];
sx q[3];
rz(-2.4637665) q[3];
sx q[3];
rz(-2.8742636) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.8159304) q[2];
sx q[2];
rz(-0.38620913) q[2];
sx q[2];
rz(-0.33829921) q[2];
rz(-2.6650688) q[3];
sx q[3];
rz(-2.3881113) q[3];
sx q[3];
rz(0.34059718) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7020096) q[0];
sx q[0];
rz(-1.5583353) q[0];
sx q[0];
rz(-0.53771341) q[0];
rz(1.733755) q[1];
sx q[1];
rz(-2.6815963) q[1];
sx q[1];
rz(-0.62526155) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.94948927) q[0];
sx q[0];
rz(-1.92983) q[0];
sx q[0];
rz(0.60199317) q[0];
x q[1];
rz(1.6444667) q[2];
sx q[2];
rz(-2.1084614) q[2];
sx q[2];
rz(1.0755838) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-9*pi/10) q[1];
sx q[1];
rz(-0.64439595) q[1];
sx q[1];
rz(2.7792536) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2042562) q[3];
sx q[3];
rz(-2.3126384) q[3];
sx q[3];
rz(0.77223611) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.9739146) q[2];
sx q[2];
rz(-2.6748071) q[2];
sx q[2];
rz(-1.7612877) q[2];
rz(-2.7050833) q[3];
sx q[3];
rz(-1.0218388) q[3];
sx q[3];
rz(2.4280587) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9647144) q[0];
sx q[0];
rz(-0.62129337) q[0];
sx q[0];
rz(-2.6253413) q[0];
rz(0.77975726) q[1];
sx q[1];
rz(-0.98194352) q[1];
sx q[1];
rz(2.0947184) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8314227) q[0];
sx q[0];
rz(-2.8084233) q[0];
sx q[0];
rz(-0.45402558) q[0];
x q[1];
rz(2.0798415) q[2];
sx q[2];
rz(-2.0915481) q[2];
sx q[2];
rz(2.6428509) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.88491625) q[1];
sx q[1];
rz(-3.0312928) q[1];
sx q[1];
rz(2.1259456) q[1];
rz(0.47886301) q[3];
sx q[3];
rz(-1.9267462) q[3];
sx q[3];
rz(-1.009481) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.20389916) q[2];
sx q[2];
rz(-0.1940618) q[2];
sx q[2];
rz(0.95721179) q[2];
rz(2.8474478) q[3];
sx q[3];
rz(-2.011994) q[3];
sx q[3];
rz(-2.8637776) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1421563) q[0];
sx q[0];
rz(-0.23040982) q[0];
sx q[0];
rz(2.9545422) q[0];
rz(-2.710178) q[1];
sx q[1];
rz(-0.43771935) q[1];
sx q[1];
rz(1.6492708) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2628415) q[0];
sx q[0];
rz(-3.0410112) q[0];
sx q[0];
rz(-2.350303) q[0];
rz(-pi) q[1];
rz(2.0772821) q[2];
sx q[2];
rz(-0.37831719) q[2];
sx q[2];
rz(-2.0997206) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.75892657) q[1];
sx q[1];
rz(-1.1622283) q[1];
sx q[1];
rz(2.9620785) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.727333) q[3];
sx q[3];
rz(-1.6633342) q[3];
sx q[3];
rz(-2.375556) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.46818587) q[2];
sx q[2];
rz(-2.2298721) q[2];
sx q[2];
rz(2.1569596) q[2];
rz(2.6268688) q[3];
sx q[3];
rz(-2.616021) q[3];
sx q[3];
rz(-1.908186) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.89727) q[0];
sx q[0];
rz(-1.4800625) q[0];
sx q[0];
rz(-0.75862128) q[0];
rz(-1.9482535) q[1];
sx q[1];
rz(-1.9494282) q[1];
sx q[1];
rz(-1.6903445) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.25989) q[0];
sx q[0];
rz(-0.30891793) q[0];
sx q[0];
rz(-2.3007459) q[0];
x q[1];
rz(-1.7733712) q[2];
sx q[2];
rz(-1.612101) q[2];
sx q[2];
rz(1.6172723) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.7566061) q[1];
sx q[1];
rz(-2.1861417) q[1];
sx q[1];
rz(-0.032757515) q[1];
x q[2];
rz(-1.162446) q[3];
sx q[3];
rz(-1.2733885) q[3];
sx q[3];
rz(1.698146) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.54185581) q[2];
sx q[2];
rz(-0.92076045) q[2];
sx q[2];
rz(-0.80121458) q[2];
rz(-0.93252212) q[3];
sx q[3];
rz(-1.2152117) q[3];
sx q[3];
rz(0.41509375) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5634609) q[0];
sx q[0];
rz(-1.7628071) q[0];
sx q[0];
rz(2.5264869) q[0];
rz(0.32147944) q[1];
sx q[1];
rz(-2.165806) q[1];
sx q[1];
rz(1.4478366) q[1];
rz(0.16531113) q[2];
sx q[2];
rz(-1.4480235) q[2];
sx q[2];
rz(1.5245246) q[2];
rz(0.88820171) q[3];
sx q[3];
rz(-0.35590469) q[3];
sx q[3];
rz(-1.1682846) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
