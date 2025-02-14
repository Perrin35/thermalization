OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.3815936) q[0];
sx q[0];
rz(-2.969709) q[0];
sx q[0];
rz(-2.5266393) q[0];
rz(-1.7847269) q[1];
sx q[1];
rz(-1.6109799) q[1];
sx q[1];
rz(-0.15712486) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2113094) q[0];
sx q[0];
rz(-2.0727949) q[0];
sx q[0];
rz(0.94663488) q[0];
rz(2.914937) q[2];
sx q[2];
rz(-1.6427411) q[2];
sx q[2];
rz(0.26357108) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.3362017) q[1];
sx q[1];
rz(-1.5861643) q[1];
sx q[1];
rz(-2.2284641) q[1];
rz(-pi) q[2];
rz(0.71977736) q[3];
sx q[3];
rz(-1.3941289) q[3];
sx q[3];
rz(2.9797785) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.5379546) q[2];
sx q[2];
rz(-0.31542888) q[2];
sx q[2];
rz(1.8559378) q[2];
rz(-1.3373318) q[3];
sx q[3];
rz(-2.1888816) q[3];
sx q[3];
rz(-3.0281236) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86785698) q[0];
sx q[0];
rz(-1.8524167) q[0];
sx q[0];
rz(2.5751172) q[0];
rz(1.6784338) q[1];
sx q[1];
rz(-0.17371829) q[1];
sx q[1];
rz(2.4991551) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.89680371) q[0];
sx q[0];
rz(-1.6224846) q[0];
sx q[0];
rz(-1.1079074) q[0];
rz(-pi) q[1];
rz(1.9011097) q[2];
sx q[2];
rz(-2.1871302) q[2];
sx q[2];
rz(-1.7834477) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.21684619) q[1];
sx q[1];
rz(-1.6213776) q[1];
sx q[1];
rz(0.84212007) q[1];
rz(-pi) q[2];
x q[2];
rz(0.94704721) q[3];
sx q[3];
rz(-0.94175856) q[3];
sx q[3];
rz(1.3082023) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.2598339) q[2];
sx q[2];
rz(-0.79462785) q[2];
sx q[2];
rz(-0.83414042) q[2];
rz(-3.0547018) q[3];
sx q[3];
rz(-1.0850685) q[3];
sx q[3];
rz(-2.0104008) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.46383816) q[0];
sx q[0];
rz(-2.0448271) q[0];
sx q[0];
rz(1.0595529) q[0];
rz(2.4067267) q[1];
sx q[1];
rz(-0.015345416) q[1];
sx q[1];
rz(2.9954092) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.138036) q[0];
sx q[0];
rz(-2.6377886) q[0];
sx q[0];
rz(-1.2362859) q[0];
rz(-pi) q[1];
rz(3.1377313) q[2];
sx q[2];
rz(-1.4465824) q[2];
sx q[2];
rz(0.079380097) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.4242658) q[1];
sx q[1];
rz(-1.4098412) q[1];
sx q[1];
rz(-0.15086727) q[1];
rz(-1.8511778) q[3];
sx q[3];
rz(-1.7903847) q[3];
sx q[3];
rz(1.4189347) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.77201468) q[2];
sx q[2];
rz(-2.7507608) q[2];
sx q[2];
rz(0.84103) q[2];
rz(-3.0817025) q[3];
sx q[3];
rz(-1.6570897) q[3];
sx q[3];
rz(-2.4976775) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.197072) q[0];
sx q[0];
rz(-2.6332025) q[0];
sx q[0];
rz(-3.0149241) q[0];
rz(-1.4177167) q[1];
sx q[1];
rz(-0.020514943) q[1];
sx q[1];
rz(0.98974481) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.58137547) q[0];
sx q[0];
rz(-0.84218431) q[0];
sx q[0];
rz(-1.5422264) q[0];
rz(-pi) q[1];
rz(0.085096882) q[2];
sx q[2];
rz(-1.1750882) q[2];
sx q[2];
rz(-2.1425785) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.77074487) q[1];
sx q[1];
rz(-0.49990067) q[1];
sx q[1];
rz(1.1558394) q[1];
rz(-pi) q[2];
rz(2.9007841) q[3];
sx q[3];
rz(-1.7205372) q[3];
sx q[3];
rz(3.0583241) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.4385472) q[2];
sx q[2];
rz(-2.788471) q[2];
sx q[2];
rz(2.9959196) q[2];
rz(-0.4438256) q[3];
sx q[3];
rz(-1.5747109) q[3];
sx q[3];
rz(-2.023196) q[3];
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
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.038427453) q[0];
sx q[0];
rz(-1.225809) q[0];
sx q[0];
rz(-2.3577754) q[0];
rz(1.8812284) q[1];
sx q[1];
rz(-0.0030219373) q[1];
sx q[1];
rz(-0.22617117) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5617666) q[0];
sx q[0];
rz(-1.7975627) q[0];
sx q[0];
rz(-1.0933601) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1532835) q[2];
sx q[2];
rz(-1.8476356) q[2];
sx q[2];
rz(-1.4275488) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.3231103) q[1];
sx q[1];
rz(-0.47814506) q[1];
sx q[1];
rz(1.7341431) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3557303) q[3];
sx q[3];
rz(-1.0556442) q[3];
sx q[3];
rz(0.011659415) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.0164531) q[2];
sx q[2];
rz(-0.21802248) q[2];
sx q[2];
rz(1.3583292) q[2];
rz(-2.6839117) q[3];
sx q[3];
rz(-1.4384559) q[3];
sx q[3];
rz(-0.55041909) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0871171) q[0];
sx q[0];
rz(-3.1396907) q[0];
sx q[0];
rz(0.052363366) q[0];
rz(0.26854435) q[1];
sx q[1];
rz(-1.152055) q[1];
sx q[1];
rz(2.9103738) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5729882) q[0];
sx q[0];
rz(-1.7460965) q[0];
sx q[0];
rz(2.8820287) q[0];
rz(-pi) q[1];
rz(0.32842689) q[2];
sx q[2];
rz(-2.3691485) q[2];
sx q[2];
rz(1.7945811) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.0655122) q[1];
sx q[1];
rz(-1.8022746) q[1];
sx q[1];
rz(-0.3409909) q[1];
rz(-pi) q[2];
rz(1.3116501) q[3];
sx q[3];
rz(-0.62302019) q[3];
sx q[3];
rz(1.3192605) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.3117567) q[2];
sx q[2];
rz(-2.629403) q[2];
sx q[2];
rz(0.83489418) q[2];
rz(3.0541776) q[3];
sx q[3];
rz(-1.089047) q[3];
sx q[3];
rz(-1.0138698) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8304623) q[0];
sx q[0];
rz(-2.1489547) q[0];
sx q[0];
rz(-0.8570689) q[0];
rz(-2.3509707) q[1];
sx q[1];
rz(-0.010684314) q[1];
sx q[1];
rz(-2.0649921) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.83504377) q[0];
sx q[0];
rz(-1.8721415) q[0];
sx q[0];
rz(-2.8623566) q[0];
rz(-2.1250379) q[2];
sx q[2];
rz(-0.86807251) q[2];
sx q[2];
rz(-1.6464485) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.5416661) q[1];
sx q[1];
rz(-1.8549207) q[1];
sx q[1];
rz(0.78943531) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4348599) q[3];
sx q[3];
rz(-1.9099209) q[3];
sx q[3];
rz(1.1992169) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.87535805) q[2];
sx q[2];
rz(-0.42576063) q[2];
sx q[2];
rz(0.88179624) q[2];
rz(-1.4074696) q[3];
sx q[3];
rz(-0.93972814) q[3];
sx q[3];
rz(0.56434977) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28500104) q[0];
sx q[0];
rz(-3.1398616) q[0];
sx q[0];
rz(-0.28975394) q[0];
rz(1.2334791) q[1];
sx q[1];
rz(-1.2954442) q[1];
sx q[1];
rz(-0.27898702) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10418071) q[0];
sx q[0];
rz(-1.54561) q[0];
sx q[0];
rz(1.7392147) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.40172949) q[2];
sx q[2];
rz(-1.3082983) q[2];
sx q[2];
rz(-2.0892734) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.853914) q[1];
sx q[1];
rz(-0.33386896) q[1];
sx q[1];
rz(-1.7823269) q[1];
x q[2];
rz(0.85106982) q[3];
sx q[3];
rz(-1.9033696) q[3];
sx q[3];
rz(-0.50902396) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.4616619) q[2];
sx q[2];
rz(-1.8671904) q[2];
sx q[2];
rz(0.28110853) q[2];
rz(-2.5680961) q[3];
sx q[3];
rz(-2.1888013) q[3];
sx q[3];
rz(-2.7219685) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9004274) q[0];
sx q[0];
rz(-0.92210162) q[0];
sx q[0];
rz(-1.7036555) q[0];
rz(-2.9519713) q[1];
sx q[1];
rz(-0.01556839) q[1];
sx q[1];
rz(1.9059034) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3251117) q[0];
sx q[0];
rz(-0.50704038) q[0];
sx q[0];
rz(-1.1757686) q[0];
rz(-pi) q[1];
rz(0.81053712) q[2];
sx q[2];
rz(-2.166368) q[2];
sx q[2];
rz(3.0308804) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.6744102) q[1];
sx q[1];
rz(-1.4100209) q[1];
sx q[1];
rz(1.2700446) q[1];
rz(-pi) q[2];
rz(0.017589324) q[3];
sx q[3];
rz(-1.5776199) q[3];
sx q[3];
rz(-1.9987035) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.8427061) q[2];
sx q[2];
rz(-1.3624531) q[2];
sx q[2];
rz(2.5617808) q[2];
rz(1.7719571) q[3];
sx q[3];
rz(-0.42391351) q[3];
sx q[3];
rz(0.4637318) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.46691) q[0];
sx q[0];
rz(-1.6413825) q[0];
sx q[0];
rz(1.7036194) q[0];
rz(-1.2314318) q[1];
sx q[1];
rz(-0.91130251) q[1];
sx q[1];
rz(-1.4968754) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9195596) q[0];
sx q[0];
rz(-2.3565533) q[0];
sx q[0];
rz(1.9505984) q[0];
rz(2.4980372) q[2];
sx q[2];
rz(-2.5179259) q[2];
sx q[2];
rz(2.8306431) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.4675967) q[1];
sx q[1];
rz(-1.5679545) q[1];
sx q[1];
rz(3.1238822) q[1];
rz(-pi) q[2];
rz(-2.1040099) q[3];
sx q[3];
rz(-0.54736558) q[3];
sx q[3];
rz(-1.7302046) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.9897495) q[2];
sx q[2];
rz(-1.5785549) q[2];
sx q[2];
rz(0.49636504) q[2];
rz(-2.1967891) q[3];
sx q[3];
rz(-3.1365518) q[3];
sx q[3];
rz(2.4387824) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4934568) q[0];
sx q[0];
rz(-1.428816) q[0];
sx q[0];
rz(-1.1270123) q[0];
rz(-1.5588749) q[1];
sx q[1];
rz(-2.8235148) q[1];
sx q[1];
rz(-2.943218) q[1];
rz(1.4737829) q[2];
sx q[2];
rz(-0.2578659) q[2];
sx q[2];
rz(0.21543287) q[2];
rz(2.583523) q[3];
sx q[3];
rz(-1.2574099) q[3];
sx q[3];
rz(-1.3539061) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
