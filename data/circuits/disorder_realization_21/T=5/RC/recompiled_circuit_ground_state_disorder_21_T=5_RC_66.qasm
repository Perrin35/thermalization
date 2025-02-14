OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.2668827) q[0];
sx q[0];
rz(-1.0615791) q[0];
sx q[0];
rz(-0.51577407) q[0];
rz(0.90142673) q[1];
sx q[1];
rz(-0.20144784) q[1];
sx q[1];
rz(2.1405061) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7667215) q[0];
sx q[0];
rz(-2.3618638) q[0];
sx q[0];
rz(-0.73260247) q[0];
rz(-pi) q[1];
rz(1.5770677) q[2];
sx q[2];
rz(-1.1561818) q[2];
sx q[2];
rz(2.8635777) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.6382003) q[1];
sx q[1];
rz(-0.93825785) q[1];
sx q[1];
rz(-1.0122385) q[1];
x q[2];
rz(1.4489861) q[3];
sx q[3];
rz(-1.7171905) q[3];
sx q[3];
rz(3.0288689) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.72656816) q[2];
sx q[2];
rz(-2.6308306) q[2];
sx q[2];
rz(1.4600352) q[2];
rz(0.36871746) q[3];
sx q[3];
rz(-1.5867149) q[3];
sx q[3];
rz(-1.6603498) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6129795) q[0];
sx q[0];
rz(-3.0280085) q[0];
sx q[0];
rz(-3.0533277) q[0];
rz(-0.9322235) q[1];
sx q[1];
rz(-2.014522) q[1];
sx q[1];
rz(-0.50311911) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4012484) q[0];
sx q[0];
rz(-0.49763864) q[0];
sx q[0];
rz(0.21792441) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2406949) q[2];
sx q[2];
rz(-1.7187727) q[2];
sx q[2];
rz(2.072538) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.14993653) q[1];
sx q[1];
rz(-2.7767065) q[1];
sx q[1];
rz(-2.8349692) q[1];
rz(-pi) q[2];
rz(-0.52946584) q[3];
sx q[3];
rz(-0.90116167) q[3];
sx q[3];
rz(0.96082276) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.8253537) q[2];
sx q[2];
rz(-2.1672858) q[2];
sx q[2];
rz(0.45315722) q[2];
rz(3.1406) q[3];
sx q[3];
rz(-1.5788635) q[3];
sx q[3];
rz(-0.7411952) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2082625) q[0];
sx q[0];
rz(-2.9587511) q[0];
sx q[0];
rz(-0.46874794) q[0];
rz(2.5543429) q[1];
sx q[1];
rz(-1.8533555) q[1];
sx q[1];
rz(-0.25150484) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55401459) q[0];
sx q[0];
rz(-1.16246) q[0];
sx q[0];
rz(0.46550444) q[0];
rz(-pi) q[1];
rz(-1.3590888) q[2];
sx q[2];
rz(-1.9518753) q[2];
sx q[2];
rz(0.43109387) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.9589825) q[1];
sx q[1];
rz(-1.512456) q[1];
sx q[1];
rz(-1.6946409) q[1];
rz(-pi) q[2];
rz(2.4819314) q[3];
sx q[3];
rz(-0.53852496) q[3];
sx q[3];
rz(-2.7646567) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.59715366) q[2];
sx q[2];
rz(-0.73828283) q[2];
sx q[2];
rz(-3.0944589) q[2];
rz(-2.2733287) q[3];
sx q[3];
rz(-1.2406113) q[3];
sx q[3];
rz(-0.49938437) q[3];
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
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.353001) q[0];
sx q[0];
rz(-2.7356739) q[0];
sx q[0];
rz(1.3141919) q[0];
rz(-2.3315673) q[1];
sx q[1];
rz(-1.6255197) q[1];
sx q[1];
rz(-2.8257418) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0994713) q[0];
sx q[0];
rz(-1.0183304) q[0];
sx q[0];
rz(-0.61037678) q[0];
rz(-pi) q[1];
x q[1];
rz(0.48607488) q[2];
sx q[2];
rz(-2.7856084) q[2];
sx q[2];
rz(-2.4497368) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.1882945) q[1];
sx q[1];
rz(-0.85894924) q[1];
sx q[1];
rz(-1.9602592) q[1];
x q[2];
rz(0.21574919) q[3];
sx q[3];
rz(-0.61601725) q[3];
sx q[3];
rz(2.8293138) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.4874068) q[2];
sx q[2];
rz(-2.5155289) q[2];
sx q[2];
rz(-0.89228863) q[2];
rz(1.4687294) q[3];
sx q[3];
rz(-1.059831) q[3];
sx q[3];
rz(1.0631801) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6650218) q[0];
sx q[0];
rz(-0.63711089) q[0];
sx q[0];
rz(0.88660216) q[0];
rz(-0.84238148) q[1];
sx q[1];
rz(-1.3096755) q[1];
sx q[1];
rz(-2.0319891) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3583325) q[0];
sx q[0];
rz(-1.9985804) q[0];
sx q[0];
rz(-0.88029998) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.38903) q[2];
sx q[2];
rz(-0.85585844) q[2];
sx q[2];
rz(0.37059298) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.2354616) q[1];
sx q[1];
rz(-1.0603728) q[1];
sx q[1];
rz(-0.6554558) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.97115406) q[3];
sx q[3];
rz(-2.1622938) q[3];
sx q[3];
rz(2.9324071) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.2716918) q[2];
sx q[2];
rz(-0.39125189) q[2];
sx q[2];
rz(-0.58763495) q[2];
rz(0.21688004) q[3];
sx q[3];
rz(-1.8606595) q[3];
sx q[3];
rz(-1.3542401) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.13123913) q[0];
sx q[0];
rz(-0.47061798) q[0];
sx q[0];
rz(-1.3838029) q[0];
rz(-2.9673987) q[1];
sx q[1];
rz(-1.1015588) q[1];
sx q[1];
rz(1.1879638) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39858946) q[0];
sx q[0];
rz(-1.849353) q[0];
sx q[0];
rz(-0.085531959) q[0];
rz(-pi) q[1];
x q[1];
rz(0.30004487) q[2];
sx q[2];
rz(-2.6159673) q[2];
sx q[2];
rz(2.2208461) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.7158051) q[1];
sx q[1];
rz(-1.5400229) q[1];
sx q[1];
rz(1.6795067) q[1];
rz(-1.8781885) q[3];
sx q[3];
rz(-2.3114021) q[3];
sx q[3];
rz(-2.8678558) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.69372988) q[2];
sx q[2];
rz(-1.8665946) q[2];
sx q[2];
rz(-2.0873783) q[2];
rz(0.55274719) q[3];
sx q[3];
rz(-0.25944969) q[3];
sx q[3];
rz(-1.118008) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.066910557) q[0];
sx q[0];
rz(-2.5801165) q[0];
sx q[0];
rz(0.018360227) q[0];
rz(-1.8439937) q[1];
sx q[1];
rz(-1.8074139) q[1];
sx q[1];
rz(2.0265354) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0783065) q[0];
sx q[0];
rz(-3.0111599) q[0];
sx q[0];
rz(3.0060769) q[0];
rz(-pi) q[1];
rz(2.0231904) q[2];
sx q[2];
rz(-1.8542737) q[2];
sx q[2];
rz(0.3800791) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.74919486) q[1];
sx q[1];
rz(-2.1304211) q[1];
sx q[1];
rz(2.1361875) q[1];
rz(-0.84172084) q[3];
sx q[3];
rz(-0.63243619) q[3];
sx q[3];
rz(3.1068592) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.6899507) q[2];
sx q[2];
rz(-0.74863282) q[2];
sx q[2];
rz(-0.79317036) q[2];
rz(-2.4801109) q[3];
sx q[3];
rz(-1.7659148) q[3];
sx q[3];
rz(-2.49559) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20872214) q[0];
sx q[0];
rz(-1.1247617) q[0];
sx q[0];
rz(-2.9710508) q[0];
rz(2.1216682) q[1];
sx q[1];
rz(-2.7248757) q[1];
sx q[1];
rz(0.65933093) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.19457283) q[0];
sx q[0];
rz(-1.5156557) q[0];
sx q[0];
rz(1.5558262) q[0];
rz(2.3691828) q[2];
sx q[2];
rz(-2.4669224) q[2];
sx q[2];
rz(1.9325369) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.48261343) q[1];
sx q[1];
rz(-0.70090997) q[1];
sx q[1];
rz(0.47710508) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.002568) q[3];
sx q[3];
rz(-2.2428277) q[3];
sx q[3];
rz(-2.6300485) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.3063804) q[2];
sx q[2];
rz(-2.7472718) q[2];
sx q[2];
rz(2.1739056) q[2];
rz(-2.5505998) q[3];
sx q[3];
rz(-1.3700181) q[3];
sx q[3];
rz(1.7041357) q[3];
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
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10385926) q[0];
sx q[0];
rz(-2.2028706) q[0];
sx q[0];
rz(-2.9528604) q[0];
rz(1.0602779) q[1];
sx q[1];
rz(-2.4791368) q[1];
sx q[1];
rz(-0.8998543) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8322139) q[0];
sx q[0];
rz(-0.77118528) q[0];
sx q[0];
rz(-0.82892771) q[0];
rz(-pi) q[1];
rz(0.84328358) q[2];
sx q[2];
rz(-1.5727864) q[2];
sx q[2];
rz(-0.28283027) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.65727216) q[1];
sx q[1];
rz(-1.5690227) q[1];
sx q[1];
rz(-0.65048762) q[1];
x q[2];
rz(2.9505355) q[3];
sx q[3];
rz(-0.85874288) q[3];
sx q[3];
rz(-2.4993757) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.4047644) q[2];
sx q[2];
rz(-1.9177723) q[2];
sx q[2];
rz(1.9723816) q[2];
rz(-2.2660008) q[3];
sx q[3];
rz(-2.4647522) q[3];
sx q[3];
rz(-0.08918795) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.635023) q[0];
sx q[0];
rz(-1.6825786) q[0];
sx q[0];
rz(-2.7154679) q[0];
rz(1.2395073) q[1];
sx q[1];
rz(-1.0303717) q[1];
sx q[1];
rz(0.32136163) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.118537) q[0];
sx q[0];
rz(-1.5154953) q[0];
sx q[0];
rz(2.1441516) q[0];
rz(1.4570792) q[2];
sx q[2];
rz(-0.4452332) q[2];
sx q[2];
rz(1.6086827) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.1758986) q[1];
sx q[1];
rz(-1.7135317) q[1];
sx q[1];
rz(-1.9121777) q[1];
rz(2.0780744) q[3];
sx q[3];
rz(-0.92417704) q[3];
sx q[3];
rz(0.77669981) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.1539803) q[2];
sx q[2];
rz(-0.93903792) q[2];
sx q[2];
rz(-0.067848094) q[2];
rz(-2.1765354) q[3];
sx q[3];
rz(-0.32876757) q[3];
sx q[3];
rz(-1.4062101) q[3];
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
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6168552) q[0];
sx q[0];
rz(-0.4586093) q[0];
sx q[0];
rz(-2.1427857) q[0];
rz(-0.85159341) q[1];
sx q[1];
rz(-2.0709745) q[1];
sx q[1];
rz(-2.4877683) q[1];
rz(0.10112496) q[2];
sx q[2];
rz(-2.58287) q[2];
sx q[2];
rz(1.9091633) q[2];
rz(0.028164955) q[3];
sx q[3];
rz(-0.71816254) q[3];
sx q[3];
rz(-0.95550334) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
