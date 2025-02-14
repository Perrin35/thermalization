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
rz(-2.9838188) q[0];
sx q[0];
rz(4.3133419) q[0];
sx q[0];
rz(11.316909) q[0];
rz(-0.14021048) q[1];
sx q[1];
rz(-1.6970716) q[1];
sx q[1];
rz(0.061847774) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0842347) q[0];
sx q[0];
rz(-1.0751197) q[0];
sx q[0];
rz(0.88229124) q[0];
rz(3.05117) q[2];
sx q[2];
rz(-1.5704463) q[2];
sx q[2];
rz(-1.4594913) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.5420336) q[1];
sx q[1];
rz(-2.2890511) q[1];
sx q[1];
rz(0.12964779) q[1];
rz(-2.2592945) q[3];
sx q[3];
rz(-1.2475444) q[3];
sx q[3];
rz(2.0417449) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.7141815) q[2];
sx q[2];
rz(-2.2990871) q[2];
sx q[2];
rz(-0.25644914) q[2];
rz(-2.9668258) q[3];
sx q[3];
rz(-1.9758965) q[3];
sx q[3];
rz(0.83785653) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8137708) q[0];
sx q[0];
rz(-2.796266) q[0];
sx q[0];
rz(-3.0178965) q[0];
rz(-2.9028614) q[1];
sx q[1];
rz(-2.3614466) q[1];
sx q[1];
rz(-0.10496584) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2793113) q[0];
sx q[0];
rz(-2.263592) q[0];
sx q[0];
rz(-0.19578085) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3388073) q[2];
sx q[2];
rz(-1.1646484) q[2];
sx q[2];
rz(-0.37877235) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.2626896) q[1];
sx q[1];
rz(-2.5715552) q[1];
sx q[1];
rz(-1.0703342) q[1];
rz(-pi) q[2];
rz(2.4363748) q[3];
sx q[3];
rz(-1.8466443) q[3];
sx q[3];
rz(-2.5666756) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.77439848) q[2];
sx q[2];
rz(-1.4906733) q[2];
sx q[2];
rz(1.4614089) q[2];
rz(-1.2857619) q[3];
sx q[3];
rz(-1.6971842) q[3];
sx q[3];
rz(2.8694966) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9937781) q[0];
sx q[0];
rz(-2.9312134) q[0];
sx q[0];
rz(-1.0149957) q[0];
rz(2.2442832) q[1];
sx q[1];
rz(-1.6474479) q[1];
sx q[1];
rz(-0.6764594) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38329002) q[0];
sx q[0];
rz(-1.5448017) q[0];
sx q[0];
rz(0.95375188) q[0];
rz(-pi) q[1];
x q[1];
rz(0.93134201) q[2];
sx q[2];
rz(-0.73340511) q[2];
sx q[2];
rz(-0.99898192) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.3520917) q[1];
sx q[1];
rz(-1.772023) q[1];
sx q[1];
rz(1.1835062) q[1];
rz(-pi) q[2];
rz(-0.018142975) q[3];
sx q[3];
rz(-2.3251495) q[3];
sx q[3];
rz(1.7545561) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.93495381) q[2];
sx q[2];
rz(-2.166344) q[2];
sx q[2];
rz(2.3812531) q[2];
rz(1.2065411) q[3];
sx q[3];
rz(-2.3307266) q[3];
sx q[3];
rz(2.2981203) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7059785) q[0];
sx q[0];
rz(-0.62788457) q[0];
sx q[0];
rz(1.5950369) q[0];
rz(-2.4226923) q[1];
sx q[1];
rz(-1.8214106) q[1];
sx q[1];
rz(0.23153201) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2255428) q[0];
sx q[0];
rz(-1.0576436) q[0];
sx q[0];
rz(1.0731359) q[0];
rz(2.3489352) q[2];
sx q[2];
rz(-1.5735627) q[2];
sx q[2];
rz(-2.1994928) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-3.0147103) q[1];
sx q[1];
rz(-1.0877123) q[1];
sx q[1];
rz(-1.4654069) q[1];
x q[2];
rz(-1.9636964) q[3];
sx q[3];
rz(-1.0024973) q[3];
sx q[3];
rz(1.6015805) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.4270758) q[2];
sx q[2];
rz(-0.38893739) q[2];
sx q[2];
rz(0.60849774) q[2];
rz(2.9454339) q[3];
sx q[3];
rz(-1.720263) q[3];
sx q[3];
rz(2.1430446) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-1.8683559) q[0];
sx q[0];
rz(-0.65160692) q[0];
sx q[0];
rz(0.77907816) q[0];
rz(-0.92998663) q[1];
sx q[1];
rz(-1.9735186) q[1];
sx q[1];
rz(-1.9680061) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.54544696) q[0];
sx q[0];
rz(-1.2327502) q[0];
sx q[0];
rz(-0.0041739504) q[0];
rz(-pi) q[1];
rz(-2.1861325) q[2];
sx q[2];
rz(-1.460607) q[2];
sx q[2];
rz(-2.2249391) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.0231578) q[1];
sx q[1];
rz(-1.506532) q[1];
sx q[1];
rz(1.8885018) q[1];
rz(-pi) q[2];
rz(0.42901943) q[3];
sx q[3];
rz(-1.3976421) q[3];
sx q[3];
rz(2.8575983) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.32213) q[2];
sx q[2];
rz(-2.9746015) q[2];
sx q[2];
rz(-2.4665311) q[2];
rz(2.2760462) q[3];
sx q[3];
rz(-1.5138488) q[3];
sx q[3];
rz(-1.2077829) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.753767) q[0];
sx q[0];
rz(-3.1238811) q[0];
sx q[0];
rz(-2.2702763) q[0];
rz(-2.9246092) q[1];
sx q[1];
rz(-1.5996409) q[1];
sx q[1];
rz(1.8035696) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9674613) q[0];
sx q[0];
rz(-2.2207157) q[0];
sx q[0];
rz(-1.6748669) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7291405) q[2];
sx q[2];
rz(-2.4096903) q[2];
sx q[2];
rz(-0.5754234) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.034778193) q[1];
sx q[1];
rz(-2.9902774) q[1];
sx q[1];
rz(-1.2447912) q[1];
rz(-0.35049482) q[3];
sx q[3];
rz(-0.51096254) q[3];
sx q[3];
rz(-0.048633752) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.013082144) q[2];
sx q[2];
rz(-1.6435813) q[2];
sx q[2];
rz(-2.3415372) q[2];
rz(-2.1498146) q[3];
sx q[3];
rz(-2.1938775) q[3];
sx q[3];
rz(0.94314027) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
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
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8868788) q[0];
sx q[0];
rz(-2.7006221) q[0];
sx q[0];
rz(-0.15175858) q[0];
rz(1.4166547) q[1];
sx q[1];
rz(-0.85985008) q[1];
sx q[1];
rz(-2.3053665) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1372414) q[0];
sx q[0];
rz(-0.81295952) q[0];
sx q[0];
rz(-1.0197958) q[0];
x q[1];
rz(-1.5791527) q[2];
sx q[2];
rz(-2.1474693) q[2];
sx q[2];
rz(1.4954665) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.6265833) q[1];
sx q[1];
rz(-1.3224416) q[1];
sx q[1];
rz(0.088661389) q[1];
rz(0.46549972) q[3];
sx q[3];
rz(-1.6878078) q[3];
sx q[3];
rz(-2.4586364) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.4862711) q[2];
sx q[2];
rz(-0.063491193) q[2];
sx q[2];
rz(0.07902321) q[2];
rz(-0.71845734) q[3];
sx q[3];
rz(-1.6301165) q[3];
sx q[3];
rz(2.2169936) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.786161) q[0];
sx q[0];
rz(-2.5434255) q[0];
sx q[0];
rz(2.2736736) q[0];
rz(2.4527841) q[1];
sx q[1];
rz(-0.30615607) q[1];
sx q[1];
rz(2.205663) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0999231) q[0];
sx q[0];
rz(-1.3893621) q[0];
sx q[0];
rz(-0.52930852) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8008119) q[2];
sx q[2];
rz(-1.182297) q[2];
sx q[2];
rz(2.4430371) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.28081255) q[1];
sx q[1];
rz(-1.6242923) q[1];
sx q[1];
rz(-1.2455275) q[1];
rz(-1.5722365) q[3];
sx q[3];
rz(-0.82732302) q[3];
sx q[3];
rz(-1.8515406) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.9800637) q[2];
sx q[2];
rz(-2.4825725) q[2];
sx q[2];
rz(2.8847983) q[2];
rz(2.1234546) q[3];
sx q[3];
rz(-2.0223821) q[3];
sx q[3];
rz(-2.170678) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37860206) q[0];
sx q[0];
rz(-2.584223) q[0];
sx q[0];
rz(-0.85365224) q[0];
rz(1.8866106) q[1];
sx q[1];
rz(-1.9957142) q[1];
sx q[1];
rz(-2.8010211) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4494394) q[0];
sx q[0];
rz(-1.3763577) q[0];
sx q[0];
rz(0.16325918) q[0];
x q[1];
rz(1.8179632) q[2];
sx q[2];
rz(-0.91193141) q[2];
sx q[2];
rz(0.62550046) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.53696991) q[1];
sx q[1];
rz(-0.83965644) q[1];
sx q[1];
rz(-2.5174058) q[1];
rz(-pi) q[2];
rz(-2.8664175) q[3];
sx q[3];
rz(-2.0180297) q[3];
sx q[3];
rz(-2.3768611) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.42369947) q[2];
sx q[2];
rz(-1.7607949) q[2];
sx q[2];
rz(2.763486) q[2];
rz(0.2374436) q[3];
sx q[3];
rz(-2.0707371) q[3];
sx q[3];
rz(2.8606991) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
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
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.045192748) q[0];
sx q[0];
rz(-1.0843596) q[0];
sx q[0];
rz(-2.4943446) q[0];
rz(1.9735533) q[1];
sx q[1];
rz(-2.2868575) q[1];
sx q[1];
rz(2.7580269) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7591568) q[0];
sx q[0];
rz(-1.54017) q[0];
sx q[0];
rz(1.5923862) q[0];
rz(0.53159376) q[2];
sx q[2];
rz(-1.2832912) q[2];
sx q[2];
rz(0.93898857) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.4216363) q[1];
sx q[1];
rz(-1.7767748) q[1];
sx q[1];
rz(0.27011807) q[1];
rz(0.10730073) q[3];
sx q[3];
rz(-1.0578053) q[3];
sx q[3];
rz(-1.9657621) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.8770404) q[2];
sx q[2];
rz(-2.2525747) q[2];
sx q[2];
rz(-1.5846579) q[2];
rz(-1.5133739) q[3];
sx q[3];
rz(-1.7729365) q[3];
sx q[3];
rz(-0.42678601) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3601111) q[0];
sx q[0];
rz(-1.7088912) q[0];
sx q[0];
rz(0.29722469) q[0];
rz(1.9602641) q[1];
sx q[1];
rz(-1.7971296) q[1];
sx q[1];
rz(2.8167579) q[1];
rz(-2.8988373) q[2];
sx q[2];
rz(-1.1052255) q[2];
sx q[2];
rz(0.33585264) q[2];
rz(0.48458002) q[3];
sx q[3];
rz(-1.5427586) q[3];
sx q[3];
rz(1.2324738) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
